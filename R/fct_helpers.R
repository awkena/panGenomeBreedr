#' Helper Functions for PanGenomeBreedr Application
#'
#' This file contains various helper functions used throughout the application.
#' Functions are organized by their purpose and usage context.
#' None of these functions are exported as they're meant for internal use only.

#-----------------------------------------------------------------------------
# VCF File Processing Functions (KASP Marker design section)
#-----------------------------------------------------------------------------
#' Extract Marker IDs and Chromosome IDs from VCF File
#'
#' This function reads a VCF file and extracts all unique marker IDs and chromosome IDs.
#'
#' @param vcf_path Character string specifying the path to the VCF file.
#'
#' @return A list with two elements:
#'   \item{vcf_matrix_markerID}{Character vector of unique marker IDs}
#'   \item{vcf_matrix_chromID}{Character vector of unique chromosome IDs}
#'
#' @details The function uses the vcfR package to read the VCF file and extracts the
#' marker IDs from the ID field and chromosome identifiers from the CHROM field in the VCF.
#' @noRd
#' @importFrom vcfR read.vcfR
marker.chr_ID <- function(vcf_path) {
  vcf_data <- vcfR::read.vcfR(vcf_path, verbose = FALSE)

  vcf_matrix <- as.matrix(vcf_data@fix) |> as.data.frame()

  vcf_matrix_chromID <- unique(vcf_matrix[["CHROM"]])
  vcf_matrix_markerID <- unique(vcf_matrix[["ID"]])

  result_list <- list(
    vcf_matrix_markerID = vcf_matrix_markerID,
    vcf_matrix_chromID = vcf_matrix_chromID
  )

  return(result_list)
}

#-----------------------------------------------------------------------------
# KASP Data Processing Functions (Marker validation section)
#-----------------------------------------------------------------------------

#' Add plates column with unique IDs
#'
#' This function adds plates column with unique ID's for kasp data file if not present.
#' This is later passed to the `nsamples_plate()` to obtain plate summary.
#'
#' @param x A data frame containing the data component from the read kasp file
#' @return A data frame with plate column added if not present
#'
#' @noRd
#' @importFrom data.table copy
plates_col <- function(x) {
  copied_data <- data.table::copy(x) # Copy the dataset

  colnames(copied_data) <- tolower(colnames(copied_data))

  # check if 'MasterPlate' and 'SNPID' exist
  if (!all(c("MasterPlate", "SNPID") %in% colnames(x))) {
    stop("Error: Columns 'MasterPlate' and 'SNPID' are missing in the input data.")
  }

  # Check if plates column exists
  if (!"plates" %in% colnames(x)) {
    x$plates <- paste(x$MasterPlate, x$SNPID, sep = "_")
  }

  return(x)
}

#' Get unique plates IDs
#'
#' This function extracts unique plates when plates column is added or already exists in Data
#'
#' @param x Results from read_kasp_csv()
#' @return A vector of unique plates IDs
#'
#' @noRd
uniq_plates <- function(x) {
  x_in <- plates_col(x)
  if (!"plates" %in% colnames(x_in)) {
    stop("Error: 'plates' column missing after applying plates_col().")
  }

  uniq_plate <- unique(x_in[["plates"]])

  return(uniq_plate) # return unique plate as a vector.
}

#' Generate a kasp marker map.
#'
#' This function parses SNPID strings and extracts chromosome and position details
#' and creates a mapping data frame
#' @param kasp_data a data frame of KASP genotype calls for one or mulitiple plates
#'
#' @returns a data frame for kasp marker map
#'
#' @noRd
#' @importFrom stringr str_match
kasp_marker_map <- function(kasp_data) {
  # Validate input
  if (is.data.frame(kasp_data) != TRUE) {
    stop("Kasp_data must be a data frame")
  }

  # Confirm required column
  if (!("SNPID" %in% colnames(kasp_data))) {
    stop("Kasp_data must have column: SNPID")
  }

  # Index required column(SNPID)
  SNP_ids <- kasp_data[["SNPID"]] |> unique()

  # Get length of SNP_ids
  SNP_ids_length <- length(SNP_ids)

  # Initialize result dataframe
  kasp_map <- data.frame(
    SNPID = character(SNP_ids_length),
    SNP_ID2 = character(SNP_ids_length),
    chr = numeric(SNP_ids_length),
    pos = integer(SNP_ids_length),
    stringsAsFactors = FALSE
  )

  # for loop to generate and insert
  for (variable in seq_len(SNP_ids_length)) {
    # Insert value into snpid column
    kasp_map[variable, "SNPID"] <- SNP_ids[variable]

    # Insert value into snpid_2
    #- Split and format first
    formatted_ids <- stringr::str_match(
      string = SNP_ids[variable],
      pattern = "([a-zA-Z]+)(\\d+)"
    ) |> as.vector()

    #- if statement to check and validate, also to crop of..
    if (SNP_ids[variable] %in% formatted_ids) {
      aa <- formatted_ids[formatted_ids != SNP_ids[variable]]

      kasp_map[variable, "SNP_ID2"] <- paste(aa[1], aa[2], sep = "|")
    } else {
      kasp_map[variable, "SNP_ID2"] <- paste(formatted_ids[1], formatted_ids[2], sep = "|") # inserted
    }

    # Insert value into chr column.
    kasp_map[variable, "chr"] <- variable

    # Insert position
    kasp_map[variable, "pos"] <- formatted_ids[length(formatted_ids)] |> as.integer()
  }

  return(kasp_map)
}

#' Extract genotype calls
#'
#' This function extracts all the genotype calls within the plates column
#'
#' @param x a list object processed by read_kasp_csv()
#' @param a an ID within the plates column
#'
#' @returns a vector of all genotype calls
#'
#' @noRd
get_calls <- function(x, a) {
  x <- plates_col(x)

  if (!"plates" %in% colnames(x)) {
    stop("Error: 'plates' column is missing after applying plates_col().")
  }

  if (!"Call" %in% colnames(x)) {
    stop("Error: 'Call' column is missing in the dataset.")
  }

  x_filtered <- x[x[["plates"]] == a, ]

  if (nrow(x_filtered) == 0) {
    return(NULL) # Return NULL instead of causing a subscript error
  }

  return(x_filtered[["Call"]]) # Returns call column as a vector
}

#' Get colnames from KASP Data
#'
#' @param x a data frame of KASP genotype calls for one or multiple plates
#'
#' @returns a character vector
#' @noRd
col_names.data <- function(x) {
  # Validate input
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("Input must be a data frame or matrix")
  }

  if (is.null(x) || nrow(x) == 0) {
    return(character(0)) # Return an empty character vector
  }

  return(colnames(x))
}




#' Convert genotype list to data frame
#'
#' This function takes a list of genotypes (typically from get_alleles()) and
#' converts it into a tabular data frame format for easier viewing and analysis.
#'
#' @param x List object with genotypes subsetted after get_alleles()
#' @return A data frame containing the genotype data with appropriate column names
#'
#' @noRd
genotypes <- function(x) {
  # Validate input
  if (!is.vector(x) && !is.list(x)) {
    stop("Input must be a vector or list")
  }

  # Unlist if needed
  x <- unlist(x)

  if (length(x) == 0) {
    warning("Input is empty")
    return(data.frame())
  }

  # Create a data frame with named columns
  result <- as.data.frame(matrix(
    nrow = 1, ncol = length(x),
    dimnames = list(NULL, names(x))
  ))

  # Assign all values at once
  result[1, ] <- x

  return(result)
}


#' Convert allele object to data frame
#'
#' This function transforms an allele object returned by get_alleles() into a
#' structured data frame format, with columns labeled as 'allele_A', 'allele_B', etc.
#'
#' @param x List object containing alleles from get_alleles()
#' @return A data frame with columns for each detected allele
#'
#' @noRd
alleles_df <- function(x) {
  # Validate input
  if (is.null(x)) {
    stop("Input cannot be NULL")
  }

  x <- unlist(x)

  if (length(x) == 0) {
    warning("No alleles found")
    return(data.frame())
  }

  # Create allele column names
  allele_cols <- paste("allele", LETTERS[1:length(x)], sep = "_")

  # Create and populate data frame
  df <- as.data.frame(matrix(
    nrow = 1,
    ncol = length(x),
    dimnames = list(NULL, allele_cols)
  ))
  df[1, ] <- x

  return(df)
}

#' Generate summary statistics for KASP genotyping plates
#'
#' This function analyzes KASP genotyping data and produces a summary report
#' on the status of each plate, including allele counts, locus type classification,
#' and overall success status.
#'
#' @param x Data frame containing KASP genotyping data
#' @param subset Column name containing plate identifiers (default: 'MasterPlate')
#' @param sep Separator character used in genotype calls (default: ':')
#' @param geno_call Column name containing genotype calls (default: 'Call')
#' @return A data frame summarizing each plate with columns for plate ID, allele count,
#'
#' @noRd
kasp_color_stat.give <- function(x,
                                 subset = "MasterPlate",
                                 sep = ":",
                                 geno_call = "Call") {
  # Get unique plates
  plates <- x[, subset]
  plates_uniq <- unique(plates)
  nplates <- length(plates_uniq)

  # Create a dynamic data frame.
  dyn_df <- as.data.frame(matrix(nrow = nplates, ncol = 4, dimnames = list(NULL, c(subset, "Allele Count", "Locus Type", "Status"))))
  # Insert and populate.
  dyn_df[, 1] <- plates_uniq

  for (i in seq_len(nplates)) {
    # Subset each master plate
    master_plate <- x[plates == plates_uniq[i], ]

    # create a color vector based on KASP geno calls
    Color <- master_plate[[geno_call]]

    # Get alleles and possible genotypes using the `get_allele()` function
    alleles_geno <- get_alleles(
      x = Color,
      sep = sep,
      data_type = "kasp"
    )
    alleles <- alleles_geno$alleles

    dyn_df[i, 2] <- length(alleles)
    # if statements.
    if (length(alleles) == 2) {
      dyn_df[i, 3] <- "Bi-allelic"
      dyn_df[i, 4] <- "Successful"
    } else if (length(alleles) == 1) {
      dyn_df[i, 3] <- "Monomorphic"
      dyn_df[i, 4] <- "Successful"
    } else if (length(alleles) == 0) {
      dyn_df[i, 3] <- "None"
      dyn_df[i, 4] <- "Failed!"
    } else if (length(alleles) > 2) {
      dyn_df[i, 3] <- "Multi-allelic"
      dyn_df[i, 4] <- "Failed!"
    }
  }
  return(dyn_df)
}


#' Update column names in sample data frame
#'
#' This function dynamically renames columns in a sample data frame to match
#' specified identifiers for MasterPlate, SNP ID, and plate information.
#'
#' @param MasterPlate New name for the master plate column
#' @param SNIPID New name for the SNP ID column
#' @param PLATE New name for the plate column
#' @param data_nsamples Data frame containing sample information to be renamed
#' @return Data frame with updated column names
#'
#' @noRd
new_colnames_nsamples <- function(MasterPlate,
                                  SNIPID,
                                  PLATE,
                                  data_nsamples) {
  # Validate input
  if (!is.data.frame(data_nsamples)) {
    stop("Input must be a data frame")
  }

  if (ncol(data_nsamples) < 4) {
    stop("Input data frame must have at least 4 columns")
  }

  result <- as.data.frame(data_nsamples) # convert to dataframe

  # Update column names, preserving the 4th column name
  original_names <- colnames(result)
  new_names <- c(MasterPlate, SNIPID, PLATE, original_names[4:length(original_names)])
  colnames(result) <- new_names

  return(result)
}


#' Split user defined height and weight
#'
#' Converts a comma-separated string into a numeric vector of height and width values.
#' The input string should contain exactly two numeric values separated by a comma.
#'
#' @param string A character string of length 1 containing two comma-separated numeric values
#'
#' @return A numeric vector of length 2, where the first element is height and the second is width
#'
#' @noRd
height_width <- function(string){
  # Check if input is a character and length 1
  if(!is.character(string) || length(string) != 1) {
    stop("Input must be a single character string")
  }

  # Split the string by comma
  splitted <- strsplit(string, split = ',')[[1]]

  # Check if exactly two values were provided
  if(length(splitted) < 2) {
    stop("Input must contain exactly two values separated by a comma")
  }

  # Convert to numeric and handle potential NA values
  h_w_vector <- as.numeric(splitted)[1:2]

  if(any(is.na(h_w_vector))) {
    stop("Both values must be convertible to numbers")
  }

  return(h_w_vector)
}





################### Writing a bulk function to carry out decision support.
#' Extract Genotype Names for User Selection
#'
#' This function extracts genotype names from a specific batch in the dataset
#' to provide options for user selection of recurrent parent (rp) and donor
#' parent (dp) in downstream analysis.
#'
#' @param data Character. Path to the CSV file containing the genotype data.
#' @param Batch Numeric. The batch number from which to extract genotype names.
#'
#' @return Character vector. A vector containing the genotype names/IDs from
#'   the specified batch, suitable for use in dropdown menus or selection lists.
#'
#' @details The function performs the following operations:
#'   1. Reads the CSV file and converts to data frame
#'   2. Examines the first 100 columns to identify the start of genotype data
#'   3. Automatically detects "batch" and "genotype" columns (case-insensitive)
#'   4. Filters data to the specified batch
#'   5. Extracts and returns genotype names for that batch
#'
#'
#' @noRd
#'
Genotypes_user <- function(data, Batch) {
  read_data <- data
  # Get first 100 columns for possible extraction.
  first_hun_col <- read_data[1, 1:100]

  store_position <- which(first_hun_col %in% c("A", "G", "C", "T", NA))[1]

  ### Use grep function to get batch and genotypes.
  #- Batch
  batch <- grep("batch", x = colnames(first_hun_col), ignore.case = TRUE, value = TRUE)
  #- Genotype
  genotype <- grep("genotype", x = colnames(first_hun_col), ignore.case = TRUE, value = TRUE)

  # Return just the data column.
  meta_data <- read_data[read_data[batch] == Batch, 1:store_position]

  # Select column for genotypes
  just_geno <- meta_data[[genotype]]

  return(just_geno) # return genotypes of selected batch.
}

#' Generate Map File and Process KASP Data
#'
#' This function takes Agriplex data and processes it to generate a map file and
#' convert KASP genotype data to numeric format. It handles SNP marker parsing,
#' removes monomorphic markers, and processes parent-offspring data.
#'
#' @param data Character. Path to the CSV file containing the raw genotype data.
#' @param Batch Numeric. The batch number to filter and process from the dataset.
#' @param sep Character. Separator used in marker IDs for parsing (e.g., "_").
#' @param sep_2 Character. Separator used in genotype data (e.g., " / " for Agriplex).
#' @param data_type Character. Type of genotyping data (e.g., "agriplex").
#' @param rp Character. Name/ID of the first parent (recurrent parent).
#' @param dp Character. Name/ID of the second parent (donor parent).
#' @param geno_vec Character vector. Vector of genotype names/IDs in the dataset.
#' @param feedback Character. One of "yes" or "no". If "yes", uses existing map file;
#'   if "no", generates new map file from marker names.
#' @param mapfile_path Character. Path to existing map file (required when feedback = "yes").
#'
#' @return A list containing two elements:
#'   \item{mapfile}{Data frame containing the marker map information}
#'   \item{proc_kasp_f}{Processed and numeric-converted KASP data}
#'   \item{par_missing_dat}{loci with any parent missing genotype}
#'   \item{genotype_error}{SNP loci with potential genotype call errors}
#'   \item{parent_het}{Heterozygous parent genotypes}
#'
#' @details The function performs the following steps:
#'   1. Reads the input CSV file and identifies genotype data columns
#'   2. Filters data by the specified batch
#'   3. Removes monomorphic markers and markers with missing parent data
#'   4. Either generates a new map file or uses an existing one
#'   5. Processes KASP data and converts to numeric format
#'
#' @noRd
#'
#'

marker_file <- function(data = NULL,
                        Batch = NULL,
                        sep = NULL,
                        sep_2 = NULL,
                        data_type = NULL,
                        rp = NULL,
                        dp = NULL,
                        geno_vec = NULL,
                        feedback = c("yes", "no"),
                        na_code = -5,
                        mapfile_path = NULL) {
  read_data <-  data |> as.data.frame()
  # Get first 100 columns for possible extraction.
  first_hun_col <- read_data[1, 1:100]

  store_position <- which(first_hun_col %in% c("A", "G", "C", "T", NA))[1]

  #- Batch
  batch <- grep("batch", x = colnames(first_hun_col), ignore.case = TRUE, value = TRUE)
  #- Genotype
  genotype <- grep("genotype", x = colnames(first_hun_col), ignore.case = TRUE, value = TRUE)

  # Return just the data column.
  return_data <- read_data[read_data[batch] == Batch, store_position:ncol(read_data)]

  # Get colnames for first.
  col_first <- colnames(return_data)[[1]]

  Prefix <- substr(x = col_first, start = 1, stop = 1) # get prefix


  dp_val <- which(geno_vec %in% c(rp, dp)) # get row index of rp and dp

  # Get rownames of data.
  rownames(return_data) <- geno_vec

  # Remove monomorphs.
  no_mono_return_dat <- parent_poly(
    x = return_data,
    rp_row = dp_val[1],
    dp_row = dp_val[2],
    sep = sep_2
  )

  # Parent Present
  no_mono_return_data <- parent_missing(
    x = no_mono_return_dat,
    rp_row = dp_val[1],
    dp_row = dp_val[2]
  )$par_present

 # Parent Missing
  par_missing_dat <- parent_missing(
    x = no_mono_return_dat,
    rp_row = dp_val[1],
    dp_row = dp_val[2]
  )$par_missing

  # Genotype error.
  genotype_error <- geno_error(x = no_mono_return_dat ,
                               rp_row = dp_val[1],
                               dp_row = dp_val[2],
                               sep = sep_2,
                               data_type = data_type)$geno_err

  # Parent Hetero.
  parent_het <- parent_het(x = no_mono_return_dat ,
                           rp_row = dp_val[1],
                           dp_row = dp_val[2],
                           sep = sep_2,
                           na_code = na_code )$par_het


  snps <- colnames(no_mono_return_data)

  # Parse snpids to generate mapfile.
  # Based on feedback
  if (feedback == "no") {
    mapfile <- parse_marker_ns(x = snps, sep = sep, prefix = Prefix) |> order_markers()

    # # View(return_data)
    # parent_missing(x = return_data  ,rp_row = dp_val[1] ,dp_row =dp_val[2] )$par_missing |> View()
    # Process and convert to numeric
    proc_kasp_f <- proc_kasp(
      x = no_mono_return_data,
      sample_id = genotype,
      map_snp_id = "snpid",
      marker_start = 1,
      kasp_map = mapfile
    ) |>
      kasp_numeric(
        rp_row = dp_val[1],
        dp_row = dp_val[2],
        sep = sep_2,
        data_type = data_type
      )
  } else if (feedback == "yes") {
    mapfile <- read.csv(file = mapfile_path, stringsAsFactors = FALSE) |> as.data.frame()

    proc_kasp_f <- proc_kasp(
      x = no_mono_return_data,
      sample_id = "Genotype", # to be selected by user
      map_snp_id = "snpid", # to be selected by user
      marker_start = 1,
      kasp_map = mapfile
    ) |>
      kasp_numeric(
        rp_row = dp_val[1],
        dp_row = dp_val[2],
        sep = sep_2,
        data_type = data_type
      )
  }


  return(list(mapfile = mapfile,
              proc_kasp_f = proc_kasp_f,
              par_missing_dat = par_missing_dat,
              genotype_error = genotype_error,
              parent_het =  parent_het

              ))
}




#' Check and Standardize Separator Format
#'
#' This function validates and standardizes separator characters used in genotype
#' data processing. It ensures the separator is a single character and converts
#' common separators to their expected formats for different genotyping platforms.
#'
#' @param sep Character. The separator character to be checked and standardized.
#'
#' @return Character. The standardized separator:
#'   \item{" / "}{If input is "/" (Agriplex format)}
#'   \item{":"}{If input is ":" (KASP format)}
#'   \item{original}{If input is any other single character}
#'
#' @details The function:
#'   1. Removes leading and trailing whitespace from the input
#'   2. Validates that the separator is exactly one character long
#'   3. Converts "/" to " / " for Agriplex data compatibility
#'   4. Preserves ":" for KASP data compatibility
#'   5. Returns other single characters unchanged
#'
#' @examples
#' check_sep("/")     # Returns " / "
#' check_sep(":")     # Returns ":"
#'
#' # This will throw an error:
#' \dontrun{
#' check_sep("//")    # Error: separator must be single character
#' }
#'
#' @noRd
#'
#
check_sep <- function(sep) {
  no_pads <- trimws(sep)
  # Ensure the lenght of sep is one
  if (nchar(no_pads) != 1) {
    stop("Confirm your separator has a length of one character.")
  }

  # Check if it tallies with any of agriplex or Kasp
  if (no_pads == "/") {
    no_pads <- " / "
  } else if (no_pads == ":") {
    no_pads <- ":"
  }

  return(no_pads)
}




#' Generate Named List from Vector
#'
#' Creates a named list where the names are derived from the input vector.
#' Each element in the resulting list is initialized as NULL and named
#' according to the corresponding element in the input vector.
#'
#' @param x A character vector or any vector that can be coerced to character.
#'   The elements will be used as names for the resulting list. Must not be
#'   NULL, empty, or contain NA values.
#'
#' @return A named list of the same length as the input vector, with each
#'   element initialized as NULL and named according to the input vector.
#'
#' @examples
#' # Basic usage with character vector
#' generate_named_list(c("apple", "banana", "cherry"))
#' # $apple
#' # NULL
#' # $banana
#' # NULL
#' # $cherry
#' # NULL
#'
#'
#' @noRd

generate_named_list <- function(x) {
  # Input validation
  if (is.null(x)) {
    stop("Argument 'x' must not be NULL", call. = FALSE)
  }

  if (length(x) == 0) {
    stop("Argument 'x' must not be empty", call. = FALSE)
  }

  if (any(is.na(x))) {
    stop("Argument 'x' must not contain NA values", call. = FALSE)
  }

  # Convert to character
  x_char <- as.character(x)

  # Check for empty strings
  if (any(x_char == "")) {
    stop("Argument 'x' must not contain empty strings", call. = FALSE)
  }

  # Create named list more efficiently
  result <- vector(mode = "list", length = length(x_char))
  names(result) <- x_char

  return(result)
}



#' Check and install required dependencies
#' @param packages Character vector of package names to check
#' @param install_missing Logical, whether to attempt installation of missing packages
#' @return Logical vector indicating which packages are available
check_dependencies <- function(packages = NULL, install_missing = TRUE) {

  if (is.null(packages)) {
    # Define your suggested packages here
    packages <- c(
      "DT", "bslib", "data.table", "fontawesome",
      "openxlsx", "reactable", "readxl", "shinyWidgets",
      "shinyalert", "shinybusy", "shinyjs", "stringr",
      "vcfR", "writexl"
    )
  }

  # Check which packages are installed
  installed <- sapply(packages, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
  })

  missing <- packages[!installed]

  if (length(missing) > 0) {
    if (install_missing) {
      message("Missing packages detected: ", paste(missing, collapse = ", "))
      message("Attempting to install missing packages...")

      tryCatch({
        utils::install.packages(missing, dependencies = TRUE)

        # Re-check after installation
        newly_installed <- sapply(missing, function(pkg) {
          requireNamespace(pkg, quietly = TRUE)
        })

        still_missing <- missing[!newly_installed]

        if (length(still_missing) > 0) {
          stop("Failed to install packages: ", paste(still_missing, collapse = ", "))
        }

        message("All packages successfully installed!")
        return(TRUE)

      }, error = function(e) {
        stop("Error installing packages: ", e$message)
      })

    } else {
      stop("Missing required packages: ", paste(missing, collapse = ", "),
           "\nPlease install them using: install.packages(c('",
           paste(missing, collapse = "', '"), "'))")
    }
  }

  return(TRUE)
}
