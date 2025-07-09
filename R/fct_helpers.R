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
#'
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
#'
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
#'
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
  allele_cols <- paste("allele", LETTERS[seq_len(x)], sep = "_")

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
height_width <- function(string) {
  # Check if input is a character and length 1
  if (!is.character(string) || length(string) != 1) {
    stop("Input must be a single character string")
  }

  # Split the string by comma
  splitted <- strsplit(string, split = ",")[[1]]

  # Check if exactly two values were provided
  if (length(splitted) < 2) {
    stop("Input must contain exactly two values separated by a comma")
  }

  # Convert to numeric and handle potential NA values
  h_w_vector <- as.numeric(splitted)[1:2]

  if (any(is.na(h_w_vector))) {
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


#' Read Map File from CSV or Excel Format
#'
#' This function reads a map file containing marker information from either
#' CSV or Excel formats and returns it as a data frame. The function automatically
#' detects the file type based on the file extension and uses the appropriate
#' reading method.
#'
#' @param filepath Character. Full path to the map file including file name and extension.
#'   Supported formats are CSV (.csv), Excel (.xlsx), and legacy Excel (.xls).
#'
#' @return A data frame containing the map file data with marker information.
#'   All data is returned as a standard data.frame regardless of input format.
#'
#' @noRd
#'
read_mapfile <- function(filepath) {
  # Get file extension
  file_ext <- tools::file_ext(filepath)

  # Read according to file type
  if (file_ext == "csv") {
    data <- read.csv(file = filepath, stringsAsFactors = FALSE)
  } else if (file_ext %in% c("xlsx", "xls")) {
    data <- readxl::read_excel(path = filepath)
  } else {
    stop("Unsupported file type. Only CSV and Excel files are allowed.")
  }

  # Ensure it's a data.frame
  return(as.data.frame(data))
}


#' Generate Mapfile, Process KASP Data, and other user decision support info
#'
#' This function takes genotype data and processes it to generate a map file and
#' convert KASP genotype data to numeric format. It handles SNP marker parsing,
#' removes monomorphic markers, and processes parent-offspring data with various
#' quality control filters.
#'
#' @param data Data frame or tibble. The input dataset containing genotype data.
#' @param Batch Numeric. The batch number to filter and process from the dataset.
#' @param batch_col Character. Name of the column containing batch information.
#' @param marker_sep Character. Separator used in marker IDs for parsing (e.g., "_").
#' @param apply_par_poly Logical. Whether to remove polymorphic parent markers. Default is TRUE.
#' @param apply_par_miss Logical. Whether to remove markers with missing parent data. Default is TRUE.
#' @param apply_geno_good Logical. Whether to apply genotype error checking. Default is TRUE.
#' @param apply_par_homo Logical. Whether to filter out heterozygous parent markers. Default is TRUE.
#' @param genotype Character. Column name for genotype identifiers.
#' @param snp_id Character. Column name for SNP-ID markers (used when feedback = "yes").
#' @param calls_sep Character. Separator used in genotype calls (e.g., " / " for Agriplex).
#' @param data_type Character. Type of genotyping data (e.g., "agriplex", "kasp").
#' @param rp Character. Name/ID of the recurrent parent.
#' @param dp Character. Name/ID of the donor parent.
#' @param Prefix Character. Prefix for generating marker names. Default is "S".
#' @param geno_vec Character vector. Vector of genotype names/IDs corresponding to rows in the dataset.
#' @param feedback Character. One of "yes" or "no". If "yes", uses existing map file;
#'   if "no", generates new map file from marker names.
#' @param na_code Numeric. Code used for missing values. Default is NA.
#' @param data_col Character vector. Names of data columns to exclude from processing.
#' @param mapfile_path Character. Path to existing map file (required when feedback = "yes").
#'
#' @return A list containing five elements:
#'   \item{mapfile}{Data frame containing the marker map information with SNP positions}
#'   \item{proc_kasp_f}{Processed and numeric-converted KASP data ready for analysis}
#'   \item{par_missing_dat}{Data frame of loci with any parent missing genotype}
#'   \item{genotype_error}{Data frame of SNP loci with potential genotype call errors}
#'   \item{parent_het}{Data frame of heterozygous parent genotypes that were identified}
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

proc_nd_map_func <- function(data = NULL,
                        Batch = NULL, # batch number
                        batch_col = NULL, # user selected batch column
                        marker_sep = NULL, # sep for marker
                        apply_par_poly = TRUE, # response for parent poly removal
                        apply_par_miss = TRUE, # response for parent missing
                        apply_geno_good = TRUE, # response for genotype good.
                        apply_par_homo = TRUE, # filter out heterozygote parent.
                        genotype = NULL, # column names for genotypes
                        snp_id = NULL, # snp-id column for markers
                        calls_sep = NULL, # genotype calls sep-- agriplex or kasp
                        data_type = NULL, # agriplex or kasp
                        rp = NULL, # recurrent parent
                        dp = NULL, # donor parent
                        Prefix = "S", # prefix if to generate marker file
                        geno_vec = NULL, # rownames for genotype calls
                        feedback = c("yes", "no"), # whether to generate mapfile or not
                        na_code = NA, # Na codes present.
                        data_col = NULL, # data colnames
                        mapfile_path = NULL) {
  if (rp == dp) {
    stop("Error: The recurrent parent must not be the same as the donor parent")
  }
  # Validate data class. Must be a datafram or tibble
  if (!(class(data)[1] %in% c("data.frame", "tbl_df", "tbl"))) {
    stop("Please check file. Must be a dataframe or tibble.")
  }

  read_data <- data

  return_data <- read_data[read_data[batch_col] == Batch, !(colnames(read_data) %in% data_col)]


  dp_val <- which(geno_vec %in% c(rp, dp)) # get row index of rp and dp

  # Get rownames of data.
  rownames(return_data) <- geno_vec

  processed_data <- return_data # rename object

  # 1. Remove monomorphic markers if requested
  if (apply_par_poly) {
    processed_data <- parent_poly(
      x = processed_data,
      rp_row = dp_val[1],
      dp_row = dp_val[2],
      sep = calls_sep
    )
  }

  # 2. Remove missing parents if requested
  if (apply_par_miss) {
    input_data <- if (!is.null(processed_data)) processed_data else return_data

    processed_data <- parent_missing(
      x = input_data,
      rp_row = dp_val[1],
      dp_row = dp_val[2]
    )$par_present

    # Parent Missing
    par_missing_dat <- parent_missing(
      x = input_data,
      rp_row = dp_val[1],
      dp_row = dp_val[2]
    )$par_missing

  } else{
    par_missing_dat <- NULL # don't compute missing parents
  }

  # 3. Apply genotype error check if requested
  if (apply_geno_good) {
    input_data <- if (!is.null(processed_data)) processed_data else return_data

    processed_data <- geno_error(
      x = input_data,
      rp_row = dp_val[1],
      dp_row = dp_val[2],
      sep = calls_sep,
      data_type = data_type
    )$geno_good

    # Genotype error.
    genotype_error <- geno_error(
      x = input_data,
      rp_row = dp_val[1],
      dp_row = dp_val[2],
      sep = calls_sep,
      data_type = data_type
    )$geno_err

  }else{
    genotype_error <- NULL # don't compute
  }

  # 4. Remove heterozygote parents.
  if (apply_par_homo) {
    input_data <- if (!is.null(processed_data)) processed_data else return_data

    # Parent Hetero.
    processed_data <- parent_het(
      x = input_data,
      rp_row = dp_val[1],
      dp_row = dp_val[2],
      sep = calls_sep,
      na_code = na_code
    )$par_hom

    # Parent Hetero.
    parent_het <- parent_het(
      x = input_data,
      rp_row = dp_val[1],
      dp_row = dp_val[2],
      sep = calls_sep,
      na_code = NA
    )$par_het
  } else {
    parent_het <- NULL
  }

  # Final output handling
  if (is.null(processed_data)) {
    final_result <- return_data
  } else {
    final_result <- processed_data
  }

  # Parse snpids to generate mapfile.
  # if user prompts for mapfile to be generated automatic
  if (feedback == "no") {
    snps <- colnames(processed_data)

    mapfile <- parse_marker_ns(x = snps, sep = marker_sep, prefix = Prefix) |> order_markers()

    # # View(return_data)
    # parent_missing(x = return_data  ,rp_row = dp_val[1] ,dp_row =dp_val[2] )$par_missing |> View()
    # Process and convert to numeric
    proc_kasp_f <- proc_kasp(
      x = processed_data,
      sample_id = genotype,
      map_snp_id = "snpid",
      marker_start = 1,
      kasp_map = mapfile
    ) |>
      kasp_numeric(
        rp_row = dp_val[1],
        dp_row = dp_val[2],
        sep = calls_sep,
        data_type = data_type
      )
  } else if (feedback == "yes") {
    mapfile <- read_mapfile(filepath = mapfile_path) |> as.data.frame()

    proc_kasp_f <- proc_kasp(
      x = processed_data,
      sample_id = genotype, # to be selected by user
      map_snp_id = snp_id, # to be selected by user
      marker_start = 1,
      kasp_map = mapfile
    ) |>
      kasp_numeric(
        rp_row = dp_val[1],
        dp_row = dp_val[2],
        sep = calls_sep,
        data_type = data_type
      )
  }
  return(list(
    mapfile = mapfile,
    proc_kasp_f = proc_kasp_f,
    par_missing_dat = par_missing_dat,
    genotype_error = genotype_error,
    parent_het = parent_het
  ))
}

#
# gg <- proc_nd_map_func(data = geno,
#                  Batch = 3 ,
#                  batch_col = 'Batch' ,
#                  marker_sep = '_',
#                  apply_par_poly = TRUE,
#                  apply_par_miss = TRUE ,
#                  apply_geno_good = TRUE,
#                  apply_par_homo = TRUE,
#                  genotype = 'Genotype' ,
#                  snp_id = 'snpid',
#                  calls_sep = ' / ',
#                  data_type = 'agriplex',
#                  rp = "BTx623a",
#                  dp = "BTx642a" ,
#                  Prefix = 'S' ,
#                  geno_vec = genotype_names,
#                  feedback = 'no' ,
#                  na_code = NA ,
#                  data_col = data_colnames,
#                  mapfile_path = NULL
#                   )



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
#' check_sep("/") # Returns " / "
#' check_sep(":") # Returns ":"
#'
#' # This will throw an error:
#' \dontrun{
#' check_sep("//") # Error: separator must be single character
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


#' Title
#'
#' @param list_loc  a list containing the start, end , chromosome and locus name
#' @param map_doc  a mapfile of the genotype calls
#' @param genofile a genotype file containing the genotype calls
#'
#' @returns a list of the subsetted individual map file and
#'
#' @noRd
#'
combi_func <- function(list_loc,
                       map_doc,
                       genofile) {
  starts <- sapply(list_loc, function(x) x["start"]) # get all starts entered
  ends <- sapply(list_loc, function(x) x["end"]) # get all ends
  chr <- sapply(list_loc, function(x) x["chr"]) # get chromosome

  # Convert mapfile columns to lower text is they haven't
  # colnames(map_doc) <- tolower(colnames(map_doc))

  # Chromosomes must be the same
  if (length(unique(chr)) != 1) {
    # Throw an error
    stop("Chromosome of focus must be the same") # one chromosome at a time
  } else {
    chr_markers <- map_doc[map_doc["chr"] == unique(chr), ] # map file only selected chr
  }

  starts_min <- min(starts) # find the mini of the subsetted
  ends_max <- max(ends) # find max of selected

  # Must be within range.
  map_range <- range(chr_markers[["pos"]])

  user_def_reg <- c(starts_min, ends_max) # user defined start and end

  # Show range and then if user enters beyond the range, flage error.
  flag_err <- all(map_range[2] > user_def_reg)
  if (isFALSE(flag_err)) {
    stop(sprintf("Start and End positions must be within the range %d - %d", map_range[1], map_range[2]))
  }

  # subset actual postions from markers and determine the length
  actual_range <- chr_markers[["pos"]][chr_markers[["pos"]] >= starts_min & chr_markers[["pos"]] <= ends_max]

  less_min <- chr_markers[["pos"]][chr_markers[["pos"]] < starts_min]

  greater_max <- chr_markers[["pos"]][chr_markers[["pos"]] > ends_max]

  if (length(actual_range) < 20) {
    len_actual <- 20 - length(actual_range) # get excess

    bulk_up <- c(less_min, greater_max)[1:len_actual] # merge max and min

    set_pos <- sort(c(actual_range, bulk_up))
  } else if (length(actual_range) >= 20) {
    set_pos <- actual_range[seq_len(20)]
  }

  # Get marker ID for selected range
  markers_in_region <- chr_markers[chr_markers[["pos"]] %in% set_pos, ][["snpid"]]

  # Get the subset data - use the same data for both x and the chr_a calculation
  chr_a <- genofile[, colnames(genofile) %in% markers_in_region]

  chr_map <- parse_marker_ns(colnames(chr_a))

  return(list(chr_a = chr_a, chr_map = chr_map)) # return map and chromosome
}



# Get data colnames.
#' Extract colnames of Agriplex file.
#'
#' @param data  a dataframe containing genotype calls
#'
#' @returns a vector of column names from the data uploaded
#' @noRd
#'
Get_dt_coln <- function(data) {
  read_data <- data |> as.data.frame()
  # Get first 100 columns for possible extraction.
  first_hun_col <- read_data[1, 1:100]

  store_position <- which(first_hun_col %in% c("A", "G", "C", "T", NA))[1]

  data_colnames <- colnames(first_hun_col)[seq_len(store_position - 1)]

  # Return
  data_colnames
}


#' Extract colnames of Agriplex file.
#'
#' @param data  a dataframe containing genotype calls
#' @param mapfile Mapfile for genotype calls
#' @param chr chromosome number of interest.
#'
#' @returns a vector of column names from the data uploaded
#' @noRd
# Function for window size manipulation.
window_size_func <- function(data, mapfile, chr) {
  # Validate if colnames of mapfile contains snpid
  snpid <- grep(x = colnames(mapfile), pattern = "id", ignore.case = TRUE, value = TRUE)

  just_1 <- data[, mapfile[mapfile["chr"] == chr, ][[snpid]]]

  return(just_1) # return
}



#' Check and Categorize Genetic Numeric Codes
#'
#' This function checks genetic numeric codes present in the recoded data and categorizes
#' them into valid genetic markers and potential issues (monomorphic, error, or missing loci).
#'
#' @param data A dataframe containing color codes to be checked against predefined
#'   genetic marker categories.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{color_code_in}{A named numeric vector of color codes present in the input data,
#'     with names corresponding to genetic marker categories}
#'   \item{na_code}{A named numeric vector containing only the problematic loci codes
#'     (monomorphic, genotype error, and missing parental loci) that were found in the input}
#' }
#' @noRd
#'

color_code_checker <- function(data , par_1 , par_2){
  # Define color codes
  color_code_names<- c( 1, 0 , 0.5 , -1 , -2, -5 ) |> sort()
  # Respective names for color codes
  col_labels <- c('Missing', 'Error','Monomorphic',par_2,
                  'Heterozygous' , par_1
  )

  # Assign names to it.
  names(col_labels) <- as.character(color_code_names)

  result <- color_code_names[color_code_names %in% data] |> sort()

  # Get values present in data
  in_data <- result[ result  %in% color_code_names]

  needed_res <- col_labels[as.character(in_data)] # subset labels


  return(needed_res)


}



#' Validate Required Column Names
#'
#' Checks if the required column names ('Batch' and 'Genotype') are present
#' in the input data frame. The function performs case-insensitive matching.
#'
#' @param data A data frame to validate column names for
#'
#' @return NULL (invisible). The function is called for its side effect of
#'   validation. If validation passes, the function completes silently.
#'
#' @details This function requires that both 'Batch' and 'Genotype' columns
#'   are present in the input data frame. Column name matching is case-insensitive,
#'   so 'batch', 'BATCH', 'Batch' are all acceptable.
#' @noRd
#'
#
check_colnames_validate <- function(data){
  must_hv_names <- c('Batch', 'Genotype') |> tolower() # required columns

  if(!(all(must_hv_names %in% tolower(colnames(data))))){

    stop(paste("Missing required columns: 'Batch', 'Genotype'"))
  }
}




