#' Helper Functions for PanGenomeBreedr Application
#'
#' This file contains various helper functions used throughout the application.
#' Functions are organized by their purpose and usage context.
#' None of these functions are exported as they're meant for internal use only.
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
    formatted_ids <- as.vector(regmatches(SNP_ids[variable], regexec("([a-zA-Z]+)(\\d+)", SNP_ids[variable]))[[1]])

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

  # Unlist if needed
  x <- unlist(x)

  if (length(x) == 0) {
    return(data.frame( Genotypes = 'None Found / Uncallable'))
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
  x <- unlist(x) # unlist object

  if (length(x) == 0) {
    # Return an empty dataframe with Non found
    return(data.frame(Alleles = 'None Found / Uncallable'))
  }

  # Create allele column names
  allele_cols <- paste("allele", LETTERS[seq_len(length(x))], sep = "_")

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
Genotypes_user <- function(data, sample_id = 'Genotype') {
  # Select column for genotypes
  just_geno <- data[[sample_id]]

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
    data <- readxl::read_excel(path = filepath,.name_repair = 'minimal')
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
#' @param marker_start Numeric. The index/column number where marker data begins.
#' @param sample_id Character. Column name containing genotype identifiers.
#' @param marker_sep Character. Separator used in marker IDs for parsing (e.g., "_").
#' @param apply_par_poly Logical. Whether to remove monomorphic parent markers. Default is TRUE.
#' @param apply_par_miss Logical. Whether to remove markers with missing parent data. Default is TRUE.
#' @param apply_geno_good Logical. Whether to apply genotype error checking. Default is TRUE.
#' @param apply_par_homo Logical. Whether to filter out heterozygous parent markers. Default is TRUE.
#' @param snp_id Character. Column name for SNP-ID markers (used when feedback = "yes").
#' @param calls_sep Character. Separator used in genotype calls (e.g., " / " or "/").
#' @param data_type Character. Type of genotyping data (e.g., "agriplex", "kasp").
#' @param rp Character. Name/ID of the recurrent parent.
#' @param dp Character. Name/ID of the donor parent.
#' @param Prefix Character. Prefix for generating marker names. Default is "S".
#' @param feedback Character. One of "yes" or "no". If "no", generates new map file from marker names;
#'   if "yes", uses existing map file. Default is "yes".
#' @param na_code Numeric or Character. Code used for missing values. Default is NA.
#' @param mapfile_path Character. Path to existing map file (required when feedback = "yes").
#'
#' @return A list containing five elements:
#' \item{mapfile}{Data frame containing the marker map information with SNP positions}
#' \item{proc_kasp_f}{Processed and numeric-converted KASP data ready for analysis}
#' \item{par_missing_dat}{Data frame of loci with any parent missing genotype}
#' \item{genotype_error}{Data frame of SNP loci with potential genotype call errors}
#' \item{parent_het}{Data frame of heterozygous parent genotypes that were identified}
#'
#' @details The function performs the following steps:
#' 1. Identifies genotype data columns starting from \code{geno_start}
#' 2. Filters and renames rows using \code{geno_vec}
#' 3. Sequentially applies QC filters (Polymorphism, Missingness, Error Checking, and Zygosity)
#' 4. Either generates a new map file based on marker name parsing or loads an existing one
#' 5. Converts genotype calls to numeric format for downstream analysis
#'
#' @noRd
proc_nd_map_func <- function(data = NULL,
                             marker_start = NULL,
                             sample_id = NULL,
                             marker_sep = NULL,
                             apply_par_poly = TRUE,
                             apply_par_miss = TRUE,
                             apply_geno_good = TRUE,
                             apply_par_homo = TRUE,
                             snp_id = NULL,
                             calls_sep = NULL,
                             data_type = NULL,
                             rp = NULL,
                             dp = NULL,
                             Prefix = "S",
                             feedback = "yes",
                             na_code = NA,
                             mapfile_path = NULL){
  if (rp == dp) {
    stop("Error: The recurrent parent must not be the same as the donor parent")
  }

  # Coerce data and set identifiers
  return_data <- as.data.frame(data)
 # rownames(return_data) <- return_data[[sample_id]]

  meta_data <- return_data[,colnames(return_data)[1: (marker_start-1)]] # Extract metadata
  processed_data <- return_data[,colnames(return_data)[marker_start: (ncol(return_data))]]


  # Remove monomorphic markers
  if (apply_par_poly) {
    processed_data <- parent_poly(
      x = processed_data,
      rp_row = rp,
      dp_row = dp,
      sep = calls_sep
    )
  }

  # Handle Missing Parents
  if (apply_par_miss) {

    miss_results <- parent_missing(
      x = processed_data,
      rp_row = rp,
      dp_row = dp
    )
    processed_data <- miss_results$par_present
    par_missing_dat <- miss_results$par_missing
  } else {
    par_missing_dat <- NULL
  }

  # Apply Genotype Error Check
  if (apply_geno_good) {

    err_results <- geno_error(
      x = processed_data,
      rp_row = rp,
      dp_row = dp,
      sep = calls_sep,
      data_type = data_type
    )
    processed_data <- err_results$geno_good
    genotype_error <- err_results$geno_err
  } else {
    genotype_error <- NULL
  }

  # Remove Heterozygote Parents
  if (apply_par_homo) {

    het_results <- parent_het(
      x = processed_data,
      rp_row = rp,
      dp_row = dp,
      sep = calls_sep,
      na_code = na_code
    )
    processed_data <- het_results$par_hom
    parent_het <- het_results$par_het
  } else {
    parent_het <- NULL
  }


  # Final Output Selection
  # If any filter was active, return the refined data else return original
  if (any(apply_par_poly, apply_par_miss, apply_geno_good, apply_par_homo)) {
    # append meta data to processed_data
    final_result <- cbind(meta_data, processed_data)

  } else {
    final_result <- return_data
  }



  # If user prompts for mapfile to be generated automatically
  if (feedback == "no") {
    # Get snp id
    snps <- colnames(final_result[,marker_start :ncol(final_result)])

    mapfile <- tryCatch({
      parse_marker_ns(x = snps, sep = marker_sep, prefix = Prefix) |>
        order_markers()
    }, error = function(e) {
      stop(paste("Failed to generate mapfile automatically:", e$message))
    })

    # Column indexing
    snpid <- safe_grep_match("snpid",choices = colnames(mapfile))
    chr <- safe_grep_match('chr',choices = colnames(mapfile))
    pos <- safe_grep_match('pos',choices = colnames(mapfile))

    # Filter the mapfile to contain only snp_ids in the processed data.
    mapfile <- mapfile[mapfile[[snpid]] %in% colnames(processed_data), ]

    # Process and numeric code the data
    proc_kasp_f <- proc_kasp(
      x = final_result,
      sample_id = sample_id,
      map_snp_id = snpid ,
      marker_start = marker_start,
      chr = chr,
      chr_pos = pos ,
      kasp_map = mapfile
    ) |>
      kasp_numeric(
        rp_row = rp,
        dp_row = dp,
        sep = calls_sep,
        data_type = data_type
      )

  } else if (feedback == "yes") {
    # User provides mapfile manually
    mapfile <- mapfile_path
    # get colnumn for mapfile
    snp_colname <- safe_grep_match(pattern = 'snp', colnames(mapfile))

    # Filter out unneeded  snip_ids
    # mapfile <- mapfile[mapfile[[snp_colname]] %in% colnames(processed_data),]
    # Keep only SNPs present in BOTH data and mapfile, in the same order
    common_snps <- intersect(mapfile[[snp_colname]], colnames(final_result))

    # subset mapfile to common_snps and preserve order of processed_data
    mapfile <- mapfile[match(common_snps, mapfile[[snp_colname]]), ]

    # Column indexing
    snpid <- safe_grep_match("snpid",choices = colnames(mapfile))
    chr <- safe_grep_match('chr',choices = colnames(mapfile))
    pos <- safe_grep_match('pos',choices = colnames(mapfile))

    # Process and numeric code the data
    proc_kasp_f <- proc_kasp(
      x = final_result,
      sample_id = sample_id,
      map_snp_id = snpid,
      marker_start = marker_start,
      chr = chr,
      chr_pos = pos,
      kasp_map = mapfile
    ) |>
      kasp_numeric(
        rp_row = rp,
        dp_row = dp,
        sep = calls_sep,
        data_type = data_type
      )
  } else {
    stop("Invalid 'feedback' value. Must be either 'yes' or 'no'.")
  }

  return(list(
    mapfile = mapfile,
    proc_kasp_f = proc_kasp_f,
    rp_name = rownames(proc_kasp_f)[rp], # recurrent parent name based on index
    dp_name =  rownames(proc_kasp_f)[dp] # donor parent name based on index
  ))
}


#' Convert SNP Positions to Genomic Windows
#'
#' @param pos_list A named list where each element is a numeric vector with 'chr' and 'pos'.
#' @return A list with 'chr', 'start', and 'end' for each trait.
#' @noRd
#'
get_trait_windows <- function(pos_list) {

  # Use lapply to iterate through each trait in your list
  lapply(pos_list, function(x) {

    # Extract original values
    chr <- unname(x["chr"])
    pos <- unname(x["pos"])

    # Add 50kb upstream and 500kb downstream of the defined position
    start <- pos - 1e6
    end   <- pos + 1.5e6

    # Return the new structure
    c(chr = chr, start = start, end = end)
  })
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
  color_code_names <- c( 1, 0 , 0.5 , -1 , -2, -5 ) |> sort()
  # Respective names for color codes
  col_labels <- c('Missing', 'Error','Monomorphic',par_2,
                  'Heterozygous' , par_1 )

  # Assign names to it.
  names(col_labels) <- as.character(color_code_names)

  result <- color_code_names[color_code_names %in% data] |> sort()

  # Get values present in data
  in_data <- result[ result  %in% color_code_names]

  needed_res <- col_labels[as.character(in_data)] # subset labels


  return(needed_res)


}








#' Read VCF file and return as data frame
#'
#' This function reads a Variant Call Format (VCF) file and converts it into
#' a data frame with appropriate column names extracted from the VCF header.
#'
#' @param vcf_file Character string specifying the path to the VCF file to read
#'
#' @return A data frame containing the VCF data with columns named according to
#'   the VCF header (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, and
#'   any sample columns)
#'
#'
#' @noRd
#'
read_vcf_as_df <- function(vcf_file) {
  # Read VCF header to get column names
  header_lines <- readLines(vcf_file)
  header <- header_lines[grep("^#CHROM", header_lines)]
  colnames <- strsplit(header, "\t")[[1]]
  colnames[1] <- "CHROM" # fix formatting
  # Read VCF data (skip header lines)
  vcf_df <- utils::read.table(vcf_file,
                              comment.char = "#", header = FALSE, sep = "\t",
                              col.names = colnames, stringsAsFactors = FALSE
  )
  return(vcf_df)
}




#' Safe Pattern Matching in Character Vector
#'
#' Safely searches for a pattern in a character vector and returns the first match.
#' If no match is found, returns the first element of the vector. If the vector
#' is empty, returns NULL.
#'
#' @param pattern A character string containing a regular expression pattern to search for.
#' @param choices A character vector to search within. Can be empty.
#' @noRd
#'
safe_grep_match <- function(pattern, choices) {
  if (length(choices) == 0) return(NULL)

  matches <- grep(pattern, choices, ignore.case = TRUE,value = TRUE)
  if (length(matches) > 0) {
    return(matches)
  }else{
    return(choices[[1]])
  }
}


#' Generate an UpSet Plot for Foreground Marker Selection
#'
#' @description
#' This function generates an UpSet plot using the \pkg{UpSetR} package to visualize
#' marker intersections from a foreground selection matrix, annotated with metadata
#' extracted from a corresponding marker information table.
#'
#' @param foreground_matrix A binary matrix or data frame of foreground selection data
#'        with markers as columns and genotypes (or individuals) as rows.
#' @param marker_info A data frame containing marker metadata.
#'        Must include columns for marker names and loci.
#' @param mainbar_y_label A character string specifying the y-axis label for the main intersection bar plot.
#' @param sets_x_label A character string specifying the x-axis label for the set size bars.
#' @param text_scale A numeric value controlling the scaling of plot text labels.
#' @param plot_type The type of metadata plot to display in the set metadata panel.
#'        Defaults to `"text"`.
#' @param assign An integer specifying the position along the axis where the metadata plot should be displayed.
#' @param column The name of the metadata column to use for annotating the set metadata panel.
#' @param colors A character string or vector of color names to apply to the metadata plot.
#'
#'@noRd
#'
run_upset_plot <- function(foreground_matrix,
                           marker_info,
                           mainbar_y_label = "Locus Intersection Size",
                           sets_x_label = "Locus Size",
                           text_scale = 1.2,
                           plot_type = "text",
                           assign = 8,
                           column = "locus",
                           colors = "firebrick2") {

  # Create metadata
  metadata <- data.frame(
    sets = marker_info[[safe_grep_match(pattern = 'mark',choices = colnames(marker_info) )]],
    locus = marker_info[[safe_grep_match(pattern = 'locus',choices = colnames(marker_info))]]
  )

  # Number of markers / sets
  nl <- ncol(foreground_matrix)

  # Prepare colors (recycle as needed)
  color_vec <- rep(colors, length.out = nl)

  # Generate UpSet plot with metadata
  plot <-  suppressWarnings(
    UpSetR::upset(
      foreground_matrix,
      nsets = nl,
      mainbar.y.label = mainbar_y_label,
      sets.x.label = sets_x_label,
      text.scale = text_scale,
      set.metadata = list(
        data = metadata,
        plots = list(
          list(
            type = plot_type,
            assign = assign,
            column = column,
            colors = color_vec
          )
        )
      )
    )
  )
  return(plot)
}





#' Compute summary result from par missing function.
#'
#' @param df a dataframe from parent missing for summary
#'
#' @returns a data frame of summarised snps
#' @noRd
par_missing_dat <- function(df){
  #Count of missing genotype calls per SNP locus (column)
  missing_per_locus <- colSums(is.na(df))

  # Count of missing genotype calls per sample (row)
  missing_per_sample <- rowSums(is.na(df))

  # Print missing calls per sample separately for reference
  missing_samples_df <- data.frame(
    Sample = rownames(df),
    Missing_Calls = missing_per_sample,
    Total_SNPs = colnames(df) |> length() ,
    row.names = NULL
  )

  return(missing_samples_df)
}




#' A padded dataframe for organizing found lines.
#'
#' This utility function takes vectors of varying lengths (SNPs present, SNPs absent, and sample names)
#' and combines them into a single data frame. It automatically pads shorter vectors
#' with NA values to ensure all columns have equal length based on the longest input.
#'
#' @param snps_present A character vector of SNP IDs that must be present.
#' @param snps_absent A character vector of SNP IDs that must be absent.
#' @param sample_name A character vector of genotype or sample names identified.
#'
#' @return A data frame with three columns: snps_present, snps_absent,
#' and sample_name, padded with NA where necessary.
#'
#' @noRd
#'
create_padded_df <- function(snps_present, snps_absent, sample_name) {

  # Determine the maximum length among all inputs
  actual_max_len <- max(
    length(snps_present),
    length(snps_absent),
    length(sample_name),
    na.rm = TRUE
  )

  # Helper to pad vectors
  pad_vector <- function(x, len) {
    if (length(x) < len) {
      return(c(x, rep(NA, len - length(x))))
    }
    return(x)
  }

  # Create and return data frame
  data.frame(
    snps_present = pad_vector(snps_present, actual_max_len),
    snps_absent  = pad_vector(snps_absent, actual_max_len),
    sample_name  = pad_vector(sample_name, actual_max_len),
    stringsAsFactors = FALSE
  )
}




#' Generate a Trait Input Row for Genomic Coordinates
#'
#' This helper function creates a stylized UI row containing inputs for
#' a trait name, chromosome number, and physical position. It is used
#' within the coordinate definition wizard to allow dynamic addition of multiple loci.
#'
#' @param id A unique integer or character string to suffix the input IDs,
#' ensuring namespace consistency.
#'
#' @return A \code{shiny::tag.list} (div) containing a \code{fluidRow} with
#' \code{textInput} and \code{numericInput} elements.
#'
#' @note The 'Remove' button is programmatically excluded for the first input
#' row (\code{id = 1}) to ensure at least one locus remains defined.
#'
#' @importFrom shiny div fluidRow column textInput numericInput actionButton
#' @noRd
makeTraitInput <- function(id, ns) {
  div(
    id = ns(paste0("trait_", id)),
    style = "border: 1px solid #dee2e6; padding: 15px; margin-bottom: 15px; border-radius: 8px; background-color: #f8f9fa; position: relative;",

    fluidRow(
      # Basic Info
      column(3, textInput(ns(paste0("name_", id)), "Trait/QTL Name:", value = paste0("QTL_", id))),
      column(2, numericInput(ns(paste0("chr_", id)), "Chr:", value = 1, min = 1)),

      # Configuration Checkboxes
      column(7,
             checkboxGroupInput(
               inputId = ns(paste0("type_", id)),
               label = "Annotation Components (Select at least one):",
               choices = c("Point Position" = "pos", "Genomic Range" = "range"),
               selected = "range",
               inline = TRUE
             )
      )
    ),

    fluidRow(
      # Conditional Point Input
      conditionalPanel(
        condition = paste0("input['", ns(paste0("type_", id)), "'].includes('pos')"),
        column(4, numericInput(ns(paste0("pos_", id)), "Exact Position (bp):", value = 1200000))
      ),

      # Conditional Range Inputs
      conditionalPanel(
        condition = paste0("input['", ns(paste0("type_", id)), "'].includes('range')"),
        column(4, numericInput(ns(paste0("start_", id)), "Start (bp):", value = 1000000)),
        column(4, numericInput(ns(paste0("end_", id)), "End (bp):", value = 1400000))
      )
    ),

    # Remove button: only for traits after the first one
    if (id > 1) {
      div(style = "position: absolute; top: 10px; right: 10px;",
          actionButton(ns(paste0("remove_", id)), "", icon = icon("trash"),
                       class = "btn-outline-danger btn-sm", title = "Remove Trait"))
    }
  )
}



#' Render a Ghost Plate Error Plot
#'
#' @description This function generates a stylized ggplot2 plot that visually
#' represents a plate layout error. It displays a "ghost" plate background 
#' with an overlayed error message.
#'
#' @param error_msg A character string containing the specific error message to display.
#' @param text_size A numeric value controlling the size of the text elements in the plot.
#' @return A ggplot2 object representing the error plot.
#' @noRd
render_ghost_plate_error <- function(error_msg, text_size = 12) {
  
  # Data setup
  ghost_data <- expand.grid(
    row = factor(rev(LETTERS[1:8])), 
    col = factor(sprintf("%02d", 1:12))
  )

  # Color Palette
  well_grey   <- "#F1F2F6"
  stroke_grey <- "#CED4DA"
  grid_grey   <- "#E9ECEF"
  text_dark   <- "#2D3436"
  error_red   <- "#D63031"

  ggplot2::ggplot(ghost_data, ggplot2::aes(x = col, y = row)) +
    # Draw the empty grid
    ggplot2::geom_tile(fill = "white", color = grid_grey, linewidth = 0.5) +
    
    # Draw the empty wells
    ggplot2::geom_point(
      size = text_size,
      shape = 21,
      fill = well_grey,
      color = stroke_grey,
      stroke = 1,
      alpha = 0.8
    ) +
    
    # Semi-transparent overlay to "de-emphasize" the background
    ggplot2::annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf,
      fill = "white",
      alpha = 0.3
    ) +
    
    # Labels and Theme setup
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(title = "PLATE LAYOUT ERROR") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = text_size,
        color = error_red
      ),
      axis.text = ggplot2::element_text(
        size = text_size,
        face = "bold",
        color = "#B2BEC3"
      )
    ) +
    
    # Central Error Message Label
    ggplot2::annotate(
      "label",
      x = 6.5,
      y = 4.5,
      label = paste0("VALIDATION FAILED:\n", error_msg),
      fill = "white",
      color = text_dark,
      size = text_size * 0.35,
      fontface = "bold",
      label.padding = ggplot2::unit(1, "lines")
    )
}





#' Render a Reactable Table with Custom Theming
#'
#' @description This helper function generates a `reactable` table with a default
#'   theme and common display options. It allows for easy and consistent styling
#'   of interactive tables throughout the application.
#'
#' @param data A data frame or tibble to be displayed in the reactable table.
#' @param theme An optional `reactableTheme` object to override the default theme.
#'   If `NULL`, a predefined default theme is used.
#'
#' @return A `reactable` HTML widget displaying the provided data.
#'
#' @noRd
 # Reactable Helper
    render_reactable <- function(data, theme = NULL) {
      if (is.null(theme)) {
        theme <- reactable::reactableTheme(
          backgroundColor = "hsl(0, 0%, 100%)",
          borderColor = "hsl(0, 0%, 89%)",
          stripedColor = "hsl(0, 0%, 97%)",
          highlightColor = "hsl(0, 0%, 96%)",
          cellPadding = "8px 12px",
          searchInputStyle = list(width = "100%", padding = "8px 12px"),
          headerStyle = list(borderWidth = "1px", padding = "8px 12px")
        )
      }
      reactable::reactable(
        data,
        defaultColDef = reactable::colDef(
          headerClass = "header",
          align = "left",
          minWidth = 100,
          width = 200
        ),
        pagination = TRUE,
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(10, 25, 50, 100),
        striped = TRUE,
        highlight = TRUE,
        bordered = TRUE,
        compact = TRUE,
        wrap = TRUE,
        theme = theme,
        style = list(maxWidth = "100%")
      )
    }






#' Plot Variant Hotspots
#'
#' @param variant_table Data frame containing basic variants
#' @param annotation_table Data frame containing annotations
#' @param region_start Numeric start position
#' @param region_end Numeric end position
#' @import ggplot2
#' @noRd
#'
plot_variant_hotspots <- function(
  variant_table,
  annotation_table = NULL,
  region_start = NULL,
  region_end = NULL
) {
  pos <- y_pos <- impact <- variant_type <- NULL
  
  plot_df <- variant_table
  if (!"variant_type" %in% colnames(plot_df)) {
    plot_df$variant_type <- sub("_.*", "", plot_df$variant_id)
  }

  # Join and Deduplicate Annotations
  impact_levels <- c("HIGH", "MODERATE", "LOW", "MODIFIER", "UNKNOWN")

  if (
    !is.null(annotation_table) &&
      nrow(annotation_table) > 0 &&
      "impact" %in% colnames(annotation_table)
  ) {
    anno_sub <- annotation_table[, c("variant_id", "impact")]
    anno_sub$impact <- factor(
      anno_sub$impact,
      levels = impact_levels,
      ordered = TRUE
    )
    anno_sub <- anno_sub[order(anno_sub$impact, na.last = TRUE), ]
    anno_sub <- anno_sub[!duplicated(anno_sub$variant_id), ]

    plot_df <- merge(plot_df, anno_sub, by = "variant_id", all.x = TRUE)
    plot_df$impact <- as.character(plot_df$impact)
    plot_df$impact[is.na(plot_df$impact)] <- "UNKNOWN"
  } else {
    plot_df$impact <- "UNKNOWN"
  }

  plot_df$impact <- factor(plot_df$impact, levels = rev(impact_levels))

  # Extract start and stop postion from table if is null
  if (is.null(region_start)) {
    region_start <- min(plot_df$pos, na.rm = TRUE)
  }
  if (is.null(region_end)) {
    region_end <- max(plot_df$pos, na.rm = TRUE)
  }

  total_vars <- nrow(plot_df)
  snp_ct <- sum(plot_df$variant_type == "SNP", na.rm = TRUE)
  indel_ct <- sum(plot_df$variant_type == "INDEL", na.rm = TRUE)

  summary_text <- sprintf(
    "Total Unique Variants: %s  |  SNPs: %s  |  INDELs: %s",
    scales::comma(total_vars),
    scales::comma(snp_ct),
    scales::comma(indel_ct)
  )

  # Assigning base elevations to each impact level
  y_mapping <- c(
    "HIGH" = 4,
    "MODERATE" = 3,
    "LOW" = 2,
    "MODIFIER" = 1,
    "UNKNOWN" = 0.5
  )
  plot_df$base_y <- y_mapping[as.character(plot_df$impact)]

  # Add micro-jitter within lanes so points don't perfectly overlap
  set.seed(42)
  plot_df$y_pos <- plot_df$base_y +
    runif(nrow(plot_df), min = -0.25, max = 0.25)

  # Color Impact levels
  impact_colors <- c(
    "HIGH" = "#e11d48",
    "MODERATE" = "#f59e0b",
    "LOW" = "#3b82f6",
    "MODIFIER" = "#94a3b8",
    "UNKNOWN" = "#cbd5e1"
  )

  type_shapes <- c("SNP" = 21, "INDEL" = 24)

  # Plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pos, y = y_pos)) +
    ggplot2::geom_hline(yintercept = c(1, 2, 3, 4), color = "#f1f5f9", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 2, color = "#cbd5e1") +
    ggplot2::geom_segment(
      ggplot2::aes(xend = pos, yend = 0),
      color = "#dcdfe3",
      alpha = 0.3,
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = impact, shape = variant_type),
      size = 3.5,
      color = "black",
      stroke = 0.5,
      alpha = 0.5
    ) +
    # Manual setting of legends and mapping values
    ggplot2::scale_fill_manual(values = impact_colors, name = "Max Impact") +
    ggplot2::scale_shape_manual(values = type_shapes, name = "Variant Type") +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::coord_cartesian(xlim = c(region_start, region_end), ylim = c(-0.1, 4.5)) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "Genomic Variant Hotspots",
      subtitle = paste0(
        plot_df$chrom[1],
        ": ",
        scales::comma(region_start),
        " - ",
        scales::comma(region_end)
      ),
      caption = summary_text,
      x = "Genomic Position (bp)"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4), order = 1),
      shape = ggplot2::guide_legend(
        override.aes = list(fill = "#fefefe14", size = 4, color = 'black'),
        order = 2
      ) # Modifyying legends
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#e2e8f0", linetype = "dashed"),
      
      legend.position = "right",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(face = "bold", size = 10, color = "#475569"),
      legend.text = ggplot2::element_text(size = 9, color = "#334155"),
      
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 16,
        color = "#0f172a",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        color = "#64748b",
        size = 12,
        hjust = 0.5,
        margin = ggplot2::margin(b = 10)
      ),
      plot.caption = ggplot2::element_text(
        face = "bold",
        color = "#1e293b",
        size = 12,
        hjust = 1,
        margin = ggplot2::margin(t = 15)
      ),
      axis.title.x = ggplot2::element_text(
        face = "bold",
        color = "#475569",
        margin = ggplot2::margin(t = 10)
      )
    )

  return(p)
}




#' Update Sidebar Navigation Buttons
#'
#' This function dynamically updates the CSS classes of sidebar navigation buttons
#' to visually indicate the active tab. It uses `shinyjs` to toggle Bootstrap 5
#' button classes between outlined and filled states.
#'
#' @param active_id Character string containing the input ID of the currently active button.
#'
#' @return No return value. Called for its side effects on the UI.
#'
#' @noRd
update_sidebar_buttons <- function(active_id) {
  buttons <- c(
    "get_db_info",
    "show_gene_cord_btn",
    "show_query_actions_btn",
    "get_pcv_sidebar_btn",
    "design_kasp_sidebar_btn"
  )
  for (btn_id in buttons) {
    if (btn_id == active_id) {
      shinyjs::removeClass(btn_id, "btn-outline-primary")
      shinyjs::addClass(btn_id, "btn-primary")
    } else {
      shinyjs::removeClass(btn_id, "btn-primary")
      shinyjs::addClass(btn_id, "btn-outline-primary")
    }
  }
}