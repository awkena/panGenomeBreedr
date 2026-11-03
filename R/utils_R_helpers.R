# Internal R helper functions for the panGenomeBreedr application
#
# Plain R helpers (data reading, parsing, QC summaries, plotting) used by the
# Shiny modules. None of them depend on a Shiny session; helpers that build UI
# or drive the app live in utils_shiny.R. None of these functions are exported.

#' Extract marker IDs and chromosome IDs from VCF file
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


#' Read map file from CSV or Excel format
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


#' Generate mapfile, process KASP data, and other user decision support info
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


#' Check and standardize separator format
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


#' Check and categorize genetic numeric codes
#'
#' This function checks genetic numeric codes present in the recoded data and categorizes
#' them into valid genetic markers and potential issues (monomorphic, error, or missing loci).
#'
#' @param data A dataframe containing color codes to be checked against predefined
#'   genetic marker categories.
#' @param par_1 parent 1
#' @param par_2 parent 2
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


#' Safe pattern matching in character vector
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


#' Generate an UpSet plot for foreground marker selection
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

