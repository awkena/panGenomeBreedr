#' Set the API endpoint URL for your private database 
#'
#' @param url Character. The full URL to your hosted API endpoint
#'   (e.g., "http://rice-genomics.myuniversity.edu:8000").
#'
#' @export
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Set a custom API endpoint URL
#' set_api_url("http://my-custom-api.com:8000")
#'
#' }
set_api_url <- function(url) {
  # Strip any accidental trailing slashes
  url <- sub("/$", "", url)

  # Set R options in current session
  options(panGenomeBreedr.api_url = url)

  Sys.setenv(PANGENOME_API_URL = url)

  message("panGenomeBreedr API endpoint successfully set to: ", url)
}





#' Get the current database API endpoint URL
#'
#' Retrieves the active API endpoint. If no custom endpoint has been set
#' by the user, it safely defaults to the official Sorghum database API.
#'
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Get the current API endpoint URL
#' get_api_url()
#'
#' }
#'
#' @returns The character string of the active API URL.
#' @export
get_api_url <- function() {
  # check if the user set a temporary custom URL in this session
  url <- getOption("panGenomeBreedr.api_url")

  # check if they set a permanent custom URL in their .Renviron
  if (is.null(url) || url == "") {
    url <- Sys.getenv("PANGENOME_API_URL")
  }

  # If nothing is set, use the sorghum AWS Server
  if (is.null(url) || url == "") {
    url <- "http://79.72.72.212:8000"
  }

  return(url)
}




#' Internal API fetcher
#' @noRd
#' @importFrom jsonlite fromJSON
.api_fetch <- function(endpoint, query = NULL, simplify = TRUE) {
  # Get url , user defined or default
  base_url <- get_api_url()
  url <- paste0(base_url, endpoint)

  # Send the request to the server
  response <- httr::GET(url, query = query)

  # Check for server errors
  if (httr::http_error(response)) {
    stop(
      paste0(
        "API Error [",
        httr::status_code(response),
        "]: Check parameters or server status."
      ),
      call. = FALSE
    )
  }

  # Extract and parse the JSON payload
  raw_json <- httr::content(response, as = "text", encoding = "UTF-8")
  data <- jsonlite::fromJSON(raw_json, simplifyDataFrame = TRUE)

  # Only force into a single dataframe if we are expecting one table
  if (simplify && is.list(data) && !is.data.frame(data)) {
    return(as.data.frame(data))
  }

  return(data)
}




#' List all tables in the database (online)
#'
#' @examples
#' \donttest{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # List all tables in the online database
#' pg_list_tables_result <- pg_list_tables()
#' print(pg_list_tables_result)
#' }
#'
#' @returns A character vector of table names.
#' @export
#' @importFrom httr GET content http_error
pg_list_tables <- function() {
  return(.api_fetch("/db/tables"))
}



#' Get variant statistics from the database (online)
#' 
#' @param include_annotations A logical value indicating whether to include
#' statistics for the annotations table. Defaults to \code{TRUE}.
#' 
#' @returns A data frame containing variant statistics.
#' 
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Get variant statistics from the online database, including annotation counts
#' pg_variant_stats_result <- pg_variant_stats(include_annotations = TRUE)
#' print(pg_variant_stats_result)
#'
#' # Get variant statistics without annotation counts
#' pg_variant_stats_no_ann <- pg_variant_stats(include_annotations = FALSE)
#' print(pg_variant_stats_no_ann)
#' }
#' @importFrom httr GET content http_error
#' @export
pg_variant_stats <- function(include_annotations = TRUE) {
  # Pass the logical parameter as a named list
  params <- list(include_annotations = include_annotations)

  return(.api_fetch("/stats/variants", query = params))
}




#' Get summary of database variant impacts by chromosome (online)
#'
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Get variant impact summary from the online database
#' pg_variant_impact_summary_result <- pg_variant_impact_summary()
#' print(pg_variant_impact_summary_result)
#' }
#'
#' @returns A data frame in wide format where each row is a chromosome and
#' columns represent the counts for each impact category.
#' 
#' @importFrom httr GET content http_error
#' @export

pg_variant_impact_summary <- function() {
  return(.api_fetch("/stats/impact"))
}




#' Get table names and row counts from the database (online)
#' 
#' @returns A data frame with two columns: \code{table} and \code{n_rows}.
#' 
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Summarize all tables in the online database
#' pg_summarize_tables_result <- pg_summarize_tables()
#' print(pg_summarize_tables_result)
#' }
#'
#' @export
#' @importFrom httr GET content http_error
pg_summarize_tables <- function() {
  return(.api_fetch("/db/summary"))
}




#' Get table schema details (online)
#'
#' @param table_name A character value specifying the name of the table.
#' Defaults to "variants". Available options are "variants", "annotations",
#' "genotypes", or "metadata".
#'
#' @returns A data frame containing column metadata.
#' 
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Get column metadata for the "variants" table
#' pg_list_table_columns_result <- pg_list_table_columns(table_name = "variants")
#' print(pg_list_table_columns_result)
#' }
#'
#' @export
#' @importFrom httr GET content http_error
pg_list_table_columns <- function(
  table_name = c("variants", "annotations", "genotypes", "metadata")
) {
  # Validate choice locally before hitting the API
  table_name <- match.arg(table_name)

  # Pass parameter as a list; httr handles the URL encoding (the '?' and '=')
  params <- list(table_name = table_name)

  return(.api_fetch("/db/columns", query = params))
}




#' Query database tables by genomic coordinates (online)
#'
#' This function connects to the public panGenomeBreedr API to retrieve data
#' from the variants, annotations, or genotypes tables based on a specific
#' chromosome and genomic range.
#'
#' @param table_name A character value specifying the target view to query. 
#'   Must be one of \code{"variants"}, \code{"annotations"}, or \code{"genotypes"}.
#' @param chrom A character value specifying the target chromosome name (e.g., \code{"Chr05"}).
#' @param start Integer. Start coordinate for the target window region.
#' @param end Integer. End coordinate for the target window region.
#' @param gene_name A character value indicating the specific Sobic gene model ID. 
#'   Utilized explicitly when subsetting the \code{"annotations"} layout.
#'
#' @returns A data frame containing the queried genomic data.
#' 
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Query variants table for a specific genomic region
#' variants_data <- pg_query_db(
#'   table_name = "variants",
#'   chrom = "Chr05",
#'   start = 75104537,
#'   end = 75106403
#' )
#' print(head(variants_data))
#'
#' # Query annotations table for a specific gene
#' annotations_data <- pg_query_db(
#'   table_name = "annotations",
#'   chrom = "Chr05",
#'   start = 75104537,
#'   end = 75106403,
#'   gene_name = "Sobic.005G213600"
#' )
#' print(head(annotations_data))
#'
#' # Query genotypes table for the same region
#' genotypes_data <- pg_query_db(table_name = "genotypes", chrom = "Chr05",
#'                               start = 75104537, end = 75106403)
#' print(head(genotypes_data))
#' }
#' 
#' @importFrom httr GET content http_error
#' @export
pg_query_db <- function(
  table_name = c("variants", "annotations", "genotypes"),
  chrom = NULL,
  start = NULL,
  end = NULL,
  gene_name = NULL
) {
  table_name <- match.arg(table_name)

  # Pack the arguments. NULLs are handled automatically by .api_fetch.
  query_params <- list(
    table_name = table_name,
    chrom = chrom,
    start = start,
    end = end,
    gene_name = gene_name
  )

  # Query result
  result <- .api_fetch("/db/query", query = query_params)

  if (table_name == 'genotypes') {
    # Reformat to include allele types and fequencies.
    query_result <- append_locus_allele_metrics(region_matrix = result)
    return(query_result)
  } else {
    return(result)
  }

}




#' Extract variants based on functional mutation impact (online)
#'
#' @param impact_level A character vector specifying the variant impact types to retain. 
#'   Allowed classifications are \code{"HIGH"}, \code{"MODERATE"}, \code{"LOW"}, and \code{"MODIFIER"}.
#' @param chrom A character value specifying the chromosome name (e.g., \code{"Chr05"}).
#' @param start Integer. Start coordinate for the genomic region.
#' @param end Integer. End coordinate for the genomic region.
#'
#' @returns A data frame containing variant information and associated functional
#'   annotations.
#' 
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Query for high-impact variants in a specific genomic region
#' high_impact_vars <- pg_query_by_impact(
#'   impact_level = "HIGH",
#'   chrom = "Chr05",
#'   start = 75104537,
#'   end = 75106403
#' )
#' print(high_impact_vars)
#'
#' }
#' 
#' @importFrom httr GET content http_error
#' @export
pg_query_by_impact <- function(
  impact_level = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
  chrom = NULL,
  start = NULL,
  end = NULL
) {
  # Squash the R vector into a single comma-separated string for the API
  impact_str <- paste(impact_level, collapse = ",")

  query_params <- list(
    impact_level = impact_str,
    chrom = chrom,
    start = start,
    end = end
  )

  return(.api_fetch("/db/impact", query = query_params))
}




#' Query variants by allele frequency (online)
#'
#' @param min_af Numeric. Minimum alternate allele frequency threshold for filtering. 
#'   Must fall between 0 and 1. Defaults to \code{0}.
#' @param max_af Numeric. Maximum alternate allele frequency threshold for filtering. 
#'   Must fall between 0 and 1. Defaults to \code{1}.
#' @param chrom A character value specifying the target chromosome name (e.g., \code{"Chr05"}).
#' @param start Integer. Start coordinate for the target window region.
#' @param end Integer. End coordinate for the target window region.
#'
#' @returns A data frame containing variant metadata and the calculated allele frequencies, 
#'   filtered by the specified thresholds.
#' 
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Query variants from the online database for a specific region
#' # and filter for common variants (MAF >= 5%)
#' common_variants_online <- pg_query_by_af(
#'   min_af = 0.05,
#'   max_af = 0.95,
#'   chrom = "Chr05",
#'   start = 75104537,
#'   end = 75106403
#' )
#'
#' # Print the head of the resulting data frame
#' print(head(common_variants_online))
#' }
#'
#' @export
pg_query_by_af <- function(
  min_af = 0,
  max_af = 1,
  chrom = NULL,
  start = NULL,
  end = NULL
) {
  if (is.null(chrom)) {
    stop("Chromosome must be specified to prevent memory overflow.")
  }

  # Get genotype matrix
  gt <- pg_query_db(
    table_name = "genotypes",
    chrom = chrom,
    start = start,
    end = end
  )

  if (is.null(gt) || nrow(gt) == 0) {
    message("No genotype data found for the specified region.")
    return(data.frame())
  }

  # Compute allele frequencies
  result <- pg_calc_af(
    gt,
    variant_id_col = 'variant_id',
    chrom_col = 'chrom',
    pos_col = 'pos'
  )

  # Filter frequencies by threshold stated
  result <- result[result$alt_af >= min_af & result$alt_af <= max_af, ]

  return(result)
}




#' Query genotypes for specific variant IDs (online)
#'
#' @param variant_ids A character vector of unique variant identifiers to look up (e.g., \code{c("SNP_Chr05_75104557")}).
#' @param variant_id_col A character value specifying the primary key column name for variant IDs. 
#'   Defaults to \code{"variant_id"}.
#' @param variants_table A character value specifying the name of the variants reference table. 
#'   Defaults to \code{"variants"}.
#' @param genotypes_table A character value specifying the name of the genotype call data table. 
#'   Defaults to \code{"genotypes"}.
#' @param meta_data A character vector of metadata columns to return alongside 
#'   the genotype matrix. This can include any field from the `variants` table 
#'   (e.g., `"chrom"`, `"pos"`, `"ref"`, `"alt"`, `"variant_type"`) as well as 
#'   dynamically calculated metrics: `"major_allele"`, `"minor_allele"`, 
#'   `"major_allele_freq"`, and `"minor_allele_freq"`. If `NULL` (default), 
#'   all available database fields and all calculated metrics are returned.
#'
#' @returns A data frame in wide format (variants x samples) with the
#'   requested metadata columns.
#' 
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Define a list of target variant IDs
#' target_markers <- c("INDEL_Chr05_75104541", "SNP_Chr05_75104557")
#'
#' # Query genotypes for these specific variant IDs from the online database,
#' # requesting specific metadata columns
#' pg_query_genotypes_result <- pg_query_genotypes(
#'   variant_ids = target_markers, 
#'   meta_data = c("chrom", "pos", "ref", "alt", "variant_type", 
#'                 "minor_allele", "major_allele", "major_allele_freq", "minor_allele_freq")
#' )
#'
#' # Print the head of the resulting data frame
#' print(head(pg_query_genotypes_result))
#' }
#'
#' @importFrom httr GET content http_error
#' @export
pg_query_genotypes <- function(
  variant_ids,
  variant_id_col = "variant_id",
  variants_table = 'variants',
  genotypes_table = 'genotypes',
  meta_data = NULL
) {
  # Prevent execution if no variant IDs are provided
  if (missing(variant_ids) || length(variant_ids) == 0) {
    warning("The 'variant_ids' vector is empty. Returning empty data frame.")
    return(data.frame())
  }

  # Predefine user request parameters
  wants_all_data <- is.null(meta_data)
  user_requested_cols <- meta_data
  calculated_metric_cols <- c(
    "major_allele",
    "minor_allele",
    "major_allele_freq",
    "minor_allele_freq"
  )

  # Formulate the metadata columns to safely request from the API
  if (!wants_all_data) {
    # Remove dynamic metrics from the API query since they don't exist in the DB table
    user_db_cols <- setdiff(user_requested_cols, calculated_metric_cols)

    # Force include columns required for table joins and allele metric calculations
    api_meta_cols <- unique(c(variant_id_col, "ref", "alt", user_db_cols))
    meta_str <- paste(api_meta_cols, collapse = ",")
  } else {
    meta_str <- NULL
  }

  # Convert R vector to a comma-separated string for the API query
  ids_str <- paste(variant_ids, collapse = ",")

  query_params <- list(
    variant_ids = ids_str,
    variant_id_col = variant_id_col,
    variants_table = variants_table,
    genotypes_table = genotypes_table,
    meta_data = meta_str
  )

  # Fetch data matrix from the API endpoint
  geno_res <- .api_fetch("/db/query_genotypes", query = query_params) 

  # Expand 
 # geno_res <- expand_matrix(geno_res)

  # Check if the API returned an empty result
  if (is.null(geno_res) || nrow(geno_res) == 0) {
    warning("No data found for the provided variant IDs.")
    return(data.frame())
  }

  # Process columns and compute metrics exactly like the offline counterpart
  if (wants_all_data) {
    # Dynamically extract metadata columns based on standard expected variant fields
    standard_metadata_fields <- c(
      variant_id_col,
      "variant_id",
      "chrom",
      "pos",
      "ref",
      "alt",
      "variant_type"
    )
    api_meta_cols <- intersect(colnames(geno_res), standard_metadata_fields)

    # Append metrics using all available metadata boundaries
    result_with_metrics <- append_locus_allele_metrics(
      geno_res,
      meta_data = api_meta_cols
    )
    return(result_with_metrics)
  } else {
    # Append metrics using the targeted metadata subset boundaries
    result_with_metrics <- append_locus_allele_metrics(
      geno_res,
      meta_data = api_meta_cols
    )

    # Isolate sample genotype columns (anything returned by the API that isn't a metadata field)
    sample_cols <- setdiff(colnames(geno_res), api_meta_cols)

    # Piece together final output layout: user-requested metadata followed by sample data
    final_cols_to_return <- c(variant_id_col,user_requested_cols, sample_cols)
    final_cols_to_return <- intersect(
      final_cols_to_return,
      colnames(result_with_metrics)
    )

    return(result_with_metrics[, final_cols_to_return, drop = FALSE])
  }
}




#' Count the distribution of variant types in the database (online)
#'
#' @param variants_table Character. The name of the table containing variant
#'   metadata. Defaults to "variants".
#'
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Count the distribution of variant types in the online database
#' pg_count_variant_types_result <- pg_count_variant_types()
#' print(pg_count_variant_types_result)
#' }
#'
#'
#' @returns A data frame with two columns: \code{variant_type} and \code{n}.
#' 
#' @importFrom httr GET content http_error
#' @export
pg_count_variant_types <- function(variants_table = "variants") {
  params <- list(variants_table = variants_table)

  return(.api_fetch("/stats/variant_types", query = params))
}




#' Summarize functional annotations and impacts for a genomic region (online)
#'
#' This function connects to the public panGenomeBreedr API to query variants within
#' a specific genomic range and returns summaries of SnpEff annotations and
#' impact levels, cross-tabulated by variant type.
#'
#' @param chrom Character. Chromosome name (e.g., "Chr05").
#' @param start Numeric. Start coordinate of the region.
#' @param end Numeric. End coordinate of the region.
#' @param annotations_table Character. Name of the annotations table. Defaults to "annotations".
#' @param variants_table Character. Name of the variants table. Defaults to "variants".
#'
#' @returns A list containing three data frames: \code{annotation_summary},
#'   \code{impact_summary}, and \code{variant_type_totals}.
#' 
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#' 
#' # Extract functional distribution summaries across a specific gene locus window
#' query_ann_summary_result <- pg_query_ann_summary(chrom = "Chr05",
#'                                               start = 75104537,
#'                                               end = 75106403)
#' print(query_ann_summary_result)
#' 
#' }
#' 
#' @export
#' @importFrom httr GET content http_error
pg_query_ann_summary <- function(
  chrom,
  start,
  end,
  annotations_table = "annotations",
  variants_table = "variants"
) {
  if (missing(chrom) || missing(start) || missing(end)) {
    stop("You must provide 'chrom', 'start', and 'end' for the region summary.")
  }

  query_params <- list(
    chrom = chrom,
    start = start,
    end = end,
    annotations_table = annotations_table,
    variants_table = variants_table
  )

  # Return as a list of dataframes
  result <- .api_fetch(
    "/stats/ann_summary",
    query = query_params,
    simplify = FALSE
  )

  if (is.null(result) || !is.data.frame(result$variant_type_totals) || nrow(result$variant_type_totals) == 0) {
    warning("No variants found in the specified region from the API.")
    return(list(
      annotation_summary = data.frame(),
      impact_summary = data.frame(),
      variant_type_totals = data.frame()
    ))
  }

  # Ensure all components are data frames before returning, as jsonlite can return empty lists
  result$annotation_summary <- as.data.frame(result$annotation_summary)
  result$impact_summary <- as.data.frame(result$impact_summary)
  result$variant_type_totals <- as.data.frame(result$variant_type_totals)

  return(result)
}




#' Get sample metadata from the database (online)
#'
#' @param query_col Character. Optional metadata column name to filter by (e.g., \code{"countryorigin"}). 
#'   If \code{NULL}, all database metadata records are returned.
#' @param query_value Character or Numeric. The explicit attribute value to match within \code{query_col}.
#'
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Retrieve all sample metadata from the online database
#' all_metadata <- pg_get_sample_metadata()
#' print(head(all_metadata))
#'
#' # Retrieve sample metadata filtered by 'countryorigin' = "Ghana"
#' ghana_metadata <- pg_get_sample_metadata(
#'   query_col = "countryorigin",
#'   query_value = "Ghana"
#' )
#' print(head(ghana_metadata))
#' }
#'
#'
#' @returns A data frame containing the sample metadata records.
#' @export
#' @importFrom httr GET content http_error
pg_get_sample_metadata <- function(query_col = NULL, query_value = NULL) {
  query_params <- list(query_col = query_col, query_value = query_value)
  return(.api_fetch("/db/metadata", query = query_params))
}



#' Filter genotypes by multiple sample metadata criteria
#'
#' @param genotype_matrix A data frame or matrix where rows are variants and
#'   columns include variant metadata followed by sample genotypes.
#' @param genotype_start_col An integer specifying the column index where sample
#'   genotype calls  begins.
#' @param filters A named list of metadata criteria to filter samples by.
#'   Names must match metadata columns, and values are the allowed entries.
#'   Example: \code{list(population = "Gates", countryorigin = c("Ethiopia", "Ghana", "Togo"))}
#'
#' @return A data frame containing the filtered genotype matrix with recalculated
#'   allele metrics for the sub-population. Returns an empty data frame if no
#'   samples match all combined criteria.
#'
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # Define filtering criteria
#' my_filters <- list(
#'   population = "Gates",
#'   countryorigin = c("Ethiopia", "Ghana", "Togo")
#' )
#'
#' # Extract the full genotype matrix first for a genomic region
#' genotype_data_region <- pg_query_db(
#'   table_name = "genotypes",
#'   chrom = "Chr05",
#'   start = 75104537,
#'   end = 75106403
#' )
#'
#' # Define the genotype_start_col.
#' genotype_start_col_val <- 11
#'
#' # Get filtered genotypes matrix
#' filtered_genotypes <- pg_query_geno_by_meta(
#'   genotype_matrix = genotype_data_region,
#'   genotype_start_col = genotype_start_col_val,
#'   filters = my_filters
#' )
#' print(head(filtered_genotypes))
#' }
#'
#' @export
#' @importFrom httr GET content http_error
pg_query_geno_by_meta <- function(genotype_matrix, genotype_start_col, filters) {

  # Validate function inputs to ensure they meet requirements.
  if (!is.data.frame(genotype_matrix) && !is.matrix(genotype_matrix)) {
    stop("`genotype_matrix` must be a data frame or matrix.", call. = FALSE)
  }
  if (nrow(genotype_matrix) == 0) {
    warning("Input `genotype_matrix` is empty.", call. = FALSE)
    return(data.frame())
  }
  if (!is.numeric(genotype_start_col) || genotype_start_col <= 1) {
    stop("`genotype_start_col` must be a numeric value greater than 1.", call. = FALSE)
  }
  if (!is.list(filters) || is.null(names(filters))) {
    stop("`filters` must be a named list. Example: list(population = 'Gates')", call. = FALSE)
  }

  # Retrieve the complete sample metadata table.
  # Note: Ensure pg_get_sample_metadata() is set up to return the full table if no args are passed.
  full_metadata <- pg_get_sample_metadata()

  if (is.null(full_metadata) || nrow(full_metadata) == 0) {
    message("No metadata could be retrieved from the database.")
    return(data.frame())
  }

  # Iteratively filter the metadata based on the provided list of criteria.
  filtered_metadata <- full_metadata
  for (col_name in names(filters)) {
    
    # Skip gracefully if the user requests a column that does not exist.
    if (!col_name %in% colnames(filtered_metadata)) {
      warning(sprintf("Column '%s' not found in metadata. Skipping this filter.", col_name), call. = FALSE)
      next
    }
    
    # Subset the metadata to keep only rows matching the current list values.
    allowed_values <- filters[[col_name]]
    filtered_metadata <- filtered_metadata[filtered_metadata[[col_name]] %in% allowed_values, , drop = FALSE]
  }

  # Stop execution if the combined filters eliminated all samples.
  if (nrow(filtered_metadata) == 0) {
    message("No samples found matching all combined metadata criteria.")
    return(data.frame())
  }

  # Extract the specific library IDs to keep.
  sample_ids_to_keep <- if ("lib" %in% colnames(filtered_metadata)) {
    unique(filtered_metadata$lib)
  } else {
    warning("Metadata is missing the 'lib' column. Cannot identify samples to keep.", call. = FALSE)
    return(data.frame())
  }

  # Subset the genotype matrix to the filtered samples.
  all_matrix_cols <- colnames(genotype_matrix)
  variant_meta_cols <- all_matrix_cols[seq_len(genotype_start_col - 1)]
  
  # Identify which of the desired sample columns are actually in the input matrix.
  available_sample_cols <- intersect(sample_ids_to_keep, all_matrix_cols)

  if (length(available_sample_cols) == 0) {
    message("None of the matching samples were found in the provided genotype matrix.")
    return(data.frame())
  }

  # Create the sub-population matrix using only available columns.
  sub_population_matrix <- genotype_matrix[, c(variant_meta_cols, available_sample_cols), drop = FALSE]

  # Recalculate allele metrics for the filtered sub-population.
  base_meta_cols <- c("variant_id", "chrom", "pos", "variant_type", "ref", "alt")
  
  # Ensure the calculation function only receives metadata it can use.
  meta_for_recalc <- intersect(base_meta_cols, variant_meta_cols)

  final_result <- append_locus_allele_metrics(
    region_matrix = sub_population_matrix,
    meta_data = meta_for_recalc
  )

  return(final_result)
}





#' Get the genomic range of a candidate gene using the Sobic ID from a GFF file (online)
#'
#' @param gene_name A character value indicating the Sobic ID of the candidate gene.
#' @param gff_path A character value specifying the path to the GFF3 file.
#' URL paths and GZ-compressed files are supported.
#'
#' @returns A list containing the chromosome, start, and end coordinates.
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#'
#' # Define path to a sorghum GFF3 file (v5.1)
#' gff_url <- "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz"
#'
#' # Retrieve coordinates for a candidate gene
#' coords <- pg_gene_coords(gene_name = "Sobic.005G213600", gff_path = gff_url)
#' print(coords)
#' }
#'
#'
#' @importFrom R.utils isUrl isGzipped gunzip
#' @export
pg_gene_coords <- function(gene_name, gff_path) {
  # We handle remote URLs by downloading the file to a temporary location first.
  if (R.utils::isUrl(gff_path)) {
    local_temp <- file.path(tempdir(), basename(gff_path))
    utils::download.file(
      gff_path,
      destfile = local_temp,
      mode = "wb",
      quiet = TRUE
    )

    # Ensure the downloaded file is removed when the function finishes to keep temp clean.
    on.exit(if (file.exists(local_temp)) unlink(local_temp), add = TRUE)
    target_path <- local_temp
  } else {
    target_path <- gff_path
  }

  # If the file is compressed, we decompress it to a temporary GFF file before reading.
  if (R.utils::isGzipped(target_path)) {
    decompressed_gff <- file.path(tempdir(), "temp_genomic_data.gff3")
    R.utils::gunzip(
      target_path,
      destname = decompressed_gff,
      overwrite = TRUE,
      remove = FALSE
    )

    # Clean up the decompressed file once the data is loaded into memory.
    on.exit(
      if (file.exists(decompressed_gff)) unlink(decompressed_gff),
      add = TRUE
    )
    read_from <- decompressed_gff
  } else {
    read_from <- target_path
  }

  # Load the GFF file. We use read.delim with a specific comment character to skip headers.
  gff <- utils::read.delim(
    read_from,
    comment.char = "#",
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  # Assign standard GFF3 column headers for logical filtering
  colnames(gff) <- c(
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes"
  )

  # Filter for gene features and extract the name attribute using regex.
  # This looks for the Name= or ID= tag within the attributes column.
  genes <- gff[gff$type == "gene", ]
  genes$extracted_id <- sub(".*(?:ID|Name)=([^;]+).*", "\\1", genes$attributes)

  # Locate the specific gene match
  gmatch <- genes[genes$extracted_id == gene_name, ]

  if (nrow(gmatch) == 0) {
    stop(paste("Gene ID", gene_name, "was not found in the provided GFF file."))
  }

  if (nrow(gmatch) > 1) {
    warning(
      "Multiple gene entries found; returning coordinates for the first occurrence."
    )
  }

  # Return the genomic range as a structured list for use in downstream database queries
  return(list(
    chrom = gmatch$seqid[1],
    start = as.integer(gmatch$start[1]),
    end = as.integer(gmatch$end[1])
  ))
}


#' Compute allele frequencies for a genotype matrix (online)
#'
#' @param gt A data frame or matrix where rows are variants and columns are samples.
#'   May also contain metadata columns (ID, Chrom, Pos).
#' @param variant_id_col Character. Name of the column containing variant IDs.
#'   Default is 'variant_id'.
#' @param chrom_col Character. Optional name of the chromosome column.
#' @param pos_col Character. Optional name of the position column.
#'
#' @returns A data frame containing the variant metadata and two calculated columns:
#' \itemize{
#'   \item \code{ref_af}: Reference allele frequency.
#'   \item \code{alt_af}: Alternate allele frequency.
#' }
#' @examples
#' \dontrun{
#' # Load the package
#' library(panGenomeBreedr)
#'
#' # First, query genotype data for a region from the online database
#' # This will be the 'gt' (genotype matrix) input for pg_calc_af
#' gt_data_region <- pg_query_db(
#'   table_name = "genotypes",
#'   chrom = "Chr05",
#'   start = 75104537,
#'   end = 75106403
#' )
#'
#' # Calculate allele frequencies for the fetched genotype data
#' af_metrics <- pg_calc_af(
#'   gt = gt_data_region,
#'   variant_id_col = "variant_id",
#'   chrom_col = "chrom",
#'   pos_col = "pos"
#' )
#' print(head(af_metrics))
#' }
#' 
#' @export
pg_calc_af <- function(
  gt,
  variant_id_col = 'variant_id',
  chrom_col = NULL,
  pos_col = NULL
) {
  # Ensure the required variant ID column is present before subsetting
  if (!variant_id_col %in% colnames(gt)) {
    stop(paste(
      "Metadata column",
      variant_id_col,
      "not found in the input data frame."
    ))
  }

  # Identify sample columns by stripping away specified metadata
  meta_cols <- c(variant_id_col, chrom_col, pos_col)
  sample_cols <- setdiff(colnames(gt), meta_cols)

  if (length(sample_cols) == 0) {
    stop("No sample columns detected. Check your column name specifications.")
  }

  # Convert the genotype portion to a matrix for high-speed vectorized operations
  gt_mat <- as.matrix(gt[, sample_cols])

  # We use a double-precision matrix to store numeric dosages for frequency math
  dosage_mat <- matrix(NA_real_, nrow = nrow(gt_mat), ncol = ncol(gt_mat))

  # Vectorized assignment for standard VCF genotype patterns
  dosage_mat[gt_mat == "0|0" | gt_mat == "0/0"] <- 0
  dosage_mat[gt_mat %in% c("0|1", "1|0", "0/1", "1/0")] <- 1
  dosage_mat[gt_mat == "1|1" | gt_mat == "1/1"] <- 2

  # Compute alternate allele frequency
  alt_af <- rowMeans(dosage_mat, na.rm = TRUE) / 2

  # Replace NaNs with NAs for downstream consistency
  alt_af[is.nan(alt_af)] <- NA

  cols_to_keep <- intersect(meta_cols, colnames(gt))
  out_meta <- gt[, cols_to_keep, drop = FALSE]

  result <- cbind(
    out_meta,
    data.frame(
      ref_af = 1 - alt_af,
      alt_af = alt_af,
      stringsAsFactors = FALSE
    )
  )

  return(result)
}


#' Filter genotypes by allele frequency (online)
#'
#' @param gt A data frame containing genotype calls queried from the genotypes table.
#' @param variant_id_col Character. The column name in \code{gt} matching unique variant identifiers. 
#'   Defaults to \code{"variant_id"}.
#' @param chrom_col Character. The column name in \code{gt} matching chromosome names. 
#'   Defaults to \code{"chrom"}.
#' @param pos_col Character. The column name in \code{gt} matching genomic positions. 
#'   Defaults to \code{"pos"}.
#' @param min_af Numeric. Minimum alternate allele frequency threshold for filtering. 
#'   Must fall between 0 and 1. Defaults to \code{0}.
#' @param max_af Numeric. Maximum alternate allele frequency threshold for filtering. 
#'   Must fall between 0 and 1. Defaults to \code{1}.
#'
#' @returns A data frame containing variant metadata and calculated frequencies
#'   for variants that passed the filter.
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#'
#' # Query genotypes table and filter to remove rare variants (MAF < 0.05)
#' filtered_vars <- pg_query_db(
#'   table_name = "genotypes",
#'   chrom = "Chr05",
#'   start = 75104537,
#'   end = 75106403
#' ) |>
#'   pg_filter_by_af(min_af = 0.05, max_af = 0.95)
#' }
#'
#' @export
pg_filter_by_af <- function(
  gt,
  variant_id_col = 'variant_id',
  chrom_col = 'chrom',
  pos_col = 'pos',
  min_af = 0,
  max_af = 1
) {
  # QC checks
  if (is.null(gt) || nrow(gt) == 0) {
    stop("The genotype matrix is empty or NULL.")
  }

  # compute allele frequencies (Using your local function)
  af_results <- pg_calc_af(
    gt = gt,
    variant_id_col = variant_id_col,
    chrom_col = chrom_col,
    pos_col = pos_col
  )

  # Filter out
  pass_idx <- which(af_results$alt_af >= min_af & af_results$alt_af <= max_af)

  result <- af_results[pass_idx, , drop = FALSE]

  # Alert the user if the thresholds were too strict, resulting in an empty set
  if (nrow(result) == 0) {
    warning(
      "No variants passed the allele frequency filter. Check your min_af/max_af thresholds."
    )
  }

  return(result)
}


#' Interactive geographic exploration of the samples (online)
#'
#' @param metadata A data frame containing sample metadata. Must include 'lat' and 'lon'
#'   columns, along with the specified coloring column.
#' @param color_by Character. The metadata column to use for point coloration.
#'   Defaults to "countryorigin".
#'
#' @returns A \code{leaflet} map object (htmlwidget) representing the interactive map.
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#'
#' # Fetch sample metadata from the database
#' meta <- pg_get_sample_metadata()
#'
#' # Explore the geographic distribution colored by genetic cluster
#' pg_map_accessions(meta, color_by = "countryorigin")
#' }
#'
#' @export
pg_map_accessions <- function(metadata, color_by = "countryorigin") {
  # Dependency chceck
  if (!requireNamespace("leaflet", quietly = TRUE)) {
    stop(
      "The 'leaflet' package is required to create interactive maps. ",
      "Please install it by running: install.packages('leaflet')",
      call. = FALSE
    )
  }
  if (!requireNamespace("tools", quietly = TRUE)) {
    stop(
      "The 'tools' package is required. Please run: install.packages('tools')",
      call. = FALSE
    )
  }

  # Check if coordinate columns exist
  if (!all(c("lat", "lon") %in% names(metadata))) {
    stop(
      "Metadata must contain 'lat' and 'lon' columns for geographic mapping."
    )
  }

  metadata$lat <- suppressWarnings(as.numeric(metadata$lat))
  metadata$lon <- suppressWarnings(as.numeric(metadata$lon))

  # Filter for complete cases now that everything is properly numeric/NA
  plot_data <- metadata[!is.na(metadata$lat) & !is.na(metadata$lon), ]

  if (nrow(plot_data) == 0) {
    stop("No samples with valid latitude and longitude found in the metadata.")
  }

  # Isolating all columns except coordinates to build the dynamic popups
  meta_cols <- setdiff(names(plot_data), c("lat", "lon"))

  # HTML popups
  html_components <- lapply(meta_cols, function(col) {
    clean_title <- tools::toTitleCase(gsub("_", " ", col))
    ifelse(
      is.na(plot_data[[col]]),
      "",
      paste0("<b>", clean_title, ":</b> ", plot_data[[col]], "<br>")
    )
  })

  popup_info <- paste0(
    "<div style='font-family: Arial, sans-serif; font-size: 12px; max-height: 250px; overflow-y: auto;'>",
    do.call(paste0, html_components),
    "</div>"
  )

  # Create a dynamic color palette
  unique_vals <- unique(plot_data[[color_by]])

  pal <- leaflet::colorFactor(
    palette = "viridis",
    domain = unique_vals,
    na.color = "transparent"
  )

  # Construct the Leaflet map
  map <- leaflet::leaflet(data = plot_data) |>
    leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron) |>
    leaflet::addCircleMarkers(
      lng = ~lon,
      lat = ~lat,
      radius = 5,
      color = ~ pal(plot_data[[color_by]]),
      stroke = FALSE,
      fillOpacity = 0.8,
      popup = popup_info,
      label = ~lib
    )

  return(map)
}