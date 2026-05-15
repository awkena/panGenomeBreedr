#' Set the API endpoint URL for pangenome database
#'
#' Allows users to connect the package to any custom-hosted database via API.
#'
#' @param url Character. The full URL to your hosted API endpoint
#'   (e.g., "http://rice-genomics.myuniversity.edu:8000").
#'
#' @export
set_api_url <- function(url) {
  # Strip any accidental trailing slashes
  url <- sub("/$", "", url)

  # Set R options in current session
  options(panGenomeBreedr.api_url = url)

  Sys.setenv(PANGENOME_API_URL = url)

  message("panGenomeBreedr API endpoint successfully set to: ", url)
}





#' Get the current API endpoint URL
#'
#' Retrieves the active API endpoint. If no custom endpoint has been set
#' by the user, it safely defaults to the official Sorghum pangenome database API.
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
    url <- "http://16.171.142.87:8000"
  }

  return(url)
}




#' Internal API Fetcher
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




#' List all tables in the connected pangenome database.
#'
#' This function connects to the public panGenomeBreedr API and retrieves the names of
#' all tables within the database.
#'
#' @returns A character vector of table names.
#' @export
#' @importFrom httr GET content http_error
pg_list_tables <- function() {
  return(.api_fetch("/db/tables"))
}



#' Get variant statistics stored in the pangenome database.
#'
#' This function connects to the public panGenomeBreedr API and calculates
#' summary statistics for variants per chromosome, including variant counts
#' and genomic ranges.
#'
#' @param include_annotations A logical value indicating whether to include
#' statistics for the annotations table. Defaults to \code{TRUE}.
#'
#' @returns A data frame containing variant statistics.
#' @export
#' @importFrom httr GET content http_error
pg_variant_stats <- function(include_annotations = TRUE) {
  # Pass the logical parameter as a named list
  params <- list(include_annotations = include_annotations)

  return(.api_fetch("/stats/variants", query = params))
}




#' Get variant statistics stored in the pangenome database based on mutation impact.
#'
#' This function connects to the public panGenomeBreedr API and summarizes the
#' distribution of mutation impacts across chromosomes.
#'
#' @returns A data frame in wide format where each row is a chromosome and
#' columns represent the counts for each impact category.
#' @export
#' @importFrom httr GET content http_error
pg_variant_impact_summary <- function() {
  return(.api_fetch("/stats/impact"))
}




#' Name and row count for each table in the pangenome database.
#'
#' This function connects to the public panGenomeBreedr API and returns a
#' summary data frame containing the table names and their respective row counts.
#'
#' @returns A data frame with two columns: \code{table} and \code{n_rows}.
#' @export
#' @importFrom httr GET content http_error
pg_summarize_tables <- function() {
  return(.api_fetch("/db/summary"))
}




#' Check the column names and types for any table in the pangenome database.
#'
#' This function connects to the public panGenomeBreedr API and retrieves
#' metadata about the columns in a specified table.
#'
#' @param table_name A character value specifying the name of the table.
#' Defaults to "variants". Available options are "variants", "annotations",
#' "genotypes", or "metadata".
#'
#' @returns A data frame containing column metadata.
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




#' Query any table in the pangenome database using chromosome and a genomic position range.
#'
#' This function connects to the public panGenomeBreedr API to retrieve data
#' from the variants, annotations, or genotypes tables based on a specific
#' chromosome and genomic range.
#'
#' @param table_name Character. One of "variants", "annotations", or "genotypes".
#' @param chrom Character. The chromosome name (e.g., "Chr05").
#' @param start,end Numeric. The genomic start and end positions.
#' @param gene_name Character. Optional Sobic ID to filter annotations (e.g., "Sobic.005G213600").
#'
#' @returns A data frame containing the queried genomic data.
#' @export
#' @importFrom httr GET content http_error
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

  return(.api_fetch("/db/query", query = query_params))
}




#' Extract variants from the annotation table based on impact type.
#'
#' This function connects to the public panGenomeBreedr API to retrieve variants
#' filtered by snpEff impact levels (HIGH, MODERATE, etc.) and optional
#' genomic coordinates.
#'
#' @param impact_level Character vector. One or more of "HIGH", "MODERATE",
#'   "LOW", "MODIFIER". Defaults to all four.
#' @param chrom Character. Optional chromosome name (e.g., "Chr05").
#' @param start,end Numeric. Optional genomic start and end positions.
#'
#' @returns A data frame containing variant information and associated functional
#'   annotations.
#' @export
#' @importFrom httr GET content http_error
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




#' Extract variants based on minimum and maximum allele frequencies within a defined region in the pangenome database.
#'
#' This function connects to the public panGenomeBreedr API to query genotypes
#' within a specific genomic range and filters the results to only include
#' variants within the specified alternate allele frequency (AF) thresholds.
#'
#' @param min_af Numeric. Minimum alternate allele frequency (0-1). Default is 0.
#' @param max_af Numeric. Maximum alternate allele frequency (0-1). Default is 1.
#' @param chrom Character. Chromosome name (e.g., "Chr05").
#' @param start,end Numeric. Genomic start and end coordinates.
#'
#' @returns A data frame containing variant metadata (ID, Chrom, Pos) and the
#'   calculated \code{ref_af} and \code{alt_af}.
#' @export
#' @importFrom httr GET content http_error
pg_query_by_af <- function(
  min_af = 0,
  max_af = 1,
  chrom = NULL,
  start = NULL,
  end = NULL
) {
  if (is.null(chrom)) {
    stop(
      "Chromosome ('chrom') must be specified for an AF query to prevent memory overflow."
    )
  }

  query_params <- list(
    min_af = min_af,
    max_af = max_af,
    chrom = chrom,
    start = start,
    end = end
  )

  return(.api_fetch("/db/query_by_af", query = query_params))
}




#' Query genotypes for one or more variant IDs from a wide-format genotype table.
#'
#' This function connects to the public panGenomeBreedr API to retrieve genomic
#' data for a specific list of variant IDs. It expands the genotype array
#' into a wide format (samples as columns).
#'
#' @param variant_ids A character vector of variant IDs to retrieve.
#' @param variant_id_col Character. Name of the ID column. Default is "variant_id".
#' @param variants_table Character. Name of the metadata table. Default is "variants".
#' @param genotypes_table Character. Name of the genotype table. Default is "genotypes".
#' @param meta_data Character vector. Specific columns to include from the
#'   variants table (e.g., "chrom", "pos", "ref", "alt"). If \code{NULL},
#'   retrieves all columns.
#'
#' @returns A data frame in wide format (variants x samples) with the
#'   requested metadata columns.
#' @export
#' @importFrom httr GET content http_error
pg_query_genotypes <- function(
  variant_ids,
  variant_id_col = "variant_id",
  variants_table = 'variants',
  genotypes_table = 'genotypes',
  meta_data = NULL
) {
  if (missing(variant_ids) || length(variant_ids) == 0) {
    warning("The 'variant_ids' vector is empty. Returning empty data frame.")
    return(data.frame())
  }

  # Convert R vectors to strings for the API
  ids_str <- paste(variant_ids, collapse = ",")
  meta_str <- if (!is.null(meta_data)) {
    paste(meta_data, collapse = ",")
  } else {
    NULL
  }

  query_params <- list(
    variant_ids = ids_str,
    variant_id_col = variant_id_col,
    variants_table = variants_table,
    genotypes_table = genotypes_table,
    meta_data = meta_str
  )

  return(.api_fetch("/db/query_genotypes", query = query_params))
}




#' Count the number of variant types in the pangenome database.
#'
#' This function connects to the public panGenomeBreedr API to perform a
#' server-side aggregation, counting the occurrences of different variant types
#' (e.g., SNP, INDEL) stored in the database.
#'
#' @param variants_table Character. The name of the table containing variant
#'   metadata. Defaults to "variants".
#'
#' @returns A data frame with two columns: \code{variant_type} and \code{n}.
#' @export
#' @importFrom httr GET content http_error
pg_count_variant_types <- function(variants_table = "variants") {
  params <- list(variants_table = variants_table)

  return(.api_fetch("/stats/variant_types", query = params))
}




#' Query the annotations table within a specified genomic region and summarize the distribution of SnpEff annotations and impact categories by variant type.
#'
#' This function connects to the public panGenomeBreedr API to query variants
#' within a specific genomic range and returns summaries of SnpEff annotations
#' and impact levels, cross-tabulated by variant type.
#'
#' @param chrom Character. Chromosome name (e.g., "Chr05").
#' @param start Numeric. Start coordinate of the region.
#' @param end Numeric. End coordinate of the region.
#' @param annotations_table Character. Name of the annotations table. Defaults to "annotations".
#' @param variants_table Character. Name of the variants table. Defaults to "variants".
#'
#' @returns A list containing three data frames: \code{annotation_summary},
#'   \code{impact_summary}, and \code{variant_type_totals}.
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
  return(.api_fetch(
    "/stats/ann_summary",
    query = query_params,
    simplify = FALSE
  ))
}




#' Retrieve sample metadata from the pangenome database.
#'
#' This function connects to the public panGenomeBreedr API to fetch accession-level
#' metadata, such as origin, race, and classification from the 'metadata' table. It supports optional
#' filtering by specific columns to easily subset populations for downstream analysis.
#'
#' @param query_col Character. The metadata column to filter by (e.g., "countryorigin").
#'   If \code{NULL}, all records are returned.
#' @param query_value Character or Numeric. The specific value to match in \code{query_col}.
#'
#' @returns A data frame containing the sample metadata records.
#' @export
#' @importFrom httr GET content http_error
pg_get_sample_metadata <- function(query_col = NULL, query_value = NULL) {
  query_params <- list(query_col = query_col, query_value = query_value)
  return(.api_fetch("/db/metadata", query = query_params))
}




#' Query genotypes filtered by sample metadata attributes in the pangenome database.
#'
#' This function connects to the public panGenomeBreedr API to retrieve genotypes for a specific
#' genomic region, but restricts the output to a subset of samples defined
#' by metadata attributes (e.g., extracting variants only for specific countries, populations, or phenotypic clusters).
#'
#' @param chrom Character. Chromosome name.
#' @param start Numeric. Genomic start range.
#' @param end Numeric. Genomic end range.
#' @param meta_col Character. Metadata column to filter samples by (e.g., "countryorigin").
#' @param meta_value Character. Value to match in \code{meta_col} (e.g., "Ghana").
#'
#' @returns A wide-format data frame (variants x filtered samples).
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#'
#' # Get genotypes for a gene, but only for samples from Ethiopia
#' eth_genotypes <- pg_query_by_metadata(
#'   chrom = "Chr05",
#'   start = 75104537,
#'   end = 75106403,
#'   meta_col = "countryorigin",
#'   meta_value = "Ethiopia"
#' )
#' }
#' @export
#' @importFrom httr GET content http_error
pg_query_by_metadata <- function(chrom, start, end, meta_col, meta_value) {
  if (
    missing(chrom) ||
      missing(start) ||
      missing(end) ||
      missing(meta_col) ||
      missing(meta_value)
  ) {
    stop(
      "All parameters (chrom, start, end, meta_col, meta_value) are required."
    )
  }

  query_params <- list(
    chrom = chrom,
    start = start,
    end = end,
    meta_col = meta_col,
    meta_value = meta_value
  )

  return(.api_fetch("/db/query_by_metadata", query = query_params))
}
