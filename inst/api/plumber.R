# plumber.R
# library(RPostgres)
# library(DBI)
# library(plumber)

# --- GLOBAL SETUP & CONNECTION MANAGEMENT ---
MASTER_CON <- NULL

ensure_con <- function() {
  if (is.null(MASTER_CON) || !DBI::dbIsValid(MASTER_CON)) {
    MASTER_CON <<- pgsql_connect(
      host = Sys.getenv("PGSQL_HOST"),
      dbname = Sys.getenv("PGSQL_DBNAME", "postgres"),
      user = Sys.getenv("PGSQL_USER", "postgres"),
      password = Sys.getenv("PGSQL_PASS")
    )
  }
  return(MASTER_CON)
}




#* @apiTitle panGenomeBreedr REST API
#* @apiDescription Public programmatic access to the Sorghum bicolor pangenome.
#* @apiContact list(name = "Israel Tetteh", email = "your.email@knust.edu.gh")

# ---  DATABASE METADATA ENDPOINTS ---

#* List all tables in the database
#* @get /db/tables
function() {
  con <- ensure_con()
  pgsql_list_tables(con)
}




#* Summarize row counts for all tables
#* @get /db/summary
function() {
  con <- ensure_con()
  pgsql_summarize_tables(con)
}




#* Get column metadata for a specific table
#* @param table_name Name of the table (variants, annotations, genotypes, metadata)
#* @get /db/columns
function(table_name = "variants") {
  con <- ensure_con()
  pgsql_list_table_columns(con, table_name)
}




# ---  VARIANT & ANNOTATION STATS ENDPOINTS ---

#* Get high-level variant statistics per chromosome
#* @param include_annotations Include annotation counts? (true/false)
#* @get /stats/variants
function(include_annotations = "true") {
  con <- ensure_con()
  inc_ann <- as.logical(include_annotations)
  pgsql_variant_stats(con, include_annotations = inc_ann)
}




#* Get variant impact summary
#* @get /stats/impact
function() {
  con <- ensure_con()
  pgsql_variant_impact_summary(con)
}




#* Count the distribution of variant types
#* @param variants_table The name of the metadata table (default: variants)
#* @get /stats/variant_types
function(variants_table = "variants") {
  con <- ensure_con()

  # Safely call your hidden SQL engine
  pgsql_count_variant_types(con, variants_table = variants_table)
}




#* Summarize genomic annotations and impacts in a specific region
#* @param chrom Chromosome name (e.g., Chr05)
#* @param start Genomic start coordinate
#* @param end Genomic end coordinate
#* @param annotations_table Table name for annotations
#* @param variants_table Table name for variants
#* @get /stats/ann_summary
function(
  chrom,
  start,
  end,
  annotations_table = "annotations",
  variants_table = "variants"
) {
  con <- ensure_con()

  # Helper to safely handle missing URL parameters
  clean_param <- function(x) {
    if (is.null(x) || x == "" || x == "NA") {
      return(NULL)
    }
    return(x)
  }

  chrom_clean <- clean_param(chrom)

  # Convert position coordinates safely to numeric
  start_num <- clean_param(start)
  if (!is.null(start_num)) {
    start_num <- as.numeric(start_num)
  }

  end_num <- clean_param(end)
  if (!is.null(end_num)) {
    end_num <- as.numeric(end_num)
  }

  # Safely call your hidden SQL engine
  pgsql_query_ann_summary(
    con = con,
    chrom = chrom_clean,
    start = start_num,
    end = end_num,
    annotations_table = annotations_table,
    variants_table = variants_table
  )
}




# ---  DATA QUERY ENDPOINTS ---

#* Query pangenome tables by genomic coordinates
#* @param table_name Which table to query
#* @param chrom Chromosome name
#* @param start Genomic start position
#* @param end Genomic end position
#* @param gene_name Specific gene name to filter annotations
#* @get /db/query
function(
  table_name = "variants",
  chrom = NULL,
  start = NULL,
  end = NULL,
  gene_name = NULL
) {
  con <- ensure_con()

  # Helper to safely handle missing URL parameters
  clean_param <- function(x) {
    if (is.null(x) || x == "" || x == "NA") {
      return(NULL)
    }
    return(x)
  }

  # Clean inputs
  chrom_clean <- clean_param(chrom)
  gene_clean <- clean_param(gene_name)

  # Convert position coordinates back to numeric if they exist
  start_num <- clean_param(start)
  if (!is.null(start_num)) {
    start_num <- as.numeric(start_num)
  }

  end_num <- clean_param(end)
  if (!is.null(end_num)) {
    end_num <- as.numeric(end_num)
  }

  # Call your hidden SQL engine
  pgsql_query_db(
    con = con,
    table_name = table_name,
    chrom = chrom_clean,
    start = start_num,
    end = end_num,
    gene_name = gene_clean
  )
}




#* Extract variants based on mutation impact
#* @param impact_level Comma-separated list of impact levels (e.g., HIGH,MODERATE)
#* @param chrom Chromosome name
#* @param start Genomic start position
#* @param end Genomic end position
#* @get /db/impact
function(
  impact_level = "HIGH,MODERATE,LOW,MODIFIER",
  chrom = NULL,
  start = NULL,
  end = NULL
) {
  con <- ensure_con()

  # Helper to safely handle missing URL parameters
  clean_param <- function(x) {
    if (is.null(x) || x == "" || x == "NA") {
      return(NULL)
    }
    return(x)
  }

  # Re-build the R vector from the comma-separated URL string
  impact_vec <- strsplit(impact_level, ",")[[1]]
  impact_vec <- trimws(impact_vec) # Remove any accidental spaces

  # Clean other inputs
  chrom_clean <- clean_param(chrom)

  start_num <- clean_param(start)
  if (!is.null(start_num)) {
    start_num <- as.numeric(start_num)
  }

  end_num <- clean_param(end)
  if (!is.null(end_num)) {
    end_num <- as.numeric(end_num)
  }

  # Safely call your hidden SQL engine
  pgsql_query_by_impact(
    con = con,
    impact_level = impact_vec,
    chrom = chrom_clean,
    start = start_num,
    end = end_num
  )
}




#* Query variants by allele frequency within a genomic region
#* @param min_af Minimum alternate allele frequency (0-1)
#* @param max_af Maximum alternate allele frequency (0-1)
#* @param chrom Chromosome name
#* @param start Genomic start position
#* @param end Genomic end position
#* @get /db/query_by_af
function(min_af = 0, max_af = 1, chrom = NULL, start = NULL, end = NULL) {
  con <- ensure_con()

  # Helper to safely handle missing URL parameters
  clean_param <- function(x) {
    if (is.null(x) || x == "" || x == "NA") {
      return(NULL)
    }
    return(x)
  }

  # Clean inputs
  chrom_clean <- clean_param(chrom)

  # Convert numeric parameters safely
  start_num <- clean_param(start)
  if (!is.null(start_num)) {
    start_num <- as.numeric(start_num)
  }

  end_num <- clean_param(end)
  if (!is.null(end_num)) {
    end_num <- as.numeric(end_num)
  }

  min_af_num <- as.numeric(min_af)
  max_af_num <- as.numeric(max_af)

  # Safely call your hidden hybrid engine
  pgsql_query_by_af(
    con = con,
    min_af = min_af_num,
    max_af = max_af_num,
    chrom = chrom_clean,
    start = start_num,
    end = end_num
  )
}




#* Query genotypes for specific variant IDs
#* @param variant_ids Comma-separated list of variant IDs
#* @param variant_id_col Column name for ID
#* @param variants_table Table name for variants
#* @param genotypes_table Table name for genotypes
#* @param meta_data Comma-separated list of metadata columns to keep
#* @get /db/query_genotypes
function(
  variant_ids = "",
  variant_id_col = "variant_id",
  variants_table = "variants",
  genotypes_table = "genotypes",
  meta_data = NULL
) {
  con <- ensure_con()

  if (variant_ids == "") {
    return(data.frame())
  }

  # Re-build R vectors from the URL strings
  ids_vec <- trimws(strsplit(variant_ids, ",")[[1]])

  meta_vec <- NULL
  if (!is.null(meta_data) && meta_data != "") {
    meta_vec <- trimws(strsplit(meta_data, ",")[[1]])
  }

  # Safely call your hidden SQL engine
  pgsql_query_genotypes(
    con = con,
    variant_ids = ids_vec,
    variant_id_col = variant_id_col,
    variants_table = variants_table,
    genotypes_table = genotypes_table,
    meta_data = meta_vec
  )
}




# ---  SAMPLE METADATA ENDPOINTS ---

#* Retrieve sample metadata
#* @param query_col The metadata column to filter by (optional)
#* @param query_value The value to match (optional)
#* @get /db/metadata
function(query_col = NULL, query_value = NULL) {
  con <- ensure_con()

  # Helper to safely handle missing URL parameters
  clean_param <- function(x) {
    if (is.null(x) || x == "" || x == "NA") {
      return(NULL)
    }
    return(x)
  }

  col_clean <- clean_param(query_col)
  val_clean <- clean_param(query_value)

  # Safely call your hidden SQL engine
  pgsql_get_sample_metadata(
    con = con,
    query_col = col_clean,
    query_value = val_clean
  )
}




#* Query genotypes filtered by sample metadata attributes
#* @param chrom Chromosome name
#* @param start Genomic start position
#* @param end Genomic end position
#* @param meta_col Metadata column to filter by
#* @param meta_value Value to match in the metadata column
#* @get /db/query_by_metadata
function(chrom, start, end, meta_col, meta_value) {
  con <- ensure_con()

  # Helper to safely handle missing URL parameters
  clean_param <- function(x) {
    if (is.null(x) || x == "" || x == "NA") {
      return(NULL)
    }
    return(x)
  }

  # Clean text inputs
  chrom_clean <- clean_param(chrom)
  col_clean <- clean_param(meta_col)
  val_clean <- clean_param(meta_value)

  # Convert position coordinates safely to numeric
  start_num <- clean_param(start)
  if (!is.null(start_num)) {
    start_num <- as.numeric(start_num)
  }

  end_num <- clean_param(end)
  if (!is.null(end_num)) {
    end_num <- as.numeric(end_num)
  }

  # Safely call your hidden SQL engine
  pgsql_query_by_metadata(
    con = con,
    chrom = chrom_clean,
    start = start_num,
    end = end_num,
    meta_col = col_clean,
    meta_value = val_clean
  )
}