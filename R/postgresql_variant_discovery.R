#' Connect to the panGenomeBreedr Database
#'
#' @param host The database endpoint. Defaults to the 'PGSQL_HOST' environment variable.
#' @param dbname The database name. Defaults to the 'PGSQL_DBNAME' environment variable or 'postgres'.
#' @param user The username. Defaults to the 'PGSQL_USER' environment variable or 'postgres'.
#' @param password The database password. Defaults to the 'PGSQL_PASS' environment variable.
#' @param port The port number. Defaults to 5432.
#'
#' @return A DBI connection object.
#' @noRd
pgsql_connect <- function(
  host = Sys.getenv("PGSQL_HOST"),
  dbname = Sys.getenv("PGSQL_DBNAME", "postgres"),
  user = Sys.getenv("PGSQL_USER", "postgres"),
  password = Sys.getenv("PGSQL_PASS"),
  port = 5432
) {
  # Prevent silent local connections
  if (host == "") {
    stop(
      "Database host not found. Please set 'PGSQL_HOST' in your .Renviron file (e.g., using usethis::edit_r_environ())."
    )
  }
  if (password == "") {
    stop(
      "Database password not found. Please set 'PGSQL_PASS' in your .Renviron file."
    )
  }

  # Attempt to establish connection
  con <- tryCatch(
    {
      DBI::dbConnect(
        RPostgres::Postgres(),
        dbname = dbname,
        host = host,
        port = port,
        user = user,
        password = password,
        sslmode = "require"
      )
    },
    error = function(e) {
      stop(paste(
        "Failed to connect to the database. Check your internet connection, VPN, or credentials.\nError:",
        e$message
      ))
    }
  )

  return(con)
}




#' List all tables in the PostgreSQL database
#'
#' This function connects to the PostgreSQL server and retrieves the names of
#' all tables within the connected database.
#'
#' @param con A \code{DBIConnection} object, as returned by \code{\link[DBI]{dbConnect}}.
#' @returns A character vector of table names.
#'
#' @details
#' This is a wrapper around \code{DBI::dbListTables}. For PostgreSQL, this
#' will typically return tables in the 'public' schema unless specified otherwise
#' in the connection.
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#' library(DBI)
#'
#' # Establish connection
#' con <- dbConnect(RPostgres::Postgres(),
#'                  dbname = "sorghum_pangenome_db",
#'                  host = "localhost",
#'                  user = "israeltawiahtetteh")
#'
#' # List tables
#' pgsql_list_tables(con)
#'
#' # Disconnect
#' dbDisconnect(con)
#' }
#'
#'
#' @import DBI
#' @importFrom RPostgres Postgres
#' @noRd
pgsql_list_tables <- function(con) {
  # Mock check for unit testing compatibility with dittodb.
  if (!inherits(con, "DBIMockConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The provided database connection is not valid or has been closed.")
    }
  }

  # Fetch all table names available in the current database schema
  tables <- DBI::dbListTables(con)

  return(tables)
}




#' Get variant statistics from the PostgreSQL database
#'
#' This function calculates summary statistics for variants per chromosome,
#' including variant counts and genomic ranges. It can optionally include
#' statistics for the annotations table.
#'
#' @param con A \code{DBIConnection} object, as returned by \code{\link[DBI]{dbConnect}}.
#' @param include_annotations A logical value indicating whether to include
#' statistics for the annotations table. Defaults to \code{TRUE}.
#'
#' @returns A data frame containing the following columns:
#' \itemize{
#'   \item \code{chrom}: Chromosome name.
#'   \item \code{n_variants}: Total number of variants on the chromosome.
#'   \item \code{min_pos}: The starting genomic position of the first variant.
#'   \item \code{max_pos}: The ending genomic position of the last variant.
#'   \item \code{n_unique_ids}: Number of unique variant IDs.
#'   \item \code{n_annotated}: (Optional) Number of variants with at least one annotation.
#' }
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#' library(DBI)
#'
#' con <- dbConnect(RPostgres::Postgres(),
#'                  dbname = "sorghum_pangenome_db",
#'                  host = "localhost",
#'                  user = "israeltawiahtetteh")
#'
#' # Get chromosome-level statistics
#' stats <- pgsql_variant_stats(con, include_annotations = TRUE)
#' print(stats)
#'
#' dbDisconnect(con)
#' }
#'
#'
#' @import DBI
#' @importFrom RPostgres Postgres
#' @noRd
pgsql_variant_stats <- function(con, include_annotations = TRUE) {
  # Mock check for unit testing compatibility with dittodb.
  if (!inherits(con, "DBIMockConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The provided database connection is not valid or has been closed.")
    }
  }

  # Aggregate basic variant metadata and counts directly on the server for efficiency
  v_query <- "
    SELECT
      chrom,
      COUNT(*) AS n_variants,
      MIN(pos) AS min_pos,
      MAX(pos) AS max_pos,
      COUNT(DISTINCT variant_id) AS n_unique_ids
    FROM variants
    GROUP BY chrom
    ORDER BY chrom
  "

  v_stats <- DBI::dbGetQuery(con, v_query)

  # If requested, cross-reference with the annotations table to see how many variants
  # are functionally characterized.
  if (isTRUE(include_annotations)) {
    a_query <- "
      SELECT chrom, COUNT(DISTINCT v.variant_id) AS n_annotated
      FROM variants v
      JOIN annotations a ON v.variant_id = a.variant_id
      GROUP BY chrom
    "

    a_counts <- DBI::dbGetQuery(con, a_query)

    # Perform a left join to keep all chromosomes found in the variants table
    v_stats <- merge(v_stats, a_counts, by = "chrom", all.x = TRUE)

    # Replace NAs with zero for chromosomes that lack any annotated variants
    v_stats$n_annotated[is.na(v_stats$n_annotated)] <- 0
  }

  return(v_stats)
}




#' Get variant impact summary from the PostgreSQL database
#'
#' This function joins the variants and annotations tables to summarize the
#' distribution of mutation impacts (e.g., HIGH, MODERATE, LOW, MODIFIER)
#' across chromosomes.
#'
#' @param con A \code{DBIConnection} object, as returned by \code{\link[DBI]{dbConnect}}.
#'
#' @returns A data frame in wide format where each row is a chromosome and
#' columns represent the counts for each impact category (e.g., \code{impact_HIGH}).
#'
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#' library(DBI)
#'
#' con <- dbConnect(RPostgres::Postgres(),
#'                  dbname = "sorghum_pangenome_db",
#'                  host = "localhost",
#'                  user = "israeltawiahtetteh")
#'
#' # Generate impact summary
#' impact_stats <- pgsql_variant_impact_summary(con)
#' head(impact_stats)
#'
#' dbDisconnect(con)
#' }
#'
#'
#' @import DBI
#' @importFrom stats reshape
#' @importFrom RPostgres Postgres
#' @noRd
pgsql_variant_impact_summary <- function(con) {
  # Mock check for unit testing compatibility with dittodb.
  if (!inherits(con, "DBIMockConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The provided database connection is not valid or has been closed.")
    }
  }

  # Join annotations with variants to group functional impact categories by chromosome
  query <- "
    SELECT v.chrom, a.impact, COUNT(*) AS n
    FROM annotations a
    JOIN variants v ON a.variant_id = v.variant_id
    GROUP BY v.chrom, a.impact
    ORDER BY v.chrom, a.impact
  "

  impact_df <- DBI::dbGetQuery(con, query)

  # Return an empty data frame if the query returns no records
  if (nrow(impact_df) == 0) {
    return(data.frame())
  }

  # Pivot the data to wide format for easier chromosome-to-chromosome comparison
  impact_wide <- stats::reshape(
    impact_df,
    idvar = "chrom",
    timevar = "impact",
    direction = "wide"
  )

  # Standardize column naming convention and fill missing combinations with zero
  colnames(impact_wide) <- gsub("n\\.", "impact_", colnames(impact_wide))
  impact_wide[is.na(impact_wide)] <- 0

  # Ensure the resulting data frame has clean indices
  rownames(impact_wide) <- NULL

  return(impact_wide)
}




#' Summarize names and row counts for each table in the PostgreSQL database
#'
#' This function iterates through all tables in the connected PostgreSQL database
#' and returns a summary data frame containing the table names and their
#' respective row counts.
#'
#' @param con A \code{DBIConnection} object, as returned by \code{\link[DBI]{dbConnect}}.
#'
#' @returns A data frame with two columns:
#' \itemize{
#'   \item \code{table}: Name of the table.
#'   \item \code{n_rows}: Total number of rows in that table.
#' }
#'
#' @details
#' For very large tables, \code{COUNT(*)} can take
#' a few seconds as PostgreSQL ensures transaction accuracy.
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#' library(DBI)
#'
#' con <- dbConnect(RPostgres::Postgres(),
#'                  dbname = "sorghum_pangenome_db",
#'                  host = "localhost",
#'                  user = "israeltawiahtetteh")
#'
#' # Summarize all tables
#' db_summary <- pgsql_summarize_tables(con)
#' print(db_summary)
#'
#' dbDisconnect(con)
#' }
#'
#'
#' @import DBI
#' @importFrom RPostgres Postgres
#' @noRd
pgsql_summarize_tables <- function(con) {
  # Mock check for unit testing compatibility with dittodb.
  if (!inherits(con, "DBIMockConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The provided database connection is not valid or has been closed.")
    }
  }

  # Pull the names of every table currently defined in the database schema
  tables <- DBI::dbListTables(con)

  # Return an empty structure immediately if no tables are found
  if (length(tables) == 0) {
    return(data.frame(
      table = character(),
      n_rows = integer(),
      stringsAsFactors = FALSE
    ))
  }

  # Looping through each table to get current size
  counts <- sapply(tables, function(tbl) {
    tbl_quoted <- DBI::dbQuoteIdentifier(con, tbl)
    query <- paste0("SELECT COUNT(*) FROM ", tbl_quoted)
    DBI::dbGetQuery(con, query)[[1]]
  })

  # Aggregate results into a clean data frame.
  # We cast row counts to numeric to prevent integer overflow errors with massive pangenome tables.
  res <- data.frame(
    table = tables,
    n_rows = as.numeric(counts),
    stringsAsFactors = FALSE
  )

  # Clear default row names for a polished output
  rownames(res) <- NULL

  return(res)
}




#' Check column names and types for any table in the PostgreSQL database
#'
#' This function retrieves metadata about the columns in a specified table,
#' including column names, data types, and nullability.
#'
#' @param con A \code{DBIConnection} object, as returned by \code{\link[DBI]{dbConnect}}.
#' @param table_name A character value specifying the name of the table.
#' Defaults to "variants", "annotations", "genotypes", or "metadata".
#'
#' @returns A data frame containing column metadata:
#' \itemize{
#'   \item \code{column_name}: The name of the column.
#'   \item \code{data_type}: The PostgreSQL data type (e.g., integer, text, text[]).
#'   \item \code{is_nullable}: Whether the column can contain NULL values.
#' }
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#' library(DBI)
#'
#' con <- dbConnect(RPostgres::Postgres(),
#'                  dbname = "sorghum_pangenome_db",
#'                  host = "localhost",
#'                  user = "israeltawiahtetteh")
#'
#' # Check the genotypes table to verify the new array type
#' pgsql_list_table_columns(con, table_name = "genotypes")
#'
#' dbDisconnect(con)
#' }
#'
#'
#' @import DBI
#' @importFrom RPostgres Postgres
#' @noRd
pgsql_list_table_columns <- function(
  con,
  table_name = c("variants", "annotations", "genotypes", "metadata")
) {
  # Ensure the provided table name matches the supported pangenome schema
  table_name <- match.arg(table_name)

  # Mock check for unit testing compatibility with dittodb.
  if (!inherits(con, "DBIMockConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The provided database connection is not valid or has been closed.")
    }
  }

  # Using interpolation to safely inject the table name into the schema query.
  query <- DBI::sqlInterpolate(
    con,
    "
    SELECT
      column_name,
      data_type,
      is_nullable
    FROM information_schema.columns
    WHERE table_name = ?tbl
    ORDER BY ordinal_position
    ",
    tbl = table_name
  )

  info <- DBI::dbGetQuery(con, query)

  # Alert the user if the query returns nothing, which usually suggests a
  # permissions issue or a missing table in the public schema.
  if (nrow(info) == 0) {
    warning(paste0(
      "Table '",
      table_name,
      "' not found in the information schema. Check your schema permissions."
    ))
  }

  return(info)
}




#' Query PostgreSQL pangenome tables using genomic coordinates
#'
#' This function retrieves data from the variants, annotations, or genotypes
#' tables based on a specific chromosome and genomic range.
#'
#' @param con A \code{DBIConnection} object.
#' @param table_name Character. One of "variants", "annotations", or "genotypes".
#' @param chrom Character. The chromosome name (e.g., "Chr05").
#' @param start,end Numeric. The genomic start and end positions.
#' @param gene_name Character. Optional Sobic ID to filter annotations (e.g., "Sobic.005G213600").
#'
#' @returns A data frame containing the queried genomic data.
#'
#' @details
#' For the 'genotypes' table, the function automatically unpacks the PostgreSQL
#' array column into individual columns named according to the sample library
#' IDs found in the metadata table.
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#' library(DBI)
#'
#' con <- dbConnect(RPostgres::Postgres(), dbname = "sorghum_pangenome_db")
#'
#' # Query a gene region in genotypes
#' gt_data <- pgsql_query_db(con, "genotypes", "Chr05", 75104537, 75106403)
#'
#' # Query annotations for a specific gene ID
#' ann_data <- pgsql_query_db(
#'   con = con,
#'   table_name = "annotations",
#'   chrom = "Chr05",
#'   gene_name = "Sobic.005G213600",
#'   start = 75104537,
#'   end = 75106403
#' )
#'
#' dbDisconnect(con)
#' }
#'
#'
#' @import DBI
#' @importFrom RPostgres Postgres
#' @noRd
pgsql_query_db <- function(
  con,
  table_name = c("variants", "annotations", "genotypes"),
  chrom = NULL,
  start = NULL,
  end = NULL,
  gene_name = NULL
) {
  # Set the target table and verify the connection is active for the session
  table_name <- match.arg(table_name)

  # Mock check for unit testing compatibility with dittodb.
  if (!inherits(con, "DBIMockConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The provided database connection is not valid or has been closed.")
    }
  }

  sql <- ""
  params <- list()

  # Logic for annotations , joining with variants to resolve genomic positions
  if (table_name == "annotations") {
    sql <- "SELECT a.*, v.chrom, v.pos FROM annotations a JOIN variants v ON a.variant_id = v.variant_id WHERE 1=1"

    if (!is.null(chrom)) {
      sql <- paste(sql, "AND v.chrom = ?ch")
      params$ch <- chrom
    }
    if (!is.null(start)) {
      sql <- paste(sql, "AND v.pos >= ?st")
      params$st <- start
    }
    if (!is.null(end)) {
      sql <- paste(sql, "AND v.pos <= ?en")
      params$en <- end
    }
    if (!is.null(gene_name)) {
      sql <- paste(sql, "AND a.gene_name = ?gn")
      params$gn <- gene_name
    }

    sql <- paste(sql, "ORDER BY v.pos")
  } else if (table_name == "genotypes") {
    # Construct query
    sql <- "SELECT g.variant_id, g.chrom, g.pos, v.variant_type, v.ref, v.alt, g.calls FROM genotypes g JOIN variants v ON g.variant_id = v.variant_id WHERE 1=1"

    if (!is.null(chrom)) {
      sql <- paste(sql, "AND g.chrom = ?ch")
      params$ch <- chrom
    }
    if (!is.null(start)) {
      sql <- paste(sql, "AND g.pos >= ?st")
      params$st <- start
    }
    if (!is.null(end)) {
      sql <- paste(sql, "AND g.pos <= ?en")
      params$en <- end
    }

    sql <- paste(sql, "ORDER BY g.pos")
  } else {
    # Default path for retrieving basic variant metadata
    sql <- "SELECT * FROM variants WHERE 1=1"

    if (!is.null(chrom)) {
      sql <- paste(sql, "AND chrom = ?ch")
      params$ch <- chrom
    }
    if (!is.null(start)) {
      sql <- paste(sql, "AND pos >= ?st")
      params$st <- start
    }
    if (!is.null(end)) {
      sql <- paste(sql, "AND pos <= ?en")
      params$en <- end
    }

    sql <- paste(sql, "ORDER BY pos")
  }

  # Execute the dynamically constructed query to prevent interpolation errors
  query <- do.call(DBI::sqlInterpolate, c(list(conn = con, sql = sql), params))
  result <- DBI::dbGetQuery(con, query)

  if (nrow(result) == 0) {
    return(result)
  }

  #  Split calls (text arrays) string into individual sample columns.
  if (table_name == "genotypes") {
    # Strip array braces and split by the comma delimiter
    cleaned_calls <- gsub("[{}]", "", result$calls)
    parsed_list <- strsplit(cleaned_calls, ",")

    # Map the resulting columns to sample names (libraries) stored in the metadata
    sample_metadata <- DBI::dbGetQuery(
      con,
      "SELECT \"lib\" FROM metadata ORDER BY array_index ASC"
    )

    # Reconstruct the matrix to hold all sample calls
    gt_matrix <- do.call(rbind, parsed_list)
    colnames(gt_matrix) <- sample_metadata$lib

    # Combine the metadata columns with the expanded genotype matrix
    info_cols <- result[, setdiff(names(result), "calls"), drop = FALSE]
    result <- cbind(
      info_cols,
      as.data.frame(gt_matrix, stringsAsFactors = FALSE)
    )
  }

  return(result)
}




#' Get the genomic range of a candidate gene using the Sobic ID from a GFF file
#'
#' This function parses a GFF3 file to extract the chromosome, start, and
#' end coordinates for a specific gene ID. It supports local files,
#' GZ-compressed files, and direct URLs.
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




#' Extract variants from PostgreSQL based on mutation impact
#'
#' This function retrieves variants from the 'annotations' table joined with
#' the 'variants' table, filtered by snpEff impact levels (HIGH, MODERATE, etc.)
#' and optional genomic coordinates.
#'
#' @param con A \code{DBIConnection} object, as returned by \code{\link[DBI]{dbConnect}}.
#' @param impact_level Character vector. One or more of "HIGH", "MODERATE",
#'   "LOW", "MODIFIER". Defaults to all four.
#' @param chrom Character. Optional chromosome name (e.g., "Chr05").
#' @param start,end Numeric. Optional genomic start and end positions.
#'
#' @returns A data frame containing variant information and associated functional
#'   annotations.
#'
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#' library(DBI)
#'
#' con <- dbConnect(RPostgres::Postgres(), dbname = "sorghum_pangenome_db")
#'
#' # Get only high-impact variants on Chromosome 5
#'  high_impact_vars <- pgsql_query_by_impact(
#'  con = con,
#'  impact_level = "HIGH",
#'  chrom = 'Chr05',
#'  start = 75104537,
#'  end = 75106403
#' )
#'
#' dbDisconnect(con)
#' }
#'
#'
#' @import DBI
#' @importFrom RPostgres Postgres
#' @noRd
pgsql_query_by_impact <- function(
  con,
  impact_level = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
  chrom = NULL,
  start = NULL,
  end = NULL
) {

  # Mock check for unit testing compatibility with dittodb.
  if (!inherits(con, "DBIMockConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The provided database connection is not valid or has been closed.")
    }
  }

  # Ensure impact levels match the expected SnpEff categories in the database.
  allowed <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  impact_level <- intersect(toupper(impact_level), allowed)

  if (length(impact_level) == 0) {
    stop(
      "Please provide valid impact levels: HIGH, MODERATE, LOW, or MODIFIER."
    )
  }

  # Format the impact levels
  impact_list <- paste(DBI::dbQuoteLiteral(con, impact_level), collapse = ",")

  # Pull variant metadata and the predicted functional effects.
  sql <- paste0(
    "
    SELECT v.*, a.allele, a.annotation, a.impact, a.gene_name,
           a.gene_id, a.feature_type, a.feature_id, a.transcript_biotype,
           a.rank, a.HGVS_c, a.HGVS_p
    FROM variants v
    JOIN annotations a ON v.variant_id = a.variant_id
    WHERE a.impact IN (",
    impact_list,
    ")"
  )

  # Apply genomic coordinate filters if a specific chromosome or window is requested.
  if (!is.null(chrom)) {
    sql <- paste0(sql, " AND v.chrom = ", DBI::dbQuoteLiteral(con, chrom))

    if (!is.null(start) && !is.null(end)) {
      # Use as.numeric to handle large genomic positions without integer overflow.
      sql <- paste0(
        sql,
        " AND v.pos BETWEEN ",
        as.numeric(start),
        " AND ",
        as.numeric(end)
      )
    }
  }

  # Order results by position to maintain genomic orientation.
  sql <- paste0(sql, " ORDER BY v.chrom, v.pos")

  # Execute and return the finalized query results.
  result <- DBI::dbGetQuery(con, sql)

  return(result)
}





#' Compute allele frequencies for a genotype matrix
#'
#' This function calculates the reference and alternate allele frequencies
#' for a variants-by-samples matrix. It handles both phased (|) and unphased (/)
#' VCF genotype strings.
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




#' Filter extracted variants based on alternate allele frequency
#'
#' Calculates allele frequencies for a genotype matrix and filters variants
#' based on a user-defined range locally on your machine. Useful for removing
#' monomorphic or rare variants from pangenome queries.
#'
#' @param gt A data frame or matrix of variants x samples, typically
#'   the output from \code{\link{pg_query_db}}.
#' @param variant_id_col Character. Column name for variant IDs. Default is 'variant_id'.
#' @param chrom_col Character. Optional column name for chromosome.
#' @param pos_col Character. Optional column name for genomic position.
#' @param min_af Numeric. Minimum alternate allele frequency threshold (0-1).
#' @param max_af Numeric. Maximum alternate allele frequency threshold (0-1).
#'
#' @returns A data frame containing variant metadata and calculated frequencies
#'   for variants that passed the filter.
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#'
#' # Query region via the API and pipe into filter to remove rare variants (MAF < 0.05)
#' # Notice: No database connection needed!
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



#' Extract variants based on allele frequencies within a genomic region
#'
#' This function queries the PostgreSQL database for genotypes within a specific
#' genomic range and filters the results to only include variants within the
#' specified alternate allele frequency (AF) thresholds.
#'
#' @param con A \code{DBIConnection} object (PostgreSQL).
#' @param min_af Numeric. Minimum alternate allele frequency (0-1). Default is 0.
#' @param max_af Numeric. Maximum alternate allele frequency (0-1). Default is 1.
#' @param chrom Character. Chromosome name (e.g., "Chr05").
#' @param start,end Numeric. Genomic start and end coordinates.
#'
#' @returns A data frame containing variant metadata (ID, Chrom, Pos) and the
#'   calculated \code{ref_af} and \code{alt_af}.
#'
#' @details
#' The function acts as a high-level wrapper that coordinates data retrieval
#' and frequency filtering. It first pulls raw genotypes using \code{\link{pgsql_query_db}}
#' and then applies the vectorized frequency logic from \code{\link{pgsql_filter_by_af}}.
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#' library(DBI)
#'
#' con <- dbConnect(RPostgres::Postgres(), dbname = "sorghum_pangenome_db")
#'
#' # Get variants on Chr05 with AF between 5% and 50%
#' common_vars <- pgsql_query_by_af(con,
#'                               min_af = 0.05,
#'                               max_af = 0.50,
#'                               chrom = "Chr05",
#'                               start = 1000000,
#'                               end = 2000000)
#'
#' dbDisconnect(con)
#' }
#'
#'
#' @import DBI
#' @importFrom RPostgres Postgres
#' @noRd
pgsql_query_by_af <- function(
  con,
  min_af = 0,
  max_af = 1,
  chrom = NULL,
  start = NULL,
  end = NULL
) {
  # Mock check for unit testing compatibility with dittodb.
  if (!inherits(con, "DBIMockConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The provided database connection is not valid or has been closed.")
    }
  }

  # Chromosome must be specified. Pulling genotypes for the entire
  # pangenome (1,676 samples) would exceed available RAM in most environments.
  if (is.null(chrom)) {
    stop(
      "Chromosome ('chrom') must be specified for an AF query to prevent memory overflow."
    )
  }

  # Retrieve the expanded genotype matrix for the requested genomic window.
  # This returns the individual sample columns required for frequency math.
  gt_data <- pgsql_query_db(
    con = con,
    table_name = "genotypes",
    chrom = chrom,
    start = start,
    end = end
  )

  # If the database returns an empty set for this region, we exit early
  if (nrow(gt_data) == 0) {
    message("No variants found in the specified genomic region.")
    return(data.frame())
  }

  # Apply the allele frequency filter using our optimized calculation engine.
  # This converts the matrix to dosages and subsets based on the AF thresholds.
  result <- pg_filter_by_af(
    gt = gt_data,
    variant_id_col = 'variant_id',
    chrom_col = 'chrom',
    pos_col = 'pos',
    min_af = min_af,
    max_af = max_af
  )

  return(result)
}




#' Query genotypes for specific variant IDs from PostgreSQL
#'
#' This function retrieves genomic data for a specific list of variant IDs.
#' It joins the 'variants' metadata with the 'genotypes' allele calls and
#' expands the genotype array into a wide format (samples as columns).
#'
#' @param con A \code{DBIConnection} object (PostgreSQL).
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
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#' library(DBI)
#'
#' con <- dbConnect(RPostgres::Postgres(), dbname = "sorghum_pangenome_db")
#'
#' # Query specific variants
#' my_ids <- c("SNP_Chr05_75104557", "INDEL_Chr05_75104541")
#' df <- pgsql_query_genotypes(con, variant_ids = my_ids)
#'
#' dbDisconnect(con)
#' }
#'
#'
#' @import DBI
#' @importFrom RPostgres Postgres
#' @noRd
pgsql_query_genotypes <- function(
  con,
  variant_ids,
  variant_id_col = "variant_id",
  variants_table = 'variants',
  genotypes_table = 'genotypes',
  meta_data = NULL
) {

  if (!inherits(con, "DBIMockConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The provided database connection is not valid or has been closed.")
    }
  }

  # Return a warning if the input vector is empty
  if (length(variant_ids) == 0) {
    warning("The 'variant_ids' vector is empty. No query was executed.")
    return(data.frame())
  }

  if (is.null(meta_data)) {
    meta_data <- DBI::dbListFields(con, variants_table)
  } else {
    meta_data <- unique(c(variant_id_col, meta_data))
  }

  v_tbl <- DBI::dbQuoteIdentifier(con, variants_table)
  g_tbl <- DBI::dbQuoteIdentifier(con, genotypes_table)
  id_col <- DBI::dbQuoteIdentifier(con, variant_id_col)
  quoted_ids <- paste(DBI::dbQuoteLiteral(con, variant_ids), collapse = ",")
  meta_cols_sql <- paste0("v.", DBI::dbQuoteIdentifier(con, meta_data), collapse = ", ")

  sql <- paste0(
    "SELECT ", meta_cols_sql, ", g.calls ",
    "FROM ", v_tbl, " v ",
    "JOIN ", g_tbl, " g ON v.", id_col, " = g.", id_col, " ",
    "WHERE v.", id_col, " IN (", quoted_ids, ") ",
    "ORDER BY v.chrom, v.pos"
  )

  res <- DBI::dbGetQuery(con, sql)

  if (nrow(res) == 0) {
    warning("No data found for the provided variant IDs.")
    return(data.frame())
  }

  cleaned_calls <- gsub("[{}]", "", res$calls)
  parsed_list <- strsplit(cleaned_calls, ",")

  sample_metadata <- DBI::dbGetQuery(
    con,
    "SELECT \"lib\" FROM metadata ORDER BY array_index ASC"
  )

  gt_matrix <- do.call(rbind, parsed_list)
  colnames(gt_matrix) <- sample_metadata$lib

  info_cols <- res[, setdiff(names(res), "calls"), drop = FALSE]
  result <- cbind(
    info_cols,
    as.data.frame(gt_matrix, stringsAsFactors = FALSE)
  )

  return(result)
}




#' Count the distribution of variant types in the PostgreSQL database
#'
#' This function performs a server-side aggregation to count the occurrences
#' of different variant types (e.g., SNP, INDEL) stored in the 'variants' table.
#'
#' @param con A \code{DBIConnection} object (PostgreSQL).
#' @param variants_table Character. The name of the table containing variant
#'   metadata. Defaults to "variants".
#'
#' @returns A data frame with two columns:
#' \itemize{
#'   \item \code{variant_type}: The category of the variant.
#'   \item \code{n}: The total count of variants for that category.
#' }
#'
#' @details
#' This function leverages PostgreSQL's \code{COUNT} and \code{GROUP BY}
#' operations. Given the size of the pangenome, this is significantly faster
#' than retrieving the data into R for counting.
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#' library(DBI)
#'
#' con <- dbConnect(RPostgres::Postgres(), dbname = "sorghum_pangenome_db")
#'
#' # Get the SNP vs INDEL breakdown
#' type_counts <- pgsql_count_variant_types(con)
#' print(type_counts)
#'
#' dbDisconnect(con)
#' }
#'
#'
#' @import DBI
#' @importFrom RPostgres Postgres
#' @noRd
pgsql_count_variant_types <- function(con, variants_table = "variants") {
  # Mock check for unit testing compatibility with dittodb.
  if (!inherits(con, "DBIMockConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The provided database connection is not valid or has been closed.")
    }
  }

  # Confirm the table exists in the schema before execution to avoid generic SQL errors.
  if (!DBI::dbExistsTable(con, variants_table)) {
    stop(sprintf("Table '%s' does not exist in the database.", variants_table))
  }

  # Verify the presence of the type column to ensure the GROUP BY operation is valid.
  cols <- DBI::dbListFields(con, variants_table)
  if (!"variant_type" %in% cols) {
    stop(sprintf(
      "The table '%s' does not have a 'variant_type' column.",
      variants_table
    ))
  }

  # Use quoted identifiers to maintain compatibility with case-sensitive PostgreSQL environments.
  tbl_quoted <- DBI::dbQuoteIdentifier(con, variants_table)

  # Grouping on the server side is vital for pangenome scale to prevent R from crashing.
  query <- paste0(
    "SELECT variant_type, COUNT(*) AS n ",
    "FROM ",
    tbl_quoted,
    " ",
    "GROUP BY variant_type ",
    "ORDER BY variant_type"
  )

  result <- DBI::dbGetQuery(con, query)

  # Coerce counts to numeric to prevent bit64 or integer overflow issues with large variant sets.
  if (nrow(result) > 0) {
    result$n <- as.numeric(result$n)
  }

  return(result)
}




#' Summarize genomic annotations and impacts in a specific region
#'
#' This function queries the PostgreSQL database for variants within a specific
#' genomic range and returns summaries of SnpEff annotations and impact levels,
#' cross-tabulated by variant type (e.g., SNP, INDEL).
#'
#' @param con A \code{DBIConnection} object (PostgreSQL).
#' @param chrom Character. Chromosome name (e.g., "Chr05").
#' @param start Numeric. Start coordinate of the region.
#' @param end Numeric. End coordinate of the region.
#' @param annotations_table Character. Name of the annotations table. Defaults to "annotations".
#' @param variants_table Character. Name of the variants table. Defaults to "variants".
#'
#' @returns A list containing three data frames:
#' \itemize{
#'   \item \code{annotation_summary}: Counts of specific SnpEff annotations per variant type.
#'   \item \code{impact_summary}: Counts of impact categories per variant type.
#'   \item \code{variant_type_totals}: Total number of SNPs and INDELs in the region.
#' }
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#' library(DBI)
#'
#' con <- dbConnect(RPostgres::Postgres(), dbname = "sorghum_pangenome_db")
#'
#' # Summarize a 2kb region on Chromosome 5
#' region_summary <- pgsql_query_ann_summary(con,
#'                                        chrom = "Chr05",
#'                                        start = 75104537,
#'                                        end = 75106403)
#'
#' # View the impact distribution
#' print(region_summary$impact_summary)
#'
#' dbDisconnect(con)
#' }
#'
#'
#' @import DBI
#' @importFrom RPostgres Postgres
#' @noRd
pgsql_query_ann_summary <- function(
  con,
  chrom = NULL,
  start = NULL,
  end = NULL,
  annotations_table = "annotations",
  variants_table = "variants"
) {
  # Mock check for unit testing compatibility with dittodb.
  if (!inherits(con, "DBIMockConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The provided database connection is not valid or has been closed.")
    }
  }

  # Defensive check to ensure the required pangenome tables exist in the schema.
  if (!DBI::dbExistsTable(con, annotations_table)) {
    stop(sprintf("Table '%s' not found.", annotations_table))
  }
  if (!DBI::dbExistsTable(con, variants_table)) {
    stop(sprintf("Table '%s' not found.", variants_table))
  }

  # Quote identifiers for safe SQL construction.
  v_tbl <- DBI::dbQuoteIdentifier(con, variants_table)
  a_tbl <- DBI::dbQuoteIdentifier(con, annotations_table)

  # Calculate the baseline variant distribution for the specified region.
  total_sql_str <- paste0(
    "SELECT variant_type, COUNT(*) AS total_variants FROM ",
    v_tbl,
    " WHERE 1=1"
  )
  params_base <- list()
  if (!is.null(chrom)) {
    total_sql_str <- paste(total_sql_str, "AND chrom = ?ch")
    params_base$ch <- chrom
  }
  if (!is.null(start)) {
    total_sql_str <- paste(total_sql_str, "AND pos >= ?st")
    params_base$st <- start
  }
  if (!is.null(end)) {
    total_sql_str <- paste(total_sql_str, "AND pos <= ?en")
    params_base$en <- end
  }
  total_sql_str <- paste(total_sql_str, "GROUP BY variant_type")

  total_sql <- do.call(
    DBI::sqlInterpolate,
    c(list(conn = con, sql = total_sql_str), params_base)
  )
  variant_type_totals <- DBI::dbGetQuery(con, total_sql)

  # Cross-tabulate specific SnpEff annotations by variant type.
  ann_sql_str <- paste0(
    "SELECT a.annotation, v.variant_type, COUNT(*) AS count FROM ",
    a_tbl,
    " a JOIN ",
    v_tbl,
    " v ON a.variant_id = v.variant_id WHERE 1=1"
  )
  params_v <- list()
  if (!is.null(chrom)) {
    ann_sql_str <- paste(ann_sql_str, "AND v.chrom = ?ch")
    params_v$ch <- chrom
  }
  if (!is.null(start)) {
    ann_sql_str <- paste(ann_sql_str, "AND v.pos >= ?st")
    params_v$st <- start
  }
  if (!is.null(end)) {
    ann_sql_str <- paste(ann_sql_str, "AND v.pos <= ?en")
    params_v$en <- end
  }
  ann_sql_str <- paste(
    ann_sql_str,
    "GROUP BY a.annotation, v.variant_type ORDER BY count DESC"
  )

  ann_sql <- do.call(
    DBI::sqlInterpolate,
    c(list(conn = con, sql = ann_sql_str), params_v)
  )
  annotation_summary <- DBI::dbGetQuery(con, ann_sql)

  # Aggregate functional impact levels to help identify high-priority regions.
  imp_sql_str <- paste0(
    "SELECT a.impact, v.variant_type, COUNT(*) AS count FROM ",
    a_tbl,
    " a JOIN ",
    v_tbl,
    " v ON a.variant_id = v.variant_id WHERE 1=1"
  )
  if (!is.null(chrom)) {
    imp_sql_str <- paste(imp_sql_str, "AND v.chrom = ?ch")
  }
  if (!is.null(start)) {
    imp_sql_str <- paste(imp_sql_str, "AND v.pos >= ?st")
  }
  if (!is.null(end)) {
    imp_sql_str <- paste(imp_sql_str, "AND v.pos <= ?en")
  }
  imp_sql_str <- paste(
    imp_sql_str,
    "GROUP BY a.impact, v.variant_type ORDER BY count DESC"
  )

  imp_sql <- do.call(
    DBI::sqlInterpolate,
    c(list(conn = con, sql = imp_sql_str), params_v)
  )
  impact_summary <- DBI::dbGetQuery(con, imp_sql)

  # Consolidate summaries into a single list for downstream reporting or plotting.
  return(list(
    annotation_summary = annotation_summary,
    impact_summary = impact_summary,
    variant_type_totals = variant_type_totals
  ))
}




#' Retrieve sample metadata from the PostgreSQL database
#'
#' This function fetches accession-level metadata, such as origin, race,
#' and classification, from the 'metadata' table. It supports optional
#' filtering by specific columns to subset populations for analysis.
#'
#' @param con A \code{DBIConnection} object (PostgreSQL).
#' @param query_col Character. The metadata column to filter by (e.g., "countryorigin").
#'   If \code{NULL}, all records are returned.
#' @param query_value Character or Numeric. The specific value to match in \code{query_col}.
#'
#' @returns A data frame containing the sample metadata records.
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#' library(DBI)
#'
#' con <- dbConnect(RPostgres::Postgres(), dbname = "sorghum_pangenome_db")
#'
#' # Fetch metadata for all accessions originating from Ghana
#' ghana_samples <- pgsql_get_sample_metadata(con,
#'                                         query_col = "countryorigin",
#'                                         query_value = "Ghana")
#'
#' dbDisconnect(con)
#' }
#'
#'
#' @import DBI
#' @importFrom RPostgres Postgres
#' @noRd
pgsql_get_sample_metadata <- function(con, query_col = NULL, query_value = NULL) {
    # Mock check for unit testing compatibility with dittodb.
  if (!inherits(con, "DBIMockConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The provided database connection is not valid or has been closed.")
    }
  }

  # Validate the existence of the filter column to prevent blind SQL failures
  if (!is.null(query_col)) {
    available_cols <- DBI::dbListFields(con, "metadata")
    if (!query_col %in% available_cols) {
      stop(sprintf("Column '%s' not found in the metadata table.", query_col))
    }
  }

  sql <- "SELECT * FROM metadata"

  #  WHERE clause query
  if (!is.null(query_col) && !is.null(query_value)) {
    sql <- paste0(
      sql,
      " WHERE ",
      DBI::dbQuoteIdentifier(con, query_col),
      " = ",
      DBI::dbQuoteLiteral(con, query_value)
    )
  }

 # sorting by array index to ensure metadata aligns with the genotype array order
  sql <- paste0(sql, " ORDER BY array_index ASC")

  return(DBI::dbGetQuery(con, sql))
}




#' Query genotypes filtered by sample metadata attributes
#'
#' This function retrieves genotypes for a genomic region but only for
#' a subset of samples defined by metadata attributes (e.g., specific
#' countries, populations, or clusters).
#'
#' @param con A \code{DBIConnection} object.
#' @param chrom Character. Chromosome name.
#' @param start,end Numeric. Genomic range.
#' @param meta_col Character. Metadata column to filter samples by (e.g., "countryorigin").
#' @param meta_value Character. Value to match in \code{meta_col} (e.g., "Ghana").
#'
#' @returns A wide-format data frame (variants x filtered samples).
#'
#' @examples
#' \dontrun{
#' # Get genotypes for a gene, but only for samples from Ethiopia
#' eth_genotypes <- pgsql_query_by_metadata(con,
#'                                       chrom = "Chr05",
#'                                       start = 75104537, end = 75106403,
#'                                       meta_col = "countryorigin",
#'                                       meta_value = "Ethiopia")
#' }
#'
#' @noRd
pgsql_query_by_metadata <- function(con, chrom, start, end, meta_col, meta_value) {
  # Get samples that match the metadata criteria using 'lib' names and 'array_index'
  sample_info <- pgsql_get_sample_metadata(con, meta_col, meta_value)

  if (nrow(sample_info) == 0) {
    stop(
      "No samples found matching the criteria: ",
      meta_col,
      " = ",
      meta_value
    )
  }

  indices <- sample_info$array_index
  sample_names <- sample_info$lib

  # Fetch the raw data for the region
  raw_data <- pgsql_query_db(
    con,
    table_name = "genotypes",
    chrom = chrom,
    start = start,
    end = end
  )

  if (nrow(raw_data) == 0) {
    return(data.frame())
  }

  result <- raw_data[, c("variant_id", "chrom", "pos", sample_names)]

  return(result)
}



#' Interactive Geographic Exploration of Sorghum Accessions
#'
#' This function generates a high-performance interactive map showing the geographic
#' distribution of sorghum lines. It dynamically generates rich, scrollable popups
#' containing all available metadata for each accession.
#'
#' @param metadata A data frame containing sample metadata. Must include 'lat' and 'lon'
#'   columns, along with the specified coloring column.
#' @param color_by Character. The metadata column to use for point coloration.
#'   Defaults to "countryorigin".
#'
#' @returns A \code{leaflet} map object (htmlwidget) representing the interactive map.
#'
#' @details
#' The function automatically filters out records with missing geographic coordinates.
#' Instead of hardcoding specific tooltip values, it dynamically reads all non-coordinate
#' columns from the provided \code{metadata} data frame and formats them into a
#' scrollable HTML popup for each point.
#'
#' \strong{Dependency Note:} To keep the core package lightweight, the \code{leaflet}
#' and \code{tools} packages are listed as "Suggested" dependencies. If they are not
#' currently installed on your system, the function will gracefully stop and prompt
#' you to install them before plotting.
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#'
#' # Fetch sample metadata from the database
#' meta <- pg_get_sample_metadata()
#'
#' # Explore the geographic distribution colored by genetic cluster
#' pg_map_accessions(meta, color_by = "kmeans_cluster")
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

  # Check if coordinate columns exist and filter for complete cases
  if (!all(c("lat", "lon") %in% names(metadata))) {
    stop("Metadata must contain 'lat' and 'lon' columns for geographic mapping.")
  }

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
      color = ~pal(plot_data[[color_by]]),
      stroke = FALSE,
      fillOpacity = 0.8,
      popup = popup_info,
      label = ~lib
    )

  return(map)
}


#' Generic Trait Association Audit (Cloud/API Version)
#'
#' This function performs an association analysis between a specific genomic
#' variant and a phenotypic trait. It pulls genotypes from the AWS API and
#' merges them with your local phenotypic data for statistical testing and plotting.
#'
#' @param variant_id Character. The specific SNP or INDEL ID to analyze.
#' @param pheno_df Data frame containing local phenotypic data. The first column
#'   must contain accession identifiers (PI numbers).
#' @param trait_col Character. The name of the column in \code{pheno_df}
#'   containing the trait scores.
#' @param trait_type Character. Either "qualitative" or "quantitative".
#'   Determines the statistical test and plot type.
#'
#' @returns A ggplot object displaying the association and statistical summary.
#'
#' @export
#' @import ggplot2
#' @import stats
pg_plot_trait_association <- function(
  variant_id,
  pheno_df,
  trait_col,
  trait_type = c("qualitative", "quantitative")
) {
  trait_type <- match.arg(trait_type)
  genotype <- NULL # Resolve R CMD check for ggplot2 global variable

  # Dependency check
  if (trait_type == "qualitative" && !requireNamespace("scales", quietly = TRUE)) {
    stop(
      "The 'scales' package is required to format qualitative trait plots. ",
      "Please install it by running: install.packages('scales')",
      call. = FALSE
    )
  }

  # 1. Pull genotype data using the API wrapper (Cloud fetch)
  gt_data <- pg_query_genotypes(variant_ids = variant_id)
  
  if (nrow(gt_data) == 0) {
    stop(paste(
      "Variant ID",
      variant_id,
      "was not found in the cloud database."
    ))
  }

  # 2. Fetch sample metadata mapping from the cloud
  pi_map <- pg_get_sample_metadata()

  # Note: To avoid creating a brand new API endpoint just for a single annotation impact, 
  # we default it to "UNKNOWN" for the API wrapper. If you want the full verdict logic, 
  # you can add a dedicated annotation endpoint later.
  impact <- "UNKNOWN"

  # Align genotype samples with phenotype PI numbers
  sample_cols <- setdiff(
    names(gt_data),
    c(
      "variant_id",
      "chrom",
      "pos",
      "ref",
      "alt",
      "qual",
      "filter",
      "variant_type"
    )
  )

  gt_long <- data.frame(
    pinumber = pi_map$pinumber[match(sample_cols, pi_map$lib)],
    genotype = as.character(unlist(gt_data[1, sample_cols])),
    stringsAsFactors = FALSE
  )

  # Ensure the phenotype data frame is ready for merging
  colnames(pheno_df)[1] <- "pinumber"
  plot_data <- merge(gt_long, pheno_df, by = "pinumber")

  # Filter out missing genotypes and phenotypes to ensure statistical integrity
  plot_data <- plot_data[
    !is.na(plot_data[[trait_col]]) &
      plot_data[[trait_col]] != "-" &
      plot_data$genotype != "./.",
  ]

  # Clean trait data based on user-defined type.
  if (trait_type == "quantitative") {
    if (!is.numeric(plot_data[[trait_col]])) {
      plot_data[[trait_col]] <- as.numeric(gsub(
        ".*[^0-9.]",
        "",
        as.character(plot_data[[trait_col]])
      ))
    }
  } else {
    # Categorical traits are forced to factors to prevent numerical misinterpretation
    plot_data[[trait_col]] <- as.factor(plot_data[[trait_col]])
  }

  # Compute population genetics parameters (MAF) for the audit
  alleles <- unlist(strsplit(plot_data$genotype, "[|/]"))
  maf <- min(table(alleles)) / length(alleles)
  p_val <- NA
  pve_val <- 0

  # Statistical tests based on the trait distribution
  if (trait_type == "quantitative") {
    p_val <- stats::kruskal.test(
      plot_data[[trait_col]] ~ plot_data$genotype
    )$p.value
    pve_val <- summary(stats::lm(
      plot_data[[trait_col]] ~ plot_data$genotype
    ))$r.squared
  } else {
    p_val <- stats::chisq.test(
      table(plot_data$genotype, plot_data[[trait_col]]),
      simulate.p.value = TRUE
    )$p.value
  }

  is_sig <- !is.na(p_val) && p_val < 0.05

  # The VALIDATE verdict is triggered by strong stats (PVE > 10%)
  verdict <- if (is_sig && pve_val > 0.1) {
    "VALIDATE"
  } else if (is_sig) {
    "CAUTION"
  } else {
    "DISCARD"
  }

  # Generate the diagnostic visualization (Local processing)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = genotype)) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette = "Set1")

  if (trait_type == "quantitative") {
    p <- p +
      ggplot2::aes(y = .data[[trait_col]], fill = genotype) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      ggplot2::geom_jitter(width = 0.2, alpha = 0.4)
  } else {
    p <- p +
      ggplot2::aes(fill = .data[[trait_col]]) +
      ggplot2::geom_bar(position = "fill", color = "white", linewidth = 0.3) +
      ggplot2::scale_y_continuous(labels = scales::percent)
  }

  # Finalize plot metadata and labels
  p <- p +
    ggplot2::labs(
      title = paste("Association Audit:", variant_id),
      subtitle = paste0(
        "P-val: ",
        format.pval(p_val, digits = 3),
        " | PVE: ",
        round(pve_val * 100, 1),
        "% | MAF: ",
        round(maf, 3),
        "\nVerdict: ",
        verdict
      ),
      x = "Genotype",
      y = if (trait_type == "quantitative") trait_col else "Frequency (%)"
    )

  return(p)
}