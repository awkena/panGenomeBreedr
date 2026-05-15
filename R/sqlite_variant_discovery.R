
#' List all tables in the SQLite database.
#' @param db_path A character value indicating the path to the SQLite database.
#' @returns A character vector of names of SQLite databases
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # List tables in the SQLite database
#' list_sqlite_tables(mini_db_path)
#'
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                              full.names = TRUE,
#'                              recursive = TRUE,
#'                              all.files = TRUE,
#'                              include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#' }
#'
#' @export
#' @import DBI
#' @importFrom RSQLite SQLite
list_sqlite_tables <- function(db_path) {

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  # assign("con", con, envir = .GlobalEnv)

  tables <- DBI::dbListTables(con)
  DBI::dbDisconnect(con)

  return(tables)

}

#' Get variants statistics stored in SQLite database
#' @param db_path A character value indicating the path to the SQLite database.
#' @param include_annotations A logical value indicating whether to include
#' statistics for the annotations table.
#' @returns A data frame of variant statistics in the SQLite databases
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # Get variant statistics in the SQLite database
#' variant_stats(mini_db_path)
#'
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                              full.names = TRUE,
#'                              recursive = TRUE,
#'                              all.files = TRUE,
#'                              include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#' }
#'
#' @export
variant_stats <- function(db_path, include_annotations = TRUE) {

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  # assign("con", con, envir = .GlobalEnv)

  # Base stats from variants table
  variant_stats <- DBI::dbGetQuery(con, "
    SELECT
      chrom,
      COUNT(*) AS n_variants,
      MIN(pos) AS min_pos,
      MAX(pos) AS max_pos,
      COUNT(DISTINCT variant_id) AS n_unique_ids
    FROM variants
    GROUP BY chrom
    ORDER BY chrom
  ")

  if (include_annotations) {
    annotation_counts <- DBI::dbGetQuery(con, "
      SELECT chrom, COUNT(DISTINCT v.variant_id) AS n_annotated
      FROM variants v
      JOIN annotations a ON v.variant_id = a.variant_id
      GROUP BY chrom
    ")
    variant_stats <- merge(variant_stats, annotation_counts, by = "chrom", all.x = TRUE)
    variant_stats$n_annotated[is.na(variant_stats$n_annotated)] <- 0
  }

  DBI::dbDisconnect(con)

  return(variant_stats)

}

#' Get variants statistics stored in SQLite database based on mutation impact.
#' @param db_path A character value indicating the path to the SQLite database.
#' @returns A data frame of variant impact summary or counts.
#'
#' @examples
#'
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # Get variant impact summary in the SQLite database
#' variant_impact_summary(mini_db_path)
#'
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                              full.names = TRUE,
#'                              recursive = TRUE,
#'                              all.files = TRUE,
#'                              include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#' }
#'
#' @export
variant_impact_summary <- function(db_path) {

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  # assign("con", con, envir = .GlobalEnv)

  query <- "
    SELECT v.chrom, a.impact, COUNT(*) AS n
    FROM annotations a
    JOIN variants v ON a.variant_id = v.variant_id
    GROUP BY v.chrom, a.impact
    ORDER BY v.chrom, a.impact
  "

  impact_df <- DBI::dbGetQuery(con, query)
  DBI::dbDisconnect(con)

  # Spread to wide format (optional)
  impact_wide <- stats::reshape(impact_df,
                                idvar = "chrom",
                                timevar = "impact",
                                direction = "wide")

  # Clean column names
  colnames(impact_wide) <- gsub("n\\.", "impact_", colnames(impact_wide))
  impact_wide[is.na(impact_wide)] <- 0
  return(impact_wide)
}

#' Name and row count for each table in SQLite database.
#' @param db_path A character value indicating the path to the SQLite database.
#' @returns A data frame row counts of tables in SQLite database.
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # Row counts for each table in the database
#' summarize_sqlite_tables(mini_db_path)
#'
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                              full.names = TRUE,
#'                              recursive = TRUE,
#'                              all.files = TRUE,
#'                              include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#' }
#'
#' @export
summarize_sqlite_tables <- function(db_path) {

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  # assign("con", con, envir = .GlobalEnv)

  tables <- DBI::dbListTables(con)
  counts <- sapply(tables, function(tbl) {
    DBI::dbGetQuery(con, paste0("SELECT COUNT(*) FROM ", tbl))[[1]]
  })

  DBI::dbDisconnect(con)
  data.frame(table = tables, n_rows = counts)

}

#' Get the genomic range of a candidate gene using the Sobic ID from a gff file.
#' @param gene_name A character value indicating the Sobic ID of candidate gene.
#' @param gff_path A character value specifying the path to gff3 file. GZ compressed
#' files and URL file paths are acceptable.
#' @returns A list object of three components indicating the chromosome, start
#' and end coordinates of candidate gene.
#' @examples
#' \donttest{
#' # example code
#'
#' # Path to GFF3 file
#' gff_path <- "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz"
#' gene_coord_gff(gene_name = "Sobic.005G213600",
#'                gff_path = gff_path)
#'
#' }
#' @export
#' @importFrom R.utils isUrl isGzipped gunzip
gene_coord_gff <- function(gene_name,
                           gff_path) {

  # Read GFF and extract gene coordinates
  if (R.utils::isUrl(gff_path)) {

    utils::download.file(gff_path,
                         destfile = file.path(tempdir(), basename(gff_path)))


    if (R.utils::isGzipped(file.path(tempdir(), basename(gff_path)))) {

      R.utils::gunzip(file.path(tempdir(), basename(gff_path)),
                      destname = file.path(tempdir(), 'my_gff.gff3'),
                      overwrite = TRUE,
                      remove = TRUE)

      gff <- utils::read.delim(file.path(tempdir(), "my_gff.gff3"),
                               comment.char = "#",
                               header = FALSE,
                               sep = "\t")
    } else {

      gff <- utils::read.delim(file.path(tempdir(), basename(gff_path)),
                               comment.char = "#",
                               header = FALSE,
                               sep = "\t")

    }

    # Ensure tempdir contents are deleted on exit
    on.exit({
      contents <- list.files(tempdir(),
                             full.names = TRUE,
                             recursive = TRUE,
                             all.files = TRUE,
                             include.dirs = TRUE)
      unlink(contents, recursive = TRUE, force = TRUE)
    }, add = TRUE)

  } else if (R.utils::isGzipped(gff_path)) {

    R.utils::gunzip(gff_path,
                    destname = file.path(dirname(gff_path), 'my_gff.gff3'),
                    overwrite = TRUE,
                    remove = TRUE)

    gff <- utils::read.delim(file.path(dirname(gff_path), 'my_gff.gff3'),
                             comment.char = "#",
                             header = FALSE,
                             sep = "\t")
  } else {

    gff <- utils::read.delim(gff_path,
                             comment.char = "#",
                             header = FALSE,
                             sep = "\t")
  }

  colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand",
                     "phase", "attributes")
  genes <- gff[gff$type == "gene", ]
  genes$gene_name <- sub(".*Name=([^;]+).*", "\\1", genes$attributes)

  gmatch <- genes[genes$gene_name == gene_name, ]
  if (nrow(gmatch) == 0) stop("Gene not found in GFF.")
  if (nrow(gmatch) > 1) warning("Multiple matches found; using the first match.")

  list(chrom = gmatch$seqid[1],
       start = as.integer(gmatch$start[1]),
       end   = as.integer(gmatch$end[1]))
}

#' Check the column names and types for any table in a SQLite database.
#' @param db_path A character value indicating the path to the SQLite database.
#' @param table_name A character value specifying the name of any table in a
#' SQLite database.
#' @returns A data frame of table column information in the SQLite database.
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # List column names and type of variable in the variants table in the database
#' list_table_columns(mini_db_path, table_name = "variants")
#'
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                              full.names = TRUE,
#'                              recursive = TRUE,
#'                              all.files = TRUE,
#'                              include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#' }
#'
#' @export
list_table_columns <- function(db_path,
                               table_name = c("variants", "annotations", "genotypes")) {

  table_name <- match.arg(table_name)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  # assign("con", con, envir = .GlobalEnv)

  info <- DBI::dbGetQuery(con, paste0("PRAGMA table_info(", table_name, ")"))
  DBI::dbDisconnect(con)

  return(info)

}


#' Query any table in your SQLite database using chromosome and a genomic position range.
#'
#' @param db_path A character value indicating the path to the SQLite database.
#' @param table_name A character value specifying the name of any table in a
#' SQLite database.
#' @param chrom A character value specifying the chromosome name.
#' @param start,end A numeric value specifying the start and end coordinates for
#' the candidate gene.
#' @param gene_name A character value indicating the Sobic ID of candidate gene,
#' if available. It used if querying the `annotations` table in the database.
#' @returns A data frame of the table queried in the SQLite database.
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # Extract snpEff annotation for variants within the candidate gene Sobic.005G213600
#' annota_region <- query_db(db_path = mini_db_path,
#'                            chrom = "Chr05",
#'                            start = 75104537,
#'                            end = 75106403,
#'                            table_name = "annotations",
#'                            gene_name = "Sobic.005G213600")
#'
#' # Extract VCF genotypes for variants within the region: 'Chr05:75104537-75106403'
#' gt_region <- query_db(db_path = mini_db_path,
#'                       chrom = "Chr05",
#'                       start = 75104537,
#'                       end = 75106403,
#'                       table_name = "genotypes")
#'
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                              full.names = TRUE,
#'                              recursive = TRUE,
#'                              all.files = TRUE,
#'                              include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#'
#' }
#'
#' @export
#' @import DBI
#' @importFrom RSQLite SQLite
query_db <- function(
  db_path,
  table_name = c("variants", "annotations", "genotypes"),
  chrom,
  start = NULL,
  end = NULL,
  gene_name = NULL
) {
  table_name <- match.arg(table_name)

  # Connect to the SQLite database
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Ensure the connection is closed when the function exits, even on failure
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Check if table exists
  all_tables <- DBI::dbListTables(con)
  if (!(table_name %in% all_tables)) {
    stop(sprintf("Table '%s' does not exist in the database.", table_name))
  }

  # ---  Annotations ---
  if (table_name == "annotations") {
    sql <- "
      SELECT a.*, v.chrom, v.pos
      FROM annotations a
      JOIN variants v ON a.variant_id = v.variant_id
      WHERE v.chrom = :chrom"

    params <- list(chrom = chrom)

    if (!is.null(start) && !is.null(end)) {
      sql <- paste(sql, "AND v.pos BETWEEN :start AND :end")
      params$start <- start
      params$end <- end
    }

    # Filter by gene_name directly in SQL
    if (!is.null(gene_name)) {
      sql <- paste(sql, "AND a.gene_name = :gene_name")
      params$gene_name <- gene_name
    }

    sql <- paste(sql, "ORDER BY v.pos")
    result <- DBI::dbGetQuery(con, sql, params = params)

    return(result)
  }

  # --- Genotypes ---
  if (table_name == "genotypes") {
    # Dynamically pull columns to support both wide and packed 'calls' formats
    gt_cols <- DBI::dbListFields(con, "genotypes")
    gt_data_cols <- setdiff(gt_cols, c("variant_id", "chrom", "pos"))
    gt_cols_sql <- paste0("g.", gt_data_cols, collapse = ", ")

    sql <- sprintf(
      "
      SELECT g.variant_id, g.chrom, g.pos, v.variant_type, v.ref, v.alt, %s
      FROM genotypes g
      JOIN variants v ON g.variant_id = v.variant_id
      WHERE g.chrom = :chrom",
      gt_cols_sql
    )

    params <- list(chrom = chrom)

    if (!is.null(start) && !is.null(end)) {
      sql <- paste(sql, "AND g.pos BETWEEN :start AND :end")
      params$start <- start
      params$end <- end
    }

    sql <- paste(sql, "ORDER BY g.pos")
    result <- DBI::dbGetQuery(con, sql, params = params)

    if (nrow(result) == 0) {
      return(result)
    }

    # Unpack the array if the database uses the packed 4-column 'calls' structure
    if ("calls" %in% gt_data_cols && length(gt_data_cols) == 1) {
      # Strip any curly braces and split by comma
      cleaned_calls <- gsub("[{}]", "", result$calls)
      parsed_list <- strsplit(cleaned_calls, ",")

      sample_names <- NULL

      # Retrieve sample names from metadata
      if ("metadata" %in% all_tables) {
        sample_metadata <- DBI::dbGetQuery(
          con,
          "SELECT lib FROM metadata ORDER BY array_index ASC"
        )

        # Safely extract the first column regardless of case sensitivity
        if (ncol(sample_metadata) > 0) {
          sample_names <- sample_metadata[[1]]
        }
      }


      if (
        is.null(sample_names) ||
          length(sample_names) != length(parsed_list[[1]])
      ) {
        sample_names <- paste0("Sample_", seq_len(length(parsed_list[[1]])))
      }

      # Reconstruct the matrix to hold all sample calls
      gt_matrix <- do.call(rbind, parsed_list)
      colnames(gt_matrix) <- sample_names

      # Combine the metadata columns with the expanded genotype matrix
      info_cols <- result[, setdiff(names(result), "calls"), drop = FALSE]
      result <- cbind(
        info_cols,
        as.data.frame(gt_matrix, stringsAsFactors = FALSE)
      )
    }
    return(result)
  }

  # --- Variants ---
  sql <- "
    SELECT * FROM variants v
    WHERE v.chrom = :chrom"

  params <- list(chrom = chrom)

  if (!is.null(start) && !is.null(end)) {
    sql <- paste(sql, "AND v.pos BETWEEN :start AND :end")
    params$start <- start
    params$end <- end
  }

  sql <- paste(sql, "ORDER BY v.pos")
  result <- DBI::dbGetQuery(con, sql, params = params)

  return(result)
}

#' Extract variants from annotation table based on impact type: LOW, MODERATE,
#' HIGH, MODIFIER.
#' @param db_path A character value indicating the path to the SQLite database.
#' @param impact_level A character value specifying the vriant impact type
#' following snpEff annotation conventions: `HIGH`, `MODERATE`, `LOW`, `MODIFIER`
#' @param chrom A character value specifying the chromosome name.
#' @param start,end A numeric value specifying the start and end coordinates for
#' the candidate gene.
#' @returns A data frame of variants showing the impact type defined in the
#' SQLite database.
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # Extract low impact variant for a region or gene
#' high_variants <- query_by_impact(db_path = mini_db_path,
#'                                impact_level = 'high',
#'                                chrom = "Chr05",
#'                                start = 75104537,
#'                                end = 75106403)
#'
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                              full.names = TRUE,
#'                              recursive = TRUE,
#'                              all.files = TRUE,
#'                              include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#' }
#'
#' @export
#'
query_by_impact <- function(db_path,
                            impact_level = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
                            chrom = NULL,
                            start = NULL,
                            end = NULL) {

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  # assign("con", con, envir = .GlobalEnv)

  # Validate impact levels
  allowed <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  impact_level <- intersect(toupper(impact_level), allowed)

  if (length(impact_level) == 0) {
    stop("No valid impact levels provided.")
  }

  impact_filter <- sprintf("a.impact IN (%s)", paste(sprintf("'%s'", impact_level),
                                                     collapse = ", "))

  region_filter <- ""
  if (!is.null(chrom)) {
    region_filter <- sprintf(" AND v.chrom = '%s'", chrom)
    if (!is.null(start) && !is.null(end)) {
      region_filter <- sprintf("%s AND v.pos BETWEEN %d AND %d", region_filter, start, end)
    }
  }

  query <- sprintf("
    SELECT v.*, a.allele, a.annotation, a.impact, a.gene_name,
           a.gene_id, a.feature_type, a.feature_id, a.transcript_biotype,
           a.rank, a.HGVS_c, a.HGVS_p
    FROM variants v
    JOIN annotations a ON v.variant_id = a.variant_id
    WHERE %s%s
    ORDER BY v.chrom, v.pos", impact_filter, region_filter)

  result <- DBI::dbGetQuery(con, query)
  DBI::dbDisconnect(con)

  return(result)

}


#' Compute allele frequencies for a VCF genotype matrix (variant x samples).
#' Chromosome and position may be included in the data.
#' @param gt A matrix or data frame of variants x samples with chromosome and
#' position meta data for each variant.
#' @param variant_id_col,chrom_col,pos_col A character value specifying the
#' column names of variant IDs, chromosome, and positions in \code{gt}.
#'
#' @returns A data frame of variants with their reference and alternate allele
#' frequencies.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # Extract variants within a candidate gene and their calculate allele
#' # frequencies
#' variant_region <- query_db(db_path = mini_db_path,
#'                               chrom = "Chr05",
#'                               start = 75104537,
#'                               end = 75106403,
#'                               table_name = "genotypes",
#'                               gene_name = "Sobic.005G213600") |>
#
#' calc_af(variant_id_col = 'variant_id',
#'        chrom_col = 'chrom',
#'        pos_col = 'pos')
#'
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                              full.names = TRUE,
#'                              recursive = TRUE,
#'                              all.files = TRUE,
#'                              include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#' }
#'
#' @export
calc_af <- function(gt,
                    variant_id_col = 'variant_id',
                    chrom_col = NULL,
                    pos_col = NULL){

  # if(!is.null(chrom)) chrom <- chrom
  # if(!is.null(pos)) pos <- pos

  # Drop variant_id, chrom, pos from sample matrix
  sample_cols <- setdiff(colnames(gt), c(variant_id_col, chrom_col, pos_col))
  gt_matrix <- gt[, sample_cols, drop = FALSE]

  # Function to convert GT string to ALT dosage (0, 1, 2)
  parse_gt <- function(gt) {
    if (is.na(gt)) return(NA)
    if (gt %in% c("0|0", "0/0")) return(0)
    if (gt %in% c("0|1", "1|0", "0/1", "1/0")) return(1)
    if (gt %in% c("1|1", "1/1")) return(2)
    return(NA)  # For unexpected/missing codes
  }

  # Convert all cells to allele dosages
  dosage_mat <- as.data.frame(apply(gt_matrix, c(1, 2), parse_gt),
                              stringsAsFactors = FALSE)

  # Compute alt allele frequency (sum of alt alleles / 2n samples)
  # Handle single variant case (update to avoid rowMeans error)
  if(nrow(dosage_mat) > 1){
    alt_af <- rowMeans(sapply(dosage_mat, as.numeric), na.rm = TRUE) / 2
  } else{
    alt_af <- mean(sapply(dosage_mat, as.numeric), na.rm = TRUE) / 2
  }


  if(is.null(chrom_col) && is.null(pos_col) ) {

    result <- data.frame(variant_id = gt[, variant_id_col], ref_af = 1 - alt_af,
                         alt_af = alt_af)

  } else {

    result <- data.frame(gt[, c(variant_id_col, chrom_col, pos_col)],
                         ref_af = 1 - alt_af,
                         alt_af = alt_af)

  }

  return(result)

}


#' Filter extracted variants based on alternate allele frequency.
#' @param gt A matrix or data frame of variants x samples with chromosome and
#' position meta data for each variant.
#' @param variant_id_col,chrom_col,pos_col A character value specifying the
#' column names of variant IDs, chromosome, and positions in \code{gt}.
#' @param min_af,max_af A numeric value indicating the minimum and maximum allele
#' frequency thresholds for filtering variants.
#'
#' @returns A data frame of filtered variants with their reference and alternate
#'  allele frequencies.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # Extract variants within a candidate gene and filter by allele frequency
#' variant_region <- query_db(db_path = mini_db_path,
#'                               chrom = "Chr05",
#'                               start = 75104537,
#'                               end = 75106403,
#'                               table_name = "genotypes",
#'                               gene_name = "Sobic.005G213600") |>
#
#' filter_by_af(variant_id_col = 'variant_id',
#'              chrom_col = 'chrom',
#'              pos_col = 'pos',
#'              min_af = 0.01,
#'              max_af = 0.99)
#'
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                              full.names = TRUE,
#'                              recursive = TRUE,
#'                              all.files = TRUE,
#'                              include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#' }
#'
#' @export
filter_by_af <- function(gt,
                         variant_id_col = 'variant_id',
                         chrom_col = 'chrom',
                         pos_col = 'pos',
                         min_af = 0,
                         max_af = 1) {

  if (nrow(gt) == 0) {
    stop("No variant found in the genotype matrix.")
  }

  result <- calc_af(gt,
                    variant_id_col = variant_id_col,
                    chrom_col = chrom_col,
                    pos_col = pos_col)

  # Filter by allele frequency
  result <- result[result$alt_af >= min_af & result$alt_af <= max_af, ]

  return(result)
}

#' Extract variants based on minimum and maximum allele frequencies within a
#' defined region in a SQLite database.
#' @param db_path A character value indicating the path to the SQLite database.
#' @param min_af,max_af A numeric value indicating the minimum and maximum allele
#' frequency thresholds for filtering variants.
#' @param chrom A character value specifying the chromosome name.
#' @param start,end A numeric value specifying the start and end coordinates for
#' the candidate gene.
#'
#' @returns A data frame of filtered variants and their allele frequencies.
#'
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # Query database to filter variants in a region based on alt allele frequency
#' filter_af <- query_by_af(db_path = mini_db_path,
#'                          min_af = 0.01,
#'                          max_af = 0.99,
#'                          chrom = "Chr05",
#'                          start = 75104537,
#'                          end = 75106403)
#'
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                              full.names = TRUE,
#'                              recursive = TRUE,
#'                              all.files = TRUE,
#'                              include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#' }
#'
#' @export
query_by_af <- function(db_path,
                        min_af = 0,
                        max_af = 1,
                        chrom = NULL,
                        start = NULL,
                        end = NULL) {

  if (is.null(chrom)) {
    stop("Chromosome must be specified.")
  }

  gt <- query_db(
    db_path = db_path,
    table_name = "genotypes",
    chrom = chrom,
    start = start,
    end = end
  )

  if (nrow(gt) == 0) {
    message("No genotype data found for specified region.")
    return(data.frame())
  }

  result <- calc_af(gt,
                    variant_id_col = 'variant_id',
                    chrom_col = 'chrom',
                    pos_col = 'pos')

  # Filter by allele frequency
  result <- result[result$alt_af >= min_af & result$alt_af <= max_af, ]

  return(result)
}


#' Query genotypes for one or more variant IDs from a wide-format genotype table.
#' @param db_path The path to your SQLite or PostgreSQL database.
#' @param variant_ids A character vector of variant IDs to query.
#' @param variant_id_col A character value specifying the column name of variant
#' IDs in the SQLite database.
#' @param variants_table,genotypes_table, A character value specifying the
#' column names of variants and genotypes tables, respectively, in the SQLite
#' database. if the `meta_data = NULL`, it returns all the metadata in
#' the variants table.
#' @param meta_data A character vector of metadata to include with the genotype data.
#' @return A data.frame of genotype data in variants x samples format with
#' metadata.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # Extract all genotypes for all samples for any set of variants
#' geno_filter <- query_genotypes(db_path = mini_db_path,
#'                                variant_ids = c("INDEL_Chr05_75104541",
#'                                                "SNP_Chr05_75104557"),
#'                                meta_data = c('chrom', 'pos', 'ref', 'alt',
#'                                              'variant_type'))
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                              full.names = TRUE,
#'                              recursive = TRUE,
#'                              all.files = TRUE,
#'                              include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#' }
#'
#' @export
#'
query_genotypes <- function(
  db_path,
  variant_ids,
  variant_id_col = "variant_id",
  variants_table = 'variants',
  genotypes_table = 'genotypes',
  meta_data = NULL
) {
  if (length(variant_ids) == 0) {
    warning("The 'variant_ids' vector is empty. No query was executed.")
    return(data.frame())
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  all_tables <- DBI::dbListTables(con)
  if (!(variants_table %in% all_tables) || !(genotypes_table %in% all_tables)) {
    stop(
      "The specified variants or genotypes table does not exist in the database."
    )
  }

  if (is.null(meta_data)) {
    meta_data <- DBI::dbListFields(con, variants_table)
  } else {
    meta_data <- unique(c(variant_id_col, meta_data))
  }

  quoted_ids <- paste(sprintf("'%s'", variant_ids), collapse = ",")
  meta_cols_sql <- paste0("v.", meta_data, collapse = ", ")

  gt_cols <- DBI::dbListFields(con, genotypes_table)
  gt_data_cols <- setdiff(gt_cols, c(variant_id_col, "chrom", "pos"))
  gt_cols_sql <- paste0("g.", gt_data_cols, collapse = ", ")

  sql <- sprintf(
    "
    SELECT %s, %s 
    FROM %s v 
    JOIN %s g ON v.%s = g.%s 
    WHERE v.%s IN (%s) 
    ORDER BY v.chrom, v.pos",
    meta_cols_sql,
    gt_cols_sql,
    variants_table,
    genotypes_table,
    variant_id_col,
    variant_id_col,
    variant_id_col,
    quoted_ids
  )

  res <- DBI::dbGetQuery(con, sql)

  if (nrow(res) == 0) {
    warning("No data found for the provided variant IDs.")
    return(data.frame())
  }

  # Unpack array if database uses packed 4-column 'calls' structure
  if ("calls" %in% gt_data_cols && length(gt_data_cols) == 1) {
    cleaned_calls <- gsub("[{}]", "", res$calls)
    parsed_list <- strsplit(cleaned_calls, ",")

    sample_names <- NULL
    if ("metadata" %in% all_tables) {
      sample_metadata <- DBI::dbGetQuery(
        con,
        "SELECT lib FROM metadata ORDER BY array_index ASC"
      )
      if (ncol(sample_metadata) > 0) sample_names <- sample_metadata[[1]]
    }

    if (
      is.null(sample_names) || length(sample_names) != length(parsed_list[[1]])
    ) {
      sample_names <- paste0("Sample_", seq_len(length(parsed_list[[1]])))
    }

    gt_matrix <- do.call(rbind, parsed_list)
    colnames(gt_matrix) <- sample_names

    info_cols <- res[, setdiff(names(res), "calls"), drop = FALSE]
    result <- cbind(
      info_cols,
      as.data.frame(gt_matrix, stringsAsFactors = FALSE)
    )
    return(result)
  }

  return(res)
}


#' Count the number of variant types in the SQLite database.
#' @param db_path The path to your SQLite or PostgreSQL database.
#' @return A data.frame of counts for variant types.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # Count the number of variant types in database
#' variant_type_count <- count_variant_types(mini_db_path)
#'
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                              full.names = TRUE,
#'                              recursive = TRUE,
#'                              all.files = TRUE,
#'                              include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#' }
#'
#' @export
#'
count_variant_types <- function(db_path) {

  table_name <- "variants"

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  # assign("con", con, envir = .GlobalEnv)

  # Check if column `variant_type` exists
  cols <- DBI::dbListFields(con, table_name)
  if (!"variant_type" %in% cols) {
    DBI::dbDisconnect(con)
    stop("The table does not have a 'variant_type' column.")
  }

  # Count by variant_type
  query <- sprintf(
    "SELECT variant_type, COUNT(*) AS n FROM %s GROUP BY variant_type ORDER BY variant_type",
    table_name
  )
  result <- DBI::dbGetQuery(con, query)
  DBI::dbDisconnect(con)

  return(result)
}

#' Query the annotations table within a specified genomic region and summarize
#' the distribution of SnpEff annotations and impact categories by variant type.
#' @param db_path A character value indicating the path to the SQLite database.
#' @param chrom A character value specifying the chromosome name.
#' @param start,end A numeric value specifying the start and end coordinates for
#' the candidate gene.
#' @param annotations_table,variants_table A character value indicating the
#' table names for snpEff annotations and variants in the SQLite database.
#'
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#'
#' # Define tempdir
#' path <- tempdir()
#'
#' # Mini SQLite database
#' mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#'
#'# Path to SQLite databases: INDEL and SNP
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#'
#' # Unzip compressed mini database and save in tempdir
#' R.utils::gunzip(mini_db,
#'                destname = mini_db_path,
#'                remove = FALSE)
#'
#' # Annotation and impact summary distribution within Chr05: 75104537-75106403
#' ann_summary <- query_ann_summary(db_path = mini_db_path,
#'                                  chrom = "Chr05",
#'                                  start = 75104537,
#'                                  end = 75106403)
#'
#' # Clean tempdir
#' contents <- list.files(tempdir(),
#'                        full.names = TRUE,
#'                        recursive = TRUE,
#'                        all.files = TRUE,
#'                        include.dirs = TRUE)
#' unlink(contents, recursive = TRUE, force = TRUE)
#' }
#'
#' @returns A list with:
#' \itemize{
#'   \item \code{annotation_summary}
#'   \item \code{impact_summary}
#'   \item \code{variant_type_totals}
#' }
#' @export
#'
query_ann_summary <- function(db_path,
                              chrom,
                              start,
                              end,
                              annotations_table = "annotations",
                              variants_table = "variants") {

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Validate required fields exist
  ann_fields <- DBI::dbListFields(con, annotations_table)
  var_fields <- DBI::dbListFields(con, variants_table)

  required_ann <- c("variant_id", "annotation", "impact")
  required_var <- c("variant_id", "chrom", "pos", "variant_type")

  if (!all(required_ann %in% ann_fields)) {
    stop("Annotations table must contain: ", paste(required_ann, collapse = ", "))
  }
  if (!all(required_var %in% var_fields)) {
    stop("Variants table must contain: ", paste(required_var, collapse = ", "))
  }

  # Query variant_type counts for region
  total_query <- sprintf("
    SELECT variant_type, COUNT(*) AS total_variants
    FROM %s
    WHERE chrom = '%s' AND pos BETWEEN %d AND %d
    GROUP BY variant_type", variants_table, chrom, start, end)

  variant_type_totals <- DBI::dbGetQuery(con, total_query)

  # Main annotation + impact query
  query <- sprintf(
    "SELECT a.annotation, a.impact, v.variant_type
     FROM %s a
     JOIN %s v ON a.variant_id = v.variant_id
     WHERE v.chrom = '%s' AND v.pos BETWEEN %d AND %d",
    annotations_table, variants_table, chrom, start, end
  )

  result <- DBI::dbGetQuery(con, query)
  DBI::dbDisconnect(con)

  if (nrow(result) == 0) {
    warning("No annotations found in the specified region.")
    return(list(annotation_summary = data.frame(),
                impact_summary = data.frame(),
                variant_type_totals = variant_type_totals))
  }

  # Annotation summary
  annotation_summary <- as.data.frame(table(result$annotation, result$variant_type))
  names(annotation_summary) <- c("annotation", "variant_type", "count")
  annotation_summary <- annotation_summary[annotation_summary$count > 0, ]
  annotation_summary <- annotation_summary[order(-annotation_summary$count), ]

  # Impact summary
  impact_summary <- as.data.frame(table(result$impact, result$variant_type))
  names(impact_summary) <- c("impact", "variant_type", "count")
  impact_summary <- impact_summary[impact_summary$count > 0, ]
  impact_summary <- impact_summary[order(-impact_summary$count), ]

  return(list(
    annotation_summary = annotation_summary,
    impact_summary = impact_summary,
    variant_type_totals = variant_type_totals
  ))
}



#' Retrieve sample metadata from the SQLite database
#'
#' This function connects to the local SQLite database and retrieves accession-level
#' metadata from the 'metadata' table. It supports optional filtering by specific
#' columns (e.g., country of origin, race, or cluster) to easily subset populations
#' for downstream analysis.
#'
#' @param db_path A character string specifying the file path to the SQLite database.
#' @param query_col Character. The metadata column to filter by (e.g., "countryorigin").
#'   If \code{NULL}, all records are returned.
#' @param query_value Character or Numeric. The specific value to match in \code{query_col}.
#'
#' @returns A data frame containing the sample metadata records, ordered by their
#'   array index to seamlessly match genotype matrices.
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#'
#' # Define tempdir and setup mini SQLite database
#' path <- tempdir()
#' mini_db <- system.file(
#'   "extdata", "mini_sorghum_variant_vcf.db.gz", 
#'   package = "panGenomeBreedr", mustWork = TRUE
#' )
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#' R.utils::gunzip(mini_db, destname = mini_db_path, remove = FALSE)
#'
#' # Fetch metadata for all accessions originating from Ethiopia
#' eth_samples <- get_sample_metadata(
#'   db_path = mini_db_path,
#'   query_col = "countryorigin",
#'   query_value = "Ethiopia"
#' )
#'
#' # Clean tempdir
#' unlink(list.files(tempdir(), full.names = TRUE, recursive = TRUE, 
#'                   include.dirs = TRUE), recursive = TRUE, force = TRUE)
#' }
#'
#' @export
#' @import DBI
#' @importFrom RSQLite SQLite
get_sample_metadata <- function(db_path, query_col = NULL, query_value = NULL) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  if (!("metadata" %in% DBI::dbListTables(con))) {
    stop("Table 'metadata' not found in the SQLite database.")
  }

  if (!is.null(query_col)) {
    available_cols <- DBI::dbListFields(con, "metadata")
    if (!query_col %in% available_cols) {
      stop(sprintf("Column '%s' not found in the metadata table.", query_col))
    }
  }

  sql <- "SELECT * FROM metadata"
  params <- list()

  if (!is.null(query_col) && !is.null(query_value)) {
    # Using parameterized SQL for safety
    sql <- paste0(sql, " WHERE ", query_col, " = :val")
    params$val <- query_value
  }

  sql <- paste0(sql, " ORDER BY array_index ASC")

  if (length(params) > 0) {
    result <- DBI::dbGetQuery(con, sql, params = params)
  } else {
    result <- DBI::dbGetQuery(con, sql)
  }
  return(result)
}


#' Query genotypes filtered by sample metadata attributes in SQLite
#'
#' This function retrieves genotypes for a specific genomic region, but restricts
#' the output to a subset of samples defined by metadata attributes (e.g., extracting
#' variants only for specific countries, populations, or phenotypic clusters).
#'
#' @param db_path A character string specifying the file path to the SQLite database.
#' @param chrom Character. Chromosome name (e.g., "Chr05").
#' @param start,end Numeric. The genomic start and end coordinates.
#' @param meta_col Character. Metadata column to filter samples by (e.g., "countryorigin").
#' @param meta_value Character. Value to match in \code{meta_col} (e.g., "Mali").
#'
#' @returns A wide-format data frame (variants x filtered samples).
#'
#' @examples
#' \dontrun{
#' library(panGenomeBreedr)
#'
#' # Define tempdir and setup mini SQLite database
#' path <- tempdir()
#' mini_db <- system.file(
#'   "extdata", "mini_sorghum_variant_vcf.db.gz", 
#'   package = "panGenomeBreedr", mustWork = TRUE
#' )
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#' R.utils::gunzip(mini_db, destname = mini_db_path, remove = FALSE)
#'
#' # Get genotypes for a specific gene region, but only for samples from Mali
#' mali_genotypes <- query_by_metadata(
#'   db_path = mini_db_path,
#'   chrom = "Chr05",
#'   start = 75104537,
#'   end = 75106403,
#'   meta_col = "countryorigin",
#'   meta_value = "Mali"
#' )
#'
#' # Clean tempdir
#' unlink(list.files(tempdir(), full.names = TRUE, recursive = TRUE, 
#'                   include.dirs = TRUE), recursive = TRUE, force = TRUE)
#' }
#'
#' @export
query_by_metadata <- function(
  db_path,
  chrom,
  start,
  end,
  meta_col,
  meta_value
) {
  sample_info <- get_sample_metadata(db_path, meta_col, meta_value)

  if (nrow(sample_info) == 0) {
    stop(paste(
      "No samples found matching the criteria:",
      meta_col,
      "=",
      meta_value
    ))
  }

  sample_names <- sample_info$lib

  # Fetch the raw data using your robust query_db function
  raw_data <- query_db(
    db_path = db_path,
    table_name = "genotypes",
    chrom = chrom,
    start = start,
    end = end
  )

  if (nrow(raw_data) == 0) {
    return(data.frame())
  }

  # Ensure we only extract samples that actually exist in the genotype matrix
  valid_samples <- intersect(sample_names, colnames(raw_data))
  result <- raw_data[, c("variant_id", "chrom", "pos", valid_samples)]

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
#' # Define tempdir and setup mini SQLite database
#' path <- tempdir()
#' mini_db <- system.file(
#'   "extdata", "mini_sorghum_variant_vcf.db.gz", 
#'   package = "panGenomeBreedr", mustWork = TRUE
#' )
#' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
#' R.utils::gunzip(mini_db, destname = mini_db_path, remove = FALSE)
#'
#' # Fetch metadata for all accessions originating from Ethiopia
#' eth_samples <- get_sample_metadata(
#'   db_path = mini_db_path,
#'   query_col = "NULL",
#'   query_value = "NULL"
#' )
#'
#' # Explore the geographic distribution colored by genetic cluster
#' query_map_accessions(meta, color_by = "kmeans_cluster")
#' }
#'
#' @export
query_map_accessions <- function(metadata, color_by = "countryorigin") {
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
    stop(
      "Metadata must contain 'lat' and 'lon' columns for geographic mapping."
    )
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
      color = ~ pal(plot_data[[color_by]]),
      stroke = FALSE,
      fillOpacity = 0.8,
      popup = popup_info,
      label = ~lib
    )

  return(map)
}


# #' Generic Trait Association Audit (Local SQLite Version)
# #'
# #' This function performs an association analysis between a specific genomic variant
# #' and a phenotypic trait. It pulls genotypes directly from the local SQLite database,
# #' merges them with provided phenotypic data, and generates a statistical summary plot
# #' using Kruskal-Wallis (quantitative) or Chi-Square (qualitative) tests.
# #'
# #' @param db_path A character string specifying the file path to the SQLite database.
# #' @param variant_id Character. The specific SNP or INDEL ID to analyze (e.g., "SNP_Chr05_75104557").
# #' @param pheno_df Data frame containing phenotypic data. The first column must contain accession identifiers (e.g., PI numbers).
# #' @param trait_col Character. The name of the column in \code{pheno_df} containing trait scores.
# #' @param trait_type Character. Either "qualitative" or "quantitative". Determines the statistical test and plot type.
# #'
# #' @returns A \code{ggplot2} object displaying the association boxplot/bar chart and a statistical summary.
# #'
# #' @examples
# #' \donttest{
# #' library(panGenomeBreedr)
# #'
# #' # Define tempdir and setup mini SQLite database
# #' path <- tempdir()
# #' mini_db <- system.file("extdata", "mini_sorghum_variant_vcf.db.gz", package = "panGenomeBreedr", mustWork = TRUE)
# #' mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
# #' R.utils::gunzip(mini_db, destname = mini_db_path, remove = FALSE)
# #'
# #' # 1. Fetch valid sample IDs to create dummy phenotype data
# #' valid_samples <- get_sample_metadata(mini_db_path)$pinumber
# #'
# #' # 2. Create a dummy quantitative phenotype dataframe
# #' set.seed(123)
# #' dummy_pheno <- data.frame(
# #'   pinumber = valid_samples,
# #'   plant_height = rnorm(length(valid_samples), mean = 150, sd = 20)
# #' )
# #'
# #' # 3. Run the association audit
# #' plot_trait_association(
# #'   db_path = mini_db_path,
# #'   variant_id = "SNP_Chr05_75104557",
# #'   pheno_df = dummy_pheno,
# #'   trait_col = "plant_height",
# #'   trait_type = "quantitative"
# #' )
# #'
# #' # Clean tempdir
# #' unlink(list.files(tempdir(), full.names = TRUE, recursive = TRUE, include.dirs = TRUE), recursive = TRUE, force = TRUE)
# #' }
# #'
# #' @export
# #' @import ggplot2
# #' @import stats
# #' @import DBI
# #' @importFrom RSQLite SQLite
# plot_trait_association <- function(
#   db_path,
#   variant_id,
#   pheno_df,
#   trait_col,
#   trait_type = c("qualitative", "quantitative")
# ) {
#   trait_type <- match.arg(trait_type)
#   genotype <- NULL

#   if (
#     trait_type == "qualitative" && !requireNamespace("scales", quietly = TRUE)
#   ) {
#     stop("The 'scales' package is required. Please install it.")
#   }

#   # 1. Pull genotype data directly from SQLite
#   gt_data <- query_genotypes(db_path = db_path, variant_ids = variant_id)

#   if (nrow(gt_data) == 0) {
#     stop(paste("Variant ID", variant_id, "was not found."))
#   }

#   # 2. Fetch functional annotations for the verdict logic
#   con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
#   ann_info <- DBI::dbGetQuery(
#     con,
#     "SELECT annotation, impact FROM annotations WHERE variant_id = :id",
#     params = list(id = variant_id)
#   )
#   DBI::dbDisconnect(con)

#   # 3. Fetch sample metadata mapping
#   pi_map <- get_sample_metadata(db_path)

#   # Align genotype samples with phenotype PI numbers
#   sample_cols <- setdiff(
#     names(gt_data),
#     c(
#       "variant_id",
#       "chrom",
#       "pos",
#       "ref",
#       "alt",
#       "qual",
#       "filter",
#       "variant_type"
#     )
#   )

#   gt_long <- data.frame(
#     pinumber = pi_map$pinumber[match(sample_cols, pi_map$lib)],
#     genotype = as.character(unlist(gt_data[1, sample_cols])),
#     stringsAsFactors = FALSE
#   )

#   colnames(pheno_df)[1] <- "pinumber"
#   plot_data <- merge(gt_long, pheno_df, by = "pinumber")

#   plot_data <- plot_data[
#     !is.na(plot_data[[trait_col]]) &
#       plot_data[[trait_col]] != "-" &
#       plot_data$genotype != "./.",
#   ]

#   if (trait_type == "quantitative") {
#     if (!is.numeric(plot_data[[trait_col]])) {
#       plot_data[[trait_col]] <- as.numeric(gsub(
#         ".*[^0-9.]",
#         "",
#         as.character(plot_data[[trait_col]])
#       ))
#     }
#   } else {
#     plot_data[[trait_col]] <- as.factor(plot_data[[trait_col]])
#   }

#   alleles <- unlist(strsplit(plot_data$genotype, "[|/]"))
#   maf <- min(table(alleles)) / length(alleles)
#   p_val <- NA
#   pve_val <- 0

#   if (trait_type == "quantitative") {
#     p_val <- stats::kruskal.test(
#       plot_data[[trait_col]] ~ plot_data$genotype
#     )$p.value
#     pve_val <- summary(stats::lm(
#       plot_data[[trait_col]] ~ plot_data$genotype
#     ))$r.squared
#   } else {
#     p_val <- stats::chisq.test(
#       table(plot_data$genotype, plot_data[[trait_col]]),
#       simulate.p.value = TRUE
#     )$p.value
#   }

#   impact <- if (nrow(ann_info) > 0) ann_info$impact[1] else "UNKNOWN"
#   is_sig <- !is.na(p_val) && p_val < 0.05

#   verdict <- if (is_sig && (pve_val > 0.1 || impact == "HIGH")) {
#     "VALIDATE"
#   } else if (is_sig) {
#     "CAUTION"
#   } else {
#     "DISCARD"
#   }

#   p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = genotype)) +
#     ggplot2::theme_minimal() +
#     ggplot2::scale_fill_brewer(palette = "Set1")

#   if (trait_type == "quantitative") {
#     p <- p +
#       ggplot2::aes(y = .data[[trait_col]], fill = genotype) +
#       ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#       ggplot2::geom_jitter(width = 0.2, alpha = 0.4)
#   } else {
#     p <- p +
#       ggplot2::aes(fill = .data[[trait_col]]) +
#       ggplot2::geom_bar(position = "fill", color = "white", linewidth = 0.3) +
#       ggplot2::scale_y_continuous(labels = scales::percent)
#   }

#   p <- p +
#     ggplot2::labs(
#       title = paste("Association Audit:", variant_id),
#       subtitle = paste0(
#         "P-val: ",
#         format.pval(p_val, digits = 3),
#         " | PVE: ",
#         round(pve_val * 100, 1),
#         "% | MAF: ",
#         round(maf, 3),
#         "\nImpact: ",
#         impact,
#         " | Verdict: ",
#         verdict
#       ),
#       x = "Genotype",
#       y = if (trait_type == "quantitative") trait_col else "Frequency (%)"
#     )

#   return(p)
# }