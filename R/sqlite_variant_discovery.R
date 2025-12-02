
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

#'Query any table in your SQLite database using chromosome and a genomic position
#' range.
#' @param db_path A character value indicating the path to the SQLite database.
#' @param table_name A character value specifying the name of any any table in a
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
#'
query_db <- function(db_path,
                     table_name = c("variants", "annotations", "genotypes"),
                     chrom,
                     start = NULL,
                     end = NULL,
                     gene_name = NULL) {

  table_name <- match.arg(table_name)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  # assign("con", con, envir = .GlobalEnv)

  all_tables <- DBI::dbListTables(con)
  if (!(table_name %in% all_tables)) {
    DBI::dbDisconnect(con)
    stop(sprintf("Table '%s' does not exist in the database.", table_name))
  }

  # Build WHERE clause
  where_clause <- sprintf("v.chrom = '%s'", chrom)
  if (!is.null(start) && !is.null(end)) {
    where_clause <- sprintf("%s AND v.pos BETWEEN %d AND %d", where_clause, start, end)
  }

  # Build query
  if (table_name == "annotations") {
    query <- sprintf("
      SELECT a.*, v.chrom, v.pos
      FROM annotations a
      JOIN variants v ON a.variant_id = v.variant_id
      WHERE %s
      ORDER BY v.pos", where_clause)

  } else if (table_name == "genotypes") {
    query <- sprintf("
      SELECT g.variant_id, g.chrom, g.pos, v.variant_type, v.ref, v.alt, %s
      FROM genotypes g
      JOIN variants v ON g.variant_id = v.variant_id
      WHERE %s
      ORDER BY g.pos",
      paste(setdiff(DBI::dbListFields(con, "genotypes"),
                                       c("variant_id", "chrom", "pos")),
            collapse = ", "), where_clause)


  } else {
    query <- sprintf("
      SELECT * FROM variants v
      WHERE %s
      ORDER BY v.pos", where_clause)
  }

  result <- DBI::dbGetQuery(con, query)
  DBI::dbDisconnect(con)

  # Extract variants within candidate gene if gene name is not NULL
  if (!is.null(gene_name) && table_name == "annotations") {

    result <- result[result$gene_name == gene_name,]

  }

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
  gt_matrix <- gt[, sample_cols]

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
  alt_af <- rowMeans(sapply(dosage_mat, as.numeric), na.rm = TRUE) / 2

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

  # Connect to database
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  # assign("con", con, envir = .GlobalEnv)

  # Build WHERE clause
  where_clause <- ""
  if (!is.null(chrom)) {
    where_clause <- sprintf("v.chrom = '%s'", chrom)
    if (!is.null(start) && !is.null(end)) {
      where_clause <- sprintf("%s AND v.pos BETWEEN %d AND %d", where_clause, start, end)
    }
  }

  query <- sprintf("
      SELECT g.*
      FROM genotypes g
      JOIN variants v ON g.variant_id = v.variant_id
      WHERE %s
      ORDER BY v.pos", where_clause)

  gt <- DBI::dbGetQuery(con, query)
  DBI::dbDisconnect(con)

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
query_genotypes <- function(db_path,
                            variant_ids,
                            variant_id_col = "variant_id",
                            variants_table = 'variants',
                            genotypes_table = 'genotypes',
                            meta_data = NULL) {

  if (length(variant_ids) == 0) {
    stop("Please provide at least one variant ID in the database.")
  }

  # Connect to the database
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Use all metadata if not specified
  if (is.null(meta_data)) {
    meta_data <- setdiff(DBI::dbListFields(con, variants_table), variant_id_col)
  }

  # Ensure join key is included in metadata
  meta_data <- unique(c(variant_id_col, meta_data))

  # Format variant IDs for SQL
  id_list <- paste(sprintf("'%s'", variant_ids), collapse = ", ")

  # Query genotype data
  geno_query <- sprintf("SELECT * FROM %s WHERE %s IN (%s)",
                        genotypes_table, variant_id_col, id_list)
  geno_df <- DBI::dbGetQuery(con, geno_query)

  # Check available metadata columns
  all_meta_cols <- DBI::dbListFields(con, variants_table)
  valid_meta_cols <- intersect(meta_data, all_meta_cols)

  if (!(variant_id_col %in% valid_meta_cols)) {
    stop(sprintf("'%s' must be included in meta_data for proper joining.", variant_id_col))
  }

  # Warn if user requested unavailable columns
  missing_cols <- setdiff(meta_data, all_meta_cols)
  if (length(missing_cols) > 0) {
    warning("The following metadata columns do not exist in the 'variants' table and will be ignored: ",
            paste(missing_cols, collapse = ", "))
  }

  # Query metadata table
  meta_query <- sprintf("SELECT %s FROM %s WHERE %s IN (%s)",
                        paste(valid_meta_cols, collapse = ", "),
                        variants_table, variant_id_col, id_list)
  meta_df <- DBI::dbGetQuery(con, meta_query)

  DBI::dbDisconnect(con)

  if (nrow(geno_df) == 0) {
    warning("No genotype data found for the specified variant IDs.")
    return(data.frame())
  }

  # Prevent duplication of shared columns except the join key
  common_cols <- intersect(names(geno_df), names(meta_df))
  cols_to_drop <- setdiff(common_cols, variant_id_col)
  geno_df <- geno_df[, !names(geno_df) %in% cols_to_drop]

  # Merge metadata and genotypes
  wide_df <- merge(meta_df, geno_df, by = variant_id_col, all.x = TRUE)

  # Drop duplicates if any (safety net)
  wide_df <- wide_df[, !duplicated(names(wide_df))]

  return(wide_df)
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

