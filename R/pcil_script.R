#' Fetch All PCIL Data Tables
#'
#' @description
#' Fetches all tables related to Pangenome-Characterized Introgression Lines (PCIL)
#' and returns them as a single list object. This function serves as a "load all" utility
#' for PCIL data, supporting both local and online database modes.
#'
#' @param con A DBI connection object to the local database. Required when `connect_db_mode = 'local'`.
#' @param connect_db_mode Character string specifying the connection mode: `'local'` or `'online'`.
#'
#' @return A list containing the following data frames:
#' \itemize{
#'   \item \code{pcil_metadata}: Metadata for PCIL families.
#'   \item \code{pcil_gene_regions}: Genomic regions of genes within the PCILs.
#'   \item \code{pcil_introgressions}: Introgression block data.
#'   \item \code{pcil_genomewide_introgressions}: Genome-wide introgression summary statistics.
#'   \item \code{pcil_inbreeding_coefficient}: Sample metadata, often used for inbreeding coefficients.
#'   \item \code{pcil_IBS_dis}: Identity-by-state (IBS) distance matrix.
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#'
#' # --- Online Mode ---
#' # Fetch all PCIL data tables from the remote database
#' pcil_data <- fetch_pcil_data(connect_db_mode = 'online')
#' }
#'
#' \dontrun{
#' # --- Offline Mode ---
#' # Requires a local database folder that contains a "pcil" subfolder;
#' # the sample dataset shipped with this package does not include one.
#' con <- connect_local_db(folder_path = "path/to/your/local_db")
#' pcil_data_local <- fetch_pcil_data(con = con, connect_db_mode = "local")
#' disconnect_local_db(con)
#' }
fetch_pcil_data <- function(
  con = NULL,
  connect_db_mode = c("local", "online")
) {
  connect_db_mode <- match.arg(connect_db_mode)

  if (connect_db_mode == "online") {
    # Fetch all PCIL tables from the API
    pcil_data <- list(
      pcil_metadata = .api_fetch("/db/pcil_metadata"),
      pcil_gene_regions = .api_fetch("/db/pcil_genes"),
      pcil_introgressions = .api_fetch("/db/pcil_introgressions"),
      pcil_genomewide_introgressions = .api_fetch("/db/pcil_genomewide"),
      pcil_inbreeding_coefficient = .api_fetch("/db/pcil_sample_metadata"),
      pcil_IBS_dis = .api_fetch("/db/pcil_ibs")
    )
  } else {
    # Fetch all PCIL tables from the local DuckDB connection
    if (
      is.null(con) || !inherits(con, "DBIConnection") || !DBI::dbIsValid(con)
    ) {
      stop("A valid 'con' object is required for 'local' mode.", call. = FALSE)
    }

    pcil_data <- list(
      pcil_metadata = DBI::dbGetQuery(con, "SELECT * FROM pcil_metadata"),
      pcil_gene_regions = DBI::dbGetQuery(con, "SELECT * FROM pcil_genes"),
      pcil_introgressions = DBI::dbGetQuery(
        con,
        "SELECT * FROM pcil_introgressions"
      ),
      pcil_genomewide_introgressions = DBI::dbGetQuery(con, "SELECT * FROM pcil_genomewide"),
      pcil_inbreeding_coefficient = DBI::dbGetQuery(con, "SELECT * FROM pcil_sample_metadata"),
      pcil_IBS_dis = DBI::dbGetQuery(con, "SELECT * FROM pcil_ibs")
    )
  }

  return(pcil_data)
}



#' Fetch Pangenome-Characterized Introgression Lines (PCIL) by Variant
#'
#' @description
#' This function identifies Pangenome-Characterized Introgression Lines (PCIL)
#' that are putative donors for a given variant of interest. It operates by
#' comparing the genotypes of recurrent parents against a pool of potential
#' donors to find families that carry the alternate allele.
#'
#' @param con A DBI connection object to the local database.
#' @param selection A character vector of variant IDs.
#' @param pcil_data A list object containing PCIL data, typically from
#'   `fetch_pcil_data()`. Must contain the `pcil_metadata` table.
#' @param connect_db_mode A character string specifying the connection mode.
#'   Can be either `'local'` (default) or `'online'`.
#'
#' @return A list containing summaries of PCIL families, genotypes, and selections.
#' \itemize{
#'   \item \code{pcil_family_summary}: A summary of PCIL families, including counts of
#'     families and lines per recurrent parent that are putative sources for the
#'     alternate allele.
#'   \item \code{geno_pi}: Genotype information for the parental lines.
#'   \item \code{pcil_summary}: Detailed information on the selected PCIL families.
#' }
#' Returns `NULL` if no matching PCIL families are found.
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#' 
#' # 1. Fetch PCIL Data
#' pcil_data <- fetch_pcil_data(connect_db_mode = 'online')
#' 
#' # 2. Define variants of interest
#' selection <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")
#' 
#' # 3. Fetch PCIL families acting as putative donors for these variants
#' results <- fetch_pcil_families_by_variant(
#'   selection = selection,
#'   pcil_data = pcil_data,
#'   connect_db_mode = 'online'
#' )
#' 
#' print(results$pcil_family_summary)
#' }
fetch_pcil_families_by_variant <- function(
  con = NULL,
  selection,
  pcil_data,
  connect_db_mode = c('local', 'online')
) {
  # Match argument for connection mode
  connect_db_mode <- match.arg(connect_db_mode)

  # Validate and extract pcil_metadata from pcil_data
  if (missing(pcil_data) || !is.list(pcil_data) || !"pcil_metadata" %in% names(pcil_data)) {
    stop("`pcil_data` must be a list from fetch_pcil_data() containing 'pcil_metadata'.", call. = FALSE)
  }
  metadata_pcil <- pcil_data$pcil_metadata

  parents_pcils <- unique(stats::na.omit(c(metadata_pcil$clan_rs_Lib, metadata_pcil$fam_rs_Lib)))
  cols <- c("variant_id", "chrom", "pos", "ref", "alt", "variant_type", parents_pcils)

  metadata <- fetch_accession_metadata(con = con, connect_db_mode = connect_db_mode)
  metadata <- metadata[metadata$lib %in% parents_pcils, ]
  metadata <- metadata[, c("pinumber", "lib", "sample")]

  variant_geno <- fetch_genotypes_by_id(
    con = con,
    connect_db_mode = connect_db_mode,
    variant_ids = selection,
    meta_data = c("variant_id", "chrom", "pos", "ref", "alt", "variant_type")
  )

  cols_present <- intersect(cols, colnames(variant_geno))
  variant_geno <- variant_geno[, cols_present]
  geno_cols <- intersect(parents_pcils, colnames(variant_geno))

  # Rename calls
  variant_geno[, geno_cols] <- lapply(variant_geno[, geno_cols], function(x) {
    x <- gsub("0|0", "Reference", x, fixed = TRUE)
    x <- gsub("1|1", "Alternate", x, fixed = TRUE)
    x <- gsub("0|1", "Heterozygous", x, fixed = TRUE)
    x <- gsub("1|0", "Heterozygous", x, fixed = TRUE)
    x <- gsub("./.", NA, x, fixed = TRUE)
    return(x)
  })

  geno_t <- as.data.frame(t(variant_geno[, geno_cols]), stringsAsFactors = FALSE)
  colnames(geno_t) <- variant_geno$variant_id
  geno_t$lib <- rownames(geno_t)
  geno_t <- geno_t[, c("lib", variant_geno$variant_id)]
  rownames(geno_t) <- NULL

  geno_t <- merge(geno_t, metadata, by = "lib", all.x = TRUE)
  geno_pi <- geno_t[!is.na(geno_t$pinumber), ]
  geno_pi <- geno_pi[, c(selection, "lib", "pinumber", "sample")]

  recurrent_parents <- unique(metadata_pcil$clan_rs_Lib)
  rp_genos <- geno_pi[geno_pi$lib %in% recurrent_parents, ]
  rp_meta <- unique(metadata_pcil[, c("clan", "clan_rs_Lib")])
  names(rp_meta) <- c("clan", "lib")
  rp_alleles <- merge(rp_genos, rp_meta, by = "lib", all.x = TRUE)

  pcil_summary_list <- list()
  for (v in selection) {
    for (i in 1:nrow(rp_genos)) {
      rp <- rp_genos[i, ]
      if (is.na(rp[[v]])) next
      putative_donors <- geno_pi[
        !is.na(geno_pi[[v]]) &
          geno_pi[[v]] != rp[[v]] &
          geno_pi[[v]] != "Heterozygous",
      ]
      if (nrow(putative_donors) >= 1) {
        pcil_rp <- metadata_pcil[metadata_pcil$clan_rs_Lib == rp$lib, ]
        selected_pcil <- pcil_rp[pcil_rp$fam_rs_Lib %in% putative_donors$lib, ]
        if (nrow(selected_pcil) >= 1) {
          selected_pcil$PEDIGREE <- paste(selected_pcil$clan, selected_pcil$family, sep = "/")
          output_selected <- selected_pcil[, c("clan", "clan_pi", "family", "family_pi", "lib_id", "sample_id")]
          output_selected$selection <- v
          pcil_summary_list[[length(pcil_summary_list) + 1]] <- output_selected
        }
      }
    }
  }

  if (length(pcil_summary_list) == 0) {
    cat("No PCIL families for", selection)
    return(NULL)
  }

  pcil_summary <- do.call(rbind, pcil_summary_list)

  agg_families <- stats::aggregate(family ~ selection + clan, data = pcil_summary, FUN = function(x) length(unique(x)))
  names(agg_families)[3] <- "families"
  agg_lines <- stats::aggregate(sample_id ~ selection + clan, data = pcil_summary, FUN = length)
  names(agg_lines)[3] <- "lines"
  pcil_family_summary <- merge(agg_families, agg_lines, by = c("selection", "clan"))

  tmp <- as.data.frame(rp_genos[, c(selection, "lib")])
  tmp$lib <- gsub("ISMB", "CSM-63", tmp$lib)
  tmp$lib <- gsub("ITIK", "IRAT204", tmp$lib)
  tmp$lib <- gsub("SAMEA12302677", "IRAT204", tmp$lib)
  tmp$lib <- gsub("ITIL", "Mota Maradi", tmp$lib)
  tmp$lib <- gsub("ISGP", "Macia", tmp$lib)
  names(tmp)[names(tmp) == "lib"] <- "clan"

  tmp_long_list <- lapply(selection, function(sel_col) {
    df_slice <- tmp[, c("clan", sel_col)]
    names(df_slice)[2] <- "rp_genotype"
    df_slice$selection <- sel_col
    df_slice
  })
  tmp_long <- do.call(rbind, tmp_long_list)

  pcil_family_summary <- merge(pcil_family_summary, tmp_long, by = c("clan", "selection"), all.x = TRUE)
  pcil_family_summary <- pcil_family_summary[, c("clan", "families", "lines", "rp_genotype", "selection")]

  pcil_family_summary$recurrent_allele <- NA
  for (v_name in selection) {
    current_variant_info <- variant_geno[variant_geno$variant_id == v_name, ]
    ref_allele <- current_variant_info$ref[1]
    alt_allele <- current_variant_info$alt[1]
    pcil_family_summary$recurrent_allele[
      pcil_family_summary$selection == v_name & pcil_family_summary$rp_genotype == "Reference"
    ] <- ref_allele
    pcil_family_summary$recurrent_allele[
      pcil_family_summary$selection == v_name & pcil_family_summary$rp_genotype == "Alternate"
    ] <- alt_allele
  }

  return(list(
    pcil_family_summary = pcil_family_summary,
    geno_pi = geno_pi,
    pcil_summary = pcil_summary
  ))
}




#' Fetch PCIL Positive Lines
#'
#' @description
#' Identifies and filters Pangenome-Characterized Introgression Lines (PCIL)
#' containing introgressions in a target gene or position window.
#' Supports multi-criteria ranking based on mean donor fraction, haploblock size,
#' and inbreeding coefficients ($F$).
#'
#' @param pcil_data A list object containing all required PCIL data tables,
#'   typically from `fetch_pcil_data()`.
#' @param variants_select_geno A character vector (gene IDs) or a data frame containing position definitions.
#' @param type Character string specifying the target type: `"gene"` or `"position"`.
#' @param sel Integer. Number of top-ranked candidate lines to select per region. Defaults to `NULL` (returns all).
#' @param donor_thresh Numeric. Minimum mean donor fraction required. Defaults to `0.75`.
#' @param block_quantile Numeric. Upper quantile threshold to exclude large introgressions. Defaults to `0.75`.
#' @param F_quantile Numeric. Lower quantile threshold for the inbreeding coefficient ($F$). Defaults to `0.25`.
#' @param window Integer. Symmetric base pair window ($+/-$) around targets when `type = "position"`. Defaults to `100`.
#' @param available_ids Data frame with columns `sample_id` and `selection` to restrict line searches per target.
#' @param result_pcil_families list object of results obtained from `fetch_pcil_families_by_variant`.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{pcil_positive}: Data frame of all matching positive introgression lines merged with genome-wide stats.
#'   \item \code{best_lines}: (Optional) Top-ranked candidate lines per region if \code{sel} is specified.
#'   \item \code{regions}: Data frame of normalized genomic target coordinates evaluated.
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#' 
#' # 1. Fetch data and families
#' pcil_data <- fetch_pcil_data(connect_db_mode = 'online')
#' selection <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")
#' 
#' results <- fetch_pcil_families_by_variant(
#'   selection = selection,
#'   pcil_data = pcil_data,
#'   connect_db_mode = 'online'
#' )
#' 
#' # 2. Get genotypes for the selected variants
#' variant_geno_sel <- fetch_genotypes_by_id(
#'   variant_ids = selection, 
#'   connect_db_mode = 'online'
#' )
#' 
#' # 3. Select PCIL positives by variant positions
#' pcil_pos_pcv <- fetch_pcil_positive(
#'   pcil_data = pcil_data,
#'   variants_select_geno = variant_geno_sel,
#'   type = "position",
#'   sel = 15,
#'   available_ids = results$pcil_summary[c("sample_id", "selection")],
#'   result_pcil_families = results,
#'   window = 0
#' )
#' }
fetch_pcil_positive <- function(
  pcil_data,
  variants_select_geno,
  type = c("gene", "position"),
  sel = NULL,
  donor_thresh = 0.75,
  block_quantile = 0.75,
  F_quantile = 0.25,
  window = 0,
  available_ids = NULL,
  result_pcil_families = NULL
) {
  type <- match.arg(type)

  if (!is.null(available_ids)) {
    required_available <- c("sample_id", "selection")
    if (!all(required_available %in% colnames(available_ids))) {
      stop(
        "available_ids must contain columns: 'sample_id', 'selection'",
        call. = FALSE
      )
    }
  }

  global_available_ids <- NULL
  if (!is.null(result_pcil_families)) {
    global_available_ids <- result_pcil_families$pcil_summary$sample_id[duplicated(
      result_pcil_families$pcil_summary$sample_id
    )]
    if (
      !is.vector(global_available_ids) && !is.character(global_available_ids)
    ) {
      stop(
        "result_pcil_families must be a list object from fetch_pcil_families_by_variant",
        call. = FALSE
      )
    }
  }

  introgressions <- pcil_data$pcil_introgressions
  genomewide_introgressions <- pcil_data$pcil_genomewide_introgressions
  gene_regions <- pcil_data$pcil_gene_regions
  inbreeding_coefficient <- pcil_data$pcil_inbreeding_coefficient

  if (type == "position") {
    variants_select_geno <- variants_select_geno[, c(
      "variant_id",
      "chrom",
      "pos"
    )]
    names(variants_select_geno) <- c("Region", "Chr", "pos")
    required <- c("Region", "Chr", "pos")
    if (!all(required %in% colnames(variants_select_geno))) {
      stop(
        "For type='position', variants_select_geno must contain: Region, Chr, pos",
        call. = FALSE
      )
    }
    message(paste0("Using +/- ", window, " bp window around positions"))
    regions <- variants_select_geno
    regions$Start <- regions$pos - window
    regions$End <- regions$pos + window
    regions <- regions[, c("Region", "Chr", "Start", "End")]
  }

  if (type == "gene") {
    if (is.character(variants_select_geno)) {
      variants_select_geno <- data.frame(
        Region = variants_select_geno,
        stringsAsFactors = FALSE
      )
    }
    if (!"Region" %in% colnames(variants_select_geno)) {
      stop(
        "Input must contain a column named 'Region' with gene IDs",
        call. = FALSE
      )
    }

    genes <- gene_regions[
      gene_regions$GeneID %in% variants_select_geno$Region,
    ]
    missing <- setdiff(variants_select_geno$Region, genes$GeneID)
    if (length(missing) > 0) {
      message(
        "The gene is not found in v5.1, please verify it is correct: ",
        paste(missing, collapse = ", ")
      )
    }

    regions <- data.frame(
      Region = genes$GeneID,
      Chr = genes$seqid,
      Start = genes$start,
      End = genes$end,
      stringsAsFactors = FALSE
    )
  }

  introgressions_in_region <- data.frame()

  for (r in seq_len(nrow(regions))) {
    tmp <- regions[r, ]
    active_introgressions <- introgressions

    if (!is.null(global_available_ids)) {
      active_introgressions <- active_introgressions[
        active_introgressions$SampleID %in% global_available_ids,
      ]
      if (nrow(active_introgressions) == 0) {
        stop(
          "None of the Sample IDs from fetch_pcil_families_by_variant are present in the introgressions table",
          call. = FALSE
        )
      }
    }

    if (!is.null(available_ids)) {
      available_ids_tofilter <- unique(available_ids$sample_id[
        available_ids$selection == tmp$Region
      ])
      if (length(available_ids_tofilter) == 0) {
        warning(
          paste("No available IDs provided for:", tmp$Region),
          call. = FALSE
        )
        next
      }
      active_introgressions <- active_introgressions[
        active_introgressions$SampleID %in% available_ids_tofilter,
      ]
    }

    tmp2 <- active_introgressions[active_introgressions$ChrLabel == tmp$Chr, ]
    tmp3 <- tmp2[
      tmp2$block_start_bp <= tmp$Start & tmp2$block_end_bp >= tmp$End,
    ]

    if (nrow(tmp3) == 0) {
      warning(paste("No PCIL (+) found for:", tmp$Region), call. = FALSE)
      next
    }

    tmp3$Region <- tmp$Region
    tmp3$Chr <- tmp$Chr
    tmp3$Start <- tmp$Start
    tmp3$End <- tmp$End
    tmp3 <- tmp3[, c(
      "SampleID",
      "Region",
      "Chr",
      "Start",
      "End",
      "ChrLabel",
      "block_start_bp",
      "block_end_bp",
      "block_len_Mb",
      "mean_donor_frac"
    )]
    introgressions_in_region <- rbind(introgressions_in_region, tmp3)
  }

  if (nrow(introgressions_in_region) == 0) {
    warning("No PCIL (+) found for any requested region", call. = FALSE)
    return(list(pcil_positive = introgressions_in_region, regions = regions))
  }

  introgressions_in_region <- merge(
    introgressions_in_region,
    genomewide_introgressions,
    by = "SampleID",
    all.x = TRUE
  )

  if (!is.null(sel)) {
    df_best_by_region <- data.frame()

    for (g in unique(introgressions_in_region$Region)) {
      tmp2 <- introgressions_in_region[introgressions_in_region$Region == g, ]
      tmp2 <- tmp2[tmp2$mean_donor_frac >= donor_thresh, ]
      if (nrow(tmp2) == 0) {
        next
      }

      if (nrow(tmp2) > 5) {
        threshold_block <- stats::quantile(
          tmp2$block_len_Mb,
          block_quantile,
          na.rm = TRUE
        )
        tmp2 <- tmp2[tmp2$block_len_Mb <= threshold_block, ]
      }

      inbreed_id_col <- if ("SampleID" %in% names(inbreeding_coefficient)) {
        "SampleID"
      } else {
        "sample_id"
      }
      tmp2 <- merge(
        tmp2,
        inbreeding_coefficient,
        by.x = "SampleID",
        by.y = inbreed_id_col,
        all.x = TRUE
      )

      if (nrow(tmp2) > 5) {
        threshold_F <- stats::quantile(tmp2$F, F_quantile, na.rm = TRUE)
        tmp2 <- tmp2[tmp2$F >= threshold_F, ]
      }

      tmp3 <- tmp2[order(tmp2$total_Mb, tmp2$total_blocks, -tmp2$F), ]
      tmp4 <- tmp3[seq_len(min(sel, nrow(tmp3))), ]
      tmp4$Rank <- seq_len(nrow(tmp4))
      df_best_by_region <- rbind(df_best_by_region, tmp4)
    }

    return(list(
      pcil_positive = introgressions_in_region,
      best_lines = df_best_by_region,
      regions = regions
    ))
  }

  return(list(pcil_positive = introgressions_in_region, regions = regions))
}




#' Fetch Negative Control PCILs
#'
#' @description
#' Identifies and ranks negative control Pangenome-Characterized Introgression
#' Lines (PCILs) for a given set of positive lines. A negative control is a
#' line that is genetically similar (low Identity-by-State distance) to a
#' positive line but lacks the specific introgression of interest.
#'
#' @param pcil_data A list object containing PCIL data tables, typically from `fetch_pcil_data()`.
#' @param pcil_positive_result A list object returned by `fetch_pcil_positive()`, which
#'   must contain the `best_lines` and `regions` data frames.
#' @param n_neg Integer. The number of top-ranked negative control lines to
#'   return for each positive line. Defaults to `10`.
#' @param available_ids An optional data frame with `sample_id` and `selection`
#'   columns to restrict the search space for negative controls.
#' @param result_pcil_families An optional list object from `fetch_pcil_families_by_variant()`
#'   to globally restrict the search space.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{pairs_best}: A data frame with the single best negative control
#'     for each positive line.
#'   \item \code{pairs_extended}: A data frame with up to `n_neg` ranked negative
#'     controls for each positive line.
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#' 
#' # 1. Fetch data and establish positive lines
#' pcil_data <- fetch_pcil_data(connect_db_mode = "online")
#' selection <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")
#' 
#' variant_geno_sel <- fetch_genotypes_by_id(variant_ids = selection, connect_db_mode = "online")
#' fam_results <- fetch_pcil_families_by_variant(
#'   selection = selection,
#'   pcil_data = pcil_data,
#'   connect_db_mode = "online"
#' )
#'
#' pcil_pos_pcv <- fetch_pcil_positive(
#'   pcil_data = pcil_data,
#'   variants_select_geno = variant_geno_sel,
#'   type = "position",
#'   sel = 15,
#'   available_ids = fam_results$pcil_summary[, c("sample_id", "selection")],
#'   result_pcil_families = fam_results,
#'   window = 0
#' )
#'
#' # 2. Fetch the top 10 negative control PCILs for the positive lines
#' pcil_neg_pcv <- fetch_pcil_negative(
#'   pcil_data = pcil_data,
#'   pcil_positive_result = pcil_pos_pcv,
#'   n_neg = 10
#' ) 
#' 
#' print(head(pcil_neg_pcv$pairs_best))
#' }
fetch_pcil_negative <- function(
  pcil_data,
  pcil_positive_result,
  n_neg = 10,
  available_ids = NULL,
  result_pcil_families = NULL
) {
  # Extract tables from pre-loaded PCIL data 
  if (missing(pcil_data) || !is.list(pcil_data)) {
    stop(
      "`pcil_data` must be a list object from fetch_pcil_data().",
      call. = FALSE
    )
  }

  required_tables <- c(
    "pcil_metadata",
    "pcil_introgressions",
    "pcil_IBS_dis",
    "pcil_genomewide_introgressions",
    "pcil_inbreeding_coefficient"
  )
  if (!all(required_tables %in% names(pcil_data))) {
    stop("`pcil_data` is missing required tables.", call. = FALSE)
  }

  pcil_metadata <- pcil_data$pcil_metadata
  introgressions <- pcil_data$pcil_introgressions
  pcil_IBS_dis <- pcil_data$pcil_IBS_dis
  pcil_genomewide_introgressions <- pcil_data$pcil_genomewide_introgressions
  pcil_inbreeding_coefficient <- pcil_data$pcil_inbreeding_coefficient

  # --- 2. Input Validation ---
  required_results <- c("best_lines", "regions")
  if (
    !is.list(pcil_positive_result) ||
      !all(required_results %in% names(pcil_positive_result))
  ) {
    stop(
      "`pcil_positive_result` must be a list containing 'best_lines' and 'regions'.",
      call. = FALSE
    )
  }

  if (
    is.null(pcil_positive_result$best_lines) ||
      nrow(pcil_positive_result$best_lines) == 0
  ) {
    warning("No 'best_lines' found in pcil_positive_result.", call. = FALSE)
    return(list(pairs_best = data.frame(), pairs_extended = data.frame()))
  }

  # Normalize SampleID column names for the metadata file
  id_col_meta <- if ("SampleID" %in% colnames(pcil_metadata)) {
    "SampleID"
  } else {
    "sample_id"
  }
  fam_col_meta <- if ("Family" %in% colnames(pcil_metadata)) {
    "Family"
  } else {
    "family"
  }

  all_pcil_pairs_extended <- list()
  all_sample_ids <- unique(pcil_metadata[[id_col_meta]])

  # Looping through each region and positive line 
  for (region_i in unique(pcil_positive_result$regions$Region)) {
    positive_lines_in_region <- pcil_positive_result$best_lines[
      pcil_positive_result$best_lines$Region == region_i,
    ]
    if (nrow(positive_lines_in_region) == 0) {
      next
    }

    # Define the pool of candidate negative lines (ones that DON'T overlap the introgression)
    lines_with_introgression <- unique(introgressions$SampleID[
      introgressions$ChrLabel == positive_lines_in_region$Chr[1] &
        introgressions$block_start_bp <= positive_lines_in_region$Start[1] &
        introgressions$block_end_bp >= positive_lines_in_region$End[1]
    ])

    candidate_pool <- setdiff(all_sample_ids, lines_with_introgression)

    # Global constraint logic
    if (!is.null(result_pcil_families)) {
      global_available_ids <- result_pcil_families$pcil_summary$sample_id[
        duplicated(result_pcil_families$pcil_summary$sample_id)
      ]
      if (length(global_available_ids) > 0) {
        candidate_pool <- intersect(candidate_pool, global_available_ids)
      }
    }

    if (!is.null(available_ids)) {
      available_for_region <- unique(available_ids$sample_id[
        available_ids$selection == region_i
      ])
      candidate_pool <- intersect(candidate_pool, available_for_region)
    }

    if (length(candidate_pool) == 0) {
      warning(
        paste("No negative candidates available for region:", region_i),
        call. = FALSE
      )
      next
    }

    # Analyze each specific positive line
    for (l in seq_len(nrow(positive_lines_in_region))) {
      focal_line <- positive_lines_in_region[l, ]
      focal_id <- focal_line$SampleID

      # Attempt to restrict to the same family
      focal_family <- pcil_metadata[[fam_col_meta]][
        pcil_metadata[[id_col_meta]] == focal_id
      ]
      if (length(focal_family) > 0) {
        focal_family <- focal_family[1]
      } else {
        focal_family <- NA
      }

      family_members <- pcil_metadata[[id_col_meta]][
        pcil_metadata[[fam_col_meta]] == focal_family &
          !is.na(pcil_metadata[[fam_col_meta]])
      ]
      family_candidates <- intersect(candidate_pool, family_members)

      current_candidates <- if (length(family_candidates) > 0) {
        family_candidates
      } else {
        candidate_pool
      }

      #  IBS EXTRACTION (Handles both Long & Wide formats) 
      is_long_format <- all(
        c("sample_1", "sample_2", "distance") %in% colnames(pcil_IBS_dis)
      )

      if (is_long_format) {
        # Long format extraction
        ibs_subset <- pcil_IBS_dis[
          pcil_IBS_dis$sample_1 == focal_id | pcil_IBS_dis$sample_2 == focal_id,
        ]

        if (nrow(ibs_subset) == 0) {
          warning(
            paste("No IBS data found for positive line:", focal_id),
            call. = FALSE
          )
          next
        }

        # Identify which column contains the candidate (the one that isn't the focal_id)
        candidates_found <- ifelse(
          ibs_subset$sample_1 == focal_id,
          ibs_subset$sample_2,
          ibs_subset$sample_1
        )
        valid_idx <- candidates_found %in% current_candidates

        if (!any(valid_idx)) {
          warning(
            paste("No overlapping candidates found in IBS data for:", focal_id),
            call. = FALSE
          )
          next
        }

        ibs_distances <- ibs_subset$distance[valid_idx]
        names(ibs_distances) <- candidates_found[valid_idx]
      } else {
        # Wide format extraction 
        id_col_ibs <- if ("SampleID" %in% colnames(pcil_IBS_dis)) {
          "SampleID"
        } else {
          "sample_id"
        }
        ibs_row <- pcil_IBS_dis[pcil_IBS_dis[[id_col_ibs]] == focal_id, ]

        if (nrow(ibs_row) == 0) {
          warning(
            paste("No IBS data found for positive line:", focal_id),
            call. = FALSE
          )
          next
        }

        valid_cand_cols <- intersect(names(ibs_row), current_candidates)
        if (length(valid_cand_cols) == 0) {
          warning(
            paste("No overlapping candidates found in IBS data for:", focal_id),
            call. = FALSE
          )
          next
        }

        ibs_candidates <- ibs_row[, valid_cand_cols, drop = FALSE]
        ibs_distances <- as.numeric(ibs_candidates[1, ])
        names(ibs_distances) <- valid_cand_cols
      }

      # Clean and Sort Distances
      ibs_distances <- ibs_distances[!is.na(ibs_distances)]
      if (length(ibs_distances) == 0) {
        warning(paste("No valid IBS values for:", focal_id), call. = FALSE)
        next
      }

      ibs_distances <- sort(ibs_distances)

      # Determine number of IBS candidates to evaluate based on n_neg
      target_n <- if (is.null(n_neg) || n_neg == 1) 3 else n_neg
      n_ibs <- min(target_n, length(ibs_distances))
      top_ibs_ids <- names(ibs_distances)[1:n_ibs]

      # Filter and rank by Genome-wide and Inbreeding metrics
      ranked_candidates <- pcil_genomewide_introgressions[
        pcil_genomewide_introgressions$SampleID %in% top_ibs_ids,
      ]

      inbreed_id_col <- if (
        "SampleID" %in% names(pcil_inbreeding_coefficient)
      ) {
        "SampleID"
      } else {
        "sample_id"
      }

      ranked_candidates <- merge(
        ranked_candidates,
        pcil_inbreeding_coefficient,
        by.x = "SampleID",
        by.y = inbreed_id_col,
        all.x = TRUE
      )

      # Sort: Smallest Mb -> Fewest Blocks -> Highest F
      ranked_candidates <- ranked_candidates[
        order(
          ranked_candidates$total_Mb,
          ranked_candidates$total_blocks,
          -ranked_candidates$F
        ),
      ]

      n_keep <- if (is.null(n_neg) || n_neg == 1) {
        1
      } else {
        min(n_neg, nrow(ranked_candidates))
      }
      best_negatives <- ranked_candidates[seq_len(n_keep), , drop = FALSE]

      if (nrow(best_negatives) == 0) {
        next
      }

      pairs_df <- data.frame(
        SampleID_Positive = rep(focal_id, nrow(best_negatives)),
        SampleID_Negative = best_negatives$SampleID,
        Region = rep(region_i, nrow(best_negatives)),
        Chr = rep(focal_line$Chr, nrow(best_negatives)),
        Start = rep(focal_line$Start, nrow(best_negatives)),
        End = rep(focal_line$End, nrow(best_negatives)),
        IBS_dis = ibs_distances[best_negatives$SampleID],
        total_Mb_neg = best_negatives$total_Mb,
        total_blocks_neg = best_negatives$total_blocks,
        F_neg = best_negatives$F,
        Family = rep(focal_family, nrow(best_negatives)),
        Pair_Rank = seq_len(nrow(best_negatives)),
        stringsAsFactors = FALSE
      )

      all_pcil_pairs_extended[[length(all_pcil_pairs_extended) + 1]] <- pairs_df
    }
  }

  if (length(all_pcil_pairs_extended) == 0) {
    return(list(pairs_best = data.frame(), pairs_extended = data.frame()))
  }

  pairs_extended <- do.call(rbind, all_pcil_pairs_extended)
  pairs_best <- pairs_extended[pairs_extended$Pair_Rank == 1, ]

  if (is.null(n_neg) || n_neg == 1) {
    return(list(pairs_best = pairs_best))
  } else {
    return(list(pairs_best = pairs_best, pairs_extended = pairs_extended))
  }
}



#' Plot PCIL Positive Introgressions
#'
#' @description
#' Visualizes the genomic regions of Pangenome-Characterized Introgression Lines (PCIL)
#' that were identified as positive for a given target. It generates a plot for each
#' target region, showing the introgression blocks for all positive lines.
#'
#' @param pcil_positive_result A list object returned by `fetch_pcil_positive`,
#'   containing `pcil_positive` and `regions` data frames.
#'
#' @return A list of ggplot objects, where each plot corresponds to a unique
#'   target region and visualizes the introgression segments of the positive PCILs.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_segment geom_vline labs theme_bw theme element_blank element_text
#'
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#' 
#' # 1. Fetch PCIL data and variant genotypes
#' pcil_data <- fetch_pcil_data(connect_db_mode = "online")
#' selection <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")
#' 
#' variant_geno_sel <- fetch_genotypes_by_id(
#'   variant_ids = selection, 
#'   connect_db_mode = "online"
#' )
#' 
#' # 2. Fetch relevant families
#' fam_results <- fetch_pcil_families_by_variant(
#'   selection = selection,
#'   pcil_data = pcil_data,
#'   connect_db_mode = "online"
#' )
#' 
#' # 3. Select positive lines
#' pcil_pos_pcv <- fetch_pcil_positive(
#'   pcil_data = pcil_data,
#'   variants_select_geno = variant_geno_sel,
#'   type = "position",
#'   sel = 15,
#'   available_ids = fam_results$pcil_summary[, c("sample_id", "selection")],
#'   result_pcil_families = fam_results,
#'   window = 0
#' )
#' 
#' # 4. Generate and display the plot
#' positive_plots <- plot_pcil_positive(pcil_positive_result = pcil_pos_pcv)
#' 
#' if (length(positive_plots) > 0) {
#'   print(positive_plots[[1]])
#' }
#' }
plot_pcil_positive <- function(pcil_positive_result) {
  Region <- block_start_bp <- block_end_bp <- SampleID <- NULL

  if (!is.list(pcil_positive_result) || !all(c("pcil_positive", "regions") %in% names(pcil_positive_result))) {
    stop("Input must be a list object returned from fetch_pcil_positive().", call. = FALSE)
  }

  pcil_pos_pcv <- pcil_positive_result
  all_pcil_pos_plots <- list()

  if (nrow(pcil_pos_pcv$pcil_positive) == 0) {
    warning("No positive PCILs to plot.", call. = FALSE)
    return(all_pcil_pos_plots)
  }

  for (region_i in unique(pcil_pos_pcv$regions$Region)) {
    region_tmp <- pcil_pos_pcv$regions[pcil_pos_pcv$regions$Region == region_i, ]

    # Filter for the current region and sort the data for plotting
    df <- pcil_pos_pcv$pcil_positive[pcil_pos_pcv$pcil_positive$Region == region_i, ]
    df <- df[order(df$block_start_bp, df$block_end_bp), ]

    # To preserve the order in the plot, convert SampleID to a factor
    # with levels in the order of appearance.
    df$SampleID <- factor(df$SampleID, levels = unique(df$SampleID))

    if (nrow(df) == 0) next

    p <- ggplot2::ggplot(df, ggplot2::aes(x = block_start_bp / 1e6, xend = block_end_bp / 1e6, y = SampleID, yend = SampleID)) +
      ggplot2::geom_segment(linewidth = 1, lineend = "round", color = "grey35") +
      ggplot2::geom_vline(xintercept = c(region_tmp$Start / 1e6, region_tmp$End / 1e6), color = "red", linewidth = 1) +
      ggplot2::labs(x = "Chromosome position (Mb)", y = "Samples with Introgressions", title = paste("All PCIL (+):", region_i)) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5))

    all_pcil_pos_plots[[region_tmp$Region]] <- p
  }

  return(all_pcil_pos_plots)
}



#' Plot Best PCIL Candidate Lines
#'
#' @description
#' Creates a genome-wide overview plot for each target region, visualizing the
#' introgression patterns of the top-ranked Pangenome-Characterized Introgression
#' Lines (PCILs) identified by `fetch_pcil_positive()`.
#'
#' @param pcil_positive_result A list object returned by `fetch_pcil_positive`,
#'   which must contain the `best_lines` and `regions` data frames.
#' @param pcil_data A list object containing all required PCIL data tables,
#'   typically from `fetch_pcil_data()`. Must contain `pcil_introgressions`.
#'
#' @return A list of ggplot objects. Each plot corresponds to a unique target
#'   region and displays the introgression segments for the best candidate lines
#'   across all chromosomes.
#'
#' @export
#' @import DBI
#' @importFrom ggplot2 ggplot aes geom_blank geom_hline geom_segment scale_color_manual geom_vline facet_wrap labs theme_minimal theme element_blank element_text element_line element_rect unit margin
#'
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#' 
#' # 1. Fetch data and find positive lines
#' pcil_data <- fetch_pcil_data(connect_db_mode = "online")
#' selection <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")
#' 
#' variant_geno_sel <- fetch_genotypes_by_id(variant_ids = selection, connect_db_mode = "online")
#' fam_results <- fetch_pcil_families_by_variant(
#'   selection = selection,
#'   pcil_data = pcil_data,
#'   connect_db_mode = "online"
#' )
#'
#' pcil_pos_pcv <- fetch_pcil_positive(
#'   pcil_data = pcil_data,
#'   variants_select_geno = variant_geno_sel,
#'   type = "position",
#'   sel = 15,
#'   available_ids = fam_results$pcil_summary[, c("sample_id", "selection")],
#'   result_pcil_families = fam_results,
#'   window = 0
#' )
#'
#' # 2. Plot the best candidate lines genome-wide
#' best_line_plots <- plot_pcil_best_lines(
#'   pcil_positive_result = pcil_pos_pcv, 
#'   pcil_data = pcil_data
#' )
#' 
#' if (length(best_line_plots) > 0) {
#'   print(best_line_plots[[1]])
#' }
#' }
plot_pcil_best_lines <- function(
  pcil_positive_result,
  pcil_data
) {
  # Define global variables to avoid CMD check notes
  Region <- SampleID <- Rank <- ChrLabel <- chr_endMb <- startMb <- endMb <- NULL
  highlight <- plot_id <- StartMb <- EndMb <- NULL

  if (
    !is.list(pcil_positive_result) ||
      !all(c("best_lines", "regions") %in% names(pcil_positive_result))
  ) {
    stop(
      "`pcil_positive_result` must be a list object from fetch_pcil_positive() containing 'best_lines' and 'regions'.",
      call. = FALSE
    )
  }
  if (
    is.null(pcil_positive_result$best_lines) || nrow(pcil_positive_result$best_lines) == 0
  ) {
    warning("No 'best_lines' found in pcil_positive_result to plot.", call. = FALSE)
    return(list())
  }

  all_pcil_best_plots <- list()

  # Fetch Introgressions Data from pcil_data 
  if (missing(pcil_data) || !is.list(pcil_data) || !"pcil_introgressions" %in% names(pcil_data)) {
    stop("`pcil_data` must be a list from fetch_pcil_data() containing 'pcil_introgressions'.", call. = FALSE)
  }
  introgressions <- pcil_data$pcil_introgressions

  chr_levels <- sprintf("Chr%02d", 1:10)
  chr_lengths <- data.frame(
    ChrLabel = factor(chr_levels, levels = chr_levels),
    chr_endMb = c(74, 78, 74, 68, 62, 62, 64, 55, 60, 61)
  )

  #  Loop Through Each Region
  for (region_i in unique(pcil_positive_result$regions$Region)) {
    best_lines_region <- pcil_positive_result$best_lines[
      pcil_positive_result$best_lines$Region == region_i,
    ]
    if (nrow(best_lines_region) == 0) {
      next
    }

    best_for_plot <- introgressions[
      introgressions$SampleID %in% best_lines_region$SampleID,
    ]
    if (nrow(best_for_plot) == 0) {
      next
    }

    best_for_plot$plot_id <- paste0(
      best_for_plot$SampleID,
      " (",
      best_for_plot$Clan,
      "/",
      best_for_plot$Family,
      ")"
    )

    sample_order <- best_lines_region$SampleID[order(best_lines_region$Rank)]

    plot_order_df <- unique(best_for_plot[, c("SampleID", "plot_id")])
    plot_order_df$rank <- match(plot_order_df$SampleID, sample_order)
    plot_order_df <- plot_order_df[order(plot_order_df$rank), ]
    plot_order <- plot_order_df$plot_id

    best_for_plot$plot_id <- factor(
      best_for_plot$plot_id,
      levels = rev(plot_order)
    )

    top_sample <- sample_order[1]
    best_for_plot$highlight <- ifelse(
      best_for_plot$SampleID == top_sample,
      "Top-Ranked",
      "Other Candidates"
    )

    # Setting plot framework
    facet_background <- expand.grid(
      ChrLabel = factor(chr_levels, levels = chr_levels),
      plot_id = levels(best_for_plot$plot_id)
    )
    facet_background <- merge(facet_background, chr_lengths, by = "ChrLabel")
    facet_background$plot_id <- factor(
      facet_background$plot_id,
      levels = levels(best_for_plot$plot_id)
    )
    best_for_plot$ChrLabel <- factor(
      best_for_plot$ChrLabel,
      levels = chr_levels
    )

    region_tmp <- pcil_positive_result$regions[
      pcil_positive_result$regions$Region == region_i,
    ]
    region_df <- data.frame(
      ChrLabel = factor(region_tmp$Chr, levels = chr_levels),
      StartMb = region_tmp$Start / 1e6,
      EndMb = region_tmp$End / 1e6
    )

    # Generate the plot
    p <- ggplot2::ggplot() +
      ggplot2::geom_blank(
        data = facet_background,
        ggplot2::aes(x = 0, y = plot_id)
      ) +
      ggplot2::geom_blank(
        data = facet_background,
        ggplot2::aes(x = chr_endMb, y = plot_id)
      ) +
      ggplot2::geom_hline(
        yintercept = seq_along(levels(best_for_plot$plot_id)),
        color = "#e9ecef",
        linewidth = 4,
        lineend = "square"
      ) +
      ggplot2::geom_segment(
        data = best_for_plot,
        ggplot2::aes(
          x = startMb,
          xend = endMb,
          y = plot_id,
          yend = plot_id,
          color = highlight
        ),
        linewidth = 4,
        lineend = "butt"
      ) +
      ggplot2::scale_color_manual(
        values = c("Top-Ranked" = "#0d6efd", "Other Candidates" = "#343a40"),
        guide = "none"
      ) +
      ggplot2::geom_vline(
        data = region_df,
        ggplot2::aes(xintercept = StartMb),
        color = "#dc3545",
        linewidth = 0.8,
        linetype = "dashed",
        alpha = 0.8
      ) +
      ggplot2::geom_vline(
        data = region_df,
        ggplot2::aes(xintercept = EndMb),
        color = "#dc3545",
        linewidth = 0.8,
        linetype = "dashed",
        alpha = 0.8
      ) +
      ggplot2::facet_wrap(
        ~ChrLabel,
        ncol = 5,
        scales = "free_x",
        drop = FALSE
      ) +
      ggplot2::labs(
        title = paste("Best PCIL Candidates for:", region_i),
        subtitle = "Genome-wide introgression patterns for top-ranked lines.",
        x = "Chromosome Position (Mb)",
        y = "PCIL Sample"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.spacing = ggplot2::unit(1.5, "lines"),
        plot.title = ggplot2::element_text(
          size = 16,
          face = "bold",
          hjust = 0.5,
          color = "#212529"
        ),
        plot.subtitle = ggplot2::element_text(
          size = 12,
          hjust = 0.5,
          color = "#6c757d",
          margin = ggplot2::margin(b = 15)
        ),
        axis.line = ggplot2::element_line(linewidth = 0.5, color = "#adb5bd"),
        axis.text = ggplot2::element_text(color = "#495057"),
        axis.title = ggplot2::element_text(face = "bold", color = "#495057"),
        strip.background = ggplot2::element_rect(
          fill = "#f8f9fa",
          color = "#dee2e6"
        ),
        strip.text = ggplot2::element_text(
          face = "bold",
          size = 11,
          color = "#212529"
        ),
        legend.position = "none"
      )

    all_pcil_best_plots[[region_i]] <- p
  }

  return(all_pcil_best_plots)
}

 

 
#' Plot PCIL Positive and Negative Pairs
#'
#' @description
#' Creates a genome-wide overview plot for each positive-negative PCIL pair
#' identified by `fetch_pcil_negative`. Each plot visualizes the introgression
#' patterns for a positive line and its corresponding negative controls.
#'
#' @param pcil_neg_results A list object returned by `fetch_pcil_negative`,
#'   which must contain the `pairs_extended` data frame.
#' @param pcil_data A list object returned by `fetch_pcil_data`, containing
#'   all necessary PCIL data tables.
#'
#' @return A list of ggplot objects. Each plot corresponds to a unique
#'   positive PCIL and displays its introgression pattern alongside its ranked
#'   negative control pairs.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_blank geom_hline geom_segment scale_color_manual geom_vline facet_wrap labs theme_classic theme element_blank element_text element_line unit
#'
#' @examples
#' \donttest{
#' library(panGenomeBreedr)
#' 
#' # 1. Run the full pipeline to identify positive and negative pairs
#' pcil_data <- fetch_pcil_data(connect_db_mode = "online")
#' selection <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")
#' 
#' variant_geno_sel <- fetch_genotypes_by_id(variant_ids = selection, connect_db_mode = "online")
#' fam_results <- fetch_pcil_families_by_variant(
#'   selection = selection,
#'   pcil_data = pcil_data,
#'   connect_db_mode = "online"
#' )
#'
#' pcil_pos_pcv <- fetch_pcil_positive(
#'   pcil_data = pcil_data,
#'   variants_select_geno = variant_geno_sel,
#'   type = "position",
#'   sel = 15,
#'   available_ids = fam_results$pcil_summary[, c("sample_id", "selection")],
#'   result_pcil_families = fam_results,
#'   window = 0
#' )
#'
#' pcil_neg_pcv <- fetch_pcil_negative(
#'   pcil_data = pcil_data,
#'   pcil_positive_result = pcil_pos_pcv,
#'   n_neg = 10
#' )
#' 
#' # 2. Generate side-by-side visual comparisons
#' pair_plots <- plot_pcil_negative_pairs(
#'   pcil_neg_results = pcil_neg_pcv, 
#'   pcil_data = pcil_data
#' )
#' 
#' if (length(pair_plots) > 0) {
#'   print(pair_plots[[1]])
#' }
#' }
plot_pcil_negative_pairs <- function(pcil_neg_results, pcil_data) {

  # --- 1. Input Validation ---
  if (!is.list(pcil_neg_results) || !"pairs_extended" %in% names(pcil_neg_results)) {
    stop("`pcil_neg_results` must be a list from fetch_pcil_negative() containing 'pairs_extended'.", call. = FALSE)
  }
  if (is.null(pcil_neg_results$pairs_extended) || nrow(pcil_neg_results$pairs_extended) == 0) {
    warning("No PCIL pairs found in `pcil_neg_results` to plot.", call. = FALSE, immediate. = TRUE)
    return(list())
  }
  required_data <- c("pcil_introgressions", "pcil_metadata")
  if (!is.list(pcil_data) || !all(required_data %in% names(pcil_data))) {
    stop("`pcil_data` must be a list from fetch_pcil_data() containing 'pcil_introgressions' and 'pcil_metadata'.", call. = FALSE)
  }

  all_pcil_pair_plots <- list()

  chr_levels <- sprintf("Chr%02d", 1:10)
  chr_lengths <- data.frame(ChrLabel = factor(chr_levels, levels = chr_levels), chr_endMb = c(74, 78, 74, 68, 62, 62, 64, 55, 60, 61))

  #  Loop through each region and positive line ---
  for (selection_i in unique(pcil_neg_results$pairs_extended$Region)) {
    pairs_variant <- pcil_neg_results$pairs_extended[pcil_neg_results$pairs_extended$Region == selection_i, ]

    for (positive_i in unique(pairs_variant$SampleID_Positive)) {
      pairs_tmp <- pairs_variant[pairs_variant$SampleID_Positive == positive_i, ]
      samples <- unique(c(positive_i, pairs_tmp$SampleID_Negative))

      # Get introgressions and merge with metadata to get Clan/Family
      best_for_plot <- pcil_data$pcil_introgressions[pcil_data$pcil_introgressions$SampleID %in% samples, ]
  best_for_plot <- merge(best_for_plot, pcil_data$pcil_metadata[, c("sample_id", "clan", "family")], by.x = "SampleID", by.y = "sample_id", all.x = TRUE)

  best_for_plot$plot_id <- paste0(best_for_plot$SampleID, " (", best_for_plot$clan, "/", best_for_plot$family, ")")

      # Order samples for plotting
      sample_order <- c(positive_i, pairs_tmp$SampleID_Negative[order(pairs_tmp$Pair_Rank)])
      plot_order_df <- unique(best_for_plot[, c("SampleID", "plot_id")])
  plot_order_df <- best_for_plot[!duplicated(best_for_plot$SampleID), c("SampleID", "plot_id")]
      plot_order_df$rank <- match(plot_order_df$SampleID, sample_order)
      plot_order_df <- plot_order_df[order(plot_order_df$rank), ]
      plot_order <- plot_order_df$plot_id
      best_for_plot$plot_id <- factor(best_for_plot$plot_id, levels = rev(plot_order))

      # Highlight positive and best negative
      top_sample <- sample_order[1]
      top_sample2 <- if (length(sample_order) > 1) sample_order[2] else NA
      best_for_plot$highlight <- ifelse(
        best_for_plot$SampleID == top_sample, "PCIL (+)",
        ifelse(!is.na(top_sample2) & best_for_plot$SampleID == top_sample2, "Best PCIL (-) pair", "Other PCIL (-) pairs")
      )

      # Create empty chromosome background
      facet_background <- expand.grid(ChrLabel = factor(chr_levels, levels = chr_levels), plot_id = levels(best_for_plot$plot_id))
      facet_background <- merge(facet_background, chr_lengths, by = "ChrLabel")
      facet_background$plot_id <- factor(facet_background$plot_id, levels = levels(best_for_plot$plot_id))
      best_for_plot$ChrLabel <- factor(best_for_plot$ChrLabel, levels = chr_levels)

      # Region of interest
      region_info <- pairs_tmp[1, ]
      region_df <- data.frame(
        ChrLabel = factor(region_info$Chr, levels = chr_levels),
        StartMb = region_info$Start / 1e6,
        EndMb = region_info$End / 1e6
      )

      # Generate Plot
      p <- ggplot2::ggplot() +
        ggplot2::geom_blank(data = facet_background, ggplot2::aes(x = 0, y = .data$plot_id)) +
        ggplot2::geom_blank(data = facet_background, ggplot2::aes(x = .data$chr_endMb, y = .data$plot_id)) +
        ggplot2::geom_hline(yintercept = seq_along(levels(best_for_plot$plot_id)), color = "grey93", linewidth = 4, lineend = "square") +
        ggplot2::geom_segment(
          data = best_for_plot,
          ggplot2::aes(x = .data$startMb, xend = .data$endMb, y = .data$plot_id, yend = .data$plot_id, color = .data$highlight),
          linewidth = 4, lineend = "butt"
        ) +
        ggplot2::scale_color_manual(
          values = c("PCIL (+)" = "#C51B8A", "Best PCIL (-) pair" = "blue4", "Other PCIL (-) pairs" = "black"),
          guide = "none"
        ) +
        ggplot2::geom_vline(data = region_df, ggplot2::aes(xintercept = .data$StartMb), color = "red", linewidth = 1) +
        ggplot2::geom_vline(data = region_df, ggplot2::aes(xintercept = .data$EndMb), color = "red", linewidth = 1) +
        ggplot2::facet_wrap(~ChrLabel, ncol = 5, scales = "free_x", drop = FALSE) +
        ggplot2::labs(
          x = "Chromosome position (Mb)",
          y = "Sample",
          title = paste("PCIL (+) and PCIL (-) pairs for", selection_i)
        ) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(linewidth = 0.5),
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = "bold", size = 11),
          panel.spacing = ggplot2::unit(1, "lines"),
          axis.text.y = ggplot2::element_text(size = 9),
          plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = "none"
        )

      plot_name <- paste(selection_i, positive_i, sep = "_")
      all_pcil_pair_plots[[plot_name]] <- p
    }
  }

  return(all_pcil_pair_plots)
}
