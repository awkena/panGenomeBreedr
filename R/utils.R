#' Parse marker names of Hapmap format with a common pattern containing chromosome
#' numbers and positions into a map file.
#' @param x A character vector containing the original marker names to be parsed.
#' @param sep A character value that serves as a unique separator between chromosome
#' positions and other components of the marker name; default value is an underscore.
#' @param prefix A character value that represents a common pattern in marker
#' names that precedes the chromosome number; the default value is `S`.
#'
#' @returns A data frame of map file consisting of the original marker names,
#' chromosome numbers and positions.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#'
#' snps <- paste0('S', 1:10, '_', 101:110)
#' map_file <- parse_marker_ns(x = snps, sep = '_', prefix = 'S')
#' }
#'
#' @details
#' The marker names to be parsed into a map file must contain the chromosome numbers
#' and their positions with a common separator, as well as a common pattern preceding
#' the marker names. For instance, `S1_101` and `S2_102` have `S` as the common
#' pattern preceding the marker names, `1` and `2` as the chromosome numbers, `101`
#' and `102` as the positions, and `_` as the separator.
#'
#' @export
#'
parse_marker_ns <- function(x,
                            sep = '_',
                            prefix = 'S') {

  # Throw up an error if all marker names are not of the Hapmap format
  #format_check <- paste0('^', prefix, '\\d+', sep, '\\d+$')
  format_check <- paste0('^', prefix, '\\d+', sep, '\\d+(\\.\\d+)?(e[+-]?\\d+)?$')

  if (any(!grepl(format_check, x))) stop('All marker names must be of the Hapmap format!')

  # Parse marker names to extract chromosome numbers and physical positions
  df <- t(as.data.frame(strsplit(x, sep, fixed = TRUE)))
  df <- as.data.frame(df)
  colnames(df) <- c('chr', 'pos')

  # Get chromosome numbers
  df$chr <- as.numeric(gsub(prefix, "", df$chr))

  # Convert marker positions to numeric values
  df$pos <- as.numeric(df$pos)

  # Column bind original SNP IDs with df
  df <- cbind(snpid = x, df)
  rownames(df) <- NULL

  return(df)

}

#' Create a heatmap to visualize and compare the genetic genetic backgrounds of
#' genotypes/lines.
#' @param x A numeric matrix with marker IDs as columns and sample IDs as row names.
#' @param map_file A data frame of map file consisting of SNP IDs and their
#' chromosome numbers and positions as columns.
#' @param snp_ids A character value indicating the column name for marker IDs
#' in \code{x}.
#' @param chr A character value indicating the column name for chromosome IDs
#' in \code{x}.
#' @param chr_pos A character value indicating the column name for chromosome
#' positions in \code{x}.
#' @param parents A character vector of length = 2 for the IDs of parents.
#' @param group_sz A positive integer value indicating the batch size for progenies
#' to include in heatmap.
#' @param pdf A logical value indicating whether to save plot as a pdf graphic
#' device when TRUE or output plot in R when FALSE.
#' @param filename A character value for path or file name for saving pdf.
#' @param legend_title A character value for specifying plot legend title.
#' @param col_mapping A character vector of length = 6 for heatmap color mapping.
#' @param col_labels A character vector of length = 6 for labels corresponding to
#' the color mapping.
#' @param alpha A numeric value between 0 and 1 for modifying the
#' opacity of colors.
#' @param text_size A numeric value for setting text size.
#' @param width A numeric value for the width of pdf device.
#' @param height A numeric value for the height of pdf device.
#' @param ... Other valid arguments that can be passed to ggplot2.
#'
#' @return A graphic object or ggplot.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#' # Create a numeric matrix of genotype scores for 10 markers and 5 samples
#' num_dat <- matrix(c(rep(1, 10), rep(0, 10),
#'                     1, 1, 0.5, 1, 1, 1, 1, 1, 0, 1,
#'                     1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
#'                     1, 1, 0, 1, 1, 1, 1, 1, 1, 0.5 ),
#'                   byrow = TRUE, ncol = 10)
#'
#' rownames(num_dat) <- c('rp', 'dp', paste0('bc1_', 1:3))
#' colnames(num_dat) <- paste0('S1', '_', c(floor(seq(1000, 10000, len = 8)),
#'                                          15000, 20000))
#'
#' # Get map file by parsing SNP IDs
#' map_file <- parse_marker_ns(colnames(num_dat))
#'
#' # Create a heatmap that compares the parents to progenies
#' cross_qc_ggplot(x = num_dat,
#'                 map_file = map_file,
#'                 snp_ids = 'snpid',
#'                 chr = 'chr',
#'                 chr_pos = 'pos',
#'                 value = 'value',
#'                 parents = c('rp', 'dp'),
#'                 group_sz = 3L,
#'                 pdf = FALSE,
#'                 legend_title = 'Heatmap_key',
#'                 alpha = 0.8,
#'                 text_size = 14)
#'
#' }
#'
#' @export
cross_qc_ggplot <- function(x,
                            map_file,
                            snp_ids = 'snpid',
                            chr = 'chr',
                            chr_pos = 'pos',
                            parents,
                            group_sz = nrow(x)-2,
                            pdf = FALSE,
                            filename = 'background_heatmap',
                            legend_title = 'Heatmap_key',
                            col_mapping,
                            col_labels,
                            alpha = 1,
                            text_size = 12,
                            width = 10.5,
                            height = 8.5,
                            ...){

  if (missing(x)) stop("Please provide a numeric matrix for the `x` argument.")
  if (!inherits(x, what = 'matrix')) stop("Argument `x` must be a matrix object.")
  if (missing(map_file)) stop("Please provide a map file for the `map_file` argument.")
  if (!inherits(map_file, what = 'data.frame')) stop("Argument `map_file` must be a data.frame object.")
  if (missing(parents)) stop("Please provide the value for the 'parents' argument.")
  if (length(parents) != 2) stop("The `parents` argument must be a character vector of length = 2.")

  snpid <- value  <- NULL # Define snpid as a global variable

  # Convert chr column to a character if it is numeric
  if (inherits(map_file[, chr], what = 'numeric')) {

    map_file[, chr] <- sprintf('Chr%02s', map_file[, chr])
    map_file[, chr] <- factor(map_file[, chr], levels = unique(map_file[, chr]))

  }

  # Subset parents data
  parent_dat <- rbind(x[which(rownames(x) == parents[1]),],
                      x[which(rownames(x) == parents[2]),])
  rownames(parent_dat) <- parents

  # Melt x into a data frame
  parent_dat <- gg_dat(num_mat = parent_dat,
                       map_file = map_file,
                       map_pos = chr_pos,
                       map_chr = chr,
                       map_snp_ids = snp_ids)

  # Subset progeny data
  progeny_dat <- x[!rownames(x) %in% unique(parent_dat$x),]

  # Create an index to split by
  grp_index <- rep(1:ceiling(nrow(progeny_dat) / group_sz), each = group_sz,
                   length.out = nrow(progeny_dat))

  # Split data frame into list of smaller data frames
  batches <- split(as.data.frame(progeny_dat), grp_index)

  batches <- lapply(batches, FUN = function(x) {

    progeny_dat <- gg_dat(num_mat = x,
                          map_file = map_file,
                          map_pos = chr_pos,
                          map_chr = chr,
                          map_snp_ids = snp_ids)

    grp <- rbind(parent_dat, progeny_dat)

    # Re-order levels of the sample ids
    grp$x <- factor(grp$x, levels = rev(c(sort(unique(parent_dat$x)),
                                          sort(unique(progeny_dat$x)))))
    grp$snpid <- factor(grp$snpid, levels = unique(grp$snpid))

    grp


  })

  nbatches <- length(batches) # Number of batches to plot

  # Create an empty list object to hold ggplots
  gg_plts <- vector(mode = 'list', length = nbatches)
  names(gg_plts) <- paste0('Batch', seq_len(nbatches))

  # Color and labels for plotting
  # Black = NAs; coral1 = RP; gold = Het; purple = DP; grey70 = Mono;
  # cornflowerblue = geno_error

  if (missing(col_mapping)) {

    col_mapping <- c('-5' = 'black',
                     '-2' = 'cornflowerblue',
                     '-1' = 'grey80',
                     '0' = 'purple2',
                     '0.5' = 'gold',
                     '1' = 'coral2')

  }

  if (missing(col_labels)) {

    col_labels <- c('-5' = "Missing",
                    '-2' = 'Error',
                    '-1' = "Monomorphic",
                    '0' = parents[2],
                    '0.5' = "Heterozygous",
                    '1' = parents[1])
  }

  for(i in seq_len(nbatches)) {

    grp <- batches[[i]]

    plt <- ggplot2::ggplot(grp, ggplot2::aes(x = snpid, y = x, fill = value)) +
      ggplot2::geom_tile(lwd = 2, linetype = 1) +
      ggplot2::scale_fill_manual( values = ggplot2::alpha(col_mapping, alpha),
                                  label = col_labels, name = legend_title) +
      ggplot2::labs(x = 'Chromosome',
                    title = paste("Heatmap for parents + progeny Batch", i)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = text_size),
                     axis.text.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_text(size = text_size, face = 'bold'),
                     axis.title.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = text_size, face = 'bold'),
                     legend.text = ggplot2::element_text(size = text_size),
                     legend.title = ggplot2::element_text(size = text_size, face = 'bold'),
                     axis.ticks.x = ggplot2::element_blank(),
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(size = text_size)) +
      ggplot2::geom_hline(yintercept = c(as.numeric(grp$x) + .5, .5),
                          col = 'white', lwd = 2.5) +
      ggplot2::facet_wrap(ggplot2::vars(chr), ncol = length(unique(grp$chr)),
                          nrow = 1, strip.position =  'bottom', scales = 'free_x')

    gg_plts[[i]] <- plt

  } # End of the loop

  if (pdf) {

    ggplot2::ggsave(filename = paste0(filename, ".pdf"),
                    plot = gridExtra::marrangeGrob(gg_plts, nrow = 1, ncol = 1),
                    device = "pdf",
                    units = "in",
                    width = width,
                    height = height)
  } else {

    return(gg_plts)

  }

}


#' Remove or filter out monomorphic loci from a data matrix or frame.
#' @param mydata A data frame of matrix consiting of markers as columns and
#' samples as rows.
#'
#' @returns A data frame of matrix object.
#'
#' @examples
#' # example code
#' library(panGenomeBreedr)
#' # Remove momomorphic loci
#' \donttest{
#' # Create a dummy marker frame data with monomorphic markers
#' mydf <- data.frame(S1_1001 = rep('A:A', 12),
#'                    S1_1011 = rep(c('T:T', 'T:C', 'C:C'), each = 4),
#'                    S2_1001 = rep(c('G:G', 'G:C', 'C:C'), each = 4),
#'                    S3_1101 = rep(c('G:G', 'C:C'), each = 6),
#'                    S1_1100 = rep(c('A:G', 'G:G', 'A:A'), each = 4),
#'                    S1_1201 = rep('C:A', 12),
#'                    row.names = sprintf('Ind_%02s', 1:12))
#'
#' # Remove monomorphic markers
#' mydf_filtered <- rm_mono(mydata = mydf)
#' }
#'
#'
#' @export

rm_mono <- function(mydata) {

  if (!inherits(mydata, what = 'data.frame') && !inherits(mydata, what = 'matrix')){

    stop('Argument `mydata` must either be a data frame or matrix.')

  }

  # Check for monomorphic loci without NAs
  mono_check <- apply(mydata, 2, function(x) length(unique(x[!is.na(x)])) > 1)

  # Remove monomorphic loci
  new_data <- mydata[,mono_check]

  return(new_data)

}

#' Calculate the proportion of recurrent parent background (RPP) fully recovered
#' in backcross progenies.
#' @param x A numeric matrix of marker genotypes for backcross progenies and recurrent
#' parent. Markers are columns and samples are rows.
#' @param map_file A data frame of map file consisting of marker IDs and their
#' chromosome numbers and positions as columns.
#' @param map_pos A character value indicating the column name for chromosome
#' positions in map file.
#' @param map_snp_ids A character value indicating the column name for SNP IDs
#' in \code{map_file}.
#' @param map_chr A character value indicating the column name for chromosome IDs
#' in \code{map_file}.
#' @param rp An integer or character value indicating the row index or sample ID
#' for the recurrent parent. if `NULL` the sample on the first row of \code{x} is
#' used as the recurrent parent.
#' @param rp_num_code A numeric value indicating the coding for the recurrent parent
#' marker background in \code{x}.
#' @param het_code A numeric value indicating the coding for heterozygous
#' marker background in \code{x}.
#' @param na_code A value indicating missing data in \code{x}.
#' @param weighted A logical value indicating whether RPP values should be weighted
#' or not.
#'
#' @details
#' The weighted RPP is computed as the sum of the product of relative distance
#' weighting of markers and genotype scores that match the recurrent parent.
#' Marker weights are obtained as half the normalized ratio of the physical
#' distance between marker `i` and marker `i+1` to the total distance for all
#' ordered markers on a chromosome.
#'
#' The computation excludes heterozygous loci and only considers chromosome regions
#' fully recovered as recurrent parent genetic background.
#'
#' @returns A data frame object comprising the RPP per chromosome and total RPP
#' of sample IDs.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#' # Create a numeric matrix of genotype scores for 10 markers and 5 samples
#' num_dat <- matrix(c(1, 1, 0.5, 0.5, 1, 1, 1, 1, 1, 1, rep(0, 10),
#'                     1, 1, 0.5, 1, 1, 1, 1, 1, 0, 1,
#'                     1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
#'                     1, 1, 0, 1, 1, 1, 1, 1, 1, 0.5 ),
#'                   byrow = TRUE, ncol = 10)
#'
#' rownames(num_dat) <- c('rp', 'dp', paste0('bc1_', 1:3))
#'
#' colnames(num_dat) <- c(paste0('S1', '_', c(floor(seq(1000, 10000, len = 5)))),
#'                             paste0('S2', '_', c(floor(seq(1000, 10000, len = 5)))))
#'
#' # Get map file by parsing SNP IDs
#' map_file <- parse_marker_ns(colnames(num_dat))
#'
#' # Calculate weighted RPP
#' rpp <- calc_rpp_bc(x = num_dat,
#'                       map_file = map_file,
#'                       map_chr = 'chr',
#'                       map_pos = 'pos',
#'                       map_snp_ids = 'snpid',
#'                       rp = 1,
#'                       rp_num_code = 1,
#'                       weighted = TRUE)
#' }
#'
#' @importFrom stats weighted.mean
#'
#' @export
#'
calc_rpp_bc <- function(x,
                        map_file,
                        map_chr = 'chr',
                        map_pos = 'pos',
                        map_snp_ids = 'snpid',
                        rp = NULL,
                        rp_num_code = 1,
                        het_code = 0.5,
                        na_code = NA,
                        weighted = TRUE) {

  if (missing(x)) stop("Provide numeric matrix of genotypes.")
  if (missing(map_file)) stop("Provide map_file of marker IDs.")

  # Resolve RP index
  rp_idx <- if (is.null(rp)) {
    1
  } else if (is.numeric(rp)) {
    rp
  } else if (is.character(rp)) {
    match(rp, rownames(x))
  }
  if (is.na(rp_idx)) stop("Recurrent parent not found.")

  # Extract RP genotype
  rp_geno <- x[rp_idx, ]

  # Handle NA code
  if (!is.na(na_code)) rp_geno[rp_geno == na_code] <- NA

  # Mask: remove loci with NA or heterozygous in RP
  keep_loci <- !is.na(rp_geno) & rp_geno != het_code
  x_filt <- x[, keep_loci, drop = FALSE]
  rp_geno <- rp_geno[keep_loci]

  # Subset map file to filtered markers
  snps_kept <- colnames(x_filt)
  map_sub <- merge(
    data.frame(snpid = snps_kept),
    map_file,
    by.x = "snpid",
    by.y = map_snp_ids,
    sort = FALSE
  )

  # Ensure SNPs in same order as in x_filt
  # map_sub <- map_sub[match(snps_kept, map_sub$snpid), ]
  map_sub <- order_markers(map_sub, chr_col = map_chr, pos_col = map_pos)

  # Split by chromosome
  chr_list <- split(seq_along(snps_kept), map_sub[[map_chr]])

  # Weighting function
  calc_weights <- function(pos) {
    n <- length(pos)
    if (n == 1) return(1)  # single-marker fallback
    d <- diff(pos)
    total <- sum(d, na.rm = TRUE)
    seg_wt <- d / total
    w <- numeric(n)
    w[1] <- seg_wt[1] / 2
    for (i in 2:(n - 1)) w[i] <- (seg_wt[i - 1] + seg_wt[i]) / 2
    w[n] <- seg_wt[n - 1] / 2
    return(w)
  }

  # Calculate RPP per chromosome
  rpp_chr <- lapply(names(chr_list), function(chr) {
    loci_idx <- chr_list[[chr]]
    geno_chr <- x_filt[, loci_idx, drop = FALSE]
    weights_chr <- if (weighted) calc_weights(map_sub[[map_pos]][loci_idx]) else NULL
    rpp_vals <- apply(geno_chr, 1, function(row) {
      if (weighted) {
        stats::weighted.mean(row == rp_num_code, weights_chr, na.rm = TRUE)
      } else {
        mean(row == rp_num_code, na.rm = TRUE)
      }
    })
    return(round(rpp_vals, 3))
  })

  names(rpp_chr) <- paste0("chr_", names(chr_list))

  # Combine into data frame
  rpp_chr_df <- as.data.frame(rpp_chr)
  rpp_chr_df$sample_id <- rownames(x_filt)

  # Compute total RPP
  if (weighted) {
    weights_all <- unlist(lapply(chr_list, function(idx) calc_weights(map_sub[[map_pos]][idx])))
    total_rpp <- apply(x_filt, 1, function(row) {
      stats::weighted.mean(row == rp_num_code, weights_all, na.rm = TRUE)
    })
  } else {
    total_rpp <- apply(x_filt, 1, function(row) {
      mean(row == rp_num_code, na.rm = TRUE)
    })
  }

  rpp_chr_df$total_rpp <- round(total_rpp, 3)

  return(rpp_chr_df[, c("sample_id", names(rpp_chr), "total_rpp")])
}


#' Compute theoretical RPP values for any specified backcross generation.
#' @param bc_gen A positive integer value indicating the backcross generation.
#' @param rpp2n A logical value indicating whether to compute theoretical RPP
#' values for the F1 all BC generations up to \code{bc_gen}.
#'
#' @returns A named vector of theoretical RPP values.
#'
#' @examples
#' # example code
#' # Calculate theoretical RPP values up to BC5
#' rpp_values <- calc_rpp_exp(bc_gen = 5, rpp2n = TRUE)
#'
#' @export

calc_rpp_exp <- function(bc_gen = 1,
                         rpp2n = FALSE) {

  # Coerce input gen into an integer
  bc_gen <- as.integer(bc_gen)

  if (bc_gen < 0) stop('Input a positve integer for the `bc_gen` argument.')

  # Loop over each generation and calculate RPP
  if (rpp2n == TRUE ) {

    if (bc_gen > 0) {

      # Create a vector to store RPP values for each generation
      rpp_values <- numeric(bc_gen + 1)

      for (n in 0:bc_gen) {
        rpp_values[n + 1] <- round(1 - 0.5^(n + 1), 4)
      }

      names(rpp_values) <- c('F1', paste0('bc', 1:bc_gen))

    } else stop('BC generation must be greater than zero.!')

  } else {

    rpp_values <- round(1 - 0.5^(bc_gen + 1), 4)
    names(rpp_values) <- paste0('bc', bc_gen)

  }

  return(rpp_values)

}


#' Visualize computed RPP values for BC progenies as a bar plot.
#' @param rpp_df A data frame containing the computed RPP values for each BC progeny.
#' @param rpp_col A character value indicating the column name of RPP values in
#' \code{rpp_df}.
#' @param rpp_sample_id A character value indicating the column name of progeny IDs
#' in \code{rpp_df}.
#' @param bc_gen An integer value indicating the BC generation for the progenies
#' \code{rpp_df}. This value is used to compute the nominal RPP values, if
#' \code{rpp_threshold} `= NULL`.
#' @param rpp_threshold A numeric value between  0 and 1 indicating the RPP
#' threshold for selecting BC progenies.
#' @param thresh_line_col A character value indicating the color of the threshold
#' line.
#' @param show_above_thresh A logical value indicating whether to subset \code{rpp_df}
#' to show only lines that have RPP values greater or equal to the desired RPP
#' threshold. Only the subset lines will be shown on the plot.
#' @param text_size A numeric value for setting text size.
#' @param text_scale_fct A numeric value for scaling text size. The default value
#' is `20\%` of \code{text_size}.
#' @param alpha A numeric value between 0 and 1 for modifying the
#' opacity of colors.
#' @param bar_width A numeric value for setting the width of plot bars.
#' @param bar_col A character value for setting the color to fill plot bars.
#' @param aspect_ratio A numeric value for setting the aspect ratio of the bar plot.
#' @param pdf A logical value indicating whether to save plot as a pdf graphic
#' device when TRUE or output plot in R when FALSE.
#' @param filename A character value for path or file name for saving pdf.
#' @param width A numeric value for the width of pdf device.
#' @param height A numeric value for the height of pdf device.
#'
#' @returns A graphical or ggplot object.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#'
#' # Observed RPP values
#' rpp_df <- data.frame(sample_id = c('rp', 'dp', paste0('bc3_', 1:8)),
#'                      rpp = c(1, 0, round(seq(0.75, 0.97, len = 8), 2)))
#'
#' # Generate bar plot for RPP values
#' rpp_barplot(rpp_df,
#'             rpp_col = 'rpp',
#'             rpp_threshold = 0.85,
#'             text_size = 18,
#'             text_scale_fct = 0.1,
#'             alpha = 0.9,
#'             bar_width = 0.5,
#'             aspect_ratio = 0.5,
#'             pdf = FALSE)
#' }
#'
#' @export
rpp_barplot <- function(rpp_df,
                        rpp_col = 'total_rpp',
                        rpp_sample_id = 'sample_id',
                        bc_gen = NULL,
                        rpp_threshold = NULL,
                        thresh_line_col = 'coral2',
                        show_above_thresh = FALSE,
                        text_size = 15,
                        text_scale_fct = 0.2,
                        alpha = 0.5,
                        bar_width = 0.5,
                        bar_col = 'cornflowerblue',
                        aspect_ratio = 0.5,
                        pdf = FALSE,
                        filename = 'rpp_barplot',
                        width = 8,
                        height = 6) {

  # Set global variables
  sample_id <- rpp <- NULL

  # Rebuild input data frame with standardized column names
  rpp_df <- data.frame(sample_id = rpp_df[, rpp_sample_id], rpp = rpp_df[, rpp_col])
  rpp_df$sample_id <- factor(rpp_df$sample_id, levels = rev(unique(rpp_df$sample_id )))

  # Use nominal RPP value as threshold if its value is NULL
  # Or set to 0.9 if both bc_gen and threshold value are NULL
  if (is.null(rpp_threshold)) {

    if (!is.null(bc_gen)) {

      rpp_threshold <- calc_rpp_exp(bc_gen)

    } else rpp_threshold <- 0.9

  }

  # Subset to show only lines that have RPP above the desired RPP threshold
  if (show_above_thresh) {

    rpp_df <- rpp_df[rpp_df[,'rpp'] >= rpp_threshold,]

  }




  plt <- ggplot2::ggplot(rpp_df, ggplot2::aes(x = rpp, y = sample_id)) +

    ggplot2::geom_bar(position = 'dodge',
                      stat = "identity",
                      width = bar_width,
                      color = ggplot2::alpha(bar_col, alpha),
                      fill = ggplot2::alpha(bar_col, alpha)) +
    ggplot2::geom_text(ggplot2::aes(x = .5, label = ifelse(rpp >= rpp_threshold, rpp * 100, "")),
                       color = 'black',
                       size = (text_size*bar_width)/(text_scale_fct*text_size)) +

    ggplot2::scale_x_continuous(limits = c(c(0, 1)), labels = function(x) paste0(x*100)) +

    ggplot2::geom_vline(xintercept = rpp_threshold, color = thresh_line_col, lwd = 2, lty = 2) +

    ggplot2::labs(title = paste('Recurrent parent genome recovery in progenies'),
                  x = "Observed RPP (%)", y = "Line") +

    ggplot2::theme(axis.text.x = ggplot2::element_text(size = text_size),
                   aspect.ratio = aspect_ratio,
                   plot.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(linewidth = 0.2),
                   plot.title = ggplot2::element_text(hjust = 0.5,
                                                      face = 'bold',
                                                      size = text_size),
                   axis.text = ggplot2::element_text(size  = text_size,
                                                     color = 'black'),
                   axis.title = ggplot2::element_text(size  = text_size,
                                                      face = 'bold',
                                                      color = 'black'))

  # Save plt as a pdf file
  if (pdf) {

    plt <- list(plt)

    ggplot2::ggsave(filename = paste0(filename, ".pdf"),
                    plot = gridExtra::marrangeGrob(plt, nrow = 1, ncol = 1),
                    device = "pdf",
                    units = "in",
                    width = width,
                    height = height)
  } else {

    return(plt)

  }


}


#' Annotate a heatmap to show introgressed loci positions.
#' @inheritParams cross_qc_ggplot
#' @param trait_pos A list object where the components correspond to the start and
#' end positions of trait loci to annotate on the heapmap.
#' @param text_scale_fct A numeric value for scaling text size. The default value
#' is `50\%` of the defined text size.
#'
#' @return A graphic object or ggplot.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#' # Create a numeric matrix of genotype scores for 10 markers and 5 samples
#' num_dat <- matrix(c(rep(1, 10), rep(0, 10),
#'                     1, 1, 0.5, 1, 1, 1, 1, 1, 0, 1,
#'                     1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
#'                     1, 1, 0, 1, 1, 1, 1, 1, 1, 0.5 ),
#'                   byrow = TRUE, ncol = 10)
#'
#' rownames(num_dat) <- c('rp', 'dp', paste0('bc1_', 1:3))
#' colnames(num_dat) <- paste0('S1', '_', c(floor(seq(1000, 10000, len = 8)),
#'                                          15000, 20000))
#'
#' # Get map file by parsing SNP IDs
#' map_file <- parse_marker_ns(colnames(num_dat))
#'
#' # Annotate a heatmap to show trait loci positions
#' cross_qc_annotate(x = num_dat,
#'                 map_file = map_file,
#'                 snp_ids = 'snpid',
#'                 chr = 'chr',
#'                 chr_pos = 'pos',
#'                 parents = c('rp', 'dp'),
#'                 trait_pos = list(loc1 = c(start = 2900, end = 4200),
#'                 loc2 = c(start = 14200, end = 15800)),
#'                 text_scale_fct = 0.5,
#'                 group_sz = 3L,
#'                 pdf = FALSE,
#'                 legend_title = 'BC1',
#'                 alpha = 0.8,
#'                 text_size = 15)
#'
#' }
#'
#' @export
cross_qc_annotate <- function(x,
                              map_file,
                              snp_ids = 'snpid',
                              chr = 'chr',
                              chr_pos = 'pos',
                              parents,
                              trait_pos,
                              group_sz = nrow(x)-2,
                              pdf = FALSE,
                              filename = 'background_heatmap',
                              legend_title = 'Heatmap_key',
                              col_mapping,
                              col_labels,
                              alpha = 0.9,
                              text_size = 12,
                              text_scale_fct = 0.5,
                              width = 10.5,
                              height = 8.5,
                              ...) {

  if (missing(x)) stop("Please provide a numeric matrix for the `x` argument.")
  if (!inherits(x, what = 'matrix')) stop("Argument `x` must be a matrix object.")
  if (missing(map_file)) stop("Please provide a map file for the `map_file` argument.")
  if (!inherits(map_file, what = 'data.frame')) stop("Argument `map_file` must be a data.frame object.")
  if (missing(parents)) stop("Please provide the value for the 'parents' argument.")
  if (length(parents) != 2) stop("The `parents` argument must be a character vector of length = 2.")
  if (missing(trait_pos)) stop("Please provide a trait positions for the `trait_pos` argument.")
  if (!inherits(trait_pos, what = 'list')) stop("Argument `trait_pos` must be a list object.")

  snpid <- value <- map_dist <- pos  <- NULL # Define global variables

  # Sort data in ascending order of chromosomes
  map_file <- order_markers(map_file, chr_col = chr, pos_col = chr_pos)

  # Calculate inter-marker distances and add it to map file
  map_dist <- lapply(split(map_file[, chr_pos], map_file[, chr]), FUN = diff)

  map_dist <- unlist(lapply(map_dist, FUN = function(x) c(0, x)))

  map_file <- cbind(map_file, map_dist = map_dist)

  # Subset parents data
  parent_dat <- rbind(x[which(rownames(x) == parents[1]),],
                      x[which(rownames(x) == parents[2]),])
  rownames(parent_dat) <- parents

  # Melt x into a data frame
  parent_dat <- gg_dat(num_mat = parent_dat,
                       map_file = map_file,
                       map_pos = chr_pos,
                       map_chr = chr,
                       map_snp_ids = snp_ids)

  # Subset progeny data
  progeny_dat <- x[!rownames(x) %in% unique(parent_dat$x),]

  # Create an index to split by
  grp_index <- rep(1:ceiling(nrow(progeny_dat) / group_sz), each = group_sz,
                   length.out = nrow(progeny_dat))

  # Split data frame into list of smaller data frames
  batches <- split(as.data.frame(progeny_dat), grp_index)

  batches <- lapply(batches, FUN = function(x) {

    progeny_dat <- gg_dat(num_mat = x,
                          map_file = map_file,
                          map_pos = chr_pos,
                          map_chr = chr,
                          map_snp_ids = snp_ids)

    grp <- rbind(parent_dat, progeny_dat)

    # Re-order levels of the sample ids
    grp$x <- factor(grp$x, levels = rev(c(sort(unique(parent_dat$x)),
                                          sort(unique(progeny_dat$x)))))
    grp$snpid <- factor(grp$snpid, levels = unique(grp$snpid))

    grp


  })

  nbatches <- length(batches) # Number of batches to plot

  # Create an empty list object to hold ggplots
  gg_plts <- vector(mode = 'list', length = nbatches)
  names(gg_plts) <- paste0('Batch', seq_len(nbatches))

  # Color and labels for plotting
  # Black = NAs; coral1 = RP; gold = Het; purple = DP; grey70 = Mono;
  # cornflowerblue = geno_error

  if (missing(col_mapping)) {

    col_mapping <- c('-5' = 'black',
                     '-2' = 'cornflowerblue',
                     '-1' = 'grey80',
                     '0' = 'purple2',
                     '0.5' = 'gold',
                     '1' = 'coral2')

  }

  if (missing(col_labels)) {

    col_labels <- c('-5' = "Missing",
                    '-2' = 'Error',
                    '-1' = "Monomorphic",
                    '0' = parents[2],
                    '0.5' = "Heterozygous",
                    '1' = parents[1])
  }


  # Annotate plot to show loci positions
  annotate_loc <- function(gg_obj, trait_pos) {

    # Get trait locus names if present, else generate names to use
    if (!is.null(names(trait_pos))) {

      pos_ns <- names(trait_pos)

    } else {

      pos_ns <- paste0('loc', 1:length(trait_pos))

    }

    gg_obj <- gg_obj

    for (i in seq_len(length(pos_ns))) {

      mid_pos <- mean(trait_pos[[i]], na.rm = TRUE)
      locus_ns <- pos_ns[i]

      gg_obj <- gg_obj +
        ggplot2::geom_vline(xintercept = trait_pos[[i]], linetype = 2,
                            lwd = 2, col = 'black') +
        ggplot2::geom_text(x = mid_pos,
                           y = length(unique(grp$x)) + 0.65,
                           label = locus_ns,
                           size = text_size/(text_scale_fct*text_size))
    }

    return(gg_obj)

  }

  for(i in seq_len(nbatches)) {

    grp <- batches[[i]]
    # Plot heat map
    plt <- ggplot2::ggplot(grp, ggplot2::aes(x = pos, y = x, fill = value,
                                             width = map_dist)) +
      ggplot2::geom_tile(lwd = 1, linetype = 1) +
      ggplot2::geom_rug(ggplot2::aes(x = pos),
                        sides = "b", length = grid::unit(0.02, "npc")) +
      ggplot2::scale_fill_manual( values = ggplot2::alpha(col_mapping, alpha),
                                  label = col_labels,
                                  name = legend_title) +

      ggplot2::xlab('Marker position (Mbp)') +
      ggplot2::scale_x_continuous(labels = function(x) paste0(x/1e6)) +
      ggplot2::expand_limits(y = c(1, length(unique(grp$x)) + 0.8)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = text_size,
                                                         face = 'bold'),
                     axis.ticks.length.y = ggplot2::unit(0.25, 'cm'),
                     axis.text.x = ggplot2::element_text(angle = 0, hjust = 0,
                                                         size = text_size),
                     axis.title.x = ggplot2::element_text(size = text_size,
                                                          face = 'bold'),
                     axis.title.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size = text_size),
                     legend.title = ggplot2::element_text(size = text_size)) +
      ggplot2::geom_hline(yintercept = c(as.numeric(grp$x) + .5, .5),
                          col = 'white', lwd = 2.5)

    # Annotate heatmap to show locus positions
    plt <- annotate_loc(gg_obj = plt, trait_pos)

    gg_plts[[i]] <- plt

  } # End of the loop

  if (pdf) {

    ggplot2::ggsave(filename = paste0(filename, ".pdf"),
                    plot = gridExtra::marrangeGrob(gg_plts, nrow = 1, ncol = 1),
                    device = "pdf",
                    units = "in",
                    width = width,
                    height = height)
  } else {

    return(gg_plts)

  }

}





#' Simulate raw SNP loci for any chromosome with or without LD.
#' @param nsnp An integer value specifying the number of SNPs to simulate.
#' @param nobs An integer value specifying the number of individuals to simulate.
#' @param chr An integer value specifying the chromosome number.
#' @param start An integer value specifying the position of the first SNP.
#' @param end An integer value specifying the position of the last SNP.
#' @param sep A separator for deriving the genotypes.
#' @param add_LD A logical value indicating whether to simulate SNPs in LD or not.
#' @param LD_range A numeric vector of `length = 2` indicating the range of LD
#' between SNPs, if \code{add_LD} = TRUE.
#'
#' @returns A data frame with the simulated SNPs as the columns and individuals
#' as rows.
#'
#' @examples
#' # example code
#' # Simulate data for 20 snp and 100 individuals on chr 1
#' geno_data <- sim_snp_dat(nsnp = 20,
#'                      nobs = 100,
#'                      chr = 1,
#'                      start = 1000,
#'                      end = 20000,
#'                      add_LD = TRUE,
#'                      LD_range = c(0.2, 1))
#'
#' @importFrom stats quantile
#' @importFrom stats rnorm
#'
#' @export
sim_snp_dat <- function(nsnp = 10L,
                        nobs = 100L,
                        chr = 1L,
                        start = 1000L,
                        end = 20000L,
                        sep = '/',
                        add_LD = FALSE,
                        LD_range = NULL) {

  # Get bi-allelic SNPs
  alleles <- replicate(nsnp, sample(c('A', 'C', 'G', 'T'), size = 2, replace = FALSE),
                       simplify = FALSE)

  # Vectorized creation of genotypes for each SNP
  genotypes <- lapply(alleles, function(x) {
    geno <- expand.grid(x, x)[-3, ]
    paste0(geno[, 1], sep, geno[, 2])
  })

  # Generate random genotypes for individuals
  dat <- as.data.frame(matrix(unlist(lapply(genotypes, sample, nobs, replace = TRUE)),
                              ncol = nsnp, byrow = FALSE))

  if (add_LD) {
    if (is.null(LD_range)) LD_range <- c(0, 1)

    # Create a Toeplitz correlation matrix
    ld_range <- seq(LD_range[2], LD_range[1], len = nsnp)
    cor_mat <- stats::toeplitz(ld_range)

    # Add a small value to the diagonal to ensure positive-definiteness
    cor_mat <- cor_mat + diag(rep(1e-6, nsnp))

    # # Simulate data with LD structure using rmvnorm
    # num_data <- mvtnorm::rmvnorm(nobs, mean = rep(0, nsnp), sigma = cor_mat)

    # Perform Cholesky decomposition on cor_mat -- faster
    chol_mat <- chol(cor_mat)

    # Generate standard normal random values
    nor_var <- matrix(stats::rnorm(nobs * nsnp), nobs, nsnp)

    # Apply Cholesky factor to introduce correlation
    num_data <- nor_var %*% chol_mat

    # Pre-compute quantile breaks -- vectorized
    quantile_breaks <- t(apply(num_data, 2, function(col) {
      stats::quantile(col, probs = c(0, 0.25, 0.75, 1))
    }))

    # Function to re-code correlated numeric matrix to a character matrix
    recode_dat <- function(column, labels, breaks) {
      cut(column,
          breaks = breaks,
          labels = labels,
          include.lowest = TRUE)
    }

    # mapply by using quantile_breaks and genotypes
    dat <- as.data.frame(mapply(recode_dat,
                                split(num_data, col(num_data)),
                                genotypes,
                                split(quantile_breaks, row(quantile_breaks))))

    # Reshuffle columns randomly
    dat <- dat[, sample(names(dat))]

  }

  # Set SNP names
  names(dat) <- paste0('S', chr, '_', floor(seq(start, end, len = nsnp)))
  rownames(dat) <- sprintf('ind_%03s', 1:nobs)

  return(dat)
}


#' Order marker IDs based on their chromosome numbers and positions in ascending
#' order.
#' @param x A data frame object for a map file or marker data file.
#' @param chr_col A character value indicating the column name for chromosome IDs
#' in \code{x}.
#' @param pos_col A character value indicating the column name for chromosome
#' positions in \code{x}.
#' @returns A data frame object of the same dimension as \code{x}.
#'
#' @examples
#' # example code
#' map_file <- data.frame(snpid = paste0('S', rep(1:2, 5), '_', 1001:1005),
#'                        chr = rep(1:2, 5),
#'                        pos = rep(1001:1005, 2))
#'
#' # Order map file
#' map_file <- order_markers(x = map_file,
#'                           chr_col = 'chr',
#'                           pos_col = 'pos')
#'
#' @export
order_markers <- function(x,
                          chr_col = 'chr',
                          pos_col = 'pos') {

  # Sort data in ascending order of chromosomes
  dat <- x[order(x[, chr_col], decreasing = FALSE),]

  # Split genotype data into chr batches
  grps <- split(dat, dat[, chr_col])

  # Function to order marker positions
  ord_pos <- function(x) {
    x[order(x[, pos_col], decreasing = FALSE),]
  }

  # Sort data in ascending order of physical positions per chromosome
  dat <- lapply(grps, FUN = ord_pos)

  # Row-bind sorted data for all chromosomes
  dat <- do.call(rbind, dat)

  return(dat)

}


#' Format marker names to comply with the Hapmap convention.
#' @param x A character vector of marker IDs to be formated.
#' @param map_file A data frame of map file consisting of marker IDs and their
#' chromosome numbers and positions as columns.
#' @param snpid_col A character value indicating the column name for marker IDs
#' in \code{map_file}.
#' @param chr_col A character value indicating the column name for chromosome IDs
#' in \code{map_file}.
#' @param pos_col A character value indicating the column name for chromosome
#' positions in \code{map_file}.
#' @param kasp_prefix A character value indicating the KASP marker prefix for the
#' crop species.
#' @param scaffold_prefix A character value indicating the KASP marker prefix for
#' the crop species.
#' @returns A character vector of markers names formated to comply with the Hapmap
#' format.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#'
#' # Marker IDs
#' snps <- c('snpSB00072', 'snpSB00106', 'snpSB00109', 'Sbv3.1_01_68028666I',
#'           'Sbv3.1_02_67884158W')
#'
#' # Map file for SNPs
#' map_file <- data.frame(snpid = snps,
#'                        chr = c(2, 5, 5, 1, 2),
#'                        pos = c(61811307, 838874, 1730282, NA, NA))
#'
#' # Format marker IDs to hapmap format
#' ns_new <- hapmap_ns_fmt(x = snps,
#'                         map_file = map_file,
#'                         snpid_col = 'snpid',
#'                         chr_col = 'chr',
#'                         pos_col = 'pos',
#'                         kasp_prefix = 'snpSB',
#'                         scaffold_prefix = 'Sbv')
#' }
#'
#' @details
#' This is an experimental function, so use it with caution. It can be used to
#' create marker IDs that can easily be parsed into a map file using the
#' `parse_marker_ns()` function. The function may not format all inputted marker
#' names in \code{x}.
#'
#' @export
#'
hapmap_ns_fmt <- function(x,
                          map_file,
                          snpid_col = 'SNP_name',
                          chr_col = 'chr',
                          pos_col = 'pos',
                          kasp_prefix = 'snpSB',
                          scaffold_prefix = 'Sbv') {

  snpid_new <- NULL

  # Format checks
  kasp_check <- paste0('^', kasp_prefix, '\\d+$')
  scaffold_check <- paste0('^', scaffold_prefix, '\\d+\\.\\d+_\\d{2}_\\d+\\w$')

  # Check if SNP ID is of the format "snpSB00072"; if TRUE, format using the
  # imported map data to the hapmap format: S(chr)(_)(pos): S1_100000
  if (any(grepl(kasp_check, x))) {

    # Get SNPs from map_dat with IDs that can not be parsed
    map_df <- map_file[grepl(kasp_check, map_file[, snpid_col]),]

    # Match SNPID in geno and df
    map_df <- merge(x = data.frame(snpid = x),
                    y = map_df,
                    by.x = 'snpid',
                    by.y = snpid_col,
                    sort = FALSE)

    # Add a new SNP name that can be parsed
    map_df$snpid_new <- paste0('S', map_df[,chr_col], '_', map_df[,pos_col])

    # Update SNP IDs in dartag data
    x[x %in% map_df$snpid] <- map_df$snpid_new

  }

  # Check if SNP ID is of the format "Sbv3.1_01_68028666I"; if TRUE,
  # extract the chr number and position from it
  # "Sbv<number>.<number>_<chr>_<pos><variant_type>
  if (any(grepl(scaffold_check, x))) {

    # Extract all SNPs with this format: "Sbv3.1_01_68028666I"
    snp_vec <- x[grepl(scaffold_check, x)]

    # Extract the chromosome number
    get_chr <- as.numeric(sub(".*_(\\d{2})_.*", "\\1", snp_vec))

    # Extract the chromosome position
    get_pos <- as.numeric(sub(".*_\\d{2}_(\\d+).*", "\\1", snp_vec))

    # Update SNP IDs in dartag data
    x[grepl(scaffold_check, x)] <- paste0('S', get_chr, '_', get_pos)
  }

  # Check if all SNP IDs in geno are of the hapmap format
  if(any(!grepl('^S\\d+_\\d+$', x))) {

    unformatted <- x[!grepl('^S\\d+_\\d+$', x)]

    warning_msg <- paste(length(unformatted), 'markers:',
                         paste0(unformatted, collapse = ', '),
                         'were not formatted!')

    warning(warning_msg)

  }

  return(x)

}

#' Find loci with unexpected homozygous genotype calls for artificial heterozygotes.
#' @param x A data matrix or frame with markers as columns and samples as rows.
#' @param rp_row,dp_row An integer or character value indicating the row index
#' or name of Parent 1 and 2.
#' @inheritParams parent_missing
#'
#' @details
#' Artificial heterozygotes are expected to show heterozygosity at polymorphic
#' loci between parents.Use this wrapper function to detect loci which did not
#' follow this prediction.
#'
#' @returns A list object with the following components:
#' 1) data frame of loci with unexpected genotype calls artificial heterozygotes
#'  if present.
#' 2) data frame of loci with expected genotype calls for artificial heterozygotes.
#'
#' @examples
#' # example code
#' library(panGenomeBreedr)
#' # Marker data
#' dat <- data.frame(snp1 = c('C:C', 'A:A', 'C:A', 'C:A'),
#'                   snp2 = c('C:C', '-:-', 'C:C', 'C:C'),
#'                   snp3 = c('T:T', 'C:C', 'C:T', 'C:T'),
#'                   snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
#'                   snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
#'                   row.names = c('rp', 'dp', 'art_het1', 'art_het2'))
#'
#' # Check for unexpected homozygous genotypes
#' homs_unexp <- find_unexp_homs(x = dat,
#'                               rp_row = 1,
#'                               dp_row = 2)$geno_unexp
#'
#' @export
#'
find_unexp_homs <- function(x,
                            rp_row,
                            dp_row,
                            na_code = NA) {

  if (missing(rp_row)) stop("Parent 1 row index number is missing!")
  if (missing(dp_row)) stop("Parent 2 row index number is missing!")

  # Replace na_code with NAs in marker data
  if (!is.na(na_code)) {

    x[x == na_code] <- NA

  }


  # Extract recurrent and donor parent data
  par_geno <- rbind(x[rp_row, ], x[dp_row, ])

  # Extract data for progenies
  prog_geno <- x[-c(rp_row, dp_row), ]

  # Check for unexpected
  col_index <- sapply(seq_len(ncol(x)), function(i) {
    any(prog_geno[, i] == par_geno[1, i]) | any(prog_geno[, i] == par_geno[2, i])
  })

  # Subset columns based on col_index
  if (any(col_index)) {

    geno_unexp <- as.data.frame(x[, col_index, drop = FALSE])
    colnames(geno_unexp) <- colnames(x)[col_index]

  } else geno_unexp <- NULL

  if (any(!col_index)) {

    geno_exp <- as.data.frame(x[, !col_index])
    colnames(geno_exp) <- colnames(x)[!col_index]

  } else geno_exp <- x

  res <- list(geno_unexp = geno_unexp,
              geno_exp = geno_exp)

  return(res)

}


#' Identify and subset InDel markers from a marker panel.
#' @inheritParams find_unexp_homs
#' @inheritParams parent_missing
#' @param sep A character used as separator for genotype calls, default is a
#' colon.
#' @param indel_sym A character value that indicates the symbols for a deletion.
#'
#' @returns A list object with the following components:
#' 1) data frame of InDel loci in marker panel if present.
#' 2) data frame of non-InDel markers in marker panel.
#'
#' @examples
#' # example code
#' library(panGenomeBreedr)
#'
#' # Marker data
#' dat <- data.frame(snp1 = c('C:C', 'A:A', 'C:A', 'C:A'),
#'                   snp2 = c('C:C', '-:-', 'C:C', 'C:C'),
#'                   snp3 = c('T:T', 'C:C', 'C:T', 'C:T'),
#'                   snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
#'                   snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
#'                   row.names = c('rp', 'dp', 'art_het1', 'art_het2'))
#'
#' # Find InDel loci
#' geno_indel <- find_indels(x = dat,
#'                           rp_row = 1,
#'                           dp_row = 2)$geno_indel
#'
#'
#' @export

find_indels  <- function(x,
                         rp_row,
                         dp_row,
                         indel_sym = '-',
                         sep = ':',
                         na_code = NA) {

  if (missing(rp_row)) stop("Parent 1 row index number is missing!")
  if (missing(dp_row)) stop("Parent 2 row index number is missing!")

  # Replace na_code with NAs in marker data
  if (!is.na(na_code)) {

    x[x == na_code] <- NA

  }

  # Extract recurrent and donor parent data
  par_dat <- rbind(x[rp_row,], x[dp_row,])

  # Pattern search for heterozygous genotype
  regexp <- paste0('^', indel_sym, '\\', sep, '\\', indel_sym, '$')

  # Col_index with logical values
  col_index <- apply(par_dat, 2, FUN = function(x) any(grepl(regexp, x)))

  # Subset columns based on col_index
  if (any(col_index)) {

    geno_indel <- as.data.frame(x[, col_index])
    colnames(geno_indel) <- colnames(x)[col_index]

  } else geno_indel <- NULL


  if (any(!col_index)) {

    geno_non_indel <- as.data.frame(x[, !col_index])
    colnames(geno_non_indel) <- colnames(x)[!col_index]

  } else geno_non_indel <- x

  res <- list(geno_indel = geno_indel,
              geno_non_indel = geno_non_indel)

  return(res)


}


#' Identify and subset loci with any parent missing genotype.
#' @inheritParams find_unexp_homs
#' @param na_code A value indicating missing data in \code{x}.
#'
#' @returns A list object with the following components:
#' 1) data frame of loci with at least one parent genotype missing, if present.
#' 2) data frame of loci with all parent genotype present.
#'
#' @examples
#' # example code
#' library(panGenomeBreedr)
#'
#' # Marker data
#' dat <- data.frame(snp1 = c('C:C', 'A:A', 'C:A', 'C:A'),
#'                   snp2 = c('C:C', NA, 'C:C', 'C:C'),
#'                   snp3 = c(NA, 'C:C', 'C:T', 'C:T'),
#'                   snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
#'                   snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
#'                   row.names = c('rp', 'dp', 'ind_1', 'ind_2'))
#'
#' # Find loci with at least one missing parent genotype
#' par_miss <- parent_missing(x = dat,
#'                           rp_row = 1,
#'                           dp_row = 2)$par_missing
#'
#' @export
parent_missing  <- function(x,
                            rp_row,
                            dp_row,
                            na_code = NA) {

  if (missing(rp_row)) stop("Parent 1 row index number is missing!")
  if (missing(dp_row)) stop("Parent 2 row index number is missing!")

  # Replace na_code with NAs in marker data
  if (!is.na(na_code)) {

    x[x == na_code] <- NA

  }

  # Subset parent marker data
  par_dat <- rbind(x[rp_row,], x[dp_row,])

  # Col_index with logical values
  col_index <- apply(par_dat, 2, FUN = function(x) anyNA(x))

  # Subset columns based on col_index
  if (any(col_index)) {

    par_missing <- as.data.frame(x[, col_index])
    colnames(par_missing) <- colnames(x)[col_index]

  } else par_missing <- NULL

  if (any(!col_index)) {

    par_present <- as.data.frame(x[, !col_index])
    colnames(par_present) <- colnames(x)[!col_index]

  } else par_present <- x

  res <- list(par_missing = par_missing,
              par_present = par_present)

  return(res)

}


#' Identify and subset loci with any heterozygous parent genotype.
#' @inheritParams find_indels
#' @inheritParams parent_missing
#'
#' @returns A list object with the following components:
#' 1) data frame of loci with at least one hetrozygous parent genotype, if present.
#' 2) data frame of loci with all homozygous parent genotype.
#'
#' @examples
#' # example code
#' library(panGenomeBreedr)
#'
#' # Marker data
#' dat <- data.frame(snp1 = c('C:A', 'A:A', 'C:A', 'C:A'),
#'                   snp2 = c('C:C', 'G:G', 'C:C', 'C:C'),
#'                   snp3 = c('C:T', 'C:C', 'C:T', 'C:T'),
#'                   snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
#'                   snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
#'                   row.names = c('rp', 'dp', 'ind_1', 'ind_2'))
#'
#' # Find loci with at least one heterozygous parent genotype
#' par_het <- parent_het(x = dat,
#'                       rp_row = 1,
#'                       dp_row = 2,
#'                       sep = ':')$par_het
#'
#' @export
parent_het  <- function(x,
                        rp_row,
                        dp_row,
                        sep = ':',
                        na_code = NA) {

  if (missing(rp_row)) stop("Parent 1 row index number is missing!")
  if (missing(dp_row)) stop("Parent 2 row index number is missing!")

  # Replace na_code with NAs in marker data
  if (!is.na(na_code)) {

    x[x == na_code] <- NA

  }

  # Subset parent marker data
  par_dat <- rbind(x[rp_row,], x[dp_row,])

  # Check if parent is heterozygous
  is_het <- function(x) {

    # Split the genotype into alleles
    alleles <- unlist(strsplit(x, sep))

    # Check if the alleles are different
    return(length(unique(alleles)) == 2)

  }

  # Col_index with logical values
  col_index <- apply(par_dat, MARGIN = c(1, 2), FUN = is_het)
  col_index <- apply(col_index, MARGIN = 2, FUN = function(x) any(x))

  # Subset columns based on col_index
  if (any(col_index)) {

    par_het <- as.data.frame(x[, col_index])
    colnames(par_het) <- colnames(x)[col_index]

  } else par_het <- NULL

  if (any(!col_index)) {

    par_hom <- as.data.frame(x[, !col_index])
    colnames(par_hom) <- colnames(x)[!col_index]

  } else par_hom <- x

  res <- list(par_het = par_het,
              par_hom = par_hom)

  return(res)

}


#' Select polymorphic loci between two parents in a marker panel.
#' @inheritParams find_indels
#'
#' @examples
#' # example code
#'
#' # Marker data
#' dat <- data.frame(snp1 = c('C:C', 'A:A', 'C:A', 'C:A'),
#'                   snp2 = c('C:C', 'C:C', 'C:C', 'C:C'),
#'                   snp3 = c('T:T', 'C:C', 'C:T', 'C:T'),
#'                   snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
#'                   snp5 = c('T:T', 'T:T', 'T:A', 'T:A'),
#'                   row.names = c('rp', 'dp', 'art_het1', 'art_het2'))
#'
#' # Find polymorphic loci
#' poly_loci <- parent_poly(x = dat,
#'                          rp_row = 1,
#'                          dp_row = 2,
#'                          sep = ':')
#'
#' @returns A data matrix or frame object of polymorphic loci between two parental
#' lines.
#'
#' @export
#'
parent_poly <- function(x,
                        rp_row,
                        dp_row,
                        sep = ':',
                        na_code = NA) {

  if (missing(rp_row)) stop("Parent 1 row index number is missing!")
  if (missing(dp_row)) stop("Parent 2 row index number is missing!")

  # Replace na_code with NAs in marker data
  if (!is.na(na_code)) {

    x[x == na_code] <- NA

  }

  # Extract recurrent and donor parent data as a character vector
  par_1 <- unlist(x[rp_row,])
  par_2 <- unlist(x[dp_row,])

  is_poly <- function(x) {

    # Get unique alleles for each parent
    al_p1 <- unique(unlist(strsplit(par_1[x], sep)))
    al_p2 <- unique(unlist(strsplit(par_2[x], sep)))

    # Check if intersect between parent alleles is less than their union
    # TRUE is polymorphic and FALSE if monomorphic
    length(intersect(al_p1, al_p2)) < length(union(al_p1, al_p2))
  }

  # Get the number of markers in x
  ncols <- ncol(x)

  # Check if parents are different at each SNP locus
  poly_loci <- sapply(seq_len(ncols), FUN = is_poly)

  # Filter loci to only polymorphic ones
  poly_genos <- x[, poly_loci, drop = FALSE]

  return(poly_genos)

}


#' Identify lines that possess favorable alleles for target loci using trait
#' predictive markers.
#' @param geno_data A data frame or matrix of marker genotype data with lines as
#' rows and markers as columns.
#' @param fore_marker_info A data frame of meta data for the trait-predictive markers
#' consisting of marker names, favorable and alternate alleles as variables.
#' @param fore_marker_col A character value specifying the column name for trait-predictive
#' markers in \code{fore_marker_info}.
#' @param fav_allele_col A character value specifying the column name for favorable
#' allele in \code{fore_marker_info}.
#' @param alt_allele_col A character value specifying the column name for alternate
#' allele in \code{fore_marker_info}.
#' @param select_type A character value of one of three options: `homo` to select
#' lines that are homozygous for the favorable allele; `hetero` to select lines
#' that are heterozygous; and `both` to select both favorable homozygotes and
#' heterozygotes.
#' @param sep A character value representing the separator used for genotype calling.
#'
#' @returns A binary data frame conversion of the original marker genotype data
#' for presence(1s) and absence(0s) of favorable alleles of target loci.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#'
#' # Marker genotype data
#' geno <- data.frame(SNP1 = c("A:A", "A:G", "G:G", "A:A"),
#'                    SNP2 = c("C:C", "C:T", "T:T", "C:T"),
#'                    SNP3 = c("G:G", "G:G", "A:G", "A:A"),
#'                    row.names = c("Line1", "Line2", "Line3", "Line4"))
#'
#' # Trait predictive markers meta data
#' marker_info <- data.frame(qtl_markers = paste0('SNP', 1:3),
#'                           fav_alleles = c('A', 'C', 'G'),
#'                           alt_alleles = c('G', 'T', 'A'))
#'
#' # Select lines where genotype is either homozygous for favorable alleles at all loci
#' foreground_select(geno_data = geno,
#'                   fore_marker_info = marker_info,
#'                   fore_marker_col = 'qtl_markers',
#'                   fav_allele_col = 'fav_alleles',
#'                   alt_allele_col = 'alt_alleles',
#'                   select_type = "homo")
#'
#' }
#'
#' @export
#'
foreground_select <- function(geno_data,
                              fore_marker_info,
                              fore_marker_col,
                              fav_allele_col,
                              alt_allele_col,
                              select_type = c("homo", "hetero", "both"),
                              sep = ":") {

  out_type <- match.arg(select_type) # Check for output type

  # Get line names
  lines <- rownames(geno_data)

  # QTL marker names
  qtl_loci <- fore_marker_info[, fore_marker_col]

  # Get favorite alleles for QTL markers
  fav_alleles <- fore_marker_info[, fav_allele_col]
  alt_alleles <- fore_marker_info[, alt_allele_col]

  # Empty matrix to hold results
  res_mat <- matrix(FALSE, nrow = nrow(geno_data), ncol = length(qtl_loci))

  # Check if QTL marker names are present in marker data
  if (!all(qtl_loci %in% colnames(geno_data))) {

    stop(paste(qtl_loci[!qtl_loci %in% colnames(geno_data)], 'not found in marker genotype data.'))

  }

  for (qtl in seq_len(length(qtl_loci))) {

    marker <- qtl_loci[qtl]

    # Subset favorite and alternate alleles
    fav_al <- fav_alleles[qtl]
    alt_al <- alt_alleles[qtl]

    if (!marker %in% colnames(geno_data)) next

    genos <- geno_data[[marker]]

    homo_fav <- paste0(fav_al, sep, fav_al)
    het1 <- paste0(fav_al, sep, alt_al)
    het2 <- paste0(alt_al, sep, fav_al)

    if (select_type == "homo") {

      res_mat[, qtl] <- genos == homo_fav

    } else if (select_type == "hetero") {

      res_mat[, qtl] <- genos == het1 | genos == het2

    } else if (select_type == "both") {

      res_mat[, qtl] <- genos == homo_fav | genos == het1 | genos == het2

    }

  }

  # Get the number of counts for loci with favorable alleles
  fav_count <- rowSums(res_mat, na.rm = TRUE)

  res_df <- data.frame(line = lines, fav_count = fav_count)

  res_mat <- as.data.frame(apply(res_mat, 2, as.integer))

  dimnames(res_mat) <- list(rownames(geno_data), qtl_loci)

  # Select lines that meet minimum threshold for favorable counts
  # selected <- res_df[res_df$fav_count >= threshold,]

  return(res_mat)

}

#' Extracts lines that have a combination of favorable alleles across target loci.
#' @param mat A binary matrix/data frame with 1s and 0s. Rows = lines, Columns = target loci.
#' @param present Character vector of column names that must be 1 (present).
#' @param absent  Character vector of column names that must be 0 (absent).
#'
#' @return Character vector of lines matching the intersection criteria.
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#'
#' # Marker genotype data
#' geno <- data.frame(SNP1 = c("A:A", "A:G", "G:G", "A:A"),
#'                    SNP2 = c("C:C", "C:T", "T:T", "C:T"),
#'                    SNP3 = c("G:G", "G:G", "A:G", "A:A"),
#'                    row.names = c("Line1", "Line2", "Line3", "Line4"))
#'
#' # Trait predictive markers meta data
#' marker_info <- data.frame(qtl_markers = paste0('SNP', 1:3),
#'                           fav_alleles = c('A', 'C', 'G'),
#'                           alt_alleles = c('G', 'T', 'A'))
#'
#' # Select lines where genotype is either homozygous for favorable alleles at all loci
#' foreground_select(geno_data = geno,
#'                   fore_marker_info = marker_info,
#'                   fore_marker_col = 'qtl_markers',
#'                   fav_allele_col = 'fav_alleles',
#'                   alt_allele_col = 'alt_alleles',
#'                   select_type = "homo") |>
#'
#'                   find_lines(present = colnames(geno))
#'
#'
#'
#' }
#'
#' @export
#'
find_lines <- function(mat,
                       present = character(),
                       absent = character()) {

  # Start with a logical vector that is TRUE for all rows
  cond <- rep(TRUE, nrow(mat))

  # If 'present' QTLs are specified, update cond to only TRUE for rows where all are 1
  if (length(present)) {

    cond <- cond & apply(mat[, present, drop = FALSE], 1, function(x) all(x == 1))

  }

  # If 'absent' QTLs are specified, update cond to only TRUE for rows where all are 0
  if (length(absent)) {

    cond <- cond & apply(mat[, absent, drop = FALSE], 1, function(x) all(x == 0))

  }

  # Return the row names (e.g., line IDs) that match the condition
  rownames(mat)[cond]
}


#' Create a heatmap to visualize and compare the genetic genetic backgrounds of
#' genotypes/lines with or without annotation for introgressed loci.
#' @inheritParams cross_qc_ggplot
#' @param trait_pos A list object where the components correspond to the chromosome number,
#' start and end positions of trait loci to annotate on the heapmap.
#' @param text_scale_fct A numeric value for scaling text size. The default value
#' is `50\%` of the defined text size.
#' @param panel_fill A character value for setting the panel background fill color.
#' @param panel_col A character value for setting the panel background border color.
#' @param label_offset A numeric value indicating the position of the trait loci
#' text labels on the heatmap. It is positioned on the donor parent by default.
#'
#' @return A ggplot graphical object.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#' # Create a numeric matrix of genotype scores for 10 markers and 5 samples
#' num_dat <- matrix(c(rep(1, 10), rep(0, 10),
#'                     1, 1, 0.5, 1, 1, 1, 1, 1, 0, 1,
#'                     1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
#'                     1, 1, 0, 1, 1, 1, 1, 1, 1, 0.5 ),
#'                   byrow = TRUE, ncol = 10)
#'
#' rownames(num_dat) <- c('rp', 'dp', paste0('bc1_', 1:3))
#' colnames(num_dat) <- paste0('S1', '_', c(floor(seq(1000, 10000, len = 8)),
#'                                          15000, 20000))
#'
#' # Get map file by parsing SNP IDs
#' map_file <- parse_marker_ns(colnames(num_dat))
#'
#' # Annotate a heatmap to show trait loci positions
#' cross_qc_heatmap(x = num_dat,
#'                 map_file = map_file,
#'                 snp_ids = 'snpid',
#'                 chr = 'chr',
#'                 chr_pos = 'pos',
#'                 parents = c('rp', 'dp'),
#'                 trait_pos = list(loc1 = c(chr = 1, start = 2900, end = 4200),
#'                 loc2 = c(chr = 1, start = 14200, end = 15800)),
#'                 text_scale_fct = 0.5,
#'                 group_sz = 3L,
#'                 pdf = FALSE,
#'                 legend_title = 'BC1',
#'                 alpha = 0.8,
#'                 text_size = 15)
#'
#' }
#'
#' @export
#'
cross_qc_heatmap <- function(x,
                             map_file,
                             snp_ids = 'snpid',
                             chr = 'chr',
                             chr_pos = 'pos',
                             parents,
                             trait_pos = NULL,
                             group_sz = nrow(x) - 2,
                             pdf = FALSE,
                             filename = 'background_heatmap',
                             legend_title = 'Heatmap_key',
                             col_mapping,
                             col_labels,
                             panel_fill = 'grey80',
                             panel_col = 'white',
                             alpha = 0.9,
                             text_size = 12,
                             text_scale_fct = 0.5,
                             width = 9,
                             height = 6.5,
                             label_offset = -1,
                             ...) {

  if (missing(x)) stop("Please provide a numeric matrix for the `x` argument.")
  if (!inherits(x, what = 'matrix')) stop("Argument `x` must be a matrix object.")
  if (missing(map_file)) stop("Please provide a map file for the `map_file` argument.")
  if (!inherits(map_file, what = 'data.frame')) stop("Argument `map_file` must be a data.frame object.")
  if (missing(parents)) stop("Please provide the value for the 'parents' argument.")
  if (length(parents) != 2) stop("The `parents` argument must be a character vector of length = 2.")

  snpid <- value <- map_dist <- pos <- cum_dist <- NULL

  if (inherits(map_file[, chr], what = 'numeric')) {
    map_file[, chr] <- sprintf('Chr%02s', map_file[, chr])
    map_file[, chr] <- factor(map_file[, chr], levels = unique(map_file[, chr]))
  }

  map_file <- order_markers(map_file, chr_col = chr, pos_col = chr_pos)
  map_dist <- lapply(split(map_file[, chr_pos], map_file[, chr]), FUN = diff)
  cum_dist <- unname(unlist(lapply(map_dist, FUN = function(x) c(0, cumsum(x)))))
  map_dist <- unname(unlist(lapply(map_dist, FUN = function(x) c(0, x))))
  map_file <- cbind(map_file, map_dist = map_dist, cum_dist)

  parent_dat <- rbind(x[which(rownames(x) == parents[1]),],
                      x[which(rownames(x) == parents[2]),])
  rownames(parent_dat) <- parents

  parent_dat <- gg_dat(num_mat = parent_dat,
                       map_file = map_file,
                       map_pos = chr_pos,
                       map_chr = chr,
                       map_snp_ids = snp_ids)

  progeny_dat <- x[!rownames(x) %in% unique(parent_dat$x),]
  grp_index <- rep(1:ceiling(nrow(progeny_dat) / group_sz), each = group_sz,
                   length.out = nrow(progeny_dat))
  batches <- split(as.data.frame(progeny_dat), grp_index)

  batches <- lapply(batches, FUN = function(x) {
    progeny_dat <- gg_dat(num_mat = x,
                          map_file = map_file,
                          map_pos = chr_pos,
                          map_chr = chr,
                          map_snp_ids = snp_ids)
    grp <- rbind(parent_dat, progeny_dat)
    grp$x <- factor(grp$x, levels = rev(c(sort(unique(parent_dat$x)),
                                          sort(unique(progeny_dat$x)))))
    grp
  })

  nbatches <- length(batches)
  gg_plts <- vector(mode = 'list', length = nbatches)
  names(gg_plts) <- paste0('Batch', seq_len(nbatches))

  if (missing(col_mapping)) {
    col_mapping <- c('-5' = 'deeppink', '-2' = 'cornflowerblue', '-1' = 'beige',
                     '0' = 'purple2', '0.5' = 'gold', '1' = 'coral2')
  }
  if (missing(col_labels)) {
    col_labels <- c('-5' = "Missing", '-2' = 'Error', '-1' = "Monomorphic",
                    '0' = parents[2], '0.5' = "Heterozygous", '1' = parents[1])
  }

  annotate_loc <- function(gg_obj, trait_pos, grp, label_offset) {

    xintercept <- chr <- mid_pos <- NULL
    if (is.null(names(trait_pos))) {

      names(trait_pos) <- paste0('loc', 1:length(trait_pos))
    }
    loc_ns <- as.list(names(trait_pos))

    trait_pos_ls <- lapply(trait_pos, FUN = function(x) {
      if (inherits(x[1], what = 'numeric')) {
        data.frame(xintercept = x[2:3], chr = sprintf('Chr%02s', x[1]))
      } else {
        data.frame(xintercept = x[2:3], chr = x[1])
      }
    })

    trait_lab <- mapply(FUN = function(x, y) {
      if (inherits(x[1], what = 'numeric')) {
        data.frame(mid_pos = mean(x[2:3], na.rm = TRUE),
                   loc_ns = y,
                   chr = sprintf('Chr%02s', x[1]))
      } else {
        data.frame(mid_pos = mean(x[2:3], na.rm = TRUE),
                   loc_ns = y,
                   chr = x[1])
      }
    }, trait_pos, loc_ns, SIMPLIFY = FALSE)

    for (i in seq_along(trait_pos)) {
      gg_obj <- gg_obj +

        ggplot2::geom_text(data = trait_lab[[i]],
                           ggplot2::aes(x = mid_pos, label = loc_ns),
                           col = "white",
                           hjust = 0.5,
                           y = length(unique(grp$x)) + label_offset,
                           size = text_size / (text_scale_fct * text_size),
                           inherit.aes = FALSE)
    }
    return(gg_obj)
  }

  for (i in seq_len(nbatches)) {
    grp <- batches[[i]]

    plt <- ggplot2::ggplot(grp, ggplot2::aes(x = pos, y = x, fill = value, width = map_dist)) +
      ggplot2::geom_tile(lwd = 1, linetype = 1) +

      ggplot2::scale_fill_manual(values = ggplot2::alpha(col_mapping, alpha),
                                 label = col_labels,
                                 name = legend_title) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5))) +
      ggplot2::xlab('Marker distance (Mbp)') +

      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  labels = function(x) paste0(x / 1e6)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::geom_hline(yintercept = c(as.numeric(grp$x) + 0.5, 0.5),
                          col = 'white', lwd = 2.5) +
      ggplot2::facet_grid(rows = NULL, cols = ggplot2::vars(chr),
                          space = "free", scales = "free_x") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = text_size, face = 'bold'),
                     axis.ticks.length.y = ggplot2::unit(0.25, 'cm'),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, size = text_size),
                     axis.title.x = ggplot2::element_text(size = text_size, face = 'bold'),
                     axis.title.y = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill = panel_fill, colour = panel_col,
                                                              linewidth = 2, linetype = "solid"),
                     legend.text = ggplot2::element_text(size = text_size),
                     legend.title = ggplot2::element_text(size = text_size),
                     strip.background = ggplot2::element_rect(fill = "grey90", colour = "grey90",
                                                              linewidth = 0.5, linetype = "solid"),
                     strip.text = ggplot2::element_text(size = text_size))

    if (!is.null(trait_pos)) {
      if (!inherits(trait_pos, what = 'list')) stop("Argument `trait_pos` must be a list object.")
      highlight_df <- do.call(rbind, lapply(trait_pos, function(loc) {
        chr_val <- if (is.numeric(loc[1])) sprintf("Chr%02s", loc[1]) else loc[1]
        pos_range <- as.numeric(loc[2:3])
        grp[as.character(grp$chr) == chr_val & grp$pos >= pos_range[1] &
              grp$pos <= pos_range[2] & grp$value == 0, ]
      }))

      plt <- plt + ggplot2::geom_tile(data = highlight_df,
                                      ggplot2::aes(x = pos, y = x, width = map_dist, height = 1),
                                      color = "black", lwd = 1.2, lty = 1,
                                      fill = NA,
                                      inherit.aes = FALSE)

      plt <- annotate_loc(gg_obj = plt, trait_pos = trait_pos, grp = grp,
                          label_offset = label_offset)
    }

    gg_plts[[i]] <- plt
  }

  if (pdf) {
    ggplot2::ggsave(filename = paste0(filename, ".pdf"),
                    plot = gridExtra::marrangeGrob(gg_plts, nrow = 1, ncol = 1),
                    device = "pdf",
                    units = "in",
                    width = width,
                    height = height)
  } else {
    return(gg_plts)
  }
}


#' Visualize genotype backgrounds with optional QTL annotations.
#' @description
#' Create a heatmap to compare the genetic backgrounds of genotypes or breeding lines
#' across chromosomes, with optional annotation of known QTL regions.
#' The function supports highlighting trait loci with vertical dashed lines and
#' top-aligned, rotated labels for clear interpretation.
#' @inheritParams cross_qc_heatmap
#' @return A ggplot graphical object.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#' # Create a numeric matrix of genotype scores for 10 markers and 5 samples
#' num_dat <- matrix(c(rep(1, 10), rep(0, 10),
#'                     1, 1, 0.5, 1, 1, 1, 1, 1, 0, 1,
#'                     1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
#'                     1, 1, 0, 1, 1, 1, 1, 1, 1, 0.5 ),
#'                   byrow = TRUE, ncol = 10)
#'
#' rownames(num_dat) <- c('rp', 'dp', paste0('bc1_', 1:3))
#' colnames(num_dat) <- paste0('S1', '_', c(floor(seq(1000, 10000, len = 8)),
#'                                          15000, 20000))
#'
#' # Get map file by parsing SNP IDs
#' map_file <- parse_marker_ns(colnames(num_dat))
#'
#' # Annotate a heatmap to show QTL positions with vertical lines
#' cross_qc_heatmap2(x = num_dat,
#'                 map_file = map_file,
#'                 snp_ids = 'snpid',
#'                 chr = 'chr',
#'                 chr_pos = 'pos',
#'                 parents = c('rp', 'dp'),
#'                 trait_pos = list(loc1 = c(chr = 1, pos = 2900),
#'                 loc2 = c(chr = 1, pos = 14200)),
#'                 text_scale_fct = 0.5,
#'                 group_sz = 3L,
#'                 pdf = FALSE,
#'                 legend_title = 'BC1',
#'                 alpha = 0.8,
#'                 text_size = 15)
#'
#' }
#'
#' @export
#'
cross_qc_heatmap2 <- function(x,
                              map_file,
                              snp_ids = 'snpid',
                              chr = 'chr',
                              chr_pos = 'pos',
                              parents,
                              trait_pos = NULL,
                              group_sz = nrow(x) - 2,
                              pdf = FALSE,
                              filename = 'background_heatmap',
                              legend_title = 'Heatmap_key',
                              col_mapping,
                              col_labels,
                              panel_fill = 'grey80',
                              panel_col = 'white',
                              alpha = 0.9,
                              text_size = 12,
                              text_scale_fct = 0.5,
                              label_offset = 0.5,
                              width = 9,
                              height = 6.5) {

  if (missing(x)) stop("Please provide a numeric matrix for the `x` argument.")
  if (!inherits(x, what = 'matrix')) stop("Argument `x` must be a matrix object.")
  if (missing(map_file)) stop("Please provide a map file for the `map_file` argument.")
  if (!inherits(map_file, what = 'data.frame')) stop("Argument `map_file` must be a data.frame object.")
  if (missing(parents)) stop("Please provide the value for the 'parents' argument.")
  if (length(parents) != 2) stop("The `parents` argument must be a character vector of length = 2.")

  snpid <- value <- map_dist <- pos <- cum_dist <- NULL

  if (inherits(map_file[, chr], what = 'numeric')) {
    map_file[, chr] <- sprintf('Chr%02s', map_file[, chr])
    map_file[, chr] <- factor(map_file[, chr], levels = unique(map_file[, chr]))
  }

  map_file <- order_markers(map_file, chr_col = chr, pos_col = chr_pos)
  map_dist <- lapply(split(map_file[, chr_pos], map_file[, chr]), FUN = diff)
  cum_dist <- unname(unlist(lapply(map_dist, FUN = function(x) c(0, cumsum(x)))))
  map_dist <- unname(unlist(lapply(map_dist, FUN = function(x) c(0, x))))
  map_file <- cbind(map_file, map_dist = map_dist, cum_dist)

  parent_dat <- rbind(x[which(rownames(x) == parents[1]),],
                      x[which(rownames(x) == parents[2]),])
  rownames(parent_dat) <- parents

  parent_dat <- gg_dat(num_mat = parent_dat,
                       map_file = map_file,
                       map_pos = chr_pos,
                       map_chr = chr,
                       map_snp_ids = snp_ids)

  progeny_dat <- x[!rownames(x) %in% unique(parent_dat$x),]
  grp_index <- rep(1:ceiling(nrow(progeny_dat) / group_sz), each = group_sz,
                   length.out = nrow(progeny_dat))
  batches <- split(as.data.frame(progeny_dat), grp_index)

  batches <- lapply(batches, FUN = function(x) {
    progeny_dat <- gg_dat(num_mat = x,
                          map_file = map_file,
                          map_pos = chr_pos,
                          map_chr = chr,
                          map_snp_ids = snp_ids)
    grp <- rbind(parent_dat, progeny_dat)
    grp$x <- factor(grp$x, levels = rev(c(sort(unique(parent_dat$x)),
                                          sort(unique(progeny_dat$x)))))
    grp
  })

  nbatches <- length(batches)
  gg_plts <- vector(mode = 'list', length = nbatches)
  names(gg_plts) <- paste0('Batch', seq_len(nbatches))

  if (missing(col_mapping)) {
    col_mapping <- c('-5' = 'deeppink', '-2' = 'cornflowerblue', '-1' = 'beige',
                     '0' = 'purple2', '0.5' = 'gold', '1' = 'coral2')
  }
  if (missing(col_labels)) {
    col_labels <- c('-5' = "Missing", '-2' = 'Error', '-1' = "Monomorphic",
                    '0' = parents[2], '0.5' = "Heterozygous", '1' = parents[1])
  }

  annotate_loc <- function(gg_obj, trait_pos, grp) {
    if (is.null(names(trait_pos))) {
      names(trait_pos) <- paste0('loc', seq_along(trait_pos))
    }

    loc_ns <- mid_pos <- NULL

    loc_names <- names(trait_pos)

    # Create trait label dataframe
    trait_lab <- mapply(function(pos, name) {
      chr_val <- if (is.numeric(pos[1])) sprintf("Chr%02s", pos[1]) else as.character(pos[1])
      mid_val <- if (length(pos) == 3) mean(as.numeric(pos[2:3]), na.rm = TRUE) else as.numeric(pos[2])
      data.frame(mid_pos = mid_val, loc_ns = name, chr = chr_val, stringsAsFactors = FALSE)
    }, trait_pos, loc_names, SIMPLIFY = FALSE)

    trait_lab_df <- do.call(rbind, trait_lab)

    # Add vertical dashed lines at each locus
    gg_obj <- gg_obj +
      ggplot2::geom_vline(data = trait_lab_df,
                          ggplot2::aes(xintercept = mid_pos),
                          linetype = "dashed", color = "black", linewidth = 1) +
      ggplot2::geom_text(data = trait_lab_df,
                         ggplot2::aes(x = mid_pos, label = loc_ns),
                         col = "black",
                         angle = 60,
                         hjust = 0,
                         y = length(unique(grp$x)) + label_offset,
                         size = text_size / (text_scale_fct * text_size),
                         inherit.aes = FALSE)
    return(gg_obj)
  }

  for (i in seq_len(nbatches)) {
    grp <- batches[[i]]

    plt <- ggplot2::ggplot(grp, ggplot2::aes(x = pos, y = x, fill = value, width = map_dist)) +
      ggplot2::geom_tile(lwd = 1, linetype = 1) +

      ggplot2::scale_fill_manual(values = ggplot2::alpha(col_mapping, alpha),
                                 label = col_labels,
                                 name = legend_title) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5))) +
      ggplot2::xlab('Chromosome') +

      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  labels = function(x) paste0(x / 1e6)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::geom_hline(yintercept = c(as.numeric(grp$x) + 0.5, 0.5),
                          col = 'white', lwd = 2.5) +
      ggplot2::coord_cartesian(clip = "off") +

      ggplot2::facet_grid(rows = NULL, cols = ggplot2::vars(chr),
                          space = "free", scales = "free_x", switch = 'x') +
      # ggplot2::facet_wrap(~ chr, nrow = 1, scales = "free_x",
      #                     strip.position = "bottom") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = text_size, face = 'bold'),
                     axis.ticks.length.y = ggplot2::unit(0.25, 'cm'),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_text(size = text_size, face = 'bold'),
                     axis.title.y = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     strip.placement = "outside",
                     plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5),
                     panel.background = ggplot2::element_rect(fill = panel_fill, colour = panel_col,
                                                              linewidth = 2, linetype = "solid"),
                     legend.text = ggplot2::element_text(size = text_size),
                     legend.title = ggplot2::element_text(size = text_size),
                     strip.background = ggplot2::element_rect(fill = "grey90", colour = "grey90",
                                                              linewidth = 0.5, linetype = "solid"),
                     strip.text = ggplot2::element_text(size = text_size))

    if (!is.null(trait_pos)) {
      if (!inherits(trait_pos, what = 'list')) {
        stop("Argument `trait_pos` must be a list object.")
      }

      plt <- annotate_loc(gg_obj = plt, trait_pos = trait_pos, grp = grp)
    }

    gg_plts[[i]] <- plt
  }

  if (pdf) {
    ggplot2::ggsave(filename = paste0(filename, ".pdf"),
                    plot = gridExtra::marrangeGrob(gg_plts, nrow = 1, ncol = 1),
                    device = "pdf",
                    units = "in",
                    width = width,
                    height = height)
  } else {
    return(gg_plts)
  }
}

