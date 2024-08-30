#' Read raw KASP results file (csv format) with one or multiple plates.
#' @param file A character value indicating the file name or path for the KASP
#' results file (csv format).
#' @param row_tags A character vector for the ordered row tags for the
#' components of the data in \code{file}.
#' @param spacing An integer value for specifying the number of empty rows between
#' data components in \code{file}.
#' @param data_type A character value indicating the file type to import;
#' currently supports either `raw` or `polished` KASP results file.
#'
#' @returns A list object of KASP results file for genotyping calls with FAM and HEX coordinates.
#'
#' @examples
#' # example code
#' library(panGenomeBreedr)
#' # Read raw bulked KASP data into R
#' \donttest{
#' path1 <-  system.file("extdata", "Genotyping_141.010_01.csv",
#'                       package = "panGenomeBreedr",
#'                       mustWork = TRUE)
#'
#' file1 <- read_kasp_csv(file = path1, data_type = 'raw')
#' # Get KASP genotyping data for plotting
#' kasp_dat <- file1$Data
#' }
#'
#' @export
#' @import utils

read_kasp_csv <- function(file,
                          row_tags = c('Statistics', 'DNA', 'SNPs','Scaling', 'Data'),
                          spacing = 2L,
                          data_type = c('raw', 'polished')) {

  # row_tags <- match.arg(row_tags)
  data_type <- match.arg(data_type)

  file <- file.path(file) # File path

  read_rows <- readLines(file, warn = FALSE)

  # Row tags for components of input data
  tags <- row_tags

  # Empty numeric vector to hold temporary results
  # Row numbers of tags in raw data file
  row_nums <- vector(mode = 'numeric', length = length(tags))

  # List components are data frames in raw kasp results file
  kasp_res <- vector(mode = 'list', length = length(tags) )
  names( row_nums) <- names(kasp_res) <- tags

  if (data_type == 'polished'){

    kasp_res[[5]] <- utils::read.csv(file, header = TRUE, na.strings = c("", NA))

  } else {

  # Get row numbers for each tag
  for (i in seq_len(length(tags))) {

    if (any(grepl(tags[i], read_rows)) == FALSE){

      stop(paste(tags[i], 'is not contained in the input file.'))

    }

    row_nums[i] <- grep(tags[i], read_rows)

  }

  for (d in seq_len(length(tags))) {

    # Read csv data of genotype calls and FAM/HEX coordinates
    dat <- utils::read.csv(file, skip = row_nums[d], header = TRUE,
                           na.strings = c("", NA))

    if (d == length(tags)) {

      # Subset data for each component
      aa <- dat

      # Remove rows with all NAs
      bb <- aa[rowSums(is.na(aa[, 2:ncol(aa)])) != (ncol(aa)-1),]

      # Remove columns with all NAs and pass data frame to kasp_res list object
      kasp_res[[d]] <- bb[, colSums(is.na(bb)) != nrow(bb)]

    } else {

      dd <- row_nums - (row_nums[d] + 1)

      # Subset data for each component
      aa <- dat[1:(dd[(d+1)] - spacing),]

      # Remove rows with all NAs
      bb <- aa[rowSums(is.na(aa[, 2:ncol(aa)])) != (ncol(aa)-1),]

      # Remove columns with all NAs
      kasp_res[[d]] <- bb[, colSums(is.na(bb)) != nrow(bb)]

    }

    }

  }

  return(kasp_res)
}


#' Get SNP or InDel alleles and possible genotypes from genotype calls in KASP data.
#' @param x A character vector of KASP genotype calls in one resaction plate.
#' @param sep A character used as separator for genotype calls, default is a
#' colon.
#' @param uncallable A character indicating `Uncallable` genotype calls, if present.
#' @param unused A character indicating `?` genotype calls, if present.
#' @param blank A character value indicating `No Template Controls (NTC)`
#' genotype calls.
#' @param others A character vector indicating other non-genotype calls in KASP
#' genotype calls, if present. These may include `'Missing', 'Bad', 'Dupe'`,
#' `'Over', 'Short'`.
#'
#'@examples
#'# example code
#' \donttest{
#' x <- panGenomeBreedr::kasp_dat$Call[1:96]
#' alleles <- get_alleles(x = x)
#' }
#'
#' @returns A list object with `length = 2` consisting of marker alleles and
#' possible genotypes in each KASP reaction plate.
#' @export


get_alleles <- function(x,
                        sep = ':',
                        blank = 'NTC',
                        uncallable = 'Uncallable',
                        unused = '?',
                        others = c('Missing', 'Bad', 'Dupe', 'Over', 'Short')
) {


  # Empty list object to hold results
  res <- vector(mode = 'list', length = 2)
  names(res) <- c('alleles', 'genotypes')

  # Create a character vector of all non-genotype calls
  non_geno_call <- c(blank, uncallable, unused, others)

  # Subset and sort genotype calls
  geno_uniq <- sort(unique(x[!x %in% non_geno_call]))

  alleles <- res[[1]] <- unique(unlist(strsplit(x = geno_uniq, split = sep)))

  homo1 <- paste(alleles[1], alleles[1], sep = sep) # homozygous genotype 1
  homo2 <- paste(alleles[2], alleles[2], sep = sep) # homozygous genotype 2
  het1 <- paste(alleles[1], alleles[2], sep = sep) # heterozygous genotype 1
  het2 <- paste(alleles[2], alleles[1], sep = sep) # heterozygous genotype 2

  res[[2]] <- c(homo1 = homo1, homo2 = homo2, het1 = het1, het2 = het2)


  return(res)
}


#' Generate pch characters for cluster plots of KASP genotype calls.
#' @param x A character vector of KASP genotype calls in one reaction plate.
#' @param sep A character used as separator for genotype calls, default is a
#' colon.
#' @param uncallable A character indicating `Uncallable` genotype calls, if present.
#' @param unused A character indicating `?` genotype calls, if present.
#' @param blank A character value indicating `No Template Controls (NTC)`
#' genotype calls.
#' @param others A character vector indicating other non-genotype calls in KASP
#' genotype calls, if present. These may include `'Missing', 'Bad', 'Dupe'`,
#' `'Over', 'Short'`.
#'
#'@examples
#'# example code
#' \donttest{
#' x <- panGenomeBreedr::kasp_dat$Call[1:96]
#' geno_pch <- kasp_pch(x = x)
#' }
#'
#' @returns A `96 x 1` data frame  of pch values for possible genotypes in
#' each KASP reaction plate.
#' @export
kasp_pch <- function(x,
                     sep = ':',
                     blank = 'NTC',
                     uncallable = 'Uncallable',
                     unused = '?',
                     others = c('Missing', 'Bad', 'Dupe', 'Over', 'Short')) {


  # Empty integer vector to hold results
  # res <- vector(mode = 'integer', length = nrow(x))

  # Get unique elements in geno call
  all_element <- unique(x)

  # Create a character vector of all possible non-genotype calls
  non_geno_call <- c(blank, uncallable, unused, others)

  # Subset and sort real genotype calls in KASP geno calls
  geno_uniq <- sort(unique(x[!x %in% non_geno_call]))

  # Combine geno_uniq with other non-genotype calls present in KASP geno call
  all_element <- c(geno_uniq, all_element[!all_element %in% geno_uniq])

  # # Combine geno_uniq and non_geno_call
  # all_calls <- c(geno_uniq, non_geno_call)
  #
  # all_element <- all_element[all_element %in% all_calls]

  if (length(geno_uniq) > 0) {

    pch_geno <- (21:(20 + length(geno_uniq))) # Assign pch values to geno_uniq
    pch_geno[pch_geno == 22] <- 24
    names(pch_geno) <- geno_uniq

  } else {

    pch_geno <- NULL

  }

  pch_non_geno <- c(22, 3, 4, 7:11) # pch for all possible non-geno calls
  names(pch_non_geno) <- c(blank, uncallable, unused, others)

  # Combine pch for real geno calls and non-geno calls present in KASP data
  # This is to be used as a loopup vector for pch values
  pch_vec <- c(pch_geno, pch_non_geno[names(pch_non_geno) %in% all_element])

  dat <- data.frame(pch = x) # Create a data frame for pch geno calls

  # Match geno calls to lookup vector for pch values for each geno call
  res <- data.frame(lapply(dat, function(i) pch_vec[i]))

  return(res)

}


#' Color-code KASP genotype calls based on LGC genomics colors for HEX and FAM.
#' @param x A data frame of KASP genotype calls for one or multiple plates.
#' @param subset A character value indicating the column name for taking subsets
#' of \code{x} for processing; default is the `MasterPlate`.
#' @param sep A character used as separator for genotype calls, default is a
#' colon.
#' @param geno_call A character indicating the column name used for genotype
#' calls in \code{x}.
#' @param uncallable A character indicating `Uncallable` genotype calls, if present.
#' @param unused A character indicating `?` genotype calls, if present.
#' @param blank A character value indicating `No Template Controls (NTC)`
#' genotype calls.
#' @param others A character vector indicating other non-genotype calls in KASP
#' genotype calls, if present. These may include `'Missing', 'Bad', 'Dupe'`,
#' `'Over', 'Short'`.
#'
#' @returns A list object with subset unit as component data frames.
#'
#' @examples
#' # example code
#' library(panGenomeBreedr)
#' \donttest{
#' dat1 <- kasp_color(x = panGenomeBreedr::kasp_dat,
#'                    subset = 'MasterPlate',
#'                    sep = ':',
#'                    geno_call = 'Call',
#'                    uncallable = 'Uncallable',
#'                    unused = '?',
#'                    blank = 'NTC')
#' }
#'
#' @details
#' This is an experimental function. The default values of some of the
#' arguments in the function are based on LGC Genomics conventions including
#' the color codes for FAM and HEX fluorescence. To make the plot color-blind
#' friendly, the heterozygotes have been coded as gold instead of green.
#' @export
#'
kasp_color <- function(x,
                       subset = 'MasterPlate',
                       sep = ':',
                       geno_call = 'Call',
                       uncallable = 'Uncallable',
                       unused = '?',
                       blank = 'NTC',
                       others = c('Missing', 'Bad', 'Dupe', 'Over', 'Short')) {

  # Get subset unit from data set
  plates <- x[, subset]
  plates_uniq <- unique(plates)

  # Number of units or plates
  nplates <- length(plates_uniq)

  # Empty list object to hold results
  res <- vector(mode = 'list', length = nplates)
  names(res) <- plates_uniq

  for (i in seq_len(nplates)) {
    # Subset each master plate
    master_plate <- x[plates == plates_uniq[i],]

    # create a color vector based on KASP geno calls
    Color <- master_plate[, geno_call]

    # Get alleles and possible genotypes using the `get_allele()` function
    alleles_geno <- get_alleles(x = Color,
                                sep = sep,
                                blank = blank,
                                uncallable = uncallable,
                                unused = unused,
                                others = others)

    alleles <- alleles_geno$alleles

    # Get pch values for all geno calls using the `kasp_pch()` function
    pch_geno <- kasp_pch(x = Color,
                         sep = sep,
                         blank = blank,
                         uncallable = uncallable,
                         unused = unused,
                         others = others)

    if (length(alleles) == 2) {

      homo1 <- alleles_geno$genotype[1] # homozygous genotype 1
      homo2 <- alleles_geno$genotype[2] # homozygous genotype 2
      homo12 <- c(homo1, homo2)
      het1 <- alleles_geno$genotype[3] # heterozygous genotype 1
      het2 <- alleles_geno$genotype[4] # heterozygous genotype 2

      # Subset homozygotes
      homs <- subset(master_plate, subset = Color == homo1 | Color == homo2)

      # Get the column index of the genotype calls
      call_col <- which(colnames(master_plate) == geno_call)

      FAM <- homs[which.max(homs$X),][, call_col] # FAM homozygote
      HEX <- homo12[which(homo12 != FAM)] # HEX homozygote

      # Recode pch symbols to be consistent with FAM and HEX genotypes
      pch_geno$pch[Color == HEX] <- 21
      pch_geno$pch[Color == FAM] <- 23
      pch_geno$pch[Color == het1 | Color == het2] <- 24

      # Assign colors based on LGC rules
      Color[Color == HEX] <- "red" # HEX homozygote
      Color[Color == FAM] <- "blue" # FAM homozygote
      Color[Color == het1] <- "yellow2" # Heterozygote
      Color[Color == het2] <- "yellow2" # Heterozygote

      # KASP FAM and HEX color coding for others
      Color[Color == uncallable] <- "darkmagenta"
      Color[Color == unused] <- "salmon4"
      Color[Color == blank] <- "black"
      Color[Color %in% others] <- "mintcream"


      master_plate <- cbind(master_plate, Color, pch_geno)

      res[[i]] <- master_plate


    } else if (length(alleles) == 1) {

      homo1 <- paste(alleles[1], alleles[1], sep = sep) # homozygous genotype 1
      Color[Color == homo1] <- "blue" # Homozygous 1

      # KASP FAM and HEX color coding for others
      Color[Color == uncallable] <- "darkmagenta"
      Color[Color == unused] <- "salmon4"
      Color[Color == blank] <- "black"
      Color[Color %in% others] <- "mintcream"

      master_plate <- cbind(master_plate, Color, pch_geno)

      res[[i]] <- master_plate

    } else if (length(alleles) == 0){

      message(cat(paste('Marker in Plate', plates_uniq[i], 'failed! \n', 'Check genotype calls.')))

      # KASP FAM and HEX color coding for others
      Color[Color == uncallable] <- "darkmagenta"
      Color[Color == unused] <- "salmon4"
      Color[Color == blank] <- "black"
      Color[Color %in% others] <- "mintcream"

      master_plate <- cbind(master_plate, Color, pch_geno)

      res[[i]] <- master_plate

    } else {

      message(cat(paste('Marker in Plate', plates_uniq[i], 'is multi-allelic! \n', 'Check genotype calls.')))

      # KASP FAM and HEX color coding for others
      Color[Color == uncallable] <- "darkmagenta"
      Color[Color == unused] <- "salmon4"
      Color[Color == blank] <- "black"
      Color[Color %in% others] <- "mintcream"

      master_plate <- cbind(master_plate, Color, pch_geno)

      res[[i]] <- master_plate
    }

  }

  return(res)

}

#' Normalize FAM and HEX fluorescence values between 0 and 1
#' @param x A numeric vector of FAM or HEX fluorescence values
#'
#' @returns A numeric vector of normalized values between 0 and 1
#'
#' @examples
#' # example code
#' library(panGenomeBreedr)
#' # Get Plate 1
#' dat1 <- panGenomeBreedr::kasp_dat[1:96,]
#' FAM_scaled <- scale_axis(dat1$X)
#'
#' @export
#'
scale_axis <- function(x) {

  norm <- (x - min(x))/(max(x) - min(x))

  return(norm)

}

#' Generate the prediction status of positive controls in a KASP assay, if present.
#' @param plate A data frame of KASP genotype calls for one plate.
#' @param geno_call A character indicating the column name used for genotype
#' calls in \code{plate}.
#' @param blank A character value indicating `No Template Controls (NTC)`
#' genotype calls.
#' @param Group_id A character value for the column ID indicating the predictions
#' of the positive controls in \code{plate}.
#' @param Group_unknown A character value representing unknown expected genotype status
#' for samples, if present. No genotype prediction can be made for such samples.
#' @returns A data frame with the prediction status of each sample added as a column.
#'
#' @details
#' The function recodes the prediction status of each sample as follows:
#' TRUE = prediction matches observed genotype call
#' FALSE = prediction does not match observed genotype call
#' Unknown = Either observed genotype call could not be made or expected genotype
#' could not be made prior to KASP genotyping or both.
#' Blank = NTC wells
#'
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#' # Get Plate 1
#'
#' dat1 <- panGenomeBreedr::beta_carotene[1:96,]
#' dat1 <- pred_status(plate = dat1,
#' geno_call = 'Call',
#' Group_id = 'Group',
#' blank = 'NTC',
#' Group_unknown = '?')
#' }
#'
#' @export
pred_status <- function(plate,
                        geno_call = 'Call',
                        Group_id = 'Group',
                        blank = 'NTC',
                        Group_unknown = NULL) {

  status <- NULL
  gg <- get_alleles(x = plate[, geno_call]) # Get possible genotypes in data

  df1 <- plate[plate[, geno_call] %in% gg$genotypes,] # Subset only samples with true calls
  df2 <- plate[!plate[, geno_call] %in% gg$genotypes,] # Subset samples with non-genotype calls

  if(!is.null(Group_unknown )) {

  # Add status for samples with true genotype calls
  df1$status <- ifelse(df1[, geno_call] == df1[, Group_id], 'True',
                       ifelse(df1[, geno_call] != df1[, Group_id] &
                                df1[, Group_id] == Group_unknown, 'Unknown', 'False'))

  # Add status for samples with non-genotype calls and NTC
  df2$status <- ifelse(df2[, geno_call] == blank, 'Blank', 'Unknown')

  } else {

    # Add status for samples with true genotype calls
    df1$status <- ifelse(df1[, geno_call] == df1[, Group_id], 'True', 'False')

    # Add status for samples with non-genotype calls and NTC
    df2$status <- ifelse(df2[, geno_call] == blank, 'Blank', 'Unknown')

}
  # rbind the two data frames df1 and df2
  df3 <- rbind(df1, df2)

  # Sort data in ascending order of chromosomes
  df3 <- df3[order(as.numeric(rownames(df3)), decreasing = FALSE),]

  return(df3)

}


#' Generate summary of predicition for positive controls in KASP genotype data,
#' if present
#' @param x A list object of KASP genotype calls processed by the `kasp_color()`
#' function.
#' @param snp_id A character value indicating the column name for SNP IDs
#' in \code{x}.
#' @param blank A character value indicating `No Template Controls (NTC)`
#' genotype calls.
#' @param Group_id A character value for the column ID indicating the predictions
#' of the positive controls in \code{x}.
#' @param geno_call A character value indicating the column name for KASP genotype
#' calls in \code{x}.
#' @param Group_unknown A character value representing unknown expected genotype status
#' for samples, if present. No genotype prediction can be made for such samples.
#'
#' @returns A list object with plates and prediction summary as components.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#' dat1 <- panGenomeBreedr::beta_carotene
#' dat1 <- kasp_color(x = beta_carotene,
#'                    subset = 'plates',
#'                    sep = ':',
#'                    geno_call = 'Call',
#'                    uncallable = 'Uncallable',
#'                    unused = '?',
#'                    blank = 'NTC')
#'
#' dat1 <- pred_summary(x = dat1,
#'                     snp_id = 'SNPID',
#'                     geno_call = 'Call',
#'                     Group_id = 'Group',
#'                     blank = 'NTC',
#'                     Group_unknown = '?')
#' dat1$summ
#' }
#'
#' @export

pred_summary <- function(x,
                         snp_id = 'SNPID',
                         Group_id = NULL,
                         blank = 'NTC',
                         Group_unknown = '?',
                         geno_call = 'Call') {

  if (!is.null(Group_id)) {

    # Get the number of plates or subset units in data input
    nplates <- length(x)

    # Get plate names
    plate_ns <- names(x)

    # Create an empty list object to hold ggplots
    res <- vector(mode = 'list', length = nplates)
    names(res) <- plate_ns

    df <- as.data.frame(matrix(data = NA, nrow = nplates, ncol = 5))
    colnames(df) <- c('plate', 'snp_id', 'false', 'true', 'unknown')

    for (i in seq_len(nplates)) {
      # Subset each plate

      plate <- x[[i]]

      df[i, 1] <- plate_ns[i]
      df[i, 2] <- plate[, snp_id][1]

      plate <- pred_status(plate = plate,
                           geno_call = geno_call,
                           Group_id = Group_id,
                           blank = blank,
                           Group_unknown = Group_unknown)

      fals <- length(plate$status[plate$status == 'False'])
      tru <- length(plate$status[plate$status == 'True'])
      unkn <- length(plate$status[plate$status == 'Unknown'])

      df[i, 3:5] <- c('False' = fals, 'True' = tru, 'Unknown' = unkn)

      res[[i]] <- plate

    }

  } else {

    stop("Provide value for the 'Group_id' argument.")

  }


  dat <- list(plates = res, summ = df)

  return(dat)

}

#' Make KASP marker genotyping QC plot.
#' @param x A list object of KASP genotype calls processed by the `kasp_color()`
#' function.
#' @param FAM A character indicating the column name for FAM fluorescence
#' coordinates in \code{x}.
#' @param HEX A character indicating the column name for HEX fluorescence
#' coordinates in \code{x}.
#' @param geno_call A character value indicating the column name for KASP genotype
#' calls in \code{x}.
#' @param color A character value indicating the column name for assigned colors
#' in \code{x}.
#' @param snp_id A character value indicating the column name for SNP IDs
#' in \code{x}.
#' @param blank A character value indicating `No Template Controls (NTC)`
#' genotype calls.
#' @param Group_id A character value for the column ID indicating the predictions
#' of the positive controls in \code{x}.
#' @param scale A logical value indicating whether to scale FAM and HEX axis
#' to values between 0 and 1.
#' @param pdf A logical value indicating whether to save plot as a pdf graphic
#' device when TRUE or output plot in R when FALSE.
#' @param width A numeric value for the width of pdf device.
#' @param height A numeric value for the height of pdf device.
#' @param filename A character value for path or file name for saving pdf.
#' @param expand_axis A numeric value to expand the axes for legend positioning.
#' @param legend.pos.x A numeric value representing the x coordinate for legend
#' placement.
#' @param legend.pos.y A numeric value representing the y coordinate for legend
#' placement.
#' @param legend.box A character value of either `horizontal` or `vertical`
#' legend placement.
#' @param legend.pos A character value for the position of legend; the
#' default value is `inside`.
#' @param alpha A numeric value between 0 and 1 for modifying the
#' opacity of colors.
#' @param text_size A numeric value for setting text size.
#' @param ... Other valid arguments that can be passed to ggplot2.
#'
#' @returns A graphic object or plot.
#'
#' @examples
#' # example code
#' library(panGenomeBreedr)
#' \donttest{
#' # Assign KASP colors to plates
#' dat1 <- kasp_color(x = panGenomeBreedr::kasp_dat,
#'                    subset = 'MasterPlate',
#'                    sep = ':',
#'                    geno_call = 'Call',
#'                    uncallable = 'Uncallable',
#'                    unused = '?',
#'                    blank = 'NTC')
#'
#' # KASP QC plot for Plate 12
#' kasp_qc_ggplot(x = dat1[12],
#'                     pdf = FALSE,
#'                     Group_id = 'Group',
#'                     scale = TRUE,
#'                     expand_axis = 0.6,
#'                     alpha = 0.5,
#'                     legend.pos.x = 0.6,
#'                     legend.pos.y = 0.8)
#' }
#'
#' @export
#' @import ggplot2
#' @import gridExtra

kasp_qc_ggplot <- function(x,
                           FAM = 'X',
                           HEX = 'Y',
                           geno_call = 'Call',
                           color = 'Color',
                           snp_id = 'SNPID',
                           blank = 'NTC',
                           Group_id = NULL,
                           scale = FALSE,
                           pdf = TRUE,
                           width = 6,
                           height = 6,
                           filename = 'kasp_qc',
                           expand_axis = 0.5,
                           legend.pos.x = 0.6,
                           legend.pos.y = 0.75,
                           legend.box = 'horizontal',
                           legend.pos = 'inside',
                           alpha = 0.5,
                           text_size = 12,
                           ...) {

  # Get the number of plates or subset units in data input
  nplates <- length(x)

  # Get plate names
  plate_ns <- names(x)

  # Create an empty list object to hold ggplots
  gg_plts <- vector(mode = 'list', length = nplates)
  names(gg_plts) <- plate_ns

  X <- NULL
  Y <- NULL
  Call <- NULL

  for (i in seq_len(nplates)) {

    plate <- x[[i]]

    snp_ns <- plate[, snp_id][1]

    title <- paste0('Plate: ', plate_ns[i], ' | ', 'SNP: ', snp_ns)

    # Generate plotting characters
    if (!is.null(Group_id)) {

      Group_pch <- plate[, Group_id] # Get predictions for positive controls
      Group_pch <- as.factor(Group_pch) # Convert to a factor
      aa <- levels(Group_pch) # Get factor levels

      if (length(aa) > 5) stop("Number of positive control groups exceeded!")

      aa <- c(aa[aa != blank], blank) # reorder factor levels
      Group_pch <- factor(Group_pch, levels = aa)

      # Create plotting characters based on the number of groups
      PCH <- vector(mode = 'integer', length = length(Group_pch))

      PCH[Group_pch != blank] <- as.integer(Group_pch[Group_pch != blank]) + 20

      PCH[PCH == 21] <- 25 # Replace all 21s with 25

      PCH[Group_pch == blank] <- 21 # Assign pch 21 to Blanks


      Group_pch <- as.character(Group_pch)
      Group_pch[Group_pch == blank] <- 'Blank'

      plate$Group_pch <- Group_pch


      # Sorting breaks for prediction labels
      grp_pch <- unique(PCH)
      grp <- unique(Group_pch)
      names(grp) <- grp_pch

      grp1 <-  sort(grp[grp != 'Blank'])
      grp2 <- sort(grp[grp == 'Blank'])

      grp <- c(grp1, grp2)
      grp_pch <- as.numeric(names(grp))

    } else {

      PCH <- 21 # If no positive controls, pch 21 for all samples

    }

    # Sorting breaks for the Observation labels
    Color <- plate[, 'Color'] # Subset color assignment column from input data

    # Get unique colors in input data
    cols <- unique(Color)

    Calls <- plate[, geno_call]

    # Replace NTC with Blank and ? with Unused
    Calls[Calls == blank] <- 'Blank'
    Calls[Calls == '?'] <- 'Unused'

    plate$Call <- Calls

    # Get unique calls in input data
    calls <- unique(Calls)

    # Assign unique colors to unique calls
    names(calls) <- cols

    # Subset and sort genotype calls
    calls2 <- sort(calls[calls != 'Blank' & calls != 'Uncallable' & calls != 'Unused'])

    # Subset and sort non-genotype calls
    calls3 <- sort(calls[calls == 'Blank' | calls == 'Uncallable' | calls == 'Unused'])

    # Re-combine unique calls
    calls <- c(calls2, calls3)
    cols <- names(calls)

    # Scale FAM and HEX values to 0 and 1s -- recommended
    if (scale == TRUE) {

      plate$X <- scale_axis(plate[, FAM]) # Scale FAM fluorescence coordinates from data input
      plate$Y <- scale_axis(plate[, HEX]) # Scale HEX fluorescence coordinates from data input

    } else {

      plate$X <- plate[, FAM] # FAM fluorescence coordinates from data input
      plate$Y <- plate[, HEX] # HEX fluorescence coordinates from data input

    }

    # Set x and y axes limits
    x_max <- round(max(plate[, FAM]), 1) + expand_axis
    y_max <- round(max(plate[, HEX]), 1) + expand_axis

    xy_max <- c(x_max, y_max)

    # Set axis limits to comparable values
    axis_max <- xy_max[which.max(xy_max)]

    # Plot title
    title <- paste0('Plate: ', plate_ns[i])

    # Plotting
    if (!is.null(Group_id)) {


        plt <- ggplot2::ggplot(plate, ggplot2::aes(x = X,
                                                   y = Y,
                                                   fill = Call,
                                                   shape = Group_pch)) +
          ggplot2::geom_point(size = 4) +
          ggplot2::scale_shape_manual(values = grp_pch,
                                      breaks = grp,
                                      name = 'Prediction') +

          ggplot2::guides(shape = ggplot2::guide_legend(order = 2)) +

          ggplot2::scale_fill_manual(values = ggplot2::alpha(cols, alpha),
                                     breaks = calls,
                                     name = 'Observation') +
          ggplot2::guides(fill = ggplot2::guide_legend(order = 0, override.aes = list(shape = 21)))

    } else {

      plt <- ggplot2::ggplot(plate, ggplot2::aes(x = X,
                                                 y = Y,
                                                 fill = Call)) +
        ggplot2::geom_point(size = 4, pch = 21) +
        ggplot2::scale_fill_manual(values = ggplot2::alpha(cols, alpha),
                                   breaks = calls,
                                   name = 'Observation')
    }

    plt <- plt + ggplot2::theme_classic() + ggplot2::xlim(c(0, axis_max)) +
      ggplot2::ylim(c(0, axis_max)) +
      ggplot2::labs(title = paste(title, '\n', 'SNP: ', snp_ns),
                    x = "FAM fluorescence", y = "HEX fluorescence") +
      ggplot2::theme(panel.background = ggplot2::element_rect(color = "black",
                                                              linewidth = 1),
                     axis.line = ggplot2::element_line(linewidth = 0.2),
                     plot.title = ggplot2::element_text(hjust = 0.5,
                                                        face = 'bold',
                                                        size = text_size),
                     axis.text = ggplot2::element_text(size  = text_size,
                                                       face = 'bold',
                                                       color = 'black'),
                     axis.title = ggplot2::element_text(size  = text_size + 4,
                                                        face = 'bold',
                                                        color = 'black'),
                     legend.text = ggplot2::element_text(size  = text_size,
                                                         color = 'black'),
                     legend.title = ggplot2::element_text(size  = text_size,
                                                          color = 'black'),
                     legend.key = ggplot2::element_blank(),
                     legend.position = legend.pos,
                     legend.position.inside = c(legend.pos.x, legend.pos.y),
                     legend.box = legend.box
      )

    gg_plts[[i]] <- plt

  }

  if (pdf == TRUE) {

    ggplot2::ggsave(filename = paste0(filename, ".pdf"),
                    plot = gridExtra::marrangeGrob(gg_plts, nrow = 1, ncol = 1),
                    device = "pdf", path = file.path(getwd()),
                    units = "in", width = width, height = height)

  } else {

    return(gg_plts)

  }

}


#' Make KASP marker genotyping QC plot overlaid with predicitons.
#' @param x A list object of KASP genotype calls processed by the `kasp_color()`
#' function.
#' @param FAM A character indicating the column name for FAM fluorescence
#' coordinates in \code{x}.
#' @param HEX A character indicating the column name for HEX fluorescence
#' coordinates in \code{x}.
#' @param geno_call A character value indicating the column name for KASP genotype
#' calls in \code{x}.
#' @param color A character value indicating the column name for assigned colors
#' in \code{x}.
#' @param PCH A character value indicating the column name for assigned PCH symbols
#' in \code{x}.
#' @param snp_id A character value indicating the column name for SNP IDs
#' in \code{x}.
#' @param blank A character value indicating `No Template Controls (NTC)`
#' genotype calls.
#' @param uncallable A character indicating `Uncallable` genotype calls, if present.
#' @param unused A character indicating `?` genotype calls, if present.
#' @param others A character vector indicating other non-genotype calls in KASP
#' genotype calls, if present. These may include `'Missing', 'Bad', 'Dupe'`,
#' `'Over', 'Short'`.
#' @param Group_id A character value for the column ID indicating the predictions
#' of the positive controls in \code{x}.
#' @param Group_unknown A character value representing unknown expected genotype status
#' for samples, if present. No genotype prediction could be made for such samples.
#' @param pred_cols A named character vector of length = 4 of colors to be used for
#' the prediction legend for positive controls, if present.
#' @param scale A logical value indicating whether to scale FAM and HEX axis
#' to values between 0 and 1.
#' @param pdf A logical value indicating whether to save plot as a pdf graphic
#' device when TRUE or output plot in R when FALSE.
#' @param width A numeric value for the width of pdf device.
#' @param height A numeric value for the height of pdf device.
#' @param filename A character value for path or file name for saving pdf.
#' @param expand_axis A numeric value to expand the axes for legend positioning.
#' @param legend.pos.x A numeric value representing the x coordinate for legend
#' placement.
#' @param legend.pos.y A numeric value representing the y coordinate for legend
#' placement.
#' @param legend.box A character value of either `horizontal` or `vertical`
#' legend placement.
#' @param legend.pos A character value for the position of legend; the
#' default value is `inside`.
#' @param alpha A numeric value between 0 and 1 for modifying the
#' opacity of colors.
#' @param text_size A numeric value for setting text size.
#' @param ... Other valid arguments that can be passed to ggplot2.
#'
#' @returns A graphic object or plot.
#'
#' @examples
#' # example code
#' library(panGenomeBreedr)
#' \donttest{
#' # Assign KASP colors to plates
#' dat1 <- kasp_color(x = panGenomeBreedr::kasp_dat,
#'                    subset = 'MasterPlate',
#'                    sep = ':',
#'                    geno_call = 'Call',
#'                    uncallable = 'Uncallable',
#'                    unused = '?',
#'                    blank = 'NTC')
#'
#' # KASP QC plot for Plate 12
#' kasp_qc_ggplot2(x = dat1[12],
#'                     pdf = FALSE,
#'                     Group_id = 'Group',
#'                     Group_unknown = NULL,
#'                     scale = TRUE,
#'                     expand_axis = 0.6,
#'                     alpha = 0.5,
#'                     legend.pos.x = 0.6,
#'                     legend.pos.y = 0.8)
#' }
#'
#' @export
#' @import ggplot2
#' @import gridExtra


kasp_qc_ggplot2 <- function(x,
                            FAM = 'X',
                            HEX = 'Y',
                            geno_call = 'Call',
                            color = 'Color',
                            PCH = 'pch',
                            snp_id = 'SNPID',
                            blank = 'NTC',
                            uncallable = 'Uncallable',
                            unused = '?',
                            others = c('Missing', 'Bad', 'Dupe', 'Over', 'Short'),
                            Group_id = NULL,
                            Group_unknown = '?',
                            pred_cols = c('Blank' = 'black', 'False' = 'red',
                                          'True' = 'blue', 'Unknown' = 'yellow2'),
                            scale = FALSE,
                            pdf = TRUE,
                            width = 6,
                            height = 6,
                            filename = 'kasp_qc',
                            expand_axis = 0.5,
                            legend.pos.x = 0.6,
                            legend.pos.y = 0.75,
                            legend.box = 'horizontal',
                            legend.pos = 'inside',
                            alpha = 0.9,
                            text_size = 12,
                            ...) {


  # Get the number of plates or subset units in data input
  nplates <- length(x)

  # Get plate names
  plate_ns <- names(x)

  # Create an empty list object to hold ggplots
  gg_plts <- vector(mode = 'list', length = nplates)
  names(gg_plts) <- plate_ns

  X <- NULL
  Y <- NULL
  Call <- NULL
  status <- NULL
  pred_cols <- pred_cols

  for (i in seq_len(nplates)) {

    plate <- x[[i]]

    snp_ns <- plate[, snp_id][i] # snp id

    title <- paste0('Plate: ', plate_ns[i], ' | ', 'SNP: ', snp_ns)

    # Extract geno calls
    Calls <- plate[, geno_call]

    # Get unique calls in input data
    calls <- unique(Calls)

    # Create a character vector of all non-genotype calls
    non_geno_call <- c(blank, uncallable, unused, others)

    # Sorting breaks for the Observation labels
    PCH_dat <- plate[, PCH] # Subset color assignment column from input data

    # Get unique colors in input data
    pch_uniq <- unique(PCH_dat)

    names(pch_uniq) <- calls
    pch_uniq <- c(sort(pch_uniq[!names(pch_uniq) %in% non_geno_call]),
                  sort(pch_uniq[names(pch_uniq) %in% non_geno_call]))

    # pch breaks
    pch_brks <- names(pch_uniq)
    pch_brks[pch_brks == blank] <- 'Blank'
    pch_brks[pch_brks == unused] <- 'Unused'
    names(pch_uniq) <- pch_brks


    # Sorting breaks for the Observation labels
    Color <- plate[, color] # Subset color assignment column from input data

    # Get unique colors in input data
    cols <- unique(Color)
    names(cols) <- calls

    cols <- c(sort(cols[!names(cols) %in% non_geno_call]),
              sort(cols[names(cols) %in% non_geno_call]))

    # col breaks
    col_brks <- names(cols)
    col_brks[col_brks == blank] <- 'Blank'
    col_brks[col_brks == unused] <- 'Unused'
    names(cols) <- col_brks


    # Generate plotting characters
    if (!is.null(Group_id)) {

      grp <- plate[, Group_id]

      if (length(unique(grp)) > 5) stop("Number of positive control groups exceeded!")

      plate <- pred_status(plate = plate,
                           geno_call = geno_call,
                           Group_id = Group_id,
                           blank = blank,
                           Group_unknown = Group_unknown)

      status_uniq <- sort(unique(plate$status))

      #pred_cols <- pred_cols[names(pred_cols) %in% status_uniq]

    }



    # Replace '?' and NTC with Blank and Unused for labeling
    plate[, geno_call][plate[, geno_call] == blank] <- 'Blank'
    plate[, geno_call][plate[, geno_call] == unused] <- 'Unused'


    # Scale FAM and HEX values to 0 and 1s -- recommended
    if (scale == TRUE) {

      plate$X <- scale_axis(plate[, FAM]) # Scale FAM fluorescence coordinates from data input
      plate$Y <- scale_axis(plate[, HEX]) # Scale HEX fluorescence coordinates from data input

    } else {

      plate$X <- plate[, FAM] # FAM fluorescence coordinates from data input
      plate$Y <- plate[, HEX] # HEX fluorescence coordinates from data input

    }

    # Set x and y axes limits
    x_max <- round(max(plate[, FAM]), 1) + expand_axis
    y_max <- round(max(plate[, HEX]), 1) + expand_axis

    xy_max <- c(x_max, y_max)

    # Set axis limits to comparable values
    axis_max <- xy_max[which.max(xy_max)]


    # Plotting
    if (!is.null(Group_id)) {


      plt <- ggplot2::ggplot(plate, ggplot2::aes(x = X,
                                                 y = Y,
                                                 fill = status,
                                                 shape = Call)) +
        ggplot2::geom_point(size = 4) +
        ggplot2::scale_shape_manual(values = pch_uniq,
                                    breaks = pch_brks,
                                    name = 'Observation') +

        ggplot2::guides(shape = ggplot2::guide_legend(order = 0)) +

        ggplot2::scale_fill_manual(values = ggplot2::alpha(pred_cols, alpha),
                                   breaks = status_uniq,
                                   name = 'Prediction') +
        ggplot2::guides(fill = ggplot2::guide_legend(order = 2, override.aes = list(shape = 21)))

    } else {

      plt <- ggplot2::ggplot(plate, ggplot2::aes(x = X,
                                                 y = Y,
                                                 fill = Call)) +
        ggplot2::geom_point(size = 4, pch = 21) +
        ggplot2::scale_fill_manual(values = ggplot2::alpha(cols, alpha),
                                   breaks = col_brks,
                                   name = 'Observation')
    }


    plt <- plt + ggplot2::theme_classic() + ggplot2::xlim(c(0, axis_max)) +
      ggplot2::ylim(c(0, axis_max)) +
      ggplot2::labs(title = paste(title, '\n', 'SNP: ', snp_ns),
                    x = "FAM fluorescence", y = "HEX fluorescence") +
      ggplot2::theme(panel.background = ggplot2::element_rect(color = "black",
                                                              linewidth = 1),
                     axis.line = ggplot2::element_line(linewidth = 0.2),
                     plot.title = ggplot2::element_text(hjust = 0.5,
                                                        face = 'bold',
                                                        size = text_size),
                     axis.text = ggplot2::element_text(size  = text_size,
                                                       face = 'bold',
                                                       color = 'black'),
                     axis.title = ggplot2::element_text(size  = text_size + 4,
                                                        face = 'bold',
                                                        color = 'black'),
                     legend.text = ggplot2::element_text(size  = text_size,
                                                         color = 'black'),
                     legend.title = ggplot2::element_text(size  = text_size,
                                                          color = 'black'),
                     legend.key = ggplot2::element_blank(),
                     legend.position = legend.pos,
                     legend.position.inside = c(legend.pos.x, legend.pos.y),
                     legend.box = legend.box
      )

    gg_plts[[i]] <- plt

  }


  if (pdf == TRUE) {

    ggplot2::ggsave(filename = paste0(filename, ".pdf"),
                    plot = gridExtra::marrangeGrob(gg_plts, nrow = 1, ncol = 1),
                    device = "pdf", path = file.path(getwd()),
                    units = "in", width = width, height = height)

  } else {

    return(gg_plts)

  }

}


#' Plot kasp genotyping plate layout.
#' @param x A list object of KASP genotype calls processed by the `kasp_color()`
#' function.
#' @param well A character value representing the column name for genotyping
#' plate wells.
#' @param color A character value indicating the column name of assigned colors
#' in \code{x}.
#' @param geno_call A character indicating the column name used for genotype
#' @param snp_id A character value indicating the column name for SNP IDs
#' in \code{x}.
#' @param pdf A logical value indicating whether to save plot as a pdf graphic
#' @param width A numeric value for the width of pdf device.
#' @param height A numeric value for the height of pdf device.
#' @param filename A character value for path or file name for saving pdf.
#' @param text_size A numeric value for text size in plot output.
#' @param ... Other valid arguments that can be passed to ggplot2.
#'
#' @returns A ggplot graphical output of plate layout.
#'
#' @examples
#' # example code
#' library(panGenomeBreedr)
#' \donttest{
#' # Assign KASP colors to plates
#' dat1 <- kasp_color(x = panGenomeBreedr::kasp_dat,
#'                    subset = 'MasterPlate',
#'                    sep = ':',
#'                    geno_call = 'Call',
#'                    uncallable = 'Uncallable',
#'                    unused = '?',
#'                    blank = 'NTC')
#' # Plot Plate 12 to see sample arrangement
#' plot_plate(dat1[12], pdf = FALSE)
#' }
#'
#' @export
#' @import ggplot2
#'
plot_plate <- function(x,
                       well = 'MasterWell',
                       color = 'Color',
                       geno_call = 'Call',
                       snp_id = 'SNPID',
                       pdf = TRUE,
                       width = 8,
                       height = 5,
                       filename = 'plate_layout',
                       text_size = 12,
                       ...){

  # Get the number of plates or subset units in data input
  nplates <- length(x)

  # Get plate names
  plate_ns <- names(x)

  gg_plts <- vector(mode = 'list', length = nplates)
  names(gg_plts) <- plate_ns

  V1 <- NULL
  V2 <- NULL

  gg_plate <- function() {

    # Make plot
    plt <- ggplot2::ggplot(pcr_plate, ggplot2::aes(x = V2, y = V1,
                                                   fill = Call)) +
      ggplot2::geom_tile(col = 'grey50',
                         alpha = 0,
                         linewidth = 0.2) +
      ggplot2::geom_point(size = text_size, alpha = 0.5, pch = 21) +
      ggplot2::theme_classic() +
      ggplot2::scale_fill_manual(values = myCol,
                                 breaks = calls,
                                 name = snp_ns) +
      ggplot2::scale_x_discrete(position = "top") +
      ggplot2::scale_y_discrete(limits = rev(levels(pcr_plate$V1))) +
      ggplot2::labs(title = title) +

      ggplot2::theme(axis.line = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5,
                                                        face = 'bold',
                                                        size = text_size),
                     axis.text = ggplot2::element_text(size  = text_size,
                                                       face = 'bold',
                                                       color = 'black')
      )

    return(plt)

  }


  for (i in seq_len(nplates)) {

  # Subset plate
  plate <- x[[i]]
  # Get well column for wells
  MasterWell <- plate[, well]
  well_freq <- table(MasterWell)

  # Stop function if well IDs are not unique.
  if( any(well_freq > 1)) {

    stop("Plate well IDs are not unique! Check the plate layout.")

  }

  # Get Color column for genotype calls
  Color <- plate[, color]
  myCol <- unique(Color)

  # Get geno call
  Call <- plate[, geno_call]

  # Replace NTC with Blank and ? with Unused
  Call[Call == 'NTC'] <- 'Blank'
  Call[Call == '?'] <- 'Unused'

  calls <- unique(Call)
  names(calls) <- myCol

  # Subset and sort genotype calls
  calls2 <- sort(calls[calls != 'Blank' & calls != 'Uncallable' & calls != 'Unused'])

  # Subset and sort non-genotype calls
  calls3 <- sort(calls[calls == 'Blank' | calls == 'Uncallable' | calls == 'Unused'])

  # Re-combine unique calls
  calls <- c(calls2, calls3)
  myCol <- names(calls)


  # Get SNP name
  snp_ns <- plate[, snp_id][1]

  # Plot title
  title <- paste0('Plate: ', plate_ns[i])

  # Split alphanumeric well IDs into alphabet and numeric
  pcr_plate <- strsplit(MasterWell, split = "(?<=[a-zA-Z])(?=[0-9])",
                        perl = TRUE)

  # Make a data frame from split well IDs
  pcr_plate <- t(as.data.frame(pcr_plate))
  pcr_plate <- as.data.frame(cbind(pcr_plate,
                                   MasterWell = MasterWell,
                                   Call = Call,
                                   Color = Color))

  # Convert Column 1 and 2 into a data frame
  pcr_plate$V1 <- as.factor(pcr_plate$V1)
  pcr_plate$V2 <- as.factor(pcr_plate$V2)
  # pcr_plate$Call <- as.factor(pcr_plate$Call)

  gg_plts[[i]] <- gg_plate()

  }

  if (pdf == TRUE) {

    ggplot2::ggsave(filename = paste0(filename, ".pdf"),
           plot = gridExtra::marrangeGrob(gg_plts, nrow = 1, ncol = 1),
           device = "pdf", path = file.path(getwd()),
           units = "in", width = width, height = height)

  } else {

    return(gg_plts)

  }

}




