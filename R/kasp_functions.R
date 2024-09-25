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
#' @param data_type A character value indicating the data source; either `kasp` or
#' `agriplex`.
#'
#'@examples
#'# example code
#' \donttest{
#' x <- panGenomeBreedr::kasp_dat$Call[1:96]
#' alleles <- get_alleles(x = x, data_type = 'kasp')
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
                        others = c('Missing', 'Bad', 'Dupe', 'Over', 'Short'),
                        data_type = c('kasp', 'agriplex')
) {

  data_type <- match.arg(data_type) # Match arguments

  # Empty list object to hold results
  res <- vector(mode = 'list', length = 2)
  names(res) <- c('alleles', 'genotypes')

  # Create a character vector of all non-genotype calls
  non_geno_call <- c(blank, uncallable, unused, others)

  # Subset and sort genotype calls
  geno_uniq <- sort(unique(x[!x %in% non_geno_call]))

  alleles <- res[[1]] <- unique(unlist(strsplit(x = geno_uniq, split = sep)))

  if (length(alleles) == 2) {

  if (data_type == 'kasp') {

    homo1 <- paste(alleles[1], alleles[1], sep = sep) # homozygous genotype 1
    homo2 <- paste(alleles[2], alleles[2], sep = sep) # homozygous genotype 2
    het1 <- paste(alleles[1], alleles[2], sep = sep) # heterozygous genotype 1
    het2 <- paste(alleles[2], alleles[1], sep = sep) # heterozygous genotype 2

  } else {

    homo1 <- paste(alleles[1]) # homozygous genotype 1
    homo2 <- paste(alleles[2]) # homozygous genotype 2
    het1 <- paste(alleles[1], alleles[2], sep = sep) # heterozygous genotype 1
    het2 <- paste(alleles[2], alleles[1], sep = sep) # heterozygous genotype 2

  }

  } else if (length(alleles == 1)) {

    if (data_type == 'kasp') {

      homo1 <- homo2 <- paste(alleles[1], alleles[1], sep = sep) # homozygous genotype 1
      het1 <- het2 <- paste(alleles[1], alleles[1], sep = sep) # heterozygous genotype 1

    } else {

      homo1 <- homo2 <- paste(alleles[1]) # homozygous genotype 1
      het1 <- het2 <- paste(alleles[1], alleles[1], sep = sep) # heterozygous genotype 1

    }

  } else if (is.null(alleles)) {

    if (data_type == 'kasp') {

      homo1 <- homo2 <- het1 <- het2 <- NULL

    } else {

      homo1 <- homo2 <- het1 <- het2 <- NULL

    }

  }

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
                                others = others,
                                data_type = 'kasp')

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
#' @param Group_unknown A character value representing unknown genotype status
#' for samples, if present. No genotype prediction can be made for such samples.
#' @returns A data frame with the prediction status of each sample added as a column.
#'
#' @details
#' The function recodes the prediction status of each sample as follows:
#' TRUE = prediction matches observed genotype call
#' FALSE = prediction does not match observed genotype call
#' Unverified = Either observed genotype call could not be made or expected genotype
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
  gg <- get_alleles(x = plate[, geno_call], data_type = 'kasp') # Get possible genotypes in data

  df1 <- plate[plate[, geno_call] %in% gg$genotypes,] # Subset only samples with true calls
  df2 <- plate[!plate[, geno_call] %in% gg$genotypes,] # Subset samples with non-genotype calls

  if(!is.null(Group_unknown )) {

  # Add status for samples with true genotype calls
  df1$status <- ifelse(df1[, geno_call] == df1[, Group_id], 'True',
                       ifelse(df1[, geno_call] != df1[, Group_id] &
                                df1[, Group_id] == Group_unknown, 'Unverified', 'False'))

  # Add status for samples with non-genotype calls and NTC
  df2$status <- ifelse(df2[, geno_call] == blank, 'Blank', 'Unverified')

  } else {

    # Add status for samples with true genotype calls
    df1$status <- ifelse(df1[, geno_call] == df1[, Group_id], 'True', 'False')

    # Add status for samples with non-genotype calls and NTC
    df2$status <- ifelse(df2[, geno_call] == blank, 'Blank', 'Unverified')

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
#' @param Group_unknown A character value representing unverified expected genotype status
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
    colnames(df) <- c('plate', 'snp_id', 'false', 'true', 'unverified')

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
      unkn <- length(plate$status[plate$status == 'Unverified'])

      df[i, 3:5] <- c('False' = fals, 'True' = tru, 'Unverified' = unkn)

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
                    device = "pdf",
                    units = "in",
                    width = width,
                    height = height)

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
#' @param Group_unknown A character value representing unknown genotype status for
#' samples, if present. No genotype prediction could be made for such samples.
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
                                          'True' = 'blue', 'Unverified' = 'yellow2'),
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
                    device = "pdf",
                    units = "in",
                    width = width,
                    height = height)

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
           device = "pdf",
           units = "in", width = width, height = height)

  } else {

    return(gg_plts)

  }

}


#' Design KASP markers based on causal variants.
#' @param vcf_file Path to the vcf file containing identified variants from
#' re-sequence data.
#' @param marker_ID Designated name of variant for marker design. Name must be
#' contained in `vcf_file`.
#' @param chr A character value representation the chromosome description in the
#' `genome_file`. Providing this helps to save memory in R.
#' @param genome_file Path to reference genome file in fasta format, either
#' compressed (.gz) or decompressed.
#' @param plot_draw A logical value indicating whether to plot 100 bp upstream
#' and downstream KASP sequence alignment to reference genome.
#' @param plot_file Path to save the sequence alignment if `plot_draw = TRUE`.
#' @param region_name A n optional character value for assigned region name.
#' @param vcf_geno_code A character vector for genotype coding for samples
#' in imported `vcf_file`. The order of the genotype codes in this vector matters:
#' `c('homo_alt_allele', 'heterozygous', 'homo_ref_allele', 'missing_genotype_call')`.
#' @param maf A numeric value between 0 and 1 representing the minor allele
#' frequency for variant subset.
#'
#' @examples
#' \donttest{
#'# Example to design a KASP marker on a substitution variant
#' path <- tempdir()
#' path1 <- "https://raw.githubusercontent.com/awkena/panGB/main/Chr02.fa.gz"
#' path2 <-  system.file("extdata", "Sobic.002G302700_SNP_snpeff.vcf",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#' ma1 <- kasp_marker_design(vcf_file = path2,
#'                           genome_file = path1,
#'                           marker_ID = "SNP_Chr02_69200443",
#'                           chr = "Chr02",
#'                           plot_draw = TRUE,
#'                           plot_file = path,
#'                           region_name = "ma1",
#'                           maf = 0.05)
#'
#' # View marker alignment output from temp folder
#' path3 <- file.path(path, list.files(path = path, "alignment_"))
#' # system(paste0('open "', path3, '"'))
#'
#' on.exit(unlink(path))
#' }

#'
#' @details
#' This function provides the intertek sequence to be used for marker development
#' for the selected casual variants. It provides all the information of allele and
#' location to fill Intertek form.It needs a vcf file of variants calls, and a genome
#' sequence of the target crop in fasta format.
#'
#'
#' @returns A data frame containing all information required for KASP marker
#' design and a DNA sequence alignment to the reference genome.
#'
#' @importFrom Biostrings fasta.index readDNAStringSet DNAStringSet
#' @importFrom GenomicRanges GRanges
#' @import VariantAnnotation
#' @import msa
#' @importFrom IRanges IRanges
#' @importFrom stats na.omit
#'
#' @export

kasp_marker_design <- function(vcf_file, #path and file name of tbi index
                               marker_ID, # what variant want to design the marker for
                               chr = NULL, # Chromosome with the variant
                               genome_file, #genome file path
                               plot_draw = TRUE, #Do you want a plot drawn? Yes or No
                               plot_file = tempdir(),
                               region_name = 'loc_1', #assigned region name
                               vcf_geno_code = c('1|1', '0|1', '0|0', '.|.'),
                               maf = 0.05) # Minor allele frequency
{

  # Extracting genomic region from vcf file
  vcf_gene <- VariantAnnotation::readVcf(file = vcf_file)

  if (vcf_gene@fixed@nrows >= 1) {

    #converting to dataframe the variant names
    variant_table <- as.data.frame(names(vcf_gene))

    #add chromsome
    variant_table$chrom <- as.character(vcf_gene@rowRanges@seqnames)

    #adding position
    variant_table$pos <- vcf_gene@rowRanges@ranges@start

    # Extracting reference allele to column
    variant_table$reference <- as.character(VariantAnnotation::ref(vcf_gene))

    # Extracting alternate allele to column
    variant_table$alternate <- as.character(unlist(VariantAnnotation::alt(vcf_gene)))

    # trying type of variation
    # variant_table$type <- NA

    tryCatch({

      variant_table$type <-  ifelse(VariantAnnotation::isInsertion(vcf_gene), "Insertion",
                                    ifelse(VariantAnnotation::isSubstitution(vcf_gene), "Substitution",
                                           ifelse(VariantAnnotation::isDeletion(vcf_gene), "Deletion", "Others")))

    })

    # getting allele frequency
    geno <- as.data.frame(VariantAnnotation::geno(vcf_gene))
    geno <- geno[geno$group_name == 'GT', -c(1:2)] # getting only the alleles

    #converting to numeric to calculate MAF
    geno <- ifelse(geno == vcf_geno_code[1], 2,
                   ifelse(geno == vcf_geno_code[2], 1,
                          ifelse(geno == vcf_geno_code[3], 0, NA)))

    # converting all to numeric
    tryCatch({
      suppressWarnings(geno <- apply(geno, MARGIN = 2, as.numeric))
    })

    # calculating frequencies
    variant_table$MAF <- round(rowSums(geno, na.rm = TRUE)/(ncol(geno)*2), digits = 3)

    # Renaming column
    names(variant_table)[1] <- c("ID")
  }

  # Reading genome by chromosome subset or whole genome
  # Indexing for genome file
  dd <-  Biostrings::fasta.index(genome_file, seqtype = "DNA")

  # Read only the sequence for the chromosome with the variant if provided,
  # else read the whole genome sequence
  if (!is.null(chr) && chr %in% dd$desc) {

    indx <-  which(dd$desc == chr)

    genome <-  Biostrings::readDNAStringSet(dd[indx, ])

  } else {

    genome <- Biostrings::readDNAStringSet(genome_file)

  }

  # genome <- Biostrings::readDNAStringSet(genome_file)

  #### extracting information from variant table
  variant <- variant_table[variant_table$ID == marker_ID,] #variant information
  pos <- as.numeric(variant$pos) #position of variant
  type <- variant$type #type of variant
  w_var <- nchar(variant$reference) # length of variant
  chromosome <- variant$chrom

  # snp variants only for adjacent highly polymorphic regions
  variants_snp_only <- variant_table[variant_table$type == 'Substitution',] # getting SNPs

  ###### Setting boundaries for markers
  ####### Upstream
  if (variant$type == 'Deletion') {
    up_bound_end <- pos # end of upstream boundary
    up_bound_start <- up_bound_end-99 #start of upstream boundary
    down_bound_start <- pos+w_var #start of downstream boundary
    down_bound_end <- down_bound_start+99 # end of downstream boundary

  } else if (variant$type == 'Substitution') {
    up_bound_end <- pos-1 # end of upstream boundary
    up_bound_start <- up_bound_end-99 #start of upstream boundary
    down_bound_start <- pos+1 #start of downstream boundary
    down_bound_end <- down_bound_start+99 # end of downstream boundary

  } else if (variant$type == 'Insertion') {
    up_bound_end <- pos # end of upstream boundary
    up_bound_start <- up_bound_end-99 #start of upstream boundary
    down_bound_start <- pos+1 #start of downstream boundary
    down_bound_end <- down_bound_start+99 # end of downstream boundary
  }

  # Process upstream and downstream variants
  up_down_variant <- function(variant_snp_dat,
                              MAF = maf,
                              start,
                              end,
                              type = c('up', 'down')) {

    type <- match.arg(type)

    # variants dataframe
    variants_100bp <- data.frame(
      ID = NA, # marker id
      reference = NA, #reference allele
      alternate = NA, #alternate allele
      pos = NA, # position
      type = NA, #type of variant
      MAF = NA, # minor allele frequency
      ontology = NA, #ontology term
      effect = NA, # effect
      gene_name = NA, #gene name
      one_letter = NA # one letter code
    )

    # variants
    if (type == 'up') {

      variants_100bp <- variant_snp_dat[variant_snp_dat$pos >= start &
                                          variant_snp_dat$pos < end,]

    } else {

      variants_100bp <- variant_snp_dat[variant_snp_dat$pos > start &
                                          variant_snp_dat$pos <= end,]
    }

    # Filtering by MAF
    variants_100bp <- variants_100bp[variants_100bp$MAF > MAF,] #MAF filtering

    #ordering by MAF
    variants_100bp <- variants_100bp[order(variants_100bp$MAF, decreasing = TRUE),]

    #selecting the two variants with higher MAF
    variants_100bp <- stats::na.omit(variants_100bp[c(1,2),])

    #ordering by pos
    variants_100bp <- variants_100bp[order(variants_100bp$pos, decreasing = FALSE),]

    return(variants_100bp)

  }

  variants_100bp_up <- up_down_variant(variant_snp_dat = variants_snp_only,
                                       MAF = maf,
                                       start = up_bound_start,
                                       end = up_bound_end,
                                       type = 'up')
  # What happens if it returns no data?
  variants_100bp_down <- up_down_variant(variant_snp_dat = variants_snp_only,
                                         MAF = maf,
                                         start = down_bound_start,
                                         end = down_bound_end,
                                         type = 'down')

  # Creating iupacode table for substitutions
  iupacode <- data.frame(one_letter = c('R', 'Y','S','W','K','M', 'R', 'Y','S','W','K','M'),
                         allele_comb1 = c('AG','CT','GC','AT','GT','AC','GA','TC','CG','TA','TG','CA'))

  # Creating collapse column of alleles for substitution
  variants_100bp_down$collapse <- paste(variants_100bp_down$reference,
                                        variants_100bp_down$alternate, sep = "")
  variants_100bp_up$collapse <- paste(variants_100bp_up$reference,
                                      variants_100bp_up$alternate, sep = "")

  # Substitution collapse allele for the the one letter iupac code
  iupac_conv <- function(iupacode = iupacode, variants_100bp) {

    if (nrow(variants_100bp) > 0) {
      i <- 1
      while (i <= nrow(variants_100bp)) {
        variants_100bp$one_letter[i] <-
          iupacode[rowSums(iupacode == variants_100bp$collapse[i]) > 0, ]$one_letter
        i <- i+1
      }
    }

    return(variants_100bp)
  }

  #upstream
  variants_100bp_up <- iupac_conv(iupacode = iupacode,
                                  variants_100bp = variants_100bp_up)
  #downstream
  variants_100bp_down <- iupac_conv(iupacode = iupacode,
                                    variants_100bp = variants_100bp_down)

  ### Getting the actual upstream and downstream sequences

  get_actual_seq <- function(variants_100bp,
                             genome = genome,
                             bound_start,
                             bound_end) {

    ####### start stop list
    start_stop_list <- c(bound_start, variants_100bp$pos - 1, bound_end)

    # creating empty sequence
    stream_sequence <- Biostrings::DNAStringSet()

    j <- 1
    for ( i in 1:(length(start_stop_list)-1)) {

      #adjust start
      if (i == 1) {

        tmp_start <- start_stop_list[i] # first iteration

      } else {

        tmp_start <- start_stop_list[i] + 1 # adding a base

      }

      #adjust stop
      if (i == length(start_stop_list)-1) {

        tmp_stop <- start_stop_list[i+1] # if is the last iteration

      } else {

        tmp_stop <- start_stop_list[i+1]-1 # removing a base

      }

      # Extracting genomic region from reference genome
      region <- GenomicRanges::GRanges(chromosome, IRanges::IRanges(tmp_start, tmp_stop))

      # Creating empty base to add
      base_to_add <- Biostrings::DNAStringSet()

      if (j <= nrow(variants_100bp)) {

        #getting one letter code for the polymorphism upstream
        base_to_add <- Biostrings::DNAStringSet(variants_100bp$one_letter[j])

        j <- j+1 #counter

      }

      # getting sequence for the region
      sequence <- Biostrings::DNAStringSet(Biostrings::getSeq(genome, region))

      # merging sequence and base to add
      tmp_seq <- c(sequence, base_to_add)

      # unlisting and making new sequence
      new_seq <- Biostrings::DNAStringSet(unlist(tmp_seq))

      #saving to upstream sequence, merging with upstream sequence
      stream_sequence <- Biostrings::DNAStringSet(unlist(c(stream_sequence, new_seq)))
    }

    return(stream_sequence)
  }

  # Upstream
  upstream_sequence <- get_actual_seq(variants_100bp = variants_100bp_up,
                                      genome = genome,
                                      bound_start = up_bound_start,
                                      bound_end = up_bound_end)

  # Downdtream
  downstream_sequence <- get_actual_seq(variants_100bp = variants_100bp_down,
                                        genome = genome,
                                        bound_start = down_bound_start,
                                        bound_end = down_bound_end)



  ## Create output sequence annotated
  if (variant$type == 'Substitution') {

    reference_allele <- variant$reference # reference allele
    alternate_allele <- variant$alternate # alternate allele

    # intertek formatted sequence
    intertek_sequence <- paste(upstream_sequence, "[",
                               variant$reference,"/", variant$alternate, "]",
                               downstream_sequence, sep = "")

  } else if (variant$type == 'Deletion') {

    #variant without anchoring base
    varc <- substr(variant$reference, start = 2, stop = nchar(variant$reference))

    #reference allele
    reference_allele <- varc

    #alternate allele
    alternate_allele <- '-'

    # intertek formatted sequence
    intertek_sequence <- paste(upstream_sequence,"[", varc, "/-", "]",
                               downstream_sequence, sep = "")

  } else if (variant$type == 'Insertion') {

    #variant without anchoring base
    varc <- substr(variant$alternate, start = 2, stop = nchar(variant$alternate))

    #reference allele
    reference_allele <- '-'

    #alternate allele
    alternate_allele <- varc

    # intertek formatted sequence
    intertek_sequence <- paste(upstream_sequence, "[-/", varc, "]",
                               downstream_sequence, sep = "")

  }

  # Creating dataframe for export
  #Empty dataframe
  marker_data <- data.frame(SNP_Name = variant$ID,
                            SNP = variant$type,
                            Marker_Name = region_name,
                            Chromosome = chromosome,
                            Chromosome_Position = variant$pos,
                            Sequence = intertek_sequence,
                            ReferenceAllele = reference_allele,
                            AlternativeAllele = alternate_allele)

  # assign("marker_data", marker_data, envir = .GlobalEnv) # creating dataframe

  # Plot alignment of upstream and downstream sequences to the reference seq
  if (plot_draw == TRUE) {

    # getting genomic range to extract
    region_control <- GenomicRanges::GRanges(chromosome, IRanges::IRanges(up_bound_start, down_bound_end))

    # getting reference sequence
    sequence_control <- VariantAnnotation::getSeq(genome, region_control)

    # merging upstream and downstream sequence
    tmp <- c(upstream_sequence, downstream_sequence)

    # merging all sequences for plotting
    sequences <- c(sequence_control, tmp)

    # converting to DNAstring set
    sequences <- Biostrings::DNAStringSet(sequences)

    # adding sequence names
    names(sequences) <- c('reference', 'upstream', 'downstream')

    # alignment
    alg <- msa::msaClustalOmega(sequences, order = 'input')

    #converting to dnastringset
    alignment <- Biostrings::DNAStringSet(alg)

    # Extract the aligned sequences
    aligned_strings <- as.character(alignment)

    # Re-implementing alignment plot in ggplot
    y <- NULL
    group <- NULL
    x <- NULL

    x <- data.frame(
      x = 1,
      y = c(0.5, 0.4, 0.3),
      seq = aligned_strings,
      group = c('reference', 'upstream', 'downstream'))

    # Define the colors you want for each group
    group_colors <- c("reference" = "black",
                      "upstream" = "red",
                      "downstream" = "blue")

    pp <- ggplot2::ggplot(data = x, ggplot2::aes(x = x, y = y, label = seq,
                                                 color = group)) +
      ggplot2::geom_text(size = 4.5, hjust = 'outward', family = 'mono',
                         key_glyph = 'rect') +

      ggplot2::labs(color = paste("Alignment for marker", marker_ID, sep = " ")) +
      ggplot2::scale_x_continuous(breaks = c(-2, 0, 2), expand = c(0, 2)) +
      ggplot2::scale_y_continuous(breaks = c(-2, 0, 2), expand = c(0, 2)) +
      ggplot2::scale_color_manual(values = group_colors,
                                  breaks = c('reference', 'upstream', 'downstream'),
                                  label = c('Reference', 'Upstream', 'Downstream')) +
      ggplot2::guides(color = ggplot2::guide_legend(label.position = "right",
                                                    title.position = 'top',
                                                    title.hjust = 0.5)) +
      ggplot2::theme(legend.title = ggplot2::element_text(size  = 12,
                                                          color = 'black'),
                     panel.background = ggplot2::element_rect('white'),
                     axis.line = ggplot2::element_blank(),
                     axis.text = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size  = 12,
                                                         color = 'black'),
                     legend.position = 'inside',
                     legend.position.inside = c(0.5, 0.6),
                     legend.direction = 'horizontal')

    # Save ggplot object as a pdf device
    ggplot2::ggsave(filename = paste0('alignment_', marker_ID, '.pdf'),
                    plot = pp,
                    device = "pdf",
                    path = file.path(plot_file),
                    units = "in",
                    width = 24,
                    height = 9)
  }

  return(marker_data)

}

#' Get a summary of the number of samples per 96-well plate in a multi-plate KASP
#' assay.
#' @param x A data frame of KASP genotype calls for one or multiple plates.
#' @param subset A character value indicating the column name for taking subsets
#' of \code{x} for processing; default is the `plates`.
#' @param snp_id A character value indicating the column name for SNP IDs
#' in \code{x}.
#' @param plate_id A character value indicating the column name for master plate
#'  having the same samples in \code{x}.
#'
#' @returns A list object with plates and a summary of number of samples
#' per plate as components.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#' dat1 <- panGenomeBreedr::beta_carotene
#' dat1 <- nsamples_plate(x = dat1,
#'                      subset = 'plates',
#'                      snp_id = 'SNPID',
#'                      plate_id = 'MasterPlate'
#'                     )$summ
#' }
#'
#' @export

nsamples_plate <- function(x,
                           subset = 'plates',
                           snp_id = 'SNPID',
                           plate_id = 'MasterPlate') {

  # Get subset unit from data set
  plates <- x[, subset]
  plates_uniq <- unique(plates)

  # Number of units or plates
  nplates <- length(plates_uniq)

  # Create an empty list object to hold ggplots
  res <- vector(mode = 'list', length = nplates)
  names(res) <- plates_uniq

  df <- as.data.frame(matrix(data = NA, nrow = nplates, ncol = 4))
  colnames(df) <- c( 'MasterPlate', 'SNPID' ,'plate', 'nobs')

  for (i in seq_len(nplates)) {

    # Subset each plate
    plate <- x[plates == plates_uniq[i],]

    df[i, 1] <- plate[, plate_id][1]
    df[i, 2] <- plate[, snp_id][1]
    df[i, 3] <- plates_uniq[i]
    df[i, 4] <- nrow(plate)

    res[[i]] <- plate

  }

  dat <- list(plates = res, summ = df)

  return(dat)

}


#' Reshape KASP data to wide format for same samples genotyped with multiple KASP
#' markers.
#' @param x A data frame of KASP genotype calls for one or multiple plates.
#' @param subset A character value indicating the column name for master plate
#'  having the same samples in \code{x}.
#' @param snp_id A character value indicating the column name for SNP IDs
#' in \code{x}.
#' @param idvar A character value indicating the column name for unique subject
#' IDs of samples in \code{x}.
#' @param geno_call A character value indicating the column name for KASP genotype
#' calls in \code{x}.
#' @param blank A character value indicating `No Template Controls (NTC)`
#' genotype calls.
#'
#' @returns A list object with reshaped data for master plates that have the same
#' subject IDs as the components. The components are data frames with Column 1 as
#' the subject IDs and the rest of the columns as the number of KASP markers assayed
#' for each master plate.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#' dat1 <- panGenomeBreedr::beta_carotene
#' plate_wide <- kasp_reshape_wide(x = dat1,
#'                                 subset = 'MasterPlate',
#'                                 snp_id = 'SNPID',
#'                                 geno_call = 'Call',
#'                                 idvar = "SubjectID",
#'                                 blank = 'NTC')
#' }
#'
#' @importFrom stats reshape
#'
#' @export
kasp_reshape_wide <- function(x,
                              subset = 'MasterPlate',
                              snp_id = 'SNPID',
                              geno_call = 'Call',
                              idvar = "SubjectID",
                              blank = 'NTC') {


  # Get subset unit from data set
  plates <- x[, subset]
  plates_uniq <- unique(plates)

  # Number of units or plates
  nplates <- length(plates_uniq)

  # Create an empty list object to hold ggplots
  res <- vector(mode = 'list', length = nplates)
  names(res) <- plates_uniq

  for (i in seq_len(nplates)) {

    # Subset master plate that has the same samples
    plate <- x[plates == plates_uniq[i], ]

    # Get relevant columns
    retain_column <- c(geno_call, idvar, snp_id)

    # Use only columns that are relevant
    plate <- plate[names(plate) %in% retain_column]

    # Remove blanks before conversion
    plate <- plate[plate[, geno_call] != blank,]

    ## long to wide (direction = "wide") requires idvar and timevar at a minimum
    plate_wide <- stats::reshape(plate,
                                 direction = "wide",
                                 idvar = "SubjectID",
                                 timevar = "SNPID")

    # Use original SNPIDs as names for the SNPs after conversion
    colnames(plate_wide)[colnames(plate_wide) != idvar] <- unique(plate[, snp_id])

    res[[i]] <- plate_wide

  }

  return(res)

}

#' Process reshaped KASP genotype data for heatmap plotting
#' @param x A data frame of genotype calls.
#' @param sample_id A string representing the column name of unique sample IDs.
#' @param marker_start An integer indicating the column index for the start of SNP calls.
#' @param kasp_map A data frame consisting of marker chromosome number and positions.
#' @param map_snp_id A character value indicating the column name for SNP IDs
#' in \code{kasp_map}.
#' @param chr A character value for the column name for the chromosome number in
#' \code{kasp_map}.
#' @param chr_pos A character value for the column name for the chromosome positions
#'  in \code{kasp_map}.
#'
#' @examples
#' \donttest{
#' # example code
#' library(panGenomeBreedr)
#' dat1 <- panGenomeBreedr::beta_carotene
#' # Reshape KASP data for each master plate for the beta_carotene data
#' plate_wide <- kasp_reshape_wide(x = dat1,
#'                                 subset = 'MasterPlate',
#'                                 snp_id = 'SNPID',
#'                                 geno_call = 'Call',
#'                                 idvar = "SubjectID",
#'                                 blank = 'NTC')
#'
#' # Get Master Plate 1
#' plate1 <- plate_wide$`SE-24-1088_P01_d1`
#'
#' # Generate a map for the beta_carotene KASP data
#' kasp_map <- data.frame(SNPID = unique(beta_carotene$SNPID),
#'                        SNPID_2 = paste0('snpSB', 1:4, '|', c(800, 803, 804, 805)),
#'                        chr = c(1:4),
#'                        pos = c(800, 803, 804, 805))
#'
#' # Process Plate1 to re-order SNPs based on chrom. and position
#' proc_plate1 <- proc_kasp(x = plate1,
#'                          kasp_map = kasp_map,
#'                          map_snp_id = "SNPID",
#'                          sample_id = "SubjectID",
#'                          marker_start = 2,
#'                          chr = 'chr',
#'                          chr_pos = 'pos')
#' }
#' @returns A list consisting of re-ordered SNPs and oredered map based on chromosome
#' number and chromosome positions.
#'
#' @details
#' This function is experimental and should only be used with reshaped KASP data.
#' The data processing involves the following steps:
#' 1. Transposing the data
#' 2. Matching marker names in KASP map file with names in KASP genotype file
#' 3. Column bind matched map for markers to the transposed data
#' 4. Sort data in ascending order of chromosome numbers
#' 5. Split data by chromosomes
#' 6. Sort data in ascending order of physical positions per chromosome
#' 7. Row bind sorted data together
#' 8. Re-transpose data to make markers columns and samples as rows for plotting.
#'
#'@export

proc_kasp <- function(x,
                      sample_id = "SubjectID",
                      marker_start = 2,
                      kasp_map,
                      map_snp_id,
                      chr = 'chr',
                      chr_pos = 'pos') {

  if(missing(kasp_map)) stop("Provide map for KASP markers used!")
  if(missing(map_snp_id)) stop("Provide column name of SNP IDs in map file!")

  # Get only genotype call data without subject/sample IDs
  if (marker_start > 1) {

    geno_mat <- x[, -c(1:marker_start-1)]

    # Set row names of genotype data to original names
    # Row names must be unique for each row

    # if there are no duplicates in the sample IDs
    if (!any(duplicated(x[, sample_id]))) {

      rownames(geno_mat) <- x[, sample_id]

    } else {

      # if there are duplicates in the sample IDs
      # Get their index numbers
      dups_ind <- which(duplicated(x[, sample_id]))

      # Show a warning message when there duplicates in sample IDs
      warnings(paste("There are duplicates in sample IDs:", dups_ind))

      # Subset sample IDs
      new_sample_id <- x[, sample_id]

      # Rename duplicated sample IDs
      new_sample_id[dups_ind] <- paste0(new_sample_id[dups_ind], dups_ind)

      rownames(geno_mat) <- new_sample_id

    }

  } else if (marker_start  == 1) {

    geno_mat <- x

  }

  # Get the number of columns in kasp map file
  ncolumns <- ncol(kasp_map)

  # Transpose imported data to make markers as rows and samples as columns
  geno_mat <- t(geno_mat)

  # Function to match SNP IDs in map file and geno file

    new_snpids <- data.frame(ids = rownames(geno_mat))

    # Match SNPID in geno_mat and kasp_map files
    df1 <- merge(x = new_snpids,
                 y = kasp_map,
                 by.x = 'ids',
                 by.y = map_snp_id,
                 sort = FALSE)

  # Column bind chr number and position to snp data
  geno_mat <- cbind(df1, geno_mat)

  # Sort data in ascending order of chromosomes
  geno_mat <- geno_mat[order(geno_mat[, chr], decreasing = FALSE),]

  # Sort data in ascending order of physical positions per chromosome
  grps <- split(geno_mat, geno_mat[, chr]) # Split genotype data into chr batches

  # Function to order marker positions
  ord_pos <- function(x) {
    x[order(x[, chr_pos], decreasing = FALSE),]
  }

  geno_mat <- lapply(grps, FUN = ord_pos)

  # Row-bind sorted data for all chromosomes
  geno_mat <- do.call(rbind, geno_mat)
  rownames(geno_mat) <- geno_mat[, 1]

  df2 <- geno_mat[, seq_len(ncolumns)]

  # geno_mat <- na.omit(geno_mat) # Remove markers with missing snp calls
  # rownames(geno_mat) <- paste0('S', geno_mat$chr,'_', geno_mat$pos)

  geno_mat <- t(geno_mat[,-seq_len(ncolumns)]) # Transpose data to make markers columns

  # if (marker_start > 1) {
  #
  #   geno_mat <- cbind(x[, c(1:marker_start-1)], geno_mat) # Add sample meta data
  #
  # }

  res <- list(ordered_geno = geno_mat, ordered_map = df2)
  return(res)
}


#' Identify SNP loci with potential genotype call errors.
#' @param x A data frame of genotype calls where all columns are SNPs and samples
#' as rows. Row names are unique sample names.
#' @param rp_row An integer or character value indication the row index or name
#' of the recurrent or Parent 1.
#' @param dp_row An integer or character value indicating the row index or name
#' of the donor or Parent 2.
#' @param sep A character used as separator for genotype calls, default is a
#' colon.
#' @param uncallable A character indicating `Uncallable` genotype calls, if present.
#' @param unused A character indicating `?` genotype calls, if present.
#' @param blank A character value indicating `No Template Controls (NTC)`
#' genotype calls.
#' @param others A character vector indicating other non-genotype calls in KASP
#' genotype calls, if present. These may include `'Missing', 'Bad', 'Dupe'`,
#' `'Over', 'Short'`.
#' @param data_type A character value indicating the data source; either `kasp` or
#' `agriplex`.
#'
#' @examples
#' # example code
#' \donttest{
#' library(panGenomeBreedr)
#' # Reshape KASP data for each master plate for the beta_carotene data
#' dat1 <- panGenomeBreedr::beta_carotene
#' plate_wide <- kasp_reshape_wide(x = dat1,
#'                                 subset = 'MasterPlate',
#'                                 snp_id = 'SNPID',
#'                                 geno_call = 'Call',
#'                                 idvar = "SubjectID",
#'                                 blank = 'NTC')
#'
#' # Get Master Plate 1
#' plate1 <- plate_wide$`SE-24-1088_P01_d1`
#'
#' # Check for genotype call  error for each SNP
#' geno_mat <- geno_error(x = plate1[,-1],
#'                   rp_row = 1,
#'                   dp_row = 7,
#'                   sep = ':',
#'                   data_type = 'kasp')
#' }
#'
#' @returns A list object with the following components:
#' 1) column index of loci with genotype call error
#' 2) data frame of loci with genotype call errors if present.
#'
#' @export
geno_error <- function(x,
                       rp_row,
                       dp_row,
                       sep = ':',
                       blank = 'NTC',
                       uncallable = 'Uncallable',
                       unused = '?',
                       others = c('Missing', 'Bad', 'Dupe', 'Over', 'Short'),
                       data_type = c('kasp', 'agriplex')
) {

  if (missing(rp_row)) stop("Parent 1 row index number is missing!")
  if (missing(dp_row)) stop("Parent 2 row index number is missing!")

  data_type <- match.arg(data_type)

  # Create a character vector of all non-genotype calls
  non_geno_call <- c(blank, uncallable, unused, others)

  # Hold column index of loci with errors
  col_index <- c()

  for (i in seq_len(ncol(x))) {

    snp <- as.vector(unname(x[, i]))

    # Subset and sort genotype calls
    snp <- snp[!snp %in% non_geno_call]

    # Get Alleles for recurrent and donor parents
    rp <- as.vector(unlist(x[rp_row,]))[i] # Recurrent parent allele
    dp <- as.vector(unlist(x[dp_row,]))[i] # Donor parent allele

    par_geno <- c(rp, dp) # Parents genotype

    # Expected genotype vector
    exp_geno <- get_alleles(x = par_geno, sep = sep, data_type = data_type)$genotypes

    #exp_geno <- par_alleles$genotypes # Expected genotypes

    geno_present <- unique(na.omit(snp)) %in% exp_geno

    if (!anyNA(par_geno) && rp == dp && any(!geno_present)) {

      col_index[i] <- i

    } else if (!anyNA(par_geno) && rp != dp && any(!geno_present)) {

      col_index[i] <- i

    }

  }

  col_index <- as.vector(na.omit(col_index))

  # Get SNPs with suspected errors and no errors
  if (length(col_index) >= 1) {

    geno_err <- as.data.frame(x[, col_index]) # Suspected error
    colnames(geno_err) <- colnames(x)[col_index]
    rownames(geno_err) <- rownames(x)

    geno_good <- as.data.frame(x[, -col_index]) # no error
    colnames(geno_good) <- colnames(x)[-col_index]
    rownames(geno_good) <- rownames(x)

  } else {

    geno_err <- NULL
    geno_good <- x

  }

  res <- list(col_index = col_index,
              geno_err = geno_err,
              geno_good = geno_good)

  return(res)

}


#' Convert processed KASP data to numeric genotypes
#' @param x A data frame of n rows of genotypes and p rows of SNPs; output from
#' `proc_kasp()` function.
#' @param rp_row An integer or character value indicating the row index or name
#' of the recurrent or Parent 1.
#' @param dp_row An integer or character value indicating the row index or name
#' of donor or Parent 2.
#' @param sep A character used as separator for genotype calls, default is a
#' colon.
#' @param uncallable A character indicating `Uncallable` genotype calls, if present.
#' @param unused A character indicating `?` genotype calls, if present.
#' @param blank A character value indicating `No Template Controls (NTC)`
#' genotype calls.
#' @param others A character vector indicating other non-genotype calls in KASP
#' genotype calls, if present. These may include `'Missing', 'Bad', 'Dupe'`,
#' `'Over', 'Short'`.
#' @param data_type A character value indicating the data source; either `kasp` or
#' `agriplex`.
#'
#' @returns A data frame of numeric codes for KASP genotype calls.
#'
#' @examples
#'
#' \donttest{
#'  # example code
#' library(panGenomeBreedr)
#'
#' # Reshape KASP data for each master plate for the beta_carotene data
#' dat1 <- panGenomeBreedr::beta_carotene
#' plate_wide <- kasp_reshape_wide(x = dat1,
#'                                 subset = 'MasterPlate',
#'                                 snp_id = 'SNPID',
#'                                 geno_call = 'Call',
#'                                 idvar = "SubjectID",
#'                                 blank = 'NTC')
#'
#' # Get Master Plate 1
#' plate1 <- plate_wide$`SE-24-1088_P01_d1`
#'
#' # Convert to numeric format for plotting
#' num_geno <- kasp_numeric(x = plate1[,-1],
#'                         rp_row = 1,
#'                         dp_row = 7,
#'                         data_type = 'kasp')
#'
#' }
#'
#' @details
#' Re-coded as 1 if homozygous for Parent 1 allele; 0 if homozygous for Parent
#' 2 allele, and 0.5 if heterozygous. If any of the parents SNP call is missing,
#' it's coded as as NA. If the Parents 1 and 2 genotypes are the same for any snp,
#' it's coded as monomorphic. Loci with progeny genotype calls different from RP
#' and DP are coded as -2.
#'
#'@export
kasp_numeric <- function(x,
                        rp_row,
                        dp_row,
                        sep = ':',
                        blank = 'NTC',
                        uncallable = 'Uncallable',
                        unused = '?',
                        others = c('Missing', 'Bad', 'Dupe', 'Over', 'Short'),
                        data_type = c('kasp', 'agriplex')
                        ) {

  if (missing(rp_row)) stop("Parent 1 row index number is missing!")
  if (missing(dp_row)) stop("Parent 2 row index number is missing!")

  data_type <- match.arg(data_type)

  # Create a character vector of all non-genotype calls
  non_geno_call <- c(blank, uncallable, unused, others)

  # Get Alleles for recurrent and donor parents
  rp <- as.vector(unlist(x[rp_row,])) # Recurrent parent allele
  dp <- as.vector(unlist(x[dp_row,])) # Donor parent allele

  # Create an empty data frame to hold numeric format of genotypes
  num_recode <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  colnames(num_recode) <- colnames(x)
  rownames(num_recode) <- rownames(x)

  for(i in seq_len(ncol(x))) {

    snp <- as.vector(unname(x[, i]))

    # Subset and sort genotype calls
    # snp <- snp[!snp %in% non_geno_call]

    par_geno <- c(rp[i], dp[i])

    # Parent genotype vector
    exp_geno <- get_alleles(x = par_geno, sep = sep, data_type = data_type)$genotypes

    # exp_geno <- par_alleles$genotypes # Expected genotypes

    # Code loci as NA if any or all parent genotype call is NA
    if (anyNA(par_geno) || all(is.na(par_geno))) {

      num_recode[,i] <- NA

    } else if (length(unique(na.omit(snp))) == 1) {

      num_recode[,i] <- -1 # if monomorphic code as -1

      # Code loci with progeny genotype calls different from RP and DP as -2
    } else if (!anyNA(par_geno) && rp[i] == dp[i] && length(unique(na.omit(snp))) > 1 ) {

      num_recode[,i] <- ifelse(snp == rp[i], -1, -2)

    } else {
      # Recode as 1 if homozygous for RP allele; 0 if homozygous for DP allele
      # 0.5 if heterozygous, -2 if other
      het1 <- exp_geno[3]
      het2 <- exp_geno[4]
      num_recode[,i] <- ifelse(snp == rp[i], 1, ifelse(snp == dp[i], 0, ifelse(snp == het1 | snp == het2, 0.5, -2)))

    }

  }

  return(num_recode)

}


#' Parse marker names with a common pattern containing chromosome numbers and
#' positions into a map file.
#' @param x A character vector containing the original marker names to be parsed.
#' @param sep A character value that serves as a unique separator between chromosome
#' positions and other components of the marker name; default value is an underscore.
#' @param prefix A character value that represents a common pattern in marker
#' names that precedes the chromome number; the default value is `S`.
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

  # Throw up errors if separator and prefix are not common to all marker names
  if (all(grepl(sep, x)) == FALSE) stop('The separator is not common to all marker names!')
  if (all(grepl(prefix, x)) == FALSE) stop('The prefix is not common to all marker names!')

  # Parse marker names to extract chromosome numbers and physical positions
  df <- t(as.data.frame(strsplit(x, sep)))
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
