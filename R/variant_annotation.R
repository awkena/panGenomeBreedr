#' Extract putative causal variants within a candidate gene from a tabix-indexed
#' snpEff annotated VCF file.
#' @param cand_gene_id A character value specifying the candidate gene ID.
#' @param gff_path A character value indicating the path to the GFF file,
#' including the complete file name.
#' @param vcf_dir A character value indicating the path to directory containing
#' snpEff annotated VCF files.
#' @param vcf_file A character value indicating the file name for snpEff
#' annotated VCF file including the .vcf.gz extension.
#' @param output_path A character value indicating the path to directory for
#' saving extracted variants.
#' @param outfile_suffix A character value indicating the file name to be used for
#' saving extracted variants.
#'
#' @returns A VCF file containing extracted variants for \code{cand_gene_id}.
#'
#' @details
#' This wrapper function operates on tabix-indexed snpEff annotated VCF files.
#' However, if a tabix-indexed VCF file is not available, it can create one from
#' the inputted VCF file.
#'
#' The file names of snpEff annotated VCF files are expected to consist of three
#' components: a common prefix, chromosome tag and a common suffix.
#'
#'
#' @examples
#' # example code
#' \donttest{
#' library(panGenomeBreedr)
#'
#' # Work from the tempdir
#' vcf_dir <- tempdir()
#'
#' # Google drive link to gff3 file
#' flink1 <- "https://drive.google.com/file/d/1XjYyJ2JLywbbniIU6oUIIxAmEBKfmHpz/view?usp=sharing"
#'
#' # Download gff3 file to tempdir()
#' gff3 <- folder_download_gd(drive_link = flink1,
#'                            output_path = vcf_dir,
#'                            is.folder = FALSE)
#'
#' # Google drive link to indel snpEff annotated vcf file on Chr05
#' flink2 <- "https://drive.google.com/file/d/1LiOeDsfIwbsCuHbw9rCJ1FLOZqICTrfs/view?usp=sharing"
#'
#' # Download indel snpEff annotated vcf file to tempdir()
#' vcf_file_indel <- folder_download_gd(drive_link = flink2,
#'                                      output_path = vcf_dir,
#'                                      is.folder = FALSE)
#' # View downloaded files in tempdir
#' # list.files(vcf_dir)
#'
#' # InDel variant extraction for lgs1 (Sobic.005G213600)
#' extract_variant(cand_gene_id = 'Sobic.005G213600',
#'                 gff_path = gff3,
#'                 vcf_dir = vcf_dir,
#'                 vcf_file = basename(vcf_file_indel),
#'                 output_path = vcf_dir,
#'                 outfile_suffix = 'lgs_variants_indel')
#'
#' # Clean tempdir after variant extraction
#' # unlink(vcf_dir, recursive = TRUE)
#' }
#'
#' @export
#'
extract_variant <- function(cand_gene_id,
                            gff_path,
                            vcf_dir,
                            vcf_file,
                            output_path = tempdir(),
                            outfile_suffix = 'variants') {

  # Conditional logic for Bioconductor dependencies
  required_pkgs <- c("rtracklayer", "GenomicRanges", "Rsamtools")

  for (pkg in required_pkgs) {

    if (!requireNamespace(pkg, quietly = TRUE)) {

      stop(sprintf("The '%s' package is required for this function.
                   Please install it via BiocManager::install('%s').", pkg, pkg),
           call. = FALSE)

    }

  }


  # Load GFF file using the rtracklayer package
  # gff_path <- file.path(gff_path)
  gff <- rtracklayer::import(gff_path)

  # Extract coordinates for the candidate gene
  # If candidate gene is in gff file, extract coordinates
  if (cand_gene_id %in% gff$Name) {

    candidate_gene <- gff[gff$type == "gene" & grepl(cand_gene_id, gff$Name), ]
    gene_coords <- as.data.frame(candidate_gene)[, c("seqnames", "start", "end")]

  } else stop('Candidate gene not found in gff file.')

  # Get chromosome
  chrom <- as.character(gene_coords$seqnames)

  # Extract genome region for candidate gene
  region <- GenomicRanges::GRanges(paste0(chrom, ":", gene_coords$start, "-",
                                          gene_coords$end))

  if (!grepl("\\.vcf\\.gz$", vcf_file, ignore.case = TRUE)) {

    stop("VCF file input is not a .gz file type.")
  }

  # Create path to snpEff annotated VCF file
  vcf_path <- file.path(vcf_dir, vcf_file)

  # Creates .tbi index if not already present
  if (!file.exists(paste0(vcf_path, ".tbi"))) {

    Rsamtools::indexTabix(vcf_path, format = "vcf")

  }

  # Extract variants with tabix
  extract <- Rsamtools::scanTabix(vcf_path, param = region)

  # Combine extracted lines into a single file
  extract <- unlist(unlist(extract))

  # Combine header and extracted variants
  if (length(extract) == 0) stop("No variants found in the specified region!")


  # Combine the header and body
  # Extract header lines from the VCF
  header_lines <- Rsamtools::headerTabix(vcf_path)$header

  output_file <- file.path(output_path, sprintf('%s_%s.vcf', cand_gene_id,
                                                outfile_suffix))

  # Write the output VCF file
  writeLines(c(header_lines, extract), output_file)

}


################################################################################

#' Get the folder or file ID from a Google Drive shareable link.
#' @param drive_link A character value indicating the shareable Google Drive link.
#' @param is.folder A logical value indicating if link is for a folder or file.
#' Set to `FALSE` if link is for a shareable file.
#'
#'#' @examples
#' # example code
#' library(panGenomeBreedr)
#' folder_link <- "https://drive.google.com/drive/folders/1BotxaUb5emlrtgo473db3gDTUCLzKi70?usp=sharing"
#' folder_id <- get_google_id(drive_link = folder_link)
#'
#' @returns A character object of Google Drive folder or file ID.
#'
get_google_id <- function(drive_link, is.folder = TRUE) {

  # Regular expression to match the file or folder ID in the link
  # Extract the file ID using base R
  if (is.folder) {

    match <- regmatches(drive_link, regexpr("folders/([a-zA-Z0-9_-]+)", drive_link))

    # Remove the surrounding "folder/" to isolate the ID
    get_id <- sub("folders/", "", match)

  } else {

    match <- regmatches(drive_link, regexpr("file/d/([a-zA-Z0-9_-]+)", drive_link))

    # Remove the surrounding "file/d/" to isolate the ID
    get_id <- sub("file/d/", "", match)

  }

  return(get_id)

}


#' Download files in a shared Google Drive folder without restrictions.
#' @param drive_link A character value indicating the shareable Google Drive link.
#' @param output_path A character value indicating the path to a directory for
#' saving downloaded files.
#' @param is.folder A logical value indicating if link is for a folder or file.
#' Set to `FALSE` if link is for a shareable file.
#'
#' @examples
#' # example code
#' \donttest{
#' library(panGenomeBreedr)
#' f_link <- "https://drive.google.com/drive/folders/1BotxaUb5emlrtgo473db3gDTUCLzKi70?usp=sharing"
#' folder_path <- folder_download_gd(drive_link = f_link)
#' }
#'
#' @returns A list or vector containing the path to directory containing downloaded
#' files from Google Drive.
#'
#' @export
#'
folder_download_gd <- function(drive_link,
                               output_path = tempdir(),
                               is.folder = TRUE) {

  # Conditional logic for Bioconductor dependencies
  pkg <- "googledrive"

    if (!requireNamespace(pkg, quietly = TRUE)) {

      stop(sprintf("The '%s' package is required for this function.
                   Please install it via utils::install.packages('%s').", pkg, pkg),
           call. = FALSE)

    }

  # No authentication required from Google Drive
  googledrive::drive_deauth()
  googledrive::drive_user()

  if (is.folder) {

    # Get folder ID from shareable link
    folder_id <- get_google_id(drive_link = drive_link, is.folder = TRUE)
    p_file <- googledrive::drive_get(googledrive::as_id(folder_id))

    # Create a directory using the same folder ID in Google Drive
    dir <- file.path(output_path, unlist(p_file[,'name']))

    if(!dir.exists(dir)) dir.create(dir)

    # List all files in the google drive folder
    files <- googledrive::drive_ls(p_file)

    # File names to be downloaded
    file_ids <- as.list(files$id)

    # Paths to save files to be downloaded
    paths <- as.list(file.path(dir, files$name))

    # Download all files to dir
    mapply(googledrive::drive_download, file_ids, paths,
           MoreArgs = list(overwrite = TRUE))

    return(paths)

  } else {

    # Get file ID from shareable link
    file_id <- get_google_id(drive_link = drive_link, is.folder = FALSE)

    p_file <-  googledrive::drive_get(googledrive::as_id(file_id))

    # Create a directory using the same folder ID in Google Drive
    dir <- file.path(output_path, unlist(p_file[,'name']))

    googledrive::drive_download(p_file, overwrite = TRUE, path = dir)

    return(dir)

  }

}
