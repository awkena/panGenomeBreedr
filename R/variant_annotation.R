#' Create a tabix bash file to run on HPC cluster
#' @param slurm_par A named character vector of `length = 6`, indicating SLURM job
#' parameters for submitting jobs to a computing cluster managed by a SLURM
#' workload manager. It must be specified in this order: number of compute nodes,
#' number of tasks, the amount of memory per CPU core, the maximum runtime
#' for the job, and the partition hardware.
#' @param  cand_gene_id A character value specifying the candidate gene ID.
#' @param tabix_output_path A character value indicating the path to directory for
#' saving tabix results output.
#' @param gff_path A character value indicating the path to the GFF file.
#' @param snpeff_snp,snpeff_indel A logical value indicating whether to run tabix
#' command on SNP or InDel snpEff annotated files.
#' @param vcf_path A character value indicating the path to directory containing
#' snpEff annotated SNP and InDel VCF files.
#' @param vcf_file_prefix A character avlue indicating the prefix for snpEff
#' annotated VCF files.
#' @param snp_vcf_suffix,indel_vcf_suffix A character value indicating the suffix
#' for snpEff annotated VCF files for SNPs or InDels.
#' @param bash_out_path  A character value indicating the path to directory for
#' saving tabix bash file.
#' @param filename A character value indicating the file name of the tabix bash
#' file.
#'
#' @returns A tabix bash file for extracting snpEff annotated variants in SNP or
#' InDel VCF files for any candidate gene.
#'
#' @examples
#' # example code
#' library(panGenomeBreedr)
#'
#' # Create a tabix bash file to extract snpEff annotation for Sobic.003G421300
#' # Path to GFF file
#' path1 <- "/pl/active/Morris_CSU/Clara_Cruet"
#' gff_path <- file.path(path1, 'resequence_pipeline_source_data/genes.gff')
#'
#' # Path to snpEff annotated vcf file
#' path2 <- "/pl/active/Morris_CSU/Sorghum_Genetic_Data/Sbicolor_v5.1"
#' vcf_path <- file.path(path2, 'VCF_update_data_nov/imputed_data_snpeff')
#'
#' create_tabix_bash(cand_gene_id = 'Sobic.003G421300',
#'                   tabix_output_path = "/scratch/alpine/awkena@colostate.edu/ssh_try/",
#'                   gff_path = gff_path,
#'                   vcf_path = vcf_path,
#'                   vcf_file_prefix = "Sorghum_d8.noduplicates.",
#'                   snp_vcf_suffix = ".snp._markernamesadded_imputed_snpeff",
#'                   indel_vcf_suffix = ".indel._markernamesadded_imputed_snpeff",
#'                   bash_out_path = tempdir(),
#'                   filename = 'tabix_test.sh')
#'
#' # Path to read bash script
#' tabix_script <- file.path(tempdir(), 'tabix_test.sh')
#'
#' # Read the bash file into R
#' tabix_content <- readLines(tabix_script)
#' print(tabix_content)
#'
#' @export
create_tabix_bash <- function(slurm_par = c(nodes = 1,
                                            ntasks = 1,
                                            mem = '1G',
                                            time ='00:05:00',
                                            partition = 'amilan',
                                            error = 'tabix_test.err'),
                              cand_gene_id = 'Sobic.003G421300',
                              tabix_output_path,
                              gff_path,
                              snpeff_snp = TRUE,
                              snpeff_indel = TRUE,
                              vcf_path,
                              vcf_file_prefix,
                              snp_vcf_suffix,
                              indel_vcf_suffix,
                              bash_out_path = tempdir(),
                              filename = 'tabix.sh'
) {

  # 1. shebang
  shebang <- "#!/bin/bash"

  # 2. Define slurm parameters
  slurm <- c(sprintf("#SBATCH --nodes=%s", slurm_par[1]),
             sprintf("#SBATCH --ntasks=%s", slurm_par[2]),
             sprintf("#SBATCH --mem=%s", slurm_par[3]),
             sprintf("#SBATCH --time=%s", slurm_par[4]),
             sprintf("#SBATCH --partition=%s", slurm_par[5]),
             sprintf("#SBATCH  --error=%s", slurm_par[6]),
             "",
             "")

  # 3. Add ReadMe text
  readme <- c("############################      README        ##########################################",
              "# This script will extract variants (SNPs and InDels) within your candidate gene.",
              "# It will extract gene region from snpEff annotations on VCF files.",
              "# It requires the candidate gene ID, gff annotation for the species, snpEff annotated",
              "# VCF files for SNPs and InDels.",
              "# Users must also have the dependency packages htslib, bcftools, and samtools installed on",
              "# their environment.",
              "# It assumes that the user has already run snpEff on variant calls.",
              "#########################################################################################",
              "",
              "",
              "set -a",
              "",
              "")

  # 4. Candidate gene input
  candidate_gene <- c("#1# .... Gene information .......#1#",
                      "# Gene name: # This must be the gene ID in gff file",
                      sprintf('Gene="%s"', cand_gene_id),
                      "",
                      "")
  # 5. tabix results output path
  tabix_res_out <- c("#2#... Output directory: it must be inside quotes ....#2#",
                     sprintf('output="%s" ', tabix_output_path),
                     "",
                     "",
                     "set +a",
                     "",
                     "")

  # 6. Load required packages
  pkgs <- c("# Loading required modules",
            "",
            "# htslib",
            "module load htslib",
            "",
            "# bcftools",
            "module load bcftools",
            "",
            "# samtools",
            "module load samtools",
            "",
            "")

  # Set environment variables
  env_var <- c("",
               "# Set environment variables",
               "export PATH=/usr/local/bin:$PATH",
               "export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH",
               "",
               "# Optional: Log loaded modules and environment variables for debugging",
               "module list",
               "env | grep -E 'PATH|LD_LIBRARY_PATH'",
               "")

  # 7. Path to gff file
  gff_dir <- c("#----------------------------------------------------------------------------------------#",
               "#~~~~~~~~Start of tabix code~~~~~~~~#",
               "#----------------------------------------------------------------------------------------#",
               "",
               "###### Getting gene information #####",
               "# Set path to gff annotation for the species ",
               sprintf('gff="%s"', gff_path),
               "echo $gff",
               "",
               "##### Extracting regions ######",
               "",
               "")

  # 8. Extracting regions from gff file
  extract_region <- c("# grepl file for gene",
                      "geneArray=($(grep 'gene.*'$Gene'' $gff))",
                      "echo ${geneArray[@]}",
                      "",
                      "#chromosome start ",
                      "Start=${geneArray[3]} ",
                      "echo $Start",
                      "",
                      "#Chromosome end",
                      "End=${geneArray[4]} ",
                      "echo $End",
                      "",
                      "# Chromosome number",
                      "Chromosome=${geneArray[0]} ",
                      "echo $Chromosome",
                      "",
                      "#region",
                      'region="$Chromosome":$Start-$End',
                      "echo $region",
                      "",
                      "")

  # 9. Path to snpEff annotated VCF files
  vcf_path <- file.path(vcf_path, vcf_file_prefix)

  # SNPS
  if(snpeff_snp) {

    snpeff_snp_vcf <- c("# Specify path to snpEff annotated  SNP VCF files",
                        "# SNPs",
                        "",
                        sprintf('genotype_snp="%s"$Chromosome"%s.vcf.gz"', vcf_path,
                                snp_vcf_suffix),
                        "echo $genotype_snp",
                        "",
                        "",
                        "#SNPS",
                        'file_snp="$Gene"_SNP_snpeff.vcf',
                        "echo $file_snp",
                        "",
                        "",
                        'tabix -h "$genotype_snp" "$region" > "$output"/"$file_snp"',
                        "")

  }  else snpeff_snp_vcf <- NULL

  # InDels
  if(snpeff_indel) {

    snpeff_indel_vcf <- c("# Specify path to snpEff annotated Indel VCF files",
                          "# INDELS",
                          "",
                          sprintf('genotype_indel="%s"$Chromosome"%s.vcf.gz"', vcf_path,
                                  indel_vcf_suffix),
                          "echo $genotype_indel",
                          "",
                          "",
                          "#Indels",
                          'file_indel="$Gene"_INDEL_snpeff.vcf',
                          "echo $file_indel",
                          "",
                          "",
                          'tabix -h "$genotype_indel" "$region" > "$output"/"$file_indel"',
                          "")

  } else snpeff_indel_vcf <- NULL

  # Write the script to the output file
  bash_script <- c(shebang, slurm, readme, candidate_gene, tabix_res_out,
                   pkgs, env_var, gff_dir, extract_region, snpeff_snp_vcf,
                   snpeff_indel_vcf)

  # Save tabix bash file in bash_out_path
  writeLines(bash_script, file.path(bash_out_path, filename))

}


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
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges
#' @import Rsamtools
#' @export
#'
extract_variant <- function(cand_gene_id,
                            gff_path,
                            vcf_dir,
                            vcf_file,
                            output_path = tempdir(),
                            outfile_suffix = 'variants') {

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
#' @importFrom googledrive drive_get as_id drive_ls drive_download drive_deauth drive_user drive_get
#' @export
#'
folder_download_gd <- function(drive_link,
                               output_path = tempdir(),
                               is.folder = TRUE) {

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
