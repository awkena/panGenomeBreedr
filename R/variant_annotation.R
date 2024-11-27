#' Create a tabix bash file to run on HPC cluster
#' @param slurm_par A named character vector of `length = 4`, indicating SLURM job
#' parameters for submitting jobs to a computing cluster managed by a SLURM
#' workload manager. It must be specified in this order: number of compute nodes,
#' number of tasks, the amount of memory per CPU core, and the maximum runtime
#' for the job.
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
#' create_tabix_bash(slurm_par = c(nodes = 1,
#'                                 ntasks = 1,
#'                                 mem_per_cpu = '1G',
#'                                 time ='00:05:00'),
#'                   cand_gene_id = 'Sobic.003G421300',
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
                                            mem_per_cpu = '1G',
                                            time ='00:05:00'),
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
             sprintf("#SBATCH --mem-per-cpu=%s", slurm_par[3]),
             sprintf("#SBATCH --time=%s", slurm_par[4]),
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
            "# module load htslib",
            "",
            "# bcftools",
            "module load bcftools",
            "",
            "# samtools",
            "module load samtools",
            "",
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
                   pkgs, gff_dir, extract_region, snpeff_snp_vcf, snpeff_indel_vcf)

  writeLines(bash_script, file.path(bash_out_path, filename))

}
