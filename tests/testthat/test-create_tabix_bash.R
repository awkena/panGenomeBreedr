test_that("create_tabix_bash works", {

  #skip_on_cran()
  path <- tempdir()
  setwd(path)

  # Create a tabix bash file to extract snpEff annotation for Sobic.003G421300
  path1 <- "/pl/active/Morris_CSU/Clara_Cruet"
  gff_path <- file.path(path1, 'resequence_pipeline_source_data/genes.gff')
  path2 <- "/pl/active/Morris_CSU/Sorghum_Genetic_Data/Sbicolor_v5.1"
  vcf_path = file.path(path2, 'VCF_update_data_nov/imputed_data_snpeff')

  expect_invisible(
  create_tabix_bash(slurm_par = c(nodes = 1,
                                  ntasks = 1,
                                  mem_per_cpu = '1G',
                                  time ='00:05:00'),
                    cand_gene_id = 'Sobic.003G421300',
                    tabix_output_path = "/scratch/alpine/awkena@colostate.edu/ssh_try/",
                    gff_path = gff_path,
                    vcf_path = vcf_path,
                    vcf_file_prefix = "Sorghum_d8.noduplicates.",
                    snp_vcf_suffix = ".snp._markernamesadded_imputed_snpeff",
                    indel_vcf_suffix = ".indel._markernamesadded_imputed_snpeff",
                    bash_out_path = path,
                    filename = 'tabix_test.sh')
  )

  # Path to read bash script
  tabix_script <- file.path(path, 'tabix_test.sh')

  # Read the bash file into R
  tabix_content <- readLines(tabix_script)
  expect_true(length(list.files(path = ".", pattern = "\\.sh$")) > 0)
  expect_equal(length(tabix_content), 106)
  on.exit(unlink(path))

})
