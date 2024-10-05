
test_that("kasp_marker_design works", {
 path1 <- "https://raw.githubusercontent.com/awkena/panGB/main/Chr02.fa.gz"

  # path1 <-  system.file("extdata", "Chr02.fa.gz",
  #                      package = "panGenomeBreedr",
  #                      mustWork = TRUE)
  path2 <-  system.file("extdata", "Sobic.002G302700_SNP_snpeff.vcf",
                        package = "panGenomeBreedr",
                        mustWork = TRUE)
  path <- tempdir()
  setwd(path)

  # Example to test new function on a substitution variant
  ma1 <- kasp_marker_design(vcf_file = path2,
                            genome_file = path1,
                            marker_ID = "SNP_Chr02_69200443",
                            chr = "Chr02",
                            plot_draw = TRUE,
                            region_name = "ma1",
                            plot_file = path,
                            maf = 0.05)


  expect_true(length(list.files(path = ".", pattern = "\\.pdf$")) > 0)
  expect_equal(dim(ma1), c(1, 8))
  on.exit(unlink(path))


})
