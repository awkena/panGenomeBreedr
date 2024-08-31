test_that("read_kasp_csv works", {

  path1 <-  system.file("extdata", "Genotyping_141.010_01.csv",
                                              package = "panGenomeBreedr",
                                              mustWork = TRUE)

                        file1 <- read_kasp_csv(file = path1, data_type = 'raw')
                        # Get KASP genotyping data for plotting
                        kasp_dat <- file1$Data
  expect_equal(dim(kasp_dat), c(768, 10) )
})

test_that("read_kasp_csv error works", {

  path1 <-  system.file("extdata", "Genotyping_141.010_01.csv",
                        package = "panGenomeBreedr",
                        mustWork = TRUE)

  expect_error(read_kasp_csv(file = path1,
                             row_tags = c('Statistics', 'DNA', 'SNPs','Scaling', 'datas'),
                             data_type = 'raw'),
               "datas is not contained in the input file.")
})
