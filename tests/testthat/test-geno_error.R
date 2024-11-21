test_that("geno_error works", {

  # Marker data with errors in snp1 and snp4
  dat <- data.frame(snp1 = c('C:A', 'A:A', 'C:A', 'C:T'),
                    snp2 = c('C:C', 'G:G', 'C:C', 'C:C'),
                    snp3 = c('C:T', 'C:C', 'C:T', 'C:T'),
                    snp4 = c('G:G', '-:-', 'G:-', 'G:A'),
                    snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
                    row.names = c('rp', 'dp', 'ind_1', 'ind_2'))

  # Check for genotype call  error for each SNP
  geno_mat <- geno_error(x = dat,
                         rp_row = 1,
                         dp_row = 2,
                         sep = ':',
                         data_type = 'kasp')

  expect_equal(dim(geno_mat$geno_good), c(4, 3))
  expect_equal(dim(geno_mat$geno_err), c(4, 2))

})
