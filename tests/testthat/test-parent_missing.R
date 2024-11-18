test_that("parent_missing works", {

  # Marker data
  dat <- data.frame(snp1 = c('C:C', 'A:A', 'C:A', 'C:A'),
                    snp2 = c('C:C', NA, 'C:C', 'C:C'),
                    snp3 = c(NA, 'C:C', 'C:T', 'C:T'),
                    snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
                    snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
                    row.names = c('rp', 'dp', 'ind_1', 'ind_2'))

  # Find loci with at least one missing parent genotype
  par_miss <- parent_missing(x = dat,
                             rp_row = 1,
                             dp_row = 2)

  expect_equal(length(par_miss), 2)
  expect_equal(dim(par_miss$par_missing), c(4, 2))
  expect_equal(dim(par_miss$par_present), c(4, 3))

})
