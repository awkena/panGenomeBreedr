test_that("parent_het works", {

  # Marker data
  dat <- data.frame(snp1 = c('C:A', 'A:A', 'C:A', 'C:A'),
                    snp2 = c('C:C', 'G:G', 'C:C', 'C:C'),
                    snp3 = c('C:T', 'C:C', 'C:T', 'C:T'),
                    snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
                    snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
                    row.names = c('rp', 'dp', 'ind_1', 'ind_2'))

  # Find loci with at least one heterozygous parent genotype
  par_het <- parent_het(x = dat,
                        rp_row = 1,
                        dp_row = 2,
                        sep = ':')

  expect_equal(length(par_het), 2)
  expect_equal(dim(par_het$par_het), c(4, 2))
  expect_equal(dim(par_het$par_hom), c(4, 3))

})
