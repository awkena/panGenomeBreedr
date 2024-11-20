test_that("get_alleles works", {

  set.seed(123)
  Call <- sample(c('A:A', 'A:T', 'T:T', '?', 'Uncallable'),
                size = 96,
                replace = TRUE,
                prob = c(.325, .25, .325, .05, .05))

  Call[c(88, 96)] <- 'NTC'

  allele_geno <- get_alleles(x = Call, data_type = 'kasp')

  expect_equal(allele_geno$alleles, c('A', 'T'))
  expect_equal(allele_geno$genotypes, c(homo1 = 'A:A',
                                        homo2 = 'T:T',
                                        het1 = 'A:T',
                                        het2 = 'T:A'))

})
