test_that("find_indels works", {

  # Marker data
  dat <- data.frame(snp1 = c('C:C', 'A:A', 'C:A', 'C:A'),
                    snp2 = c('C:C', '-:-', 'C:C', 'C:C'),
                    snp3 = c('T:T', 'C:C', 'C:T', 'C:T'),
                    snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
                    snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
                    row.names = c('rp', 'dp', 'art_het1', 'art_het2'))

  # Check for unexpected homozygous genotypes
  geno_indel <- find_indels(x = dat,
                            rp_row = 1,
                            dp_row = 2)

  expect_equal(length(geno_indel), 2)
  expect_equal(dim(geno_indel$geno_indel), c(4, 2))
  expect_equal(dim(geno_indel$geno_non_indel), c(4, 3))



})
