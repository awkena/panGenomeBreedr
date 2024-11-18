test_that("find_unexp_homs works", {

  # Marker data
  dat <- data.frame(snp1 = c('C:C', 'A:A', 'C:A', 'C:A'),
                    snp2 = c('C:C', '-:-', 'C:C', 'C:C'),
                    snp3 = c('T:T', 'C:C', 'C:T', 'C:T'),
                    snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
                    snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
                    row.names = c('rp', 'dp', 'art_het1', 'art_het2'))

  # Check for unexpected homozygous genotypes
  homs_unexp <- find_unexp_homs(x = dat,
                                rp_row = 1,
                                dp_row = 2)

  expect_equal(length(homs_unexp), 2)
  expect_equal(dim(homs_unexp$geno_unexp), c(4, 2))
  expect_equal(dim(homs_unexp$geno_exp), c(4, 3))
})
