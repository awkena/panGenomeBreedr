test_that("parent_poly works", {

  # Marker data
  dat <- data.frame(snp1 = c('C:C', 'A:A', 'C:A', 'C:A'),
                    snp2 = c('C:C', 'C:C', 'C:C', 'C:C'),
                    snp3 = c('T:T', 'C:C', 'C:T', 'C:T'),
                    snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
                    snp5 = c('T:T', 'T:T', 'T:A', 'T:A'),
                    row.names = c('rp', 'dp', 'art_het1', 'art_het2'))

  # Find polymorphic loci
  poly_loci <- parent_poly(x = dat,
                           rp_row = 1,
                           dp_row = 2,
                           sep = ':')

  expect_equal(dim(poly_loci), c(4, 3))
  expect_equal(names(poly_loci), c('snp1', 'snp3', 'snp4'))

})
