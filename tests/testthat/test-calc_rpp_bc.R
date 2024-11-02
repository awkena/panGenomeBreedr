test_that("calc_rpp_bc works", {

  # Create a numeric matrix of genotype scores for 10 markers and 5 samples
  num_dat <- matrix(c(rep(1, 10), rep(0, 10),
                      1, 1, 0.5, 1, 1, 1, 1, 1, 0, 1,
                      1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
                      1, 1, 0, 1, 1, 1, 1, 1, 1, 0.5 ),
                    byrow = TRUE, ncol = 10)

  rownames(num_dat) <- c('rp', 'dp', paste0('bc1_', 1:3))
  colnames(num_dat) <- paste0('S1', '_', c(floor(seq(1000, 10000, len = 8)),
                                           15000, 20000))

  # Get map file by parsing SNP IDs
  map_file <- parse_marker_ns(colnames(num_dat))

  # Calculate weighted RPP
  rpp <- calc_rpp_bc(x = num_dat,
                     map_file = map_file,
                     map_chr = 'chr',
                     map_pos = 'pos',
                     map_snp_ids = 'snpid',
                     rp = 1,
                     rp_num_code = 1,
                     weighted = TRUE)

  expect_equal(dim(rpp), c(5, 2))

})
