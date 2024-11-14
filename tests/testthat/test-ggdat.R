test_that("ggdat works", {

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

  # Convert num_geno to a long format data frame
  df <- gg_dat(num_mat = num_dat,
               map_file = map_file)

  expect_equal(dim(df), c(50, 5))
})
