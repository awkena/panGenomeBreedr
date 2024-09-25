test_that("parse_marker_ns works", {

  snps <- paste0('S', 1:10, '_', 101:110)
  map_file <- parse_marker_ns(x = snps, sep = '_', prefix = 'S')

  expect_equal(dim(map_file), c(10, 3))

})
