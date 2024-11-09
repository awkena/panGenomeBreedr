test_that("order_markers works", {

  map_file <- data.frame(snpid = paste0('S', rep(1:2, 5), '_', 1001:1005),
                         chr = rep(1:2, 5),
                         pos = rep(1001:1005, 2))

  # Order map file
  map_file <- order_markers(x = map_file,
                            chr_col = 'chr',
                            pos_col = 'pos')

  expect_equal(dim(map_file), c(10, 3))
  expect_equal(map_file$chr, rep(1:2, each = 5))
  expect_equal(map_file$pos, rep(1001:1005, 2))

})
