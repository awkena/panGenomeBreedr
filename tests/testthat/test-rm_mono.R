test_that("rm_mono works", {
  # Create a dummy marker frame data with monomorphic markers
  mydf <- data.frame(S1_1001 = rep('A:A', 12),
                     S1_1011 = rep(c('T:T', 'T:C', 'C:C'), each = 4),
                     S2_1001 = rep(c('G:G', 'G:C', 'C:C'), each = 4),
                     S3_1101 = rep(c('G:G', 'C:C'), each = 6),
                     S1_1100 = rep(c('A:G', 'G:G', 'A:A'), each = 4),
                     S1_1201 = rep('C:A', 12),
                     row.names = sprintf('Ind_%02s', 1:12))

  # Remove monomorphic markers
  mydf_filtered <- rm_mono(mydata = mydf)

  expect_equal(dim(mydf_filtered), c(12, 4))
})
