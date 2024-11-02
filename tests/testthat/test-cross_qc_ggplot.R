test_that("cross_qc_ggplot works", {

  #skip_on_cran()
  path <- tempdir()
  setwd(path)

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

  # Create a heatmap that compares the parents to progenies
  expect_invisible(cross_qc_ggplot(x = num_dat,
                                   map_file = map_file,
                                   snp_ids = 'snpid',
                                   chr = 'chr',
                                   chr_pos = 'pos',
                                   value = 'value',
                                   parents = c('rp', 'dp'),
                                   group_sz = 3L,
                                   pdf = TRUE,
                                   legend_title = 'Heatmap_key',
                                   alpha = 0.8,
                                   text_size = 14,
                                   width = 12,
                                   height = 10))

  expect_true(length(list.files(path = ".", pattern = "\\.pdf$")) > 0)
  on.exit(unlink(path))
})
