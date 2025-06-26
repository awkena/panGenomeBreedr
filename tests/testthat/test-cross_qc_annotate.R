test_that("cross_qc_heatmap2() creates valid ggplot objects with trait annotations", {
  # skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")

  # Example genotype matrix
  num_dat <- matrix(c(rep(1, 10), rep(0, 10),
                      1, 1, 0.5, 1, 1, 1, 1, 1, 0, 1,
                      1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
                      1, 1, 0, 1, 1, 1, 1, 1, 1, 0.5),
                    byrow = TRUE, ncol = 10)
  rownames(num_dat) <- c("rp", "dp", paste0("bc1_", 1:3))
  colnames(num_dat) <- paste0("S1_", c(1000, 2000, 3000, 4000, 5000,
                                       6000, 7000, 8000, 15000, 20000))

  # Map file
  map_file <- parse_marker_ns(colnames(num_dat))

  # Locus positions to annotate
  loci <- list(
    "QTL1" = c("1", 2900, 4200),
    "QTL2" = c("1", 14200, 15800)
  )

  # Run function
  plt_list <- cross_qc_annotate(
    x = num_dat,
    map_file = map_file,
    snp_ids = "snpid",
    chr = "chr",
    chr_pos = "pos",
    parents = c("rp", "dp"),
    trait_pos = loci,
    pdf = FALSE
  )

  # Basic checks
  expect_type(plt_list, "list")
  expect_s3_class(plt_list[[1]], "gg")

  # Check that expected annotations are present
  gdata <- ggplot_build(plt_list[[1]])
  vline_layers <- which(vapply(gdata$data, function(d) "xintercept" %in% names(d), logical(1)))
  expect_gt(length(vline_layers), 0)

  xintercepts <- unlist(lapply(gdata$data[vline_layers], function(d) d$xintercept))
  expected_lines <- c(2900, 4200, 14200, 15800)
  expect_true(all(expected_lines %in% round(xintercepts)))
})
