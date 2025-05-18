test_that("cross_qc_heatmap() creates valid ggplot objects with trait annotations", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")

  # Simulated genotype matrix (5 lines x 10 SNPs)
  mat <- matrix(c(rep(1, 10),
                  rep(0, 10),
                  c(1,1,0.5,1,1,1,1,1,0,1),
                  c(1,1,0,1,1,1,1,1,1,1),
                  c(1,1,0,1,1,1,1,1,1,0.5)),
                byrow = TRUE, ncol = 10)

  rownames(mat) <- c("rp", "dp", "bc1_1", "bc1_2", "bc1_3")
  colnames(mat) <- paste0("S1_", c(floor(seq(1000, 10000, len = 8)), 15000, 20000))

  # Map file (parsed SNP info)
  map_file <- parse_marker_ns(colnames(mat))

  # Trait position annotations
  trait_loci <- list(
    loc1 = c(chr = 1, start = 2900, end = 4200),
    loc2 = c(chr = 1, start = 14200, end = 15800)
  )

  # Generate heatmap
  gg_list <- cross_qc_heatmap(
    x = mat,
    map_file = map_file,
    snp_ids = "snpid",
    chr = "chr",
    chr_pos = "pos",
    parents = c("rp", "dp"),
    trait_pos = trait_loci,
    text_scale_fct = 0.5,
    group_sz = 2L,
    pdf = FALSE,
    legend_title = "Genotypes",
    alpha = 0.8,
    text_size = 10
  )

  # Check output: should be a list of ggplot objects
  expect_true(is.list(gg_list))
  expect_true(all(sapply(gg_list, inherits, "gg")))
  expect_length(gg_list, ceiling(3 / 2))  # 3 progeny in group_sz = 2 → 2 batches

})
