test_that("introgression_check() calls each progeny at a region", {
  num_dat <- matrix(c(rep(1, 6), rep(0, 6),
                      1, 1, 0, 0, 1, 1,      # donor at the region
                      1, 1, 1, 1, 1, 1,      # recurrent
                      1, 1, 0.5, 0.5, 1, 1,  # heterozygous
                      1, 1, -5, -5, 1, 1),   # missing at the region
                    byrow = TRUE, ncol = 6)
  rownames(num_dat) <- c("rp", "dp", paste0("bc_", 1:4))
  colnames(num_dat) <- paste0("S1_", c(1000, 2000, 3000, 3500, 6000, 8000))
  map_file <- parse_marker_ns(colnames(num_dat))
  qtl <- list(qtl1 = c(chr = 1, start = 2500, end = 4000))

  res <- introgression_check(num_dat, map_file, parents = c("rp", "dp"),
                             trait_pos = qtl, expected = "donor")

  expect_equal(res$sample, paste0("bc_", 1:4))
  expect_equal(res$markers_from, rep("region", 4))
  expect_equal(res$n_markers, rep(2, 4))
  expect_equal(res$n_informative, c(2, 2, 2, 0))
  expect_equal(res$call, c("Donor", "Recurrent", "Heterozygous", "No informative markers"))
  expect_equal(res$pass, c(TRUE, FALSE, FALSE, NA))
  expect_equal(res$region[1], "2,500-4,000")

  het <- introgression_check(num_dat, map_file, parents = c("rp", "dp"),
                             trait_pos = qtl, expected = "het")
  expect_equal(het$pass, c(FALSE, FALSE, TRUE, NA))

  carrier <- introgression_check(num_dat, map_file, parents = c("rp", "dp"),
                                 trait_pos = qtl, expected = "carrier")
  expect_equal(carrier$pass, c(TRUE, FALSE, TRUE, NA))
})

test_that("introgression_check() uses the nearest markers for a position or an empty region", {
  num_dat <- matrix(c(rep(1, 4), rep(0, 4), 1, 0, 0, 1),
                    byrow = TRUE, ncol = 4)
  rownames(num_dat) <- c("rp", "dp", "bc_1")
  colnames(num_dat) <- paste0("S1_", c(1000, 2000, 3000, 9000))
  map_file <- parse_marker_ns(colnames(num_dat))

  pt <- introgression_check(num_dat, map_file, parents = c("rp", "dp"),
                            trait_pos = list(q = c(chr = 1, pos = 2500)),
                            n_nearest = 2)
  expect_equal(pt$markers_from, "nearest")
  expect_equal(pt$n_markers, 2)
  expect_equal(pt$call, "Donor")

  empty <- introgression_check(num_dat, map_file, parents = c("rp", "dp"),
                               trait_pos = list(q = c(chr = 1, start = 5000, end = 6000)),
                               n_nearest = 1)
  expect_equal(empty$markers_from, "nearest")
  expect_equal(empty$n_markers, 1)

  other_chr <- introgression_check(num_dat, map_file, parents = c("rp", "dp"),
                                   trait_pos = list(q = c(chr = 2, pos = 100)))
  expect_equal(other_chr$call, "No informative markers")
  expect_true(is.na(other_chr$pass))
})

test_that("introgression_check() checks its inputs", {
  num_dat <- matrix(c(1, 0, 1), ncol = 1, dimnames = list(c("rp", "dp", "bc"), "S1_100"))
  map_file <- parse_marker_ns(colnames(num_dat))
  expect_error(introgression_check(num_dat, map_file, parents = c("rp", "dp"),
                                   trait_pos = c(chr = 1, pos = 100)), "named list")
  expect_error(introgression_check(num_dat, map_file, snp_ids = "id", parents = c("rp", "dp"),
                                   trait_pos = list(q = c(chr = 1, pos = 100))), "columns")
  expect_error(introgression_check(num_dat, map_file, parents = c("rp", "dp"),
                                   trait_pos = list(q = c(chr = 1))), "start")
})
