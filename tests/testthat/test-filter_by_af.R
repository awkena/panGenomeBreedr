test_that("filter_by_af() filters variants correctly by allele frequency thresholds", {

  # Mock genotype matrix
  gt_df <- data.frame(
    variant_id = c("v1", "v2", "v3", "v4"),
    chrom = c("Chr01", "Chr01", "Chr02", "Chr02"),
    pos = c(100, 200, 300, 400),
    s1 = c("0/0", "0/1", "1/1", "0/0"),
    s2 = c("0/1", "1/1", "0/0", "0/1"),
    s3 = c("1/1", "0/0", "0/1", "1/1"),
    stringsAsFactors = FALSE)

  # Full range: should retain all
  all_retained <- filter_by_af(gt_df, min_af = 0, max_af = 1)
  expect_equal(nrow(all_retained), 4)

  # Moderate AF filter: 0.3–0.7 should include v1 and v2
  mod_range <- filter_by_af(gt_df, min_af = 0.3, max_af = 0.7)
  expect_equal(nrow(mod_range), 4)
  expect_equal(sort(mod_range$variant_id), c("v1", "v2", "v3", "v4"))

  # Strict filter: AF = 0.5 only
  exact_half <- filter_by_af(gt_df, min_af = 0.5, max_af = 0.5)
  expect_equal(nrow(exact_half), 4)
  expect_equal(exact_half$variant_id, c("v1", "v2", "v3", "v4"))

  # Too strict: no variants pass
  no_pass <- filter_by_af(gt_df, min_af = 0.9, max_af = 1.0)
  expect_equal(nrow(no_pass), 0)

  # Expect error for empty input
  expect_error(filter_by_af(gt_df[0, ]), "No variant found in the genotype matrix.")
})
