gt_region <- fetch_table_region(
  con = con_test,
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
)

test_that("filter_by_allele_frequency() works correctly with a given range", {
  min_freq <- 0.05
  max_freq <- 0.95

  filtered_vars <- filter_by_allele_frequency(
    gt = gt_region,
    min_af = min_freq,
    max_af = max_freq
  )

  #  Check output structure and that rows were filtered
  expect_s3_class(filtered_vars, "data.frame")
  expect_true(nrow(filtered_vars) < nrow(gt_region))

  #  Check for expected columns
  expected_cols <- c("variant_id", "chrom", "pos", "ref_af", "alt_af")
  expect_true(all(expected_cols %in% colnames(filtered_vars)))

  #  Check that filtering logic is correct
  expect_true(all(filtered_vars$alt_af >= min_freq & filtered_vars$alt_af <= max_freq))
})

test_that("filter_by_allele_frequency() works with default min/max values", {
  # With default min_af=0 and max_af=1, it should return all variants
  filtered_vars <- filter_by_allele_frequency(gt = gt_region)
  all_freqs <- calculate_allele_frequencies(gt = gt_region)

  expect_equal(nrow(filtered_vars), nrow(all_freqs))
})

test_that("filter_by_allele_frequency() errors on empty input", {
  empty_gt <- gt_region[0, , drop = FALSE]

  expect_error(
    filter_by_allele_frequency(gt = empty_gt),
    "No variant found in the genotype matrix.",
    fixed = TRUE
  )
})
