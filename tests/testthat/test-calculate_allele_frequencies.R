gt_region <- fetch_table_region(
  con = con_test,
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
)

test_that("calculate_allele_frequencies() works with full metadata", {
  af_metrics <- calculate_allele_frequencies(
    gt = gt_region,
    variant_id_col = "variant_id",
    chrom_col = "chrom",
    pos_col = "pos"
  )

  # Check output structure and dimensions
  expect_s3_class(af_metrics, "data.frame")
  expect_equal(nrow(af_metrics), nrow(gt_region))

  # Check for expected columns
  expected_cols <- c("variant_id", "chrom", "pos", "ref_af", "alt_af")
  expect_true(all(expected_cols %in% colnames(af_metrics)))

  # Check that frequencies sum to 1 (within a small tolerance)
  expect_true(all(abs(af_metrics$ref_af + af_metrics$alt_af - 1) < 1e-6))

  # Check that frequencies are valid probabilities
  expect_true(all(af_metrics$alt_af >= 0 & af_metrics$alt_af <= 1))
})

test_that("calculate_allele_frequencies() works with minimal metadata (ID only)", {
  af_metrics <- calculate_allele_frequencies(
    gt = gt_region,
    variant_id_col = "variant_id",
    chrom_col = NULL, 
    pos_col = NULL   
  )

  expect_s3_class(af_metrics, "data.frame")
  expect_equal(colnames(af_metrics), c("variant_id", "ref_af", "alt_af"))
})

test_that("calculate_allele_frequencies() handles single-row data frames", {
  single_variant_gt <- gt_region[1, , drop = FALSE]

  af_metrics <- calculate_allele_frequencies(gt = single_variant_gt)

  expect_s3_class(af_metrics, "data.frame")
  expect_equal(nrow(af_metrics), 1)
  expect_true(abs(af_metrics$ref_af + af_metrics$alt_af - 1) < 1e-6)
})

test_that("calculate_allele_frequencies() handles NA genotypes correctly", {
  # Create a small, manual data frame with a missing genotype call
  manual_gt <- data.frame(
    variant_id = "test_snp_1",
    sample1 = "0/0", # 0 alt alleles
    sample2 = "0/1", # 1 alt allele
    sample3 = "1/1", # 2 alt alleles
    sample4 = NA     # NA should be ignored
  )

  af_metrics <- calculate_allele_frequencies(gt = manual_gt, variant_id_col = "variant_id")

  # Expected alt_af = (0 + 1 + 2) / (2 * 3 non-NA samples) = 3 / 6 = 0.5
  expect_equal(af_metrics$alt_af, 0.5)
  expect_equal(af_metrics$ref_af, 0.5)
})
