test_that("pg_calc_af correctly computes frequencies for phased and unphased genotypes", {
  # Build a matrix containing a mix of standard VCF notation formats
  dummy_gt <- data.frame(
    variant_id = c("v1", "v2", "v3"),
    chrom = c("Chr1", "Chr1", "Chr1"),
    pos = c(10, 20, 30),
    sample_A = c("0/0", "0/1", "1/1"), # Unphased
    sample_B = c("0|0", "1|0", "1|1"), # Phased
    stringsAsFactors = FALSE
  )

  res <- pg_calc_af(
    gt = dummy_gt,
    variant_id_col = "variant_id",
    chrom_col = "chrom",
    pos_col = "pos"
  )

  # Validate output structure
  expect_s3_class(res, "data.frame")
  expect_true(all(c("ref_af", "alt_af") %in% colnames(res)))

  # Validate exact mathematical outcomes
  # v1: Hom Ref (0 + 0 dosage = 0 alt alleles out of 4)
  expect_equal(res$alt_af[1], 0.0)
  expect_equal(res$ref_af[1], 1.0)

  # v2: Heterozygous (1 + 1 dosage = 2 alt alleles out of 4)
  expect_equal(res$alt_af[2], 0.5)
  expect_equal(res$ref_af[2], 0.5)

  # v3: Hom Alt (2 + 2 dosage = 4 alt alleles out of 4)
  expect_equal(res$alt_af[3], 1.0)
  expect_equal(res$ref_af[3], 0.0)
})


test_that("pg_calc_af handles missing genotypes gracefully", {
  # Simulate a matrix where some samples dropped out or failed QC
  dummy_gt_na <- data.frame(
    variant_id = c("v1", "v2"),
    sample_A = c("0/0", NA),
    sample_B = c("./.", "1/1"), # Unrecognized strings should become NA during matrix conversion
    sample_C = c("1/1", "1/1"),
    stringsAsFactors = FALSE
  )

  res <- pg_calc_af(dummy_gt_na)

  # v1: sample_B is missing. A is 0, C is 2. Mean dosage = 1. alt_af = 1 / 2 = 0.5
  expect_equal(res$alt_af[1], 0.5)

  # v2: sample_A is missing. B is 2, C is 2. Mean dosage = 2. alt_af = 2 / 2 = 1.0
  expect_equal(res$alt_af[2], 1.0)
})


test_that("pg_calc_af enforces required metadata and sample columns", {
  dummy_gt <- data.frame(
    variant_id = "v1",
    chrom = "Chr1",
    sample_A = "0/0"
  )

  # Trap missing ID column
  expect_error(
    pg_calc_af(dummy_gt, variant_id_col = "wrong_id_name"),
    "not found in the input data frame"
  )

  # Trap missing sample columns 
  expect_error(
    pg_calc_af(
      dummy_gt[, c("variant_id", "chrom")],
      variant_id_col = "variant_id",
      chrom_col = "chrom"
    ),
    "No sample columns detected"
  )
})