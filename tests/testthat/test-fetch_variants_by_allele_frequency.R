test_that("fetch_variants_by_allele_frequency() works correctly", {
  common_vars <- fetch_variants_by_allele_frequency(
    con = con_test,
    min_af = 0.05,
    max_af = 0.95,
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )

  expect_s3_class(common_vars, "data.frame")
  expect_gt(nrow(common_vars), 0)
  expect_true(all(common_vars$alt_af >= 0.05 & common_vars$alt_af <= 0.95))
  expect_true(all(c("variant_id", "chrom", "pos", "ref_af", "alt_af") %in% colnames(common_vars)))
})

test_that("fetch_variants_by_allele_frequency() validates inputs", {
  # Chromosome is required
  expect_error(
    fetch_variants_by_allele_frequency(con = con_test),
    "Chromosome must be specified.",
    fixed = TRUE
  )

  # Connection is required for local mode (error comes from internal fetch_table_region call)
  expect_error(
    fetch_variants_by_allele_frequency(con = NULL, chrom = "Chr05"),
    "A valid 'con' object is required for 'local' mode. Please run connect_local_db() first.",
    fixed = TRUE
  )
})

test_that("fetch_variants_by_allele_frequency() handles no results gracefully", {
  # Test a region with no variants, expecting a message
  expect_message(
    no_vars <- fetch_variants_by_allele_frequency(con = con_test, chrom = "Chr01", start = 1, end = 10),
    "No genotype data found for the specified region."
  )
  expect_s3_class(no_vars, "data.frame")
  expect_equal(nrow(no_vars), 0)

  # Test a filter that returns no variants, expecting a warning
  expect_warning(
    rare_vars <- fetch_variants_by_allele_frequency(con = con_test, min_af = 0.99, chrom = "Chr05"),
    "No variants passed the allele frequency filter. Check your min_af/max_af thresholds."
  )
  expect_s3_class(rare_vars, "data.frame")
  expect_equal(nrow(rare_vars), 0)
})
