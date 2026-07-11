test_that("pg_filter_by_af enforces input quality control checks", {
  # Check that passing a NULL object immediately stops execution
  expect_error(
    pg_filter_by_af(NULL),
    "The genotype matrix is empty or NULL."
  )

  # Check that passing an empty data frame triggers the same safeguard
  expect_error(
    pg_filter_by_af(data.frame()),
    "The genotype matrix is empty or NULL."
  )
})


test_that("pg_filter_by_af correctly subsets based on thresholds", {
  # Mock the internal calculator so we strictly control the frequencies
  local_mocked_bindings(
    pg_calc_af = function(gt, variant_id_col, chrom_col, pos_col) {
      data.frame(
        variant_id = c("v1", "v2", "v3", "v4"),
        alt_af = c(0.01, 0.05, 0.50, 0.99),
        stringsAsFactors = FALSE
      )
    }
  )

  # Dummy input 
  dummy_gt <- data.frame(variant_id = c("v1", "v2", "v3", "v4"))

  # Test common variant filtering (MAF >= 0.05 and <= 0.95)
  res <- pg_filter_by_af(dummy_gt, min_af = 0.05, max_af = 0.95)

  # Validate that the extreme frequencies were dropped and only v2 and v3 remain
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(c("v2", "v3") %in% res$variant_id))
  expect_false("v1" %in% res$variant_id) # 0.01 is out
  expect_false("v4" %in% res$variant_id) # 0.99 is out
})


test_that("pg_filter_by_af warns when thresholds eliminate all variants", {
  local_mocked_bindings(
    pg_calc_af = function(...) {
      data.frame(variant_id = c("v1", "v2"), alt_af = c(0.01, 0.99))
    }
  )

  dummy_gt <- data.frame(variant_id = c("v1", "v2"))

  # Apply incredibly strict thresholds that neither variant can pass
  expect_warning(
    res <- pg_filter_by_af(dummy_gt, min_af = 0.40, max_af = 0.60),
    "No variants passed the allele frequency filter."
  )

  # Verify it gracefully returns an empty data frame instead of breaking
  expect_equal(nrow(res), 0)
})