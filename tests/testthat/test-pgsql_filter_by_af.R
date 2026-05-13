test_that("pgsql_filter_by_af correctly subsets variants based on frequency", {
  # 1. Create a "Test Population"
  test_data <- data.frame(
    variant_id = c("v_fixed_ref", "v_seg", "v_fixed_alt"),
    chrom = rep("Chr03", 3),
    pos = c(1, 2, 3),
    S1 = c("0/0", "0/1", "1/1"),
    S2 = c("0/0", "0/0", "1/1"),
    stringsAsFactors = FALSE
  )

  # Test: Common Variants Only (0.05 to 0.95)
  # Only v_seg should pass (AF = 0.25)
  res_common <- pgsql_filter_by_af(test_data, min_af = 0.1, max_af = 0.9)
  expect_equal(nrow(res_common), 1)
  expect_equal(res_common$variant_id, "v_seg")

  #  Test: Inclusive Range (0 to 1)
  # All should pass
  res_all <- pgsql_filter_by_af(test_data, min_af = 0, max_af = 1)
  expect_equal(nrow(res_all), 3)

  # Test: Fixed Alternate Only
  res_fixed <- pgsql_filter_by_af(test_data, min_af = 0.9, max_af = 1)
  expect_equal(nrow(res_fixed), 1)
  expect_equal(res_fixed$variant_id, "v_fixed_alt")
})

test_that("pgsql_filter_by_af handles empty inputs and no-passers", {
  #  Test NULL/Empty input error
  expect_error(pgsql_filter_by_af(NULL), "genotype matrix is empty")

  #  Test Warning when nothing passes
  test_data <- data.frame(
    variant_id = "v1",
    chrom = "Chr01",
    pos = 10,
    S1 = "0/0",
    S2 = "0/0",
    stringsAsFactors = FALSE
  )

  # AF is 0, filter is 0.5-1.0
  expect_warning(
    pgsql_filter_by_af(test_data, min_af = 0.5, max_af = 1),
    "No variants passed"
  )
})