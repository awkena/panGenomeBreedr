test_that("query_by_af() queries database and filters by allele frequency", {
  # Skip check if the global test connection wasn't initialized cleanly
  skip_if_not(
    exists("con_test"),
    message = "Global test connection 'con_test' not found."
  )

  # Test full range query within your target Chr05 window (should extract all regional markers)
  result <- query_by_af(
    con = con_test,
    min_af = 0,
    max_af = 1,
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )

  expect_true(is.data.frame(result))
  expect_true(nrow(result) >= 0)
  if (nrow(result) > 0) {
    expect_true(all(
      c("variant_id", "chrom", "pos", "ref_af", "alt_af") %in% colnames(result)
    ))
    expect_true(all(result$chrom == "Chr05"))
  }

  # Test with an impossible frequency threshold window (should return a clean empty dataframe)
  filtered_empty <- query_by_af(
    con = con_test,
    min_af = 1.1,
    max_af = 1.2,
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )

  expect_true(is.data.frame(filtered_empty))
  expect_equal(nrow(filtered_empty), 0)

  # Test with a nonexistent chromosome region to check graceful empty handling
  empty_region <- query_by_af(
    con = con_test,
    min_af = 0,
    max_af = 1,
    chrom = "ChrNonexistent",
    start = 1,
    end = 1000
  )

  expect_true(is.data.frame(empty_region))
  expect_equal(nrow(empty_region), 0)
})