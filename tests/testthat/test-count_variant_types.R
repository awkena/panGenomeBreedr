test_that("count_variant_types() returns correct counts of variant types", {
  # Skip check if the global test connection wasn't initialized cleanly
  skip_if_not(
    exists("con_test"),
    message = "Global test connection 'con_test' not found."
  )

  # Run against your live mini Parquet dataset 
  result <- count_variant_types(con_test)

  # Validate output structure
  expect_true(is.data.frame(result))
  expect_named(result, c("variant_type", "n"))
  expect_true(nrow(result) >= 1)

  # Assert that reported type counts make logical sense
  expect_type(result$n, "double") 
  expect_true(all(result$n > 0))
  expect_true("SNP" %in% result$variant_type)

  # Test error handling for missing variant_type column 
  con_bad <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")

  # Write a malformed mock variants table lacking 'variant_type'
  bad_variants <- data.frame(
    variant_id = c("v1", "v2"),
    chrom = c("Chr05", "Chr05"),
    pos = c(75104537, 75104600)
  )
  DBI::dbWriteTable(con_bad, "variants", bad_variants)

  # Register cleanup to disconnect immediately if an unexpected error occurs
  on.exit(DBI::dbDisconnect(con_bad, shutdown = TRUE), add = TRUE)

  # Run expectation check against the bad connection object
expect_error(
  count_variant_types(con_bad),
  "The 'variants' table does not have a 'variant_type' column"
)

  # Explicitly clean up the temporary bad memory connection instance
  DBI::dbDisconnect(con_bad, shutdown = TRUE)
})