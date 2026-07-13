test_that("pg_count_variant_types queries correct endpoint and passes custom parameters", {
  # Mock the internal fetcher to intercept the network call
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Assert the function targets the correct API route
      expect_equal(endpoint, "/stats/variant_types")

      # Assert the custom table name parameter was passed correctly
      expect_equal(query$variants_table, "custom_variants_table")

      # Return a dummy response mimicking a variant type distribution table
      data.frame(
        variant_type = c("SNP", "INDEL", "MNP"),
        n = c(5000, 150, 20),
        stringsAsFactors = FALSE
      )
    }
  )

  # Execute the function with a custom table name
  res <- pg_count_variant_types(variants_table = "custom_variants_table")

  # Validate the output structure matches the promised schema
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
  expect_true(all(c("variant_type", "n") %in% colnames(res)))
  expect_equal(res$variant_type[1], "SNP")
  expect_equal(res$n[2], 150)
})


test_that("pg_count_variant_types handles default arguments correctly", {
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Assert that the default table name is "variants" when the user provides nothing
      expect_equal(query$variants_table, "variants")

      data.frame() # Return empty data frame for simplicity
    }
  )

  # Execute the function without specifying any arguments
  res <- pg_count_variant_types()
  expect_s3_class(res, "data.frame")
})