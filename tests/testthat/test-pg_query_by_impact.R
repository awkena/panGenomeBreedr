test_that("pg_query_by_impact formats impact levels correctly and queries endpoint", {
  # Mock the internal fetcher
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Assert the API is routing to the correct endpoint
      expect_equal(endpoint, "/db/impact")

      # Assert the string collapsing worked perfectly
      expect_equal(query$impact_level, "HIGH,MODERATE")

      # Assert coordinate parameters are passed
      expect_equal(query$chrom, "Chr05")
      expect_equal(query$start, 100)

      # Return dummy data representing variants
      data.frame(
        variant_id = c("SNP_1", "INDEL_2"),
        impact = c("HIGH", "MODERATE"),
        stringsAsFactors = FALSE
      )
    }
  )

  # Execute the function with a vector of impact levels
  res <- pg_query_by_impact(
    impact_level = c("HIGH", "MODERATE"),
    chrom = "Chr05",
    start = 100,
    end = 500
  )

  # Validate output
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_equal(res$impact[1], "HIGH")
})


test_that("pg_query_by_impact defaults to all impact levels if none provided", {
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # The default argument is a vector of all 4 types, so it should collapse all 4
      expect_equal(query$impact_level, "HIGH,MODERATE,LOW,MODIFIER")
      data.frame() # Return empty df for simplicity
    }
  )

  res <- pg_query_by_impact(chrom = "Chr05", start = 100, end = 500)
  expect_s3_class(res, "data.frame")
})