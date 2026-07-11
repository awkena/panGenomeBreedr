test_that("pg_get_sample_metadata retrieves all records when no arguments are provided", {
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Assert the API is routing to the correct endpoint
      expect_equal(endpoint, "/db/metadata")

      # Assert that empty arguments translate to NULL in the query list
      expect_null(query$query_col)
      expect_null(query$query_value)

      # Return a dummy full metadata table
      data.frame(
        lib = c("S1", "S2", "S3"),
        countryorigin = c("Ghana", "Ethiopia", "Togo"),
        stringsAsFactors = FALSE
      )
    }
  )

  res <- pg_get_sample_metadata()

  # Validate output
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
})

test_that("pg_get_sample_metadata passes filtering parameters correctly", {
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Assert the filtering parameters were packaged correctly
      expect_equal(query$query_col, "countryorigin")
      expect_equal(query$query_value, "Ghana")

      # Return a dummy filtered metadata table
      data.frame(
        lib = c("S1", "S4"),
        countryorigin = c("Ghana", "Ghana"),
        stringsAsFactors = FALSE
      )
    }
  )

  # Execute the function with explicit filters
  res <- pg_get_sample_metadata(
    query_col = "countryorigin",
    query_value = "Ghana"
  )

  # Validate output reflects the mocked subset
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(res$countryorigin == "Ghana"))
})