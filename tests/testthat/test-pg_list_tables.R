test_that("pg_list_tables successfully retrieves and returns table names", {
  # Mock the internal fetcher to intercept the network call
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Strict Assertion: Ensure the function targets the correct API route
      expect_equal(endpoint, "/db/tables")

      # Return a dummy character vector mimicking the database response
      c("variants", "annotations", "genotypes", "metadata")
    }
  )

  # Execute the function
  res <- pg_list_tables()

  # Validate the output structure matches expectations
  expect_type(res, "character")
  expect_length(res, 4)
  expect_equal(res[1], "variants")
  expect_true("genotypes" %in% res)
})