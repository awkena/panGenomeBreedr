test_that("pg_summarize_tables queries correct endpoint and returns summary schema", {
  # Mock the internal fetcher
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Assert the function targets the correct route
      expect_equal(endpoint, "/db/summary")

      # Return a dummy response mimicking a database table summary
      data.frame(
        table = c("variants", "annotations", "genotypes", "metadata"),
        n_rows = c(50000, 150000, 50000, 1676),
        stringsAsFactors = FALSE
      )
    }
  )

  # Execute the function
  res <- pg_summarize_tables()

  # Validate the output structure
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 4)
  expect_true(all(c("table", "n_rows") %in% colnames(res)))
  expect_equal(res$table[1], "variants")
  expect_equal(res$n_rows[4], 1676)
})