test_that("pg_list_table_columns handles arguments and queries correct endpoint", {
  # Mock the internal fetcher
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Route validation
      expect_equal(endpoint, "/db/columns")

      # Parameter validation
      expect_true(!is.null(query$table_name))

      # Return a dummy response mimicking a schema table
      data.frame(
        column_name = c("variant_id", "chrom", "pos"),
        data_type = c("character", "character", "integer"),
        stringsAsFactors = FALSE
      )
    }
  )

  # Test successful execution with a valid argument
  res <- pg_list_table_columns(table_name = "variants")

  # Validate output structure
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
  expect_true("data_type" %in% colnames(res))
  expect_equal(res$column_name[1], "variant_id")

  # Test default argument fallback 
  res_default <- pg_list_table_columns()
  expect_s3_class(res_default, "data.frame")

  # Test error handling for invalid table names
  expect_error(
    pg_list_table_columns(table_name = "invalid_table"),
    "'arg' should be one of"
  )
})