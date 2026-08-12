test_that("count_variant_types() works correctly in local mode", {
  counts <- count_variant_types(con = con_test)

  # Check output structure
  expect_s3_class(counts, "data.frame")
  expect_equal(colnames(counts), c("variant_type", "n"))
  expect_gt(nrow(counts), 0)

  #  Check content
  expect_true(all(c("SNP", "INDEL") %in% counts$variant_type))
  expect_true(is.numeric(counts$n))
  expect_true(all(counts$n > 0))
})

test_that("count_variant_types() handles errors correctly", {
  # Error on invalid connection
  expect_error(
    count_variant_types(con = NULL),
    "A valid 'con' object is required for 'local' mode. Please run connect_local_db() first.",
    fixed = TRUE
  )

  # Error on invalid table name (DBI will throw an error)
  expect_error(
    count_variant_types(con = con_test, variants_table = "non_existent_table")
  )

  # Error on invalid connect_db_mode
  expect_error(count_variant_types(con = con_test, connect_db_mode = "invalid"))
})
