test_that("list_columns() works for 'variants' table in local mode", {
  schema <- list_columns(con = con_test, table_name = "variants", connect_db_mode = 'local')

  # Check that the output is a data frame
  expect_s3_class(schema, "data.frame")
  expect_true(all(c("column_name", "column_type", "null", "key") %in% colnames(schema)))
  expect_true(all(c("variant_id", "chrom", "pos", "ref", "alt") %in% schema$column_name))
})

test_that("list_columns() works for all expected tables in local mode", {
  expected_tables <- c("variants", "annotations", "genotypes", "metadata")
  for (tbl in expected_tables) {
    schema <- list_columns(con = con_test, table_name = tbl, connect_db_mode = 'local')
    expect_s3_class(schema, "data.frame")
    expect_gt(nrow(schema), 0, label = paste("Row count for table:", tbl))
    expect_true(all(c("column_name", "column_type") %in% colnames(schema)), label = paste("Columns for table:", tbl))
  }
})

test_that("list_columns() errors correctly without a valid connection", {
  expect_error(
    list_columns(con = NULL, table_name = "variants", connect_db_mode = 'local'),
    "A valid 'con' object is required for 'local' mode. Please run connect_local_db() first.",
    fixed = TRUE
  )
})

test_that("list_columns() validates table_name argument", {
  expect_error(
    list_columns(con = con_test, table_name = "invalid_table", connect_db_mode = 'local'),
    "'arg' should be one of \"variants\", \"annotations\", \"genotypes\", \"metadata\""
  )
})

test_that("list_columns() validates connect_db_mode argument", {
  expect_error(
    list_columns(con = con_test, table_name = "variants", connect_db_mode = "invalid_mode"),
    "'arg' should be one of \"local\", \"online\""
  )
})
