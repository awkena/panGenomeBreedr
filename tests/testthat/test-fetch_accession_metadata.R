test_that("fetch_accession_metadata() retrieves all metadata correctly", {
  all_meta <- fetch_accession_metadata(con = con_test)

  expect_s3_class(all_meta, "data.frame")
  expect_gt(nrow(all_meta), 0)
  expect_true(all(c("lib", "sample", "countryorigin", "lat", "lon") %in% colnames(all_meta)))
  # Check that the data is ordered by array_index
  expect_equal(all_meta$array_index, sort(all_meta$array_index))
})

test_that("fetch_accession_metadata() filters by column and value correctly", {
  ethiopia_meta <- fetch_accession_metadata(
    con = con_test,
    query_col = "countryorigin",
    query_value = "Ethiopia"
  )

  expect_s3_class(ethiopia_meta, "data.frame")
  expect_gt(nrow(ethiopia_meta), 0)
  expect_true(all(ethiopia_meta$countryorigin == "Ethiopia"))

  # Check that it returns fewer rows than the full set
  all_meta_rows <- nrow(fetch_accession_metadata(con = con_test))
  expect_lt(nrow(ethiopia_meta), all_meta_rows)
})

test_that("fetch_accession_metadata() returns an empty data.frame for non-matching filters", {
  no_match <- fetch_accession_metadata(
    con = con_test,
    query_col = "countryorigin",
    query_value = "Atlantis" # A value not in the dataset
  )
  expect_s3_class(no_match, "data.frame")
  expect_equal(nrow(no_match), 0)
})

test_that("fetch_accession_metadata() handles errors correctly", {
  # Invalid connection
  expect_error(fetch_accession_metadata(con = NULL), "A valid 'con' object is required", fixed = TRUE)

  # Invalid query column
  expect_error(fetch_accession_metadata(con = con_test, query_col = "invalid_column"), "Column 'invalid_column' not found", fixed = TRUE)

  # Invalid connect_db_mode
  expect_error(fetch_accession_metadata(con = con_test, connect_db_mode = "invalid"))
})
