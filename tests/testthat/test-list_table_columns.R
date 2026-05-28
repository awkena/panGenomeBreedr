test_that("list_table_columns() returns correct column info for allowed tables", {
  skip_if_not(
    exists("con_test"),
    message = "Global test connection 'con_test' not found."
  )

  # Test variants table view column presence
  variant_info <- list_table_columns(con_test, "variants")
  expect_true(is.data.frame(variant_info))
  expect_true("variant_id" %in% variant_info$column_name)
  expect_true("chrom" %in% variant_info$column_name)

  # Test annotations table view column presence
  annot_info <- list_table_columns(con_test, "annotations")
  expect_true(is.data.frame(annot_info))
  expect_true("variant_id" %in% annot_info$column_name)

  # Test genotypes table view column presence
  geno_info <- list_table_columns(con_test, "genotypes")
  expect_true(is.data.frame(geno_info))
  expect_true("variant_id" %in% geno_info$column_name)

  #  Check constraint check mapping error handles correctly
  expect_error(
    list_table_columns(con_test, table_name = "nonexistent"),
    "'arg' should be one of"
  )
})