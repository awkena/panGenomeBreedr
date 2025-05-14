test_that("list_table_columns() returns correct column info for allowed tables", {
  skip_on_cran()
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  # Create temporary SQLite database
  db_path <- tempfile(fileext = ".db")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Write mock tables with known structure
  DBI::dbWriteTable(con, "variants",
                    data.frame(variant_id = character(), chrom = character(), pos = integer()))
  DBI::dbWriteTable(con, "annotations",
                    data.frame(variant_id = character(), effect = character()))
  DBI::dbWriteTable(con, "genotypes",
                    data.frame(sample_id = character(), variant_id = character(), genotype = character()))

  DBI::dbDisconnect(con)

  # Test "variants" table structure
  variant_info <- list_table_columns(db_path, "variants")
  expect_true(is.data.frame(variant_info))
  expect_equal(variant_info$name, c("variant_id", "chrom", "pos"))
  expect_equal(variant_info$type, c("TEXT", "TEXT", "INTEGER"))

  # Test "annotations" table structure
  annot_info <- list_table_columns(db_path, "annotations")
  expect_equal(annot_info$name, c("variant_id", "effect"))
  expect_equal(annot_info$type, c("TEXT", "TEXT"))

  # Test "genotypes" table structure
  geno_info <- list_table_columns(db_path, "genotypes")
  expect_equal(geno_info$name, c("sample_id", "variant_id", "genotype"))
  expect_equal(geno_info$type, rep("TEXT", 3))

  # Check error on invalid table name (not in match.arg)
  expect_error(list_table_columns(db_path, table_name = "nonexistent"),
               "'arg' should be one of")

  # Cleanup
  unlink(db_path)
})
