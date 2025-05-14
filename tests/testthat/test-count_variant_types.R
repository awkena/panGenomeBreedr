test_that("count_variant_types() returns correct counts of variant types", {
  skip_on_cran()
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  # Create temporary SQLite DB
  db_path <- tempfile(fileext = ".db")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Create variants table with variant_type
  variants <- data.frame(
    variant_id = paste0("v", 1:6),
    chrom = rep("Chr01", 6),
    pos = seq(100, 600, by = 100),
    ref = rep("A", 6),
    alt = rep("T", 6),
    variant_type = c("SNP", "SNP", "INDEL", "SNP", "INDEL", "SNP"),
    stringsAsFactors = FALSE)

  DBI::dbWriteTable(con, "variants", variants)

  DBI::dbDisconnect(con)

  # Run function
  result <- count_variant_types(db_path)

  # Validate output
  expect_true(is.data.frame(result))
  expect_named(result, c("variant_type", "n"))
  expect_equal(nrow(result), 2)
  expect_equal(result$n[result$variant_type == "SNP"], 4)
  expect_equal(result$n[result$variant_type == "INDEL"], 2)

  # Test error for missing 'variant_type' column
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  DBI::dbExecute(con, "DROP TABLE variants")
  DBI::dbWriteTable(con, "variants", variants[, -6])  # omit variant_type
  DBI::dbDisconnect(con)

  expect_error(count_variant_types(db_path),
               "The table does not have a 'variant_type' column")

  # Cleanup
  unlink(db_path)
})
