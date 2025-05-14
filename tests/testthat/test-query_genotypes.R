test_that("query_genotypes() returns sample x variant genotype matrix", {
  skip_on_cran()
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  # Create temp DB
  db_path <- tempfile(fileext = ".db")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Mock genotype table (wide format: one row per variant)
  genotypes <- data.frame(
    variant_id = c("v1", "v2"),
    chrom = c("Chr01", "Chr01"),
    pos = c(100, 200),
    sample1 = c("0/0", "0/1"),
    sample2 = c("1/1", "0/0"),
    sample3 = c("0/1", "1/1"),
    stringsAsFactors = FALSE)

  DBI::dbWriteTable(con, "genotypes", genotypes)

  DBI::dbDisconnect(con)

  # Query both variants
  result <- query_genotypes(db_path, variant_ids = c("v1", "v2"))

  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 3)  # 1 for sample + 2 variant columns
  expect_equal(nrow(result), 3)  # 3 samples
  expect_named(result, c("sample", "v1", "v2"))

  # Check transposition
  expect_equal(result$sample, c("sample1", "sample2", "sample3"))
  expect_equal(result$v1[1], "0/0")  # sample1, v1

  # Test with a single variant ID
  single_result <- query_genotypes(db_path, variant_ids = "v2")
  expect_equal(ncol(single_result), 2)  # sample + v2
  expect_equal(nrow(single_result), 3)
  expect_named(single_result, c("sample", "v2"))

  # Error if variant_ids is empty
  expect_error(query_genotypes(db_path, variant_ids = character(0)),
               "Please provide at least one variant ID")

  # Cleanup
  unlink(db_path)
})
