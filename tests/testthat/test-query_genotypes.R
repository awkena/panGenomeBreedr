test_that("query_genotypes() retrieves genotypes and metadata from SQLite", {
  skip_on_cran()
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  # Setup temp SQLite DB
  db_path <- tempfile(fileext = ".db")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Create mock variants table
  variants <- data.frame(
    variant_id = c("v1", "v2", "v3"),
    chrom = c("Chr01", "Chr01", "Chr02"),
    pos = c(100, 200, 300),
    ref = c("A", "T", "C"),
    alt = c("G", "A", "T"),
    variant_type = c("SNP", "SNP", "INDEL"),
    extra = c("x", "y", "z")  # included for testing metadata filtering
  )
  DBI::dbWriteTable(con, "variants", variants)

  # Create mock genotypes table
  genotypes <- data.frame(
    variant_id = c("v1", "v2", "v3"),
    sample1 = c("0/0", "0/1", "1/1"),
    sample2 = c("1/1", "0/0", "0/1"),
    stringsAsFactors = FALSE
  )
  DBI::dbWriteTable(con, "genotypes", genotypes)

  DBI::dbDisconnect(con)

  # ---- TEST 1: Full query with valid metadata ----
  result <- query_genotypes(db_path,
                            variant_ids = c("v1", "v3"),
                            meta_data = c("chrom", "pos", "ref", "alt", "variant_type"))

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_true(all(c("chrom", "pos", "ref", "alt", "variant_type", "sample1", "sample2") %in% names(result)))

  # ---- TEST 2: Default metadata (meta_data = NULL) ----
  result_default <- query_genotypes(db_path,
                                    variant_ids = c("v2", "v3"))
  expect_true("chrom" %in% names(result_default))
  expect_equal(nrow(result_default), 2)

  # ---- TEST 3: Missing metadata column (warns but proceeds) ----
  expect_warning({
    result_partial <- query_genotypes(db_path,
                                      variant_ids = c("v1"),
                                      meta_data = c("chrom", "nonexistent_col"))
  }, "do not exist in the 'variants' table")

  expect_true("chrom" %in% names(result_partial))
  expect_false("nonexistent_col" %in% names(result_partial))

  # ---- TEST 4: No variant IDs provided (error) ----
  expect_error(query_genotypes(db_path, variant_ids = character(0)),
               "at least one variant ID")

  # ---- TEST 5: Nonexistent variant ID (warns and returns empty) ----
  expect_warning({
    result_empty <- query_genotypes(db_path, variant_ids = "v999")
  }, "No genotype data found")

  expect_true(is.data.frame(result_empty))
  expect_equal(nrow(result_empty), 0)

  # Cleanup
  unlink(db_path)
})
