test_that("query_by_af() queries database and filters by allele frequency", {
  skip_on_cran()
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  # Create temporary SQLite DB
  db_path <- tempfile(fileext = ".db")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Mock variants table
  variants <- data.frame(
    variant_id = c("v1", "v2", "v3", "v4"),
    chrom = c("Chr01", "Chr01", "Chr02", "Chr02"),
    pos = c(100, 200, 300, 400),
    ref = c("A", "T", "C", "G"),
    alt = c("G", "A", "T", "C"),
    variant_type = c("SNP", "SNP", "SNP", "SNP")
  )
  DBI::dbWriteTable(con, "variants", variants)

  # Mock genotypes table
  genotypes <- data.frame(
    variant_id = c("v1", "v2", "v3", "v4"),
    chrom = c("Chr01", "Chr01", "Chr02", "Chr02"),
    pos = c(100, 200, 300, 400),
    s1 = c("0/0", "0/1", "1/1", "0/0"),
    s2 = c("0/1", "1/1", "0/0", "0/1"),
    s3 = c("1/1", "0/0", "0/1", "1/1"),
    stringsAsFactors = FALSE
  )
  DBI::dbWriteTable(con, "genotypes", genotypes)

  DBI::dbDisconnect(con)

  # Test query in Chr01 within 0–1 AF (should include v1 and v2)
  result <- query_by_af(db_path,
                        min_af = 0,
                        max_af = 1,
                        chrom = "Chr01",
                        start = 50,
                        end = 250)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_true(all(result$variant_id %in% c("v1", "v2")))

  # Test with stricter threshold — only v1 (AF = 0.5)
  filtered <- query_by_af(db_path,
                          min_af = 0.49,
                          max_af = 0.51,
                          chrom = "Chr01",
                          start = 50,
                          end = 250)
  expect_equal(nrow(filtered), 2)
  expect_equal(filtered$variant_id, c("v1", "v2"))

  # Test with region that returns no results
  empty_result <- query_by_af(db_path,
                              min_af = 0.1,
                              max_af = 0.2,
                              chrom = "Chr03",
                              start = 1,
                              end = 1000)
  expect_true(is.data.frame(empty_result))
  expect_equal(nrow(empty_result), 0)

  # Cleanup
  unlink(db_path)
})
