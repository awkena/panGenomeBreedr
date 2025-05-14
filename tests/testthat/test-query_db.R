test_that("query_db() retrieves correct records from SQLite database", {
  skip_on_cran()
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  # Setup temporary SQLite DB
  db_path <- tempfile(fileext = ".db")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Create variants table
  variants_df <- data.frame(
    variant_id = c("v1", "v2", "v3", "v4"),
    chrom = c("Chr01", "Chr01", "Chr02", "Chr02"),
    pos = c(100, 200, 300, 400),
    ref = c("A", "T", "C", "G"),
    alt = c("G", "A", "T", "C"),
    variant_type = c("SNP", "SNP", "INDEL", "SNP")
  )
  DBI::dbWriteTable(con, "variants", variants_df)

  # Create annotations table
  annotations_df <- data.frame(
    variant_id = c("v1", "v2", "v3"),
    gene_name = c("Sobic.001G000100", "Sobic.001G000100", "Sobic.001G000200"),
    impact = c("HIGH", "MODERATE", "LOW")
  )
  DBI::dbWriteTable(con, "annotations", annotations_df)

  # Create genotypes table
  genotypes_df <- data.frame(
    variant_id = c("v1", "v2", "v3", "v4"),
    chrom = c("Chr01", "Chr01", "Chr02", "Chr02"),
    pos = c(100, 200, 300, 400),
    sample1 = c("0/1", "1/1", "0/0", "0/1"),
    sample2 = c("0/0", "0/1", "1/1", "0/0")
  )
  DBI::dbWriteTable(con, "genotypes", genotypes_df)

  DBI::dbDisconnect(con)

  # Test query from variants
  result_variants <- query_db(db_path, table_name = "variants", chrom = "Chr01",
                              start = 50, end = 150)
  expect_true(nrow(result_variants) == 1)
  expect_equal(result_variants$variant_id, "v1")

  # Test query from annotations (gene + region)
  result_annot <- query_db(db_path, table_name = "annotations", chrom = "Chr01",
                           start = 50, end = 250, gene_name = "Sobic.001G000100")

  expect_equal(nrow(result_annot), 2)
  expect_true(all(result_annot$gene_name == "Sobic.001G000100"))

  # Test query from genotypes
  result_geno <- query_db(db_path, table_name = "genotypes", chrom = "Chr02",
                          start = 300, end = 400)

  expect_equal(nrow(result_geno), 2)
  expect_true(all(c("sample1", "sample2", "ref", "alt", "variant_type") %in% colnames(result_geno)))

  # Check for table not existing
  expect_error(query_db(db_path, table_name = "nonexistent", chrom = "Chr01",
                        start = 1, end = 100), "'arg' should be one of")
  # Clean up
  unlink(db_path)
})
