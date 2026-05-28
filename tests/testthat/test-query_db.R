test_that("query_db() retrieves correct records from the database pipeline", {
  # Skip check if the global test connection wasn't initialized cleanly
  skip_if_not(
    exists("con_test"),
    message = "Global test connection 'con_test' not found."
  )

  # Test coordinate-window query from variants table
  result_variants <- query_db(
    con = con_test,
    table_name = "variants",
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )

  expect_true(is.data.frame(result_variants))
  if (nrow(result_variants) > 0) {
    expect_true(all(result_variants$chrom == "Chr05"))
    expect_true("variant_id" %in% colnames(result_variants))
  }

  # Test query from annotations table within the locus region
  result_annot <- query_db(
    con = con_test,
    table_name = "annotations",
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )

  expect_true(is.data.frame(result_annot))
  if (nrow(result_annot) > 0) {
    expect_true("variant_id" %in% colnames(result_annot))
    expect_true("annotation" %in% colnames(result_annot))
  }

  # Test multi-column query from genotypes table (Validates cross-table metadata appends)
  result_geno <- query_db(
    con = con_test,
    table_name = "genotypes",
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )

  expect_true(is.data.frame(result_geno))
  if (nrow(result_geno) > 0) {
    # Ensure wide genotype properties are pulled out alongside the underlying matrix matrix
    expect_true(all(
      c("variant_id", "chrom", "pos", "ref", "alt", "variant_type") %in%
        colnames(result_geno)
    ))
  }

  # Check for strict error handling when calling invalid tables (Verifies match.arg constraints)
  expect_error(
    query_db(
      con = con_test,
      table_name = "nonexistent",
      chrom = "Chr05",
      start = 75104537,
      end = 75106403
    ),
    "'arg' should be one of"
  )
})