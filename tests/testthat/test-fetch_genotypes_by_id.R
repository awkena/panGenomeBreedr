test_that("fetch_genotypes_by_id() works with default meta_data (all)", {
  target_ids <- c("INDEL_Chr05_75104541", "SNP_Chr05_75104557")
  
  genotypes <- fetch_genotypes_by_id(
    con = con_test,
    variant_ids = target_ids,
    meta_data = NULL 
  )
  
  # Check output structure and dimensions
  expect_s3_class(genotypes, "data.frame")
  expect_equal(nrow(genotypes), 2)
  
  # Check that all requested variants are present
  expect_true(all(target_ids %in% genotypes$variant_id))
  
  # Check for expected metadata, calculated metrics, and sample columns
  expect_true(all(c("chrom", "pos", "ref", "alt", "variant_type") %in% colnames(genotypes)))
  expect_true(all(c("major_allele", "minor_allele_freq") %in% colnames(genotypes)))
  expect_gt(ncol(genotypes), 10) # Ensure sample columns are present
})

test_that("fetch_genotypes_by_id() works with specific meta_data", {
  target_ids <- c("INDEL_Chr05_75104541", "SNP_Chr05_75104557")
  
  genotypes <- fetch_genotypes_by_id(
    con = con_test,
    variant_ids = target_ids,
    meta_data = c("chrom", "pos", "minor_allele_freq")
  )
  
  # Check output structure
  expect_s3_class(genotypes, "data.frame")
  expect_equal(nrow(genotypes), 2)
  
  # Check that only requested metadata and metrics are present
  expected_cols <- c("variant_id", "chrom", "pos", "minor_allele_freq")
  expect_true(all(expected_cols %in% colnames(genotypes)))
  
  # Check that non-requested columns are absent
  expect_false("variant_type" %in% colnames(genotypes))
  expect_false("major_allele" %in% colnames(genotypes))
})

test_that("fetch_genotypes_by_id() automatically detects chromosome", {
  target_ids <- c("INDEL_Chr05_75104541", "SNP_Chr05_75104557")
  
  # Automatically detect chcromosome number if not provided by user
  result_auto_chrom <- fetch_genotypes_by_id(con = con_test, variant_ids = target_ids, chrom = NULL)
  
  # Explicitly provide 'Chr05'
  result_explicit_chrom <- fetch_genotypes_by_id(con = con_test, variant_ids = target_ids, chrom = "Chr05")
  
  expect_equal(nrow(result_auto_chrom), 2)
  expect_equal(result_auto_chrom, result_explicit_chrom)
})

test_that("fetch_genotypes_by_id() handles edge cases and errors", {
  # Empty variant_ids vector
  expect_warning(
    empty_res <- fetch_genotypes_by_id(con = con_test, variant_ids = character(0)),
    "The 'variant_ids' vector is empty. No query was executed."
  )
  expect_s3_class(empty_res, "data.frame")
  expect_equal(nrow(empty_res), 0)
  
  # Non-existent variant_ids
  expect_warning(
    no_match_res <- fetch_genotypes_by_id(con = con_test, variant_ids = "SNP_Chr01_12345"),
    "No data found for the provided variant IDs."
  )
  expect_s3_class(no_match_res, "data.frame")
  expect_equal(nrow(no_match_res), 0)
  
  # Invalid connection
  expect_error(
    fetch_genotypes_by_id(con = NULL, variant_ids = "SNP_Chr05_75104557"),
    "A valid 'con' object is required for 'local' mode. Please run connect_local_db() first.",
    fixed = TRUE
  )
})
