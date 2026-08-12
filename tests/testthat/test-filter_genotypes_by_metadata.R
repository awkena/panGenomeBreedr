gt_matrix <- fetch_table_region(
  con = con_test,
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
)

# In the test data, the first 10 columns are metadata, so genotypes start at column 11.
genotype_start_col_test <- 11

test_that("filter_genotypes_by_metadata() works with a single filter", {
  filters <- list(countryorigin = "Ethiopia")

  filtered_gt <- filter_genotypes_by_metadata(
    con = con_test,
    genotype_matrix = gt_matrix,
    genotype_start_col = genotype_start_col_test,
    filters = filters
  )

  # Check output structure
  expect_s3_class(filtered_gt, "data.frame")
  expect_gt(nrow(filtered_gt), 0)

  # Check that sample columns were filtered
  expect_lt(ncol(filtered_gt), ncol(gt_matrix))

  # Check that allele frequencies were recalculated and are different
  original_maf <- gt_matrix$minor_allele_freq
  recalculated_maf <- filtered_gt$minor_allele_freq
  expect_false(identical(original_maf, recalculated_maf))
})

test_that("filter_genotypes_by_metadata() works with multiple filters", {
  filters <- list(countryorigin = "Ethiopia", population = "Gates")

  filtered_gt <- filter_genotypes_by_metadata(
    con = con_test,
    genotype_matrix = gt_matrix,
    genotype_start_col = genotype_start_col_test,
    filters = filters
  )

  expect_s3_class(filtered_gt, "data.frame")
  expect_gt(nrow(filtered_gt), 0)
  expect_lt(ncol(filtered_gt), ncol(gt_matrix))
})

test_that("filter_genotypes_by_metadata() handles no matching samples", {
  filters <- list(countryorigin = "Atlantis") # A value not in the dataset

  expect_message(
    no_match <- filter_genotypes_by_metadata(
      con = con_test, genotype_matrix = gt_matrix,
      genotype_start_col = genotype_start_col_test, filters = filters
    ),
    "No samples found matching all combined metadata criteria."
  )

  expect_s3_class(no_match, "data.frame")
  expect_equal(nrow(no_match), 0)
})

test_that("filter_genotypes_by_metadata() handles invalid inputs", {
  expect_error(filter_genotypes_by_metadata(con = NULL, genotype_matrix = gt_matrix, genotype_start_col = 11, filters = list(a = 1)), "A valid 'con' object is required")
  expect_error(filter_genotypes_by_metadata(con = con_test, genotype_matrix = "not a df", genotype_start_col = 11, filters = list(a = 1)), "`genotype_matrix` must be a data frame or matrix.")
  expect_error(filter_genotypes_by_metadata(con = con_test, genotype_matrix = gt_matrix, genotype_start_col = 1, filters = list(a = 1)), "`genotype_start_col` must be a numeric value greater than 1.")
  expect_error(filter_genotypes_by_metadata(con = con_test, genotype_matrix = gt_matrix, genotype_start_col = 11, filters = "not a list"), "`filters` must be a named list.")
})

test_that("filter_genotypes_by_metadata() warns on invalid filter column", {
  filters <- list(countryorigin = "Ethiopia", invalid_col = "test")

  expect_warning(
    filtered_gt <- filter_genotypes_by_metadata(
      con = con_test, genotype_matrix = gt_matrix,
      genotype_start_col = genotype_start_col_test, filters = filters
    ),
    "Column 'invalid_col' not found in metadata. Skipping this filter."
  )

  # Should still filter by the valid column and produce a result
  expect_gt(nrow(filtered_gt), 0)
  expect_lt(ncol(filtered_gt), ncol(gt_matrix))
})
