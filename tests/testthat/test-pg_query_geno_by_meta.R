test_that("pg_query_geno_by_meta strictly validates input parameters", {
  dummy_gt <- data.frame(variant_id = "v1", S1 = "0/0")

  # Check invalid genotype matrix
  expect_error(
    pg_query_geno_by_meta(
      genotype_matrix = "not_a_matrix",
      genotype_start_col = 2,
      filters = list(pop = "A")
    ),
    "must be a data frame or matrix"
  )

  # Check empty genotype matrix
  expect_warning(
    res_empty <- pg_query_geno_by_meta(data.frame(), 2, list(pop = "A")),
    "Input `genotype_matrix` is empty"
  )
  expect_equal(nrow(res_empty), 0)

  # Check invalid start column
  expect_error(
    pg_query_geno_by_meta(dummy_gt, genotype_start_col = 1, list(pop = "A")),
    "numeric value greater than 1"
  )

  # Check invalid filter list (unnamed)
  expect_error(
    pg_query_geno_by_meta(dummy_gt, 2, filters = list("A", "B")),
    "must be a named list"
  )
})


test_that("pg_query_geno_by_meta handles exhaustive filters and missing sample columns", {
  dummy_gt <- data.frame(
    variant_id = "v1",
    S1 = "0/0",
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    pg_get_sample_metadata = function(...) {
      data.frame(
        lib = c("S1", "S2"),
        country = c("Ghana", "Togo"),
        stringsAsFactors = FALSE
      )
    }
  )

  # Test filter that eliminates all samples from the metadata
  expect_message(
    res_no_match <- pg_query_geno_by_meta(dummy_gt, 2, list(country = "Mali")),
    "No samples found matching all combined metadata criteria"
  )
  expect_equal(nrow(res_no_match), 0)

  # Test filter that keeps a sample (S2) which doesn't exist in dummy_gt matrix
  expect_message(
    res_missing_col <- pg_query_geno_by_meta(
      dummy_gt,
      2,
      list(country = "Togo")
    ),
    "None of the matching samples were found in the provided genotype matrix"
  )
  expect_equal(nrow(res_missing_col), 0)
})


test_that("pg_query_geno_by_meta correctly subsets samples and routes to metrics calculator", {
  # Build dummy genotype matrix with 3 samples
  dummy_gt <- data.frame(
    variant_id = c("v1", "v2"),
    ref = c("A", "C"),
    alt = c("T", "G"),
    S1 = c("0/0", "1/1"),
    S2 = c("0/1", "0/0"),
    S3 = c("1/1", "0/1"),
    stringsAsFactors = FALSE
  )

  # Mock metadata fetch to classify the 3 samples
  local_mocked_bindings(
    pg_get_sample_metadata = function(...) {
      data.frame(
        lib = c("S1", "S2", "S3"),
        population = c("PopA", "PopB", "PopA"),
        country = c("Ghana", "Ghana", "Togo"),
        stringsAsFactors = FALSE
      )
    },

    # Mock the downstream calculator to attach a flag proving it executed
    append_locus_allele_metrics = function(region_matrix, meta_data) {
      region_matrix$recalculated <- TRUE
      region_matrix
    }
  )

  # Execute with multiple filter criteria
  res <- pg_query_geno_by_meta(
    genotype_matrix = dummy_gt,
    genotype_start_col = 4,
    filters = list(population = "PopA", country = "Ghana")
  )

  # Validate output structure
  cols <- colnames(res)
  expect_true("variant_id" %in% cols)
  expect_true("S1" %in% cols)
  expect_false("S2" %in% cols)
  expect_false("S3" %in% cols)

  # Validate the routing hit the metrics calculator
  expect_true(res$recalculated[1])
})


test_that("pg_query_geno_by_meta warns and gracefully skips invalid filter columns", {
  dummy_gt <- data.frame(
    variant_id = "v1",
    S1 = "0/0",
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    pg_get_sample_metadata = function(...) {
      data.frame(lib = "S1", country = "Ghana", stringsAsFactors = FALSE)
    },
    append_locus_allele_metrics = function(region_matrix, meta_data) {
      region_matrix
    }
  )

  # Pass a filter column ("climate") that doesn't exist in the metadata
  expect_warning(
    res <- pg_query_geno_by_meta(
      dummy_gt,
      2,
      list(climate = "Arid", country = "Ghana")
    ),
    "Column 'climate' not found in metadata. Skipping this filter."
  )

  # Verify it still successfully applied the valid filter and returned the data
  expect_true("S1" %in% colnames(res))
})