test_that("pg_query_genotypes safely handles missing or empty variant inputs", {
  # Check that providing an empty vector triggers the warning and returns an empty data frame
  expect_warning(
    res_empty <- pg_query_genotypes(c()),
    "vector is empty"
  )
  expect_s3_class(res_empty, "data.frame")
  expect_equal(nrow(res_empty), 0)

  # Check that providing no arguments at all triggers the same safeguard
  expect_warning(
    res_missing <- pg_query_genotypes(),
    "vector is empty"
  )
  expect_equal(nrow(res_missing), 0)
})

test_that("pg_query_genotypes handles empty database returns gracefully", {
  # Intercept the API to simulate the database finding no matches for the requested IDs
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      data.frame()
    }
  )

  # Verify the function alerts the user and halts further processing safely
  expect_warning(
    res <- pg_query_genotypes(c("INVALID_SNP")),
    "No data found for the provided variant IDs"
  )
  expect_equal(nrow(res), 0)
})

test_that("pg_query_genotypes processes default requests (all metadata)", {
  # Mock the network fetch and the downstream allele metrics calculator
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Verify the ID vector was squashed correctly
      expect_equal(query$variant_ids, "V1,V2")
      # Verify meta_data is NULL when requesting all data
      expect_null(query$meta_data)

      data.frame(
        variant_id = c("V1", "V2"),
        chrom = c("Chr1", "Chr1"),
        ref = c("A", "G"),
        alt = c("T", "C"),
        Sample_1 = c("0/0", "0/1"),
        stringsAsFactors = FALSE
      )
    },

    append_locus_allele_metrics = function(region_matrix, meta_data) {
      # Simulate the calculator appending a metric
      region_matrix$major_allele_freq <- c(1.0, 0.5)
      region_matrix
    }
  )

  res <- pg_query_genotypes(c("V1", "V2"))

  # Validate that the metrics were appended and the sample columns were retained
  expect_equal(nrow(res), 2)
  expect_true("major_allele_freq" %in% colnames(res))
  expect_true("Sample_1" %in% colnames(res))
})

test_that("pg_query_genotypes correctly parses and routes subsetted metadata requests", {
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      expect_equal(query$meta_data, "variant_id,ref,alt,chrom")

      data.frame(
        variant_id = "V1",
        ref = "A",
        alt = "T",
        chrom = "Chr1",
        Sample_1 = "0/1",
        stringsAsFactors = FALSE
      )
    },

    append_locus_allele_metrics = function(region_matrix, meta_data) {
      # Simulate the metric calculation
      region_matrix$minor_allele_freq <- 0.5
      region_matrix
    }
  )

  # Request a mix of a standard database column and a calculated metric column
  res <- pg_query_genotypes(
    variant_ids = "V1",
    meta_data = c("chrom", "minor_allele_freq")
  )

  # Validate column reconstruction
  cols <- colnames(res)
  expect_true(all(
    c("variant_id", "chrom", "minor_allele_freq", "Sample_1") %in% cols
  ))
  expect_false("ref" %in% cols)
  expect_false("alt" %in% cols)

  expect_equal(res$minor_allele_freq[1], 0.5)
})