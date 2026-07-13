library(testthat)

test_that("pg_query_db formats queries correctly and bypasses metrics for variants", {
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Assert the API is routing to the correct endpoint
      expect_equal(endpoint, "/db/query")

      # Assert the parameters are bundled correctly into the query list
      expect_equal(query$table_name, "variants")
      expect_equal(query$chrom, "Chr05")
      expect_equal(query$start, 100)
      expect_equal(query$end, 500)
      expect_equal(query$gene_name, "Sobic.005G123")

      # Return a dummy variant data frame
      data.frame(variant_id = c("v1", "v2"), chrom = c("Chr05", "Chr05"))
    },

    # If the function accidentally calls the metric calculator, force a failure
    append_locus_allele_metrics = function(region_matrix, meta_data = NULL) {
      stop("Metrics calculator should not be called for the variants table.")
    }
  )

  res <- pg_query_db(
    table_name = "variants",
    chrom = "Chr05",
    start = 100,
    end = 500,
    gene_name = "Sobic.005G123"
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
})


test_that("pg_query_db correctly routes genotypes through the allele metrics calculator", {
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Return a raw, "unprocessed" dummy genotype matrix
      data.frame(variant_id = "v1", ref = "A", alt = "T", sample_1 = "0/1")
    },

    append_locus_allele_metrics = function(region_matrix, meta_data = NULL) {
      # Add a flag to the matrix to prove this function was successfully executed
      region_matrix$metrics_were_calculated <- TRUE
      return(region_matrix)
    }
  )

  # Trigger the genotype specific routing
  res <- pg_query_db(table_name = "genotypes", chrom = "Chr05")

  # Verify the returned object contains the flag from the metrics calculator
  expect_true(res$metrics_were_calculated)
})


test_that("pg_query_db enforces valid table names", {
  # Verify that the internal match.arg() successfully blocks bad table requests
  expect_error(
    pg_query_db(table_name = "invalid_table_name"),
    "'arg' should be one of"
  )
})