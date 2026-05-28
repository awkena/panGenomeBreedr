test_that("query_genotypes() retrieves genotypes and metadata from the database", {
  # Skip check if the global test connection wasn't initialized cleanly
  skip_if_not(
    exists("con_test"),
    message = "Global test connection 'con_test' not found."
  )

  # Inspect available variant IDs present in your real mini database view slice
  # to fetch valid target keys dynamically for testing
  sample_variants <- DBI::dbGetQuery(
    con_test,
    "SELECT variant_id FROM variants LIMIT 2"
  )$variant_id

  skip_if(
    length(sample_variants) < 2,
    message = "Insufficient variant data rows in mini dataset to execute test."
  )

  v1 <- sample_variants[1]
  v2 <- sample_variants[2]

  # Full query with explicit valid metadata parameters 
  result <- query_genotypes(
    con = con_test,
    variant_ids = c(v1, v2),
    meta_data = c("chrom", "pos", "ref", "alt", "variant_type")
  )

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_true(all(
    c("variant_id", "chrom", "pos", "ref", "alt", "variant_type") %in%
      names(result)
  ))

  # Default metadata evaluation path (meta_data = NULL) 
  result_default <- query_genotypes(con = con_test, variant_ids = c(v1, v2))

  expect_true(is.data.frame(result_default))
  expect_true("chrom" %in% names(result_default))
  expect_equal(nrow(result_default), 2)

  # No variant IDs provided (triggers defensive warning & returns empty df) 
  expect_warning(
    {
      result_empty_input <- query_genotypes(
        con = con_test,
        variant_ids = character(0)
      )
    },
    "The 'variant_ids' vector is empty"
  )

  expect_true(is.data.frame(result_empty_input))
  expect_equal(nrow(result_empty_input), 0)

  #  Nonexistent variant ID handling (warns and returns empty df) 
  expect_warning(
    {
      result_empty_locus <- query_genotypes(
        con = con_test,
        variant_ids = "NONEXISTENT_MARKER_ID_999"
      )
    },
    "No data found for the provided variant IDs"
  )

  expect_true(is.data.frame(result_empty_locus))
  expect_equal(nrow(result_empty_locus), 0)
})