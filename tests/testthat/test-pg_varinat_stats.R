test_that("pg_variant_stats queries the correct endpoint and passes parameters", {
  # --- TEST 1: include_annotations = TRUE (Default) ---
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Assert the endpoint is correct
      expect_equal(endpoint, "/stats/variants")

      # Assert the parameter was passed as TRUE
      expect_true(query$include_annotations)

      # Return a dummy response simulating the database
      data.frame(
        chrom = c("Chr01", "Chr02"),
        variant_count = c(500, 300),
        annotation_count = c(1500, 900)
      )
    }
  )

  res_true <- pg_variant_stats(include_annotations = TRUE)
  expect_s3_class(res_true, "data.frame")
  expect_equal(nrow(res_true), 2)
  expect_true("annotation_count" %in% colnames(res_true))

  # --- TEST 2: include_annotations = FALSE ---
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Assert the parameter was passed as FALSE this time
      expect_false(query$include_annotations)

      # Return a dummy response without the annotation columns
      data.frame(
        chrom = c("Chr01", "Chr02"),
        variant_count = c(500, 300)
      )
    }
  )

  res_false <- pg_variant_stats(include_annotations = FALSE)
  expect_s3_class(res_false, "data.frame")
  expect_false("annotation_count" %in% colnames(res_false))
})