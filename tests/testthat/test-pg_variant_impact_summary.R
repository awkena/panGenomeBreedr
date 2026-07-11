test_that("pg_variant_impact_summary queries correct endpoint and returns expected schema", {
  # Mock the internal fetcher
  local_mocked_bindings(
    .api_fetch = function(endpoint, query = NULL, simplify = TRUE) {
      # Assert the function targets the correct route
      expect_equal(endpoint, "/stats/impact")

      # Return a dummy response mimicking the wide-format impact summary
      data.frame(
        chrom = c("Chr01", "Chr02"),
        HIGH = c(15, 8),
        MODERATE = c(120, 95),
        LOW = c(300, 250),
        MODIFIER = c(4500, 3200),
        stringsAsFactors = FALSE
      )
    }
  )

  # Execute the function
  res <- pg_variant_impact_summary()

  # Validate the output structure
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(
    c("chrom", "HIGH", "MODERATE", "LOW", "MODIFIER") %in% colnames(res)
  ))
  expect_equal(res$HIGH[1], 15)
})