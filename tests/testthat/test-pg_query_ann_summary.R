test_that("pg_query_ann_summary enforces required coordinate arguments", {
  # Verify that omission of any critical positional arguments triggers the stop condition
  expect_error(
    pg_query_ann_summary(start = 75104537, end = 75106403),
    "You must provide 'chrom', 'start', and 'end'"
  )

  expect_error(
    pg_query_ann_summary(chrom = "Chr05", start = 75104537),
    "You must provide 'chrom', 'start', and 'end'"
  )
})


test_that("pg_query_ann_summary successfully parses standard region summaries", {
  local_mocked_bindings(
    .api_fetch = function(endpoint, query, simplify) {
      # Assert the API is routing to the correct endpoint
      expect_equal(endpoint, "/stats/ann_summary")

      # Assert the query parameters are packaged correctly
      expect_equal(query$chrom, "Chr05")
      expect_equal(query$start, 100)
      expect_equal(query$end, 500)
      expect_equal(query$annotations_table, "annotations")
      expect_equal(query$variants_table, "variants")

      # Return a well-formed dummy response
      list(
        annotation_summary = data.frame(
          feature = c("gene", "transcript"),
          count = c(12, 15)
        ),
        impact_summary = data.frame(
          impact = c("HIGH", "MODERATE"),
          count = c(2, 45)
        ),
        variant_type_totals = data.frame(
          type = c("SNP", "INDEL"),
          count = c(150, 12)
        )
      )
    }
  )

  res <- pg_query_ann_summary(chrom = "Chr05", start = 100, end = 500)

  # Validate that the returned object contains the expected structured data frames
  expect_type(res, "list")
  expect_named(
    res,
    c("annotation_summary", "impact_summary", "variant_type_totals")
  )
  expect_s3_class(res$impact_summary, "data.frame")
  expect_equal(nrow(res$variant_type_totals), 2)
  expect_equal(res$impact_summary$impact[1], "HIGH")
})


test_that("pg_query_ann_summary safely handles empty API returns", {
  local_mocked_bindings(
    .api_fetch = function(endpoint, query, simplify) {
      # Return an empty list, simulating what jsonlite does for empty JSON arrays
      list()
    }
  )

  # Check that it throws the expected warning
  expect_warning(
    res <- pg_query_ann_summary(chrom = "Chr05", start = 1, end = 100),
    "No variants found in the specified region from the API."
  )

  # Validate that the function gracefully falls back to empty data frames
  expect_type(res, "list")
  expect_named(
    res,
    c("annotation_summary", "impact_summary", "variant_type_totals")
  )
  expect_s3_class(res$annotation_summary, "data.frame")
  expect_equal(nrow(res$variant_type_totals), 0)
})