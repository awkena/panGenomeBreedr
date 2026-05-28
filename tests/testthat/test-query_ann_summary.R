test_that("query_ann_summary() returns correct summaries by variant type", {
  # Skip check if the global test connection wasn't initialized cleanly
  skip_if_not(
    exists("con_test"),
    message = "Global test connection 'con_test' not found."
  )

  # Run query against the real extracted window on Chr05
  summary <- query_ann_summary(
    con = con_test,
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )

  # Check structure of output list
  expect_type(summary, "list")
  expect_named(
    summary,
    c("annotation_summary", "impact_summary", "variant_type_totals")
  )

  #  Structural Checks on Annotation Summary
  expect_s3_class(summary$annotation_summary, "data.frame")
  expect_true(all(
    c("annotation", "variant_type", "count") %in%
      colnames(summary$annotation_summary)
  ))
  expect_true(nrow(summary$annotation_summary) >= 0)

  # Structural Checks on Impact Summary
  expect_s3_class(summary$impact_summary, "data.frame")
  expect_true(all(
    c("impact", "variant_type", "count") %in% colnames(summary$impact_summary)
  ))
  expect_true(nrow(summary$impact_summary) >= 0)

  # Structural Checks on Global Variant Type Totals
  expect_s3_class(summary$variant_type_totals, "data.frame")
  expect_true(all(
    c("variant_type", "total_variants") %in%
      colnames(summary$variant_type_totals)
  ))
  expect_true(nrow(summary$variant_type_totals) > 0) # Asserts that real data was found and summarized
})