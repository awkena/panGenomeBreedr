test_that("summarize_annotations() works correctly in local mode for a valid region", {
  # Use a known region with data
  summary_list <- summarize_annotations(
    con = con_test,
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )

  # 1. Check that the output is a list
  expect_type(summary_list, "list")
  # 2. Check that the list contains the three expected data frames
  expect_named(summary_list, c("annotation_summary", "impact_summary", "variant_type_totals"))

  # 3. Check each data frame's structure and content
  # annotation_summary
  expect_s3_class(summary_list$annotation_summary, "data.frame")
  expect_true(all(c("annotation", "variant_type", "count") %in% colnames(summary_list$annotation_summary)))
  expect_gt(nrow(summary_list$annotation_summary), 0)

  # impact_summary
  expect_s3_class(summary_list$impact_summary, "data.frame")
  expect_true(all(c("impact", "variant_type", "count") %in% colnames(summary_list$impact_summary)))
  expect_gt(nrow(summary_list$impact_summary), 0)

  # variant_type_totals
  expect_s3_class(summary_list$variant_type_totals, "data.frame")
  expect_true(all(c("variant_type", "total_variants") %in% colnames(summary_list$variant_type_totals)))
  expect_gt(nrow(summary_list$variant_type_totals), 0)
})

test_that("summarize_annotations() errors correctly for missing required arguments", {
  # Missing chrom
  expect_error(
    summarize_annotations(con = con_test, start = 1, end = 100),
    "You must provide 'chrom', 'start', and 'end' for the region summary.",
    fixed = TRUE
  )
  # Missing start
  expect_error(
    summarize_annotations(con = con_test, chrom = "Chr01", end = 100),
    "You must provide 'chrom', 'start', and 'end' for the region summary.",
    fixed = TRUE
  )
})

test_that("summarize_annotations() errors correctly without a valid connection", {
  expect_error(
    summarize_annotations(con = NULL, chrom = "Chr05", start = 1, end = 100),
    "A valid 'con' object is required for 'local' mode. Please run connect_local_db() first.",
    fixed = TRUE
  )
})

test_that("summarize_annotations() returns empty data frames and warns for regions with no variants", {
  expect_warning(
    empty_summary <- summarize_annotations(con = con_test, chrom = "Chr01", start = 1, end = 10),
    "No variants found in the specified region.",
    fixed = TRUE
  )
  expect_equal(nrow(empty_summary$annotation_summary), 0)
  expect_equal(nrow(empty_summary$impact_summary), 0)
  expect_equal(nrow(empty_summary$variant_type_totals), 0)
})

test_that("summarize_annotations() validates connect_db_mode argument", {
  expect_error(
    summarize_annotations(con = con_test, chrom = "Chr05", start = 1, end = 100, connect_db_mode = "invalid_mode")
  )
})