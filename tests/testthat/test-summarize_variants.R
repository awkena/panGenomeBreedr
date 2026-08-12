test_that("summarize_variants() works in local mode with annotations (default)", {
  stats <- summarize_variants(con = con_test, connect_db_mode = 'local')

  expect_s3_class(stats, "data.frame")
  expect_gt(nrow(stats), 0)

  expected_cols <- c("chrom", "n_variants", "min_pos", "max_pos", "n_unique_ids", "n_annotated")
  expect_true(all(expected_cols %in% colnames(stats)))

  expect_true(is.numeric(stats$n_variants))
  expect_true(is.numeric(stats$n_annotated))
  expect_false(any(is.na(stats$n_annotated)))
})

test_that("summarize_variants() works in local mode without annotations", {
  stats <- summarize_variants(con = con_test, include_annotations = FALSE)

  expect_s3_class(stats, "data.frame")
  expect_gt(nrow(stats), 0)

  expected_cols <- c("chrom", "n_variants", "min_pos", "max_pos", "n_unique_ids")
  expect_true(all(expected_cols %in% colnames(stats)))
  expect_false("n_annotated" %in% colnames(stats))
})

test_that("summarize_variants() errors correctly without a valid connection", {
  expect_error(
    summarize_variants(con = NULL, connect_db_mode = 'local'),
    "A valid 'con' object is required for 'local' mode. Please run connect_local_db() first.",
    fixed = TRUE
  )
})

test_that("summarize_variants() validates connect_db_mode argument", {
  expect_error(summarize_variants(con = con_test, connect_db_mode = "invalid_mode"))
})
