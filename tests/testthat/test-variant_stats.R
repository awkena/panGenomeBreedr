test_that("variant_stats() returns correct statistics with and without annotations", {
  # Skip check if the global test connection wasn't initialized cleanly
  skip_if_not(
    exists("con_test"),
    message = "Global test connection 'con_test' not found."
  )

  #  Run with include_annotations = TRUE 
  stats_with_annot <- variant_stats(con = con_test, include_annotations = TRUE)

  expect_true(is.data.frame(stats_with_annot))
  expect_true(nrow(stats_with_annot) >= 1) # Asserts it captured rows for active chromosomes (e.g., Chr05)

  # Ensure all summary calculation headers are present, including the annotation count
  required_headers_annot <- c(
    "chrom",
    "n_variants",
    "min_pos",
    "max_pos",
    "n_unique_ids",
    "n_annotated"
  )
  expect_true(all(required_headers_annot %in% names(stats_with_annot)))

  # Ensure positions and counts make arithmetic sense
  expect_true(all(stats_with_annot$n_variants >= 0))
  expect_true(all(stats_with_annot$n_annotated >= 0))
  expect_true(all(stats_with_annot$max_pos >= stats_with_annot$min_pos))

  # Run with include_annotations
  stats_no_annot <- variant_stats(con = con_test, include_annotations = FALSE)

  expect_true(is.data.frame(stats_no_annot))
  expect_equal(nrow(stats_no_annot), nrow(stats_with_annot))

  # Ensure core headers are retained, but annotation summary columns are cleanly excluded
  required_headers_no_annot <- c(
    "chrom",
    "n_variants",
    "min_pos",
    "max_pos",
    "n_unique_ids"
  )
  expect_true(all(required_headers_no_annot %in% names(stats_no_annot)))
  expect_false("n_annotated" %in% names(stats_no_annot))
})