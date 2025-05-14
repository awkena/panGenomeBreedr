test_that("variant_stats() returns correct statistics with and without annotations", {
  skip_on_cran()
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  # Create temporary SQLite database
  db_path <- tempfile(fileext = ".db")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Create dummy 'variants' table
  variants_df <- data.frame(
    variant_id = paste0("var", 1:6),
    chrom = c("Chr01", "Chr01", "Chr01", "Chr02", "Chr02", "Chr03"),
    pos = c(100, 200, 300, 400, 500, 600))
  DBI::dbWriteTable(con, "variants", variants_df)

  # Create dummy 'annotations' table (partial overlap with variants)
  annotations_df <- data.frame(
    variant_id = c("var1", "var3", "var4", "var6"),
    annotation = c("missense", "synonymous", "nonsense", "intron"))
  DBI::dbWriteTable(con, "annotations", annotations_df)

  DBI::dbDisconnect(con)

  # Run variant_stats() with include_annotations = TRUE
  stats_with_annot <- variant_stats(db_path, include_annotations = TRUE)

  expect_true(is.data.frame(stats_with_annot))
  expect_equal(nrow(stats_with_annot), 3)
  expect_true(all(c("chrom", "n_variants", "min_pos", "max_pos", "n_unique_ids", "n_annotated") %in% names(stats_with_annot)))
  expect_equal(stats_with_annot$n_annotated[stats_with_annot$chrom == "Chr01"], 2)

  # Run variant_stats() with include_annotations = FALSE
  stats_no_annot <- variant_stats(db_path, include_annotations = FALSE)

  expect_true(is.data.frame(stats_no_annot))
  expect_equal(nrow(stats_no_annot), 3)
  expect_true(all(c("chrom", "n_variants", "min_pos", "max_pos", "n_unique_ids") %in% names(stats_no_annot)))
  expect_false("n_annotated" %in% names(stats_no_annot))

  # Cleanup
  unlink(db_path)
})
