test_that("variant_impact_summary() returns impact summary in wide format", {
  skip_on_cran()
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  # Create temporary SQLite database
  db_path <- tempfile(fileext = ".db")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Write mock variants table
  variants_df <- data.frame(
    variant_id = paste0("var", 1:6),
    chrom = c("Chr01", "Chr01", "Chr02", "Chr02", "Chr03", "Chr03"),
    pos = seq(100, 600, by = 100))
  DBI::dbWriteTable(con, "variants", variants_df)

  # Write mock annotations table
  annotations_df <- data.frame(
    variant_id = c("var1", "var2", "var3", "var4", "var5"),
    impact = c("MODERATE", "HIGH", "LOW", "MODERATE", "HIGH"))
  DBI::dbWriteTable(con, "annotations", annotations_df)

  DBI::dbDisconnect(con)

  # Run the function
  summary_df <- variant_impact_summary(db_path)

  # Check structure
  expect_true(is.data.frame(summary_df))
  expect_true("chrom" %in% names(summary_df))
  expect_true(any(grepl("impact_", names(summary_df))))

  # Check row count = unique chromosomes in joined data
  expect_equal(nrow(summary_df), 3)  # Chr01, Chr02, Chr03

  # Check counts in output
  chr01_row <- summary_df[summary_df$chrom == "Chr01", ]
  expect_equal(chr01_row$impact_MODERATE, 1)
  expect_equal(chr01_row$impact_HIGH, 1)

  chr03_row <- summary_df[summary_df$chrom == "Chr03", ]
  expect_equal(chr03_row$impact_HIGH, 1)
  expect_equal(chr03_row$impact_MODERATE, 0)  # should be filled as 0

  # Cleanup
  unlink(db_path)
})
