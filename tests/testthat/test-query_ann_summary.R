test_that("query_ann_summary() returns correct summaries by variant type", {
  skip_on_cran()
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  # Setup temporary SQLite DB
  db_path <- tempfile(fileext = ".db")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Mock variants table
  variants <- data.frame(
    variant_id = c("v1", "v2", "v3"),
    chrom = c("Chr01", "Chr01", "Chr01"),
    pos = c(100, 150, 180),
    variant_type = c("SNP", "SNP", "INDEL"))

  DBI::dbWriteTable(con, "variants", variants)

  # Mock annotations table
  annotations <- data.frame(
    variant_id = c("v1", "v2", "v3"),
    annotation = c("missense_variant", "synonymous_variant", "stop_gained"),
    impact = c("HIGH", "LOW", "HIGH"))

  DBI::dbWriteTable(con, "annotations", annotations)

  DBI::dbDisconnect(con)

  # Run query
  summary <- query_ann_summary(db_path,
                               chrom = "Chr01",
                               start = 100,
                               end = 180)

  # Check structure of output
  expect_type(summary, "list")
  expect_named(summary, c("annotation_summary", "impact_summary", "variant_type_totals"))

  # Check counts
  expect_equal(nrow(summary$annotation_summary), 3)
  expect_true("missense_variant" %in% summary$annotation_summary$annotation)
  expect_true("stop_gained" %in% summary$annotation_summary$annotation)

  expect_equal(nrow(summary$impact_summary), 3)
  expect_equal(summary$impact_summary$count[summary$impact_summary$impact == "HIGH"],
               c(1,1))

  expect_equal(nrow(summary$variant_type_totals), 2)
  expect_equal(summary$variant_type_totals$total_variants[summary$variant_type_totals$variant_type == "SNP"], 2)

  # Cleanup
  unlink(db_path)
})
