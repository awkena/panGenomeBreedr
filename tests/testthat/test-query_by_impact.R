test_that("query_by_impact() filters variants by impact and region", {
  skip_on_cran()
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  # Create temp DB
  db_path <- tempfile(fileext = ".db")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Mock variants table
  variants_df <- data.frame(
    variant_id = c("v1", "v2", "v3", "v4"),
    chrom = c("Chr01", "Chr01", "Chr02", "Chr02"),
    pos = c(100, 200, 300, 400),
    ref = c("A", "T", "C", "G"),
    alt = c("G", "A", "T", "C"),
    variant_type = c("SNP", "SNP", "INDEL", "SNP"))

  DBI::dbWriteTable(con, "variants", variants_df)

  # Mock annotations table
  annotations_df <- data.frame(
    variant_id = c("v1", "v2", "v3", "v4"),
    allele = c("G", "A", "T", "C"),
    annotation = c("missense", "synonymous", "stop_gained", "intron"),
    impact = c("HIGH", "LOW", "HIGH", "MODIFIER"),
    gene_name = c("Sobic.001G000100", "Sobic.001G000100", "Sobic.001G000200",
                  "Sobic.001G000300"),
    gene_id = paste0("gene", 1:4),
    feature_type = rep("gene", 4),
    feature_id = paste0("feat", 1:4),
    transcript_biotype = rep("protein_coding", 4),
    rank = c(1, 1, 1, 1),
    HGVS_c = c("c.100A>G", "c.200T>A", "c.300C>T", "c.400G>C"),
    HGVS_p = c("p.Lys34Glu", "p.Leu50=", "p.Arg100*", "p.Val120Ala"))

   DBI::dbWriteTable(con, "annotations", annotations_df)

  DBI::dbDisconnect(con)

  # Query by impact only
  result_high <- query_by_impact(db_path, impact_level = "HIGH")
  expect_equal(nrow(result_high), 2)
  expect_true(all(result_high$impact == "HIGH"))

  # Query by impact (lowercase input)
  result_case <- query_by_impact(db_path, impact_level = "high")
  expect_equal(nrow(result_case), 2)

  # Query with region filtering
  result_region <- query_by_impact(db_path, impact_level = "HIGH", chrom = "Chr02",
                                   start = 250, end = 350)

  expect_equal(nrow(result_region), 1)
  expect_equal(result_region$chrom, "Chr02")
  expect_equal(result_region$pos, 300)

  # Query with invalid impact level
  expect_error(query_by_impact(db_path, impact_level = "nonsense"),
               "No valid impact levels provided.")

  # Cleanup
  unlink(db_path)
})
