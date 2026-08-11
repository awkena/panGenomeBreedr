test_that("fetch_table_region() validates required arguments", {
  # Chromosome is required
  expect_error(
    fetch_table_region(con = con_test, table_name = "variants"),
    "A valid chromosome ('chrom') is required to query locus data.",
    fixed = TRUE
  )
  # Connection is required for local mode
  expect_error(
    fetch_table_region(con = NULL, table_name = "variants", chrom = "Chr05"),
    "A valid 'con' object is required for 'local' mode. Please run connect_local_db() first.",
    fixed = TRUE
  )
})

test_that("fetch_table_region() works for 'variants' table", {
  variants <- fetch_table_region(
    con = con_test,
    table_name = "variants",
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )
  expect_s3_class(variants, "data.frame")
  expect_gt(nrow(variants), 0)
  expect_true(all(variants$chrom == "Chr05"))
  expect_true(all(variants$pos >= 75104537 & variants$pos <= 75106403))
})

test_that("fetch_table_region() works for 'annotations' table", {
  annots <- fetch_table_region(
    con = con_test,
    table_name = "annotations",
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )
  expect_s3_class(annots, "data.frame")
  expect_gt(nrow(annots), 0)
  expect_true("gene_name" %in% colnames(annots))
})

test_that("fetch_table_region() for 'annotations' filters by gene_name", {
  annots <- fetch_table_region(
    con = con_test,
    table_name = "annotations",
    chrom = "Chr05",
    gene_name = "Sobic.005G213600"
  )
  expect_s3_class(annots, "data.frame")
  expect_gt(nrow(annots), 0)
  expect_true(all(annots$gene_name == "Sobic.005G213600"))
})

test_that("fetch_table_region() works for 'genotypes' table", {
  genos <- fetch_table_region(
    con = con_test,
    table_name = "genotypes",
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )
  expect_s3_class(genos, "data.frame")
  expect_gt(nrow(genos), 0)
  # Check for standard and calculated allele metric columns
  expect_true(all(c("variant_id", "chrom", "pos", "major_allele_freq", "minor_allele_freq") %in% colnames(genos)))
  # Check that sample columns are present
  expect_gt(ncol(genos), 10)
})

test_that("fetch_table_region() returns an empty data.frame for regions with no data", {
  no_data <- fetch_table_region(
    con = con_test,
    table_name = "variants",
    chrom = "Chr01",
    start = 1,
    end = 10
  )
  expect_s3_class(no_data, "data.frame")
  expect_equal(nrow(no_data), 0)
})
