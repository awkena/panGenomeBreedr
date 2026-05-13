test_that("pgsql_calc_af computes accurate frequencies for mixed genotypes", {
  # A Genotype sample"
  test_gt <- data.frame(
    variant_id = c("v1", "v2", "v3"),
    chrom = c("Chr03", "Chr03", "Chr03"),
    pos = c(100, 200, 300),
    SampleA = c("0/0", "./.", "0|0"),
    SampleB = c("0/1", "0/1", "0|1"),
    SampleC = c("1/1", "1/1", "1|1"),
    stringsAsFactors = FALSE
  )

  # Execute
  res <- pgsql_calc_af(test_gt, chrom_col = "chrom", pos_col = "pos")

  # 2. Assertions
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)

  expect_equal(res$alt_af[1], 0.5)
  expect_equal(res$ref_af[1], 0.5)
 
  expect_equal(res$alt_af[2], 0.75)

  expect_equal(res$alt_af[3], 0.5)

  #  Metadata Preservation
  expect_true(all(c("variant_id", "chrom", "pos") %in% colnames(res)))
})

test_that("pgsql_calc_af handles 100% missing data gracefully", {
  missing_gt <- data.frame(
    variant_id = "v_all_na",
    SampleA = "./.",
    SampleB = "./.",
    stringsAsFactors = FALSE
  )

  res <- pgsql_calc_af(missing_gt)

  # Should return NA rather than NaN or crashing
  expect_true(is.na(res$alt_af[1]))
  expect_true(is.na(res$ref_af[1]))
})