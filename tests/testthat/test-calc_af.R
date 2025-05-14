test_that("calc_af() computes allele frequencies correctly from VCF-style genotype strings", {
  # Mock genotype data
  gt_df <- data.frame(
    variant_id = c("v1", "v2", "v3", "v4"),
    chrom = c("Chr01", "Chr01", "Chr02", "Chr02"),
    pos = c(100, 200, 300, 400),
    s1 = c("0/0", "0/1", "1/1", NA),
    s2 = c("0/1", "1/1", "0/0", "0/1"),
    s3 = c("1/1", "0/0", "0/1", "1/1"),
    stringsAsFactors = FALSE)

  # Run calc_af() with chrom and pos included
  af_df <- calc_af(gt_df,
                   variant_id_col = "variant_id",
                   chrom_col = "chrom",
                   pos_col = "pos")

  # Check output structure
  expect_true(is.data.frame(af_df))
  expect_named(af_df, c("variant_id", "chrom", "pos", "ref_af", "alt_af"))
  expect_equal(nrow(af_df), 4)

  # Check known values
  # For v1: genotypes = 0, 1, 2 → total alt = (0+1+2) = 3 → AF = 3 / (2*3) = 0.5
  expect_equal(round(af_df$alt_af[af_df$variant_id == "v1"], 2), 0.5)
  expect_equal(round(af_df$ref_af[af_df$variant_id == "v1"], 2), 0.5)

  # For v4: genotypes = NA, 1, 2 → alt = 3, n = 2 → AF = 3 / 4 = 0.75
  expect_equal(round(af_df$alt_af[af_df$variant_id == "v4"], 2), 0.75)

  # Run again without chrom/pos columns
  af_simple <- calc_af(gt_df[, c("variant_id", "s1", "s2", "s3")],
                       variant_id_col = "variant_id")
  expect_named(af_simple, c("variant_id", "ref_af", "alt_af"))

  # Confirm rows match
  expect_equal(nrow(af_simple), 4)
})
