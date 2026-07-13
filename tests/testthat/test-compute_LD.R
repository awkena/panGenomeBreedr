test_that("compute_LD works in all-vs-all mode", {
  # Create mock genotype data
  mock_geno_df <- data.frame(
    variant_id = c("SNP_1", "SNP_2", "SNP_3", "SNP_4_mono"),
    chrom = "Chr01",
    pos = c(100, 200, 300, 400),
    ref = "A",
    alt = "T",
    sample1 = c("0|0", "0|1", "1|1", "0|0"),
    sample2 = c("0|1", "1|1", "0|0", "0|0"),
    sample3 = c("1|1", "0|0", "0|1", "0|0"),
    stringsAsFactors = FALSE
  )

  ld_results <- compute_LD(df = mock_geno_df, genotype_start_col = 6)

  # Check output type and columns
  expect_s3_class(ld_results, "data.frame")
  expect_named(ld_results, c("variant_1", "position_1", "variant_type_1", "variant_2", "position_2", "variant_type_2", "distance_bp", "R2", "D_prime"))

  # Check that monomorphic variant is excluded
  expect_false("SNP_4_mono" %in% ld_results$variant_1)
  expect_false("SNP_4_mono" %in% ld_results$variant_2)

  # Check number of pairs (3 polymorphic variants -> 3 pairs)
  expect_equal(nrow(ld_results), 3)

  # Check values for a specific pair (SNP_1 vs SNP_2)
  pair12 <- ld_results[ld_results$variant_1 == "SNP_1" & ld_results$variant_2 == "SNP_2", ]
  expect_equal(pair12$distance_bp, 100)
  expect_true(pair12$R2 >= 0 && pair12$R2 <= 1)
  expect_true(pair12$D_prime >= 0 && pair12$D_prime <= 1)
})

test_that("compute_LD works in targeted mode", {
  mock_geno_df <- data.frame(
    variant_id = c("SNP_1", "SNP_2", "SNP_3"),
    chrom = "Chr01",
    pos = c(100, 200, 300),
    sample1 = c("0|0", "0|1", "1|1"),
    sample2 = c("0|1", "1|1", "0|0"),
    sample3 = c("1|1", "0|0", "0|1"),
    stringsAsFactors = FALSE
  )

  target_variants <- c("SNP_1")
  ld_results <- compute_LD(df = mock_geno_df, target_variant_ids = target_variants, genotype_start_col = 4)

  # Check that only target variant is in variant_1
  expect_true(all(ld_results$variant_1 == "SNP_1"))
  expect_equal(nrow(ld_results), 2) # SNP_1 vs SNP_2 and SNP_1 vs SNP_3
})


test_that("compute_LD handles errors correctly", {
  mock_geno_df <- data.frame(
    variant_id = c("SNP_1", "SNP_2"),
    chrom = "Chr01",
    position = c(100, 200), # incorrect pos column name for the default
    sample1 = c("0|0", "0|1"),
    stringsAsFactors = FALSE
  )
  # No 'pos' column
  expect_error(compute_LD(mock_geno_df[, -3], genotype_start_col = 4),
               "Could not find a physical position column")

  # Invalid target variant
  expect_error(compute_LD(mock_geno_df, target_variant_ids = "SNP_99", genotype_start_col = 4),
               "The following anchor variant ID\\(s\\) were not found")

  # No polymorphic pairs
  mono_df <- data.frame(
    variant_id = c("SNP_1", "SNP_2"),
    pos = c(100, 200),
    sample1 = c("0|0", "1|1"),
    sample2 = c("0|0", "1|1"),
    stringsAsFactors = FALSE
  )
  expect_error(compute_LD(mono_df, genotype_start_col = 3),
               "No valid polymorphic pairs were successfully calculated.")
})