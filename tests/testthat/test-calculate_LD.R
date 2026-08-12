# mock genotype dataframe for testing
mock_geno_df <- data.frame(
  variant_id = c("SNP_1", "SNP_2", "INDEL_3", "SNP_4_mono"),
  pos = c(100, 200, 300, 400),
  sample1 = c("0|0", "0/1", "1|1", "0|0"),
  sample2 = c("0/1", "1/1", "0|0", "0|0"),
  sample3 = c("1|1", "0|0", "0/1", "0|0"),
  stringsAsFactors = FALSE
)

test_that("calculate_LD works correctly in All-vs-All mode", {

  ld_results <- calculate_LD(
    df = mock_geno_df,
    target_variant_ids = NULL,
    genotype_start_col = 3
  )

  # Check output structure
  expect_s3_class(ld_results, "data.frame")
  expected_cols <- c(
    "variant_1", "position_1", "variant_type_1",
    "variant_2", "position_2", "variant_type_2",
    "distance_bp", "R2", "D_prime"
  )
  expect_equal(colnames(ld_results), expected_cols)

  # Check number of pairs (3 polymorphic variants: (3*2)/2 = 3 pairs)
  expect_equal(nrow(ld_results), 3)

  # Monomorphic variants should be excluded
  expect_false("SNP_4_mono" %in% ld_results$variant_1)
  expect_false("SNP_4_mono" %in% ld_results$variant_2)

  # Check variant types are correctly identified
  expect_true("SNP" %in% c(ld_results$variant_type_1, ld_results$variant_type_2))
  expect_true("INDEL" %in% c(ld_results$variant_type_1, ld_results$variant_type_2))

  # Check a specific calculation (SNP_1 vs SNP_2)
  snp1_vs_snp2 <- ld_results[ld_results$variant_1 == "SNP_1" & ld_results$variant_2 == "SNP_2", ]
  expect_equal(snp1_vs_snp2$R2, round(1/9, 5))
  expect_equal(snp1_vs_snp2$D_prime, round(1/3, 5))
  expect_equal(snp1_vs_snp2$distance_bp, 100)
})

test_that("calculate_LD works correctly in Targeted mode", {

  target <- "SNP_1"
  ld_results <- calculate_LD(
    df = mock_geno_df,
    target_variant_ids = target,
    genotype_start_col = 3
  )

  # Check output structure
  expect_s3_class(ld_results, "data.frame")

  # Should compute LD against 2 other polymorphic variants
  expect_equal(nrow(ld_results), 2)

  # All `variant_1` should be the target
  expect_true(all(ld_results$variant_1 == target))

  # Target should not be in `variant_2`
  expect_false(target %in% ld_results$variant_2)
})

test_that("calculate_LD handles multiple targets correctly", {
  targets <- c("SNP_1", "INDEL_3")
  ld_results <- calculate_LD(
    df = mock_geno_df,
    target_variant_ids = targets,
    genotype_start_col = 3
  )

  # Each of the 2 targets is compared with 2 other polymorphic variants
  expect_equal(nrow(ld_results), 4)
  expect_true(all(unique(ld_results$variant_1) %in% targets))
})


test_that("calculate_LD handles errors and warnings", {

  # Error on missing position column
  df_no_pos <- mock_geno_df
  colnames(df_no_pos)[2] <- "location"
  expect_error(
    calculate_LD(df_no_pos, genotype_start_col = 3),
    "Could not find a physical position column"
  )

  # Error on missing target variant
  expect_error(
    calculate_LD(mock_geno_df, target_variant_ids = "SNP_DOES_NOT_EXIST", genotype_start_col = 3),
    "The following anchor variant ID\\(s\\) were not found in the dataframe: SNP_DOES_NOT_EXIST"
  )

  # Returns empty data.frame and warns when only monomorphic target is provided
  expect_warning(
    result <- calculate_LD(mock_geno_df, target_variant_ids = "SNP_4_mono", genotype_start_col = 3),
    "Anchor variant SNP_4_mono is monomorphic. Skipping."
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), 9)
})
