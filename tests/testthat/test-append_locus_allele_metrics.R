test_that("append_locus_allele_metrics calculates frequencies correctly", {
  # Mock genotype matrix
  mock_matrix <- data.frame(
    variant_id = "VAR1",
    chrom = "Chr01",
    pos = 123,
    ref = "A",
    alt = "G",
    variant_type = "SNP",
    sample1 = "0|0", # ref/ref
    sample2 = "0|1", # het
    sample3 = "1|1", # alt/alt
    sample4 = "1|1", # alt/alt
    stringsAsFactors = FALSE
  )

  result <- append_locus_allele_metrics(mock_matrix, meta_data = c("variant_id", "chrom", "pos", "ref", "alt", "variant_type"))

  # Expected calculations:
  # count_0 (ref 'A') = 2 (s1) + 1 (s2) = 3
  # count_1 (alt 'G') = 1 (s2) + 2 (s3) + 2 (s4) = 5
  # total = 8
  # minor_freq = 3/8 = 0.375
  # major_freq = 5/8 = 0.625
  # minor_allele = 'A' (ref)
  # major_allele = 'G' (alt)

  expect_equal(result$minor_allele_freq, 0.375)
  expect_equal(result$major_allele_freq, 0.625)
  expect_equal(result$minor_allele, "A")
  expect_equal(result$major_allele, "G")
})

test_that("append_locus_allele_metrics handles column ordering", {
  mock_matrix <- data.frame(
    variant_id = "VAR1",
    chrom = "Chr01",
    pos = 123,
    ref = "A",
    alt = "G",
    variant_type = "SNP",
    sample1 = "0|0",
    sample2 = "0|1",
    stringsAsFactors = FALSE
  )

  result <- append_locus_allele_metrics(mock_matrix, meta_data = c("variant_id", "chrom", "pos", "ref", "alt", "variant_type"))

  expected_order <- c(
    "variant_id", "chrom", "pos", "ref", "alt", "variant_type",
    "major_allele", "minor_allele", "major_allele_freq", "minor_allele_freq",
    "sample1", "sample2"
  )
  expect_equal(names(result), expected_order)
})

test_that("append_locus_allele_metrics handles monomorphic cases", {
  # All reference
  mock_matrix_ref <- data.frame(
    variant_id = "VAR1", chrom = "Chr01", pos = 123, ref = "A", alt = "G", variant_type = "SNP",
    sample1 = "0|0", sample2 = "0|0",
    stringsAsFactors = FALSE
  )
  result_ref <- append_locus_allele_metrics(mock_matrix_ref, meta_data = c("variant_id", "chrom", "pos", "ref", "alt", "variant_type"))
  expect_equal(result_ref$minor_allele_freq, 0)
  expect_equal(result_ref$major_allele_freq, 1)
  expect_equal(result_ref$major_allele, "A")

  # All alternate
  mock_matrix_alt <- data.frame(
    variant_id = "VAR1", chrom = "Chr01", pos = 123, ref = "A", alt = "G", variant_type = "SNP",
    sample1 = "1|1", sample2 = "1|1",
    stringsAsFactors = FALSE
  )
  result_alt <- append_locus_allele_metrics(mock_matrix_alt, meta_data = c("variant_id", "chrom", "pos", "ref", "alt", "variant_type"))
  expect_equal(result_alt$minor_allele_freq, 0)
  expect_equal(result_alt$major_allele_freq, 1)
  expect_equal(result_alt$major_allele, "G")
})