test_that("foreground_select() correctly identifies lines with favorable alleles", {
  # Example genotype data
  geno <- data.frame(
    SNP1 = c("A:A", "A:G", "G:G", "A:A"),
    SNP2 = c("C:C", "C:T", "T:T", "C:T"),
    SNP3 = c("G:G", "G:G", "A:G", "A:A"),
    row.names = c("Line1", "Line2", "Line3", "Line4"),
    stringsAsFactors = FALSE)

  # Meta info for trait-predictive markers
  marker_info <- data.frame(
    qtl_markers = paste0("SNP", 1:3),
    fav_alleles = c("A", "C", "G"),
    alt_alleles = c("G", "T", "A"),
    stringsAsFactors = FALSE)

  # Test homozygous selection
  homo_res <- foreground_select(
    geno_data = geno,
    fore_marker_info = marker_info,
    fore_marker_col = "qtl_markers",
    fav_allele_col = "fav_alleles",
    alt_allele_col = "alt_alleles",
    select_type = "homo")

  expect_true(is.data.frame(homo_res))
  expect_equal(dim(homo_res), c(4, 3))
  expect_equal(homo_res["Line1", "SNP1"], 1)  # A:A
  expect_equal(homo_res["Line2", "SNP1"], 0)  # A:G
  expect_equal(homo_res["Line4", "SNP3"], 0)  # A:A is not G:G

  # Test heterozygous selection
  hetero_res <- foreground_select(
    geno_data = geno,
    fore_marker_info = marker_info,
    fore_marker_col = "qtl_markers",
    fav_allele_col = "fav_alleles",
    alt_allele_col = "alt_alleles",
    select_type = "hetero")

  expect_equal(hetero_res["Line2", "SNP1"], 1)  # A:G
  expect_equal(hetero_res["Line3", "SNP3"], 1)  # A:G

  # Test both selection
  both_res <- foreground_select(
    geno_data = geno,
    fore_marker_info = marker_info,
    fore_marker_col = "qtl_markers",
    fav_allele_col = "fav_alleles",
    alt_allele_col = "alt_alleles",
    select_type = "both")

  expect_equal(both_res["Line1", "SNP1"], 1)  # A:A
  expect_equal(both_res["Line2", "SNP1"], 1)  # A:G
  expect_equal(both_res["Line3", "SNP1"], 0)  # G:G

  # Check for missing marker column
  bad_marker_info <- marker_info
  bad_marker_info$qtl_markers[1] <- "MISSING_SNP"

  expect_error(
    foreground_select(
      geno_data = geno,
      fore_marker_info = bad_marker_info,
      fore_marker_col = "qtl_markers",
      fav_allele_col = "fav_alleles",
      alt_allele_col = "alt_alleles",
      select_type = "homo"),
    "not found in marker genotype data"
  )
})
