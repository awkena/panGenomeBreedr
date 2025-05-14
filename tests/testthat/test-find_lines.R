test_that("find_lines() identifies lines based on presence and absence of loci", {
  # Mock binary matrix from foreground_select
  mat <- data.frame(
    SNP1 = c(1, 0, 0, 1),
    SNP2 = c(1, 0, 0, 0),
    SNP3 = c(1, 1, 0, 0),
    row.names = c("Line1", "Line2", "Line3", "Line4"))

  # Lines with all loci present
  result_all_present <- find_lines(mat, present = c("SNP1", "SNP2", "SNP3"))
  expect_equal(result_all_present, "Line1")

  # Lines with SNP3 absent (0)
  result_absent <- find_lines(mat, absent = c("SNP3"))
  expect_equal(result_absent, c("Line3", "Line4"))

  # Combination: SNP1 and SNP2 present, SNP3 absent
  result_combined <- find_lines(mat, present = c("SNP1", "SNP2"), absent = c("SNP3"))
  expect_equal(result_combined, character())

  # No filter: should return all
  result_all <- find_lines(mat)
  expect_equal(result_all, rownames(mat))

  # Case with no matching lines
  result_none <- find_lines(mat, present = c("SNP1", "SNP2", "SNP3"), absent = c("SNP2"))
  expect_length(result_none, 0)
  expect_type(result_none, "character")
})
