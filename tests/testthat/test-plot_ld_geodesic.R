test_that("plot_ld_geodesic returns correct structure", {
  # Mock data
  ld_df <- data.frame(
    variant_1 = "SNP_1",
    variant_2 = "SNP_2",
    R2 = 0.9,
    D_prime = 0.95
  )
  query_db_geno <- data.frame(
    variant_id = c("SNP_1", "SNP_2"),
    pos = c(1000, 2000)
  )
  query_db_annot <- data.frame(
    variant_id = c("SNP_1", "SNP_2"),
    impact = c("MODERATE", "LOW")
  )

  # Run function
  result <- plot_ld_geodesic(
    ld_df = ld_df,
    query_db_geno = query_db_geno,
    query_db_annot = query_db_annot,
    threshold = 0.1,
    block_threshold = 0.8
  )

  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("plot", "table"))
  expect_s3_class(result$plot, "ggplot")
  expect_s3_class(result$table, "data.frame")
})

test_that("plot_ld_geodesic identifies haploblocks correctly", {
  # Mock data
  ld_df <- data.frame(
    variant_1 = c("SNP_1", "SNP_1", "SNP_3"),
    variant_2 = c("SNP_2", "SNP_3", "SNP_4"),
    R2 = c(0.9, 0.5, 0.2) # Only SNP_1 and SNP_2 form a block
  )
  query_db_geno <- data.frame(
    variant_id = c("SNP_1", "SNP_2", "SNP_3", "SNP_4"),
    pos = c(1000, 2000, 3000, 4000)
  )
  query_db_annot <- data.frame(
    variant_id = c("SNP_1", "SNP_2", "SNP_3", "SNP_4"),
    impact = c("HIGH", "MODERATE", "LOW", "MODIFIER")
  )

  result <- plot_ld_geodesic(
    ld_df = ld_df,
    query_db_geno = query_db_geno,
    query_db_annot = query_db_annot,
    threshold = 0.1,
    block_threshold = 0.8 # R2 >= 0.8 for a block
  )

  haplo_table <- result$table
  expect_equal(ncol(haplo_table), 2)
  expect_true("Block_1" %in% colnames(haplo_table))
  expect_equal(sum(!is.na(haplo_table$Block_1)), 2)
  expect_true(all(c("SNP_1", "SNP_2") %in% haplo_table$Block_1))
  expect_true(all(c("HIGH", "MODERATE") %in% haplo_table$Block_1_Impact_Level))
})

test_that("plot_ld_geodesic handles no variants meeting threshold", {
  ld_df <- data.frame(variant_1 = "SNP_1", variant_2 = "SNP_2", R2 = 0.1)
  query_db_geno <- data.frame(variant_id = c("SNP_1", "SNP_2"), pos = c(1000, 2000))
  query_db_annot <- data.frame(variant_id = c("SNP_1", "SNP_2"), impact = "LOW")

  expect_error(
    plot_ld_geodesic(
      ld_df = ld_df,
      query_db_geno = query_db_geno,
      query_db_annot = query_db_annot,
      threshold = 0.5 # Higher than any R2 value
    ),
    "No variants meet the current threshold parameters."
  )
})

test_that("plot_ld_geodesic handles no haploblocks", {
    ld_df <- data.frame(variant_1 = "SNP_1", variant_2 = "SNP_2", R2 = 0.7)
    query_db_geno <- data.frame(variant_id = c("SNP_1", "SNP_2"), pos = c(1000, 2000))
    query_db_annot <- data.frame(variant_id = c("SNP_1", "SNP_2"), impact = "LOW")

    result <- plot_ld_geodesic(
        ld_df = ld_df,
        query_db_geno = query_db_geno,
        query_db_annot = query_db_annot,
        threshold = 0.5,
        block_threshold = 0.9 # Higher than any R2 value
    )

    expect_s3_class(result$table, "data.frame")
    expect_equal(nrow(result$table), 0)
    expect_equal(ncol(result$table), 0)
})