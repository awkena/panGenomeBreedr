

# --- Helper: Create mock variant data ---
mock_variants <- data.frame(
  chrom.x = rep("Chr05", 5),
  pos.x = c(1000, 2000, 3000, 4000, 5000),
  variant_id = c("v1", "v2", "v3", "v4", "v5"),
  variant_type = c("SNP", "SNP", "INDEL", "SNP", "INDEL"),
  impact = c("HIGH", "MODERATE", "LOW", "MODIFIER", "HIGH"),
  stringsAsFactors = FALSE
)

# Mock data with custom column names to test tidy evaluation
mock_custom_cols <- data.frame(
  my_chr = rep("Chr01", 3),
  my_pos = c(150, 250, 350),
  variant_id = c("vA", "vB", "vC"),
  variant_type = c("SNP", "INDEL", "SNP"),
  impact = c("LOW", "LOW", "MODERATE"),
  stringsAsFactors = FALSE
)

# --- Tests ---

test_that("plot_variant_hotspot returns a valid ggplot object with expected layers", {

  plt <- plot_variant_hotspot(mock_variants)

  # Check if the output is successfully built as a ggplot object
  expect_s3_class(plt, c("gg", "ggplot"))

  # Check that the correct geometries were added to the plot in the right order
  expect_true(inherits(plt$layers[[1]]$geom, "GeomHline"))
  expect_true(inherits(plt$layers[[2]]$geom, "GeomSegment"))
  expect_true(inherits(plt$layers[[3]]$geom, "GeomPoint"))

})

test_that("plot_variant_hotspot correctly stratifies y-coordinates based on Impact", {

  plt <- plot_variant_hotspot(mock_variants)

  # The modified data frame is stored directly in the plot object
  plot_data <- plt$data

  # Verify base_y mapping exactly matches the function's internal logic
  expect_equal(plot_data$base_y[plot_data$impact == "HIGH"], c(4, 4))
  expect_equal(plot_data$base_y[plot_data$impact == "MODERATE"], 3)
  expect_equal(plot_data$base_y[plot_data$impact == "LOW"], 2)
  expect_equal(plot_data$base_y[plot_data$impact == "MODIFIER"], 1)

  # Check that plot_y contains the jitter (should not be exactly equal to base_y)
  expect_true(all(plot_data$plot_y != plot_data$base_y))
  expect_true(all(abs(plot_data$plot_y - plot_data$base_y) <= 0.3))

})

test_that("plot_variant_hotspot correctly calculates summary statistics for the title", {

  plt <- plot_variant_hotspot(mock_variants)

  # 5 variants total: 3 SNPs, 2 INDELs
  expected_title <- "Total Unique Variants: 5 | SNPs: 3 | INDELS: 2"

  # Check the generated title in the plot labels
  expect_equal(plt$labels$title, expected_title)

  # Check subtitle for correct chromosome and comma-formatted boundaries
  expected_subtitle <- "Chr05: 1,000 - 5,000"
  expect_equal(plt$labels$subtitle, expected_subtitle)

})

test_that("plot_variant_hotspot successfully handles custom column names via tidy evaluation", {

  # Function should not fail when using non-standard column names
  expect_silent(
    plt <- plot_variant_hotspot(var_df = mock_custom_cols,
                                chrom_col = "my_chr",
                                pos_col = "my_pos")
  )

  expect_s3_class(plt, "ggplot")

  # Verify the subtitle picked up the custom chromosome name
  expected_subtitle <- "Chr01: 150 - 350"
  expect_equal(plt$labels$subtitle, expected_subtitle)

})

test_that("plot_variant_hotspot respects manual boundary overrides", {

  plt <- plot_variant_hotspot(mock_variants, start_pos = 0, end_pos = 10000)

  # Check subtitle reflects manual boundaries
  expected_subtitle <- "Chr05: 0 - 10,000"
  expect_equal(plt$labels$subtitle, expected_subtitle)

  # Check the underlying x-axis scale limits
  scale_x <- plt$scales$get_scales("x")
  expect_equal(scale_x$limits, c(0, 10000))

})
