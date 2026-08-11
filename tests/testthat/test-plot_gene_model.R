
# --- Helper: Create mock gene data frames ---
# Simulates the output of gene_coord_gff() for a positive strand gene
mock_gene_pos <- data.frame(
  ID = "Sobic.005G213600",
  Feature = c("gene", "five_prime_UTR", "CDS", "three_prime_UTR"),
  Chromosome = "Chr05",
  Start = c(75104500, 75104500, 75104601, 75106001),
  End = c(75106500, 75104600, 75106000, 75106500),
  Strand = "+",
  stringsAsFactors = FALSE
)

# Simulates a negative strand gene to test the arrow direction logic
mock_gene_neg <- data.frame(
  ID = "Sobic.DUP",
  Feature = c("gene", "CDS"),
  Chromosome = "Chr05",
  Start = c(1000, 1200),
  End = c(2000, 1800),
  Strand = "-",
  stringsAsFactors = FALSE
)


# --- Tests ---

test_that("plot_gene_model returns a valid ggplot object with expected layers", {

  p <- plot_gene_model(mock_gene_pos)

  # 1. Check if the output is successfully built as a ggplot object
  expect_s3_class(p, c("gg", "ggplot"))

  # 2. Check that the correct geometries were added to the plot
  # Layer 1 should be the geom_rect (features)
  expect_true(inherits(p$layers[[1]]$geom, "GeomRect"))
  # Layer 2 should be the geom_segment (backbone/arrow)
  expect_true(inherits(p$layers[[2]]$geom, "GeomSegment"))

})


test_that("plot_gene_model accurately calculates structural y-heights based on feature type", {

  p <- plot_gene_model(mock_gene_pos)

  # Extract the data specifically passed to the geom_rect layer
  # (Since ggplot() was called empty, p$data is empty, data lives in the layer)
  layer_data <- p$layers[[1]]$data

  # Check CDS heights
  cds_data <- layer_data[layer_data$Feature == "CDS", ]
  expect_equal(cds_data$ymin, -0.15)
  expect_equal(cds_data$ymax, 0.15)

  # Check UTR heights
  utr_data <- layer_data[grepl("UTR", layer_data$Feature), ]
  expect_true(all(utr_data$ymin == -0.075))
  expect_true(all(utr_data$ymax == 0.075))

  # Check backbone (gene) heights
  gene_data <- layer_data[layer_data$Feature == "gene", ]
  expect_equal(gene_data$ymin, -0.01)
  expect_equal(gene_data$ymax, 0.01)

})


test_that("plot_gene_model maps the arrow direction correctly based on strand (+)", {

  p_pos <- plot_gene_model(mock_gene_pos)

  # Extract the data specifically passed to the geom_segment layer
  backbone_data <- p_pos$layers[[2]]$data

  # For a positive strand, x_start should be Start, and x_end should be End (pointing right)
  expect_equal(backbone_data$x_start, backbone_data$Start)
  expect_equal(backbone_data$x_end, backbone_data$End)

})


test_that("plot_gene_model maps the arrow direction correctly based on strand (-)", {

  p_neg <- plot_gene_model(mock_gene_neg)

  # Extract the data specifically passed to the geom_segment layer
  backbone_data <- p_neg$layers[[2]]$data

  # For a negative strand, x_start should be End, and x_end should be Start (pointing left)
  expect_equal(backbone_data$x_start, backbone_data$End)
  expect_equal(backbone_data$x_end, backbone_data$Start)

})


test_that("plot_gene_model handles data frames with missing features gracefully", {

  # Test with just a gene backbone and CDS, no UTRs
  minimal_gene <- mock_gene_pos[mock_gene_pos$Feature %in% c("gene", "CDS"), ]

  # Function should not fail
  expect_silent(p <- plot_gene_model(minimal_gene))
  expect_s3_class(p, "ggplot")

})
