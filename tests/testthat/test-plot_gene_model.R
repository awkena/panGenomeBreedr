
# --- Helper: Create mock gene data frames ---
mock_gene_pos <- data.frame(
  ID = c(
    "Sobic.005G213600",
    "Sobic.005G213600.1",
    "Sobic.005G213600.1",
    "Sobic.005G213600.1",
    "Sobic.005G213600.1"
  ),
  Feature = c("gene", "mRNA", "five_prime_UTR", "CDS", "three_prime_UTR"),
  Chromosome = "Chr05",
  Start = c(75104500, 75104500, 75104500, 75104601, 75106001),
  End = c(75106500, 75106500, 75104600, 75106000, 75106500),
  Strand = "+",
  stringsAsFactors = FALSE
)

# Simulates a negative strand gene to test the arrow direction logic
mock_gene_neg <- data.frame(
  ID = c("Sobic.DUP", "Sobic.DUP.1", "Sobic.DUP.1"),
  Feature = c("gene", "mRNA", "CDS"),
  Chromosome = "Chr05",
  Start = c(1000, 1000, 1200),
  End = c(2000, 2000, 1800),
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


test_that("plot_gene_model drops the gene summary row when a transcript is present", {

  p <- plot_gene_model(mock_gene_pos)
  layer_data <- p$layers[[1]]$data

  # The gene-level backbone should not be drawn as its own row...
  expect_equal(nrow(layer_data[layer_data$Feature == "gene", ]), 0)
  # ...only the transcript's own ID should appear as a row
  expect_equal(unique(layer_data$ID), "Sobic.005G213600.1")

})


test_that("plot_gene_model falls back to the gene row when no transcripts exist", {

  # A gene with no annotated mRNA/transcript children at all
  gene_only <- data.frame(
    ID = "Sobic.LONELY",
    Feature = "gene",
    Chromosome = "Chr01",
    Start = 100,
    End = 500,
    Strand = "+",
    stringsAsFactors = FALSE
  )

  expect_silent(p <- plot_gene_model(gene_only))

  layer_data <- p$layers[[1]]$data
  expect_equal(nrow(layer_data), 1)
  expect_equal(layer_data$Feature, "gene")

})


test_that("plot_gene_model accurately calculates structural y-heights based on feature type", {

  # row_height = 1 pins the base ratios (0.15 / 0.075 / 0.01) so the expected
  # values below are exact, independent of whatever row_height defaults to
  p <- plot_gene_model(mock_gene_pos, row_height = 1)

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

  # Check transcript backbone (mRNA) heights, on its own row
  mrna_data <- layer_data[layer_data$Feature == "mRNA", ]
  expect_equal(mrna_data$ymin, -0.01)
  expect_equal(mrna_data$ymax, 0.01)

})


test_that("plot_gene_model scales bar thickness proportionally with row_height", {

  p_full <- plot_gene_model(mock_gene_pos, row_height = 1)
  p_half <- plot_gene_model(mock_gene_pos, row_height = 0.5)

  cds_full <- p_full$layers[[1]]$data
  cds_half <- p_half$layers[[1]]$data
  cds_full <- cds_full[cds_full$Feature == "CDS", ]
  cds_half <- cds_half[cds_half$Feature == "CDS", ]

  # Halving row_height should exactly halve bar thickness, so proportions
  # (and thus visual chunkiness) stay constant as spacing is tuned
  expect_equal(cds_half$ymax, cds_full$ymax * 0.5)
  expect_equal(cds_half$ymin, cds_full$ymin * 0.5)

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

  # Test with just the transcript backbone and CDS, no UTRs
  minimal_gene <- mock_gene_pos[mock_gene_pos$Feature %in% c("mRNA", "CDS"), ]

  # Function should not fail
  expect_silent(p <- plot_gene_model(minimal_gene))
  expect_s3_class(p, "ggplot")

})


test_that("plot_gene_model stacks multiple transcripts on separate rows", {

  mock_multi <- rbind(
    mock_gene_pos,
    data.frame(
      ID = c("Sobic.005G213600.2", "Sobic.005G213600.2"),
      Feature = c("mRNA", "CDS"),
      Chromosome = "Chr05",
      Start = c(75104550, 75104650),
      End = c(75106500, 75106000),
      Strand = "+",
      stringsAsFactors = FALSE
    )
  )

  p <- plot_gene_model(mock_multi)
  layer_data <- p$layers[[1]]$data

  # Two distinct transcripts should end up on two distinct rows
  rows_by_id <- unique(layer_data[, c("ID", "y_row")])
  expect_equal(nrow(rows_by_id), 2)
  expect_equal(length(unique(rows_by_id$y_row)), 2)

})
