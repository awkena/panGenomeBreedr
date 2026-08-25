mock_pcil_data <- list(
  pcil_metadata = data.frame(
    sample_id = c("LineA", "LineB", "LineC"),
    family = c("Fam1", "Fam1", "Fam2"),
    clan = c("Clan1", "Clan1", "Clan2"),
    stringsAsFactors = FALSE
  ),
  pcil_gene_regions = data.frame(
    GeneID = "TestGene",
    seqid = "Chr03",
    start = 1000,
    end = 2000,
    stringsAsFactors = FALSE
  ),
  pcil_introgressions = data.frame(
    SampleID = c("LineA", "LineB", "LineC"),
    ChrLabel = c("Chr03", "Chr03", "Chr03"),
    block_start_bp = c(1500, 5000, 10000), # LineA overlaps the gene/position!
    block_end_bp = c(1600, 6000, 12000),
    block_len_Mb = c(0.001, 0.001, 0.002),
    mean_donor_frac = c(0.9, 0.8, 0.95),
    stringsAsFactors = FALSE
  ),
  pcil_genomewide_introgressions = data.frame(
    SampleID = c("LineA", "LineB", "LineC"),
    total_blocks = c(10, 15, 5),
    total_Mb = c(20.5, 30.0, 10.0),
    stringsAsFactors = FALSE
  ),
  pcil_inbreeding_coefficient = data.frame(
    sample_id = c("LineA", "LineB", "LineC"),
    F = c(0.85, 0.75, 0.90),
    stringsAsFactors = FALSE
  ),
  pcil_IBS_dis = data.frame(
    sample_1 = c("LineA", "LineA", "LineB"),
    sample_2 = c("LineB", "LineC", "LineC"),
    distance = c(0.1, 0.4, 0.2), # LineA is closely related to LineB (dist = 0.1)
    stringsAsFactors = FALSE
  )
)

mock_variants_geno <- data.frame(
  variant_id = "SNP_Test",
  chrom = "Chr03",
  pos = 1550,
  stringsAsFactors = FALSE
)


# ==============================================================================
# TEST fetch_pcil_positive()
# ==============================================================================
test_that("fetch_pcil_positive filters correct overlapping lines", {
  res <- fetch_pcil_positive(
    pcil_data = mock_pcil_data,
    variants_select_geno = mock_variants_geno,
    type = "position",
    sel = 5,
    window = 0
  )

  # Check 1: Did it return the expected list structure?
  expect_type(res, "list")
  expect_named(res, c("pcil_positive", "best_lines", "regions"))

  # Check 2: Did it correctly identify LineA as the positive carrier?
  # (Because pos 1550 falls inside LineA's block_start 1500 and block_end 1600)
  expect_equal(nrow(res$best_lines), 1)
  expect_equal(res$best_lines$SampleID[1], "LineA")

  # Check 3: Did it merge the inbreeding coefficient (F) properly?
  expect_true("F" %in% colnames(res$best_lines))
  expect_equal(res$best_lines$F[1], 0.85)
})


# ==============================================================================
#  TEST fetch_pcil_negative()
# ==============================================================================
test_that("fetch_pcil_negative identifies genetically similar non-carriers", {
  # Create a mock positive result pretending LineA is our focal positive line
  mock_pos_result <- list(
    best_lines = data.frame(
      SampleID = "LineA",
      Region = "SNP_Test",
      Chr = "Chr03",
      Start = 1550,
      End = 1550,
      stringsAsFactors = FALSE
    ),
    regions = data.frame(
      Region = "SNP_Test",
      Chr = "Chr03",
      Start = 1550,
      End = 1550,
      stringsAsFactors = FALSE
    )
  )

  # Run the negative selection function
  res_neg <- fetch_pcil_negative(
    pcil_data = mock_pcil_data,
    pcil_positive_result = mock_pos_result,
    n_neg = 1
  )

  # Check 1: Output structure
  expect_type(res_neg, "list")
  expect_true("pairs_best" %in% names(res_neg))

  expect_equal(nrow(res_neg$pairs_best), 1)
  expect_equal(res_neg$pairs_best$SampleID_Negative[1], "LineB")
  expect_equal(res_neg$pairs_best$IBS_dis[1], 0.1)
  expect_equal(res_neg$pairs_best$Family[1], "Fam1")
})


# ==============================================================================
# TEST PLOTTING FUNCTIONS 
# ==============================================================================
test_that("Plotting functions generate ggplot objects without error", {
  # Mock positive result
  mock_pos_result <- list(
    pcil_positive = data.frame(
      SampleID = "LineA",
      Region = "SNP_Test",
      Chr = "Chr03",
      Start = 1550,
      End = 1550,
      ChrLabel = "Chr03",
      block_start_bp = 1500,
      block_end_bp = 1600,
      stringsAsFactors = FALSE
    ),
    best_lines = data.frame(
      SampleID = "LineA",
      Region = "SNP_Test",
      Chr = "Chr03",
      Start = 1550,
      End = 1550,
      Rank = 1,
      stringsAsFactors = FALSE
    ),
    regions = data.frame(
      Region = "SNP_Test",
      Chr = "Chr03",
      Start = 1550,
      End = 1550,
      stringsAsFactors = FALSE
    )
  )

  # Mock pcil_introgressions with a plot_id for the best_lines logic
  mock_data_for_plot <- mock_pcil_data
  mock_data_for_plot$pcil_introgressions$Clan <- c("Clan1", "Clan1", "Clan2")
  mock_data_for_plot$pcil_introgressions$Family <- c("Fam1", "Fam1", "Fam2")
  mock_data_for_plot$pcil_introgressions$startMb <- mock_data_for_plot$pcil_introgressions$block_start_bp /
    1e6
  mock_data_for_plot$pcil_introgressions$endMb <- mock_data_for_plot$pcil_introgressions$block_end_bp /
    1e6

  p1 <- plot_pcil_positive(mock_pos_result)
  p2 <- plot_pcil_best_lines(mock_pos_result, mock_data_for_plot)

  # Validate they output a list of plots that inherit the ggplot class
  expect_type(p1, "list")
  expect_s3_class(p1[[1]], "ggplot")

  expect_type(p2, "list")
  expect_s3_class(p2[[1]], "ggplot")
})