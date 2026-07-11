test_that("marker.chr_ID correctly parses the SNP snpeff VCF fixture", {
  # Locate the VCF file within your package infrastructure
  vcf_path <- system.file(
    "extdata",
    "Sobic.002G302700_SNP_snpeff.vcf",
    package = "panGenomeBreedr",
    mustWork = TRUE
  )

  # Access internal function using triple colon
  # Note: ensure vcfR is available in your test environment
  skip_if_not_installed("vcfR")

  # Execution
  result <- panGenomeBreedr:::marker.chr_ID(vcf_path)

  # Assertions
  expect_type(result, "list")
  expect_named(result, c("vcf_matrix_markerID", "vcf_matrix_chromID"))

  # Verify that we extracted character vectors
  expect_type(result$vcf_matrix_markerID, "character")
  expect_type(result$vcf_matrix_chromID, "character")

  # Logic check: ensure at least one chromosome was identified
  expect_gt(length(result$vcf_matrix_chromID), 0)
})


test_that("plates_col adds a plates column only when missing", {
  # Scenario A: Column is missing - should create it
  df_missing <- data.frame(
    MasterPlate = c("P1", "P1"),
    SNPID = c("SNP01", "SNP02"),
    stringsAsFactors = FALSE
  )

  result_created <- panGenomeBreedr:::plates_col(df_missing)

  expect_true("plates" %in% colnames(result_created))
  expect_equal(result_created$plates, c("P1_SNP01", "P1_SNP02"))

  # Scenario B: Column already exists - should remain untouched
  df_existing <- data.frame(
    MasterPlate = c("P1"),
    SNPID = c("SNP01"),
    plates = c("PRE-EXISTING"),
    stringsAsFactors = FALSE
  )

  result_preserved <- panGenomeBreedr:::plates_col(df_existing)
  expect_equal(result_preserved$plates, "PRE-EXISTING")
})


test_that("plates_col stops when critical columns are missing", {
  # Create a malformed data frame missing 'SNPID'
  df_bad <- data.frame(
    MasterPlate = c("P1"),
    OtherCol = c("Data")
  )

  expect_error(
    panGenomeBreedr:::plates_col(df_bad),
    "Error: Columns 'MasterPlate' and 'SNPID' are missing"
  )
})


test_that("uniq_plates correctly extracts unique identifiers from KASP data", {
  # Scenario: Data needs the plates column to be created
  df_input <- data.frame(
    MasterPlate = c("P1", "P1", "P2"),
    SNPID = c("S1", "S2", "S1"),
    stringsAsFactors = FALSE
  )

  # Access internal function
  result <- panGenomeBreedr:::uniq_plates(df_input)

  # Assertions
  expect_type(result, "character")
  expect_length(result, 3) # "P1_S1", "P1_S2", "P2_S1"
  expect_true(all(c("P1_S1", "P1_S2", "P2_S1") %in% result))
})


test_that("uniq_plates handles pre-existing plates column", {
  # Scenario: plates column already exists, should ignore MasterPlate/SNPID concatenation
  df_input <- data.frame(
    MasterPlate = c("P1", "P2"),
    SNPID = c("S1", "S1"),
    plates = c("CUSTOM_PLATE_A", "CUSTOM_PLATE_B"),
    stringsAsFactors = FALSE
  )

  result <- panGenomeBreedr:::uniq_plates(df_input)

  expect_equal(result, c("CUSTOM_PLATE_A", "CUSTOM_PLATE_B"))
})


test_that("uniq_plates throws an error if plates column cannot be generated", {
  # Scenario: Missing mandatory columns (SNPID) so plates_col() should fail
  df_bad <- data.frame(MasterPlate = c("P1"))

  expect_error(
    panGenomeBreedr:::uniq_plates(df_bad),
    "Error: Columns 'MasterPlate' and 'SNPID' are missing"
  )
})


test_that("kasp_marker_map correctly extracts SNP details from SNPID strings", {
  # Create a mock data frame with various ID formats
  df_kasp <- data.frame(
    SNPID = c("Chr01_100", "Chr02_500", "Chr01_100"), # Chr01_100 repeated
    stringsAsFactors = FALSE
  )

  # Access internal function
  result <- panGenomeBreedr:::kasp_marker_map(df_kasp)

  # Assertions
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2) # Should be unique SNPIDs (Chr01_100, Chr02_500)

  # Verify data types and structure
  expect_true(all(c("SNPID", "SNP_ID2", "chr", "pos") %in% colnames(result)))
  expect_type(result$pos, "integer")

  # Check one specific extraction logic
  expect_equal(result$pos[1], 100L)
  expect_equal(result$SNP_ID2[1], "Chr|01")
})

test_that("kasp_marker_map throws errors for invalid input", {
  # Missing SNPID column
  expect_error(
    panGenomeBreedr:::kasp_marker_map(data.frame(OtherCol = 1)),
    "Kasp_data must have column: SNPID"
  )

  # Not a data frame
  expect_error(
    panGenomeBreedr:::kasp_marker_map("not a dataframe"),
    "Kasp_data must be a data frame"
  )
})


test_that("get_calls extracts specific genotype calls for a given plate", {
  # Mock data with a 'plates' and 'Call' column
  df_kasp <- data.frame(
    MasterPlate = c("P1", "P1", "P2"),
    SNPID = c("S1", "S2", "S1"),
    Call = c("AA", "AB", "BB"),
    stringsAsFactors = FALSE
  )

  # Access internal function
  result <- panGenomeBreedr:::get_calls(df_kasp, "P1_S1")

  # Check extraction
  expect_equal(result, "AA")

  # Check NULL return for non-existent plate
  null_result <- panGenomeBreedr:::get_calls(df_kasp, "NON_EXISTENT")
  expect_null(null_result)
})


test_that("col_names.data handles column name extraction and empty inputs", {
  # Test with a standard data frame
  df_data <- data.frame(SNPID = 1, Call = 2)
  expect_equal(panGenomeBreedr:::col_names.data(df_data), c("SNPID", "Call"))

  # Test with empty data frame
  expect_equal(panGenomeBreedr:::col_names.data(data.frame()), character(0))

  # Test error handling for invalid input
  expect_error(
    panGenomeBreedr:::col_names.data("not_a_df"),
    "Input must be a data frame or matrix"
  )
})


test_that("genotypes converts list to correctly named data frame", {
  # Test standard input
  input_list <- list(S1 = "AA", S2 = "AB", S3 = "BB")
  result <- panGenomeBreedr:::genotypes(input_list)

  expect_s3_class(result, "data.frame")
  expect_equal(colnames(result), c("S1", "S2", "S3"))
  expect_equal(as.list(result[1, ]), input_list)

  # Test empty input handling
  empty_result <- panGenomeBreedr:::genotypes(list())
  expect_equal(empty_result$Genotypes, "None Found / Uncallable")
})


test_that("alleles_df assigns correctly formatted column names", {
  # Test standard input
  input_list <- list(a = "C", b = "T")
  result <- panGenomeBreedr:::alleles_df(input_list)

  expect_s3_class(result, "data.frame")
  expect_equal(colnames(result), c("allele_A", "allele_B"))

  # Test empty input handling
  empty_result <- panGenomeBreedr:::alleles_df(list())
  expect_equal(empty_result$Alleles, "None Found / Uncallable")
})


test_that("kasp_color_stat.give correctly classifies plate success status", {
  # Create a mock dataframe covering all 4 classification scenarios
  df_kasp <- data.frame(
    MasterPlate = c("P1", "P1", "P2", "P3", "P4", "P4", "P4"),
    Call = c("A:A", "A:T", "A:A", "Uncallable", "A:A", "T:T", "G:G"), # Simulating counts
    stringsAsFactors = FALSE
  )

  result <- panGenomeBreedr:::kasp_color_stat.give(df_kasp)

  # Assertions
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)

  # Check Status classifications
  # P1 (2 alleles) -> Bi-allelic / Successful
  # P2 (1 allele)  -> Monomorphic / Successful
  # P3 (0 alleles) -> None / Failed!
  # P4 (3 alleles) -> Multi-allelic / Failed!

  expect_equal(result$Status[result$MasterPlate == "P1"], "Successful")
  expect_equal(result$Status[result$MasterPlate == "P3"], "Failed!")
})


test_that("kasp_color_stat.give handles invalid input gracefully", {
  # Missing columns should trigger an error if the function performs internal checks
  expect_error(
    panGenomeBreedr:::kasp_color_stat.give(data.frame(WrongCol = 1)),
    "undefined columns selected"
  )
})


test_that("new_colnames_nsamples renames columns while preserving extra data", {
  # Mock data frame with 4 columns
  df_input <- data.frame(a = 1, b = 2, c = 3, d = 4)

  result <- panGenomeBreedr:::new_colnames_nsamples(
    "Plate",
    "SNP",
    "ID",
    df_input
  )

  # Assertions
  expect_equal(colnames(result), c("Plate", "SNP", "ID", "d"))
  expect_error(panGenomeBreedr:::new_colnames_nsamples(
    "P",
    "S",
    "I",
    data.frame(a = 1)
  ))
})




test_that("Genotypes_user extracts specified sample IDs", {
  df_input <- data.frame(
    Genotype = c("S1", "S2", "S3"),
    Value = c(1, 2, 3),
    stringsAsFactors = FALSE
  )

  # Extract IDs
  result <- panGenomeBreedr:::Genotypes_user(df_input, sample_id = "Genotype")

  expect_equal(result, c("S1", "S2", "S3"))
  expect_type(result, "character")
})


test_that("read_mapfile handles CSV and Excel file types", {
  # Create a temporary CSV file
  tmp_csv <- tempfile(fileext = ".csv")
  write.csv(data.frame(a = 1:2), tmp_csv, row.names = FALSE)

  # Test CSV reading
  res_csv <- panGenomeBreedr:::read_mapfile(tmp_csv)
  expect_s3_class(res_csv, "data.frame")
  expect_equal(colnames(res_csv), "a")

  # Test unsupported file error
  tmp_txt <- tempfile(fileext = ".txt")
  writeLines("test", tmp_txt)
  expect_error(panGenomeBreedr:::read_mapfile(tmp_txt), "Unsupported file type")

  # Clean up
  unlink(c(tmp_csv, tmp_txt))
})



test_that("proc_nd_map_func stops if parents are identical", {
  expect_error(
    panGenomeBreedr:::proc_nd_map_func(rp = 1, dp = 1),
    "recurrent parent must not be the same"
  )
})



test_that("check_sep standardizes separators for different platforms", {
  # Agriplex conversion
  expect_equal(panGenomeBreedr:::check_sep("/"), " / ")

  # KASP preservation
  expect_equal(panGenomeBreedr:::check_sep(":"), ":")

  # Custom separator preservation
  expect_equal(panGenomeBreedr:::check_sep("-"), "-")

  # Validation error
  expect_error(panGenomeBreedr:::check_sep("//"), "length of one character")
})





test_that("safe_grep_match returns expected column matches", {
  choices <- c("SNP_ID", "Chr", "Pos")

  # Match "snp" -> "SNP_ID"
  expect_equal(panGenomeBreedr:::safe_grep_match("snp", choices), "SNP_ID")

  # Fallback to first element if pattern not found
  expect_equal(panGenomeBreedr:::safe_grep_match("missing", choices), "SNP_ID")

  # Handle empty input
  expect_null(panGenomeBreedr:::safe_grep_match("snp", character(0)))
})



test_that("read_vcf_as_df reads standard VCF files", {
  # Create a minimal VCF fixture
  tmp_vcf <- tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
      "1\t100\t.\tA\tT\t.\tPASS\t."
    ),
    tmp_vcf
  )

  res <- panGenomeBreedr:::read_vcf_as_df(tmp_vcf)

  expect_s3_class(res, "data.frame")
  expect_true("CHROM" %in% colnames(res))
  expect_equal(nrow(res), 1)

  unlink(tmp_vcf)
})


test_that("par_missing_dat calculates missingness correctly", {
  # Create data with some NAs
  df_miss <- data.frame(
    S1 = c(NA, "A"),
    S2 = c("B", NA),
    row.names = c("Samp1", "Samp2")
  )

  result <- panGenomeBreedr:::par_missing_dat(df_miss)

  expect_equal(result$Missing_Calls, c(1, 1))
  expect_equal(result$Total_SNPs, c(2, 2))
})





test_that("run_upset_plot returns an upset object", {
  # Minimal mock inputs
  mat <- data.frame(S1 = c(1, 0), S2 = c(0, 1))
  info <- data.frame(marker = c("S1", "S2"), locus = c("L1", "L2"))

  # Run plot generator
  plot <- panGenomeBreedr:::run_upset_plot(mat, info)

  expect_s3_class(plot, "upset")
})


test_that("makeTraitInput generates consistent Shiny tag structures", {
  # Mock a Shiny namespace function
  ns <- function(x) x

  # Generate row for id 1 (should not have a trash button)
  res_1 <- panGenomeBreedr:::makeTraitInput(1, ns)

  # Generate row for id 2 (should have a trash button)
  res_2 <- panGenomeBreedr:::makeTraitInput(2, ns)

  # Verify structure
  expect_s3_class(res_1, "shiny.tag")
  expect_true(grepl("trait_1", as.character(res_1)))

  # Ensure only ID > 1 has the actionButton (trash icon)
  expect_false(grepl("remove_1", as.character(res_1)))
  expect_true(grepl("remove_2", as.character(res_2)))
})

test_that("render_ghost_plate_error produces a valid ggplot object", {
  plot <- panGenomeBreedr:::render_ghost_plate_error("Invalid file format")

  # Verify ggplot integrity
  expect_s3_class(plot, "ggplot")
  expect_equal(plot$labels$title, "PLATE LAYOUT ERROR")
})



test_that("plot_variant_hotspots generates a ggplot object", {
  # Create minimal variant and annotation tables
  variants <- data.frame(
    variant_id = "SNP_100",
    chrom = "Chr1",
    pos = 100,
    variant_type = "SNP"
  )
  annos <- data.frame(
    variant_id = "SNP_100",
    impact = "HIGH"
  )

  # Run plot generator
  p <- panGenomeBreedr:::plot_variant_hotspots(variants, annos)

  # Verify ggplot structure
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Genomic Variant Hotspots")
})


