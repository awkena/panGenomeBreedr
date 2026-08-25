mock_marker_snp <- list(
  marker_data = data.frame(
    SNP_Name = "SNP_Chr03_11361160",
    SNP = "Substitution",
    Marker_Name = "example_snp",
    Chromosome = "Chr03",
    Chromosome_Position = 11361160,
    Sequence = "ACTG...[G/A]...CTGAAA",
    ReferenceAllele = "G",
    AlternativeAllele = "A",
    stringsAsFactors = FALSE
  )
)

# A standard INDEL to test the 'I' suffix
mock_marker_indel <- list(
  marker_data = data.frame(
    SNP_Name = "INDEL_Chr05_12345",
    SNP = "Insertion",
    Marker_Name = "example_indel",
    Chromosome = "Chr05",
    Chromosome_Position = 12345,
    Sequence = "ACTG...[-/T]...CTGAAA",
    ReferenceAllele = "-",
    AlternativeAllele = "T",
    stringsAsFactors = FALSE
  )
)

# A SNP with a non-standard IUPAC pair to test the 'X' suffix fallback
mock_marker_non_iupac <- list(
  marker_data = data.frame(
    SNP_Name = "SNP_Chr02_54321",
    SNP = "Substitution",
    Marker_Name = "example_non_iupac",
    Chromosome = "Chr02",
    Chromosome_Position = 54321,
    Sequence = "ACTG...[C/N]...CTGAAA",
    ReferenceAllele = "C",
    AlternativeAllele = "N",
    stringsAsFactors = FALSE
  )
)

# A combination of multiple markers to test batch processing
mock_marker_multiple <- list(
  marker_data = rbind(
    mock_marker_snp$marker_data,
    mock_marker_indel$marker_data,
    mock_marker_non_iupac$marker_data
  )
)

# Common parameters for the function
genome_version <- "Sbv5.1"
region_name <- "qTest"
trait <- "Test Trait"
owner <- "Test Owner"


test_that("make_intertek_table correctly processes a single SNP", {
  result <- make_intertek_table(
    marker_data = mock_marker_snp,
    genome_version = genome_version,
    region_name = region_name,
    trait = trait,
    owner = owner
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$`SNP Name*`, "Sbv5.1_03_11361160_R") # G/A -> AG -> R
  expect_equal(result$`SNP*`, "G/A")
  expect_equal(result$Trait, trait)
  expect_equal(result$`Gene/QTL`, region_name)
  expect_equal(result$`Owner*`, owner)
  expect_equal(result$`Chromosome Number`, "03")
})

test_that("make_intertek_table correctly processes a single INDEL", {
  result <- make_intertek_table(
    marker_data = mock_marker_indel,
    genome_version = genome_version,
    region_name = region_name,
    trait = trait,
    owner = owner
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$`SNP Name*`, "Sbv5.1_05_12345_I")
  expect_equal(result$`SNP*`, "-/T")
})

test_that("make_intertek_table handles non-IUPAC SNP pairs correctly", {
  result <- make_intertek_table(
    marker_data = mock_marker_non_iupac,
    genome_version = genome_version,
    region_name = region_name,
    trait = trait,
    owner = owner
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$`SNP Name*`, "Sbv5.1_02_54321_X") # C/N is not a standard IUPAC pair
})

test_that("make_intertek_table correctly processes multiple markers", {
  result <- make_intertek_table(
    marker_data = mock_marker_multiple,
    genome_version = genome_version,
    region_name = region_name,
    trait = trait,
    owner = owner
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_equal(result$`SNP Name*`, c("Sbv5.1_03_11361160_R", "Sbv5.1_05_12345_I", "Sbv5.1_02_54321_X"))
  expect_equal(result$`SNP*`, c("G/A", "-/T", "C/N"))
})



# Test for error handling
test_that("make_intertek_table throws errors for invalid input structures", {
  # Not a list
  expect_error(
    make_intertek_table(marker_data = data.frame()),
    "`marker_data` must be a list containing a `\\$marker_data` data.frame."
  )

  # List without $marker_data element
  expect_error(
    make_intertek_table(marker_data = list(wrong_name = data.frame())),
    "`marker_data` must be a list containing a `\\$marker_data` data.frame."
  )

  # Empty data frame
  expect_error(
    make_intertek_table(marker_data = list(marker_data = data.frame())),
    "`marker_data\\$marker_data` must be a non-empty data.frame."
  )
})

test_that("make_intertek_table throws error for missing required columns", {
  bad_data <- mock_marker_snp$marker_data
  bad_data$SNP_Name <- NULL # Remove a required column

  expect_error(
    make_intertek_table(marker_data = list(marker_data = bad_data)),
    "Missing required columns in `marker_data\\$marker_data`: SNP_Name"
  )
})