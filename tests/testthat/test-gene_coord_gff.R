library(testthat)

# --- Helper: Create a minimal mock GFF3 file for robust local testing ---
mock_gff_lines <- c(
  "##gff-version 3",
  "Chr05\tphytozomev13\tgene\t75104500\t75106500\t.\t+\t.\tID=Sobic.005G213600;Name=Sobic.005G213600;locusName=Sobic.005G213600",
  "Chr05\tphytozomev13\tmRNA\t75104500\t75106500\t.\t+\t.\tID=Sobic.005G213600.1;Name=Sobic.005G213600.1;Parent=Sobic.005G213600",
  "Chr05\tphytozomev13\tfive_prime_UTR\t75104500\t75104600\t.\t+\t.\tID=Sobic.005G213600.1.UTR5.1;Parent=Sobic.005G213600.1",
  "Chr05\tphytozomev13\tCDS\t75104601\t75106000\t.\t+\t0\tID=Sobic.005G213600.1.CDS.1;Parent=Sobic.005G213600.1",
  "Chr05\tphytozomev13\tthree_prime_UTR\t75106001\t75106500\t.\t+\t.\tID=Sobic.005G213600.1.UTR3.1;Parent=Sobic.005G213600.1",
  # Add a duplicate gene to test the warning trigger
  "Chr05\tphytozomev13\tgene\t1000\t2000\t.\t-\t.\tID=Sobic.DUP;Name=Sobic.DUP",
  "Chr05\tphytozomev13\tgene\t3000\t4000\t.\t-\t.\tID=Sobic.DUP_2;Name=Sobic.DUP"
)

# Wrapper to generate a fresh file per test
setup_mock_gff <- function() {
  temp_mock_gff <- tempfile(pattern = "mock_sorghum_", fileext = ".gff3")
  writeLines(mock_gff_lines, temp_mock_gff)
  return(temp_mock_gff)
}

# --- Tests ---

test_that("gene_coord_gff correctly parses local uncompressed GFF structure", {

  mock_gff <- setup_mock_gff()
  on.exit(unlink(mock_gff)) # Clean up this specific file after test

  res <- gene_coord_gff(gene_name = "Sobic.005G213600", gff_path = mock_gff)

  # Check class and dimensions
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 5) # gene, mRNA, 5'UTR, CDS, 3'UTR

  # Check expected column names
  expected_cols <- c("ID", "Feature", "Chromosome", "Start", "End", "Strand")
  expect_true(all(expected_cols %in% colnames(res)))

  # Check data types
  expect_type(res$Start, "integer")
  expect_type(res$End, "integer")
  expect_type(res$ID, "character")

  # Check coordinate accuracy for the parent gene
  gene_row <- res[res$Feature == "gene", ]
  expect_equal(gene_row$Start, 75104500)
  expect_equal(gene_row$End, 75106500)
  expect_equal(gene_row$Strand, "+")
})


test_that("gene_coord_gff handles missing genes appropriately", {

  mock_gff <- setup_mock_gff()
  on.exit(unlink(mock_gff))

  expect_error(
    gene_coord_gff(gene_name = "Sobic.NOT_REAL", gff_path = mock_gff),
    "Gene not found in GFF."
  )

})


test_that("gene_coord_gff warns when multiple gene copies exist", {

  mock_gff <- setup_mock_gff()
  on.exit(unlink(mock_gff))

  expect_warning(
    gene_coord_gff(gene_name = "Sobic.DUP", gff_path = mock_gff),
    "Multiple gene matches found; using the first match."
  )

})


test_that("gene_coord_gff correctly processes online gzipped files", {

  # Skip on CRAN to prevent build failures if GitHub servers timeout
  skip_on_cran()

  # Skip if no internet connection is available locally
  skip_if_offline()

  gff_url <- "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz"

  # This tests the URL downloading and R.utils::gunzip flow
  res <- gene_coord_gff(gene_name = "Sobic.005G213600", gff_path = gff_url)

  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0)
  expect_true(all(c("gene", "mRNA", "CDS") %in% res$Feature))

})
