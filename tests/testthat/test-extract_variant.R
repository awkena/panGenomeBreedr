test_that("extract_variant() extracts variants from tabix-indexed VCF for a gene", {
  skip_on_cran()
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("GenomicRanges")

  # Create temp directory
  dir <- tempdir()

  # ---- Step 1: Create mock GFF with one gene ----
  gff_path <- file.path(dir, "mock.gff3")
  gff_contents <- paste(
    "##gff-version 3",
    "Chr01\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=Sobic.001G000100",
    sep = "\n")

  writeLines(gff_contents, con = gff_path)

  # ---- Step 2: Create mock VCF ----
  vcf_path <- file.path(dir, "mock_chr01.vcf")
  vcf_lines <- c(
    "##fileformat=VCFv4.2",
    "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "Chr01\t150\t.\tA\tG\t.\t.\tANN=G|missense_variant|MODERATE|Sobic.001G000100",
    "Chr01\t180\t.\tC\tT\t.\t.\tANN=T|synonymous_variant|LOW|Sobic.001G000100")

  writeLines(vcf_lines, con = vcf_path)

  # ---- Step 3: Compress and index VCF ----
  vcf_bgz_path <- file.path(dir, "mock_chr01.vcf.gz")
  Rsamtools::bgzip(vcf_path, dest = vcf_bgz_path, overwrite = TRUE)
  Rsamtools::indexTabix(vcf_bgz_path, format = "vcf")

  # ---- Step 4: Run extract_variant() ----
  expect_silent({
    extract_variant(
      cand_gene_id = "Sobic.001G000100",
      gff_path = gff_path,
      vcf_dir = dir,
      vcf_file = basename(vcf_bgz_path),
      output_path = dir,
      outfile_suffix = "test"
    )
  })

  # ---- Step 5: Check that output file exists and is valid ----
  out_file <- file.path(dir, "Sobic.001G000100_test.vcf")
  expect_true(file.exists(out_file))

  contents <- readLines(out_file)
  expect_true(any(grepl("missense_variant", contents)))
  expect_true(any(grepl("synonymous_variant", contents)))
  expect_true(any(grepl("^#", contents)))  # header lines included

  # ---- Step 6: Clean up ----
  unlink(c(gff_path, vcf_path, vcf_bgz_path, paste0(vcf_bgz_path, ".tbi"), out_file))
})
