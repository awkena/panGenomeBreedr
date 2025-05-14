test_that("gene_coord_gff() extracts correct coordinates from a local GFF3 file", {
  skip_on_cran()
  skip_if_not_installed("R.utils")

  # Mock GFF3 data
  gff_data <- c("Chr01\tPhytozome\tgene\t100\t500\t.\t+\t.\tID=Sobic.001G000100;Name=Sobic.001G000100",
                "Chr01\tPhytozome\tgene\t600\t900\t.\t-\t.\tID=Sobic.001G000200;Name=Sobic.001G000200",
                "Chr01\tPhytozome\tmRNA\t100\t500\t.\t+\t.\tID=mRNA1;Parent=Sobic.001G000100")

  # Save to a temp plain-text file
  gff_path <- tempfile(fileext = ".gff3")
  writeLines(gff_data, con = gff_path)

  # Also make a GZ version for gzipped test
  gff_path_gz <- paste0(gff_path, ".gz")
  R.utils::gzip(gff_path, destname = gff_path_gz, overwrite = TRUE,
                remove = FALSE)

  # Test plain GFF
  result <- gene_coord_gff("Sobic.001G000100", gff_path)
  expect_type(result, "list")
  expect_named(result, c("chrom", "start", "end"))
  expect_equal(result$chrom, "Chr01")
  expect_equal(result$start, 100)
  expect_equal(result$end, 500)

  # Test gzipped GFF
  result_gz <- gene_coord_gff("Sobic.001G000200", gff_path_gz)
  expect_equal(result_gz$chrom, "Chr01")
  expect_equal(result_gz$start, 600)
  expect_equal(result_gz$end, 900)

  # Test for gene not found
  expect_error(gene_coord_gff("Sobic.UNKNOWN", gff_path), "Gene not found")

  # Clean up
  unlink(gff_path)
  unlink(gff_path_gz)
})
