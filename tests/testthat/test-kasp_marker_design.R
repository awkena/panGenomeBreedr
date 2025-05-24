test_that("kasp_marker_design works correctly for substitution variant", {

  skip_if_offline()  # Skip if there's no internet access

  # Skip test if any required Bioconductor package is missing
  bioc_pkgs <- c("Biostrings", "BSgenome", "GenomicRanges", "IRanges", "msa")
  missing_pkgs <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    skip(paste("Missing Bioconductor packages:",
               paste(missing_pkgs, collapse = ", ")))
  }

  # Define input paths
  path1 <- "https://raw.githubusercontent.com/awkena/panGB/main/Chr02.fa.gz"
  path2 <- system.file("extdata", "Sobic.002G302700_SNP_snpeff.vcf",
                       package = "panGenomeBreedr", mustWork = TRUE)

  tmp_dir <- tempdir()
  plot_path <- file.path(tmp_dir, "plots")
  dir.create(plot_path, showWarnings = FALSE)

  # Run marker design
  marker <- kasp_marker_design(vcf_file = path2,
                               variant_id_col = "ID",
                               chrom_col = "CHROM",
                               pos_col = "POS",
                               ref_al_col = "REF",
                               alt_al_col = "ALT",
                               genome_file = path1,
                               geno_start = 10,
                               marker_ID = "SNP_Chr02_69200443",
                               chr = "Chr02",
                               plot_draw = TRUE,
                               plot_file = plot_path,
                               region_name = "ma1",
                               maf = 0.05)

  # === ASSERTIONS ===
  expect_s3_class(marker, "data.frame")
  expect_named(marker, c("SNP_Name", "SNP", "Marker_Name",
                         "Chromosome", "Chromosome_Position",
                         "Sequence", "ReferenceAllele", "AlternativeAllele"))
  expect_equal(nrow(marker), 1)

  pdf_files <- list.files(plot_path, pattern = "\\.pdf$", full.names = TRUE)
  expect_true(length(pdf_files) > 0)
  expect_true(file.exists(pdf_files[1]))

  # Clean up
  unlink(plot_path, recursive = TRUE)
})
