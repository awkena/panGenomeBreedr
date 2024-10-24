test_that("ggdat works", {

  # Set path to the directory where your data is located
  path1 <-  system.file("extdata", "agriplex_dat.csv",
                        package = "panGenomeBreedr",
                        mustWork = TRUE)

  # Import raw Agriplex data file
  geno <- read.csv(file = path1, header = TRUE, colClasses = c("character"))

  # Parse snp ids to generate a map file
  snps <- colnames(geno)[-c(1:6)] # Get snp ids
  map_file <- parse_marker_ns(x = snps, sep = '_', prefix = 'S')

  # Process genotype data to re-order SNPs based on chromosome and positions
  stg5 <- proc_kasp(x = geno[geno$Batch == 3,], # stg5 NILs
                    kasp_map = map_file,
                    map_snp_id = "snpid",
                    sample_id = "Genotype",
                    marker_start = 7,
                    chr = 'chr',
                    chr_pos = 'pos')

  map_file <- stg5$ordered_map # Ordered map
  stg5 <- stg5$ordered_geno # Ordered geno

  # Convert to numeric format for plotting
  num_geno <- kasp_numeric(x = stg5,
                           rp_row = 1, # Recurrent parent row ID
                           dp_row = 3, # Donor parent row ID
                           sep = ' / ',
                           data_type = 'agriplex')

  # Convert num_geno to a long format data frame
  df <- gg_dat(num_mat = num_geno,
               map_file = map_file)

  expect_equal(dim(df), c(21789, 5))
})
