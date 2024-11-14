test_that("proc_kasp works", {

  # Load the beta_carotene KASP data
  dat1 <- panGenomeBreedr::beta_carotene
  plate_wide <- kasp_reshape_wide(x = dat1,
                                  subset = 'MasterPlate',
                                  snp_id = 'SNPID',
                                  geno_call = 'Call',
                                  idvar = "SubjectID",
                                  blank = 'NTC')

  # Get Master Plate 1
  plate1 <- plate_wide$`SE-24-1088_P01_d1`

  # Generate a map for the beta_carotene KASP data
  kasp_map <- data.frame(SNPID = unique(beta_carotene$SNPID),
                         SNPID_2 = paste0('snpSB', 1:4, '|', c(800, 803, 804, 805)),
                         chr = c(1:4),
                         pos = c(800, 803, 804, 805))

  # Process Plate1 to re-order SNPs based on chrom. and position
  proc_plate1 <- proc_kasp(x = plate1,
                           kasp_map = kasp_map,
                           map_snp_id = "SNPID",
                           sample_id = "SubjectID",
                           marker_start = 2,
                           chr = 'chr',
                           chr_pos = 'pos')


  expect_equal(dim(proc_plate1), c(94, 4))

})
