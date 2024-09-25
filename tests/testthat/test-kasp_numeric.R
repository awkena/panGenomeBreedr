test_that("kasp_numeric works", {

  # Reshape KASP data for each master plate for the beta_carotene data
  dat1 <- panGenomeBreedr::beta_carotene
  plate_wide <- kasp_reshape_wide(x = dat1,
                                  subset = 'MasterPlate',
                                  snp_id = 'SNPID',
                                  geno_call = 'Call',
                                  idvar = "SubjectID",
                                  blank = 'NTC')

  # Get Master Plate 1
  plate1 <- plate_wide$`SE-24-1088_P01_d1`
  rownames(plate1) <- plate1[,1]

  # Convert to numeric format for plotting
  num_geno <- kasp_numeric(x = plate1[,-1],
                          rp_row = 1,
                          dp_row = 7,
                          data_type = 'kasp')

  expect_equal(dim(num_geno), c(94, 4))



})
