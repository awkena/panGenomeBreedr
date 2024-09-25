test_that("geno_error works", {
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

  # Check for genotype call  error for each SNP
  geno_mat <- geno_error(x = plate1[,-1],
                    rp_row = 1,
                    dp_row = 7,
                    sep = ':',
                    data_type = 'kasp')

  expect_equal(dim(geno_mat$geno_good), c(94, 3))
  expect_equal(dim(geno_mat$geno_err), c(94, 1))

})
