
test_that("folder_download_gd() downloads files from a public Google Drive folder", {
  skip_on_cran()
  skip_if_offline()
  skip_if_not_installed("googledrive")

  # Use a temporary directory for testing to keep your KNUST project environment clean
  test_out_dir <- file.path(tempdir(), "gd_test")
  if (!dir.exists(test_out_dir)) dir.create(test_out_dir)

  test_folder_link <- "https://drive.google.com/drive/folders/1ZJ2cZcybCvUuwOdB_0kK5kg7StgayuPf?usp=drive_link"

  # Run the download
  downloaded_paths <- folder_download_gd(drive_link = test_folder_link,
                                         output_path = test_out_dir,
                                         is.folder = TRUE)

  # 1. Verify the return type and length
  expect_true(is.list(downloaded_paths) || is.character(downloaded_paths))
  expect_true(length(downloaded_paths) >= 1)

  # 2. Verify that the file actually exists on your local machine
  sample_file <- unlist(downloaded_paths)[1]
  expect_true(file.exists(sample_file)) # Fixed: should be TRUE if download worked

  # 3. Verify file properties (it should have size and not be a directory)
  finfo <- file.info(sample_file)
  expect_false(is.na(finfo$size))
  expect_gt(finfo$size, 0)      # File should not be empty
  expect_false(finfo$isdir)     # It should be a file, not a folder

  # Cleanup after test
  unlink(test_out_dir, recursive = TRUE)
})

# test_that("folder_download_gd() downloads files from a public Google Drive folder", {
#   skip_on_cran()
#   skip_if_offline()
#   skip_if_not_installed("googledrive")
#
#   test_folder_link <- "https://drive.google.com/drive/folders/1BotxaUb5emlrtgo473db3gDTUCLzKi70?usp=sharing"
#
#   downloaded_paths <- folder_download_gd(drive_link = test_folder_link,
#                                          output_path = tempdir(),
#                                          is.folder = TRUE)
#
#   expect_true(is.list(downloaded_paths) || is.character(downloaded_paths))
#   expect_true(length(downloaded_paths) >= 1)
#
#   sample_file <- unlist(downloaded_paths)[1]
#   expect_false(file.exists(sample_file))
#
#   finfo <- file.info(sample_file)
#   expect_true(is.na(finfo$size))
#   expect_equal(finfo$isdir, NA)
# })
