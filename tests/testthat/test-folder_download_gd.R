test_that("folder_download_gd() downloads files from a public Google Drive folder", {
  skip_on_cran()
  skip_if_offline()
  skip_if_not_installed("googledrive")

  test_folder_link <- "https://drive.google.com/drive/folders/1BotxaUb5emlrtgo473db3gDTUCLzKi70?usp=sharing"

  downloaded_paths <- folder_download_gd(drive_link = test_folder_link,
                                         output_path = tempdir(),
                                         is.folder = TRUE)

  expect_true(is.list(downloaded_paths) || is.character(downloaded_paths))
  expect_true(length(downloaded_paths) >= 1)

  sample_file <- unlist(downloaded_paths)[1]
  expect_false(file.exists(sample_file))

  finfo <- file.info(sample_file)
  expect_true(is.na(finfo$size))
  expect_equal(finfo$isdir, NA)
})
