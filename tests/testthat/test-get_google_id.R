test_that("get_google_id() extracts correct IDs from Google Drive links", {
  # Folder link
  folder_link <- "https://drive.google.com/drive/folders/1BotxaUb5emlrtgo473db3gDTUCLzKi70?usp=sharing"
  folder_id <- get_google_id(drive_link = folder_link, is.folder = TRUE)
  expect_equal(folder_id, "1BotxaUb5emlrtgo473db3gDTUCLzKi70")

  # File link
  file_link <- "https://drive.google.com/file/d/1XjYyJ2JLywbbniIU6oUIIxAmEBKfmHpz/view?usp=sharing"
  file_id <- get_google_id(drive_link = file_link, is.folder = FALSE)
  expect_equal(file_id, "1XjYyJ2JLywbbniIU6oUIIxAmEBKfmHpz")

  # Broken folder link (missing ID)
  broken_folder <- "https://drive.google.com/drive/folders/"
  broken_id <- get_google_id(drive_link = broken_folder, is.folder = TRUE)
  expect_equal(broken_id, character())  # returns empty string

  # Broken file link (non-matching pattern)
  bad_file <- "https://example.com/file/xyz"
  bad_id <- get_google_id(drive_link = bad_file, is.folder = FALSE)
  expect_equal(bad_id, character())  # no match
})
