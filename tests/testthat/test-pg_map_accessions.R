test_that("pg_map_accessions enforces coordinate column requirements", {
  # Provide a data frame entirely missing the required mapping columns
  bad_meta <- data.frame(
    lib = c("S1", "S2"),
    countryorigin = c("Ghana", "Togo")
  )

  expect_error(
    pg_map_accessions(bad_meta),
    "Metadata must contain 'lat' and 'lon' columns"
  )
})


test_that("pg_map_accessions handles uncoercible coordinates and empty subsets", {
  # Provide a data frame where the coordinates are corrupted or explicitly NA
  na_meta <- data.frame(
    lib = c("S1", "S2"),
    lat = c("missing_data", NA),
    lon = c("unknown", NA),
    countryorigin = c("Ghana", "Togo")
  )

  # The function should attempt to coerce them to numeric, fail (generating NAs),
  # filter out the NAs, realize the dataframe is now empty, and throw this error.
  expect_error(
    pg_map_accessions(na_meta),
    "No samples with valid latitude and longitude found"
  )
})


test_that("pg_map_accessions successfully builds a leaflet htmlwidget", {
  # Safely abort this specific test block if the testing environment lacks leaflet
  skip_if_not_installed("leaflet")
  skip_if_not_installed("tools")

  # Build a clean, valid metadata table with realistic coordinate values
  valid_meta <- data.frame(
    lib = c("S1", "S2", "S3"),
    lat = c(5.6037, 9.4005, 6.1333),
    lon = c(-0.1870, -0.8393, 1.2167),
    countryorigin = c("Ghana", "Ghana", "Togo"),
    stringsAsFactors = FALSE
  )

  # Execute the mapping function
  map_obj <- pg_map_accessions(valid_meta, color_by = "countryorigin")

  # Validate that the resulting object is a properly structured htmlwidget
  expect_s3_class(map_obj, "leaflet")
  expect_s3_class(map_obj, "htmlwidget")

  # Ensure the leaflet object contains the expected internal data structures
  expect_true("x" %in% names(map_obj))
})