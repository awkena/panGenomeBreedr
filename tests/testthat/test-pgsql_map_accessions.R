test_that("pgsql_map_accessions returns a valid interactive leaflet object", {
  # Skip the test gracefully if the suggested packages are not present
  skip_if_not_installed("leaflet")
  skip_if_not_installed("tools")

  # Mock data with missing coordinates
  mock_metadata <- data.frame(
    lib = c("Line_A", "Line_B", "Line_C", "Line_D"),
    lat = c(5.6, 6.7, NA, 7.8),
    lon = c(-0.2, -1.5, 2.0, NA),
    countryorigin = c("Ghana", "Ghana", "Togo", "Benin"),
    elevation = c(200, 350, 150, 400),
    kmeans_cluster = c(1, 2, 1, 2),
    stringsAsFactors = FALSE
  )

  # Execute
  p_inter <- pgsql_map_accessions(mock_metadata)

  # Core Assertions
  expect_s3_class(p_inter, "leaflet")
  expect_s3_class(p_inter, "htmlwidget")

  # TRACE INDEXING EXPLAINED (Leaflet Edition):
  # Leaflet stores its building steps in the 'x$calls' list.
  # We extract the 'addCircleMarkers' call to inspect the data that was actually plotted.
  marker_calls <- Filter(
    function(call) call$method == "addCircleMarkers",
    p_inter$x$calls
  )

  # Ensure the marker layer was successfully added
  expect_length(marker_calls, 1)

  # Extract the latitudes passed to the map layer
  # The args list for addCircleMarkers has the latitudes at index 1
  passed_lats <- marker_calls[[1]]$args[[1]]

  # We check to verify only 2 complete points (Line_A and Line_B) were passed to the map
  expect_equal(length(passed_lats), 2)
})

test_that("pgsql_map_accessions handles custom coloring without noise", {
  skip_if_not_installed("leaflet")
  skip_if_not_installed("tools")

  mock_metadata <- data.frame(
    lib = c("S1", "S2"),
    lat = c(5.6, 6.7),
    lon = c(-0.2, -1.5),
    countryorigin = c("Ghana", "Ghana"),
    elevation = c(200, 350),
    kmeans_cluster = c("A", "B"),
    stringsAsFactors = FALSE
  )

  # The function is robust and should execute flawlessly without throwing errors
  expect_error(
    pgsql_map_accessions(mock_metadata, color_by = "kmeans_cluster"),
    NA
  )
})

test_that("pgsql_map_accessions triggers graceful stops for bad inputs", {
  skip_if_not_installed("leaflet")
  skip_if_not_installed("tools")

  # Missing coordinate columns
  bad_metadata <- data.frame(
    lib = c("S1", "S2"),
    countryorigin = c("Ghana", "USA")
  )
  expect_error(
    pgsql_map_accessions(bad_metadata),
    "Metadata must contain 'lat' and 'lon' columns"
  )

  # All coordinates are missing
  empty_coords <- data.frame(
    lib = c("S1", "S2"),
    lat = c(NA, NA),
    lon = c(NA, NA),
    countryorigin = c("Ghana", "USA")
  )
  expect_error(
    pgsql_map_accessions(empty_coords),
    "No samples with valid latitude and longitude found"
  )
})