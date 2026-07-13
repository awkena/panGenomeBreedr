test_that("disconnect_local_db disconnects an active connection and respects quiet mode", {
  mini_folder <- system.file(
    "extdata",
    "mini_curated_sorghum_variant_resource",
    package = "panGenomeBreedr",
    mustWork = FALSE
  )
  skip_if(
    mini_folder == "",
    message = "Mini extdata database folder not found."
  )

  # Establish connection for testing normal verbose path
  con_verbose <- connect_local_db(
    folder_path = mini_folder,
    max_memory = "2GB",
    quiet = TRUE
  )

  expect_message(
    disconnect_local_db(con_verbose, quiet = FALSE),
    "Successfully disconnected from the local database"
  )
  expect_false(dbIsValid(con_verbose))

  # Establish connection for testing quiet mode execution path
  con_quiet <- connect_local_db(
    folder_path = mini_folder,
    max_memory = "2GB",
    quiet = TRUE
  )

  expect_silent(
    disconnect_local_db(con_quiet, quiet = TRUE)
  )
  expect_false(dbIsValid(con_quiet))
})


test_that("disconnect_local_db handles invalid or missing connections gracefully", {
  mini_folder <- system.file(
    "extdata",
    "mini_curated_sorghum_variant_resource",
    package = "panGenomeBreedr",
    mustWork = FALSE
  )
  skip_if(
    mini_folder == "",
    message = "Mini extdata database folder not found."
  )

  # Test behavior with a connection that has already been severed
  con_dead <- connect_local_db(
    folder_path = mini_folder,
    max_memory = "2GB",
    quiet = TRUE
  )
  dbDisconnect(con_dead, shutdown = TRUE)

  expect_message(
    disconnect_local_db(con_dead, quiet = FALSE),
    "The database connection is already closed or invalid"
  )

  # Test behavior when the 'con' object is completely missing from the call
  expect_message(
    disconnect_local_db(quiet = FALSE),
    "The database connection is already closed or invalid"
  )

  # Test behavior when passed an invalid object type entirely
  expect_message(
    disconnect_local_db("not_a_connection", quiet = FALSE),
    "The database connection is already closed or invalid"
  )
})