test_that("connect_local_db successfully connects and initializes views with extdata", {
  # Locate the mini dataset bundled with the package
  mini_folder <- system.file(
    "extdata",
    "mini_curated_sorghum_variant_resource",
    package = "panGenomeBreedr",
    mustWork = FALSE
  )

  # Skip gracefully if the package assets are not built in this environment
  skip_if(
    mini_folder == "",
    message = "Mini extdata database folder not found."
  )

  # Execute successful initialization path
  expect_message(
    con <- connect_local_db(
      folder_path = mini_folder,
      max_memory = "2GB",
      quiet = FALSE
    ),
    "Successfully connected to the local offline database"
  )

  # Ensure an active DBI Connection handle is returned
  expect_s4_class(con, "DBIConnection")
  expect_true(dbIsValid(con))

  # Validate that DuckDB mounted all 4 required views cleanly
  mounted_tables <- dbListTables(con)
  expect_true(all(
    c("genotypes", "variants", "metadata", "annotations") %in% mounted_tables
  ))

  # Disconnect to release session locks
  dbDisconnect(con, shutdown = TRUE)
})

test_that("connect_local_db respects the quiet parameter", {
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

  # Ensure no console messages are triggered when quiet = TRUE
  expect_silent(
    con <- connect_local_db(
      folder_path = mini_folder,
      max_memory = "2GB",
      quiet = TRUE
    )
  )

  dbDisconnect(con, shutdown = TRUE)
})

test_that("connect_local_db stops early if a core Parquet asset is missing", {
  # Build a temporary sandbox directory with a missing component
  corrupted_folder <- file.path(tempdir(), "broken_db_locus")
  dir.create(corrupted_folder, showWarnings = FALSE)
  on.exit(unlink(corrupted_folder, recursive = TRUE), add = TRUE)

  # Write an empty file for only one of the requirements to trip the file checker loop
  writeLines("", file.path(corrupted_folder, "genotypes.parquet"))

  # Ensure it throws the expected critical file omission warning
  expect_error(
    connect_local_db(folder_path = corrupted_folder),
    "Critical file missing: variants.parquet"
  )
})