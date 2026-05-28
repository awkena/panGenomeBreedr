# ==============================================================================
# Global Test Setup for panGenomeBreedr (DuckDB Backend)
# ==============================================================================

mini_folder <- system.file(
  "extdata/mini_curated_sorghum_variant_resource",
  package = "panGenomeBreedr",
  mustWork = TRUE
)

# Connect to the mini parquet database
con_test <- connect_local_db(
  folder_path = mini_folder,
  max_memory = "1GB",
  quiet = TRUE
)

# Disconnect and release RAM after use.
withr::defer(
  disconnect_local_db(con_test, quiet = TRUE),
  teardown_env()
)
