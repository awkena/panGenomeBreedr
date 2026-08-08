
# --- Helper: Create Mock Data ---
mock_gene_df <- data.frame(
  ID = "Sobic.005G213600",
  Feature = c("gene", "five_prime_UTR", "CDS", "three_prime_UTR"),
  Chromosome = "Chr05",
  Start = c(75104500, 75104500, 75104601, 75106001),
  End = c(75106500, 75104600, 75106000, 75106500),
  Strand = "+",
  stringsAsFactors = FALSE
)

mock_pg_annota <- data.frame(
  variant_id = c("INDEL_Chr05_75104881", "SNP_Chr05_75105000"),
  impact = c("HIGH", "LOW"),
  stringsAsFactors = FALSE
)

mock_pg_gt <- data.frame(
  variant_id = c("INDEL_Chr05_75104881", "SNP_Chr05_75105000"),
  chrom.x = c("Chr05", "Chr05"),
  pos.x = c(75104881, 75105000),
  variant_type = c("INDEL", "SNP"),
  stringsAsFactors = FALSE
)


# --- Tests ---

test_that("hotspot_overlay_plot generates a patchwork object in online mode", {

  # Mock the internal dependencies to isolate the wrapper's logic
  local_mocked_bindings(
    gene_coord_gff = function(gene_name, gff_path) mock_gene_df,
    pg_query_db = function(table_name, chrom, start, end, gene_name = NULL) {
      if (table_name == "annotations") return(mock_pg_annota)
      if (table_name == "variants") return(mock_pg_gt)
    }
  )

  p <- hotspot_overlay_plot(
    gene_name = "Sobic.005G213600",
    gff_path = "dummy.gff3",
    connect_db_mode = "online"
  )

  # 1. Check if the output is successfully built as a patchwork object
  expect_s3_class(p, c("patchwork", "gg", "ggplot"))

  # 2. Check that the object contains exactly 2 stacked plots
  # patchwork standard indexing treats [[1]] as top, [[2]] as bottom
  expect_s3_class(p[[1]], "ggplot")
  expect_s3_class(p[[2]], "ggplot")
})


test_that("hotspot_overlay_plot processes local DB mode correctly", {

  # Mock the local database connection and query functions
  local_mocked_bindings(
    gene_coord_gff = function(gene_name, gff_path) mock_gene_df,
    connect_local_db = function(folder_path) return("mock_connection"),
    query_db = function(con, table_name, chrom, start, end, gene_name = NULL) {
      if (table_name == "annotations") return(mock_pg_annota)
      if (table_name == "variants") return(mock_pg_gt)
    }
  )

  p <- hotspot_overlay_plot(
    gene_name = "Sobic.005G213600",
    gff_path = "dummy.gff3",
    connect_db_mode = "local",
    local_db_path = "fake/path/to/db"
  )

  expect_s3_class(p, c("patchwork", "gg", "ggplot"))
})


test_that("hotspot_overlay_plot successfully adds selected_variants annotation layers", {

  local_mocked_bindings(
    gene_coord_gff = function(gene_name, gff_path) mock_gene_df,
    pg_query_db = function(table_name, chrom, start, end, gene_name = NULL) {
      if (table_name == "annotations") return(mock_pg_annota)
      if (table_name == "variants") return(mock_pg_gt)
    }
  )

  # Generate annotated plot
  p_annotated <- hotspot_overlay_plot(
    gene_name = "Sobic.005G213600",
    gff_path = "dummy.gff3",
    connect_db_mode = "online",
    selected_variants = c("INDEL_Chr05_75104881")
  )

  # Extract geom classes across both plots in the patchwork to ensure the layers were added
  all_geoms <- c(
    sapply(p_annotated[[1]]$layers, function(l) class(l$geom)[1]),
    sapply(p_annotated[[2]]$layers, function(l) class(l$geom)[1])
  )

  expect_true("GeomTextRepel" %in% all_geoms)
  expect_true("GeomVline" %in% all_geoms)
})


test_that("hotspot_overlay_plot throws an error if local mode is selected without a path", {

  # Because the wrapper now fails fast, we do not need to mock gene_coord_gff here
  expect_error(
    hotspot_overlay_plot(
      gene_name = "Sobic.005G213600",
      gff_path = "dummy.gff3",
      connect_db_mode = "local",
      local_db_path = NULL
    ),
    "Please provide the path to the local DB."
  )

})
