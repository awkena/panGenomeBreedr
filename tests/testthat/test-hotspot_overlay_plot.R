
# --- Helper: Create Mock Data ---
# Shaped like the real gene_coord_gff() output: a gene backbone row plus a
# transcript row (with its own ID) and that transcript's CDS/UTR children.
mock_gene_df <- data.frame(
  ID = c(
    "Sobic.005G213600",
    "Sobic.005G213600.1",
    "Sobic.005G213600.1",
    "Sobic.005G213600.1",
    "Sobic.005G213600.1"
  ),
  Feature = c("gene", "mRNA", "five_prime_UTR", "CDS", "three_prime_UTR"),
  Chromosome = "Chr05",
  Start = c(75104500, 75104500, 75104500, 75104601, 75106001),
  End = c(75106500, 75106500, 75104600, 75106000, 75106500),
  Strand = "+",
  stringsAsFactors = FALSE
)

# annotations_df and genotypes_df as hotspot_overlay_plot() merges them
# itself (by variant_id); chrom/pos overlap between the two, so the merge
# suffixes them to chrom.x/pos.x -- matching plot_variant_hotspot()'s
# defaults.
mock_annota_df <- data.frame(
  variant_id = c("INDEL_Chr05_75104881", "SNP_Chr05_75105000"),
  chrom = c("Chr05", "Chr05"),
  pos = c(75104881, 75105000),
  impact = c("HIGH", "LOW"),
  stringsAsFactors = FALSE
)

mock_geno_df <- data.frame(
  variant_id = c("INDEL_Chr05_75104881", "SNP_Chr05_75105000"),
  chrom = c("Chr05", "Chr05"),
  pos = c(75104881, 75105000),
  variant_type = c("INDEL", "SNP"),
  minor_allele_freq = c(0.12, 0.2),
  stringsAsFactors = FALSE
)


# --- Tests ---

test_that("hotspot_overlay_plot generates a patchwork object", {

  # Mock gene_coord_gff() to isolate the wrapper's own logic from the GFF
  # parser; annotations_df/genotypes_df are passed directly, as the real
  # caller (mod_variant_discovery.R) already fetches them beforehand.
  local_mocked_bindings(
    gene_coord_gff = function(gene_name, gff_path) mock_gene_df
  )

  p <- hotspot_overlay_plot(
    gene_name = "Sobic.005G213600",
    gff_path = "dummy.gff3",
    annotations_df = mock_annota_df,
    genotypes_df = mock_geno_df
  )

  # 1. Check if the output is successfully built as a patchwork object
  expect_s3_class(p, c("patchwork", "gg", "ggplot"))

  # 2. Check that the object contains exactly 2 stacked plots
  # patchwork standard indexing treats [[1]] as top, [[2]] as bottom
  expect_s3_class(p[[1]], "ggplot")
  expect_s3_class(p[[2]], "ggplot")
})


test_that("hotspot_overlay_plot successfully adds selected_variants annotation layers", {

  local_mocked_bindings(
    gene_coord_gff = function(gene_name, gff_path) mock_gene_df
  )

  # Generate annotated plot
  p_annotated <- hotspot_overlay_plot(
    gene_name = "Sobic.005G213600",
    gff_path = "dummy.gff3",
    annotations_df = mock_annota_df,
    genotypes_df = mock_geno_df,
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
