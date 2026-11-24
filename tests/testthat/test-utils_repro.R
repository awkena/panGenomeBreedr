test_that("repro_call writes one argument per line", {
  code <- panGenomeBreedr:::repro_call(
    "kasp_color",
    list(
      x = panGenomeBreedr:::repro_sym("kasp_data"),
      others = c("Missing", "Bad"),
      assign_cols = c(FAM = "blue", HEX = "gold"),
      threshold = 0.95,
      Group_id = NULL,
      scale = TRUE
    ),
    assign = "kasp_coded"
  )

  expect_equal(code, c(
    "kasp_coded <- kasp_color(",
    "  x = kasp_data,",
    "  others = c(\"Missing\", \"Bad\"),",
    "  assign_cols = c(FAM = \"blue\", HEX = \"gold\"),",
    "  threshold = 0.95,",
    "  Group_id = NULL,",
    "  scale = TRUE",
    ")"
  ))

  # The generated code parses and recreates the values
  args <- eval(parse(text = sub("^kasp_coded <- kasp_color", "list", paste(code, collapse = "\n"))),
               envir = list(kasp_data = "data"))
  expect_equal(args$assign_cols, c(FAM = "blue", HEX = "gold"))
  expect_null(args$Group_id)
})

test_that("repro_file_name keeps bundle paths inside data/", {
  expect_equal(panGenomeBreedr:::repro_file_name("../../etc/passwd"), "passwd")
  expect_equal(panGenomeBreedr:::repro_file_name("my file (1).csv"), "my_file__1_.csv")
  expect_equal(panGenomeBreedr:::repro_data_path("a.csv"), file.path("data", "a.csv"))
})

test_that("repro_recorder collects upstream steps in order and reports missing ones", {
  rec <- panGenomeBreedr:::repro_recorder()
  step <- panGenomeBreedr:::repro_step
  rec$register("read", function() step("Read", "a <- 1"))
  rec$register("color", function() step("Colour", "b <- a + 1", needs = "read"))
  rec$register("qc", function() step("QC", "c <- b + 1", needs = "color"))
  rec$register("pred", function() stop("not run yet"))

  res <- rec$collect("qc")
  expect_equal(names(res$steps), c("read", "color", "qc"))
  expect_length(res$missing, 0)

  res2 <- rec$collect("pred")
  expect_equal(res2$missing, "pred")
  expect_null(rec$collect("unknown")$steps$unknown)
})

test_that("KASP bundle reruns outside the app and matches direct calls", {
  skip_if_not_installed("panGenomeBreedr")
  skip_if_not(requireNamespace("yaml", quietly = TRUE))

  path <- system.file("extdata", "Genotyping_141.010_01.csv",
                      package = "panGenomeBreedr")
  read_step <- panGenomeBreedr:::repro_kasp_read_step(list(
    path = path,
    name = "Genotyping_141.010_01.csv",
    datatype = "raw",
    row_tags = c("Statistics", "DNA", "SNPs", "Scaling", "Data"),
    spacing = 2
  ))
  color_args <- list(
    uncallable = "Uncallable", unused = "?", blank = "NTC",
    others = c("Missing", "Bad", "Dupe", "Over", "Short"), sep = ":",
    subset = "plates", geno_call = "Call",
    assign_cols = c(FAM = "blue", HEX = "gold", het = "forestgreen")
  )
  pred_args <- list(
    Group_unknown = "?", blank = "NTC", snp_id = "SNPID",
    Group_id = "Group", geno_call = "Call", rate_out = TRUE
  )
  steps <- list(
    read = read_step,
    nsamples = panGenomeBreedr:::repro_kasp_nsamples_step(list(
      subset = "plates", snp_id = "SNPID", plate_id = "MasterPlate"
    )),
    color = panGenomeBreedr:::repro_kasp_color_step(color_args),
    pred = panGenomeBreedr:::repro_kasp_pred_step(pred_args)
  )

  zip_file <- tempfile(fileext = ".zip")
  panGenomeBreedr:::repro_bundle(zip_file, steps, "Test analysis",
                                 script_name = "marker_validation.R")
  out <- tempfile("bundle_")
  utils::unzip(zip_file, exdir = out)

  expect_setequal(
    list.files(out, recursive = TRUE),
    c("marker_validation.R", "settings.yml", "sessionInfo.txt", "README.txt",
      "data/Genotyping_141.010_01.csv")
  )
  settings <- yaml::read_yaml(file.path(out, "settings.yml"))
  expect_equal(settings$steps$color$geno_call, "Call")

  # Run the script in its own folder and environment
  env <- new.env()
  withr::with_dir(out, sys.source(file.path(out, "marker_validation.R"), envir = env))
  expect_true(file.exists(file.path(out, "results", "prediction_summary.csv")))
  expect_true(file.exists(file.path(out, "results", "samples_per_plate.csv")))

  # Same numbers as calling the functions directly
  kasp <- read_kasp_csv(path)
  kd <- panGenomeBreedr:::plates_col(kasp$Data)
  direct <- do.call(pred_summary,
                    c(list(x = do.call(kasp_color, c(list(x = kd), color_args))),
                      pred_args))
  expect_equal(env$pred$summ, direct$summ)
})

test_that("KASP marker design script passes the app's settings through", {
  vcf <- system.file("extdata", "Sobic.002G302700_SNP_snpeff.vcf",
                     package = "panGenomeBreedr")
  genome <- tempfile(fileext = ".gz")
  writeLines("placeholder", genome)

  design <- list(
    variant_id_col = "ID", chrom_col = "CHROM", pos_col = "POS",
    ref_al_col = "REF", alt_al_col = "ALT", geno_start = 10,
    marker_IDs = "SNP_Chr02_69200443", chr = "Chr02",
    region_name = "ma1", maf = 0.05
  )
  intertek <- list(genome_version = "Sbv5.1", region_name = "ma1",
                   trait = "Example trait", owner = "panGenomeBreedr example")
  step <- panGenomeBreedr:::repro_md_step(list(
    source = "vcf",
    variant_file = vcf, variant_name = "Sobic.002G302700_SNP_snpeff.vcf",
    genome_file = genome, genome_name = "Chr02.fa.gz",
    genome_url = "https://example.org/Chr02.fa.gz",
    design = design, intertek = intertek
  ))

  zip_file <- tempfile(fileext = ".zip")
  panGenomeBreedr:::repro_bundle(zip_file, list(design = step),
                                 "Test design", script_name = "kasp_marker_design.R")
  out <- tempfile("bundle_")
  utils::unzip(zip_file, exdir = out)
  expect_true(all(file.exists(file.path(out, c(
    "kasp_marker_design.R", "data/Sobic.002G302700_SNP_snpeff.vcf",
    "data/Chr02.fa.gz"
  )))))

  # Run the script with the two panGB calls replaced by recorders
  calls <- list()
  env <- new.env()
  env$library <- function(...) invisible(NULL)
  env$kasp_marker_design <- function(...) {
    calls$design <<- list(...)
    list(marker_data = data.frame(SNP_Name = "m1"), plot = list())
  }
  env$make_intertek_table <- function(...) {
    calls$intertek <<- list(...)
    data.frame(SNP_Name = "m1")
  }
  withr::with_dir(out, sys.source(file.path(out, "kasp_marker_design.R"), envir = env))

  expect_equal(calls$design[names(design)], design)
  expect_equal(calls$design$vcf_file, "data/Sobic.002G302700_SNP_snpeff.vcf")
  expect_null(calls$design$gt_df)
  expect_equal(calls$design$genome_file, "data/Chr02.fa.gz")
  expect_equal(calls$intertek[names(intertek)], intertek)
  expect_true(file.exists(file.path(out, "results", "KASP_markers_for_Intertek.csv")))
})

test_that("Decision support script matches the app's Process Data and RPP", {
  skip_if_not_installed("panGenomeBreedr")
  path <- system.file("extdata", "agriplex_dat.csv", package = "panGenomeBreedr")
  dat <- panGenomeBreedr:::read_mapfile(path)
  dat <- dat[dat[["Batch"]] == "3", ]

  for (flt in list(c(TRUE, TRUE, TRUE, TRUE), c(FALSE, FALSE, FALSE, FALSE),
                   c(TRUE, FALSE, TRUE, FALSE))) {
    res <- panGenomeBreedr:::proc_nd_map_func(
      data = dat, marker_start = 7, sample_id = "Genotype", marker_sep = "_",
      apply_par_poly = flt[1], apply_par_miss = flt[2],
      apply_geno_good = flt[3], apply_par_homo = flt[4],
      calls_sep = " / ", data_type = "agriplex", rp = 1, dp = 3,
      Prefix = "S", feedback = "no"
    )
    process <- panGenomeBreedr:::repro_ds_process_step(list(
      data_name = "agriplex_dat.csv", data_path = path,
      cluster_by = "Batch", cluster_focus = "3", marker_start = 7,
      sample_id = "Genotype", calls_sep = " / ", data_type = "agriplex",
      rp = 1, dp = 3, rp_name = res$rp_name, dp_name = res$dp_name,
      apply_par_poly = flt[1], apply_par_miss = flt[2],
      apply_geno_good = flt[3], apply_par_homo = flt[4],
      feedback = "no", marker_sep = "_", prefix = "S",
      map_cols = c(snpid = "snpid", chr = "chr", pos = "pos")
    ))
    rpp_args <- list(map_chr = "chr", map_pos = "pos", map_snp_ids = "snpid",
                     rp = res$rp_name, na_code = -5, weighted = TRUE)
    rpp <- panGenomeBreedr:::repro_ds_rpp_step(list(
      rpp = rpp_args, plot = NULL, sample_col = "sample_id", rpp_col = "total_rpp",
      threshold = 0.9, threshold_note = "", parents = c(res$rp_name, res$dp_name),
      replicates = character()
    ))

    zip_file <- tempfile(fileext = ".zip")
    panGenomeBreedr:::repro_bundle(zip_file, list(process = process, rpp = rpp),
                                   "Test", script_name = "mabc.R")
    out <- tempfile("bundle_")
    utils::unzip(zip_file, exdir = out)
    env <- new.env()
    withr::with_dir(out, sys.source(file.path(out, "mabc.R"), envir = env))

    expect_identical(env$geno_num, res$proc_kasp_f)
    expect_equal(env$map, res$mapfile, ignore_attr = TRUE)
    expect_equal(env$rpp, do.call(calc_rpp_bc, c(list(x = res$proc_kasp_f,
                                                      map_file = res$mapfile),
                                                 rpp_args)))
  }
})

test_that("Variant discovery script runs every step with the app's settings", {
  vd <- function(f, ...) get(f, envir = asNamespace("panGenomeBreedr"))(...)
  filters <- list(
    genotypes = list(list(column = "minor_allele_freq", value = c(0.1, 0.5))),
    variants = list(),
    annotations = list(list(column = "impact", value = c("HIGH", "MODERATE")))
  )
  steps <- list(
    connect = vd("repro_vd_connect_step", list(mode = "local", db_path = "/data/my db")),
    region = vd("repro_vd_region_step", list(method = "coords", mode = "local",
                                            chrom = "Chr05", start = 100, end = 900)),
    filters = vd("repro_vd_filters_step", list(
      mode = "local", filters = filters,
      sample_filters = list(country = "Ethiopia"),
      plot_filters = filters[c("variants", "annotations", "genotypes")]
    )),
    hotspot = vd("repro_vd_hotspot_step", list(method = "coords", text_sz = 2.5,
                                              plot_annotations = TRUE)),
    ld = vd("repro_vd_ld_step", list(target_ids = NULL, metric = "R2", threshold = 0.2,
                                    block_threshold = 0.8, filled = FALSE,
                                    palette_option = "magma")),
    pcil_data = vd("repro_vd_pcil_data_step", list(mode = "local")),
    pcil_families = vd("repro_vd_pcil_families_step", list(mode = "local", selection = "v2")),
    pcil_positive = vd("repro_vd_pcil_positive_step", list(
      type = "position", variant_ids = c("v2", "v3"), sel = 5, donor_thresh = 0.75,
      block_quantile = 0.75, F_quantile = 0.25, window = 1000,
      use_families = TRUE, use_family_ids = TRUE, global_ids = NULL
    )),
    pcil_negative = vd("repro_vd_pcil_negative_step", list(
      n_neg = 3, use_families = TRUE, use_family_ids = TRUE, global_ids = c("L1", "L2")
    )),
    kasp = vd("repro_vd_kasp_step", list(
      mode = "local", candidates = NULL, meta_cols = c("chrom", "pos"),
      genome_name = "Chr05.fa.gz",
      design = list(variant_id_col = "variant_id", chrom_col = "chrom", pos_col = "pos",
                    ref_al_col = "ref", alt_al_col = "alt", geno_start = 10,
                    marker_IDs = "v2", region_name = "r1", maf = 0.05),
      intertek = list(genome_version = "v5", region_name = "r1", trait = "t", owner = "o")
    )),
    disconnect = vd("repro_vd_disconnect_step")
  )
  script <- vd("repro_script", steps, "Test", script_name = "variant_discovery.R",
               notes = panGenomeBreedr:::repro_vd_notes)
  expect_silent(parse(text = script))
  expect_match(script, "db_path <- \"/data/my db\"", fixed = TRUE)
  expect_false(grepl("data/", gsub("/data/my db", "", script, fixed = TRUE), fixed = TRUE))

  # Fake database tables
  variants <- data.frame(variant_id = paste0("v", 1:4), chrom = "Chr05",
                         pos = c(100, 300, 500, 700), minor_allele_freq = c(0.05, 0.2, 0.3, 0.4))
  annotations <- data.frame(variant_id = paste0("v", c(1:4, 2)),
                            impact = c("LOW", "HIGH", "MODERATE", "MODIFIER", "LOW"))
  genotypes <- cbind(variants, S1 = 0, S2 = 2, S3 = 1)
  calls <- list()
  rec <- function(name, value) {
    force(value)
    function(...) {
      calls[[name]] <<- c(calls[[name]], list(list(...)))
      value
    }
  }
  plot <- ggplot2::ggplot()
  env <- new.env()
  env$connect_local_db <- rec("connect", "CON")
  env$disconnect_local_db <- rec("disconnect", NULL)
  env$fetch_table_region <- function(con, table_name, chrom, start, end) {
    calls$region <<- c(calls$region, list(list(con, table_name, chrom, start, end)))
    get(table_name)
  }
  env$summarize_annotations <- rec("summary", list(annotation_summary = data.frame(a = 1),
                                                   impact_summary = data.frame(a = 1),
                                                   variant_type_totals = data.frame(a = 1)))
  env$fetch_accession_metadata <- rec("meta", data.frame(lib = c("S1", "S3"),
                                                         country = c("Ethiopia", "Sudan")))
  env$plot_variant_hotspot <- rec("hotspot", plot)
  env$calculate_LD <- rec("ld", data.frame(a = 1))
  env$plot_ld_geodesic <- rec("ld_plot", list(plot = plot, table = data.frame(b = 1)))
  env$fetch_pcil_data <- rec("pcil_data", "PCIL")
  env$fetch_pcil_families_by_variant <- rec("families", list(
    pcil_family_summary = data.frame(a = 1),
    pcil_summary = data.frame(sample_id = "L1", selection = "v2", other = 1)
  ))
  env$fetch_pcil_positive <- rec("positive", list(best_lines = data.frame(a = 1),
                                                  pcil_positive = data.frame(a = 1)))
  env$plot_pcil_positive <- rec("pos_plot", list(plot))
  env$plot_pcil_best_lines <- rec("best_plot", list(plot))
  env$fetch_pcil_negative <- rec("negative", list(pairs_best = data.frame(a = 1),
                                                  pairs_extended = data.frame(a = 1)))
  env$plot_pcil_pairs <- rec("pairs_plot", list(plot))
  env$fetch_genotypes_by_id <- rec("geno_by_id", genotypes)
  env$kasp_marker_design <- rec("design", list(plot = list(plot)))
  env$make_intertek_table <- rec("intertek", data.frame(SNP_Name = "m1"))
  env$library <- function(...) invisible(NULL)

  out <- tempfile("vd_")
  dir.create(out)
  writeLines(script, file.path(out, "variant_discovery.R"))
  writeLines("fasta", file.path(out, "Chr05.fa.gz"))
  withr::with_dir(out, sys.source("variant_discovery.R", envir = env))

  expect_equal(calls$connect[[1]]$folder_path, "/data/my db")
  expect_equal(calls$region[[1]], list("CON", "variants", "Chr05", 100, 900))
  # Row filters, sample filters and "Refine plot"
  expect_equal(env$genotypes_shown$variant_id, c("v2", "v3", "v4"))
  expect_equal(setdiff(names(env$genotypes_shown), names(variants)), c("S1"))
  expect_equal(env$annotations_shown$variant_id, c("v2", "v3"))
  expect_equal(env$selected_variants, c("v2", "v3"))
  # Hotspot: most severe impact per selected variant
  hs <- calls$hotspot[[1]]$var_df
  expect_equal(hs$variant_id, c("v2", "v3"))
  expect_equal(hs$impact, c("HIGH", "MODERATE"))
  expect_identical(calls$ld[[1]]$df, env$genotypes_shown)
  expect_equal(calls$ld_plot[[1]]$palette_option, "magma")
  # PCIL
  expect_equal(calls$families[[1]]$selection, "v2")
  pos <- calls$positive[[1]]
  expect_equal(pos$variants_select_geno$variant_id, c("v2", "v3"))
  expect_equal(pos$available_ids, data.frame(sample_id = "L1", selection = "v2"))
  expect_equal(pos$window, 1000)
  expect_equal(calls$negative[[1]]$global_available_ids, c("L1", "L2"))
  expect_equal(calls$negative[[1]]$n_neg, 3)
  # KASP
  expect_equal(calls$geno_by_id[[1]]$variant_ids, paste0("v", 1:4))
  expect_equal(calls$design[[1]]$genome_file, "Chr05.fa.gz")
  expect_equal(calls$design[[1]]$chr, "Chr05")
  expect_equal(calls$intertek[[1]]$owner, "o")
  expect_length(calls$disconnect, 1)
  expect_true(all(file.exists(file.path(out, "results", c(
    "region_variants.csv", "filtered_genotypes.csv", "hotspot_plot.pdf",
    "ld_geodesic_plot.pdf", "pcil_positive_best_lines.csv", "pcil_negative_pairs.pdf",
    "KASP_markers_for_Intertek.csv"
  )))))
})

test_that("Variant discovery script for an online database and a gene", {
  vd <- function(f, ...) get(f, envir = asNamespace("panGenomeBreedr"))(...)
  steps <- list(
    connect = vd("repro_vd_connect_step", list(mode = "online", api_url = "https://api.example.org")),
    region = vd("repro_vd_region_step", list(method = "gene", mode = "online",
                                            gene_name = "Sobic.005G213600",
                                            gff_source = "file", gff_name = "my genes.gff3")),
    filters = vd("repro_vd_filters_step", list(mode = "online", filters = list())),
    hotspot = vd("repro_vd_hotspot_step", list(method = "gene", gene_name = "Sobic.005G213600",
                                              text_sz = 3))
  )
  script <- vd("repro_script", steps, "Test", script_name = "variant_discovery.R",
               notes = panGenomeBreedr:::repro_vd_notes)
  expect_silent(parse(text = script))
  expect_match(script, "set_api_url(\n  url = \"https://api.example.org\"\n)", fixed = TRUE)
  expect_match(script, "connect_db_mode = \"online\"", fixed = TRUE)
  expect_false(grepl("con = con", script, fixed = TRUE))
  expect_match(script, "gff_file <- \"my_genes.gff3\"", fixed = TRUE)
  expect_match(script, "panGenomeBreedr:::hotspot_overlay_plot(", fixed = TRUE)
  expect_match(script, "selected_variants <- NULL", fixed = TRUE)
})
