# test-mod_variant_discovery.R
# UI test
test_that("mod_variant_discovery_ui produces expected UI elements", {
  ui <- mod_variant_discovery_ui("test")
  expect_true(inherits(ui, "shiny.tag.list"))

  html <- as.character(ui)
  expect_true(grepl('id="test-connect_btn"', html))
  expect_true(grepl('id="test-disconnect_btn"', html))
  expect_true(grepl('id="test-db_path"', html))
  expect_true(grepl('id="test-status_badge"', html))
  expect_true(grepl('id="test-nav_id"', html))
  expect_true(grepl('id="test-query_dbase_btn"', html))
  expect_true(grepl('id="test-get_pcv_btn"', html))
})

# Server test for basic operations
test_that("mod_variant_discovery_server handles basic operations", {
  # Mock base R functions first
  with_mocked_bindings(
    file.exists = function(...) TRUE,
    .package = "base",
    {
      # Mock your custom functions (these don't need .package since they're in your namespace)
      with_mocked_bindings(
        list_sqlite_tables = function(...) c("variants", "annotations", "genotypes"),
        variant_impact_summary = function(...) data.frame(impact = c("HIGH", "MODERATE"), count = c(10, 20)),
        summarize_sqlite_tables = function(...) data.frame(table_name = c("variants", "annotations"), row_count = c(100, 200)),
        count_variant_types = function(...) data.frame(variant_type = c("SNP", "INDEL"), count = c(80, 20)),
        variant_stats = function(...) data.frame(stat = c("Total Variants", "SNPs"), value = c(100, 80)),
        query_db = function(...) data.frame(variant_id = 1:5, pos = 1000:1004),
        query_ann_summary = function(...) list(
          annotation_summary = data.frame(annotation = c("missense", "synonymous"), count = c(15, 25)),
          impact_summary = data.frame(impact = c("HIGH", "MODERATE"), count = c(10, 20)),
          variant_type_totals = data.frame(type = c("SNP", "INDEL"), count = c(80, 20))
        ),
        query_by_impact = function(...) data.frame(variant_id = 1:3, impact = rep("HIGH", 3)),
        query_genotypes = function(...) data.frame(variant_id = 1:3, pos = 1000:1002, sample1 = c(0, 1, 2), sample2 = c(1, 1, 0)),
        filter_by_af = function(...) data.frame(variant_id = 1:2, af = c(0.1, 0.2)),
        gene_coord_gff = function(...) list(chrom = "Chr1", start = 1000, end = 2000),
        show_toast = function(...) NULL,
        {
          # Mock DBI functions with namespace specification
          with_mocked_bindings(
            dbConnect = function(...) NULL,
            dbDisconnect = function(...) NULL,
            .package = "DBI",
            {
              testServer(mod_variant_discovery_server, {
                session$setInputs(db_path = "test.db", connect_btn = 1)
                expect_true(rv$connected)
                expect_equal(rv$tables, c("variants", "annotations", "genotypes"))

                session$setInputs(get_info = 1)
                expect_true(!is.null(rv$variant_impact))
                expect_true(!is.null(rv$sqlite_summary))
                expect_true(!is.null(rv$variant_count))
                expect_true(!is.null(rv$variant_stats))

                session$setInputs(get_cord = 1)
                session$setInputs(gene_name = "Sobic.001G001000", input_method = "url",
                                  gff_url = "test.gff", submit = 1)
                expect_equal(values$result$chrom, "Chr1")

                session$setInputs(query_database = "q_entire", table_name = "annotations", query_dbase_btn = 1)
                expect_true(!is.null(values$query_db_val))

                session$setInputs(impact_level = "HIGH", af_range = c(0.05, 0.95), get_pcv_btn = 1)
                expect_true(!is.null(values$query_geno_react))

                session$setInputs(disconnect_btn = 1)
                expect_false(rv$connected)
              })
            }
          )
        }
      )
    }
  )
})

# Server test for handling failed connections
test_that("mod_variant_discovery_server handles connection errors", {
  with_mocked_bindings(
    file.exists = function(...) FALSE,
    .package = "base",
    {
      with_mocked_bindings(
        show_toast = function(...) NULL,
        {
          testServer(mod_variant_discovery_server, {
            session$setInputs(db_path = "invalid.db", connect_btn = 1)
            expect_false(rv$connected)
          })
        }
      )
    }
  )
})

# KASP marker design test
test_that("mod_variant_discovery_server handles KASP marker design", {
  with_mocked_bindings(
    kasp_marker_design_alt = function(...) list(
      marker_data = data.frame(marker_id = "test", sequence = "ACGT"),
      plot = ggplot2::ggplot()
    ),
    show_toast = function(...) NULL,
    {
      testServer(mod_variant_discovery_server, {
        values$query_geno_react <- data.frame(
          variant_id = c("var1", "var2"),
          chrom = c("Chr1", "Chr1"),
          pos = c(1000, 2000),
          ref = c("A", "C"),
          alt = c("G", "T"),
          variant_type = c("SNP", "SNP"),
          sample1 = c(0, 1),
          sample2 = c(1, 2)
        )

        session$setInputs(
          push_1 = 1,
          modal_genome_file = list(datapath = "test.fa"),
          modal_marker_ID = "var1",
          modal_reg_name = "test region",
          modal_maf = 0.05,
          modal_draw_plot = TRUE,
          modal_run_but = 1
        )

        expect_true(!is.null(kasp_design_result()))
      })
    }
  )
})
