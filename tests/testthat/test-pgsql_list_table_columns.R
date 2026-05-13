test_that("pgsql_list_table_columns returns schema metadata for all tables", {
  skip_if_not_installed("dittodb")
  skip_if_not_installed("RPostgres")
  skip("Mock SQL needs re-recording")

  dittodb::with_mock_path("fixtures", {
    dittodb::with_mock_db({
      # Connect with recorded credentials
      con <- DBI::dbConnect(
        RPostgres::Postgres(),
        dbname = 'sorghum_pangenome_db',
        host = 'localhost',
        port = 5432,
        user = 'israeltawiahtetteh'
      )

      # Test variants table
      cols_var <- pgsql_list_table_columns(con, table_name = "variants")
      expect_s3_class(cols_var, "data.frame")
      expect_true(all(
        c("variant_id", "chrom", "pos") %in% cols_var$column_name
      ))

      #  Test annotations table
      cols_ann <- pgsql_list_table_columns(con, table_name = "annotations")
      expect_true(all(
        c("variant_id", "impact", "gene_name") %in% cols_ann$column_name
      ))

      # Test genotypes table
      cols_geno <- pgsql_list_table_columns(con, table_name = "genotypes")
      expect_true("variant_id" %in% cols_geno$column_name)
      # Check for specific columns if they exist in your schema (e.g., sample names or the calls array)
      expect_true(any(c("chrom", "pos") %in% cols_geno$column_name))

      # Test metadata table
      cols_meta <- pgsql_list_table_columns(con, table_name = "metadata")
      expect_true(all(c("lib", "pinumber") %in% cols_meta$column_name))

      # Test error handling for invalid table name
      expect_error(pgsql_list_table_columns(con, table_name = "non_existent"))

      DBI::dbDisconnect(con)
    })
  })
})
