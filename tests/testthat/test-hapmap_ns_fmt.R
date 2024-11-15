test_that("hapmap_ns_fmt works", {

  # Marker IDs
  snps <- c('snpSB00072', 'snpSB00106', 'snpSB00109', 'Sbv3.1_01_68028666I',
            'Sbv3.1_02_67884158W')

  # Map file for SNPs
  map_file <- data.frame(snpid = snps,
                         chr = c(2, 5, 5, 1, 2),
                         pos = c(61811307, 838874, 1730282, NA, NA))

  # Format marker IDs to hapmap format
  ns_new <- hapmap_ns_fmt(x = snps,
                          map_file = map_file,
                          snpid_col = 'snpid',
                          chr_col = 'chr',
                          pos_col = 'pos',
                          kasp_prefix = 'snpSB',
                          scaffold_prefix = 'Sbv')

  expect_true(all(grepl('^S\\d+_\\d+$', ns_new)), TRUE) # expect TRUE

})
