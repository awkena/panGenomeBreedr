test_that("kasp_color works", {

  set.seed(123)
  dat1 <- data.frame(DaughterPlate = rep('S000940416', time = 96),
                     MasterPlate = rep('SE-24-0391_P01_d2', time = 96),
                     MasterWell = paste0(rep(LETTERS[1:8], times = 12),
                                         rep(sprintf("%02d", 1:12), each = 8)),
                     Call = sample(c('A:A', 'A:T', 'T:T', '?', 'Uncallable'),
                                   size = 96,
                                   replace = TRUE,
                                   prob = c(.325, .25, .325, .05, .05)),
                     X = runif(96),
                     Y = runif(96),
                     SNPID = rep('snpSB00720', time = 96))

  dat1$Call[dat1$MasterWell == 'H11' | dat1$MasterWell == 'H12'] <- 'NTC'

  dat2 <- kasp_color(dat1)

  expect_equal(dim(dat2[[1]]), c(96, 9))
})
