test_that("pred_summary works", {

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
                     SNPID = rep('snpSB00720', time = 96),
                     Group = sample(c('A:A', 'A:T', 'T:T', '?'), size = 96,
                                    replace = TRUE, prob = c(.325, .25, .325, .1)))

  dat1$Group[dat1$MasterWell == 'H11' | dat1$MasterWell == 'H12'] <- 'NTC'

  dat1$Call[dat1$MasterWell == 'H11' | dat1$MasterWell == 'H12'] <- 'NTC'

  # Subset based on MasterPlate x SNPID combination -- 8 plates
  dat1 <- kasp_color(x = dat1,
                     subset = 'MasterPlate',
                     sep = ':',
                     geno_call = 'Call',
                     uncallable = 'Uncallable',
                     unused = '?',
                     blank = 'NTC')

  dat1 <- pred_summary(x = dat1,
                      geno_call = 'Call',
                      Group_id = 'Group',
                      blank = 'NTC',
                      Group_unknown = '?')

  expect_equal(dim(dat1$summ), c(1, 5))
  expect_equal(length(dat1$plates), 1)
})
