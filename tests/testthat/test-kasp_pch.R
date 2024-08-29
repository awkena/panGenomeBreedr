test_that("kasp_pch works", {

  Call <- sample(c('A:A', 'A:T', 'T:T', '?', 'Uncallable'),
                 size = 96,
                 replace = TRUE,
                 prob = c(.325, .25, .325, .05, .05))

  geno_pch <- kasp_pch(x = Call)

  expect_equal(dim(geno_pch), c(96, 1))

})
