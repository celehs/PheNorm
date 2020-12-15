
test_that("PheNorm_Prob probs has right dimension",{
  set.seed(1234)
  fit.dat <- read.csv("https://raw.githubusercontent.com/celehs/PheNorm/master/data-raw/data.csv")
  fit.phenorm=PheNorm.Prob("ICD", "utl", fit.dat, nm.X = NULL,
                           corrupt.rate=0.3, train.size=nrow(fit.dat))
  expect_equal(dim(fit.dat)[1], length(fit.phenorm$probs))
})

test_that("PheNorm_Prob works",{
  set.seed(1234)
  fit.dat <- read.csv("https://raw.githubusercontent.com/celehs/PheNorm/master/data-raw/data.csv")
  fit.phenorm=PheNorm.Prob("ICD", "utl", fit.dat, nm.X = NULL,
                           corrupt.rate=0.3, train.size=nrow(fit.dat))
  expect_equal(round(fit.phenorm$probs[1],6), 0.466247)

})
