library(tRNAscanImport)

context("dot bracket annotations")
test_that("dot bracket annotations:",{
  db <- c(">>>>>>>....<<<<<<<")
  actual <- getBasePairing(db)[[1]]
  expect_equal(actual$pos, 1:18)
  expect_equal(actual$forward, actual[order(actual$pos, decreasing = TRUE),]$reverse)
  db <- c(">>>>>>>[[[<<<<<<<]]]{{{{}}}}(((())))")
  actual <- getBasePairing(db)[[1]]
  expect_equal(sum(actual$forward),sum(actual$reverse))
  expect_equal(sum(actual[1:7,]$forward),sum(actual[11:17,]$reverse))
  expect_equal(sum(actual[8:10,]$forward),sum(actual[18:20,]$reverse))
  expect_equal(sum(actual[21:24,]$forward),sum(actual[25:28,]$reverse))
  expect_equal(sum(actual[29:32,]$forward),sum(actual[33:36,]$reverse))
})