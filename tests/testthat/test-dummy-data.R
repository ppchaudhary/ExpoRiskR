test_that("generate_dummy_exporisk returns expected components", {
  d <- generate_dummy_exporisk()

  expect_true(is.list(d))
  expect_true(all(c("microbiome","metabolome","exposures","meta") %in% names(d)))

  expect_true(is.matrix(d$microbiome))
  expect_true(is.matrix(d$metabolome))
  expect_true(is.matrix(d$exposures))
  expect_true(is.data.frame(d$meta))

  expect_equal(nrow(d$microbiome), nrow(d$metabolome))
  expect_equal(nrow(d$microbiome), nrow(d$exposures))
  expect_equal(nrow(d$microbiome), nrow(d$meta))
})
