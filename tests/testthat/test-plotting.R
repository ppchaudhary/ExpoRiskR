test_that("at least one exported plotting function exists", {
  plots <- grep("^plot_", ls("package:ExpoRiskR"), value = TRUE)
  expect_true(length(plots) >= 1)
})
