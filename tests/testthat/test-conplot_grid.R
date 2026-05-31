test_conplot_grid_data <- function() {
  epsilon <- 10^seq(-4, -1, length.out = 8)
  R0 <- seq(1.05, 4, length.out = 10)
  prob <- outer(
    epsilon, R0,
    function(epsilon, R0) plogis(2 * log(R0) + log10(epsilon) + 3)
  )
  list(epsilon = epsilon, R0 = R0, prob = prob)
}

with_test_pdf <- function(code) {
  pdf.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf.file)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

test_that("plot_conplot_grid returns metadata invisibly for valid grids", {
  grid <- test_conplot_grid_data()

  with_test_pdf({
    result <- withVisible(
      plot_conplot_grid(grid$epsilon, grid$R0, grid$prob, label.cex = 0.7)
    )
  })

  expect_false(result$visible)
  expect_type(result$value, "list")
  expect_equal(result$value$orientation, "rows: epsilon; columns: R0")
  expect_true(all(is.finite(result$value$levels)))
})

test_that("plot_conplot_grid validates probability matrix dimensions", {
  grid <- test_conplot_grid_data()

  expect_error(
    plot_conplot_grid(grid$epsilon, grid$R0, grid$prob[-1, ]),
    "prob must have dimensions length\\(epsilon\\) by length\\(R0\\)"
  )
})

test_that("plot_conplot_grid validates label size", {
  grid <- test_conplot_grid_data()

  expect_error(
    plot_conplot_grid(grid$epsilon, grid$R0, grid$prob, label.cex = 0),
    "label.cex must be a positive finite scalar"
  )
})

test_that("plot_conplot_grid validates probability range", {
  grid <- test_conplot_grid_data()
  grid$prob[1, 1] <- 1.01

  expect_error(
    plot_conplot_grid(grid$epsilon, grid$R0, grid$prob),
    "prob values must be in \\[0, 1\\]"
  )
})

test_that("plot_conplot_grid plots with a colour legend", {
  grid <- test_conplot_grid_data()

  with_test_pdf({
    result <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      colour.legend = TRUE
    )
  })

  expect_true(result$colour.legend)
  expect_true(result$filled)
})

test_that("plot_conplot_grid plots contours without filled background", {
  grid <- test_conplot_grid_data()

  with_test_pdf({
    result <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      filled = FALSE
    )
  })

  expect_false(result$filled)
  expect_true(all(is.finite(result$levels)))
})
