test_conplot_grid_data <- function() {
  epsilon <- 10^seq(-5, -1.7, length.out = 30)
  R0 <- 2^seq(0, 5, length.out = 35)
  prob <- outer(
    epsilon, R0,
    function(epsilon, R0) plogis(3 * log(R0) + 180 * epsilon - 8)
  )
  list(epsilon = epsilon, R0 = R0, prob = prob)
}

test_conplot_disease_data <- function() {
  data.frame(
    epsilon = c(0.001, 0.003),
    R0 = c(4, 10),
    label = c("A", "B")
  )
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
      plot_conplot_grid(
        grid$epsilon, grid$R0, grid$prob,
        label.cex = 0.7,
        show.diseases = FALSE
      )
    )
  })

  expect_false(result$visible)
  expect_type(result$value, "list")
  expect_equal(result$value$orientation, "rows: epsilon; columns: R0")
  expect_equal(result$value$axes, "x: epsilon; y: R0")
  expect_true(all(is.finite(result$value$levels)))
  expect_equal(result$value$fill.breaks, sort(unique(c(0, result$value$high.levels, 1))))
  expect_equal(result$value$label.bg.cex.mult, 4)
  expect_equal(result$value$manual.label.point.cex, 0.7 * 4)
})

test_that("plot_conplot_grid preserves manuscript-style defaults", {
  grid <- test_conplot_grid_data()

  with_test_pdf({
    result <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      disease.data = test_conplot_disease_data()
    )
  })

  expect_equal(result$xlim, c(0, 0.02))
  expect_equal(result$ylim, c(1, 32))
  expect_equal(result$log, "y")
  expect_equal(result$label.cex, 0.5)
  expect_equal(result$manual.label.point.cex, 2)
  expect_true(result$show.manual.labels)
  expect_true(result$show.quadratic)
  expect_true(result$show.local.minimum)
  expect_true(result$diseases.plotted)
})

test_that("manual label background circles scale with label size", {
  grid <- test_conplot_grid_data()

  with_test_pdf({
    small <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      label.cex = 0.5,
      show.diseases = FALSE
    )
  })
  with_test_pdf({
    large <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      label.cex = 1,
      show.diseases = FALSE
    )
  })

  expect_equal(unique(small$manual.labels$label.bg.cex), 2)
  expect_equal(unique(large$manual.labels$label.bg.cex), 4)
  expect_equal(large$manual.label.point.cex / small$manual.label.point.cex, 2)
})

test_that("math labels are prepared through earnmisc nice_text", {
  non.tikz <- burnout:::conplot_nice_label(
    "basic reproduction number ($\\Rn$)",
    tikz.mode = FALSE,
    warn = FALSE
  )
  tikz <- burnout:::conplot_nice_label(
    "basic reproduction number ($\\Rn$)",
    tikz.mode = TRUE,
    warn = FALSE
  )

  expect_s3_class(non.tikz, "latexexpression")
  expect_true(grepl("\\mathcal R_0", tikz, fixed = TRUE))
})

test_that("plot_conplot_grid accepts explicit tikz label mode without compiling LaTeX", {
  grid <- test_conplot_grid_data()

  with_test_pdf({
    result <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      use.tikz = TRUE,
      show.diseases = FALSE
    )
  })

  expect_true(result$use.tikz)
})

test_that("plot_conplot_grid honours caller-level use.tikz for labels", {
  grid <- test_conplot_grid_data()

  with_test_pdf({
    result <- local({
      use.tikz <- TRUE
      plot_conplot_grid(
        grid$epsilon, grid$R0, grid$prob,
        show.diseases = FALSE
      )
    })
  })

  expect_true(result$use.tikz)
})

test_that("plot_conplot_grid works on a normal PDF device", {
  grid <- test_conplot_grid_data()

  with_test_pdf({
    result <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      use.tikz = NULL,
      show.diseases = FALSE
    )
  })

  expect_null(result$use.tikz)
  expect_true(result$show.manual.labels)
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
      colour.legend = TRUE,
      show.diseases = FALSE
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
      filled = FALSE,
      show.diseases = FALSE
    )
  })

  expect_false(result$filled)
  expect_true(all(is.finite(result$levels)))
})

test_that("plot_conplot_grid can disable overlays and annotations", {
  grid <- test_conplot_grid_data()

  with_test_pdf({
    result <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      show.manual.labels = FALSE,
      show.overlays = FALSE,
      show.diseases = FALSE
    )
  })

  expect_false(result$show.manual.labels)
  expect_false(result$show.quadratic)
  expect_false(result$show.local.minimum)
  expect_false(result$diseases.plotted)
})
