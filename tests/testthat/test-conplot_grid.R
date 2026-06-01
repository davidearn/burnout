test_conplot_grid_data <- function() {
  epsilon <- 10^seq(-5, -1.7, length.out = 30)
  R0 <- 2^seq(0, 5, length.out = 35)
  prob <- outer(
    epsilon, R0,
    function(epsilon, R0) plogis(3 * (log(R0) - log(3))^2 + 180 * epsilon - 8)
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
  expect_true(result$show.local.minimum.label)
  expect_equal(result$local.minimum.label$col, "darkred")
  expect_equal(result$local.minimum.label$cex, result$label.cex)
  expect_equal(result$local.minimum.label$offset.lines, 1)
  expect_true(is.finite(result$local.minimum.label$offset.inches))
  expect_true(all(is.finite(result$local.minimum.label$position)))
  expect_true(all(is.finite(result$local.minimum.label$curve.position)))
  expect_equal(result$local.minimum.label$position[["epsilon"]], 0.01)
  expect_equal(result$local.minimum.label$curve.position[["epsilon"]], 0.01)
  expect_gt(
    result$local.minimum.label$position[["R0"]],
    result$local.minimum.label$curve.position[["R0"]]
  )
  expect_true(is.finite(result$local.minimum.label$srt))
  expect_true(result$show.n.legend)
  expect_equal(result$n.legend$label, "$n = 10^6$")
  expect_equal(result$n.legend$position, "topright")
  expect_equal(result$n.legend$col, "black")
  expect_equal(result$n.legend$cex, 1)
  expect_equal(result$n.legend$bty, "n")
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
  tikz.x <- burnout:::conplot_nice_label(
    "mean infectious period / mean lifetime ($\\varepsilon$)",
    tikz.mode = TRUE,
    warn = FALSE
  )
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
  tikz.n <- burnout:::conplot_nice_label(
    "$n = 10^6$",
    tikz.mode = TRUE,
    warn = FALSE
  )

  expect_true(grepl("\\varepsilon", tikz.x, fixed = TRUE))
  expect_s3_class(non.tikz, "latexexpression")
  expect_true(grepl("\\mathcal R_0", tikz, fixed = TRUE))
  expect_true(grepl("n = 10^6", tikz.n, fixed = TRUE))
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
  expect_true(grepl("\\varepsilon", result$labels$xlab, fixed = TRUE))
  expect_true(grepl("n = 10^6", result$n.legend$plotting.label, fixed = TRUE))
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

test_that("manual labels tolerate sparse NA rows and columns", {
  grid <- test_conplot_grid_data()
  grid$prob[, 1:3] <- NA_real_
  grid$prob[1:2, ] <- NA_real_

  with_test_pdf({
    result <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      show.diseases = FALSE,
      show.local.minimum = FALSE
    )
  })

  expect_true(result$show.manual.labels)
  expect_true(anyNA(result$manual.labels$x) || anyNA(result$manual.labels$y))
})

test_that("local-minimum label can be configured and disabled", {
  grid <- test_conplot_grid_data()

  with_test_pdf({
    default <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      show.diseases = FALSE
    )
  })
  with_test_pdf({
    custom <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      show.diseases = FALSE,
      local.minimum.label.col = "blue",
      local.minimum.label.cex = 0.8,
      local.minimum.label.position = c(0.012, 2.2),
      local.minimum.label.offset.lines = 0
    )
  })
  with_test_pdf({
    raised <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      show.diseases = FALSE,
      local.minimum.label.offset.lines = 2
    )
  })
  with_test_pdf({
    inherited.colour <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      show.diseases = FALSE,
      local.minimum.col = "purple",
      local.minimum.label.col = NULL
    )
  })
  with_test_pdf({
    disabled <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      show.diseases = FALSE,
      show.local.minimum.label = FALSE
    )
  })

  expect_equal(default$local.minimum.label$col, "darkred")
  expect_equal(default$local.minimum.label$label, "minimum persistence probability")
  expect_equal(default$local.minimum.label$position[["epsilon"]], 0.01)
  expect_equal(default$local.minimum.label$curve.position[["epsilon"]], 0.01)
  expect_true(is.finite(default$local.minimum.label$position[["R0"]]))
  expect_gt(
    default$local.minimum.label$position[["R0"]],
    default$local.minimum.label$curve.position[["R0"]]
  )
  expect_equal(default$local.minimum.label$offset.lines, 1)
  expect_gt(default$local.minimum.label$offset.inches, 0)
  expect_true(is.finite(default$local.minimum.label$srt))

  expect_equal(custom$local.minimum.label$col, "blue")
  expect_equal(custom$local.minimum.label$cex, 0.8)
  expect_equal(custom$local.minimum.label$offset.lines, 0)
  expect_equal(unname(custom$local.minimum.label$position), c(0.012, 2.2))
  expect_equal(unname(custom$local.minimum.label$curve.position), c(0.012, 2.2))
  expect_gt(
    raised$local.minimum.label$position[["R0"]],
    default$local.minimum.label$position[["R0"]]
  )
  expect_gt(
    raised$local.minimum.label$offset.inches,
    default$local.minimum.label$offset.inches
  )
  expect_equal(inherited.colour$local.minimum.label$col, "purple")

  expect_false(disabled$show.local.minimum.label)
  expect_null(disabled$local.minimum.label)
})

test_that("local-minimum label requires the local-minimum curve", {
  grid <- test_conplot_grid_data()

  expect_error(
    plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      show.local.minimum = FALSE,
      show.local.minimum.label = TRUE
    ),
    "show.local.minimum.label = TRUE requires show.local.minimum = TRUE"
  )
})

test_that("local-minimum curve tolerates sparse edge profiles", {
  grid <- test_conplot_grid_data()
  grid$prob[, 1:3] <- NA_real_
  grid$prob[1:2, ] <- NA_real_

  with_test_pdf({
    result <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      show.diseases = FALSE,
      show.manual.labels = FALSE,
      show.quadratic = FALSE,
      show.local.minimum = TRUE,
      show.local.minimum.label = FALSE
    )
  })

  expect_true(result$show.local.minimum)
  expect_true(any(is.finite(result$local.minimum$y)))
  expect_null(result$local.minimum.label)
})

test_that("local-minimum overlay is skipped when no finite curve exists", {
  grid <- test_conplot_grid_data()
  grid$prob <- outer(
    grid$epsilon, grid$R0,
    function(epsilon, R0) plogis(3 * log(R0) + 180 * epsilon - 8)
  )

  with_test_pdf({
    result <- NULL
    expect_warning(
      result <- plot_conplot_grid(
        grid$epsilon, grid$R0, grid$prob,
        show.diseases = FALSE,
        show.manual.labels = FALSE,
        show.quadratic = FALSE
      ),
      "local-minimum curve has no finite points"
    )
  })

  expect_true(result$show.local.minimum)
  expect_null(result$local.minimum.label)
  expect_false(any(is.finite(result$local.minimum$y)))
})

test_that("talk figure helper returns full conplot metadata invisibly", {
  script.candidates <- c(
    file.path("sandbox", "tracked", "conplot_funs.R"),
    file.path("..", "..", "sandbox", "tracked", "conplot_funs.R")
  )
  script.file <- script.candidates[file.exists(script.candidates)]
  skip_if(!length(script.file), "conplot talk helper is unavailable")
  script.file <- script.file[[1L]]

  script.env <- new.env(parent = globalenv())
  sys.source(script.file, envir = script.env)

  grid <- test_conplot_grid_data()
  plot.args <- list(
    show.diseases = FALSE,
    show.n.legend = FALSE,
    filled = FALSE
  )

  with_test_pdf({
    result <- withVisible(
      script.env$make_conplot_figure(
        grid = list(epsvals = grid$epsilon, Rvals = grid$R0, prob = grid$prob),
        device.settings = list(output.file = NULL, width = 6, height = 6),
        plot.args = plot.args,
        output.mode = "device"
      )
    )
  })

  expect_false(result$visible)
  expect_type(result$value, "list")
  expect_true(all(is.finite(result$value$local.minimum.label$position)))
  expect_gt(
    result$value$local.minimum.label$position[["R0"]],
    result$value$local.minimum.label$curve.position[["R0"]]
  )
  expect_false(result$value$filled)
})

test_that("talk figure script selects model-specific grid defaults", {
  sirs.candidates <- c(
    file.path("sandbox", "tracked", "conplot_talk_figure_sirs.R"),
    file.path("..", "..", "sandbox", "tracked", "conplot_talk_figure_sirs.R")
  )
  sir.candidates <- c(
    file.path("sandbox", "tracked", "conplot_talk_figure_sir.R"),
    file.path("..", "..", "sandbox", "tracked", "conplot_talk_figure_sir.R")
  )
  sirs.file <- sirs.candidates[file.exists(sirs.candidates)]
  sir.file <- sir.candidates[file.exists(sir.candidates)]
  skip_if(!length(sirs.file), "SIRS conplot talk script is unavailable")
  skip_if(!length(sir.file), "SIR conplot talk script is unavailable")
  sirs.file <- sirs.file[[1L]]
  sir.file <- sir.file[[1L]]

  old.options <- options(burnout.conplot_talk_figure.skip_main = TRUE)
  on.exit(options(old.options), add = TRUE)

  sirs.env <- new.env(parent = globalenv())
  sys.source(sirs.file, envir = sirs.env)

  expect_equal(sirs.env$prob.file.name, "prob.RData.eta=0.01")
  expect_equal(sirs.env$output.basename, "conplot_talk_figure_sirs")
  expect_false(sirs.env$plot.args$show.quadratic)
  expect_true(sirs.env$plot.args$show.local.minimum)
  expect_equal(sirs.env$plot.args$local.minimum.xlow, 0)
  expect_false(sirs.env$plot.args$show.local.minimum.label)

  sir.env <- new.env(parent = globalenv())
  sys.source(sir.file, envir = sir.env)

  expect_equal(sir.env$prob.file.name, "prob_sir.RData")
  expect_equal(sir.env$output.basename, "conplot_talk_figure_sir")
  expect_true(sir.env$plot.args$show.quadratic)
  expect_true(sir.env$plot.args$show.local.minimum)
  expect_true(sir.env$plot.args$show.local.minimum.label)
  expect_false("local.minimum.xlow" %in% names(sir.env$plot.args))
})

test_that("sample-size legend can be configured and disabled", {
  grid <- test_conplot_grid_data()

  with_test_pdf({
    custom <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      show.diseases = FALSE,
      n.legend.label = "$n = 10^5$",
      n.legend.position = c(0.002, 20),
      n.legend.cex = 0.8,
      n.legend.col = "blue",
      n.legend.bty = "o"
    )
  })
  with_test_pdf({
    disabled <- plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      show.diseases = FALSE,
      show.n.legend = FALSE
    )
  })

  expect_equal(custom$n.legend$label, "$n = 10^5$")
  expect_equal(unname(custom$n.legend$position), c(0.002, 20))
  expect_equal(custom$n.legend$cex, 0.8)
  expect_equal(custom$n.legend$col, "blue")
  expect_equal(custom$n.legend$bty, "o")

  expect_false(disabled$show.n.legend)
  expect_null(disabled$n.legend)
})

test_that("sample-size legend position is validated", {
  grid <- test_conplot_grid_data()

  expect_error(
    plot_conplot_grid(
      grid$epsilon, grid$R0, grid$prob,
      n.legend.position = "upperright"
    ),
    "n.legend.position must be a valid legend keyword"
  )
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
      show.n.legend = FALSE,
      show.diseases = FALSE
    )
  })

  expect_false(result$show.manual.labels)
  expect_false(result$show.quadratic)
  expect_false(result$show.local.minimum)
  expect_false(result$show.local.minimum.label)
  expect_false(result$show.n.legend)
  expect_false(result$diseases.plotted)
})
