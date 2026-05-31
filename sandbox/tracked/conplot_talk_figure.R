## Recreate the manuscript-style conplot from the preserved probability grid.
##
## Run from the package root:
##
##     Rscript sandbox/tracked/conplot_talk_figure.R
##
## Input:
##
##     ../sources/conplot_standalone/prob.RData
##
## Output by default:
##
##     sandbox/conplot_talk_figure.pdf
##
## Set use.tikz.output to TRUE below to write a standalone tikz .tex file
## instead. The tikz path uses earnmisc::tikz_open(), so plot labels are
## prepared automatically for tikz by earnmisc::nice_text().

## Slide/talk settings. The defaults below keep the manuscript-style appearance.
## Increase label.cex or set colour.legend to TRUE for presentation variants.
use.tikz.output <- FALSE

plot.settings <- list(
    output.file = file.path(
        "sandbox",
        if (use.tikz.output) "conplot_talk_figure.tex" else "conplot_talk_figure.pdf"
    ),
    width = 6,
    height = 6,
    ##label.cex = 0.5,
    label.cex = 1,
    label.bg.cex.mult = 4,
    colour.legend = FALSE,
    use.tikz = NULL,
    cex.lab = NULL,
    cex.axis = NULL,
    cex.main = NULL,
    show.diseases = TRUE,
    show.overlays = TRUE,
    show.manual.labels = TRUE
)

load_burnout_package <- function(package.root) {
    if (requireNamespace("pkgload", quietly = TRUE)) {
        message("Loading local development burnout with pkgload::load_all(\".\").")
        pkgload::load_all(package.root, quiet = TRUE)
    } else {
        message("pkgload is not available; falling back to library(burnout).")
        suppressPackageStartupMessages(library(burnout))
    }

    if (!exists("plot_conplot_grid", envir = asNamespace("burnout"), inherits = FALSE)) {
        stop(
            "The loaded burnout package does not export plot_conplot_grid(). ",
            "Install or load the current development checkout.",
            call. = FALSE
        )
    }
}

find_prob_file <- function(package.root) {
    script.arg <- commandArgs(trailingOnly = FALSE)
    script.arg <- script.arg[grepl("^--file=", script.arg)]
    script.file <- if (length(script.arg)) {
        sub("^--file=", "", script.arg[[1L]])
    } else {
        file.path(package.root, "sandbox", "tracked", "conplot_talk_figure.R")
    }
    script.dir <- dirname(normalizePath(script.file, mustWork = FALSE))

    prob.candidates <- unique(c(
        file.path(package.root, "../sources/conplot_standalone/prob.RData"),
        file.path(package.root, "../../sources/conplot_standalone/prob.RData"),
        file.path(script.dir, "../../../sources/conplot_standalone/prob.RData"),
        file.path(script.dir, "../../sources/conplot_standalone/prob.RData")
    ))

    exists <- file.exists(prob.candidates)
    if (!any(exists)) {
        stop(
            "Could not find prob.RData. Tried:\n",
            paste("  -", prob.candidates, collapse = "\n"),
            call. = FALSE
        )
    }

    normalizePath(prob.candidates[exists][[1L]], mustWork = TRUE)
}

load_prob_grid <- function(prob.file) {
    grid.env <- new.env(parent = emptyenv())
    load(prob.file, envir = grid.env)

    required.names <- c("epsvals", "Rvals", "prob")
    missing.names <- setdiff(required.names, ls(grid.env, all.names = TRUE))
    if (length(missing.names)) {
        stop(
            "prob.RData is missing required object(s): ",
            paste(missing.names, collapse = ", "),
            call. = FALSE
        )
    }

    epsvals <- get("epsvals", envir = grid.env)
    Rvals <- get("Rvals", envir = grid.env)
    prob <- get("prob", envir = grid.env)

    if (!is.numeric(epsvals) || !is.vector(epsvals)) {
        stop("epsvals must be a numeric vector.", call. = FALSE)
    }
    if (!is.numeric(Rvals) || !is.vector(Rvals)) {
        stop("Rvals must be a numeric vector.", call. = FALSE)
    }
    if (!is.matrix(prob) || !is.numeric(prob)) {
        stop("prob must be a numeric matrix.", call. = FALSE)
    }

    expected.dim <- c(length(epsvals), length(Rvals))
    if (!identical(dim(prob), expected.dim)) {
        stop(
            "prob dimensions must be length(epsvals) by length(Rvals); ",
            "expected ", paste(expected.dim, collapse = " x "),
            " but got ", paste(dim(prob), collapse = " x "), ".",
            call. = FALSE
        )
    }

    list(epsvals = epsvals, Rvals = Rvals, prob = prob)
}

open_conplot_device <- function(settings, use.tikz.output) {
    output.file <- settings$output.file
    dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)

    if (use.tikz.output) {
        earnmisc::tikz_open(
            file = output.file,
            width = settings$width,
            height = settings$height,
            standAlone = TRUE
        )
    } else {
        grDevices::pdf(
            output.file,
            width = settings$width,
            height = settings$height,
            onefile = TRUE
        )
    }
}

make_conplot_figure <- function(grid, settings, use.tikz.output) {
    open_conplot_device(settings, use.tikz.output)
    device.open <- TRUE
    on.exit({
        if (device.open) {
            grDevices::dev.off()
        }
    }, add = TRUE)

    burnout::plot_conplot_grid(
        epsilon = grid$epsvals,
        R0 = grid$Rvals,
        prob = grid$prob,
        label.cex = settings$label.cex,
        label.bg.cex.mult = settings$label.bg.cex.mult,
        colour.legend = settings$colour.legend,
        use.tikz = settings$use.tikz,
        cex.lab = settings$cex.lab,
        cex.axis = settings$cex.axis,
        cex.main = settings$cex.main,
        show.diseases = settings$show.diseases,
        show.overlays = settings$show.overlays,
        show.manual.labels = settings$show.manual.labels
    )

    grDevices::dev.off()
    device.open <- FALSE
}

package.root <- normalizePath(".", mustWork = TRUE)
description.file <- file.path(package.root, "DESCRIPTION")
if (!file.exists(description.file)) {
    stop("Run this script from the burnout package root.", call. = FALSE)
}

package.name <- read.dcf(description.file)[1L, "Package"]
if (!identical(unname(package.name), "burnout")) {
    stop("Run this script from the burnout package root.", call. = FALSE)
}

load_burnout_package(package.root)

prob.file <- find_prob_file(package.root)
grid <- load_prob_grid(prob.file)

message(
    "Loaded ", basename(prob.file), " with prob dimensions ",
    paste(dim(grid$prob), collapse = " x "), "."
)

plot.settings$output.file <- file.path(package.root, plot.settings$output.file)
make_conplot_figure(grid, plot.settings, use.tikz.output)

message("Wrote ", normalizePath(plot.settings$output.file, mustWork = TRUE))
