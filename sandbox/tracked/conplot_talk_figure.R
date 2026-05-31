## Create a slide-friendly conplot figure from the preserved probability grid.
##
## Run from the package root:
##
##     Rscript sandbox/tracked/conplot_talk_figure.R
##
## Input:
##
##     ../sources/conplot_standalone/prob.RData
##
## Output:
##
##     sandbox/conplot_talk_figure.pdf

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

make_conplot_pdf <- function(grid, output.file) {
    dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)

    grDevices::pdf(output.file, width = 11, height = 6.5, onefile = TRUE)
    pdf.open <- TRUE
    on.exit({
        if (pdf.open) {
            grDevices::dev.off()
        }
    }, add = TRUE)

    graphics::par(mar = c(5.2, 5.4, 4.1, 1.2), cex = 1.05)

    burnout::plot_conplot_grid(
        epsilon = grid$epsvals,
        R0 = grid$Rvals,
        prob = grid$prob,
        label.cex = 1.05,
        colour.legend = TRUE,
        filled = TRUE,
        main = "Persistence probability after a major outbreak",
        xlab = expression(R[0]),
        ylab = expression(epsilon),
        contour.lwd = 1.1,
        cex.lab = 1.25,
        cex.axis = 1.05,
        cex.main = 1.15
    )

    grDevices::dev.off()
    pdf.open <- FALSE
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

output.file <- file.path(package.root, "sandbox", "conplot_talk_figure.pdf")
make_conplot_pdf(grid, output.file)

message("Wrote ", normalizePath(output.file, mustWork = TRUE))
