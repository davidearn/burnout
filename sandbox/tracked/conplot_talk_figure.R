## Recreate the manuscript-style conplot from the preserved probability grid.
##
## Run from the package root:
##
##     Rscript sandbox/tracked/conplot_talk_figure.R
##
## Input:
##
##     ../sources/conplot_standalone/prob_sir.RData
##     ../sources/conplot_standalone/prob_sirs.RData
##     ../sources/conplot_standalone/prob.RData.eta=0.01
##
## In an interactive RStudio session, the default draws to the active graphics
## device. With Rscript, the default writes:
##
##     sandbox/conplot_talk_figure.pdf
##
## Set output.mode to "tikz" to write:
##
##     sandbox/conplot_talk_figure.tex
##
## Supported output.mode values are "device", "pdf", and "tikz". The tikz
## path uses earnmisc::tikz_open(), so plot labels are prepared automatically
## for tikz by earnmisc::nice_text().

## Slide/talk settings. The defaults below keep the manuscript-style appearance.
## Increase label.cex or set colour.legend to TRUE for presentation variants.
output.mode <- if (interactive()) "device" else "pdf"
##output.mode <- "tikz"

conplot.model <- Sys.getenv("BURNOUT_CONPLOT_MODEL", unset = "sirs")
##conplot.model <- "sir"
##conplot.model <- "sirs"

prob.file.name <- Sys.getenv("BURNOUT_CONPLOT_FILE", unset = "")
if (!nzchar(prob.file.name)) {
    prob.file.name <- switch(
        conplot.model,
        sir = "prob_sir.RData",
        ##sirs = "prob_sirs.RData",
        sirs = "prob.RData.eta=0.01",
        stop("conplot.model must be \"sir\" or \"sirs\".", call. = FALSE)
    )
}

device.settings <- list(
    output.file = NULL,
    width = 6,
    height = 6
)

plot.args <- list(
    label.cex = 1,
    label.bg.cex.mult = 4,
    colour.legend = FALSE,
    use.tikz = NULL,
    cex.lab = NULL,
    cex.axis = NULL,
    cex.main = NULL,
    show.diseases = TRUE,
    show.manual.labels = TRUE,
    show.quadratic = TRUE,
    show.local.minimum = TRUE,
    show.local.minimum.label = TRUE,
    local.minimum.label = "minimum persistence probability",
    local.minimum.label.cex = NULL,
    local.minimum.label.col = NULL,
    local.minimum.label.position = NULL,
    local.minimum.label.offset.lines = 1,
    show.n.legend = TRUE,
    n.legend.label = "$n = 10^6$",
    n.legend.position = "topright",
    n.legend.cex = NULL,
    ##n.legend.col = "black",
    n.legend.col = "grey95",
    n.legend.bty = "n"
)

apply_conplot_model_defaults <- function(plot.args, conplot.model) {
    if (!conplot.model %in% c("sir", "sirs")) {
        stop("conplot.model must be \"sir\" or \"sirs\".", call. = FALSE)
    }
    if (identical(conplot.model, "sirs")) {
        plot.args$show.quadratic <- FALSE
        plot.args$show.local.minimum <- TRUE
        plot.args$show.local.minimum.label <- FALSE
    }
    plot.args
}

plot.args <- apply_conplot_model_defaults(plot.args, conplot.model)

resolve_output_mode <- function(output.mode) {
    if (!is.character(output.mode) || length(output.mode) != 1L || is.na(output.mode)) {
        stop("output.mode must be one of \"device\", \"pdf\", or \"tikz\".", call. = FALSE)
    }
    if (!output.mode %in% c("device", "pdf", "tikz")) {
        stop("output.mode must be one of \"device\", \"pdf\", or \"tikz\".", call. = FALSE)
    }
    output.mode
}

resolve_output_file <- function(package.root, output.mode, output.file = NULL) {
    output.mode <- resolve_output_mode(output.mode)
    if (identical(output.mode, "device")) {
        return(NULL)
    }
    if (is.null(output.file)) {
        output.file <- file.path(
            "sandbox",
            if (identical(output.mode, "tikz")) {
                "conplot_talk_figure.tex"
            } else {
                "conplot_talk_figure.pdf"
            }
        )
    }
    normalizePath(file.path(package.root, output.file), mustWork = FALSE)
}

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

find_prob_file <- function(package.root, prob.file.name) {
    if (!is.character(prob.file.name) || length(prob.file.name) != 1L ||
        is.na(prob.file.name) || !nzchar(prob.file.name)) {
        stop("prob.file.name must be a non-empty character scalar.", call. = FALSE)
    }

    script.arg <- commandArgs(trailingOnly = FALSE)
    script.arg <- script.arg[grepl("^--file=", script.arg)]
    script.file <- if (length(script.arg)) {
        sub("^--file=", "", script.arg[[1L]])
    } else {
        file.path(package.root, "sandbox", "tracked", "conplot_talk_figure.R")
    }
    script.dir <- dirname(normalizePath(script.file, mustWork = FALSE))

    prob.candidates <- unique(c(
        file.path(package.root, "../sources/conplot_standalone", prob.file.name),
        file.path(package.root, "../../sources/conplot_standalone", prob.file.name),
        file.path(script.dir, "../../../sources/conplot_standalone", prob.file.name),
        file.path(script.dir, "../../sources/conplot_standalone", prob.file.name)
    ))

    exists <- file.exists(prob.candidates)
    if (!any(exists)) {
        stop(
            "Could not find ", prob.file.name, ". Tried:\n",
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

open_conplot_device <- function(device.settings, output.mode) {
    output.mode <- resolve_output_mode(output.mode)
    if (identical(output.mode, "device")) {
        return(FALSE)
    }

    output.file <- device.settings$output.file
    if (is.null(output.file)) {
        stop("device.settings$output.file must be set for pdf or tikz output.", call. = FALSE)
    }
    dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)

    if (identical(output.mode, "tikz")) {
        earnmisc::tikz_open(
            file = output.file,
            width = device.settings$width,
            height = device.settings$height,
            standAlone = TRUE
        )
    } else {
        grDevices::pdf(
            output.file,
            width = device.settings$width,
            height = device.settings$height,
            onefile = TRUE
        )
    }
    TRUE
}

make_conplot_figure <- function(grid, device.settings, plot.args, output.mode) {
    output.mode <- resolve_output_mode(output.mode)
    device.open <- open_conplot_device(device.settings, output.mode)
    on.exit({
        if (device.open) {
            grDevices::dev.off()
        }
    }, add = TRUE)

    protected.args <- c("epsilon", "R0", "prob")
    protected.present <- intersect(protected.args, names(plot.args))
    if (length(protected.present)) {
        stop(
            "plot.args must not contain core grid argument(s): ",
            paste(protected.present, collapse = ", "),
            call. = FALSE
        )
    }

    plot.info <- do.call(
        burnout::plot_conplot_grid,
        c(list(
        epsilon = grid$epsvals,
        R0 = grid$Rvals,
        prob = grid$prob
        ), plot.args)
    )

    if (device.open) {
        grDevices::dev.off()
        device.open <- FALSE
    }

    invisible(plot.info)
}

main <- function(output.mode.value) {
    output.mode <- output.mode.value
    output.mode <- resolve_output_mode(output.mode)
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

    prob.file <- find_prob_file(package.root, prob.file.name)
    grid <- load_prob_grid(prob.file)

    message(
        "Loaded ", basename(prob.file), " for ", conplot.model,
        " conplot with prob dimensions ",
        paste(dim(grid$prob), collapse = " x "), "."
    )

    device.settings$output.file <- resolve_output_file(
        package.root,
        output.mode,
        device.settings$output.file
    )
    plot.info <- make_conplot_figure(
        grid = grid,
        device.settings = device.settings,
        plot.args = plot.args,
        output.mode = output.mode
    )

    if (identical(output.mode, "device")) {
        message("Drew conplot on the active graphics device.")
    } else {
        message("Wrote ", normalizePath(device.settings$output.file, mustWork = TRUE))
    }

    invisible(plot.info)
}

if (!isTRUE(getOption("burnout.conplot_talk_figure.skip_main", FALSE))) {
    plot.info <- main(output.mode)
}
