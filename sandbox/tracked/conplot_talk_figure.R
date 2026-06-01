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
        sirs = "prob_sirs.RData",
        stop("conplot.model must be \"sir\" or \"sirs\".", call. = FALSE)
    )
}

plot.settings <- list(
    output.file = NULL,
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
    show.manual.labels = TRUE,
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

apply_conplot_model_defaults <- function(settings, conplot.model) {
    if (!conplot.model %in% c("sir", "sirs")) {
        stop("conplot.model must be \"sir\" or \"sirs\".", call. = FALSE)
    }
    if (identical(conplot.model, "sirs")) {
        settings$show.overlays <- FALSE
        settings$show.local.minimum.label <- FALSE
    }
    settings
}

plot.settings <- apply_conplot_model_defaults(plot.settings, conplot.model)

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

open_conplot_device <- function(settings, output.mode) {
    output.mode <- resolve_output_mode(output.mode)
    if (identical(output.mode, "device")) {
        return(FALSE)
    }

    output.file <- settings$output.file
    if (is.null(output.file)) {
        stop("settings$output.file must be set for pdf or tikz output.", call. = FALSE)
    }
    dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)

    if (identical(output.mode, "tikz")) {
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
    TRUE
}

make_conplot_figure <- function(grid, settings, output.mode) {
    output.mode <- resolve_output_mode(output.mode)
    device.open <- open_conplot_device(settings, output.mode)
    on.exit({
        if (device.open) {
            grDevices::dev.off()
        }
    }, add = TRUE)

    plot.info <- burnout::plot_conplot_grid(
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
        show.manual.labels = settings$show.manual.labels,
        show.local.minimum.label = settings$show.local.minimum.label,
        local.minimum.label = settings$local.minimum.label,
        local.minimum.label.cex = settings$local.minimum.label.cex,
        local.minimum.label.col = settings$local.minimum.label.col,
        local.minimum.label.position = settings$local.minimum.label.position,
        local.minimum.label.offset.lines = settings$local.minimum.label.offset.lines,
        show.n.legend = settings$show.n.legend,
        n.legend.label = settings$n.legend.label,
        n.legend.position = settings$n.legend.position,
        n.legend.cex = settings$n.legend.cex,
        n.legend.col = settings$n.legend.col,
        n.legend.bty = settings$n.legend.bty
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

    plot.settings$output.file <- resolve_output_file(
        package.root,
        output.mode,
        plot.settings$output.file
    )
    plot.info <- make_conplot_figure(grid, plot.settings, output.mode)

    if (identical(output.mode, "device")) {
        message("Drew conplot on the active graphics device.")
    } else {
        message("Wrote ", normalizePath(plot.settings$output.file, mustWork = TRUE))
    }

    invisible(plot.info)
}

if (!isTRUE(getOption("burnout.conplot_talk_figure.skip_main", FALSE))) {
    plot.info <- main(output.mode)
}
