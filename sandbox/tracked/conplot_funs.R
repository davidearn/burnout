## Shared helpers for the burnout conplot talk scripts.
##
## Source this file from scripts run at the package root. The helpers load a
## preserved probability grid and pass it to burnout::plot_conplot_grid().
## They do not compute probability grids.

resolve_output_mode <- function(output.mode) {
    if (!is.character(output.mode) || length(output.mode) != 1L || is.na(output.mode)) {
        stop("output.mode must be one of \"device\", \"pdf\", or \"tikz\".", call. = FALSE)
    }
    if (!output.mode %in% c("device", "pdf", "tikz")) {
        stop("output.mode must be one of \"device\", \"pdf\", or \"tikz\".", call. = FALSE)
    }
    output.mode
}

is_absolute_path <- function(path) {
    grepl("^(/|[A-Za-z]:[/\\\\])", path)
}

resolve_output_file <- function(package.root,
                                output.mode,
                                output.file = NULL,
                                output.basename = "conplot_talk_figure") {
    output.mode <- resolve_output_mode(output.mode)
    if (identical(output.mode, "device")) {
        return(NULL)
    }
    if (!is.character(output.basename) || length(output.basename) != 1L ||
        is.na(output.basename) || !nzchar(output.basename)) {
        stop("output.basename must be a non-empty character scalar.", call. = FALSE)
    }

    if (is.null(output.file)) {
        output.file <- file.path(
            "sandbox",
            paste0(output.basename, if (identical(output.mode, "tikz")) ".tex" else ".pdf")
        )
    }
    if (!is.character(output.file) || length(output.file) != 1L ||
        is.na(output.file) || !nzchar(output.file)) {
        stop("output.file must be NULL or a non-empty character scalar.", call. = FALSE)
    }
    if (!is_absolute_path(output.file)) {
        output.file <- file.path(package.root, output.file)
    }
    normalizePath(output.file, mustWork = FALSE)
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
        file.path(package.root, "sandbox", "tracked", "conplot_funs.R")
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
            basename(prob.file), " is missing required object(s): ",
            paste(missing.names, collapse = ", "),
            call. = FALSE
        )
    }

    epsvals <- get("epsvals", envir = grid.env)
    Rvals <- get("Rvals", envir = grid.env)
    prob <- get("prob", envir = grid.env)
    eta <- if (exists("eta", envir = grid.env, inherits = FALSE)) {
        get("eta", envir = grid.env)
    } else {
        NULL
    }

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

    list(
        epsvals = epsvals,
        Rvals = Rvals,
        prob = prob,
        eta = eta,
        prob.file = prob.file
    )
}

resolve_device_scalar <- function(device.settings, name, default) {
    value <- device.settings[[name]]
    if (is.null(value)) {
        return(default)
    }
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
        stop("device.settings$", name, " must be a positive finite scalar.", call. = FALSE)
    }
    value
}

open_conplot_device <- function(device.settings, output.mode) {
    output.mode <- resolve_output_mode(output.mode)
    if (identical(output.mode, "device")) {
        return(list(opened = FALSE, tikz.info = NULL))
    }

    output.file <- device.settings$output.file
    if (is.null(output.file)) {
        stop("device.settings$output.file must be set for pdf or tikz output.", call. = FALSE)
    }
    dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)

    width <- resolve_device_scalar(device.settings, "width", 6)
    height <- resolve_device_scalar(device.settings, "height", 6)

    if (identical(output.mode, "tikz")) {
        tikz.info <- earnmisc::tikz_open(
            file = output.file,
            width = width,
            height = height,
            standAlone = TRUE
        )
        return(list(opened = TRUE, tikz.info = tikz.info))
    }

    grDevices::pdf(
        output.file,
        width = width,
        height = height,
        onefile = TRUE
    )
    list(opened = TRUE, tikz.info = NULL)
}

compile_conplot_tikz <- function(device.info, output.mode, compile.tikz) {
    if (!identical(output.mode, "tikz") || !isTRUE(compile.tikz)) {
        return(invisible(NULL))
    }
    if (is.null(device.info$tikz.info)) {
        stop("Internal error: tikz output requested without tikz device metadata.", call. = FALSE)
    }
    earnmisc::tikz_compile(device.info$tikz.info)
    invisible(NULL)
}

make_conplot_figure <- function(grid,
                                device.settings,
                                plot.args,
                                output.mode,
                                compile.tikz = TRUE) {
    output.mode <- resolve_output_mode(output.mode)
    device.info <- open_conplot_device(device.settings, output.mode)
    device.open <- isTRUE(device.info$opened)
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
    compile_conplot_tikz(device.info, output.mode, compile.tikz)

    invisible(plot.info)
}

run_conplot_script <- function(conplot.model,
                               prob.file.name,
                               output.basename,
                               device.settings,
                               plot.args,
                               output.mode,
                               compile.tikz = TRUE) {
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

    if (!conplot.model %in% c("sir", "sirs")) {
        stop("conplot.model must be \"sir\" or \"sirs\".", call. = FALSE)
    }

    load_burnout_package(package.root)

    prob.file <- find_prob_file(package.root, prob.file.name)
    grid <- load_prob_grid(prob.file)

    eta.text <- if (is.null(grid$eta)) "" else paste0("; eta = ", paste(grid$eta, collapse = ", "))
    message(
        "Loaded ", basename(prob.file), " for ", conplot.model,
        " conplot with prob dimensions ",
        paste(dim(grid$prob), collapse = " x "), eta.text, "."
    )

    device.settings$output.file <- resolve_output_file(
        package.root = package.root,
        output.mode = output.mode,
        output.file = device.settings$output.file,
        output.basename = output.basename
    )

    plot.info <- make_conplot_figure(
        grid = grid,
        device.settings = device.settings,
        plot.args = plot.args,
        output.mode = output.mode,
        compile.tikz = compile.tikz
    )

    if (identical(output.mode, "device")) {
        message("Drew conplot on the active graphics device.")
    } else {
        message("Wrote ", normalizePath(device.settings$output.file, mustWork = TRUE))
    }

    invisible(plot.info)
}
