##' Plot a precomputed conplot probability grid
##'
##' Plot a contour-style figure from precomputed burnout or persistence
##' probabilities. This function does not compute probabilities; it only
##' plots an existing grid, such as the kind of grid used for the historical
##' conplot figure. Future helpers may compute compatible grids directly.
##'
##' The probability matrix must have rows corresponding to `epsilon` values
##' and columns corresponding to `R0` values. This is the orientation used by
##' the historical `prob.RData` reference grid.
##'
##' @param epsilon finite numeric vector of epsilon values for the vertical
##'   axis. Values must be strictly increasing.
##' @param R0 finite numeric vector of reproduction-number values for the
##'   horizontal axis. Values must be strictly increasing.
##' @param prob numeric probability matrix with `length(epsilon)` rows and
##'   `length(R0)` columns. Values must be in `[0, 1]`; `NA` values are
##'   allowed and are left blank in the plot.
##' @param label.cex positive finite scalar controlling contour-label size.
##' @param colour.legend logical scalar. If `TRUE`, add a right-hand colour
##'   ramp legend for the filled background.
##' @param filled logical scalar. If `TRUE`, draw a filled image-like
##'   background before contour lines.
##' @param levels numeric vector of contour levels. If `NULL`, the union of
##'   `high.levels` and `low.levels` is used.
##' @param high.levels,low.levels numeric vectors used to construct contour
##'   levels when `levels = NULL`.
##' @param log character string passed to base graphics for logarithmic axes.
##'   The default `"y"` uses a logarithmic epsilon axis.
##' @param xlim,ylim finite numeric vectors of length two giving plot limits.
##' @param xlab,ylab,main plot labels passed to `plot`.
##' @param fill.colours optional vector of colours for the filled background.
##'   If `NULL`, a default slide-friendly colour ramp is used.
##' @param contour.col,contour.lwd colour and line width for contour lines.
##' @param ... additional arguments passed to `plot`.
##'
##' @return Invisibly returns a list with plotting metadata, including the
##'   contour levels, fill breaks, axis ranges, and probability-matrix
##'   orientation.
##'
##' @export
##'
##' @examples
##' epsilon <- 10^seq(-4, -1, length.out = 20)
##' R0 <- seq(1.05, 4, length.out = 25)
##' prob <- outer(
##'     epsilon, R0,
##'     function(epsilon, R0) plogis(2 * log(R0) + log10(epsilon) + 3)
##' )
##' plot_conplot_grid(epsilon, R0, prob, label.cex = 0.7)
plot_conplot_grid <- function(epsilon,
                              R0,
                              prob,
                              label.cex = 0.8,
                              colour.legend = FALSE,
                              filled = TRUE,
                              levels = NULL,
                              high.levels = c(seq(0.1, 0.9, by = 0.1), 0.95),
                              low.levels = 10^(-c(1, 2, 4, 8, 12)),
                              log = "y",
                              xlim = range(R0),
                              ylim = range(epsilon),
                              xlab = expression(R[0]),
                              ylab = expression(epsilon),
                              main = "Persistence probability",
                              fill.colours = NULL,
                              contour.col = "black",
                              contour.lwd = 1,
                              ...) {

    epsilon <- validate_conplot_grid_vector(epsilon, "epsilon")
    R0 <- validate_conplot_grid_vector(R0, "R0")
    prob <- validate_conplot_prob_matrix(prob, epsilon, R0)
    label.cex <- validate_conplot_positive_scalar(label.cex, "label.cex")
    colour.legend <- validate_conplot_logical_scalar(colour.legend, "colour.legend")
    filled <- validate_conplot_logical_scalar(filled, "filled")
    log <- validate_conplot_log(log, epsilon = epsilon, R0 = R0)
    xlim <- validate_conplot_range(xlim, "xlim")
    ylim <- validate_conplot_range(ylim, "ylim")

    if (is.null(levels)) {
        high.levels <- validate_conplot_levels(high.levels, "high.levels")
        low.levels <- validate_conplot_levels(low.levels, "low.levels")
        levels <- sort(unique(c(high.levels, low.levels)))
    } else {
        levels <- validate_conplot_levels(levels, "levels")
        high.levels <- validate_conplot_levels(high.levels, "high.levels")
        low.levels <- validate_conplot_levels(low.levels, "low.levels")
    }

    fill.breaks <- sort(unique(c(0, levels, 1)))
    if (length(fill.breaks) < 2L) {
        stop("fill breaks must contain at least two distinct values.", call. = FALSE)
    }
    if (is.null(fill.colours)) {
        fill.colours <- grDevices::hcl.colors(
            n = length(fill.breaks) - 1L,
            palette = "YlOrRd",
            rev = TRUE
        )
    } else {
        fill.colours <- validate_conplot_colours(
            fill.colours,
            expected.length = length(fill.breaks) - 1L
        )
    }

    z <- t(prob)
    zlim <- range(prob, na.rm = TRUE)

    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)

    if (colour.legend) {
        graphics::par(fig = c(0, 0.82, 0, 1))
    }

    graphics::plot.default(
        x = R0,
        y = rep(ylim[1L], length(R0)),
        type = "n",
        xlim = xlim,
        ylim = ylim,
        log = log,
        xlab = xlab,
        ylab = ylab,
        main = main,
        ...
    )

    if (filled) {
        graphics::image(
            x = R0,
            y = epsilon,
            z = z,
            breaks = fill.breaks,
            col = fill.colours,
            add = TRUE
        )
    }

    graphics::contour(
        x = R0,
        y = epsilon,
        z = z,
        levels = levels,
        add = TRUE,
        labcex = label.cex,
        col = contour.col,
        lwd = contour.lwd
    )
    graphics::box()

    if (colour.legend) {
        draw_conplot_colour_legend(
            fill.breaks = fill.breaks,
            fill.colours = fill.colours
        )
    }

    out <- list(
        levels = levels,
        fill.breaks = fill.breaks,
        xlim = xlim,
        ylim = ylim,
        zlim = zlim,
        log = log,
        filled = filled,
        colour.legend = colour.legend,
        label.cex = label.cex,
        orientation = "rows: epsilon; columns: R0"
    )

    invisible(out)
}

##' Validate a conplot grid vector
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return A numeric vector.
##' @noRd
validate_conplot_grid_vector <- function(x, name) {
    if (!is.numeric(x) || !is.vector(x) || length(x) < 2L) {
        stop(name, " must be a numeric vector with length at least 2.", call. = FALSE)
    }
    if (anyNA(x) || any(!is.finite(x))) {
        stop(name, " must contain only finite values.", call. = FALSE)
    }
    if (any(diff(x) <= 0)) {
        stop(name, " must be strictly increasing.", call. = FALSE)
    }
    x
}

##' Validate a conplot probability matrix
##'
##' @param prob probability matrix.
##' @param epsilon epsilon grid vector.
##' @param R0 reproduction-number grid vector.
##'
##' @return A numeric matrix.
##' @noRd
validate_conplot_prob_matrix <- function(prob, epsilon, R0) {
    if (!is.matrix(prob) || !is.numeric(prob)) {
        stop("prob must be a numeric matrix.", call. = FALSE)
    }
    expected.dim <- c(length(epsilon), length(R0))
    if (!identical(dim(prob), expected.dim)) {
        stop(
            "prob must have dimensions length(epsilon) by length(R0); ",
            "expected ", paste(expected.dim, collapse = " x "),
            " but got ", paste(dim(prob), collapse = " x "), ".",
            call. = FALSE
        )
    }
    if (any(is.nan(prob))) {
        stop("prob may contain NA values but not NaN values.", call. = FALSE)
    }
    observed <- prob[!is.na(prob)]
    if (length(observed) == 0L) {
        stop("prob must contain at least one non-NA value.", call. = FALSE)
    }
    if (any(!is.finite(observed))) {
        stop("all non-NA prob values must be finite.", call. = FALSE)
    }
    if (any(observed < 0 | observed > 1)) {
        stop("all non-NA prob values must be in [0, 1].", call. = FALSE)
    }
    prob
}

##' Validate a positive scalar
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return A numeric scalar.
##' @noRd
validate_conplot_positive_scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        !is.finite(x) || x <= 0) {
        stop(name, " must be a positive finite scalar.", call. = FALSE)
    }
    x
}

##' Validate a logical scalar
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return A logical scalar.
##' @noRd
validate_conplot_logical_scalar <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop(name, " must be a scalar logical value.", call. = FALSE)
    }
    x
}

##' Validate contour levels
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return A sorted numeric vector.
##' @noRd
validate_conplot_levels <- function(x, name) {
    if (!is.numeric(x) || length(x) == 0L || anyNA(x) || any(!is.finite(x))) {
        stop(name, " must be a non-empty finite numeric vector.", call. = FALSE)
    }
    if (any(x < 0 | x > 1)) {
        stop(name, " values must be in [0, 1].", call. = FALSE)
    }
    sort(unique(x))
}

##' Validate plot limits
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return A finite numeric vector of length two.
##' @noRd
validate_conplot_range <- function(x, name) {
    if (!is.numeric(x) || length(x) != 2L || anyNA(x) || any(!is.finite(x))) {
        stop(name, " must be a finite numeric vector of length 2.", call. = FALSE)
    }
    if (x[1L] >= x[2L]) {
        stop(name, " must be increasing.", call. = FALSE)
    }
    x
}

##' Validate axis log setting
##'
##' @param log base-graphics log setting.
##' @param epsilon epsilon grid vector.
##' @param R0 reproduction-number grid vector.
##'
##' @return A character scalar.
##' @noRd
validate_conplot_log <- function(log, epsilon, R0) {
    if (!is.character(log) || length(log) != 1L || is.na(log)) {
        stop("log must be a character scalar.", call. = FALSE)
    }
    if (!log %in% c("", "x", "y", "xy")) {
        stop('log must be one of "", "x", "y", or "xy".', call. = FALSE)
    }
    if (grepl("x", log, fixed = TRUE) && any(R0 <= 0)) {
        stop("R0 values must be positive when log includes the x axis.", call. = FALSE)
    }
    if (grepl("y", log, fixed = TRUE) && any(epsilon <= 0)) {
        stop("epsilon values must be positive when log includes the y axis.", call. = FALSE)
    }
    log
}

##' Validate conplot fill colours
##'
##' @param fill.colours vector of colours.
##' @param expected.length required colour-vector length.
##'
##' @return A character vector of colours.
##' @noRd
validate_conplot_colours <- function(fill.colours, expected.length) {
    if (!is.character(fill.colours) || length(fill.colours) != expected.length ||
        anyNA(fill.colours)) {
        stop(
            "fill.colours must be a character vector with length ",
            expected.length, ".",
            call. = FALSE
        )
    }
    fill.colours
}

##' Draw a right-hand colour legend for a conplot grid
##'
##' @param fill.breaks numeric break points.
##' @param fill.colours colours corresponding to break intervals.
##'
##' @return `NULL`, invisibly.
##' @noRd
draw_conplot_colour_legend <- function(fill.breaks, fill.colours) {
    graphics::par(fig = c(0.84, 1, 0, 1), new = TRUE, mar = c(5.1, 0.5, 4.1, 3.5))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = range(fill.breaks), xaxs = "i", yaxs = "i")

    graphics::rect(
        xleft = 0,
        ybottom = fill.breaks[-length(fill.breaks)],
        xright = 1,
        ytop = fill.breaks[-1L],
        col = fill.colours,
        border = NA
    )
    graphics::box()

    axis.at <- pretty(fill.breaks)
    axis.at <- axis.at[axis.at >= min(fill.breaks) & axis.at <= max(fill.breaks)]
    graphics::axis(side = 4, at = axis.at, las = 1)
    graphics::mtext("probability", side = 4, line = 2.4)

    invisible(NULL)
}
