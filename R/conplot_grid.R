##' Plot a precomputed conplot probability grid
##'
##' Plot a contour-style figure from precomputed burnout or persistence
##' probabilities. This function does not compute probabilities; it only
##' plots an existing grid, such as the historical `prob.RData` reference
##' grid. The default graphical choices are designed to reproduce the
##' original conplot manuscript figure closely, while allowing slide-oriented
##' adjustments.
##'
##' The probability matrix must have rows corresponding to `epsilon` values
##' and columns corresponding to `R0` values. By default, `epsilon` is plotted
##' on the horizontal axis and `R0` on a logarithmic vertical axis, matching
##' the original standalone conplot.
##'
##' @param epsilon finite numeric vector of epsilon values. Values must be
##'   strictly increasing.
##' @param R0 finite numeric vector of reproduction-number values. Values
##'   must be strictly increasing.
##' @param prob numeric probability matrix with `length(epsilon)` rows and
##'   `length(R0)` columns. Values must be in `[0, 1]`; `NA` values are
##'   allowed and are left blank in the plot.
##' @param label.cex positive finite scalar controlling the size of manual
##'   contour labels and any contour labels drawn by `contour()`. The default
##'   is `0.7` in interactive sessions and `0.5` otherwise,
##'   preserving manuscript-style scripted output while making RStudio-style
##'   inspection easier.
##' @param label.bg.cex.mult positive finite scalar. Manual contour-label
##'   background circles are drawn with point size
##'   `label.cex * label.bg.cex.mult` unless `manual.label.point.cex` is
##'   supplied explicitly. The default reproduces the original manuscript
##'   circle size when `label.cex = 0.5`.
##' @param use.tikz optional logical scalar passed to
##'   [earnmisc::nice_text()] for TeX-like labels. If `NULL`, a scalar logical
##'   caller-level `use.tikz` variable is honoured when present; otherwise
##'   [earnmisc::nice_text()] detects an active tikz graphics device. This
##'   means tikz and ordinary non-tikz devices both work without manual label
##'   conversion.
##' @param label.warn logical scalar passed to [earnmisc::nice_text()]. If
##'   `TRUE`, warn about unsupported TeX macros or missing optional conversion
##'   support.
##' @param colour.legend logical scalar. If `TRUE`, add a right-hand colour
##'   ramp legend for the filled background.
##' @param filled logical scalar. If `TRUE`, draw the original greyscale
##'   filled background before contour lines.
##' @param levels numeric vector of contour levels. If `NULL`, the union of
##'   `high.levels` and `low.levels` is used.
##' @param high.levels,low.levels numeric vectors used to construct contour
##'   levels when `levels = NULL`. `high.levels` also define the filled
##'   greyscale break points by default.
##' @param log character string passed to base graphics for logarithmic axes.
##'   The default `"y"` uses the original logarithmic `R0` axis.
##' @param xlim,ylim finite numeric vectors of length two giving plot limits.
##'   Defaults match the original manuscript figure.
##' @param xlab,ylab,main plot labels. Character labels are processed with
##'   [earnmisc::nice_text()]. The default epsilon label uses
##'   `"$\\varepsilon$"` and the default reproduction-number label uses the
##'   `"$\\Rn$"` macro from `earnmisc`.
##' @param fill.colours optional vector of colours for the filled background.
##'   If `NULL`, the original high-level greyscale ramp is used.
##' @param fill.breaks optional vector of filled-background break points. If
##'   `NULL`, `c(0, high.levels, 1)` is used, as in the original figure.
##' @param grey.adjust numeric scalar controlling the original greyscale ramp
##'   when `fill.colours = NULL`.
##' @param contour.col,contour.lwd,contour.lty colour, line width, and line
##'   type for contour lines.
##' @param contour.labels labels passed to `contour()`. The default `""`
##'   suppresses automatic labels, preserving the original manual label
##'   placement.
##' @param show.manual.labels logical scalar. If `TRUE`, add the original
##'   hand-positioned high- and low-level labels. Labels whose profiles
##'   contain too few finite values, for example because an edge row or
##'   column is all `NA`, are omitted rather than causing the plot to fail.
##' @param manual.label.x epsilon value for high-level labels. If `NULL`, use
##'   `xlim[2] - 0.002`, matching the standalone code.
##' @param manual.label.y `R0` value for low-level labels.
##' @param manual.label.point.cex optional positive finite scalar for the
##'   filled points behind manual labels. If `NULL`, it is computed from
##'   `label.cex * label.bg.cex.mult`.
##' @param show.overlays logical scalar controlling both original dark-red
##'   overlay curves unless `show.quadratic` or `show.local.minimum` is
##'   supplied explicitly.
##' @param show.quadratic logical scalar. If `TRUE`, draw the original dotted
##'   quadratic curve.
##' @param quadratic.coefficients named or unnamed numeric vector with
##'   intercept, slope, and curvature for the quadratic overlay.
##' @param quadratic.from,quadratic.to finite scalars giving the horizontal
##'   range for the quadratic overlay.
##' @param quadratic.col,quadratic.lty,quadratic.lwd graphical parameters for
##'   the quadratic overlay.
##' @param show.local.minimum logical scalar. If `TRUE`, draw the original
##'   local-minimum curve estimated directly from the precomputed grid.
##' @param local.minimum.xlow,local.minimum.ylow,local.minimum.yhigh numeric
##'   limits used when estimating the local-minimum curve.
##' @param local.minimum.col,local.minimum.lwd graphical parameters for the
##'   local-minimum curve.
##' @param show.local.minimum.label logical scalar. If `TRUE`, label the
##'   local-minimum curve. The default follows `show.local.minimum`.
##' @param local.minimum.label character label for the local-minimum curve.
##'   The label is processed with [earnmisc::nice_text()].
##' @param local.minimum.label.col optional label colour. If `NULL`, the
##'   label uses `local.minimum.col`.
##' @param local.minimum.label.cex optional positive finite scalar for the
##'   local-minimum label. If `NULL`, the value of `label.cex` is used.
##' @param local.minimum.label.position optional finite numeric vector of
##'   length two giving the epsilon and `R0` coordinates used as the
##'   local-minimum label anchor before applying the vertical offset. If
##'   `NULL`, the midpoint of the quadratic overlay range is used with the
##'   corresponding quadratic `R0` value.
##' @param local.minimum.label.offset.lines finite numeric scalar giving the
##'   vertical label offset in text-line heights. Positive values place the
##'   label visually above the curve; the default avoids overlap with the
##'   dark-red curves on the logarithmic `R0` axis.
##' @param show.n.legend logical scalar. If `TRUE`, draw the original
##'   top-right sample-size annotation.
##' @param n.legend.label character label for the sample-size annotation.
##'   The default is `"$n = 10^6$"` and is processed with
##'   [earnmisc::nice_text()].
##' @param n.legend.position character legend position such as `"topright"`,
##'   or a finite numeric vector of length two giving coordinates. The
##'   default preserves the original top-right placement.
##' @param n.legend.cex optional positive finite scalar for the sample-size
##'   annotation. If `NULL`, the base legend default `1` is used.
##' @param n.legend.col,n.legend.bty colour and box type passed to
##'   `legend()` for the sample-size annotation.
##' @param show.diseases logical scalar. If `TRUE`, add the original disease
##'   point annotations when `disease.data` is supplied or the suggested
##'   `sirr` package is available.
##' @param disease.data optional data frame with at least `epsilon`, `R0`, and
##'   `label` columns. If `NULL`, compatible data are loaded from
##'   `sirr::diseaseParameters` when available.
##' @param disease.names character vector of disease identifiers to show when
##'   using `sirr::diseaseParameters`.
##' @param disease.cex,disease.point.cex positive finite scalars controlling
##'   disease-label and disease-point size.
##' @param disease.point.bg,disease.box.col colours for disease points and the
##'   white-out boxes behind labels.
##' @param cex.lab,cex.axis,cex.main optional positive finite scalars passed to
##'   base graphics for label, axis, and title sizes.
##' @param ... additional arguments passed to `plot()`.
##'
##' @return Invisibly returns a list with plotting metadata, including
##'   contour levels, fill breaks, axis ranges, overlay settings, annotation
##'   settings, probability-matrix orientation, and resolved label metadata.
##'   In particular, `local.minimum.label$position` contains the actual
##'   offset label coordinates used for plotting, and
##'   `local.minimum.label$curve.position` contains the unoffset curve-anchor
##'   coordinates.
##'
##' @export
##'
##' @examples
##' epsilon <- 10^seq(-5, -1.7, length.out = 30)
##' R0 <- 2^seq(0, 5, length.out = 35)
##' prob <- outer(
##'     epsilon, R0,
##'     function(epsilon, R0) plogis(3 * log(R0) + 180 * epsilon - 8)
##' )
##' plot_conplot_grid(
##'     epsilon, R0, prob,
##'     show.diseases = FALSE,
##'     show.local.minimum = FALSE
##' )
plot_conplot_grid <- function(epsilon,
                              R0,
                              prob,
                              label.cex = conplot_default_label_cex(),
                              label.bg.cex.mult = 4,
                              use.tikz = NULL,
                              label.warn = TRUE,
                              colour.legend = FALSE,
                              filled = TRUE,
                              levels = NULL,
                              high.levels = c(seq(0.1, 0.9, by = 0.1), 0.95),
                              low.levels = 10^(-c(1, 2, 4, 8, 12)),
                              log = "y",
                              xlim = c(0, 0.02),
                              ylim = c(1, 32),
                              xlab = "mean infectious period / mean lifetime ($\\varepsilon$)",
                              ylab = "basic reproduction number ($\\Rn$)",
                              main = "contours of persistence probability",
                              fill.colours = NULL,
                              fill.breaks = NULL,
                              grey.adjust = 0.7,
                              contour.col = "black",
                              contour.lwd = 2,
                              contour.lty = 1,
                              contour.labels = "",
                              show.manual.labels = TRUE,
                              manual.label.x = NULL,
                              manual.label.y = 2,
                              manual.label.point.cex = NULL,
                              show.overlays = TRUE,
                              show.quadratic = show.overlays,
                              quadratic.coefficients = c(
                                  intercept = 2.572629848,
                                  slope = -27.71866282,
                                  curvature = -109.8
                              ),
                              quadratic.from = xlim[1L],
                              quadratic.to = xlim[2L],
                              quadratic.col = "darkred",
                              quadratic.lty = "dotted",
                              quadratic.lwd = 1,
                              show.local.minimum = show.overlays,
                              local.minimum.xlow = 0.004,
                              local.minimum.ylow = 1.5,
                              local.minimum.yhigh = 5,
                              local.minimum.col = "darkred",
                              local.minimum.lwd = 3,
                              show.local.minimum.label = show.local.minimum,
                              local.minimum.label = "minimum persistence probability",
                              local.minimum.label.col = NULL,
                              local.minimum.label.cex = NULL,
                              local.minimum.label.position = NULL,
                              local.minimum.label.offset.lines = 1,
                              show.n.legend = TRUE,
                              n.legend.label = "$n = 10^6$",
                              n.legend.position = "topright",
                              n.legend.cex = NULL,
                              n.legend.col = "black",
                              n.legend.bty = "n",
                              show.diseases = TRUE,
                              disease.data = NULL,
                              disease.names = c(
                                  "measles", "pertussis", "mumps", "rubella",
                                  "chickenpox", "influenza", "smallpox",
                                  "scarletfever", "HIV", "covid19",
                                  "covid19delta", "ebola", "pneumonicplague"
                              ),
                              disease.cex = 0.9,
                              disease.point.cex = 1,
                              disease.point.bg = "darkred",
                              disease.box.col = "white",
                              cex.lab = NULL,
                              cex.axis = NULL,
                              cex.main = NULL,
                              ...) {

    epsilon <- validate_conplot_grid_vector(epsilon, "epsilon")
    R0 <- validate_conplot_grid_vector(R0, "R0")
    prob <- validate_conplot_prob_matrix(prob, epsilon, R0)
    label.cex <- validate_conplot_positive_scalar(label.cex, "label.cex")
    label.bg.cex.mult <- validate_conplot_positive_scalar(
        label.bg.cex.mult, "label.bg.cex.mult"
    )
    label.warn <- validate_conplot_logical_scalar(label.warn, "label.warn")
    use.tikz <- resolve_conplot_use_tikz(
        use.tikz = use.tikz,
        envir = parent.frame(),
        warn = label.warn
    )
    colour.legend <- validate_conplot_logical_scalar(colour.legend, "colour.legend")
    filled <- validate_conplot_logical_scalar(filled, "filled")
    show.manual.labels <- validate_conplot_logical_scalar(
        show.manual.labels, "show.manual.labels"
    )
    show.overlays <- validate_conplot_logical_scalar(show.overlays, "show.overlays")
    show.quadratic <- validate_conplot_logical_scalar(show.quadratic, "show.quadratic")
    show.local.minimum <- validate_conplot_logical_scalar(
        show.local.minimum, "show.local.minimum"
    )
    show.local.minimum.label <- validate_conplot_logical_scalar(
        show.local.minimum.label, "show.local.minimum.label"
    )
    if (show.local.minimum.label && !show.local.minimum) {
        stop(
            "show.local.minimum.label = TRUE requires show.local.minimum = TRUE.",
            call. = FALSE
        )
    }
    show.n.legend <- validate_conplot_logical_scalar(show.n.legend, "show.n.legend")
    show.diseases <- validate_conplot_logical_scalar(show.diseases, "show.diseases")
    log <- validate_conplot_log(log, epsilon = epsilon, R0 = R0)
    xlim <- validate_conplot_range(xlim, "xlim")
    ylim <- validate_conplot_range(ylim, "ylim")

    high.levels <- validate_conplot_levels(high.levels, "high.levels")
    low.levels <- validate_conplot_levels(low.levels, "low.levels")
    if (is.null(levels)) {
        levels <- sort(unique(c(low.levels, high.levels)))
    } else {
        levels <- validate_conplot_levels(levels, "levels")
    }

    if (is.null(fill.breaks)) {
        fill.breaks <- sort(unique(c(0, high.levels, 1)))
    } else {
        fill.breaks <- validate_conplot_levels(fill.breaks, "fill.breaks")
    }
    if (length(fill.breaks) < 2L) {
        stop("fill.breaks must contain at least two distinct values.", call. = FALSE)
    }
    if (is.unsorted(fill.breaks, strictly = TRUE)) {
        stop("fill.breaks must be strictly increasing.", call. = FALSE)
    }

    grey.adjust <- validate_conplot_probability_scalar(grey.adjust, "grey.adjust")
    if (is.null(fill.colours)) {
        fill.colours <- rev(grDevices::grey(1 - grey.adjust + fill.breaks[-1L] * grey.adjust))
    } else {
        fill.colours <- validate_conplot_colours(
            fill.colours,
            expected.length = length(fill.breaks) - 1L
        )
    }

    contour.lwd <- validate_conplot_positive_scalar(contour.lwd, "contour.lwd")
    if (is.null(manual.label.point.cex)) {
        manual.label.point.cex <- label.cex * label.bg.cex.mult
    } else {
        manual.label.point.cex <- validate_conplot_positive_scalar(
            manual.label.point.cex, "manual.label.point.cex"
        )
    }
    manual.label.y <- validate_conplot_finite_scalar(manual.label.y, "manual.label.y")
    if (is.null(manual.label.x)) {
        manual.label.x <- xlim[2L] - 0.002
    }
    manual.label.x <- validate_conplot_finite_scalar(manual.label.x, "manual.label.x")

    quadratic.coefficients <- validate_conplot_quadratic_coefficients(
        quadratic.coefficients
    )
    quadratic.from <- validate_conplot_finite_scalar(quadratic.from, "quadratic.from")
    quadratic.to <- validate_conplot_finite_scalar(quadratic.to, "quadratic.to")
    if (quadratic.from >= quadratic.to) {
        stop("quadratic.from must be less than quadratic.to.", call. = FALSE)
    }
    quadratic.lwd <- validate_conplot_positive_scalar(quadratic.lwd, "quadratic.lwd")

    local.minimum.xlow <- validate_conplot_finite_scalar(
        local.minimum.xlow, "local.minimum.xlow"
    )
    local.minimum.ylow <- validate_conplot_finite_scalar(
        local.minimum.ylow, "local.minimum.ylow"
    )
    local.minimum.yhigh <- validate_conplot_finite_scalar(
        local.minimum.yhigh, "local.minimum.yhigh"
    )
    local.minimum.lwd <- validate_conplot_positive_scalar(
        local.minimum.lwd, "local.minimum.lwd"
    )
    local.minimum.label.col <- validate_conplot_optional_colour_scalar(
        local.minimum.label.col, "local.minimum.label.col"
    )
    if (is.null(local.minimum.label.col)) {
        local.minimum.label.col <- local.minimum.col
    }
    local.minimum.label.cex <- validate_conplot_optional_positive_scalar(
        local.minimum.label.cex, "local.minimum.label.cex"
    )
    if (is.null(local.minimum.label.cex)) {
        local.minimum.label.cex <- label.cex
    }
    local.minimum.label.position <- validate_conplot_optional_position(
        local.minimum.label.position,
        "local.minimum.label.position"
    )
    local.minimum.label.offset.lines <- validate_conplot_finite_scalar(
        local.minimum.label.offset.lines,
        "local.minimum.label.offset.lines"
    )
    n.legend.position <- validate_conplot_legend_position(
        n.legend.position, "n.legend.position"
    )
    n.legend.cex <- validate_conplot_optional_positive_scalar(
        n.legend.cex, "n.legend.cex"
    )
    if (is.null(n.legend.cex)) {
        n.legend.cex <- 1
    }
    n.legend.col <- validate_conplot_colour_scalar(n.legend.col, "n.legend.col")
    n.legend.bty <- validate_conplot_character_scalar(n.legend.bty, "n.legend.bty")

    disease.cex <- validate_conplot_positive_scalar(disease.cex, "disease.cex")
    disease.point.cex <- validate_conplot_positive_scalar(
        disease.point.cex, "disease.point.cex"
    )
    cex.lab <- validate_conplot_optional_positive_scalar(cex.lab, "cex.lab")
    cex.axis <- validate_conplot_optional_positive_scalar(cex.axis, "cex.axis")
    cex.main <- validate_conplot_optional_positive_scalar(cex.main, "cex.main")

    local.minimum.label.source <- local.minimum.label
    n.legend.label.source <- n.legend.label

    xlab <- conplot_nice_label(xlab, tikz.mode = use.tikz, warn = label.warn)
    ylab <- conplot_nice_label(ylab, tikz.mode = use.tikz, warn = label.warn)
    main <- conplot_nice_label(main, tikz.mode = use.tikz, warn = label.warn)
    local.minimum.label <- conplot_nice_label(
        local.minimum.label,
        tikz.mode = use.tikz,
        warn = label.warn
    )
    n.legend.label <- conplot_nice_label(
        n.legend.label,
        tikz.mode = use.tikz,
        warn = label.warn
    )

    disease.data <- resolve_conplot_disease_data(
        disease.data = disease.data,
        disease.names = disease.names,
        show.diseases = show.diseases
    )

    zlim <- range(prob, na.rm = TRUE)
    axis.at <- 2^(0:5)
    axis.at <- axis.at[axis.at >= ylim[1L] & axis.at <= ylim[2L]]

    plot.args <- list(
        x = NA,
        y = NA,
        las = 1,
        log = log,
        xlim = xlim,
        ylim = ylim,
        xaxs = "i",
        yaxs = "i",
        yaxt = "n",
        xlab = xlab,
        ylab = ylab,
        main = main
    )
    plot.args <- c(plot.args, list(...))
    if (!is.null(cex.lab)) {
        plot.args$cex.lab <- cex.lab
    }
    if (!is.null(cex.axis)) {
        plot.args$cex.axis <- cex.axis
    }
    if (!is.null(cex.main)) {
        plot.args$cex.main <- cex.main
    }

    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)

    if (colour.legend) {
        graphics::par(fig = c(0, 0.82, 0, 1))
    }

    do.call(graphics::plot.default, plot.args)

    if (!filled) {
        graphics::abline(
            h = 2^(1:4),
            v = seq(0.005, 0.015, by = 0.005),
            col = grDevices::grey(0.95),
            lwd = contour.lwd
        )
    }

    axis.args <- list(side = 2, at = axis.at, labels = axis.at, las = 1)
    if (!is.null(cex.axis)) {
        axis.args$cex.axis <- cex.axis
    }
    do.call(graphics::axis, axis.args)

    if (filled) {
        graphics::image(
            x = epsilon,
            y = R0,
            z = prob,
            col = fill.colours,
            breaks = fill.breaks,
            log = log,
            add = TRUE
        )
    }

    graphics::contour(
        x = epsilon,
        y = R0,
        z = prob,
        log = log,
        levels = levels,
        labels = contour.labels,
        add = TRUE,
        col = contour.col,
        lwd = contour.lwd,
        lty = contour.lty,
        labcex = label.cex
    )

    manual.labels <- NULL
    if (show.manual.labels) {
        manual.labels <- draw_conplot_manual_labels(
            epsilon = epsilon,
            R0 = R0,
            prob = prob,
            high.levels = high.levels,
            low.levels = low.levels,
            x = manual.label.x,
            y = manual.label.y,
            grey.adjust = grey.adjust,
            point.cex = manual.label.point.cex,
            label.cex = label.cex,
            use.tikz = use.tikz,
            label.warn = label.warn
        )
    }

    quadratic.points <- NULL
    if (show.quadratic) {
        quadratic.points <- draw_conplot_quadratic(
            coefficients = quadratic.coefficients,
            from = quadratic.from,
            to = quadratic.to,
            col = quadratic.col,
            lty = quadratic.lty,
            lwd = quadratic.lwd
        )
    }

    local.minimum <- NULL
    local.minimum.label.metadata <- NULL
    if (show.local.minimum) {
        local.minimum <- conplot_local_min_curve(
            epsilon,
            R0,
            prob,
            xlow = local.minimum.xlow,
            ylow = local.minimum.ylow,
            yhigh = local.minimum.yhigh
        )
        local.minimum.available <- any(is.finite(local.minimum$x) & is.finite(local.minimum$y))
        if (local.minimum.available) {
            graphics::lines(
                local.minimum$x,
                local.minimum$y,
                col = local.minimum.col,
                lwd = local.minimum.lwd
            )
        } else {
            warning(
                "local-minimum curve has no finite points; skipping the ",
                "local-minimum overlay and label.",
                call. = FALSE
            )
        }
        if (show.local.minimum.label && local.minimum.available) {
            local.minimum.label.metadata <- draw_conplot_local_minimum_label(
                label = local.minimum.label,
                label.source = local.minimum.label.source,
                coefficients = quadratic.coefficients,
                from = quadratic.from,
                to = quadratic.to,
                position = local.minimum.label.position,
                offset.lines = local.minimum.label.offset.lines,
                col = local.minimum.label.col,
                cex = local.minimum.label.cex,
                log = log
            )
        }
    }

    diseases.plotted <- FALSE
    if (show.diseases && !is.null(disease.data)) {
        draw_conplot_diseases(
            disease.data = disease.data,
            point.bg = disease.point.bg,
            point.cex = disease.point.cex,
            label.cex = disease.cex,
            box.col = disease.box.col
        )
        diseases.plotted <- TRUE
    }

    n.legend.metadata <- NULL
    if (show.n.legend) {
        n.legend.metadata <- draw_conplot_n_legend(
            label = n.legend.label,
            label.source = n.legend.label.source,
            position = n.legend.position,
            col = n.legend.col,
            cex = n.legend.cex,
            bty = n.legend.bty
        )
    }

    graphics::box()

    if (colour.legend) {
        draw_conplot_colour_legend(
            fill.breaks = fill.breaks,
            fill.colours = fill.colours
        )
    }

    out <- list(
        levels = levels,
        high.levels = high.levels,
        low.levels = low.levels,
        fill.breaks = fill.breaks,
        fill.colours = fill.colours,
        xlim = xlim,
        ylim = ylim,
        zlim = zlim,
        log = log,
        filled = filled,
        colour.legend = colour.legend,
        label.cex = label.cex,
        label.bg.cex.mult = label.bg.cex.mult,
        manual.label.point.cex = manual.label.point.cex,
        use.tikz = use.tikz,
        show.manual.labels = show.manual.labels,
        manual.labels = manual.labels,
        show.quadratic = show.quadratic,
        quadratic.points = quadratic.points,
        show.local.minimum = show.local.minimum,
        local.minimum = local.minimum,
        show.local.minimum.label = show.local.minimum.label,
        local.minimum.label = local.minimum.label.metadata,
        show.n.legend = show.n.legend,
        n.legend = n.legend.metadata,
        show.diseases = show.diseases,
        diseases.plotted = diseases.plotted,
        labels = list(
            xlab = xlab,
            ylab = ylab,
            main = main
        ),
        orientation = "rows: epsilon; columns: R0",
        axes = "x: epsilon; y: R0"
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

##' Validate a finite scalar
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return A numeric scalar.
##' @noRd
validate_conplot_finite_scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
        stop(name, " must be a finite numeric scalar.", call. = FALSE)
    }
    x
}

##' Validate a positive scalar
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return A numeric scalar.
##' @noRd
validate_conplot_positive_scalar <- function(x, name) {
    x <- validate_conplot_finite_scalar(x, name)
    if (x <= 0) {
        stop(name, " must be a positive finite scalar.", call. = FALSE)
    }
    x
}

##' Validate an optional positive scalar
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return `NULL` or a numeric scalar.
##' @noRd
validate_conplot_optional_positive_scalar <- function(x, name) {
    if (is.null(x)) {
        return(NULL)
    }
    validate_conplot_positive_scalar(x, name)
}

##' Validate an optional logical scalar
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return `NULL` or a logical scalar.
##' @noRd
validate_conplot_optional_logical_scalar <- function(x, name) {
    if (is.null(x)) {
        return(NULL)
    }
    validate_conplot_logical_scalar(x, name)
}

##' Resolve tikz label mode for conplot labels
##'
##' @param use.tikz explicit tikz mode or `NULL`.
##' @param envir caller environment.
##' @param warn logical scalar controlling warnings.
##'
##' @return `NULL` or a scalar logical value.
##' @noRd
resolve_conplot_use_tikz <- function(use.tikz, envir, warn = TRUE) {
    warn <- validate_conplot_logical_scalar(warn, "warn")

    if (!is.null(use.tikz)) {
        return(validate_conplot_logical_scalar(use.tikz, "use.tikz"))
    }

    if (exists("use.tikz", envir = envir, inherits = FALSE)) {
        caller.use.tikz <- get("use.tikz", envir = envir, inherits = FALSE)
        if (is.logical(caller.use.tikz) &&
            length(caller.use.tikz) == 1L &&
            !is.na(caller.use.tikz)) {
            return(caller.use.tikz)
        }
        if (warn) {
            warning(
                "Ignoring caller-level `use.tikz` because it is not a scalar logical value.",
                call. = FALSE
            )
        }
    }

    NULL
}

##' Validate a probability scalar
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return A numeric scalar.
##' @noRd
validate_conplot_probability_scalar <- function(x, name) {
    x <- validate_conplot_finite_scalar(x, name)
    if (x < 0 || x > 1) {
        stop(name, " must be in [0, 1].", call. = FALSE)
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
    if (grepl("x", log, fixed = TRUE) && any(epsilon <= 0)) {
        stop("epsilon values must be positive when log includes the x axis.", call. = FALSE)
    }
    if (grepl("y", log, fixed = TRUE) && any(R0 <= 0)) {
        stop("R0 values must be positive when log includes the y axis.", call. = FALSE)
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

##' Validate a character scalar
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return A character scalar.
##' @noRd
validate_conplot_character_scalar <- function(x, name) {
    if (!is.character(x) || length(x) != 1L || is.na(x)) {
        stop(name, " must be a character scalar.", call. = FALSE)
    }
    x
}

##' Validate a colour scalar
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return A character scalar.
##' @noRd
validate_conplot_colour_scalar <- function(x, name) {
    validate_conplot_character_scalar(x, name)
}

##' Validate an optional colour scalar
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return `NULL` or a character scalar.
##' @noRd
validate_conplot_optional_colour_scalar <- function(x, name) {
    if (is.null(x)) {
        return(NULL)
    }
    if (!is.character(x) || length(x) != 1L || is.na(x)) {
        stop(name, " must be NULL or a character scalar.", call. = FALSE)
    }
    x
}

##' Validate an optional two-coordinate position
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return `NULL` or a numeric vector with `epsilon` and `R0` entries.
##' @noRd
validate_conplot_optional_position <- function(x, name) {
    if (is.null(x)) {
        return(NULL)
    }
    if (!is.numeric(x) || length(x) != 2L || anyNA(x) || any(!is.finite(x))) {
        stop(name, " must be NULL or a finite numeric vector of length 2.", call. = FALSE)
    }
    names(x) <- c("epsilon", "R0")
    x
}

##' Validate a conplot legend position
##'
##' @param x object to validate.
##' @param name argument name for error messages.
##'
##' @return A character scalar or numeric vector with `x` and `y` entries.
##' @noRd
validate_conplot_legend_position <- function(x, name) {
    if (is.character(x) && length(x) == 1L && !is.na(x)) {
        valid.positions <- c(
            "bottomright", "bottom", "bottomleft", "left", "topleft",
            "top", "topright", "right", "center"
        )
        if (!x %in% valid.positions) {
            stop(
                name, " must be a valid legend keyword or a finite numeric ",
                "vector of length 2.",
                call. = FALSE
            )
        }
        return(x)
    }
    if (is.numeric(x) && length(x) == 2L && !anyNA(x) && all(is.finite(x))) {
        names(x) <- c("x", "y")
        return(x)
    }
    stop(
        name, " must be a valid legend keyword or a finite numeric vector of length 2.",
        call. = FALSE
    )
}

##' Validate quadratic coefficients
##'
##' @param x coefficient vector.
##'
##' @return A numeric vector with intercept, slope, and curvature.
##' @noRd
validate_conplot_quadratic_coefficients <- function(x) {
    if (!is.numeric(x) || length(x) != 3L || anyNA(x) || any(!is.finite(x))) {
        stop(
            "quadratic.coefficients must be a finite numeric vector of length 3.",
            call. = FALSE
        )
    }
    names(x) <- c("intercept", "slope", "curvature")
    x
}

##' Default manual conplot label size
##'
##' @return A positive numeric scalar.
##' @noRd
conplot_default_label_cex <- function() {
    if (interactive()) {
        return(0.7)
    }
    0.5
}

##' Prepare a conplot text label
##'
##' @param x label object.
##' @param tikz.mode optional logical tikz mode.
##' @param warn logical scalar controlling warnings.
##'
##' @return A label suitable for base graphics.
##' @noRd
conplot_nice_label <- function(x, tikz.mode = NULL, warn = TRUE) {
    if (!is.character(x)) {
        return(x)
    }

    if (is.null(tikz.mode)) {
        return(earnmisc::nice_text(x, warn = warn))
    }

    earnmisc::nice_text(x, use.tikz = tikz.mode, warn = warn)
}

##' Interpolate a grid column profile at an epsilon value
##'
##' @param epsilon epsilon grid.
##' @param R0 reproduction-number grid.
##' @param prob probability grid.
##' @param x epsilon value.
##'
##' @return Numeric vector of probabilities over `R0`.
##' @noRd
conplot_profile_at_epsilon <- function(epsilon, R0, prob, x) {
    vapply(
        seq_along(R0),
        function(j) conplot_approx_or_na(epsilon, prob[, j], xout = x),
        numeric(1L)
    )
}

##' Interpolate a grid row profile at an R0 value
##'
##' @param epsilon epsilon grid.
##' @param R0 reproduction-number grid.
##' @param prob probability grid.
##' @param y reproduction-number value.
##'
##' @return Numeric vector of probabilities over `epsilon`.
##' @noRd
conplot_profile_at_R0 <- function(epsilon, R0, prob, y) {
    vapply(
        seq_along(epsilon),
        function(i) conplot_approx_or_na(R0, prob[i, ], xout = y),
        numeric(1L)
    )
}

##' Interpolate a profile when enough finite values are available
##'
##' @param x coordinate vector.
##' @param y profile values.
##' @param xout coordinate to interpolate.
##'
##' @return Interpolated numeric scalar, or `NA_real_` when interpolation is
##'   unsupported.
##' @noRd
conplot_approx_or_na <- function(x, y, xout) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 2L) {
        return(NA_real_)
    }
    stats::approx(x[ok], y[ok], xout = xout, rule = 2)$y
}

##' Find a contour crossing by interpolation
##'
##' @param x coordinate vector.
##' @param z profile values.
##' @param level contour level.
##'
##' @return A scalar coordinate or `NA_real_`.
##' @noRd
conplot_find_crossing <- function(x, z, level) {
    ok <- is.finite(x) & is.finite(z)
    x <- x[ok]
    z <- z[ok]
    if (length(x) < 2L) {
        return(NA_real_)
    }
    d <- z - level
    exact <- which(d == 0)
    if (length(exact)) {
        return(x[exact[[1L]]])
    }
    crossings <- which(d[-length(d)] * d[-1L] < 0)
    if (!length(crossings)) {
        return(NA_real_)
    }
    j <- crossings[[1L]]
    stats::approx(z[j:(j + 1L)], x[j:(j + 1L)], xout = level)$y
}

##' Draw original manual conplot labels
##'
##' @param epsilon epsilon grid.
##' @param R0 reproduction-number grid.
##' @param prob probability grid.
##' @param high.levels high probability levels.
##' @param low.levels low probability levels.
##' @param x epsilon value for high-level labels.
##' @param y reproduction number for low-level labels.
##' @param grey.adjust greyscale adjustment.
##' @param point.cex point size.
##' @param label.cex label size.
##' @param use.tikz optional logical tikz mode.
##' @param label.warn logical scalar controlling warnings.
##'
##' @return Data frame of manual label positions.
##' @noRd
draw_conplot_manual_labels <- function(epsilon,
                                       R0,
                                       prob,
                                       high.levels,
                                       low.levels,
                                       x,
                                       y,
                                       grey.adjust,
                                       point.cex,
                                       label.cex,
                                       use.tikz,
                                       label.warn) {
    high.profile <- conplot_profile_at_epsilon(epsilon, R0, prob, x)
    high.y <- vapply(
        high.levels,
        function(level) conplot_find_crossing(R0, high.profile, level),
        numeric(1L)
    )
    high.ok <- is.finite(high.y)
    high.cols <- rev(grDevices::grey(
        1 - grey.adjust + pmin(high.levels + 0.1, 1) * grey.adjust
    ))

    if (any(high.ok)) {
        graphics::points(
            x = rep(x, sum(high.ok)),
            y = high.y[high.ok],
            pch = 19,
            cex = point.cex,
            col = high.cols[high.ok]
        )
        graphics::text(
            x = rep(x, sum(high.ok)),
            y = high.y[high.ok],
            labels = high.levels[high.ok],
            col = "black",
            cex = label.cex
        )
    }

    low.profile <- conplot_profile_at_R0(epsilon, R0, prob, y)
    low.x <- vapply(
        low.levels,
        function(level) conplot_find_crossing(epsilon, low.profile, level),
        numeric(1L)
    )
    low.ok <- is.finite(low.x)

    if (any(low.ok)) {
        graphics::points(
            x = low.x[low.ok],
            y = rep(y, sum(low.ok)),
            pch = 19,
            col = "white",
            cex = point.cex
        )
        graphics::text(
            x = low.x[low.ok],
            y = rep(y, sum(low.ok)),
            labels = conplot_low_level_labels(
                low.levels[low.ok],
                use.tikz = use.tikz,
                label.warn = label.warn
            ),
            col = "black",
            cex = label.cex
        )
    }

    rbind(
        data.frame(
            type = "high",
            level = high.levels,
            x = x,
            y = high.y,
            label.cex = label.cex,
            label.bg.cex = point.cex
        ),
        data.frame(
            type = "low",
            level = low.levels,
            x = low.x,
            y = y,
            label.cex = label.cex,
            label.bg.cex = point.cex
        )
    )
}

##' Build low-level label expressions
##'
##' @param x low-level values.
##' @param use.tikz optional logical tikz mode.
##' @param label.warn logical scalar controlling warnings.
##'
##' @return Character vector or expression vector.
##' @noRd
conplot_low_level_labels <- function(x, use.tikz, label.warn) {
    powers <- as.integer(round(log10(x)))
    conplot_nice_label(
        paste0("$10^{", powers, "}$"),
        tikz.mode = use.tikz,
        warn = label.warn
    )
}

##' Draw the original quadratic overlay
##'
##' @param coefficients quadratic coefficients.
##' @param from,to plotting range.
##' @param col,lty,lwd graphical parameters.
##'
##' @return Data frame of curve coordinates.
##' @noRd
draw_conplot_quadratic <- function(coefficients, from, to, col, lty, lwd) {
    x <- seq(from, to, length.out = 200L)
    y <- conplot_quadratic_value(coefficients, x)
    graphics::lines(x, y, col = col, lty = lty, lwd = lwd)
    data.frame(x = x, y = y)
}

##' Evaluate the quadratic conplot overlay
##'
##' @param coefficients quadratic coefficients.
##' @param x epsilon values.
##'
##' @return Numeric `R0` values on the quadratic overlay.
##' @noRd
conplot_quadratic_value <- function(coefficients, x) {
    coefficients[["intercept"]] +
        coefficients[["slope"]] * x +
        coefficients[["curvature"]] * x^2
}

##' Draw the local-minimum-curve label
##'
##' @param label label to draw.
##' @param label.source original label before device-specific conversion.
##' @param coefficients quadratic coefficients used for placement and slope.
##' @param from,to horizontal range for the quadratic overlay.
##' @param position optional explicit curve-anchor position.
##' @param offset.lines vertical offset in text-line heights.
##' @param col,cex graphical parameters for the label.
##' @param log base-graphics log setting.
##'
##' @return List of label metadata.
##' @noRd
draw_conplot_local_minimum_label <- function(label,
                                             label.source,
                                             coefficients,
                                             from,
                                             to,
                                             position,
                                             offset.lines,
                                             col,
                                             cex,
                                             log) {
    if (is.null(position)) {
        x <- mean(c(from, to))
        position <- c(
            epsilon = x,
            R0 = conplot_quadratic_value(coefficients, x)
        )
    }
    curve.position <- position
    offset <- conplot_offset_position_above(
        position = curve.position,
        cex = cex,
        offset.lines = offset.lines,
        log = log
    )
    label.position <- offset$position

    srt <- conplot_quadratic_label_angle(
        coefficients = coefficients,
        x = curve.position[["epsilon"]],
        from = from,
        to = to,
        log = log
    )

    graphics::text(
        x = label.position[["epsilon"]],
        y = label.position[["R0"]],
        labels = label,
        col = col,
        cex = cex,
        srt = srt,
        adj = c(0.5, 0.5),
        xpd = NA
    )

    list(
        label = label.source,
        plotting.label = label,
        position = label.position,
        curve.position = curve.position,
        offset.lines = offset.lines,
        offset.inches = offset$offset.inches,
        col = col,
        cex = cex,
        srt = srt
    )
}

##' Offset a label position visually upward
##'
##' @param position numeric vector with `epsilon` and `R0` entries.
##' @param cex text expansion factor.
##' @param offset.lines offset in text-line heights.
##' @param log base-graphics log setting.
##'
##' @return List with the offset position and physical offset in inches.
##' @noRd
conplot_offset_position_above <- function(position, cex, offset.lines, log) {
    offset.inches <- offset.lines * graphics::strheight(
        "M",
        units = "inches",
        cex = cex
    )
    if (offset.inches == 0) {
        return(list(position = position, offset.inches = offset.inches))
    }

    user.y <- conplot_log_transform(position[["R0"]], axis = "y", log = log)
    inch.y <- graphics::grconvertY(user.y, from = "user", to = "inches")
    shifted.user.y <- graphics::grconvertY(
        inch.y + offset.inches,
        from = "inches",
        to = "user"
    )
    shifted.y <- conplot_inverse_log_transform(
        shifted.user.y,
        axis = "y",
        log = log
    )

    out <- position
    out[["R0"]] <- shifted.y
    list(position = out, offset.inches = offset.inches)
}

##' Draw the sample-size legend
##'
##' @param label label to draw.
##' @param label.source original label before device-specific conversion.
##' @param position legend keyword or coordinates.
##' @param col,cex,bty graphical parameters for the legend text.
##'
##' @return List of legend metadata.
##' @noRd
draw_conplot_n_legend <- function(label,
                                  label.source,
                                  position,
                                  col,
                                  cex,
                                  bty) {
    legend.args <- list(
        legend = label,
        text.col = col,
        cex = cex,
        bty = bty
    )
    if (is.character(position)) {
        legend.args$x <- position
    } else {
        legend.args$x <- position[["x"]]
        legend.args$y <- position[["y"]]
    }

    do.call(graphics::legend, legend.args)

    list(
        label = label.source,
        plotting.label = label,
        position = position,
        col = col,
        cex = cex,
        bty = bty
    )
}

##' Compute visual rotation for the quadratic conplot overlay
##'
##' @param coefficients quadratic coefficients.
##' @param x epsilon value at the label centre.
##' @param from,to horizontal range for the quadratic overlay.
##' @param log base-graphics log setting.
##'
##' @return Rotation angle in degrees for `text(srt = ...)`.
##' @noRd
conplot_quadratic_label_angle <- function(coefficients, x, from, to, log) {
    span <- to - from
    dx <- span * 0.01
    x1 <- max(from, x - dx)
    x2 <- min(to, x + dx)
    if (x1 >= x2) {
        return(0)
    }

    y1 <- conplot_quadratic_value(coefficients, x1)
    y2 <- conplot_quadratic_value(coefficients, x2)
    if (!all(is.finite(c(x1, x2, y1, y2))) ||
        (grepl("x", log, fixed = TRUE) && any(c(x1, x2) <= 0)) ||
        (grepl("y", log, fixed = TRUE) && any(c(y1, y2) <= 0))) {
        return(0)
    }

    plot.x <- conplot_log_transform(c(x1, x2), axis = "x", log = log)
    plot.y <- conplot_log_transform(c(y1, y2), axis = "y", log = log)
    if (!all(is.finite(c(plot.x, plot.y)))) {
        return(0)
    }
    inch.x <- graphics::grconvertX(plot.x, from = "user", to = "inches")
    inch.y <- graphics::grconvertY(plot.y, from = "user", to = "inches")
    as.numeric(atan2(diff(inch.y), diff(inch.x)) * 180 / pi)
}

##' Transform coordinates to base-graphics user coordinates for log axes
##'
##' @param x coordinate values.
##' @param axis `"x"` or `"y"`.
##' @param log base-graphics log setting.
##'
##' @return Numeric coordinates in the active plot's user coordinate system.
##' @noRd
conplot_log_transform <- function(x, axis, log) {
    if (grepl(axis, log, fixed = TRUE)) {
        return(log10(x))
    }
    x
}

##' Transform base-graphics user coordinates back from log axes
##'
##' @param x coordinate values in active plot user coordinates.
##' @param axis `"x"` or `"y"`.
##' @param log base-graphics log setting.
##'
##' @return Numeric coordinates on the original plotting scale.
##' @noRd
conplot_inverse_log_transform <- function(x, axis, log) {
    if (grepl(axis, log, fixed = TRUE)) {
        return(10^x)
    }
    x
}

##' Estimate the original local-minimum curve
##'
##' @param x epsilon grid.
##' @param y reproduction-number grid.
##' @param z probability grid.
##' @param xlow,ylow,yhigh limits for accepting local minima.
##'
##' @return Data frame of local-minimum coordinates.
##' @noRd
conplot_local_min_curve <- function(x, y, z, xlow = 0.004, ylow = 1, yhigh = Inf) {
    if (length(x) != dim(z)[1L] || length(y) != dim(z)[2L]) {
        stop("local-minimum grid dimensions are inconsistent.", call. = FALSE)
    }
    ymin <- rep(NA_real_, length(x))
    yminj <- rep(NA_integer_, length(x))

    for (i in seq_along(x)) {
        if (sum(z[i, ] == 0, na.rm = TRUE) > 10L || x[i] < xlow) {
            next
        }
        zi <- ifelse(y >= ylow & y <= yhigh, z[i, ], NA_real_)
        if (all(is.na(zi))) {
            next
        }
        j <- which.min(zi)
        finite.index <- which(!is.na(zi))
        if (!length(finite.index) || j == min(finite.index)) {
            next
        }
        yminj[i] <- j
        ymin[i] <- y[j]
    }

    data.frame(x = x, y = ymin, xi = seq_along(x), yj = yminj)
}

##' Resolve disease annotations
##'
##' @param disease.data optional disease data.
##' @param disease.names disease identifiers to retain.
##' @param show.diseases whether annotations are requested.
##'
##' @return A data frame or `NULL`.
##' @noRd
resolve_conplot_disease_data <- function(disease.data, disease.names, show.diseases) {
    if (!show.diseases) {
        return(NULL)
    }
    if (is.null(disease.data)) {
        disease.data <- load_conplot_sirr_disease_data()
        if (is.null(disease.data)) {
            return(NULL)
        }
    }

    required.names <- c("epsilon", "R0", "label")
    missing.names <- setdiff(required.names, names(disease.data))
    if (length(missing.names)) {
        stop(
            "disease.data is missing required column(s): ",
            paste(missing.names, collapse = ", "),
            call. = FALSE
        )
    }

    if ("disease" %in% names(disease.data)) {
        disease.data <- disease.data[disease.data$disease %in% disease.names, , drop = FALSE]
        disease.data <- adjust_conplot_disease_data(disease.data)
    }

    disease.data
}

##' Load disease annotations from sirr when available
##'
##' @return A data frame or `NULL`.
##' @noRd
load_conplot_sirr_disease_data <- function() {
    if (!requireNamespace("sirr", quietly = TRUE)) {
        return(NULL)
    }
    env <- new.env(parent = emptyenv())
    utils::data("diseaseParameters", package = "sirr", envir = env)
    if (!exists("diseaseParameters", envir = env, inherits = FALSE)) {
        return(NULL)
    }
    get("diseaseParameters", envir = env, inherits = FALSE)
}

##' Apply original disease-label nudges
##'
##' @param disease.data disease data frame.
##'
##' @return Adjusted disease data frame.
##' @noRd
adjust_conplot_disease_data <- function(disease.data) {
    disease.data[disease.data$disease == "measles", "R0"] <- 18
    disease.data[disease.data$disease == "pertussis", "R0"] <- 16
    disease.data[disease.data$disease == "covid19delta", "R0"] <- 7.4
    disease.data
}

##' Draw disease annotations
##'
##' @param disease.data disease data frame.
##' @param point.bg,box.col colours.
##' @param point.cex,label.cex size controls.
##'
##' @return `NULL`, invisibly.
##' @noRd
draw_conplot_diseases <- function(disease.data,
                                  point.bg,
                                  point.cex,
                                  label.cex,
                                  box.col) {
    if (nrow(disease.data) == 0L) {
        return(invisible(NULL))
    }
    non.hiv <- disease.data$label != "HIV"
    if (any(non.hiv)) {
        conplot_white_out_string(
            x = disease.data$epsilon[non.hiv],
            y = disease.data$R0[non.hiv],
            string = disease.data$label[non.hiv],
            col = box.col,
            cex = label.cex
        )
    }
    graphics::points(
        x = disease.data$epsilon,
        y = disease.data$R0,
        pch = 21,
        bg = point.bg,
        cex = point.cex
    )
    graphics::text(
        x = disease.data$epsilon,
        y = disease.data$R0,
        labels = disease.data$label,
        pos = 4,
        xpd = NA,
        cex = label.cex
    )
    invisible(NULL)
}

##' White out a label background
##'
##' @param x,y label anchor coordinates.
##' @param string label text.
##' @param col fill colour.
##' @param cex text size.
##'
##' @return `NULL`, invisibly.
##' @noRd
conplot_white_out_string <- function(x, y, string, col = "white", cex = 1) {
    sw <- graphics::strwidth(string) * cex
    sh <- graphics::strheight(string) * cex
    xfrsz <- 0.0001
    yfrsz <- 0.1
    xshift <- 0.0004

    graphics::rect(
        xleft = xshift + x - xfrsz,
        ybottom = y * 2^(-(sh / 2 + yfrsz)),
        xright = xshift + x + sw + xfrsz,
        ytop = y * 2^(sh / 2 + yfrsz),
        col = col,
        border = NA
    )
    invisible(NULL)
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
