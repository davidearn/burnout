## Recreate the manuscript-style SIRS conplot from the preserved probability grid.
##
## Run from the package root:
##
##     Rscript sandbox/tracked/conplot_talk_figure_sirs.R
##
## Input:
##
##     ../sources/conplot_standalone/prob.RData.eta=0.01
##
## In an interactive RStudio session, the default draws to the active graphics
## device. With Rscript, the default writes:
##
##     sandbox/conplot_talk_figure_sirs.pdf
##
## Set output.mode to "tikz" to write and compile:
##
##     sandbox/conplot_talk_figure_sirs.tex

helper.candidates <- c(
    file.path("sandbox", "tracked", "conplot_funs.R"),
    file.path("..", "..", "sandbox", "tracked", "conplot_funs.R")
)
helper.file <- helper.candidates[file.exists(helper.candidates)]
if (!length(helper.file)) {
    stop("Could not find sandbox/tracked/conplot_funs.R.", call. = FALSE)
}
source(helper.file[[1L]], local = TRUE)

output.mode <- if (interactive()) "device" else "pdf"
##output.mode <- "tikz"

conplot.model <- "sirs"
prob.file.name <- "prob.RData.eta=0.01"
output.basename <- "conplot_talk_figure_sirs"

device.settings <- list(
    output.file = NULL,
    width = 6,
    height = 6
)

## Edit this list to pass any plot_conplot_grid() argument through to the
## plotting helper. The SIRS grid has data-derived local minima below the
## original SIR xlow cutoff, so local.minimum.xlow is set to zero.
plot.args <- list(
    label.cex = 1,
    label.bg.cex.mult = 4,
    colour.legend = FALSE,
    low.levels = c(0),
    use.tikz = NULL,
    cex.lab = NULL,
    cex.axis = NULL,
    cex.main = NULL,
    show.diseases = TRUE,
    show.manual.labels = TRUE,
    show.quadratic = FALSE,
    show.local.minimum = TRUE,
    local.minimum.xlow = 0,
    show.local.minimum.label = FALSE,
    local.minimum.label = "minimum persistence probability",
    local.minimum.label.cex = NULL,
    local.minimum.label.col = NULL,
    local.minimum.label.position = NULL,
    local.minimum.label.offset.lines = 1,
    show.n.legend = TRUE,
    n.legend.label = "$n = 10^6$",
    n.legend.position = "topright",
    n.legend.cex = NULL,
    n.legend.col = "grey95",
    n.legend.bty = "n"
)

if (!isTRUE(getOption("burnout.conplot_talk_figure.skip_main", FALSE))) {
    plot.info <- run_conplot_script(
        conplot.model = conplot.model,
        prob.file.name = prob.file.name,
        output.basename = output.basename,
        device.settings = device.settings,
        plot.args = plot.args,
        output.mode = output.mode
    )
}
