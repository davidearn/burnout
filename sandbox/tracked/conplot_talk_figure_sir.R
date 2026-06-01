## Recreate the manuscript-style SIR conplot from the preserved probability grid.
##
## Run from the package root:
##
##     Rscript sandbox/tracked/conplot_talk_figure_sir.R
##
## Input:
##
##     ../sources/conplot_standalone/prob_sir.RData
##
## In an interactive RStudio session, the default draws to the active graphics
## device. With Rscript, the default writes:
##
##     sandbox/conplot_talk_figure_sir.pdf
##
## Set output.mode to "tikz" to write and compile:
##
##     sandbox/conplot_talk_figure_sir.tex

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

conplot.model <- "sir"
prob.file.name <- "prob_sir.RData"
output.basename <- "conplot_talk_figure_sir"

device.settings <- list(
    output.file = NULL,
    width = 6,
    height = 6
)

## Edit this list to pass any plot_conplot_grid() argument through to the
## plotting helper.
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
