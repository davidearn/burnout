## Compatibility wrapper for the split conplot talk scripts.
##
## Prefer running one of:
##
##     Rscript sandbox/tracked/conplot_talk_figure_sir.R
##     Rscript sandbox/tracked/conplot_talk_figure_sirs.R
##
## This wrapper keeps the older command useful by running the SIRS script.

message(
    "conplot_talk_figure.R is a compatibility wrapper; ",
    "use conplot_talk_figure_sir.R or conplot_talk_figure_sirs.R directly."
)

script.candidates <- c(
    file.path("sandbox", "tracked", "conplot_talk_figure_sirs.R"),
    file.path("..", "..", "sandbox", "tracked", "conplot_talk_figure_sirs.R")
)
script.file <- script.candidates[file.exists(script.candidates)]
if (!length(script.file)) {
    stop("Could not find sandbox/tracked/conplot_talk_figure_sirs.R.", call. = FALSE)
}
source(script.file[[1L]], local = TRUE)
