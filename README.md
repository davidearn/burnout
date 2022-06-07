
The `burnout` package requires R ≥ 4.2.0 in order for the documentation to be displayed as desired.  As mentioned [here](https://cran.r-project.org/doc/manuals/r-devel/NEWS.html), 4.2.0 automatically supports [KaTeX](https://katex.org/docs/support_table.html) and [MathJax](https://www.mathjax.org/), which I use liberally in [roxygen](https://roxygen2.r-lib.org/) documentation.

# generating graphs from Parsons _et al_ (2022)

Generate the analytical curves in Figure 3 via
```
example(plot_P1)
```

Generate Figures B1 and B2 via
```
example(plot_x_in)
```
Note that the green/blue curves (i.e., crude/better approximations of
peak prevalence `ymax`) in Todd's plots correspond to dotted/solid
curves in this plot.

The function `compare_funs` creates a data frame with values of two
functions of `R0` and `epsilon` at a grid of points.  It returns an
object of class `compare_funs`, for which there is a `plot` method.
So, for example, the following compares the exact and approximate
expressions for Kendalls's q:
```
cf <- compare_funs()
plot(cf)
```
To compare the crude and better approximations of `xin`, use
```
cxin <- compare_funs(x_in, x_in_crude, Rmax=20)
plot(cxin)
```
In the above plot (with `Rmax = 20`) the bottom two panels (for `epsilon = 0.01` and
`0.1`) correspond to Figure B1 in the paper.

To compare van Herwaarden's (1997) approximation to ours:
```
cp <- compare_funs(P1_prob, vanH_prob)
plot(cp)
```
_Note, however, that the above works only because of a `try` catch inside `vanH_prob`.  We actually need to resolve this._

A better test is to compare the burnout probabilities (conditional on
not fizzling), as opposed to the persistence probabiities, because
`P1` is derived from the burnout probability identically in all cases.
The formulae that are different are the components from burning out
conditional on not fizzling:
```
cb <- compare_funs(burnout_prob, burnout_prob_vanH)
plot(cb)
```
