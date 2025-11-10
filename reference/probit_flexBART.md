# Probit flexBART for binary outcomes

Fit a BART model of a binary responses using the the Albert and Chib
(1993) data augmentation for probit models.

## Usage

``` r
probit_flexBART(formula, train_data, test_data = NULL,...)
```

## Arguments

- formula:

  an object of class [`formula`](https://rdrr.io/r/stats/formula.html)
  (or one that can be coerced to the class): a symbolic description of
  the model to be fitted. The details of model specification are given
  under 'Details'.

- train_data:

  an object of class `data.frame` containing data used to train the
  model. As usual, rows (resp. columns) correspond to observations
  (resp. variables)

- test_data:

  an optional object of class `data.frame` containing test-set (i.e.,
  out-of-sample) data. Default is `NULL`.

- ...:

  Additional arguments for setting prior hyperparameters (e.g., number
  of trees, \\\mu\_{0}\\, \\\tau\\, etc.) and MCMC control parameters
  (e.g., number of chains, iterations, etc.) See
  [`flexBART`](https://skdeshpande91.github.io/flexBART/reference/flexBART.md)
  for details about additional arguments.

## Details

For a binary response \\Y\\, \\p\\ predictors \\X\_{1}, \ldots,
X\_{p}\\, `probit_flexBART()` models \\P(Y=1 \vert X = x) = \Phi(f(x)\\,
where \\\Phi\\ is the standard normal cumulative distribution function.
`probit_flexBART()` combines the Albert & Chib (1993) data augmentation
strategy for probit regression with the usual Bayesian backfitting used
to fit (VC)BART models.

### The formula argument

**Currently, `probit_flexBART()` only supports fitting single ensemble
probit BART models.** So, the only valid formula will look something
like `Y~bart(.)` or `Y ~ bart(x1+x2)`. As with
[`flexBART()`](https://skdeshpande91.github.io/flexBART/reference/flexBART.md),
you must include the string “bart” on the right-hand side of the formula
object.

### Prior specification

`probit_flexBART()` approximates the function \\f(x)\\ with an ensemble
of binary regression trees. It also specifies independent priors on the
trees in the ensemble that are essentially identical to those deployed
by
[`flexBART()`](https://skdeshpande91.github.io/flexBART/reference/flexBART.md).
That is, the tree structure is generated using a branching process in
which the probability that a node at depth \\d\\ is non-terminal is
\\\alpha \times (1 + d)^{-\beta}\\. Then, decision rules are drawn
sequentially from the root down to each leaf. Finally, independent
\\N(\mu_0, \tau^2)\\ priors are specified for the outputs in each leaf.

With this specification, the marginal prior of any evaluation of the
regression function \\f(x)\\ is \\N(M \times \mu\_{0}, \tau^{2} \times
M)\\, where \\M\\ is the number of trees in the ensemble. Thus, for each
\\x\\, the induced prior for \\P(Y = 1 \vert X = x)\\ places 95%
probability on the interval \\\[\Phi(M \times \mu\_{0} - 2 \times \tau
\times \sqrt{M}), \Phi(M \times \mu\_{0} + 2 \times \tau \times
\sqrt{M})\]\\. By default, `probit_flexBART()` sets \\\tau =
1/\sqrt{M}\\ and \\\Phi^{-1}(\overline{y})/M\\ (i.e.,
`qnorm(mean(Y))/M`). Use the `mu0_vec` and `tau_vec` arguments to set
other hyperparameter values.

## Value

An object of `class` “flexBART” (essentially a list) containing

- dinfo:

  Essentially a list containing information about the input and output
  variables. Used by
  [`predict.flexBART`](https://skdeshpande91.github.io/flexBART/reference/predict.flexBART.md).

- trees:

  A list (or length `nd`) of character vectors (of lenght `M`)
  containing textual representations of the regression trees. These
  strings are parsed by
  [`predict.flexBART`](https://skdeshpande91.github.io/flexBART/reference/predict.flexBART.md)
  to reconstruct the C++ representations of the sampled trees.

- scaling_info:

  Essentially a list containing information for re-scaling raw MCMC
  output to the original outcome scale. Used by
  [`predict.flexBART`](https://skdeshpande91.github.io/flexBART/reference/predict.flexBART.md).

- M:

  A copy of the argument `M_vec`. Used by
  [`predict.flexBART()`](https://skdeshpande91.github.io/flexBART/reference/predict.flexBART.md).

- cov_ensm:

  An \\p \times R\\ binary matrix encoding whose (j,r)-element is 1 if
  trees in the ensemble for \\\beta\_{r}(X)\\ can split on \\X\_{j}\\.

- prob.train.mean:

  Vector containing posterior mean of \\P(Y=1 \vert X = x)\\ for the
  training data.

- prob.train:

  Matrix with `nd` rows and `length(Y_train)` columns containing
  posterior samples of \\P(Y=1 \vert X = x)\\ for the training data.
  Each row corresponds to a posterior sample of the regression
  functionand each column corresponds to a training observation. Only
  returned if `save_samples = TRUE`.

- prob.test.mean:

  Vector containing posterior mean of \\P(Y=1 \vert X = x)\\ on testing
  data, if testing data is provided.

- prob.test:

  If testing data was supplied, matrix containing posterior samples of
  the regression function evaluated on the testing data. Structure is
  similar to that of `prob.train`. Only returned if testing data is
  passed and `save_samples = TRUE`.

- varcounts:

  Matrix that counts the number of times a variable was used in a
  decision rule in each MCMC iteration. Structure is similar to that of
  `prob.train`, with rows corresponding to MCMC iteration and columns
  corresponding to predictors

- timing:

  Vector of runtimes for each chain

## See also

[`flexBART()`](https://skdeshpande91.github.io/flexBART/reference/flexBART.md)
for continuous outcomes.

## References

Albert, J.H. and Chib, S. (1993) Bayesian analysis of binary and
polychotomous data. *Journal of the American Statistical Association*.
**88**(422):669–679.
[doi:10.2307/2290350](https://doi.org/10.2307/2290350) .
