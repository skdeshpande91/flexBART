# Predicting new observations with previously fitted flexBART model

`predict.flexBART()` can take the output of
[`flexBART`](https://skdeshpande91.github.io/flexBART/reference/flexBART.md)
and
[`probit_flexBART`](https://skdeshpande91.github.io/flexBART/reference/probit_flexBART.md)
and use it to make predictions at new inputs.

## Usage

``` r
# S3 method for class 'flexBART'
predict(object, newdata, ...)
```

## Arguments

- object:

  object of class inheriting from “flexBART”.

- newdata:

  Data frame in which to look for variables with which to predict.
  Cannot be omitted.

- ...:

  Additional optional arguments governing whether to print progress
  (i.e., `verbose` and `print_every`).

## Details

Make predictions at new inputs based on the output of
[`flexBART`](https://skdeshpande91.github.io/flexBART/reference/flexBART.md).
When the training and/or testing dataset is large, it is recommended to
run
[`flexBART`](https://skdeshpande91.github.io/flexBART/reference/flexBART.md)
with the arguments `save_samples = FALSE` and `save_trees = TRUE`. Then,
to access posterior samples of the function evaluations, pass the fitted
“flexBART” object to `predict.flexBART()` along with the appropriate
`newdata` argument. If `fit` were produced by
[`probit_flexBART`](https://skdeshpande91.github.io/flexBART/reference/probit_flexBART.md),
then the function outputs draws of the fitted probabilities.

## Value

When there is only one ensemble, a matrix containing posterior samples
of the regression function evaluated at the supplied inputs. For models
with multiple ensembles, a list with three elements: (i) a matrix
containing posterior samples of the regression function evaluation and
(ii) an array containing evaluations of the identified slopes; and (iii)
an array containing evaluations of all slopes on the standardized scale.
