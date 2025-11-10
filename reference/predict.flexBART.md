# Predicting new observations with previously fitted flexBART model

`predict.flexBART()` can take the output of
[`flexBART()`](https://skdeshpande91.github.io/flexBART/reference/flexBART.md)
and
[`probit_flexBART()`](https://skdeshpande91.github.io/flexBART/reference/probit_flexBART.md)
and use it to make predictions at new inputs.

## Usage

``` r
# S3 method for class 'flexBART'
predict(object, newdata, ...)
```

## Arguments

- object:

  object of class inhereting from “flexBART”.

- newdata:

  Data frame in whch to look for variables with which to predict. Cannot
  be omitted.

- ...:

  Additional optional arguments governing whether to print progress
  (i.e., `verbose` and `print_every`).

## Details

Make predictions at new inputs based on the output of
[`flexBART()`](https://skdeshpande91.github.io/flexBART/reference/flexBART.md).
Useful when the testing dataset is quite large. If `fit` were produced
by
[`probit_flexBART()`](https://skdeshpande91.github.io/flexBART/reference/probit_flexBART.md),
then the function outputs draws of the fitted probabilities.

## Value

When there is only one ensemble, a matrix containing posterior samples
of the regression function evaluated at the supplied inputs. For models
with multiple ensembles, a list with three elements: (i) a matrix
containing posterior samples of the regression function evaluation and
(ii) an array containing evaluations of the identified slopes; and (iii)
an array containing evaluations of all slopes on the standardized scale.
