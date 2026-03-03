# flexBART 2.0.5

## Major changes
  * Refactor internal function to detect nesting structure, which previously threw an error if it encountered a categorical variable that was incompletely observed (i.e., not every level appeared in the supplied data)
  * Fixes error in how heteroskedastic BART saves samples of the residual standard deviation function with test data.


# flexBART 2.0.0

## Major changes
  * Add support for varying coefficient models (i.e., multipe additve ensembles)
  * Remove unnecessary passes over the full training data during backfitting and computing conditional posterior of leaf nodes
  * Introduce pre-processing functions
  * Add support for nested categorical features
  * Introduce formula interface

# flexBART 1.0.0

## Major changes

  * Initial package
