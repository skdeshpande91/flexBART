# flexBART
  <!-- badges: start -->
  [![R-CMD-check](https://github.com/skdeshpande91/flexBART/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/skdeshpande91/flexBART/actions/workflows/R-CMD-check.yaml)
  [![CRAN status](https://www.r-pkg.org/badges/version/flexBART)](https://CRAN.R-project.org/package=flexBART)
  <!-- badges: end -->

Welcome to version 2.0 of the **flexBART** package!
**flexBART** (>= 2.0.0) is a new implementation of BART that is designed to fit flexible varying coefficient models using ensembles of binary regression trees.
In addition to the flexible priors for categorical decision rules introduced in earlier versions, this new version introduces a formula interface and implements a lot of data pre-processing that (hopefully) makes it easier than ever to fit BART models.

## Installation & Basic Usage

It is highly recommended that you install R version 4.0.0 or later before installing **flexBART**.
Before installing **flexBART**, ensure that you have set up an appropriate C++ toolchain for your system.

  * For macOS: we recommend using the [macrtools](https://github.com/coatless-mac/macrtools) package
  * For Windows: we recommend using Rtools, which can be downloaded [here](https://cran.r-project.org/bin/windows/Rtools/). Please make sure you download the version of Rtools that corresponds to your R version (e.g., RTools45 for R version 4.5.x) 
  * For Linux: we recommend following [these instructions](https://github.com/stan-dev/rstan/wiki/Configuring-C-Toolchain-for-Linux) from the Stan development team.

Once your C++ toolchain is configured, you can install **flexBART** using `devtools::install_github`:
```
devtools::install_github(repo = "skdeshpande91/flexBART")
```

### Basic Usage

Starting in version 2.0.0, **flexBART** features a formula interface and allows users to pass their data as `data.frame` or `tibble` objects.
So, given a data frame `train_data` containing named columns for an outcome (e.g., `Y`) and predictors, you can fit a simple BART model to predict `Y` using all the predictors by running
```
flexBART(formula = Y ~ bart(.), train_data = train_data)
```

**flexBART** also supports fitting varying coefficient models of the form

$$
Y = \beta_{0}(X) + \beta_{1}(X)Z_{1} + \cdots + \beta_{R}Z_{R} + \sigma \epsilon; \epsilon \sim N(0,1),
$$

where each coefficient function $\beta_{r}(X)$ is approximated with its own tree ensemble.
To fit such a model in **flexBART**, you can use a formula like `Y ~ bart(.) + Z1 * bart(.) + Z2 * bart(.)`, including a separate `bart()` for each coefficient function.

The formula interface also provides fine control over the predictor variables used in each ensemble.
To allow an ensemble to only split on a few variables (e.g., `X1`, `X2`, and `X3`), you would specify `bart(X1 + X2 + X3)` and to allow an ensemble to split on all variables except `X1` and `X2`, you would specify `bart(.-X1-X2)`. 
Note that when it detects multiple ensembles in the formula, **flexBART** will *not* include any of the $Z_{r}$'s as splitting variables when it expands the `.`
So, to include, say, a piece-wise linear function, $X_{1} * \beta_{1}(X_{1}),$ you would need to specify `X1 * bart(X1)` in the formula argument. 

See the package articles at [the package website](https://skdeshpande91.github.io/flexBART) for more details.

<!--
## Pre-processing

Like earlier version (e.g., 1.2.0 and earlier), the latest version of **flexBART** assumes that all continuous predictors are re-scaled to the interval [-1,1] and represents the distinct values of categorical predictors with non-negative integers.
But unlike those earlier versions, which required users to perform such re-scaling and conversion themselves, **flexBART** now automates the pre-processing.


### Manually specifying cutpoints for numeric predictors

Internally, **flexBART** treats all predictors passed as a `factor` or `character` as categorical.
It then checks whether each numerical predictor is discrete (e.g., age measured in years) or whether it is continuous by looking at the number of pairwise differences between consecutive values.
Decision rules based on numerical predictors take the form $\{X_{j} < c\}.$

If **flexBART** detects that $X_{j}$ is continuous, it will rescale the supplied values of $X_{j}$ to the interval [-1,1] and allow regression trees to select the cutpoint $c$ uniformly from that interval.
**flexBART** adds $0.1\text{sd}(X_{j})$ to the maximum value of $X_{j}$ and subtracts $0.1\text{sd}(X_{j})$ from the minimum value of $X_{j}$ before re-scaling the predictor.
If testing data is provided, **flexBART** determines the min, max, and standard deviation of $X_{j}$ using both the training and testing data.

If, on the other hand, **flexBART** determines that $X_{j}$ is discrete, it will not re-scale the predictor and instead forces regression trees to select the cutpoint $c$ from the unique values of $X_{j}.$
If testing data is provided, **flexBART** determines the unique values of $X_{j}$ using both the training and testing data.

### Priors for categorical predictors

In **flexBART** ensembles, decision rules based on categorical predictors take the form $\{X_{j} \in \mathcal{C}\}$ where $\mathcal{C}$ is a random subset of the discrete values that $X_{j}$ can assume.
This is in stark contrast to most other implementations of BART, which one-hot encode categorical predictors.
Please see [Deshpande (2024)](https://doi.org/10.1080/10618600.2024.2431072) for arguments against the use of one-hot encoding with BART.

Internally, **flexBART** determines the set of available values of $X_{j}$ by looking at the `levels()` of all predictors saved as `factor()` variables.
As a result, **flexBART** is able to make predictions at new values of a categorical predictor not present in the training data *so long as these values are included as levels of that predictor.*


**flexBART** also includes support for network-structured categorical predictors (e.g., spatial areas with known adjacency structure). 
To force the "cutset" $\mathcal{C}$ to correspond to connect components of these networks, you should provide the corresponding adjacency matrices via the `adjacency_list` argument.
This argument should be a named list with one element per network-structured predictor.
Each element should be a binary or weighted adjacency matrix whose row and column names correspond to the levels of the predictor.
**flexBART** implements four different priors over decision rues for network-structured predictors. See the documentation and Section 3.2 of [Deshpande (2024)](https://doi.org/10.1080/10618600.2024.2431072) for details about these priors.
-->


<!--
The **flexBART** package overcomes this limitation by utilizing a new prior for decision trees.
Like other implementations, the drawing a decision tree from the new prior is accomplished by first simulating a branching process and then randomly drawing decision rules for each non-terminal (i.e. non-leaf) node of the tree.
Decision rules are drawn in two steps. First a variable is selected (either uniformly at random or from a vector of probabilities that itself is given a Dirichlet prior; see Linero (2018) for details).
If this variable is a continuous predictor, a single cut-point is drawn uniformly from the set of available cut-points.
However, if the splitting variable is a categorical predictor, the new prior assigns each available levels of that predictor randomly to the left or right branch.
In this way, the prior over decision tree now allows multiple levels of a categorical variable to be assigned to both the left and right branches of an internal node.

Building on this, the package also provides support for *structured* categorical predictors for which there are *a priori* preferences about which levels of the predictor ought to be clustered together.
These preferences are operationalized with a network whose vertices correspond to the levels of the predictor and whose edges encode co-clustering preference.
An example would be spatial regions where each region is represented by a vertex in a network and an edge is drawn between vertices whose corresponding regions are geographically/spatially adjacent. We would like to get the decision tree prior to respect the supplied adjacency information.

## Installation and basic usage

The package source files are contained in the sub-directory flexBART.
To install, you can either download that directory and build and install the package from the command line.
Alternatively, you can install using `devtools::install_github`:
```
devtools::install_github(repo = "skdeshpande91/flexBART", subdir = "flexBART")
```

The `examples` subdirectory contains some case studies showing how to use flexBART with network data and comparing the run-time to the **BART** package. 


## Notes about the implementation

**flexBART** is a re-implementation of BART that tries to produce more thoughtful splits on categorical variables.
While **flexBART** is largely based on the design principles in **BART** package, it contains a couple of improvements designed to make the code more readable and faster.
By far the most salient is that in the main MCMC loop we no longer perform any tree traversals; that is, we do not loop over all of the observations and trace each observation's path from the root node of a tree to one of the leafs.
Instead, we keep track of the partition of observations induced by each tree and update them as needed in the Metropolis-Hastings step.
In the code, we represent the partition as a `std::map<int,std::vector<int>>` where the key is the node id of the leaf and the value is a vector holding the observation.
-->
