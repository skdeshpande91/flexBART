# flexBART


## Motivation
Default implementations of Bayesian Additive Regression Trees (BART) represent categorical predictors using several binary indicators, one for each level of each categorical predictor. Axis-aligned decision rules are well-defined with these indicators; they send one level of a categorical predictor to the left and all others levels to the right (or vice versa).
Regression trees built with these rules partition the set of all levels of a categorical predictor by recursively removing one level at a time. 
Unfortunately, most partitions of the levels cannot be built with this ``remove one at a time'' strategy, meaning that default implementations of BART are extremely limited in their ability to ``borrow strength'' across groups of levels.

The **flexBART** package overcomes this limitation by utilizing a new prior for decision trees.
Like other implementations, the drawing a decision tree from the new prior is accomplished by first simulating a branching process and then randomly drawing decision rules for each non-terminal (i.e. non-leaf) node of the tree.
Decision rules are drawn in two steps. First a variable is selected (either uniformly at random or from a vector of probabilities that itself is given a Dirichlet prior; see Linero (2018) for details).
If this variable is a continuous predictor, a single cut-point is drawn uniformly from the set of available cut-points.
However, if the splitting variable is a categorical predictor, the new prior assigns each available levels of that predictor randomly to the left or right branch.
In this way, the prior over decision tree now allows multiple levels of a categorical variable to be assigned to both the left and right branches of an internal node.

Building on this, the package also provides support for *structured* categorical predictors for which there are a prior preferences about which levels of the predictor ought to be clustered together.
These preferences are operationalized with a network whose vertices correspond to the levels of the predictor and whose edges encode co-clustering preference.
An example would be spatial regions where each region is represented by a vertex in a network and an edge is drawn between vertices whose corresponding regions are geographically/spatially adjacent. We would like to get the decision tree prior to respect the supplied adjacency information.

## Installation and basic usage


## Notes about the implementation

**flexBART** is a re-implementation of BART that tries to produce more thoughtful splits on categorical variables.
While **flexBART** is largely based on the design principles in **BART** package, it contains a couple of improvements designed to make the code more readable and faster.
By far the most salient is that in the main MCMC loop we no longer perform any tree traversals; that is, we do not loop over all of the observations and trace each observation's path from the root node of a tree to one of the leafs.
Instead, we keep track of the partition of observations induced by each tree and update them as needed in the Metropolis-Hastings step.
In the code, we represent the partition as a `std::map<int,std::vector<int>>` where the key is the node id of the leaf and the value is a vector holding the observation.
