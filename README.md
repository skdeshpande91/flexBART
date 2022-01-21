# flexBART


## Motivation
The original BART model works pretty well already, especially when you have continuous predictors.
But the way BART treats categorical variables -- creating a dummy variable for each level -- is a bit unsatisfying.

To elaborate on this, consider running BART with a singel categorical predictor X with levels {0,1,2,3}. 
By default, BART will create 3 dummary variables X1, X2, and X3, one each for the levels 1, 2, and 3.
It then assembles decision trees where decisions will be of the form ``If Xj = 0, go to left; else go to the right."
This decision rule corresponds to partitioning the set {0,1,2,3} into {0,2,3} (left) and {1}. 
More generally, default implementations of BART produce decision trees that correspond to recursively removing single elements from a set of levels.
Unfortunately, not every partition of the set of levels can be recovered using such a recursive partitioning strategy.

To wit, we cannot form the partition {0,1}, {2,3} using the prior process used in default implementations of BART.
More formally, I conjecture that the default prior of BART lacks *full support* in the sense that the tree prior places 0 prior probability on certain partitions of the set of levels of a categorical predictor.

From a practical standpoint, this is probably much ado about nothing. After all, BART doesn't just use a single tree; it builds an entire ensemble of them.
And if one really has a response surface where there is some similiarities across categorical levels 0 & 1 and some different similarity across levels 2 & 3, BART'll probably tease it out.

Nevertheless, it is still interesting to see if we can introduce some simple modifications of BART's tree prior that do have full support.

Another downside of the default approach (i.e. creating dummy variables and splitting on individual values of them) occurs when you have some prior structural knowledge about the levels.
This is most salient in spatial settings where a categorical variable can be used to represent some geographical unit (e.g. state, census tract, block) and there is a natural adjacency structure between the units.
It might be desirable to ensure that when BART splits on such a structure categorical variable, the partition respects the assumed adjacency structure.

Finally, BART relies on axis-aligned splits of continuous variables. Several authors have reported that building decision trees with random combinations of features can achieve some marginal improvements in predictive accuracy. It is not clear whether a similar phenomenon happens with BART.

## Towards a more flexible BART

`flexBART` is a re-implementation of BART that tries to produce more thoughtful splits on categorical variables and enable investigation of using non-axis aligned splitting rules based on random combinations of continuous predictors.
While `flexBART` is largely based on the design principles in `BART` package but contains a couple of improvements related to speed. 
For instance, in the main MCMC loop of `BART::wbart()`, multiple tree traversals are carried out, which can be expensive when there are lots of observations. In `flexBART`, no tree traversals are performed. 

