# flexBART


## Motivation
The original BART model works pretty well already, especially when you have continuous predictors.
But the way BART treats categorical variables -- creating a dummy variable for each level -- is a little bit unsatisfying.


Consider running BART with a single categorical predictor X that has levels {0, 1, 2, 3}.
By default, BART creates 3 dummy variables, one each for the levels 1, 2, and 3. Let's call these variables X1, X2, and X3 where X1 = 1 if and only if X = 1 and X1 = 0 otherwise.
BART then assembles decision trees where each decision rule be of the form ``If Xj = 0, go to the left; else go to the right.''
By tracing the paths from the root node to each leaf, we can construct a partition over the set of levels of X. 
Each step of this construction involves peeling away a single level from a larger set of levels.
It is not difficult to see that we cannot produce the partition {0,2}, {1,3} using this construction: the very first decision rule (at the root node) would produce a singleton that cannot be merged with any other subsequent set of levels.
This argument suggests the following conjecture: the default implementation of the BART prior lacks *full support* in the sense that it places 0 prior probability on certain partitions of the set of levels of a categorical predictor.

In this simple example, the default BART would never allow a tree in the ensemble to create two groups of observations, one containing all observations with X = 0 or 2 and one containing all observations with X = 1 or 3.
In other words, the default BART is precluded from trying out a wide class of ways to ``borrow strength'' across categorical levels.
It would be interesting (if only to me) to see whether a more flexible BART that was allowed to try out a richer set of strategies for pooling observations across categorical levels offered any real practical benefits.


Once we allow categorical decision rules that cluster multiple levels together simultaneously, the next natural extension is to permit decision rules that obey certain structure.
For instance, in analyzing aggregated spatial data (e.g. at the census tract or block level), we might use a categorical variable to encode the spatial unit. 
If we have a known adjacency structure, it would be nice to restrict our decision trees to partition spatial units in a spatially contiguous fashion.
This new implementation of BART allows for this possibility.

Now that we're thinking about more exotic decision rules, why stop with categorical variables? BART traditionally relies on axis-aligned decision rules of continuous variables.
But several authors have, at various points since the introduction of random forests, reported moderate gains when they allow decision rules based on random combinations of features.
Is the same true for BART?
This new implementation also allows for so-called ``random combination'' decision rules.
But early tests haven't been super promising (more to come!)

### Towards a more flexible BART

`flexBART` is a re-implementation of BART that tries to produce more thoughtful splits on categorical variables and enable investigation of using non-axis aligned splitting rules based on random combinations of continuous predictors.
While `flexBART` is largely based on the design principles in `BART` package, it contains a couple of improvements designed to make the code more readable and faster.
By far the most salient is that in the main MCMC loop we no longer perform any tree traversals; that is, we do not loop over all of the observations and trace each observation's path from the root node of a tree to one of the leafs.
Instead, we keep track of the partition of observations induced by each tree and update them as needed in the Metropolis-Hastings step.
In the code, we represent the partition as a `std::map<int,std::vector<int>>` where the key is the node id of the leaf and the value is a vector holding the observation.

### Example

The file `scripts/grid_example.R` is a small spatial example. 
In this example, we have 100 spatial units arranged in a 10 x 10 grid (nodes in the graph represent spatial units and edges represent adjacency relations).
We then defined a piecewise constant function mu on this grid by first creating 5 clusters of spatial units and assigning a constant value within each cluster.
We generated a single observation at each node by adding standard normal noise to value of mu.


We then randomly selected 10 nodes to be our testing data and training flexBART and BART on the remaining 90 datapoints. In the figure below, the nodes in the test set are colored in gray.
Both flexBART and BART were aware of the full adjacency structure in training; that is, although they did not encounter data from the testing nodes, we did not modify the underlying adjacency matrix of the nodes. 
We compared how well flexBART (which can exploit adjacency information) and regular BART (which cannot) were able to (i) recover F evaluate at each spatial unit in the training dataset and (ii) interpolate the value of F at each spatial unit in the testing dataset. 

![](https://github.com/skdeshpande91/flexBART/blob/main/figures/grid_example.png")

In the figure, the color scale runs from about -9 (dark blue) to 0 (beige/yellow-ish) to about 9 (dark red). 
The same scale is used for all panels.

flexBART (w/o adjacency) uses decision rules that partition levels of the categorical variable uniformly, while flexBART (w/ adjacency) uses decision rules that partition levels in a way that respects the underlying adjacency structure.
Regular BART does pretty poorly on this example. Essentially it pools together all of the node and makes almost the same prediction at each. While this behavior is consistent with all trees in the ensemble being equal to a root node (i.e. a stump), we found that the chain visited non-trivial trees. 



