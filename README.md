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
As a sort of preview, here are a few prior draws from this new regression tree prior.

When we put them all together, we end up with some neat looking spatial functions! 


Of course, we can extend this to create a cute spatio-temporal process using an ensemble of trees that can partition space and time. 
Here are a few prior draws when we allow ourselves to split on space & time (space once again being encoded by census tract id). 



Now that we're thinking about more exotic decision rules, why stop with categorical variables? BART traditionally relies on axis-aligned decision rules of continuous variables.
But several authors have, at various points since the introduction of random forests, reported moderate gains when they allow decision rules based on random combinations of features.
Is the same true for BART?
This new implementation also allows for so-called ``random combination'' decision rules.

### Towards a more flexible BART

`flexBART` is a re-implementation of BART that tries to produce more thoughtful splits on categorical variables and enable investigation of using non-axis aligned splitting rules based on random combinations of continuous predictors.
While `flexBART` is largely based on the design principles in `BART` package, it contains a couple of improvements designed to make the code more readable and faster.
Perhaps the biggest change is that in the main MCMC loop we no longer perform any tree traversals. 

