# Under the Hood

Most pertinent details of **DiDInt.jl** can be found in the [DID-INT paper](https://doi.org/10.48550/arXiv.2412.14447), however, some details that might be of particular interest to users of the **DiDInt.jl** package are elaborated upon here.

## Four Step Estimation Procedure

The **DID-INT** estimator is computed in four steps:

1. The means of the outcome of interest (conditional on covariates) are computed for each period for each state (denoted as $\lambda$ in the DID-INT paper)
2. Within state differences between these $\lambda$ values are computed and are labelled as coming from treated or untreated states
3. A (weighted) regression of differences on an intercept term and a dummy variable indicating treated/untreated status is then used to compute subaggregate ATTs for a given subset of the data (more details on the choice of subset in the aggregation section)
4. The (weighted) mean of subaggregate ATTs is taken in order to compute the aggregate ATT

More explicitly:

1. 


## Weighting

## Aggregation

## Randomization Inference

## Period Grid Construction & Date Matching Procedure