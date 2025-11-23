# DiDInt.jl

Welcome to the DiDInt.jl documentation!

**DiDInt.jl** is a package for implementing intersection difference-in-differences (DID-INT), a method developed by [Karim & Webb (2025)](https://doi.org/10.48550/arXiv.2412.14447). **DiDInt.jl** allows for unbiased estimation of the average effect of treatment on the treated (ATT) in cases when the common causal covariates (CCC) assumption is violated. That is, it allows for unbiased estimation of the ATT when the effects of covariates on the outcome of interest vary by state, time, or both. Further, **DiDInt.jl** works for both staggered adoption (assuming that once treated, always treated) and common adoption scenarios.

**DiDInt.jl** introduces two functions: one for estimation (`didint()`, and one for plotting `didint_plot()`). The details of these functions can be found [here](functions.md).

Some further details (a look under the hood) can be found on this [page](details.md). Topics covered on this page include everything ranging from weighting options, aggregation methods, computation of p-values via the randomization inference procedure described in [MacKinnon and Webb (2020)](https://doi.org/10.1016%2Fj.jeconom.2020.04.024), and how **DiDInt.jl** processes dates and builds periods. 

Examples can be found [here](examples.md)!