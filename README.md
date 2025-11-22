# DiDInt.jl
Intersection difference-in-differences!

**DiDInt.jl** introduces two functions:
- `didint()`, which is used for estimating the average effect of treatment on the treated (ATT) while accounting for various violations of the common causal covariates (CCC) assumption
- `didint_plot()` which produces datasets that can easily be used to make parallel trends plots that help determine which CCC violation should be accounted for when using `didint()`. `didint_plot` can also produce datasets that can easily be used to make event study plots.

To learn more about the common causal covariates assumption and intersection difference-in-differences, please see [Karim & Webb (2025)](https://arxiv.org/abs/2412.14447).

A Stata wrapper for this package is available here: **[didintjl](https://github.com/ebjamieson97/didintjl)**.

## Installation 
```
using Pkg
Pkg.add(url="https://github.com/ebjamieson97/DiDInt.jl")
```

## `didint()` - Estimate ATT

The function `didint()` estimates the ATT while adjusting for covariates that may vary by state, time, or both.
```julia
didint(
       outcome::Union{AbstractString, Symbol},
       state::Union{AbstractString, Symbol},
       time::Union{AbstractString, Symbol},
       data::DataFrame;
       gvar::Union{AbstractString, Symbol, Nothing} = nothing,
       treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing,
       treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date, Nothing} = nothing,
       date_format::Union{AbstractString, Nothing} = nothing,
       covariates::Union{T, Vector{T}} where T <: Union{AbstractString, Symbol} = nothing,
       ccc::AbstractString = "int",
       agg::AbstractString = "cohort",
       weighting::AbstractString = "both",
       ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
       freq::Union{AbstractString, Nothing} = nothing,
       freq_multiplier::Number = 1,
       start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
       end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
       nperm::Number = 999,
       verbose::Bool = true,
       seed::Number = rand(1:1000000),
       use_pre_controls::Bool = false,
       notyet::Union{Nothing, Bool} = nothing,
       hc::Union{AbstractString, Number} = "hc3"
      )
``` 

### Example
```julia
using DiDInt
# Placeholder for now, sorry
```

### Parameters

#### Required Parameters
- **outcome**
  Name of the column which identifies the outcome of interest.
- **state**
  Name of the column which identifies the state membership of the observation.
- **time**
  Name of the column which identifies the date of the observation.
- **data** 
  The DataFrame to be used for the analysis.

#### Treatment Specification
- **gvar**  
  Name of the column which indicates time of first treatment for each state.
- **treated_states**  
  A vector of strings (or a single string) noting the treated state(s).
- **treatment_times**  
  A vector (or single entry) denoting the associated treatment times of the treated_states. The order should match `treated_states` (i.e., the first treated state corresponds to the first treatment time, and so on).

#### Model Specification
- **ccc**  
  Specify which version of DID-INT should be used. Options: `"hom"`, `"time"`, `"state"`, `"add"`, `"int"`.
- **agg**  
  Enter the aggregation method as a string. Options: `"cohort"`, `"simple"`, `"state"`, `"sgt"`, `"none"`.
- **weighting**  
  Specify which weighting method should be used. Options: `"both"`, `"att"`, `"diff"`, `"none"`.
- **covariates**  
  A vector of covariates entered as either strings or symbols (or a single covariate string or symbol), or `nothing` (default).
- **notyet**  
  Determine if pre-treatment periods from treated states should be used as controls.
- **ref**  
  A dictionary specifying which category in a categorical variable should be used as the reference (baseline) category.

#### Date Processing & Period Grid Construction
- **date_format**  
  Date format (e.g., `"yyyy"` or `"yyyy-mm-dd"`) to be used when parsing string dates from the time column, or `start_date`, `end_date`, and `treatment_times` arguments.
- **freq**  
  A string indicating the desired timeframe of a period for the analysis for staggered adoption scenarios. Options: `"year"`, `"month"`, `"week"`, `"day"`.
- **freq_multiplier**  
  An integer by which the `freq` argument should be multiplied in a staggered adoption scenario (e.g., if a two-year period is desired, set `freq = "year"` and `freq_multiplier = 2`).
- **start_date**  
  Any data prior to this date is dropped, and serves as the starting date for the period grid construction if activated.
- **end_date**  
  Any data after this date is dropped, and serves as the end date for the period grid construction if activated.

#### Inference
- **nperm**  
  The number of unique permutations (not including the initial assignment of treatment times) to be considered when performing the randomization inference.
- **verbose**  
  A boolean option for displaying progress of the randomization procedure.
- **seed**  
  An integer to set the random seed for the randomization inference procedure.
- **hc**  
  Specify which heteroskedasticity-consistent covariance matrix estimator (HCCME) should be used. Options: `0`, `1`, `2`, `3`, `4` (or `"hc0"`, `"hc1"`, `"hc2"`, `"hc3"`, `"hc4"`).

### Returns
A DataFrame of results including the estimate of the ATT as well as standard errors and p-values.

## `didint_plot()` - Prepare Data for Visualization

The `didint_plot()` function produces a dataset in a long format that can easily be used for plotting parallel trends or event study plots.
```julia
didint_plot(
            outcome::Union{AbstractString, Symbol},
            state::Union{AbstractString, Symbol},
            time::Union{AbstractString, Symbol},
            data::DataFrame;
            gvar::Union{AbstractString, Symbol, Nothing} = nothing,
            treated_states::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Nothing} = nothing,
            treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date, Nothing} = nothing,
            date_format::Union{AbstractString, Nothing} = nothing,
            covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing} = nothing,
            ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,
            ccc::Union{AbstractString, Vector{<:AbstractString}} = "all",
            event::Bool = false,
            weights::Bool = true,
            ci::Number = 0.95,
            freq::Union{AbstractString, Nothing} = nothing,
            freq_multiplier::Number = 1,
            start_date::Union{AbstractString, Number, Date, Nothing} = nothing,
            end_date::Union{AbstractString, Number, Date, Nothing} = nothing,
            hc::Union{AbstractString, Number} = "hc3"
           )
```

### Example
```julia
using DiDInt
# Placeholder for now, sorry
```

### Parameters

#### Required Parameters
- **outcome**  
  Input the name of the column which identifies the outcome of interest.
- **state**  
  Input the name of the column which identifies the state membership of the observation.
- **time**  
  Input the name of the column which identifies the date of the observation.
- **data**  
  The DataFrame to be used for the analysis.

#### Treatment Specification
- **gvar**  
  Name of the column which indicates time of first treatment for each state.
- **treatment_times**  
  A vector (or single entry) denoting the associated treatment times of the treated_states. The order should match `treated_states` (i.e., the first treated state corresponds to the first treatment time, and so on).

#### Model Specifications
- **ccc**  
  Specify which versions of DID-INT should be used. Options are either `"all"`, or any combination of: `"none"`, `"hom"`, `"time"`, `"state"`, `"add"`, `"int"`.
- **covariates**  
  A vector of covariates entered as either strings or symbols (or a single covariate string or symbol), or `nothing` (default).
- **ref**  
  A dictionary specifying which category in a categorical variable should be used as the reference (baseline) category.

#### Event Study Plot Options
- **event**  
  Specify if data should be prepared for an event study plot as opposed to a parallel trends plot.
- **treated_states**  
  A vector of strings (or a single string) noting the treated state(s).
- **weights**  
  Whether to use weighted means when computing event study estimates. If `true`, estimates are computed as weighted averages of state-level means for each period relative to treatment; if `false`, uses simple unweighted averages.
- **ci**  
  Define the size of confidence bands for the event study plot.
- **hc**  
  Specify which heteroskedasticity-consistent covariance matrix estimator (HCCME) should be used. Options: `0`, `1`, `2`, `3`, `4` (or `"hc0"`, `"hc1"`, `"hc2"`, `"hc3"`, `"hc4"`).

#### Date Processing & Period Grid Construction
- **date_format**  
  Date format (e.g., `"yyyy"` or `"yyyy-mm-dd"`) to be used when parsing string dates from the time column, or `start_date`, `end_date`, and `treatment_times` arguments.
- **freq**  
  A string indicating the desired timeframe of a period for the analysis for staggered adoption scenarios. Options: `"year"`, `"month"`, `"week"`, `"day"`.
- **freq_multiplier**  
  An integer by which the `freq` argument should be multiplied in a staggered adoption scenario (e.g., if a two-year period is desired, set `freq = "year"` and `freq_multiplier = 2`).
- **start_date**  
  Any data prior to this date is dropped, and serves as the starting date for the period grid construction if activated.
- **end_date**  
  Any data after this date is dropped, and serves as the end date for the period grid construction if activated.

### Returns
A DataFrame of means and means residualized by the specified covariates for each of the specified common causal covariates (CCC) violations by period for each state, or, a DataFrame of means of the treated states by periods before/after treatment (again, residualized by the specified covariates and for each of the specified CCC violations).