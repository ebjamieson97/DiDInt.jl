# DiDInt.jl
Intersection difference-in-differences!

## Installation 
```
using Pkg
Pkg.add(url="https://github.com/ebjamieson97/DiDInt.jl")
```

## Usage

The core function didint() estimates the ATT while adjusting for covariates that may vary by state, time, or both. Its interface is defined as follows:
```
didint(outcome::AbstractString, state::AbstractString,  
       time::AbstractString, data::DataFrame,  
       treated_states::Union{Vector{<:AbstractString}, AbstractString},  
       treatment_times::Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date};  
       date_format::Union{AbstractString, Nothing} = nothing,  
       covariates::Union{Vector{<:AbstractString}, AbstractString, Nothing} = nothing,  
       ccc::AbstractString = "int", agg::AbstractString = "state",  
       ref::Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing,  
       freq::Union{AbstractString, Nothing} = nothing, freq_multiplier::Number = 1,  
       autoadjust::Bool = false, nperm::Number = 1000, verbose::Bool = true)
``` 

### Parameters

- **outcome** (AbstractString):  
  Name of the column identifying the outcome of interest.

- **state** (AbstractString):  
  Name of the column identifying the state membership.

- **time** (AbstractString):  
  Name of the column identifying the date of the observation.

- **data** (DataFrame):  
  The DataFrame containing the data for analysis.

- **treated_states** (Union{Vector{<:AbstractString}, AbstractString}):  
  A vector (or single string) identifying the treated state(s).

- **treatment_times** (Union{T, Vector{T}} where T <: Union{AbstractString, Number, Date}):  
  A vector (or single entry) of treatment times. The order should match that of treated_states (i.e. the first treated state corresponds to the first treatment time, and so on).

- **covariates** (Union{Vector{<:AbstractString}, AbstractString, Nothing}):  
  Optional covariates as a vector (or single string) of column names, or nothing (default).

- **ccc** (AbstractString = "int"):  
  Specify which version of DID-INT to use. Options: "hom", "time", "state", "add", "int" (default).

- **agg** (AbstractString = "state"):  
  Weighting method for aggregation. Options include "cohort", "simple", "state", and "unweighted".

- **ref** (Union{Dict{<:AbstractString, <:AbstractString}, Nothing} = nothing):  
  Dictionary specifying the reference category for any categorical variable.

- **freq** (Union{AbstractString, Nothing} = nothing):  
  Timeframe for periods in staggered adoption scenarios (e.g., "year", "month", "week", "day").

- **freq_multiplier** (Number = 1):  
  A multiplier for the freq argument (e.g., set to 2 for a two-year period).

- **autoadjust** (Bool = false):  
  Automatically determine the period length in staggered adoption settings.

- **nperm** (Number = 1000):  
  Number of unique permutations to use for randomization inference.

- **verbose** (Bool = true):  
  Display progress during computation.
