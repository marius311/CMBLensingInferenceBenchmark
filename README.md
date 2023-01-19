# CMBLensingInferenceTestProblem.jl

A small wrapper around [CMBLensing.jl](https://github.com/marius311/CMBLensing.jl) which gives you a "black-box" posterior on which to test various inference algorithms. 

See [README.ipynb](README.ipynb) for complete example and a canonical solution to compare against. 

## Usage

**Requires Julia 1.9 beta**

Clone this repo, instantiate the package environment, then open and run README.ipynb. 

## Note

If instead you want to install CMBLensingInferenceTestProblem.jl itself into a parent environment, make sure to also install the exact versions of CMBLensing, JLD2, and Memoization which are pinned in the repo's Manifest.toml. To see these versions, run `pkg> st` from this repo, then in the parent environment run e.g. `pkg> add CMBLensing#xxxxxxx` with the appropriate commit hash.

(can get rid of this note once I get all the necessary changes merged upstream)