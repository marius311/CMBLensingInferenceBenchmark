# CMBLensingInferenceTestProblem.jl

A small wrapper around [CMBLensing.jl](https://github.com/marius311/CMBLensing.jl) which gives you a "black-box" posterior on which to test various inference algorithms. 

See [README.ipynb](README.ipynb) for complete example and a canonical HMC solution to compare against.

## Requirements

* Julia 1.9+

## Usage (working copy)

```
git clone https://github.com/marius311/CMBLensingInferenceTestProblem
cd CMBLensingInferenceTestProblem
julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'
```

Then open and run [README.ipynb](README.ipynb).

## Usage (as a dependency)

From the Julia environment for your project:

```
pkg> add https://github.com/marius311/CMBLensingInferenceTestProblem
```