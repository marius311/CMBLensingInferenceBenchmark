
module CMBLensingInferenceTestProblem

using CMBLensing
using ComponentArrays
using LogDensityProblems
using LinearAlgebra

export load_cmb_lensing_problem

struct CMBLensingLogDensityProblem
    ds
    Ωstart
    Ωtrue
    Λmass
end

function load_cmb_lensing_problem(;storage, T, Nside)

    (;f,  ϕ, ds) = load_sim(; storage, T, Nside, pol=:P, θpix=3, μKarcminT=1, beamFWHM=1, seed=0)
    (;f°, ϕ°) = mix(ds; f, ϕ)
    θ = ComponentVector(r=0.2, Aϕ=1)
    Ωtrue = LenseBasis(FieldTuple(;f°, ϕ°, θ=log.(θ)))

    Ωstart = let
        (;f, ϕ) = MAP_joint(ds, progress=false)
        (;f°, ϕ°) = mix(ds; f, ϕ)
        LenseBasis(FieldTuple(;f°, ϕ°, θ=log.(θ)))
    end

    Λmass = let
        (;Cf, Cϕ, Nϕ, D, G, B̂, M̂, Cn̂) = ds
        Diagonal(FieldTuple(
            f° = diag(pinv(D)^2 * (pinv(Cf) + B̂'*M̂'*pinv(Cn̂)*M̂*B̂)),
            ϕ° = diag(pinv(G)^2 * (pinv(Cϕ) + pinv(Nϕ))),
            θ = ComponentVector(r=1., Aϕ=1)
        ))
    end

    return CMBLensingLogDensityProblem(ds, Ωstart, Ωtrue, Λmass)
    
end

(prob::CMBLensingLogDensityProblem)(Ω::FieldTuple) = logpdf(Mixed(prob.ds); Ω.f°, Ω.ϕ°, θ=exp.(Ω.θ))
LogDensityProblems.logdensity(prob::CMBLensingLogDensityProblem, Ω::FieldTuple) = prob(Ω)
LogDensityProblems.capabilities(prob::CMBLensingLogDensityProblem) = 0
LogDensityProblems.dimension(prob::CMBLensingLogDensityProblem) = length(prob.Ωstart)

function to_from_vec(Ω::FieldTuple)
    to_vec(Ω) = Ω[:]
    from_vec(vec) = first(promote(vec, Ω))
    return to_vec, from_vec
end

end