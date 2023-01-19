
@time using CMBLensing, ComponentArrays, LinearAlgebra, IterativeSolvers, 
    LogDensityProblems, AdvancedHMC, Zygote, BenchmarkTools, PyPlot, ProgressMeter, Random, NamedTupleTools

##
storage = Array
(;f, ϕ, ds) = @time load_sim(; storage, T=Float64, Nside=64, pol=:P, θpix=3, μKarcminT=1, beamFWHM=1, seed=0)
(f°, ϕ°) = mix(ds; f, ϕ)
θ = ComponentVector(r=0.1, Aϕ=1)
##
Ω₀ = LenseBasis(FieldTuple(select(mix(ds; select(simulate(Xoshiro(1), ds), (:f, :ϕ))...), (:f°, :ϕ°))))
Ωvec₀ = Ω₀[:]
##
(;Cf, Cϕ, Nϕ, D, G, B̂, M̂, Cn̂) = ds;
Λmass = Diagonal(FieldTuple(
    f° = diag(pinv(D)^2 * (pinv(Cf) + B̂'*M̂'*pinv(Cn̂)*M̂*B̂)),
    ϕ° = diag(pinv(G)^2 * (pinv(Cϕ) + pinv(Nϕ)))
))
##
struct CMBProblem
    ds
    Ω₀
    Λmass
end
# function LogDensityProblems.logdensity(p::CMBProblem, zvec)
#     (;ds, Ω₀) = p
#     (;f°, ϕ°, θ°) = copy(first(promote(zvec, Ω₀)))
#     logpdf(Mixed(ds); f°, ϕ°, θ=exp.(θ°))
# end
function LogDensityProblems.logdensity(p::CMBProblem, zvec)
    (;ds, Ω₀) = p
    Ω = copy(first(promote(zvec, Ω₀)))
    logpdf(Mixed(ds); Ω...)
end
LogDensityProblems.dimension(p::CMBProblem) = length(p.Ω₀)
LogDensityProblems.capabilities(::Type{<:CMBProblem}) = LogDensityProblems.LogDensityOrder{0}()
##


U = Ω -> logpdf(Mixed(ds); Ω...)
Ω = Ω₀
chain = []
@showprogress for i=1:100
    Ω, = state = hmc_step(U, Ω, Λmass; symp_kwargs=[(N=25, ϵ=0.005)], progress=false, always_accept=false)
    push!(chain, state)
end

figure() do 
    plot(getindex.(chain, 2))
end
plot([Ω₀.f°, Ω.f°])
plot([Ω₀.ϕ°, Ω.ϕ°, ϕ°])



##

p = CMBProblem(ds, Ω₀)

n_samples, n_adapts = 2_000, 1_000

# Define a Hamiltonian system
metric = DiagEuclideanMetric(LogDensityProblems.dimension(p))
hamiltonian = Hamiltonian(metric, p, CMBLensing.Zygote)

find_good_stepsize(hamiltonian, Ωvec₀)

integrator = Leapfrog(initial_ϵ)


blackbox_logpdf(Ωvec₀)

gradient(blackbox_logpdf, Ωvec₀)[1]




gradient(z -> logpdf(ds; z...), z)[1]

f = FlatMap(rand(10,10))
ft = FieldTuple(;f, θ)

ft2 = ft .+ ft # broadcasting works
ft3 = Diagonal(ft2) * ft # diagonal operators
ft3'ft3 # inner prod
(;f, θ) = ft2 # unpacking components

ft2.θ

promote(ft2, ft2[:])[2].θ