using BenchmarkTools
using TaylorIntegration
using DiffEqCallbacks
using ParameterizedFunctions
using Parameters

function prob_setup(g, t)
    p0, q0 = ics_separate(g)
    z0 = ics(g)

    p = PhysicalParameters(B=0.15)
    prob1 = DynamicalODEProblem(ṗ, q̇, p0[1], q0[1], t, p)
    prob2 = ODEProblem(ż, z0[1], t, p)
    @unpack A,B,D = p
    prob3 = ODEProblem(h_eqs, Array(z0[1]), t, (A,B,D))

    return prob1, prob2, prob3
end

h_eqs = @ode_def HamiltEqs begin
  dp₀ = -A * q₀ - 3 * B / √2 * (q₂^2 - q₀^2) - D * q₀ * (q₀^2 + q₂^2)
  dp₂ = -q₂ * (A + 3 * √2 * B * q₀ + D * (q₀^2 + q₂^2))
  dq₀ = A * p₀
  dq₂ = A * p₂
end A B D

function g(resid, u)
    resid[1] = H((u[1],u[2]),(u[3],u[4]), PhysicalParameters(B=0.15)) - E
    resid[2:4] .= 0
end

const cb = ManifoldProjection(g, nlopts=Dict(:ftol=>1e-13))
const E = 0.1

function run_suite(g, t)
    p1, p2, p3 = prob_setup(g, t)
    suite = BenchmarkGroup()

    suite["ODEProblem"] = BenchmarkGroup()
    suite["ODEProblem"]["Vern9"] = @benchmarkable solve($p2, Vern9(), abstol=1e-14, reltol=1e-14)
    suite["ODEProblem"]["Vern9+ManifoldProjection"] = @benchmarkable solve($p3, Vern9(), callback=cb, abstol=1e-14, reltol=1e-14)
    suite["ODEProblem"]["TaylorMethod"] = @benchmarkable solve($p3, TaylorMethod(50), abstol=1e-20)

    suite["DynamicalODEProblem"] = BenchmarkGroup()
    suite["DynamicalODEProblem"]["DPRKN12"] =  @benchmarkable solve($p1, DPRKN12(), abstol=1e-14, reltol=1e-14)
    suite["DynamicalODEProblem"]["KahanLi8"] =  @benchmarkable solve($p1, KahanLi8(), dt=1e-2)
    suite["DynamicalODEProblem"]["SofSpa10"] =  @benchmarkable solve($p1, SofSpa10(), dt=1e-2)

    tune!(suite)
    run(suite)
end

energy_err(sol) = map(
    i->H((sol[1,i], sol[2,i]), (sol[3,i], sol[4,i]), PhysicalParameters(B=0.15)) - E,
    axes(sol,2))

function E_conservation(g, t; extended=false)
    p = prob_setup(g, t, extended=extended)

    suite = Dict{String,NamedTuple}()
end
# 
# h=load()
#
# results = run_suite(h, 10.)
#
# results["ODEProblem"]["Vern9"].times
#
# p1, p2 = prob_setup(g, 10.)
#
# sol = solve(p1, DPRKN12(), abstol=1e-14, reltol=1e-14)
#
# i=12
#
# energy_err(sol)
