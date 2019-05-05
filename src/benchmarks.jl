using BenchmarkTools
using TaylorIntegration
using DiffEqCallbacks
using ParameterizedFunctions
using Parameters

function prob_setup(g, t; rescaling=false)
    p0, q0 = ics_separate(g)
    z0 = ics(g)

    p = PhysicalParameters(B=0.15)
    if !rescaling
        prob1 = DynamicalODEProblem(ṗ, q̇, p0[1], q0[1], t, p)
        prob2 = ODEProblem(ż, z0[1], t, p)
    else
        d0 = 1e-9
        prob1 = λproblem(ż, (z0[1], z0[1].+d0/√4), t, p)
        prob2 = λproblem(ṗ, q̇, (p0[1], p0[1].+d0/√2), (q0[1], q0[1].+d0/√2), t, p)
    end
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

function timings(g, t; rescaling=false)
    p1, p2, p3 = prob_setup(g, t, rescaling=rescaling)
    suite = BenchmarkGroup()

    suite["ODEProblem"] = BenchmarkGroup()
    suite["ODEProblem"]["Vern9"] = @benchmarkable solve($p2, Vern9(), abstol=1e-14, reltol=1e-14)
    !rescaling && suite["ODEProblem"]["Vern9+ManifoldProjection"] = @benchmarkable solve($p3, Vern9(), callback=cb, abstol=1e-14, reltol=1e-14)
    !rescaling && suite["ODEProblem"]["TaylorMethod"] = @benchmarkable solve($p3, TaylorMethod(50), abstol=1e-20)

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

function E_conservation(g, t; rescaling=false, short=true)
    p = prob_setup(g, t, rescaling=rescaling)

    suite = Dict{String,NamedTuple}()
    GC.gc()
    sol, t = @timed solve(p[1], Vern9(), abstol=1e-14, reltol=1e-14)
    suite["Vern9"] = (t=t, err=energy_err(sol))
    if !rescaling && short
        GC.gc()
        sol, t = @timed solve(p[3], Vern9(), callback=cb, abstol=1e-14, reltol=1e-14)
        suite["Vern9+ManifoldProjection"] = (t=t, err=energy_err(sol))
        GC.gc()
        sol, t = @timed solve(p[3], TaylorMethod(50), abstol=1e-20)
        suite["TaylorMethod"] = (t=t, err=energy_err(sol))
    end
    GC.gc()
    sol, t = @timed solve(p[2], DPRKN12(), abstol=1e-14, reltol=1e-14)
    suite["DPRKN12"] = (t=t, err=energy_err(sol))
    GC.gc()
    sol, t = @timed solve(p[2], KahanLi8(), dt=1e-2)
    suite["KahanLi8"] = (t=t, err=energy_err(sol))
    GC.gc()
    sol, t = @timed solve(p[2], SofSpa10(), dt=1e-2)
    suite["SofSpa10"] = (t=t, err=energy_err(sol))

    return suite
end

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
