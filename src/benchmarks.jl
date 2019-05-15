using BenchmarkTools
using TaylorIntegration
using DiffEqCallbacks
using ParameterizedFunctions
using RecursiveArrayTools
using Parameters
using Plots: plot, plot!, scatter, scatter!, histogram, histogram!, savefig,
    pgfplots
using Colors

pgfplots()

function prob_setup(g, t; rescaling=false, Ttr=1e3, τ=5.)
    p0, q0 = ics_separate(g, E, PoincareRand(n=5), 0.15)
    z0 = ics(g)

    p = PhysicalParameters(B=0.15)
    if !rescaling
        prob1 = DynamicalODEProblem(ṗ, q̇, p0[1], q0[1], t, p)
        prob2 = ODEProblem(ż, z0[1], t, p)
    else
        alg = TimeRescaling(T=t, Ttr=Ttr, τ=τ)
        prob1 = λproblem(p0[1], q0[1], alg, params=p)
        prob2 = λproblem(z0[1], alg, params=p)
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
    if !rescaling
        suite["ODEProblem"]["Vern9+ManifoldProjection"] = @benchmarkable solve($p3, Vern9(), callback=cb, abstol=1e-14, reltol=1e-14)
        suite["ODEProblem"]["TaylorMethod"] = @benchmarkable solve($p3, TaylorMethod(50), abstol=1e-20)
    end

    suite["DynamicalODEProblem"] = BenchmarkGroup()
    suite["DynamicalODEProblem"]["DPRKN12"] =  @benchmarkable solve($p1, DPRKN12(), abstol=1e-14, reltol=1e-14)
    suite["DynamicalODEProblem"]["KahanLi8"] =  @benchmarkable solve($p1, KahanLi8(), dt=1e-2)
    suite["DynamicalODEProblem"]["SofSpa10"] =  @benchmarkable solve($p1, SofSpa10(), dt=1e-2)

    tune!(suite)
    run(suite)
end

energy_err(sol,offset=0) = DiffEqArray(map(
    i->H((sol[1,i], sol[2,i]),
         (sol[3+offset,i], sol[4+offset,i]),
         PhysicalParameters(B=0.15)
        ) - E,
    axes(sol,2)), sol.t)

function E_conservation(g, t; rescaling=false, short=true)
    if short
        Ttr = 0.
        saveat = Float64[]
    else
        Ttr = 1e3
        saveat = 1.
    end
    p = prob_setup(g, t, rescaling=rescaling, Ttr=Ttr)

    suite = Dict{String,NamedTuple}()
    GC.gc()
    sol, t = @timed solve(p[2], Vern9(), abstol=1e-14, reltol=1e-14, saveat=saveat)
    suite["Vern9"] = (t=t, err=energy_err(sol))
    if !rescaling
        GC.gc()
        sol, t = @timed solve(p[3], Vern9(), callback=cb, abstol=1e-14, reltol=1e-14)
        suite["Vern9+ManifoldProjection"] = (t=t, err=energy_err(sol))
        if short
            GC.gc()
            sol, t = @timed solve(p[3], TaylorMethod(50), abstol=1e-20)
            suite["TaylorMethod"] = (t=t, err=energy_err(sol))
        end
    end
    offset = rescaling ? 2 : 0
    GC.gc()
    sol, t = @timed solve(p[1], DPRKN12(), abstol=1e-14, reltol=1e-14)
    suite["DPRKN12"] = (t=t, err=energy_err(sol, offset))
    GC.gc()
    sol, t = @timed solve(p[1], KahanLi8(), dt=1e-2, saveat=saveat)
    suite["KahanLi8"] = (t=t, err=energy_err(sol, offset))
    GC.gc()
    sol, t = @timed solve(p[1], SofSpa10(), dt=1e-2, saveat=saveat)
    suite["SofSpa10"] = (t=t, err=energy_err(sol, offset))

    return suite
end

function collect_results(g, t; rescaling=false, short=true)
    E_results = E_conservation(g, t, rescaling=rescaling, short=short)
    if short
        t_results = timings(g, t, rescaling=rescaling)
        integ1 = keys(t_results["ODEProblem"])
        integ2 = keys(t_results["DynamicalODEProblem"])
        solvers = Iterators.flatten((integ1,integ2))

        # times are in nanoseconds
        ts1 = [minimum(t_results["ODEProblem"][i].times) for i in integ1]
        ts2 = [minimum(t_results["DynamicalODEProblem"][i].times) for i in integ2]
        ts = Iterators.flatten((ts1,ts2))
    else
        solvers = keys(E_results)
        ts = [E_results[i].t for i in solvers]
    end
    E_err = [E_results[i].err for i in solvers]

    return solvers, ts, E_err
end

function short_benchmark(g; rescaling=false)
    solvers, ts, E_errs = collect_results(g, 100., rescaling=rescaling, short=true)
    p1 = plot(background_color=bg,
              xlabel="t",
              ylabel="Energy error",
              framestyle=:grid,
              tex_output_standalone=true
        )
    for (E_err,solver) in zip(E_errs, solvers)
        plot!(p1, E_err, label=solver)
    end
    p2 = plot(background_color=bg,
              xlabel="Solver",
              ylabel="Computation time (ms)",
              framestyle=:grid,
              tex_output_standalone=true,
              legend=false
        )
    scatter!(p2, string.(solvers), ts./1e6, m=10, color=ac, markerstrokealpha=0,
        markerstrokewidth=0)

    return p1, p2
end

function long_benchmark(g; rescaling=false)
    solvers, ts, E_errs = collect_results(g, 1e4, rescaling=rescaling, short=false)
    p1 = plot(background_color=bg,
              xlabel="t",
              ylabel="Energy error",
              framestyle=:grid,
              tex_output_standalone=true
        )
    for (E_err,solver) in zip(E_errs, solvers)
        plot!(p1, E_err, label=solver)
    end
    p2 = plot(background_color=bg,
              xlabel="Solver",
              ylabel="Computation time (s)",
              framestyle=:grid,
              tex_output_standalone=true,
              legend=false
        )
    scatter!(p2, string.(solvers), ts, m=10, color=ac, markerstrokewidth=0,
        markerstrokealpha=0)

    return p1, p2
end
