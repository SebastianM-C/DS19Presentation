module DS19Presentation

export load, slide4, slide5, slide6, slide7, slide8, slide9_10, slide11_12,
    slide13_14_15_16, slide17_18, slide19_20, slide21_22, slide23_24, slide25_26,
    slide_a1_2, slide_a3_4,
    pgfplots, save_animation, animate

using DataDeps
using MD5
using Serialization
using StorageGraphs

# ENV["DATADEPS_ALWAYS_ACCEPT"] = true

function __init__()
    include("$(@__DIR__)/data.jl")
end

using NuclearSurfaceVibrations
using .Classical
using .Visualizations
using .InitialConditions: depchain, initial_conditions!
using .Classical.ParallelTrajectories
using StaticArrays
using OrdinaryDiffEq
using LaTeXStrings
using IntervalArithmetic
using StatsBase
using AbstractPlotting: Node, save, hbox, vbox

include("benchmarks.jl")

# open("src/data.jl", "w") do fh
#     registration = generate(Figshare(), "https://figshare.com/articles/DS19_test_database/8078699")
#     print(fh, registration)
# end

const bg = colorant"#FAFAFA"
const ac = colorant"#EB811B"

function load()
    deserialize(datadep"DS19 presentation materials test/graph.jls")
end

function ics_separate(g, E=0.1, ic_alg=PoincareRand(n=5), B=0.15)
    p = PhysicalParameters(B=B)
    q0, p0 = initial_conditions!(g, E, alg=ic_alg, params=p)
    p₀ = [SVector{2}(p0[i, :]) for i ∈ axes(p0, 1)]
    q₀ = [SVector{2}(q0[i, :]) for i ∈ axes(q0, 1)]

    return p₀, q₀
end

function ics(g, E=0.1, ic_alg=PoincareRand(n=5), B=0.15)
    p₀, q₀ = ics_separate(g, E, ic_alg, B)
    z0 = [vcat(p₀[i], q₀[i]) for i ∈ axes(q₀, 1)]
end

function get_sol(g, i=1, ic_alg=PoincareRand(n=5); B=0.15, E=0.1)
    p = PhysicalParameters(B=B)
    z0 = ics(g, E, ic_alg, B)
    prob = ODEProblem(ż, z0[i], 40., p)
    sol = solve(prob, Vern9(), abstol=1e-14, reltol=1e-14, maxiters=1e9)
end

function get_sim(g, ic_alg=PoincareRand(n=5); B=0.15, E=0.1)
    p = PhysicalParameters(B=B)
    q0, p0 = initial_conditions!(g, E, alg=ic_alg, params=p)
    poincaremap(q0, p0, params=p, t=40)
end

function get_psol(g, i=1, ic_alg=PoincareRand(n=5); B=0.15, E=0.1)
    p = PhysicalParameters(B=B)
    z0 = ics(g, E, ic_alg, B)
    d0 = 1e-3
    pprob = parallel_problem(ż, (z0[i], z0[i].+d0/√4), 40., p)
    psol = solve(pprob, Vern9(), abstol=1e-14, reltol=1e-14, maxiters=1e9)
end

function slide4(g; saveimage=false)
    sol = get_sol(g)
    t = Node(0.)

    surface_sc = animate_solution(sol, t)
    if saveimage
        path = joinpath("assets", "nucleus.png")
        save(path, surface_sc)
    end
    return surface_sc
end

function slide5(g; saveimage=false, savevideo=false)
    sol = get_sol(g)
    t = Node(0.)

    surface_sc = animate_solution(sol, t)
    section_sc = θϕ_sections(sol, t, surface_sc.limits[])
    sc = vbox(surface_sc, section_sc, sizes=[0.7, 0.3])
    if saveimage
        path = joinpath("assets", "nucleus-with-sections.png")
        save(path, sc)
    end
    if savevideo
        path = joinpath("assets", "nucleus-with-sections.mp4")
        save_animation(sc, t, (0, 10), path)
    end
    return t, sc
end

function slide6(g; saveimage=false, savevideo=false)
    i = 3
    sol = get_sol(g, i, B=0.5, E=0.3)
    sim = get_sim(g, B=0.5, E=0.3)
    t = Node(0.)

    surface_sc = animate_solution(sol, t)
    line_sc = path_animation3D(sol, t)
    plot_slice!(line_sc, sim[i])

    sc = vbox(surface_sc, line_sc, sizes=[0.5,0.5])
    if saveimage
        path = joinpath("assets", "nucleus-with-poincare.png")
        save(path, sc)
    end
    if savevideo
        path = joinpath("assets", "nucleus-with-poincare.mp4")
        save_animation(sc, t, (0, 40), path)
    end
    return t, sc
end

function slide7(g; saveimage=false, savevideo=false)
    psol = get_psol(g, 3, B=0.5, E=0.3)
    t = Node(0.)

    parallel_sc = parallel_paths(psol, t, linecolor=ac)
    psc = paths_distance(psol, t)
    lpsc = paths_distance_log(psol, t)
    dsc = hbox(psc, lpsc)
    sc = vbox(parallel_sc, dsc, sizes=[0.5, 0.5])

    if saveimage
        path = joinpath("assets", "dist-with-log.png")
        save(path, sc)
    end
    if savevideo
        path = joinpath("assets", "dist-with-log.mp4")
        save_animation(sc, t, (0, 40), path)
    end
    return t, sc
end

function slide8(g; saveimage=false)
    sc = poincare_explorer(g, 120., DynSys(T=1e5), PoincareRand(n=500), t=1e4,
        params=PhysicalParameters(B=0.5), markersize=0.08)
    if saveimage
        path = joinpath("assets", "poincare.png")
        save(path, sc)
    end
    return sc
end

function slide9_10(g)
    p1, p2 = short_benchmark(g)
    savefig(p1, "assets/short-benchmark-E.tex")
    savefig(p2, "assets/short-benchmark-t.tex")

    return p1, p2
end

function slide11_12(g)
    p1, p2 = long_benchmark(g)
    savefig(p1, "assets/long-benchmark-E.tex")
    savefig(p2, "assets/long-benchmark-t.tex")

    return p1, p2
end

function slide_a1_2(g)
    p1, p2 = short_benchmark(g, rescaling=true)
    savefig(p1, "assets/short-benchmark-rescaling-E.tex")
    savefig(p2, "assets/short-benchmark-rescaling-t.tex")

    return p1, p2
end

function slide_a3_4(g)
    p1, p2 = long_benchmark(g, rescaling=true)
    savefig(p1, "assets/long-benchmark-rescaling-E.tex")
    savefig(p2, "assets/long-benchmark-rescaling-t.tex")

    return p1, p2
end

function slide13_14_15_16(g)
    E = 120.
    p = PhysicalParameters(B=0.5)
    ic_alg = PoincareRand(n=500)
    ic_dep = InitialConditions.depchain(p,E,ic_alg)
    for T in 10. .^(4:7)
        λs = g[:λ, ic_dep..., (λ_alg=DynSys(T=T),)][1]
        tpow = Int(log10(T))
        plt = histogram(λs, nbins = 50,
            label=latexstring("T=10^$tpow"), lw=0, framestyle=:grid,
            color=colorant"#6699CC", background_color=bg, markerstrokealpha=0,
            tex_output_standalone=true)
        savefig(plt, "assets/hist-lambda-$tpow.tex")
    end
end

function slide17_18(g)
    E = 120.
    p = PhysicalParameters(B=0.5)
    ic_alg = PoincareRand(n=500)
    λhist, shist = selected_hist(g, E, DynSys(T=1e5), ic_alg, params=p)
    plt = plot(λhist, label=L"T=10^5", tex_output_standalone=true,
        lw=0, framestyle=:grid, background_color=bg, color=colorant"#6699CC")
    savefig(plt, "assets/hist-lambda-selected-before.tex")

    plot!(plt, shist, label="chaotic", lw=0, color=colorant"#DDCC77")
    savefig(plt, "assets/hist-lambda-selected-after.tex")
    return plt
end

function slide19_20(g)
    p = PhysicalParameters(B=0.5)
    ic_alg = PoincareRand(n=500)
    p1 = mean_over_ic(g, DynSys(T=1e5), ic_alg, params=p, Einterval=10:10:3000,
        framestyle=:grid, background_color=bg, color=colorant"#6699CC", lw=2.5,
        markerstrokewidth=0, markerstrokealpha=0, tex_output_standalone=true)
    savefig(p1, "assets/mean-over-ic.tex")
    p2 = mean_over_ic(g, DynSys(T=1e5), ic_alg, params=p, Einterval=0.01:0.01:10,
        framestyle=:grid, background_color=bg, color=colorant"#6699CC", lw=2.5,
        markerstrokewidth=0, markerstrokealpha=0, tex_output_standalone=true)
    savefig(p2, "assets/mean-over-ic-low.tex")
    return p1, p2
end

function slide21_22(g)
    p = PhysicalParameters(B=0.5)
    ic_alg = PoincareRand(n=500)
    p1 = mean_over_E(g, DynSys(T=1e5), 10:10:3000, ic_alg=ic_alg,
        Binterval=0.1:0.02:0.6, legend=false, markersize=9,
        framestyle=:grid, background_color=bg, color=colorant"#6699CC",
        markerstrokewidth=0, markerstrokealpha=0, tex_output_standalone=true)
    savefig(p1, "assets/mean-over-E.tex")
    p2 = mean_over_E(g, DynSys(T=1e5), 0.01:0.01:10, ic_alg=ic_alg,
        Binterval=0.1:0.02:0.6, legend=false, markersize=9,
        framestyle=:grid, background_color=bg, color=colorant"#6699CC",
        markerstrokewidth=0, markerstrokealpha=0, tex_output_standalone=true)
    savefig(p2, "assets/mean-over-E-low.tex")
    return p1, p2
end

function slide23_24(g)
    p = PhysicalParameters(B=0.5)
    ic_alg = PoincareRand(n=500)
    # p1 = mean_over_ic(g, DInftyAlgorithm(T=1e5), ic_alg, params=p,
    #     Einterval=10:10:3000, framestyle=:grid, background_color=bg,
    #     color=colorant"#6699CC", lw=2.5, markerstrokewidth=0, markerstrokealpha=0,
    #     tex_output_standalone=true)
    # savefig(p1, "assets/mean-over-ic-dinf.tex")
    p2 = mean_over_ic(g, DInftyAlgorithm(T=1e5), ic_alg, params=p,
        Einterval=0.01:0.01:10, framestyle=:grid, background_color=bg,
        color=colorant"#6699CC", lw=2.5, markerstrokewidth=0, markerstrokealpha=0,
        tex_output_standalone=true)
    savefig(p2, "assets/mean-over-ic-low-dinf.tex")
    return p2
end

function slide25_26(g)
    p = PhysicalParameters(B=0.5)
    ic_alg = PoincareRand(n=500)
    # p1 = mean_over_ic(g, DynSys(T=1e5), DInftyAlgorithm(T=1e5), ic_alg, params=p,
    #     Einterval=10:10:3000, framestyle=:grid, background_color=bg,
    #     color=colorant"#6699CC", lw=2.5, markerstrokewidth=0, markerstrokealpha=0,
    #     tex_output_standalone=true)
    # savefig(p1, "assets/mean-over-ic-Gamma.tex")
    p2 = mean_over_ic(g, DynSys(T=1e5), DInftyAlgorithm(T=1e5), ic_alg, params=p,
        Einterval=0.01:0.02:10, framestyle=:grid, background_color=bg,
        color=colorant"#6699CC", lw=2.5, markerstrokewidth=0, markerstrokealpha=0,
        tex_output_standalone=true)
    savefig(p2, "assets/mean-over-ic-low-Gamma.tex")
    return p2
end

end # module
