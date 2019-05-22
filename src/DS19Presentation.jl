module DS19Presentation

export load, slide4, slide5, slide6, slide8_9, slide10_11, slide13, slide15,
    slide17_18_19_20, slide22_23, slide24_25, slide27, slide28, slide30_31,
    slide_a2, slide_a3, slide_a4_5, slide_a6_7, slide_a8, slide_a9,
    pgfplots, save_animation, animate

# Fix for https://github.com/JuliaIO/ImageMagick.jl/issues/140
using ImageMagick

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
#     registration = generate(Figshare(), "https://figshare.com/articles/DS19_presentation_materials/8146106")
#     print(fh, registration)
# end

const bg = colorant"#FAFAFA"
const ac = colorant"#EB811B"

function load()
    deserialize(datadep"DS19 presentation materials/graph.jls")
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

"""
    slide5(g; saveimage=false, savevideo=false)

This function returns the time node and the scene for the
simulation of the nuclear surface and its sections.
Call as `t, sc = slide5(g);` and then display the scene
using `display(sc)` and use `animate(t, (0, 40))` to animate
it.
"""
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

function nucleus_poincare(g, B, E; i=1)
    sol = get_sol(g, i, B=B, E=E)
    sim = get_sim(g, B=B, E=E)
    t = Node(0.)

    surface_sc = animate_solution(sol, t)
    line_sc = path_animation3D(sol, t)
    plot_slice!(line_sc, sim[i])

    sc = vbox(surface_sc, line_sc, sizes=[0.5,0.5])
    return t, sc
end

"""
    slide6(g; saveimage=false, savevideo=false)

This function returns the time node and the scene for the
simulation of the nuclear surface and the motion in phase space.
Call as `t, sc = slide6(g);` and then display the scene
using `display(sc)` and use `animate(t, (0, 40))` to animate
it.
"""
function slide6(g; saveimage=false, savevideo=false)
    t, sc = nucleus_poincare(g, 0.5, 0.3, i=3)
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

"""
    slide8_9(g)

This function returns the two plots for short time benchmarks.
Call as `p1, p2 = slide8_9(g)`.
"""
function slide8_9(g)
    p1, p2 = short_benchmark(g)
    savefig(p1, "assets/short-benchmark-E.tex")
    savefig(p2, "assets/short-benchmark-t.tex")

    return p1, p2
end

"""
    slide10_11(g)

This function returns the two plots for long time benchmarks.
Run twice in a row for meaningfull results, as the first time
will include precompilation times. See the julia documentation for
details. Call as `p1, p2 = slide10_11(g)`.
"""
function slide10_11(g)
    p1, p2 = long_benchmark(g)
    savefig(p1, "assets/long-benchmark-E.tex")
    savefig(p2, "assets/long-benchmark-t.tex")

    return p1, p2
end

"""
    slide13(g; saveimage=false, savevideo=false)

This function returns the time node and the scene for the
simulation of the distance between two close trajectories.
Call as `t, sc = slide13(g);` and then display the scene
using `display(sc)` and use `animate(t, (0, 40))` to animate
it.
"""
function slide13(g; saveimage=false, savevideo=false)
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

"""
    slide15(g; saveimage=false)

This function returns the interactive Poincaré section.
Call as `slide15(g)`.
"""
function slide15(g; saveimage=false)
    sc = poincare_explorer(g, 120., DynSys(T=1e5), PoincareRand(n=500), t=1e4,
        params=PhysicalParameters(B=0.5), markersize=0.08)
    if saveimage
        path = joinpath("assets", "poincare.png")
        save(path, sc)
    end
    return sc
end

"""
    slide17_18_19_20(g)

This function produces the .tex files for the histogram evolution
of λ with respect to the integration time in the `assets` folder.
"""
function slide17_18_19_20(g)
    E = 120.
    p = PhysicalParameters(B=0.5)
    ic_alg = PoincareRand(n=500)
    ic_dep = InitialConditions.depchain(p,E,ic_alg)
    for T in 10. .^(4:7)
        λs = g[:λ, ic_dep..., (λ_alg=DynSys(T=T),)][1]
        tpow = Int(log10(T))
        plt = histogram(λs, nbins = 50, xlabel=L"\lambda", ylabel=L"N",
            label=latexstring("T=10^$tpow"), lw=0, framestyle=:grid,
            color=colorant"#6699CC", background_color=bg, markerstrokealpha=0,
            tex_output_standalone=true)
        savefig(plt, "assets/hist-lambda-$tpow.tex")
    end
end

"""
    slide22_23(g)

This functuion returns the histogram with the selected chaotic part.
"""
function slide22_23(g)
    E = 120.
    p = PhysicalParameters(B=0.5)
    ic_alg = PoincareRand(n=500)
    λhist, shist = selected_hist(g, E, DynSys(T=1e5), ic_alg, params=p)
    plt = plot(λhist, xlabel=L"\lambda", ylabel=L"N", label=L"T=10^5",
        tex_output_standalone=true, lw=0, framestyle=:grid, background_color=bg,
        color=colorant"#6699CC")
    savefig(plt, "assets/hist-lambda-selected-before.tex")

    plot!(plt, shist, label="chaotic", lw=0, color=colorant"#DDCC77")
    savefig(plt, "assets/hist-lambda-selected-after.tex")
    return plt
end

"""
    slide24_25(g)

This function returns the plot of the maximal Lyapunov exponent averaged
over initial conditions.
Call as `p1, p2 = slide24_25(g)`.
"""
function slide24_25(g)
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

"""
    slide27(g)

This function returns the plot of d∞ averaged over initial conditions.
Call as `slide27(g)`.
"""
function slide27(g)
    p = PhysicalParameters(B=0.5)
    ic_alg = PoincareRand(n=500)
    plt = mean_over_ic(g, DInftyAlgorithm(T=1e5), ic_alg, params=p,
        Einterval=0.01:0.01:10, framestyle=:grid, background_color=bg,
        color=colorant"#6699CC", lw=2.5, markerstrokewidth=0, markerstrokealpha=0,
        tex_output_standalone=true)
    savefig(plt, "assets/mean-over-ic-dinf.tex")
    return plt
end

"""
    slide28(g)

This function returns the plot of Γ averaged over initial conditions.
Call as `slide28(g)`.
"""
function slide28(g)
    p = PhysicalParameters(B=0.5)
    ic_alg = PoincareRand(n=500)
    plt = mean_over_ic(g, DynSys(T=1e5), DInftyAlgorithm(T=1e5), ic_alg, params=p,
        Einterval=0.01:0.02:10, framestyle=:grid, background_color=bg,
        color=colorant"#6699CC", lw=2.5, markerstrokewidth=0, markerstrokealpha=0,
        tex_output_standalone=true)
    savefig(plt, "assets/mean-over-ic-Gamma.tex")
    return plt
end

"""
    slide30_31(g)

This function returns the plot of the maximal Lyapunov exponent averaged
over the energy.
Call as `p1, p2 = slide30_31(g)`.
"""
function slide30_31(g)
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

function slide_a2()
    x = y = range(-5, 5, length=100)
    plt = contour(x,y,
        (x,y)->V((x,y),PhysicalParameters(B=0.5)),
        xlabel=L"q_0", ylabel=L"q_2", colorbar_title=L"V",
        levels=50, framestyle=:grid, tex_output_standalone=true,
        background_color=bg)
    savefig(plt, "assets/potential-contour.tex")

    return plt
end

function slide_a3(g; saveimage=false, savevideo=false)
    t, sc = nucleus_poincare(g, 0.4, 0.1, i=1)
    if saveimage
        path = joinpath("assets", "nucleus-with-poincare-regular.png")
        save(path, sc)
    end
    if savevideo
        path = joinpath("assets", "nucleus-with-poincare-regular.mp4")
        save_animation(sc, t, (0, 40), path)
    end
    return t, sc
end

function slide_a4_5(g)
    p1, p2 = short_benchmark(g, rescaling=true)
    savefig(p1, "assets/short-benchmark-rescaling-E.tex")
    savefig(p2, "assets/short-benchmark-rescaling-t.tex")

    return p1, p2
end

function slide_a6_7(g)
    p1, p2 = long_benchmark(g, rescaling=true)
    savefig(p1, "assets/long-benchmark-rescaling-E.tex")
    savefig(p2, "assets/long-benchmark-rescaling-t.tex")

    return p1, p2
end

function slide_a8(g)
    p = PhysicalParameters(B=0.5)
    ic_alg = PoincareRand(n=500)
    plt = mean_over_ic(g, DynSys(T=1e5), ic_alg, params=p, Einterval=10:10:3000,
        framestyle=:grid, background_color=bg, color=colorant"#6699CC", lw=2.5,
        markerstrokewidth=0, markerstrokealpha=0, tex_output_standalone=true,
        label=latexstring("T=10^5"), legend=:topright)
    mean_over_ic(g, DynSys(T=1e6), ic_alg, params=p, Einterval=10:10:3000,
        plt=plt, color=colorant"#DDCC77", lw=2.5, markerstrokewidth=0,
        markerstrokealpha=0, label=latexstring("T=10^6"), legend=:topright)
    savefig(plt, "assets/mean-over-ic-T-comparison.tex")

    return plt
end

function slide_a9(g)
    p = PhysicalParameters(B=0.5)
    ic_alg = PoincareRand(n=500)
    plt = mean_over_ic(g, DynSys(T=1e5), ic_alg, params=p, Einterval=10:10:3000,
        framestyle=:grid, background_color=bg, color=colorant"#6699CC", lw=2.5,
        markerstrokewidth=0, markerstrokealpha=0, tex_output_standalone=true,
        label=latexstring("B=$(p.B)"), legend=:topright)

    p = PhysicalParameters(B=0.6)
    mean_over_ic(g, DynSys(T=1e6), ic_alg, params=p, Einterval=10:10:3000,
        plt=plt, color=colorant"#DDCC77", lw=2.5, markerstrokewidth=0,
        markerstrokealpha=0, label=latexstring("B=$(p.B)"), legend=:topright)
    savefig(plt, "assets/mean-over-ic-B-difference.tex")

    return plt
end

end # module
