module DS19Presentation

export load, slide4, slide5

using DataDeps
using MD5
using Serialization
using StorageGraphs

ENV["DATADEPS_ALWAYS_ACCEPT"] = true

function __init__()
    include("$(@__DIR__)/data.jl")
end

using NuclearSurfaceVibrations
using .Classical
using .Visualizations
using .InitialConditions: depchain, initial_conditions!
using StaticArrays
using OrdinaryDiffEq
using Makie

# open("src/data.jl", "w") do fh
#     registration = generate(Figshare(), "https://figshare.com/articles/DS19_test_database/8078699")
#     print(fh, registration)
# end

function load()
    deserialize(datadep"DS19 test database/graph.jls")
end

function get_sol(g, E, ic_alg, p)
    q0, p0 = initial_conditions!(g, E, alg=ic_alg, params=p)
    p₀ = [SVector{2}(p0[i, :]) for i ∈ axes(p0, 1)]
    q₀ = [SVector{2}(q0[i, :]) for i ∈ axes(q0, 1)]
    z0 = [vcat(p₀[i], q₀[i]) for i ∈ axes(q₀, 1)]
    prob = ODEProblem(ż, z0[1], 100., p)
    sol = solve(prob, Vern9(), abstol=1e-14, reltol=1e-14, maxiters=1e9)
end

function slide4(g; saveimage=false)
    E = 0.1
    p = PhysicalParameters(B=0.15)
    ic_alg = PoincareRand(n=10)

    sol = get_sol(g, E, ic_alg, p)
    t = Node(0.)
    surface_sc = animate_solution(sol, t)
    if saveimage
        save("assets/nucleus.png", surface_sc)
    end
    return surface_sc
end

function slide5(g; saveimage=false, savevideo=false)
    E = 0.1
    p = PhysicalParameters(B=0.15)
    ic_alg = PoincareRand(n=10)

    sol = get_sol(g, E, ic_alg, p)
    t = Node(0.)
    surface_sc = animate_solution(sol, t)
    section_sc = θϕ_sections(sol, t, surface_sc.limits[])
    sc = vbox(surface_sc, section_sc, sizes=[0.7, 0.3])
    if saveimage
        save("assets/nucleus-with-sections.png", sc)
    end
    if savevideo
        save_animation(sc, t, (0, 40), "nucleus-with-sections.webm")
    end
    return t, surface_sc
end

end # module
