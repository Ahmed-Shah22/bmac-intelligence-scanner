# BAC-Explorer v12 — reactive cone, working
# Tests: single reaction (6 rays), independent pair (36 rays), repressed (32 rays)
#
# To run: julia +1.12, then include this file. Server starts on port 8135.

# ==== Load environment ====
using Pkg
Pkg.activate("C:\\\\BindingAndCatalysis.jl")
using BindingAndCatalysis, LinearAlgebra, SparseArrays, Polyhedra
using Bonito, WGLMakie, Makie
Base.muladd(a::Real, b::BindingAndCatalysis.ExactLogExpr, c::BindingAndCatalysis.ExactLogExpr) = (a * b) + c
WGLMakie.activate!()

# ==== extract_rays, figure, and app code follow ====
# Paste the working extract_rays function, figure setup, and App() do block here.
function extract_rays(N_matrix, L_matrix)
    n_species = size(N_matrix, 2)
    x_sym = [Symbol("x$i") for i in 1:n_species]
    q_sym = [Symbol("t$i") for i in 1:size(L_matrix, 1)]
    K_sym = [Symbol("K$i") for i in 1:size(N_matrix, 1)]
    model = Bnc(N=N_matrix, L=L_matrix, x_sym=x_sym, q_sym=q_sym, K_sym=K_sym)
        rays_out = Vector{Vector{Float64}}()
            for rgm in get_regimes(model; singular=false)
                    try
            poly = get_polyhedron(rgm; chart=:qK)
            for ray in collect(Polyhedra.rays(Polyhedra.vrep(poly)))
                            v = ray.a
                if any(x -> abs(x) > 1e-9, v) && length(v) >= 3
                    push!(rays_out, v[1:3])
                end
            end
        catch
            continue
        end
    end
    return model, rays_out
end

fig = Figure(size=(700, 500))
ax = LScene(fig[1,1]; show_axis=true)
xs_obs = Observable(Float64[])
ys_obs = Observable(Float64[])
zs_obs = Observable(Float64[])
lines!(ax, xs_obs, ys_obs, zs_obs; color=:steelblue, linewidth=2)

app = App() do
    n_field = Bonito.TextField("1,1,-1")
        l_field = Bonito.TextField("1,0,1 ; 0,1,1")
    btn = Bonito.Button("Compute R-index + Cone")
    result = Bonito.Observable("Click to compute. Cone updates with each circuit.")
        on(btn) do _
        try
            parse_row(s) = parse.(Int, split(strip(s), r"[,\s]+"))
            n_rows = split(strip(n_field.value[]), r"[;\n]+")
            l_rows = split(strip(l_field.value[]), r"[;\n]+")
            N_matrix = reduce(vcat, [reshape(parse_row(r), 1, :) for r in n_rows if !isempty(strip(r))])
            L_matrix = reduce(vcat, [reshape(parse_row(r), 1, :) for r in l_rows if !isempty(strip(r))])
                        model, rays = extract_rays(N_matrix, L_matrix)
            volumes = get_volumes(model; reltol=0.01, time_limit=15.0)
            r_total = sum(v.mean for v in volumes)
            n_pos = count(v -> v.mean > 0, volumes)
            result[] = "Total R-index: " * string(round(r_total, digits=4)) *
                       "\nPositive regimes: " * string(n_pos) *
                       "\nCone rays: " * string(length(rays))
                                   xs = Float64[]; ys = Float64[]; zs = Float64[]
            for r in rays
                push!(xs, 0.0); push!(ys, 0.0); push!(zs, 0.0)
                push!(xs, r[1]); push!(ys, r[2]); push!(zs, r[3])
            end
            xs_obs[] = xs; ys_obs[] = ys; zs_obs[] = zs
        catch e
            result[] = "Error: " * string(e)
        end
    end
    DOM.div(
            DOM.h1("BAC-Explorer v12"),
        DOM.p("N matrix:"), n_field,
        DOM.p("L matrix:"), l_field,
        btn,
        DOM.hr(),
        DOM.h2("Result"),
                DOM.pre(result),
                        DOM.hr(),
                                DOM.h2("Regime Cone"),
        fig,
    )
end

server = Bonito.Server(app, "127.0.0.1", 8135)

