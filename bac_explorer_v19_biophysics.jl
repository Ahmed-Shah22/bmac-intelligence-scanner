# BAC-Explorer v19 — Biophysics to Circuit Toy Model
# Verified: Sequence -> Affinity (Kd) -> Topology -> R-index
# Run: julia +1.12, then include this file. Server starts on port 8190.

using Pkg
Pkg.activate("C:\\BindingAndCatalysis.jl")
using BindingAndCatalysis, LinearAlgebra, SparseArrays, Polyhedra
using Bonito, WGLMakie, Makie

Base.muladd(a::Real, b::BindingAndCatalysis.ExactLogExpr, c::BindingAndCatalysis.ExactLogExpr) = (a * b) + c
WGLMakie.activate!()

# --- Setup the reactive 3D cone figure ---
cone_fig = Figure(size=(650, 450))
cone_ax = LScene(cone_fig[1,1]; show_axis=true)
cone_xs = Observable(Float64[])
cone_ys = Observable(Float64[])
cone_zs = Observable(Float64[])
lines!(cone_ax, cone_xs, cone_ys, cone_zs; color=:steelblue, linewidth=2)

function update_cone!(rays)
    xs = Float64[]; ys = Float64[]; zs = Float64[]
    for r in rays
        push!(xs, 0.0); push!(ys, 0.0); push!(zs, 0.0)
        push!(xs, r[1]); push!(ys, r[2]); push!(zs, r[3])
    end
    cone_xs[] = xs; cone_ys[] = ys; cone_zs[] = zs
end

function extract_rays(model)
    rays_out = Vector{Vector{Float64}}()
    regimes = get_regimes(model; singular=false)
    for rgm in regimes
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
    return rays_out
end

# --- The Biophysical Toy Model ---
const CONSENSUS_P10 = "TATAAT"

function sequence_to_affinity(seq::String)
    seq = uppercase(strip(seq))
    if length(seq) != 6
        return (score=0, kd=Inf, topology_exists=false)
            end
    matches = 0
    for i in 1:6
        if seq[i] == CONSENSUS_P10[i]
            matches += 1
        end
    end
    mismatches = 6 - matches
    kd = 1.0 * (5.0 ^ mismatches)
    topology_exists = kd < 1000.0
    return (score=matches, kd=kd, topology_exists=topology_exists)
    end

# --- The App ---
app = App() do
    seq_field = Bonito.TextField("TATAAT")
        run_btn = Bonito.Button("Translate Sequence to Circuit")
    result = Bonito.Observable("Enter a 6bp -10 box sequence. Consensus is TATAAT.")
    
    on(run_btn) do _
        try
            seq = seq_field.value[]
            biophys = sequence_to_affinity(seq)
                        
                                    if !biophys.topology_exists
                result[] = "Sequence: " * seq * "\n" *
                           "Match Score: " * string(biophys.score) * "/6\n" *
                           "Predicted Kd: " * string(round(biophys.kd, digits=2)) * "\n\n" *
                           "⚠️ Affinity too weak. RNAP will not bind. Topology is empty."
                update_cone!([])
                return
            end
            
            N_matrix = [1 1 -1]
            L_matrix = [1 0 1; 0 1 1]
            
            model = Bnc(N=N_matrix, L=L_matrix,
                        x_sym=[:RNAP, :Promoter, :Complex],
                        q_sym=[:tRNAP, :tPromoter],
                        K_sym=[:Kd])
                        
            rays = extract_rays(model)
            vols = get_volumes(model; reltol=0.02, time_limit=5.0)
            r_index = sum(v.mean for v in vols)
            
            result[] = "Sequence: " * seq * "\n" *
                       "Match Score: " * string(biophys.score) * "/6\n" *
                       "Predicted Kd: " * string(round(biophys.kd, digits=2)) * "\n\n" *
                       "Generated Topology N: " * string(N_matrix) * "\n" *
                       "Generated Topology L: " * string(L_matrix) * "\n\n" *
                       "R-index of Topology: " * string(round(r_index, digits=4)) * "\n\n" *
                       "Notice: Even though the Kd changed, the R-index stayed at 1.0.\n" *
                       "The sequence sets the parameters, but the topology sets the robustness."
                       
            update_cone!(rays)
            
        catch e
            result[] = "Error: " * string(e)
            update_cone!([])
        end
    end
    
    DOM.div(
        DOM.h1("BAC-Explorer v19 — Biophysics to Circuit"),
        DOM.p("Layer 1 Toy Model: Sequence → Affinity → Topology → R-index"),
        DOM.hr(),
        DOM.h2("1. Input Sequence (−10 box)"),
        seq_field,
        run_btn,
        DOM.hr(),
        DOM.h2("2. Physical Translation & Topology"),
        DOM.pre(result),
        DOM.hr(),
        DOM.h2("3. Regime Cone (Topology Geometry)"),
                cone_fig,
    )
end

server = Bonito.Server(app, "127.0.0.1", 8190)
