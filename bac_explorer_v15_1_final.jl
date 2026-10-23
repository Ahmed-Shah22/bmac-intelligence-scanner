using Pkg
Pkg.activate("C:\\BindingAndCatalysis.jl")
using BindingAndCatalysis, LinearAlgebra, SparseArrays, Polyhedra
using Bonito, WGLMakie, Makie
Base.muladd(a::Real, b::BindingAndCatalysis.ExactLogExpr, c::BindingAndCatalysis.ExactLogExpr) = (a * b) + c
WGLMakie.activate!()

const PROMOTER_DB = Dict(
    "TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC" => ("J23100", "constitutive"),
    "TTTACAGCTAGCTCAGTCCTAGGTATTATGCTAGC" => ("J23101", "constitutive"),
    "TTGACAGCTAGCTCAGTCCTAGGTACTGTGCTAGC" => ("J23102", "constitutive"),
    "CTGATAGCTAGCTCAGTCCTAGGGATTATGCTAGC" => ("J23103", "constitutive"),
    "TTGACAGCTAGCTCAGTCCTAGGTATTGTGCTAGC" => ("J23104", "constitutive"),
    "TTTACGGCTAGCTCAGTCCTAGGTACTATGCTAGC" => ("J23105", "constitutive"),
    "TTTACGGCTAGCTCAGTCCTAGGTATAGTGCTAGC" => ("J23106", "constitutive"),
    "TAACCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAAGGATCC" => ("lacUV5", "repressed"),
    "TCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATACTGAGCAC" => ("pTet", "repressed"),
    "CTCCATACCCGTTTTTTTGGGCTAGC" => ("pBAD", "activated"),
)

const TOPOLOGY_N = Dict(
    "constitutive" => ("1,1,-1", "1,0,1 ; 0,1,1"),
    "repressed"    => ("1,1,-1,0,0 ; 0,1,0,1,-1", "1,0,1,0,0 ; 0,1,1,0,1 ; 0,0,0,1,1"),
    "activated"    => ("1,1,-1,0,0 ; 0,0,1,1,-1", "1,0,1,0,0 ; 0,1,1,1,1 ; 0,0,0,1,1"),
    )
    
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

cone_fig = Figure(size=(650, 450))
cone_ax = LScene(cone_fig[1,1]; show_axis=true)
cone_xs = Observable(Float64[])
cone_ys = Observable(Float64[])
cone_zs = Observable(Float64[])
lines!(cone_ax, cone_xs, cone_ys, cone_zs; color=:steelblue, linewidth=2)

app = App() do
    seq_field = Bonito.TextField("TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC")
        analyze_btn = Bonito.Button("Analyze sequence")
    download_btn = Bonito.Button("Download results")
    result = Bonito.Observable("Paste a promoter sequence, then click Analyze.")
    interpretation = Bonito.Observable("")
    save_status = Bonito.Observable("")
    last_result = Ref("")

    on(analyze_btn) do _
        save_status[] = ""
        try
            seq = uppercase(strip(seq_field.value[]))
            if isempty(seq)
                result[] = "Empty sequence."
                interpretation[] = ""
                cone_xs[] = Float64[]; cone_ys[] = Float64[]; cone_zs[] = Float64[]
                return
            end
            if !haskey(PROMOTER_DB, seq)
                result[] = "Unknown sequence."
                interpretation[] = "Not in the reference dataset. A learned sequence-to-topology model is the next research step."
                last_result[] = ""
                cone_xs[] = Float64[]; cone_ys[] = Float64[]; cone_zs[] = Float64[]
                return
            end
            name, topology = PROMOTER_DB[seq]
            N_str, L_str = TOPOLOGY_N[topology]
            parse_row(s) = parse.(Int, split(strip(s), r"[,\s]+"))
            n_rows = split(strip(N_str), r"[;\n]+")
            l_rows = split(strip(L_str), r"[;\n]+")
            N_matrix = reduce(vcat, [reshape(parse_row(r), 1, :) for r in n_rows if !isempty(strip(r))])
            L_matrix = reduce(vcat, [reshape(parse_row(r), 1, :) for r in l_rows if !isempty(strip(r))])
                        model, rays = extract_rays(N_matrix, L_matrix)
            volumes = get_volumes(model; reltol=0.01, time_limit=15.0)
            r_total = sum(v.mean for v in volumes)
            n_pos = count(v -> v.mean > 0, volumes)

            result[] = "Recognized: " * name *
                                   "\nTopology: " * topology *
                       "\nReactions: " * string(size(N_matrix, 1)) *
                       "\nSpecies: " * string(size(N_matrix, 2)) *
                       "\n\nR-index total: " * string(round(r_total, digits=4)) *
                       "\nPositive regimes: " * string(n_pos) *
                       "\nCone rays: " * string(length(rays))

            interp = r_total > 0.75 ? "Robust design. Most of parameter space produces the desired behavior." :
                     r_total > 0.4  ? "Moderate robustness. About half of parameter space produces the desired behavior." :
                     r_total > 0.1  ? "Fragile design. Only a small fraction of parameter space produces the desired behavior." :
                                      "Very fragile. The desired behavior occupies an extremely narrow region."
                                                  interpretation[] = interp
            last_result[] = "Sequence: " * seq * "\nName: " * name * "\nTopology: " * topology *
                            "\nR-index: " * string(round(r_total, digits=4)) * "\nRegimes: " * string(n_pos) *
                            "\nInterpretation: " * interp * "\n"

            xs = Float64[]; ys = Float64[]; zs = Float64[]
            for r in rays
                push!(xs, 0.0); push!(ys, 0.0); push!(zs, 0.0)
                push!(xs, r[1]); push!(ys, r[2]); push!(zs, r[3])
            end
            cone_xs[] = xs; cone_ys[] = ys; cone_zs[] = zs
        catch e
            result[] = "Error: " * string(e)
            interpretation[] = ""
            cone_xs[] = Float64[]; cone_ys[] = Float64[]; cone_zs[] = Float64[]
        end
    end

    on(download_btn) do _
        if isempty(last_result[])
            save_status[] = "Nothing to download yet."
                    else
            path = "C:\\bmac-julia\\bac_result_$(rand(1000:9999)).txt"
                        write(path, last_result[])
            save_status[] = "Saved to: " * path
        end
    end

    DOM.div(
        DOM.h1("BAC-Explorer v15.1"),
        DOM.hr(),
        DOM.h2("1. Input"), seq_field, analyze_btn,
        DOM.hr(),
        DOM.h2("2. Analysis"), DOM.pre(result),
        DOM.hr(),
        DOM.h2("3. Interpretation"), DOM.p(interpretation),
        DOM.hr(),
        DOM.h2("4. Export"), download_btn, DOM.p(save_status),
        DOM.hr(),
        DOM.h2("5. Regime cone"),
        cone_fig,
    )
end

server = Bonito.Server(app, "127.0.0.1", 8155)

