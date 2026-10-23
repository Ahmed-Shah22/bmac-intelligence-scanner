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

function build_catalysis_model()
    N = [ 1  1  0  0  0  0 -1  0  0  0 ;
          0  0  1  0  1  0  0 -1  0  0 ;
          0  0  1  1  0  0  0  0 -1  0 ;
          0  0  0  0  1  1  0  0  0 -1 ]
    x_sym = [:I, :A, :Astar, :B, :Bstar, :E, :C1, :C2, :C3, :C4]
    q_sym = [:tI, :tA, :tAstar, :tB, :tBstar, :tE]
    K_sym = [:KA1, :KA2, :KB1, :KB2]
    model = Bnc(N=N, x_sym=x_sym, q_sym=q_sym, K_sym=K_sym)
    Pi_mat = diagm(ones(Int, 4))
    Gamma_mat = [  1  -1   0   0 ;
                  -1   1   0   0 ;
                   0   0   1  -1 ;
                   0   0  -1   1 ]
    update_catalysis!(model;
        Π = Pi_mat, Γ = Gamma_mat,
        x_picked = [:C1, :C2, :C3, :C4],
        q_picked = [:tAstar, :tA, :tBstar, :tB],
        w_sym = [:tA, :tB],
        k_sym = [:kA1, :kA2, :kB1, :kB2])
    return model
end

function binding_order_matrix(rgm)
    binding_rgm = get_binding_regime(rgm)
    is_singular(binding_rgm) ? get_H_numerically(binding_rgm) : get_H(binding_rgm)
end

function binding_order_scalar(rgm, output_sym::Symbol, input_sym::Symbol)
    model = get_binding_network(rgm)
    H = binding_order_matrix(rgm)
    return H[locate_sym_x(model, output_sym), locate_sym_qK(model, input_sym)]
end

function steady_state_invariance(rgm; atol=1e-6)
    model = get_binding_network(rgm)
    tI_idx = locate_sym_wKk(model, :tI)
    active_idx = locate_sym_x(model, :Astar)
    total_active_idx = locate_sym_qcat(model, :tAstar)
    H = get_H(rgm)
    F, _ = get_qcat_F_F0(rgm)
    return (;
        A_active = abs(H[active_idx, tI_idx]) <= atol,
        tA_active = abs(F[total_active_idx, tI_idx]) <= atol)
end

function input_responsiveness(rgm; threshold=0.1)
    flux_drive = binding_order_scalar(rgm, :C1, :tI) - binding_order_scalar(rgm, :C2, :tI)
    active_gate = binding_order_scalar(rgm, :Astar, :tAstar)
    return (;
        A_active = (flux_drive > threshold) && (active_gate > threshold),
        tA_active = flux_drive > threshold)
end

function adaptation_signature(rgm; threshold=0.1)
    inv = steady_state_invariance(rgm)
    resp = input_responsiveness(rgm; threshold)
    return (;
        A_active = inv.A_active && resp.A_active,
        tA_active = inv.tA_active && resp.tA_active)
end

function parse_row(s)
    parse.(Int, split(strip(s), r"[,\s]+"))
end

function build_binding_model(N_str, L_str)
    N_matrix = reduce(vcat, [reshape(parse_row(r), 1, :) for r in split(strip(N_str), r"[;\n]+") if !isempty(strip(r))])
    L_matrix = reduce(vcat, [reshape(parse_row(r), 1, :) for r in split(strip(L_str), r"[;\n]+") if !isempty(strip(r))])
    return Bnc(N=N_matrix, L=L_matrix,
               x_sym=[Symbol("x$i") for i in 1:size(N_matrix,2)],
               q_sym=[Symbol("t$i") for i in 1:size(L_matrix,1)],
               K_sym=[Symbol("K$i") for i in 1:size(N_matrix,1)])
end

function extract_rays(model, is_bnc)
    rays_out = Vector{Vector{Float64}}()
    regimes = is_bnc ? get_bnc_regimes(model; singular=false) : get_regimes(model; singular=false)
    chart = is_bnc ? :wKk : :qK
    for rgm in regimes
        try
            poly = get_polyhedron(rgm; chart=chart)
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
    return regimes, rays_out
end

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

app = App() do
    seq_field = Bonito.TextField("TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC")
        n_field = Bonito.TextField("1,1,-1,0,0 ; 0,1,0,1,-1")
    l_field = Bonito.TextField("1,0,1,0,0 ; 0,1,1,0,1 ; 0,0,0,1,1")
    seq_btn = Bonito.Button("1. Look up sequence")
    bind_btn = Bonito.Button("2. Run binding circuit")
    cat_btn = Bonito.Button("3. Run catalysis circuit")
    dl_btn = Bonito.Button("Download last result")
    result = Bonito.Observable("Click a button to run.")
        save_status = Bonito.Observable("")
    last_result = Ref("")

    on(seq_btn) do _
        save_status[] = ""
        try
            seq = uppercase(strip(seq_field.value[]))
            if !haskey(PROMOTER_DB, seq)
                result[] = "Unknown sequence."
                update_cone!([])
                return
            end
            name, topology = PROMOTER_DB[seq]
            N_str, L_str = TOPOLOGY_N[topology]
            model = build_binding_model(N_str, L_str)
            regimes, rays = extract_rays(model, false)
                        vols = get_volumes(model; reltol=0.01, time_limit=15.0)
            r = sum(v.mean for v in vols)
            n_pos = count(v -> v.mean > 0, vols)
            txt = "Sequence: " * name * "\nTopology: " * topology *
                  "\nR-index: " * string(round(r, digits=4)) *
                  "\nRegimes: " * string(n_pos) *
                  "\nRays: " * string(length(rays))
            result[] = txt
            last_result[] = txt
            update_cone!(rays)
        catch e
            result[] = "Error: " * string(e)
            update_cone!([])
        end
    end

    on(bind_btn) do _
        save_status[] = ""
        try
            model = build_binding_model(n_field.value[], l_field.value[])
            regimes, rays = extract_rays(model, false)
            vols = get_volumes(model; reltol=0.01, time_limit=15.0)
            r = sum(v.mean for v in vols)
            n_pos = count(v -> v.mean > 0, vols)
            txt = "Custom binding circuit\nR-index: " * string(round(r, digits=4)) *
                  "\nRegimes: " * string(n_pos) *
                  "\nRays: " * string(length(rays))
            result[] = txt
            last_result[] = txt
            update_cone!(rays)
        catch e
            result[] = "Error: " * string(e)
            update_cone!([])
        end
    end

    on(cat_btn) do _
        save_status[] = ""
        try
            model = build_catalysis_model()
                        full_dim = get_bnc_regimes(model; singular=false)
            stable = filter(is_stable, full_dim)
            valid_A = filter(r -> adaptation_signature(r).A_active, stable)
            valid_tA = filter(r -> adaptation_signature(r).tA_active, stable)
            shared = intersect(valid_A, valid_tA)
            r_A = sum(v.mean for v in get_volumes(valid_A; reltol=0.01, time_limit=30.0))
            r_tA = sum(v.mean for v in get_volumes(valid_tA; reltol=0.01, time_limit=30.0))
            r_sh = sum(v.mean for v in get_volumes(shared; reltol=0.01, time_limit=30.0))
            txt = "Competitive adaptation (4 binding + 4 catalysis)\n" *
                  "Full-dim: " * string(length(full_dim)) *
                                    "  Stable: " * string(length(stable)) *
                  "\nValid A*: " * string(length(valid_A)) *
                  "  Valid tA*: " * string(length(valid_tA)) *
                  "  Shared: " * string(length(shared)) *
                  "\n\nR(A*) = " * string(round(r_A, digits=6)) *
                  "\nR(tA*) = " * string(round(r_tA, digits=6)) *
                  "\nR(shared) = " * string(round(r_sh, digits=6)) *
                  "\n\nPaper: 0.019, 0.033, 0.014"
            result[] = txt
            last_result[] = txt
            rays = Vector{Vector{Float64}}()
                        for rgm in union(valid_A, valid_tA)
                try
                    poly = get_polyhedron(rgm; chart=:wKk)
                    for ray in collect(Polyhedra.rays(Polyhedra.vrep(poly)))
                                            v = ray.a
                        if any(x -> abs(x) > 1e-9, v) && length(v) >= 3
                            push!(rays, v[1:3])
                        end
                    end
                catch
                    continue
                end
            end
            update_cone!(rays)
        catch e
            result[] = "Error: " * string(e)
            update_cone!([])
        end
    end

    on(dl_btn) do _
        if isempty(last_result[])
            save_status[] = "Nothing to download yet."
        else
            path = "C:\\bmac-julia\\bac_result_$(rand(1000:9999)).txt"
            write(path, last_result[])
            save_status[] = "Saved to: " * path
        end
    end

    DOM.div(
        DOM.h1("BAC-Explorer v18"),
        DOM.p("Three buttons. No dropdown."),
        DOM.hr(),
        DOM.h2("A. Sequence lookup"),
        seq_field,
        seq_btn,
        DOM.hr(),
        DOM.h2("B. Custom binding circuit"),
        DOM.p("N:"), n_field,
        DOM.p("L:"), l_field,
        bind_btn,
        DOM.hr(),
        DOM.h2("C. Catalysis circuit (competitive adaptation)"),
        cat_btn,
        DOM.hr(),
        DOM.h2("Result"),
        DOM.pre(result),
        DOM.h2("Download"),
        dl_btn, DOM.p(save_status),
        DOM.hr(),
        DOM.h2("Regime cone"),
        cone_fig,
    )
end

server = Bonito.Server(app, "127.0.0.1", 8180)

