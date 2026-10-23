using Pkg
Pkg.activate("C:\\BindingAndCatalysis.jl")
using BindingAndCatalysis, LinearAlgebra, SparseArrays, Bonito
Base.muladd(a::Real, b::BindingAndCatalysis.ExactLogExpr, c::BindingAndCatalysis.ExactLogExpr) = (a * b) + c

app = App() do
    n_field = Bonito.TextField("1,1,-1")
        l_field = Bonito.TextField("1,0,1 ; 0,1,1")
    btn = Bonito.Button("Compute R-index")
    result = Bonito.Observable("Enter N and L. Commas, spaces, semicolons all work.")
    on(btn) do _
        try
            parse_row(s) = parse.(Int, split(strip(s), r"[,\s]+"))
                        n_rows = split(strip(n_field.value[]), r"[;\n]+")
            l_rows = split(strip(l_field.value[]), r"[;\n]+")
            N_matrix = reduce(vcat, [reshape(parse_row(r), 1, :) for r in n_rows if !isempty(strip(r))])
            L_matrix = reduce(vcat, [reshape(parse_row(r), 1, :) for r in l_rows if !isempty(strip(r))])
            n_species = size(N_matrix, 2)
            x_sym = [Symbol("x$i") for i in 1:n_species]
            q_sym = [Symbol("t$i") for i in 1:size(L_matrix, 1)]
            K_sym = [Symbol("K$i") for i in 1:size(N_matrix, 1)]
            model = Bnc(N=N_matrix, L=L_matrix, x_sym=x_sym, q_sym=q_sym, K_sym=K_sym)
            volumes = get_volumes(model; reltol=0.01, time_limit=15.0)
            r_total = sum(v.mean for v in volumes)
            n_pos = count(v -> v.mean > 0, volumes)
            detail = join([string(round(v.mean, digits=4)) for v in volumes], ", ")
            result[] = "Total R-index: " * string(round(r_total, digits=4)) *
                       "\nPositive regimes: " * string(n_pos) *
                       "\nVolumes: " * detail
        catch e
            result[] = "Error: " * string(e)
        end
    end
    DOM.div(
        DOM.h1("BAC-Explorer v11"),
        DOM.p("Accepts commas, spaces, or semicolons."),
        DOM.p("N matrix:"), n_field,
        DOM.p("L matrix:"), l_field,
        btn,
        DOM.hr(),
        DOM.h2("Result"),
        DOM.pre(result),
    )
end

server = Bonito.Server(app, "127.0.0.1", 8130)
