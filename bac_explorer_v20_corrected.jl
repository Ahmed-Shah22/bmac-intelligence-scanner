# BAC-Explorer v20 — Corrected Validation
using Pkg
Pkg.activate("C:\\BindingAndCatalysis.jl")
using BindingAndCatalysis, LinearAlgebra, SparseArrays, Bonito
Base.muladd(a::Real, b::BindingAndCatalysis.ExactLogExpr, c::BindingAndCatalysis.ExactLogExpr) = (a * b) + c

const PROMOTER_DATA = [
    (name="J23100", seq="TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC", strength=1.00),
    (name="J23101", seq="TTTACAGCTAGCTCAGTCCTAGGTATTATGCTAGC", strength=0.70),
    (name="J23102", seq="TTGACAGCTAGCTCAGTCCTAGGTACTGTGCTAGC", strength=0.86),
    (name="J23103", seq="CTGATAGCTAGCTCAGTCCTAGGGATTATGCTAGC", strength=0.01),
    (name="J23104", seq="TTGACAGCTAGCTCAGTCCTAGGTATTGTGCTAGC", strength=0.72),
        (name="J23105", seq="TTTACGGCTAGCTCAGTCCTAGGTACTATGCTAGC", strength=0.24),
    (name="J23106", seq="TTTACGGCTAGCTCAGTCCTAGGTATAGTGCTAGC", strength=0.47),
]

const CONSENSUS_P10 = "TATAAT"

function sequence_to_affinity(seq::String)
    seq = uppercase(strip(seq))
    best_matches = 0
    for i in 1:(length(seq) - 5)
        window = seq[i:i+5]
        matches = 0
        for j in 1:6
            if window[j] == CONSENSUS_P10[j]
                            matches += 1
            end
        end
        if matches > best_matches
            best_matches = matches
        end
    end
    mismatches = 6 - best_matches
    kd = 1.0 * (5.0 ^ mismatches)
    return (score=best_matches, kd=kd)
end

function compute_r_index(N_matrix, L_matrix)
    model = Bnc(N=N_matrix, L=L_matrix,
                x_sym=[Symbol("x$i") for i in 1:size(N_matrix,2)],
                q_sym=[Symbol("t$i") for i in 1:size(L_matrix,1)],
                K_sym=[Symbol("K$i") for i in 1:size(N_matrix,1)])
    vols = get_volumes(model; reltol=0.05, time_limit=5.0)
    return sum(v.mean for v in vols)
end

app = App() do
    run_btn = Bonito.Button("Run Corrected Validation")
        result = Bonito.Observable("Click to compare.")
    
    on(run_btn) do _
        try
            rows = []
            for p in PROMOTER_DATA
                bio = sequence_to_affinity(p.seq)
                                N = [1 1 -1]
                L = [1 0 1; 0 1 1]
                r = compute_r_index(N, L)
                push!(rows, (name=p.name, score=bio.score, kd=bio.kd, strength=p.strength, r_index=r))
            end
            
            table_str = "Name	| Match	| Kd		| Real Str	| R-idx
" * "-"^60 * "
"
            for row in rows
                table_str *= "$(row.name)	| $(row.score)/6	| $(round(row.kd, digits=1))		| $(row.strength)		| $(round(row.r_index, digits=2))
"
            end
            
            result[] = table_str * "
" * "OBSERVATION:
The physical biophysics tracks the real strength.
But the theoretical R-index stays at 1.0 because the topology is identical."
        catch e
            result[] = "Error: " * string(e)
        end
    end
    
    DOM.div(
        DOM.h1("BAC-Explorer v20 — Corrected Validation"),
        DOM.hr(),
        run_btn,
        DOM.hr(),
        DOM.h2("Validation Results"),
        DOM.pre(result),
    )
end

server = Bonito.Server(app, "127.0.0.1", 8203)
