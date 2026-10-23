using BindingAndCatalysis

repressed_model = Bnc(
    N = [ 1   1  -1   0   0;
          0   1   0   1  -1],
    L = [ 1   0   1   0   0;
          0   1   1   0   1;
          0   0   0   1   1],
    x_sym = [:RNAP, :P, :C_P, :Repressor, :C_R],
    q_sym = [:tRNAP, :tP, :tRep],
    K_sym = [:K_RNAP, :K_Rep],
)

println("Conservation laws:")
show_conservation(repressed_model)
println()

regimes_rep = get_regimes(repressed_model; singular=false)
println("Number of full-dimensional regimes: ", length(regimes_rep))
println()

volumes_rep = get_volumes(repressed_model; reltol=0.01, time_limit=30.0)
println("Per-regime volumes:")
for (i, v) in enumerate(volumes_rep)
    println("  regime ", i, ": ", v)
end
