function classify_promoter(name::String, regulator::String)
    if regulator == "Constitutive"
        return (topology = :constitutive, N = [1  1 -1], L = [1  0  1; 0  1  1], x_sym = [:RNAP, :P, :C_P], q_sym = [:tRNAP, :tP], K_sym = [:K_P])
    elseif regulator in ["LacI", "TetR"]
        return (topology = :repressed_competitive, N = [1  1 -1  0  0; 0  1  0  1 -1], L = [1  0  1  0  0; 0  1  1  0  1; 0  0  0  1  1], x_sym = [:RNAP, :P, :C_P, :Repressor, :C_R], q_sym = [:tRNAP, :tP, :tRep], K_sym = [:K_RNAP, :K_Rep])
    elseif regulator == "AraC"
        return (topology = :activated_sequential, N = [1  1 -1  0  0; 0  0  1  1 -1], L = [1  0  1  0  0; 0  1  1  1  1; 0  0  0  1  1], x_sym = [:AraC, :P, :C_A, :RNAP, :C_RNAP], q_sym = [:tAraC, :tP, :tRNAP], K_sym = [:K_AraC, :K_RNAP])
    else
        error("Unknown regulator: $regulator")
    end
end
