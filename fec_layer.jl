using Pkg
Pkg.activate("C:\\\\BindingAndCatalysis.jl")
using BindingAndCatalysis, LinearAlgebra, SparseArrays, Bonito, DifferentialEquations
Base.muladd(a::Real, b::BindingAndCatalysis.ExactLogExpr, c::BindingAndCatalysis.ExactLogExpr) = (a * b) + c

function fec_dynamics!(du, u, p, t)
    K_atp, n_hill = p
    # Glucose is constant (external supply). We only model ATP/ADP dynamics.
    v_glycolysis = 1.0 * (u[1]^n_hill / (K_atp^n_hill + u[1]^n_hill))
    du[1] = v_glycolysis - 0.5*u[1] # ATP production and consumption
    du[2] = 0.5*u[1] - v_glycolysis # ADP balances ATP
end

app = App() do
    n_slider = Bonito.Slider(1:1:5, value=2, label="ATP Regulation Exponent (n)")
    run_btn = Bonito.Button("Run FEC Simulation")
    result = Bonito.Observable("Adjust the exponent and click to simulate.")
    on(run_btn) do _
        try
            n_val = n_slider.value[]
            p = (0.5, n_val)
            u0 = [1.0, 0.1]
            tspan = (0.0, 50.0)
            prob = ODEProblem(fec_dynamics!, u0, tspan, p)
            sol = solve(prob, Tsit5(), saveat=0.1)
            atp_traj = sol[1, :]
            # Proper oscillation detection: count peaks in the last half of the simulation
            tail = atp_traj[end÷2:end]
            peaks = 0
            for i in 2:length(tail)-1
                if tail[i] > tail[i-1] && tail[i] > tail[i+1]
                    peaks += 1
                end
            end
            is_oscillating = peaks > 3
            status_str = is_oscillating ? "YES (Glycolytic Oscillation)" : "NO (Stable Steady State)"
            out = "── FEC: DYNAMIC METABOLISM ──\n\n"
            out *= "Regulation Exponent (n): " * string(n_val) * "\n"
            out *= "Final ATP Concentration: " * string(round(atp_traj[end], digits=2)) * "\n"
            out *= "Oscillation Detected: " * status_str * "\n"
            out *= "----------------------------------------\n\n"
            out *= "The FEC theory states that the exponent 'n' is the control knob.\n"
            out *= "At high n, the system exhibits limit-cycle oscillations.\n"
            out *= "At low n, it relaxes to a stable steady state."
            result[] = out
        catch e
            result[] = "Error: " * string(e)
        end
    end
    DOM.div(
        DOM.style("body { background: radial-gradient(circle at top left, #2A0A0A, #0B0B12); margin: 0; padding: 0; } .canvas { min-height: 100vh; padding: 80px 20px; display: flex; justify-content: center; font-family: 'Helvetica Neue', sans-serif; color: #E8D8B1; } .frame { max-width: 680px; width: 100%; } .header-sigil { display: flex; justify-content: space-between; align-items: center; border-bottom: 1px solid #5A3A22; padding-bottom: 20px; margin-bottom: 60px; font-size: 24px; } .sigil-group { display: flex; gap: 20px; color: #D4AF37; } .caltech-mark { color: #FF6C0C; font-weight: bold; } .control-mark { color: #00A86B; } h1.title { font-family: 'Times New Roman'; font-size: 52px; font-weight: 400; line-height: 1.05; margin: 0 0 32px 0; color: #FDFBF7; } .subtitle { font-size: 18px; line-height: 1.7; color: #B8A99A; margin: 0 0 60px 0; font-weight: 300; } .art-btn { background-color: #D4AF37 !important; color: #0B0B12 !important; border: none !important; padding: 16px 40px !important; font-size: 12px !important; letter-spacing: 4px !important; text-transform: uppercase !important; cursor: pointer !important; font-weight: bold !important; transition: all 0.4s ease !important; } .art-btn:hover { background-color: #00A86B !important; color: #FDFBF7 !important; } .placard { background: transparent !important; border: 1px solid #3A2A1A !important; padding: 32px !important; font-family: 'Courier New', Courier, monospace !important; font-size: 14px !important; line-height: 1.8 !important; color: #D4AF37 !important; margin-top: 40px !important; box-shadow: inset 0 0 20px rgba(0,0,0,0.5); white-space: pre-wrap; } input[type=\"range\"] { -webkit-appearance: none; width: 100%; height: 4px; background: #D4AF37; outline: none; opacity: 0.7; transition: .2s; margin-bottom: 30px; } input[type=\"range\"]:hover { opacity: 1; } input[type=\"range\"]::-webkit-slider-thumb { -webkit-appearance: none; appearance: none; width: 16px; height: 16px; background: #00A86B; cursor: pointer; border-radius: 50%; }"),
        DOM.div(class="canvas", DOM.div(class="frame", DOM.div(class="header-sigil", DOM.div(class="sigil-group", DOM.span("❖"), DOM.span("∞", class="control-mark"), DOM.span("☰")), DOM.div("BMAC LAB // FEC LAYER", style="font-size: 11px; letter-spacing: 4px; color: #8A7A6A;"), DOM.div(class="sigil-group", DOM.span("✦", class="caltech-mark"), DOM.span("◈"), DOM.span("❁"))), DOM.h1("Flux Exponent Control.", class="title"), DOM.p("The R-index analyzes structural robustness. FEC simulates dynamic metabolism. Here, the binding network regulates the *exponent* of the flux. Adjust the ATP regulation exponent and observe the transition from stable steady state to glycolytic oscillation.", class="subtitle"), DOM.h4("REGULATION EXPONENT (n):", style="font-size: 11px; letter-spacing: 2px; color: #8A7A6A;"), n_slider, run_btn, DOM.div(DOM.pre(result, class="placard"))))
    )
end
server = Bonito.Server(app, "127.0.0.1", 8242)
