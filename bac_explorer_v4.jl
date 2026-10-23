using Bonito, WGLMakie, Makie
WGLMakie.activate!()

const CIRCUITS = Dict(
    "Competitive Adaptation" => [
        "R(free A-star adaptation) = 0.0194 +/- 0.47%",
        "R(total tA-star adaptation) = 0.0329 +/- 0.36%",
        "R(shared, both) = 0.0141 +/- 0.55%",
    ],
    "Michaelis-Menten" => [
        "R(linear, low substrate) = 0.3335 +/- 0.88%",
        "R(saturated, high substrate) = 0.3321 +/- 0.88%",
        "R(third regime) = 0.3345 +/- 0.87%",
        "R(equality regime) = 0.0000",
    ],
    "Hill (n=2, sequential)" => [
        "R(C2 ~ tS^2 * tE) = 0.2624 +/- 1.04%",
        "R(C2 ~ tS^2 / tE) = 0.2917 +/- 0.97%",
        "R(C2 ~ tS * tE) = 0.0824 +/- 2.07%",
        "R(C2 ~ tS) = 0.0623 +/- 2.40%",
        "R(C2 ~ tE, saturated) = 0.3013 +/- 0.94%",
    ],
    "Promoter: Constitutive" => [
        "Topology: 1 binding reaction (RNAP + P <-> C_P)",
        "Total regimes: 3",
        "Positive-volume regimes: 3 (each ~ 0.333)",
        "Sum of volumes: 1.000",
        "Interpretation: R-index for 'promoter bound' = 1/3",
        "Applies to: J23100, J23101, J23102, J23103, J23104, J23105, J23106",
    ],
    "Promoter: Repressed (competitive)" => [
        "Topology: 2 competing reactions (RNAP and Repressor bind same promoter)",
        "Total regimes: 12",
        "Positive-volume regimes: 8",
        "Largest regime: 0.244",
        "Smallest positive regime: 0.023",
        "Sum of volumes: 1.000",
        "Applies to: lacUV5 (LacI), pTet (TetR)",
    ],
    "Promoter: Activated (sequential)" => [
        "Topology: 2 sequential reactions (AraC binds, then RNAP)",
        "Total regimes: 12",
        "Positive-volume regimes: 6",
        "Half of parameter space collapses to zero (requires exact equalities)",
        "Largest regime: 0.150",
        "Sum of positive volumes: 0.500",
        "Applies to: pBAD (AraC)",
    ],
    "Promoter + RBS (combined)" => [
        "Topology: 2 independent binding reactions",
        "Total regimes: 16",
        "Positive-volume regimes: 9 (each ~ 0.111)",
        "Sum of volumes: 1.000",
        "Interpretation: independence reproduces (1/3) x (1/3) = 1/9",
    ],
)

const RAY_MATRIX = [
    -1.0  0.0 -1.0  0.0 -1.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0;
     0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0 -1.0  0.0  0.0  0.0;
     1.0  1.0  1.0  0.0  1.0  1.0  1.0  0.0 -1.0  0.0 -1.0  0.0;
     1.0  0.0  1.0  0.0  1.0  1.0  1.0  0.0 -1.0  0.0  0.0  0.0;
     0.0  1.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0;
     0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0  0.0;
     0.0  0.0  1.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
     0.0  0.0 -1.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0;
     1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0;
     1.0  0.0  1.0  0.0  1.0  0.0  0.0  0.0 -1.0  0.0  0.0  0.0;
    -0.0 -0.0 -0.0 -1.0 -0.0 -0.0 -0.0  0.0 -0.0  0.0 -1.0  0.0;
]

function build_cone_figure()
    fig = Figure(size=(700, 500))
    ax = LScene(fig[1, 1]; show_axis=true)
    xs = Float64[]; ys = Float64[]; zs = Float64[]
    for i in 1:size(RAY_MATRIX, 1)
        push!(xs, 0.0); push!(ys, 0.0); push!(zs, 0.0)
        push!(xs, RAY_MATRIX[i, 1])
        push!(ys, RAY_MATRIX[i, 2])
        push!(zs, RAY_MATRIX[i, 3])
    end
    lines!(ax, xs, ys, zs; color=:steelblue, linewidth=2)
    return fig
end

cone_fig = build_cone_figure()

app = App() do
    dropdown = Bonito.Dropdown(collect(keys(CIRCUITS)))
    content = map(dropdown.value) do name
        DOM.div(
            DOM.h2(name),
            DOM.ul([DOM.li(row) for row in CIRCUITS[name]])
        )
    end
    DOM.div(
        DOM.h1("BAC-Explorer v4"),
        DOM.p("Biocircuit Realizability Analysis -- with Phase 1 promoter dataset"),
        DOM.p("Select a circuit or promoter topology:"),
        dropdown,
        DOM.hr(),
        DOM.h2("R-index results"),
        content,
        DOM.hr(),
        DOM.h2("Regime cone (Competitive Adaptation reference, 11 rays)"),
        DOM.p("Interactive: drag to rotate. Scroll to zoom."),
        cone_fig
    )
end

server = Bonito.Server(app, "127.0.0.1", 8095)
