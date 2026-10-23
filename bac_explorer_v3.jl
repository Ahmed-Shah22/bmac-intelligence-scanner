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
        DOM.h1("BAC-Explorer v3"),
        DOM.p("Select a circuit:"),
        dropdown,
        DOM.hr(),
        DOM.h2("R-index results"),
        content,
        DOM.hr(),
        DOM.h2("Regime cone (Competitive Adaptation, 11 rays)"),
        DOM.p("Interactive: drag to rotate. Scroll to zoom."),
        cone_fig
    )
end

server = Bonito.Server(app, "127.0.0.1", 8093)
