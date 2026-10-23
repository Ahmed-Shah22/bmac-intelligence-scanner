using Bonito

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

app = App() do
    DOM.div(
        DOM.h1("BAC-Explorer v2"),
        DOM.p("Biocircuit Realizability Analysis"),
        DOM.hr(),
        DOM.h2("Competitive Adaptation Circuit"),
        DOM.ul([DOM.li(row) for row in CIRCUITS["Competitive Adaptation"]]),
        DOM.hr(),
        DOM.h2("Michaelis-Menten"),
        DOM.ul([DOM.li(row) for row in CIRCUITS["Michaelis-Menten"]]),
        DOM.hr(),
        DOM.h2("Hill (n=2, sequential binding)"),
        DOM.ul([DOM.li(row) for row in CIRCUITS["Hill (n=2, sequential)"]]),
    )
end

server = Bonito.Server(app, "127.0.0.1", 8092)
