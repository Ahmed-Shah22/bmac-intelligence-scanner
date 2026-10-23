using Bonito

app = App() do
    DOM.div(
        DOM.h1("BAC-Explorer"),
        DOM.p("Biocircuit Realizability Analysis"),
        DOM.hr(),
        DOM.h2("Competitive Adaptation Circuit"),
        DOM.p("R-index results (Monte Carlo, reltol 1%):"),
        DOM.ul(
            DOM.li("R(free A⋆ adaptation) = 0.0194  ± 0.47%"),
            DOM.li("R(total tA⋆ adaptation) = 0.0329 ± 0.36%"),
            DOM.li("R(shared, both) = 0.0141 ± 0.55%"),
        ),
        DOM.hr(),
        DOM.p("Flagship regime: index 28, volume 0.00459"),
        DOM.p("Interpretation: 1.4% of parameter space produces robust adaptation."),
    )
end

server = Bonito.Server(app, "127.0.0.1", 8090)