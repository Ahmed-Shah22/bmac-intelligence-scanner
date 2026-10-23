using Bonito

app = App() do
    n_field = Bonito.TextField("")
    btn = Bonito.Button("Compute")
    result = Bonito.Observable("Enter a row of N, then click Compute.")

    on(btn.value) do _
        raw = n_field.value[]
        if isempty(strip(raw))
            result[] = "Please enter a row of the N matrix."
            return
        end
        try
            vals = parse.(Int, split(strip(raw)))
            result[] = "Parsed N = " * string(vals) * "  (length " * string(length(vals)) * ")"
        catch e
            result[] = "Parse error: " * string(e)
        end
    end

    DOM.div(
        DOM.h1("BAC-Explorer v5"),
        DOM.p("Step 1: User-input form"),
        DOM.hr(),
        DOM.p("Enter a row of the N matrix as space-separated integers."),
        DOM.p("Example: 1 1 -1  (reactant +1, reactant +1, product -1)"),
        n_field,
        btn,
        DOM.hr(),
        DOM.h2("Result"),
        DOM.p(result),
    )
end

server = Bonito.Server(app, "127.0.0.1", 8098)
