using Bonito

app = App() do
    n_field = Bonito.TextField("")
    btn = Bonito.Button("Parse N")
    result = Bonito.Observable("Enter a matrix, then click Parse N.")

    on(btn.value) do _
        raw = n_field.value[]
        if isempty(strip(raw))
            result[] = "Please enter an N matrix."
            return
        end
        try
            rows = split(strip(raw), ';')
            matrix = [parse.(Int, split(strip(r))) for r in rows]
            n_rows = length(matrix)
            n_cols = length(matrix[1])
            display = join(["[" * join(row, " ") * "]" for row in matrix], "\n")
            result[] = "Parsed N matrix: " * string(n_rows) * "x" * string(n_cols) * "\n\n" * display
        catch e
            result[] = "Parse error: " * string(e)
        end
    end

    DOM.div(
        DOM.h1("BAC-Explorer v6"),
        DOM.p("Step 2: Parse a full N matrix"),
        DOM.hr(),
        DOM.p("Enter the N matrix. Separate rows with semicolons."),
        DOM.p("Example 1 reaction: 1 1 -1"),
        DOM.p("Example 2 reactions: 1 1 -1 ; 0 0 1"),
        n_field,
        btn,
        DOM.hr(),
        DOM.h2("Result"),
        DOM.pre(result),
    )
end

server = Bonito.Server(app, "127.0.0.1", 8099)
