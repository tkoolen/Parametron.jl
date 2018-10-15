# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

using Documenter, Parametron

makedocs(
    modules = [Parametron, Parametron.Functions],
    checkdocs = :exports,
    format = :html,
    root = @__DIR__,
    sitename ="Parametron.jl",
    authors = "Twan Koolen and contributors.",
    pages = [
        "Home" => "index.md",
        "API" => [
            "api/functions.md",
            "api/parameter.md",
            "api/lazyexpression.md",
            "api/model.md",
            "api/debug.md",
        ]
    ],
    html_prettyurls = parse(Bool, get(ENV, "CI", "false")),
)

deploydocs(
    repo = "github.com/tkoolen/Parametron.jl.git",
    target = "build"
)
