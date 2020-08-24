using Documenter, Robust32s

makedocs(
    modules = [Robust32s],
    sitename = "Robust32s",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages  = Any[
        "Overview" => "index.md",
        "How It Works" => "howitworks.md"
    ]   
)
        
deploydocs(
    repo = "github.com/JeffreySarnoff/Robust32s.jl.git",
    target = "build"
)
