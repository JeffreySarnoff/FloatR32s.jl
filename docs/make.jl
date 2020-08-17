using Documenter, Robust32s

makedocs(
    pages  = Any[
        "Overview" => "index.md",
    ]   
)
        
deploydocs(
    repo = "github.com/JeffreySarnoff/Robust32s.jl.git",
)
