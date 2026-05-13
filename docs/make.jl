using Documenter
using SchroedingerEquation
using LinearAlgebra

makedocs(
    sitename = "SchroedingerEquation.jl",
    format = Documenter.HTML(),
    modules = [SchroedingerEquation],
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Basis & Grids" => "basis.md",
            "Hamiltonians" => "hamiltonian.md",
            "Solvers" => "solvers.md",
        ],
        "Python Support" => "python.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/physarief78/SchroedingerEquation.jl.git",
)
