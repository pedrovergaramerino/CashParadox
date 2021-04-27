push!(LOAD_PATH,"../src/")
using Documenter, CashParadox

makedocs(modules = [CashParadox], sitename = "CashParadox.jl")

deploydocs(repo = "github.com/pedrovergaramerino/CashParadox.jl.git", devbranch = "main")
