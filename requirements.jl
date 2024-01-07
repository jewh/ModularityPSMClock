# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# installs the required packages for the project

using Pkg 

dependencies = [
    "CSV",
    "DataFrames",
    "StaticArrays",
    "LinearAlgebra",
    "Distributions",
    "PyPlot",
    "Plots",
    "Random",
    "CairoMakie",
    "Glob",
    "ClusterManagers",
    "Distributed",
    "Dictionaries",
    "Statistics",
    "BenchmarkTools"
]

Pkg.add(dependencies)
