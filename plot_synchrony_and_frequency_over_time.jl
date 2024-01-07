# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# plot synchrony and frequency over time for specific simulations

using CairoMakie, DataFrames, CSV
include("kuramoto.jl")

tissue = TissueParams(
    11.0, # cell diameter in μm .
    8.71, # coefficient for intercellular force .
    1.67, # cell advection velocity in anterior .
    3.0, # cell advection speed in posterior .
    0.7, # % x coordinate dividing the PSM into anterior and posterior .
    300.0 + 60.0 + 25.0, # length of the tissue in μm. NOTE - this is dynamically modified if abs(uₐ) > 0
    0., # x-coordinate of the most anterior cell .
    0.4, # length scale of the cell motility gradient. Inf64 to set motility constant across tissue .
    1.0, # intrinsic cell motility in the posterior of the tissue. This is from Uriu et al. 2017, not 2021! (1.0 in latter) .
    3.0, # coefficient for shape of motility gradient (needs to be 1.0 for no gradient) .
    0.026, # phase noise for cell directionality (/min) .
    20.0, # coefficient for boundary force .
    1.0, # length-scale of boundary force .
    300.0, # x coord of centre of PSM torus .
    60.0 + 25., # y coord of centre of PSM torus .
    25., # z coord of centre of PSM torus .
    25.0, # radius of PSM cylinder . NOTE - this is dynamically modified if abs(mᵣ) > 0
    60.0, # radius of PSM torus . 
    0.0015, # desired density (from C code). NOTE - this is dynamically modified if abs(md) > 0
    2 * 883 + 555, # number of cells in tissue
    2 * 883 + 555, # total number of cells that have been in the tissue
    187.5 - 15., # duration of growth phase (mins)
    15., # duration of M phase (mins)
    100., # anterior limit of cell addition .
    0.75, # coupling strength of cell cycle to clock
    1.1, # increase in μm in cell diameter at posterior vs anterior
    325., # length intercept. NOTE - this is NOT dynamically modified if abs(uₐ) > 0
    0.0, # rate at which xₐ moves to the posterior.
    0.002123, # density intercept . NOTE - this is NOT dynamically modified if abs(md) > 0
    0.0, # rate at which density ρ₀ changes 
    27.6, # radius intercept. NOTE - this is NOT dynamically modified if abs(mᵣ) > 0
    0.0, # rate at which radius changes
    253.6, # time at which the PSM length begins to decrease (Warning: will yield errors for infs or nans)
    300.0 + 60.0 + 25.0, # length of lattice in x
    2 * 60. + 2 * 25., # length of lattice in y
    2 * 25., # length of lattice in z
    0.0, # anterior limit of cell division
    9.2, # initial diameter of cells before they start to shrink
    0.0, # rate of cell dimaeter decrease
    11.0, # side of the lattice cubes. This should be equal to the size of the largest cell in the simulation.
    536.4 # time at which we stop shrinking the cells in the PSM
    )

phase = PhaseParams(
    0.07,
    0.0013,
    0.2094,
    0.66, # 0.66 by default
    3.07,
    2101,
    )

fname = raw"data\20221201_frequency_tracking_tissue_simulation_randomnewcells_seed0.csv"

df = DataFrame(CSV.File(fname, delim=","))

# only get the anterior left-hand side
df = df[df[!, Symbol("Position X Reference Frame")] .< 25 + 11, :]
df = df[df[!, Symbol("Position X Reference Frame")] .>= 25, :]
df = df[df[!, Symbol("Position Y Reference Frame")] .< 85, :]

# now get a trace for mean frequency and the synchrony at anterior over time
T = sort(collect(Set(df[!, :Time_Mins])))
r = Vector{Float64}(undef, length(T))
f = Vector{Float64}(undef, length(T))

for (i, t) in enumerate(T)

    dt = df[df[!, :Time_Mins] .== t, :]

    r[i] = phase_order(dt[!, :Phase_sim])
    f[i] = Statistics.mean(dt[!, Symbol("dθ/dt")])

end


fig = Figure()

ax1 = Axis(fig[1,1], ylabel=L"r")
ax2 = Axis(fig[2,1], xlabel = "t (min)", ylabel="dθ/dt (min⁻¹)")

lines!(ax1, T, r; color=Makie.wong_colors()[1])
xlims!(ax1, 0, 1000)
ylims!(ax1, 0, 1.1)

lines!(ax2, T, f; color=Makie.wong_colors()[1])
xlims!(ax2, 0, 1000)
hlines!(ax2, [avgω(tissue.xa, phase, tissue)], linestyle=:dash, color=:black)
# text!(ax2, 900, 0.1415, text=L"\frac{1}{d_{c}} \int_{x_{a}}^{x_{a} + d_{c}} \omega(x) dx")

save("figures/20231215_traces_of_synchrony_and_frequency.svg", fig)
save("figures/20231215_traces_of_synchrony_and_frequency.pdf", fig)


fig
