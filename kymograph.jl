# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# file for plotting kymographs of phase over time for specific simulations

# Imports
using CairoMakie, LaTeXStrings
include("experiments.jl") # code for simulating
include("save_simulations.jl") # load code for saving data as .csv
include("simulation_analysis.jl") # code for plotting


translate_phase(x::Float64) = 0.5 * (1 - cos(x))

"""
Transform the data into an X x T matrix where X no. x points, T no. time points.

NOTE: This will not work with a simulation where tissue length changes.
"""
function kymograph(df, phase::PhaseParams, tissue::TissueParams; type="frequency")

    @assert type ∈ ("frequency", "synchrony", "phase") "'$(type)' is not a valid argument for 'type'."

    dx = tissue.dc # step for going along the x axis (μm)

    xmin, xmax = tissue.xa, tissue.xa + tissue.L
    xs = xmin:dx:xmax # x range

    ts = unique(df[!, :Time_Mins])

    kmg = Matrix{Float64}(undef, length(ts), length(xs)-1) # matrix which is our kymograph

    for (j, t) in enumerate(ts)

        du = df[df[!, :Time_Mins] .== t, :]


        for i = 2:length(xs)
            dt = du[du[!, Symbol("Position X Reference Frame")] .>= xs[i-1], :]
            dt = dt[dt[!, Symbol("Position X Reference Frame")] .< xs[i], :]

            if type == "frequency"

                kmg[j, i-1] = mean(dt[!, Symbol("dθ/dt")])

            elseif type == "synchrony"

                kmg[j, i-1] = phase_order(dt[!, :Phase_sim])

            elseif type == "phase"

                kmg[j, i-1] = mean(translate_phase.(dt[!, :Phase_sim]))

            end


        end

    end

    return kmg, ts, xs

end



# load in parameters
begin
    # params for default simulation
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
        # 2100,
        2601
        )

    tissue.xa += 25.
    tissue.Xc += 25.
    tissue.Yc += 25.
    tissue.Zc += 25.
    tissue.Lx += 25.
    tissue.Ly += 25.
    tissue.Lz += 25.
end

# code for comparing adding posterior, DP + LV, control kymographs

fig = Figure()

file = raw"data\20221201_frequency_tracking_tissue_simulation_randomnewcells_seed0.csv"

df = DataFrame(CSV.File(file))

# now select out one half of the PSM. Use LHS for ease.
df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

kmg, ts, xs = kymograph(df, phase, tissue; type="synchrony")

dt = ts[2] - ts[1]
dx = xs[2] - xs[1]

ax1 = Axis(fig[1,1], aspect=AxisAspect((dx * length(xs))/(length(ts)*dt)))
ax1.ylabel = L"t \; \text{(min)}"
ax1.xlabel = L"x \; \text{(μm)}"
hm1 = heatmap!(ax1, xs .- tissue.xa, ts, kmg', colorrange=(0,1), colormap = :viridis)

# now posterior-only addition

file = raw"data\20221201_frequency_tracking_adding_posterior_zero_cells_simulation_seed200.csv"

df = DataFrame(CSV.File(file))

# now select out one half of the PSM. Use LHS for ease.
df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

kmg, ts, xs = kymograph(df, phase, tissue; type="synchrony")

dt = ts[2] - ts[1]
dx = xs[2] - xs[1]

ax2 = Axis(fig[1, 2], aspect=AxisAspect((dx * length(xs))/(length(ts)*dt)))
ax2.ylabel = L"t \; \text{(min)}"
ax2.xlabel = L"x \; \text{(μm)}"
heatmap!(ax2, xs .- tissue.xa, ts, kmg', colorrange=(0,1), colormap = :viridis)

# Now DP + LV

file = raw"data\20221201_frequency_tracking_adding_posterior_lateral_zero_cells_simulation_seed300.csv"

df = DataFrame(CSV.File(file))

# now select out one half of the PSM. Use LHS for ease.
df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

kmg, ts, xs = kymograph(df, phase, tissue; type="synchrony")

dt = ts[2] - ts[1]
dx = xs[2] - xs[1]

ax3 = Axis(fig[1, 3], aspect=AxisAspect((dx * length(xs))/(length(ts)*dt)))
ax3.ylabel = L"t \; \text{(min)}"
ax3.xlabel = L"x \; \text{(μm)}"
heatmap!(ax3, xs .- tissue.xa, ts, kmg', colorrange=(0,1), colormap = :viridis)

cbar = Colorbar(fig[1, 4], hm1, label=L"r")
cbar.tellheight = true

ax1.title = "Random"
ax2.title = "DP"
ax3.title = "DP + LV"

# colgap!(fig.layout, 2)
# rowgap!(fig.layout, 10)
# rowsize!(fig.layout, 1, Aspect(1, 1))

fig

save("figures/20230907_celladdition_kymographs.svg", fig)


# # code for making left-right figure for posterior addition of cells
# begin
#     file = raw"data\posterior_simulations\20220719_frequency_tracking_adding_posterior_zero_cells_simulation_seed30_deterministic.csv"

#     df = DataFrame(CSV.File(file))



#     # now select out one half of the PSM. Use LHS for ease.
#     df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

#     Plots.theme(
#         :vibrant; 
#         framestyle=:box
#         )

#     kmg, ts, xs = kymograph(df, phase, tissue; type="synchrony")

#     rfigl = Plots.heatmap(ts, xs[1:end-1], kmg', color=:viridis, colorbartitle="r", clims=(0, 1), colorbar=false)
#     Plots.xlims!((0, 1000))
#     Plots.xlabel!("Time (mins)")
#     Plots.ylabel!("x (μm)")


#     # Now plot the synchrony in anterior over time
#     times = sort(unique(df[!, :Time_Mins]))
#     rs = Vector{Float64}(undef, length(times))

#     for (t, time) in enumerate(times)
#         dt = df[df[!, :Time_Mins] .== time, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .>= tissue.xa, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .< tissue.xa + tissue.dc, :]
#         rs[t] = phase_order(dt[!, :Phase_sim])
#     end

#     rfigl2 = Plots.plot(times, rs, label=false, color="#e37")
#     Plots.xlims!((0, 1000))
#     Plots.ylims!((0, 1))
#     Plots.ylabel!("r")
#     Plots.xlabel!("Time (mins)")

#     # Now right hand side
#     df = DataFrame(CSV.File(file))

#     # now select out one half of the PSM. Use LHS for ease.
#     df = df[df[!, Symbol("Position Y Reference Frame")] .>= tissue.Yc, :]

#     kmg, ts, xs = kymograph(df, phase, tissue; type="synchrony")

#     rfigr = Plots.heatmap(ts, xs[1:end-1], kmg', color=:viridis, colorbartitle="r", clims=(0, 1), colorbar=false)
#     Plots.xlims!((0, 1000))
#     Plots.xlabel!("Time (mins)")
#     Plots.ylabel!("x (μm)")


#     # Now plot the synchrony in anterior over time
#     times = sort(unique(df[!, :Time_Mins]))
#     rs = Vector{Float64}(undef, length(times))

#     for (t, time) in enumerate(times)
#         dt = df[df[!, :Time_Mins] .== time, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .>= tissue.xa, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .< tissue.xa + tissue.dc, :]
#         rs[t] = phase_order(dt[!, :Phase_sim])
#     end

#     rfigr2 = Plots.plot(times, rs, label=false, color="#e37")
#     Plots.xlims!((0, 1000))
#     Plots.ylims!((0, 1))
#     Plots.ylabel!("r")
#     Plots.xlabel!("Time (mins)")

#     # Now create a colorbar to add to the multipanel plot

#     cbar = Plots.scatter(
#         [0, 0], [0, 1], zcolor=[0, 1], clims=(0, 1), xlims=(1.05, 1.1), label="", c=:viridis, 
#         colorbartitle="r", framestyle=:none
#         )

#     l = @layout [Plots.grid(2, 2) a{0.035w}]

#     bigfig = Plots.plot(rfigl, rfigr, rfigl2, rfigr2, cbar, layout=l, size=(1000, 650), leftmargin=2.5mm)

#     Plots.savefig(bigfig, "figures/20220324_frequency_tracking_adding_posterior_zero_cells_simulation_seed30_deterministic_LR_kymograph.svg")
# end

# # code for comparing adding posterior, DP + LV, control kymographs, for frequency
# begin

#     file = raw"data\20220324_frequency_tracking_tissue_simulation_randomnewcells_seed0_deterministic.csv"

#     df = DataFrame(CSV.File(file))

#     # now select out one half of the PSM. Use LHS for ease.
#     df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

#     kmg, ts, xs = kymograph(df, phase, tissue; type="frequency")

#     for (i, x) in enumerate(xs[1:end-1])

#         kmg[:, i] .-= avgω(x, phase, tissue)

#     end

#     ffigl = Plots.heatmap(ts, xs[1:end-1], kmg', color=:balance, colorbartitle="r", clims=(-0.02, 0.02), colorbar=false)
#     Plots.xlims!((0, 1000))
#     Plots.xlabel!("Time (mins)")
#     Plots.ylabel!("x (μm)")
#     Plots.title!("Random")

#     # Now plot the frequency in anterior over time
#     times = sort(unique(df[!, :Time_Mins]))
#     fs = Vector{Float64}(undef, length(times))

#     for (t, time) in enumerate(times)
#         dt = df[df[!, :Time_Mins] .== time, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .>= tissue.xa, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .< tissue.xa + tissue.dc, :]
#         fs[t] = mean(dt[!, Symbol("dθ/dt")])
#     end

#     ffigl2 = Plots.plot(times, fs, label=false, color="#e37")
#     Plots.xlims!((0, 1000))
#     Plots.ylims!((0.14, 0.16))
#     Plots.ylabel!("dθ/dt")
#     Plots.xlabel!("Time (mins)")

#     # now posterior-only addition

#     file = raw"data\posterior_simulations\20220719_frequency_tracking_adding_posterior_zero_cells_simulation_seed30_deterministic.csv"

#     df = DataFrame(CSV.File(file))

#     # now select out one half of the PSM. Use LHS for ease.
#     df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

#     kmg, ts, xs = kymograph(df, phase, tissue; type="frequency")

#     for (i, x) in enumerate(xs[1:end-1])

#         kmg[:, i] .-= avgω(x, phase, tissue)

#     end

#     ffigm = Plots.heatmap(ts, xs[1:end-1], kmg', color=:balance, colorbartitle="r", clims=(-0.02, 0.02), colorbar=false)
#     Plots.xlims!((0, 1000))
#     Plots.xlabel!("Time (mins)")
#     Plots.ylabel!("x (μm)")
#     Plots.title!("Dorso-Posterior")

#     # Now plot the frequency in anterior over time
#     times = sort(unique(df[!, :Time_Mins]))
#     fs = Vector{Float64}(undef, length(times))

#     for (t, time) in enumerate(times)
#         dt = df[df[!, :Time_Mins] .== time, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .>= tissue.xa, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .< tissue.xa + tissue.dc, :]
#         fs[t] = mean(dt[!, Symbol("dθ/dt")])
#     end

#     ffigm2 = Plots.plot(times, fs, label=false, color="#e37")
#     Plots.xlims!((0, 1000))
#     Plots.ylims!((0.14, 0.16))
#     Plots.ylabel!("dθ/dt")
#     Plots.xlabel!("Time (mins)")

#     # Now DP + LV

#     file = raw"data\posterior_simulations\lateral_posterior_simulations\20220324_frequency_tracking_adding_posterior_lateral_zero_cells_simulation_seed40_deterministic.csv"

#     df = DataFrame(CSV.File(file))

#     # now select out one half of the PSM. Use LHS for ease.
#     df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

#     kmg, ts, xs = kymograph(df, phase, tissue; type="frequency")

#     for (i, x) in enumerate(xs[1:end-1])

#         kmg[:, i] .-= avgω(x, phase, tissue)

#     end

#     ffigr = Plots.heatmap(ts, xs[1:end-1], kmg', color=:balance, colorbartitle="r", clims=(-0.02, 0.02), colorbar=false)
#     Plots.xlims!((0, 1000))
#     Plots.xlabel!("Time (mins)")
#     Plots.ylabel!("x (μm)")
#     Plots.title!("DP + LV")

#     # Now plot the frequency in anterior over time
#     times = sort(unique(df[!, :Time_Mins]))
#     fs = Vector{Float64}(undef, length(times))

#     for (t, time) in enumerate(times)
#         dt = df[df[!, :Time_Mins] .== time, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .>= tissue.xa, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .< tissue.xa + tissue.dc, :]
#         fs[t] = mean(dt[!, Symbol("dθ/dt")])
#     end

#     ffigr2 = Plots.plot(times, fs, label=false, color="#e37")
#     Plots.xlims!((0, 1000))
#     Plots.ylims!((0.14, 0.16))
#     Plots.ylabel!("dθ/dt")
#     Plots.xlabel!("Time (mins)")

#     # Now create a colorbar to add to the multipanel plot

#     cbartitle = L"""$\frac{d \theta}{dt} - \frac{1}{d_{c}} \int_{x_{a}}^{x_{a} + d_{c}} \omega(x) dx$"""

#     cbartitle = L"\Delta\omega"

#     cbar = Plots.scatter(
#         [0, 0], [0, 1], zcolor=[0, 0], clims=(-0.02, 0.02), xlims=(1, 1.1), label="", c=:balance, 
#         colorbartitle=cbartitle, framestyle=:none
#         )

#     l = @layout [Plots.grid(2, 3) a{0.035w}]

#     bigfig = Plots.plot(ffigl, ffigm, ffigr, ffigl2, ffigm2, ffigr2, cbar, layout=l, size=(1000, 650), leftmargin=2.5mm)

#     Plots.savefig(bigfig, "figures/20220717_comparingposition_deterministic_kymograph_frequency.svg")
# end

# # code for comparing adding posterior, DP + LV, control kymographs, for frequency
# # now without subtracting avgω, so we can see the frequency gradient
# begin

#     dc = 0.02
#     clims = (frequency(tissue.xa, phase, tissue) - dc, frequency(tissue.xa + tissue.L, phase, tissue) + dc)

#     file = raw"data\20220324_frequency_tracking_tissue_simulation_randomnewcells_seed0_deterministic.csv"

#     df = DataFrame(CSV.File(file))

#     # now select out one half of the PSM. Use LHS for ease.
#     df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

#     kmg, ts, xs = kymograph(df, phase, tissue; type="frequency")

#     ffigl = Plots.heatmap(ts, xs[1:end-1], kmg', color=:plasma, colorbartitle="r", clims=clims, colorbar=false)
#     Plots.xlims!((0, 1000))
#     Plots.xlabel!("Time (mins)")
#     Plots.ylabel!("x (μm)")
#     Plots.title!("Random")

#     # Now plot the frequency in anterior over time
#     times = sort(unique(df[!, :Time_Mins]))
#     fs = Vector{Float64}(undef, length(times))

#     for (t, time) in enumerate(times)
#         dt = df[df[!, :Time_Mins] .== time, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .>= tissue.xa, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .< tissue.xa + tissue.dc, :]
#         fs[t] = mean(dt[!, Symbol("dθ/dt")])
#     end

#     ffigl2 = Plots.plot(times, fs, label=false, color="#e37")
#     Plots.xlims!((0, 1000))
#     Plots.ylims!((0.14, 0.16))
#     Plots.ylabel!("dθ/dt")
#     Plots.xlabel!("Time (mins)")

#     # now posterior-only addition

#     file = raw"data\posterior_simulations\20220719_frequency_tracking_adding_posterior_zero_cells_simulation_seed30_deterministic.csv"

#     df = DataFrame(CSV.File(file))

#     # now select out one half of the PSM. Use LHS for ease.
#     df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

#     kmg, ts, xs = kymograph(df, phase, tissue; type="frequency")

#     ffigm = Plots.heatmap(ts, xs[1:end-1], kmg', color=:plasma, colorbartitle="r", clims=clims, colorbar=false)
#     Plots.xlims!((0, 1000))
#     Plots.xlabel!("Time (mins)")
#     Plots.ylabel!("x (μm)")
#     Plots.title!("Dorso-Posterior")

#     # Now plot the frequency in anterior over time
#     times = sort(unique(df[!, :Time_Mins]))
#     fs = Vector{Float64}(undef, length(times))

#     for (t, time) in enumerate(times)
#         dt = df[df[!, :Time_Mins] .== time, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .>= tissue.xa, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .< tissue.xa + tissue.dc, :]
#         fs[t] = mean(dt[!, Symbol("dθ/dt")])
#     end

#     ffigm2 = Plots.plot(times, fs, label=false, color="#e37")
#     Plots.xlims!((0, 1000))
#     Plots.ylims!((0.14, 0.16))
#     Plots.ylabel!("dθ/dt")
#     Plots.xlabel!("Time (mins)")

#     # Now DP + LV

#     file = raw"data\posterior_simulations\lateral_posterior_simulations\20220324_frequency_tracking_adding_posterior_lateral_zero_cells_simulation_seed40_deterministic.csv"

#     df = DataFrame(CSV.File(file))

#     # now select out one half of the PSM. Use LHS for ease.
#     df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

#     kmg, ts, xs = kymograph(df, phase, tissue; type="frequency")

#     ffigr = Plots.heatmap(ts, xs[1:end-1], kmg', color=:plasma, colorbartitle="r", clims=clims, colorbar=false)
#     Plots.xlims!((0, 1000))
#     Plots.xlabel!("Time (mins)")
#     Plots.ylabel!("x (μm)")
#     Plots.title!("DP + LV")

#     # Now plot the frequency in anterior over time
#     times = sort(unique(df[!, :Time_Mins]))
#     fs = Vector{Float64}(undef, length(times))

#     for (t, time) in enumerate(times)
#         dt = df[df[!, :Time_Mins] .== time, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .>= tissue.xa, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .< tissue.xa + tissue.dc, :]
#         fs[t] = mean(dt[!, Symbol("dθ/dt")])
#     end

#     ffigr2 = Plots.plot(times, fs, label=false, color="#e37")
#     Plots.xlims!((0, 1000))
#     Plots.ylims!((0.14, 0.16))
#     Plots.ylabel!("dθ/dt")
#     Plots.xlabel!("Time (mins)")

#     # Now create a colorbar to add to the multipanel plot

#     cbartitle = "dθ/dt"

#     cbar = Plots.scatter(
#         [0, 0], [0, 1], zcolor=[0, 0], clims=clims, xlims=(1, 1.1), label="", c=:plasma, 
#         colorbartitle=cbartitle, framestyle=:none
#         )

#     l = @layout [Plots.grid(2, 3) a{0.035w}]

#     bigfig = Plots.plot(ffigl, ffigm, ffigr, ffigl2, ffigm2, ffigr2, cbar, layout=l, size=(1000, 650), leftmargin=2.5mm)

#     Plots.savefig(bigfig, "figures/20220717_comparingposition_deterministic_kymograph_frequencygradient.svg")
# end


# # code for comparing synchrony across different densities
# begin

#     dc = 0.02
#     clims = (frequency(tissue.xa, phase, tissue) - dc, frequency(tissue.xa + tissue.L, phase, tissue) + dc)

#     file = raw"data\20220324_frequency_tracking_tissue_simulation_randomnewcells_seed0_deterministic.csv"

#     df = DataFrame(CSV.File(file))

#     # now select out one half of the PSM. Use LHS for ease.
#     df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

#     kmg, ts, xs = kymograph(df, phase, tissue; type="frequency")

#     ffigl = Plots.heatmap(ts, xs[1:end-1], kmg', color=:plasma, colorbartitle="r", clims=clims, colorbar=false)
#     Plots.xlims!((0, 1000))
#     Plots.xlabel!("Time (mins)")
#     Plots.ylabel!("x (μm)")
#     Plots.title!("Random")

#     # Now plot the frequency in anterior over time
#     times = sort(unique(df[!, :Time_Mins]))
#     fs = Vector{Float64}(undef, length(times))

#     for (t, time) in enumerate(times)
#         dt = df[df[!, :Time_Mins] .== time, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .>= tissue.xa, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .< tissue.xa + tissue.dc, :]
#         fs[t] = mean(dt[!, Symbol("dθ/dt")])
#     end

#     ffigl2 = Plots.plot(times, fs, label=false, color="#e37")
#     Plots.xlims!((0, 1000))
#     Plots.ylims!((0.14, 0.16))
#     Plots.ylabel!("dθ/dt")
#     Plots.xlabel!("Time (mins)")

#     # now posterior-only addition

#     file = raw"data\posterior_simulations\20220719_frequency_tracking_adding_posterior_zero_cells_simulation_seed30_deterministic.csv"

#     df = DataFrame(CSV.File(file))

#     # now select out one half of the PSM. Use LHS for ease.
#     df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

#     kmg, ts, xs = kymograph(df, phase, tissue; type="frequency")

#     ffigm = Plots.heatmap(ts, xs[1:end-1], kmg', color=:plasma, colorbartitle="r", clims=clims, colorbar=false)
#     Plots.xlims!((0, 1000))
#     Plots.xlabel!("Time (mins)")
#     Plots.ylabel!("x (μm)")
#     Plots.title!("Dorso-Posterior")

#     # Now plot the frequency in anterior over time
#     times = sort(unique(df[!, :Time_Mins]))
#     fs = Vector{Float64}(undef, length(times))

#     for (t, time) in enumerate(times)
#         dt = df[df[!, :Time_Mins] .== time, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .>= tissue.xa, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .< tissue.xa + tissue.dc, :]
#         fs[t] = mean(dt[!, Symbol("dθ/dt")])
#     end

#     ffigm2 = Plots.plot(times, fs, label=false, color="#e37")
#     Plots.xlims!((0, 1000))
#     Plots.ylims!((0.14, 0.16))
#     Plots.ylabel!("dθ/dt")
#     Plots.xlabel!("Time (mins)")

#     # Now DP + LV

#     file = raw"data\posterior_simulations\lateral_posterior_simulations\20220324_frequency_tracking_adding_posterior_lateral_zero_cells_simulation_seed40_deterministic.csv"

#     df = DataFrame(CSV.File(file))

#     # now select out one half of the PSM. Use LHS for ease.
#     df = df[df[!, Symbol("Position Y Reference Frame")] .<= tissue.Yc, :]

#     kmg, ts, xs = kymograph(df, phase, tissue; type="frequency")

#     ffigr = Plots.heatmap(ts, xs[1:end-1], kmg', color=:plasma, colorbartitle="r", clims=clims, colorbar=false)
#     Plots.xlims!((0, 1000))
#     Plots.xlabel!("Time (mins)")
#     Plots.ylabel!("x (μm)")
#     Plots.title!("DP + LV")

#     # Now plot the frequency in anterior over time
#     times = sort(unique(df[!, :Time_Mins]))
#     fs = Vector{Float64}(undef, length(times))

#     for (t, time) in enumerate(times)
#         dt = df[df[!, :Time_Mins] .== time, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .>= tissue.xa, :]
#         dt = dt[dt[!, Symbol("Position X Reference Frame")] .< tissue.xa + tissue.dc, :]
#         fs[t] = mean(dt[!, Symbol("dθ/dt")])
#     end

#     ffigr2 = Plots.plot(times, fs, label=false, color="#e37")
#     Plots.xlims!((0, 1000))
#     Plots.ylims!((0.14, 0.16))
#     Plots.ylabel!("dθ/dt")
#     Plots.xlabel!("Time (mins)")

#     # Now create a colorbar to add to the multipanel plot

#     cbartitle = "dθ/dt"

#     cbar = Plots.scatter(
#         [0, 0], [0, 1], zcolor=[0, 0], clims=clims, xlims=(1, 1.1), label="", c=:plasma, 
#         colorbartitle=cbartitle, framestyle=:none
#         )

#     l = @layout [Plots.grid(2, 3) a{0.035w}]

#     bigfig = Plots.plot(ffigl, ffigm, ffigr, ffigl2, ffigm2, ffigr2, cbar, layout=l, size=(1000, 650), leftmargin=2.5mm)

#     Plots.savefig(bigfig, "figures/20220717_comparingposition_deterministic_kymograph_frequencygradient.svg")
# end
