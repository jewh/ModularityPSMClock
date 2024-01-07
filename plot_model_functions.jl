# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Plot the basic functions in the model for cell movement and phase dynamics, i.e. frequency, advection velocity, rate of motility

include("experiments.jl")

using CairoMakie

# define params
begin
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
        2100,
        )
end

# set theme
publication_theme = Theme(
    fontsize=26, 
    font="CMU Serif",
    Axis=(
        xlabelsize=20, 
        ylabelsize=20,
        xgridstyle=:dash, 
        ygridstyle=:dash,
        xtickalign=1, 
        ytickalign=1, 
        yticksize=10, 
        xticksize=10,
        xlabelpadding=3,
        titlesize=30,
        titlealign=:left
        ),
    Legend=(
        framecolor=(:black, 0.5), 
        bgcolor=(:white)
        ),
    Colorbar=(
        ticksize=16, 
        tickalign=1, 
        spinewidth=0.5
        ),
    BoxPlot = Attributes(
        mediancolor=:black,
        whiskerwidth=0.3,
        strokecolor=:black,
        strokewidth=1,
        colormap=(:Egypt, 0.7), #(colormap, opacity)
        colorrange=(1, length(to_colormap(:Egypt))), # needs to start at 0 to deal with boolean variables for colour (1 otherwise)
        outliercolor=(:black, 0.5)
        ),
)

set_theme!(publication_theme)

# first plot a figure showing the model functions
fig = Figure(backgroundcolor=:grey90, resolution=(550, 550));

ax1 = Axis(
    fig[1,1],
    # aspect = AxisAspect(400 / 60),
    xlabel = L"(x - x_{a}) / L",
    ylabel = L"\omega(x) \; \text{(min^{-1})}",
    yticks = [0.14, 0.21]
    )

# xlims!(0, 1)
ylims!(0.14, 0.21)

ax2 = Axis(
    fig[2,1],
    # aspect = AxisAspect(400 / 60),
    xlabel = L"(x - x_{a}) / L",
    ylabel = L"v_{d}(x) \; \text{(μm⋅min^{-1})}",
    yticks = [0, 1.67]
    )

# xlims!(0, 1)
ylims!(0, 1.67)

ax3 = Axis(
    fig[3,1],
    # aspect = AxisAspect(400 / 60),
    xlabel = L"(x - x_{a}) / L",
    ylabel = L"v_{0}(x) \; \text{(μm⋅min^{-1})}",
    yticks = [0, 1]
    )

# xlims!(0, 1)
ylims!(0, 1)

# now plot

X = tissue.xa:(tissue.xa+tissue.L)
nx = [(x - tissue.xa) / tissue.L for x in X]

lines!(ax1, nx, [frequency(x, phase, tissue) for x in X])
lines!(ax2, nx, [-advection_velocity([x, 0, 0], tissue)[1] for x in X])
lines!(ax3, nx, [intrinsic_motion([x, 0, 0], tissue) for x in X])

save("figures/model_functions.svg", fig)




# plot a variety of motility profiles

fig = Figure(backgroundcolor=:grey90, resolution=(550, 550));

ax1 = Axis(
    fig[1,1],
    # aspect = AxisAspect(400 / 60),
    xlabel = L"(x - x_{a}) / L",
    ylabel = L"v_{0}(x)",
    yticks = [0, 1]
    )

X = tissue.xa:(tissue.xa+tissue.L)
nx = [(x - tissue.xa) / tissue.L for x in X]

for xv = 0:0.1:1

    for h = 0:0.5:5

        tissue.h = h
        tissue.Xᵥ = xv

        lines!(ax1, nx, [intrinsic_motion([x, 0, 0], tissue) for x in X])

    end

end

save("figures/motility space.svg", fig)



# plot the motility profile and how it changes
begin
    fig = Figure(resolution=(707, 500))

    tissue.h = 3.0
    tissue.Xᵥ = 0.4

    x = 0:tissue.L

    nx = [(x - tissue.xa) / tissue.L for x in x]

    # just the native function

    ax1 = Axis(
        fig[1,1],
        # aspect = AxisAspect(400 / 60),
        xlabel = L"(x - x_{a}) / L",
        ylabel = L"v_{0}(x)",
        yticks = [0, 1],
        xticks = [1]
        )

    ylims!(0, 1)
    xlims!(0, 1)

    y = [intrinsic_motion([x, 0, 0], tissue) for x in x]

    vlines!(ax1, [1 - tissue.Xᵥ], color=:black, linestyle=:dash)
    lines!(ax1, nx, y)

    text!(ax1, 0.65, 0.2; text=L"1 - X_{v}")

    fig

    ax2 = Axis(
        fig[1,2],
        # aspect = AxisAspect(400 / 60),
        xlabel = L"(x - x_{a}) / L",
        ylabel = L"v_{0}(x)",
        yticks = [0, 1],
        xticks = [1]
        )

        
    ylims!(0, 1)
    xlims!(0, 1)

    tissue.Xᵥ = 0.4

    y = [intrinsic_motion([x, 0, 0], tissue) for x in x]

    lines!(ax2, nx, y)
    vlines!(ax2, [1 - tissue.Xᵥ]; color=Cycled(1), linestyle=:dash)

    tissue.Xᵥ = 0.6

    y = [intrinsic_motion([x, 0, 0], tissue) for x in x]

    lines!(ax2, nx, y)
    vlines!(ax2, [1 - tissue.Xᵥ]; color=Cycled(2), linestyle=:dash)
    text!(ax2, 0.22, 0.8; text=L"\uparrow X_{v}\;\leftarrow", color=Makie.wong_colors()[2])

    tissue.Xᵥ = 0.2

    y = [intrinsic_motion([x, 0, 0], tissue) for x in x]

    lines!(ax2, nx, y)
    vlines!(ax2, [1 - tissue.Xᵥ]; color=Cycled(3), linestyle=:dash)
    text!(ax2, 0.8, 0.2; text=L"\rightarrow \; \downarrow X_{v}", color=Makie.wong_colors()[3])

    fig


    ax3 = Axis(
        fig[2,1],
        # aspect = AxisAspect(400 / 60),
        xlabel = L"(x - x_{a}) / L",
        ylabel = L"v_{0}(x)",
        yticks = [0, 1],
        xticks = [1]
        )

        
    ylims!(0, 1)
    xlims!(0, 1)

    tissue.Xᵥ = 0.4
    tissue.h = 3.0

    y = [intrinsic_motion([x, 0, 0], tissue) for x in x]

    lines!(ax3, nx, y)

    tissue.h = 5.0

    y = [intrinsic_motion([x, 0, 0], tissue) for x in x]

    lines!(ax3, nx, y)
    text!(ax3, 0.4, 0.5; text=L"\uparrow h", color=Makie.wong_colors()[2])

    tissue.h = 1.0

    y = [intrinsic_motion([x, 0, 0], tissue) for x in x]

    lines!(ax3, nx, y)
    text!(ax3, 0.8, 0.5; text=L"\downarrow h", color=Makie.wong_colors()[3])

    save("figures/model_intrinsicmotion.svg", fig)

    fig
end

# plot the trace of random directionality at the posterior and anterior

begin

    tissue.Dφ = 0.026

    seed = 0

    Random.seed!(seed) # 0, 6 good. 10 very good

    xstart1 = tissue.xa + 50
    xstart2 = tissue.L - 50

    x1 = SVector{3, Float64}(xstart1, 0, 0)
    x2 = SVector{3, Float64}(xstart2, 0, 0)
    n = SVector{3, Float64}(1, 0, 0)

    cell1 = FCell(x1, n, 0, 0, 0, 0)

    cell2 = FCell(x2, n, 0, 0, 0, 0)

    N = 10000

    M1 = Matrix{Float64}(undef, 2, N)
    M2 = Matrix{Float64}(undef, 2, N)

    dt = 0.01

    for t = 1:N

        M1[1, t] = x1[1] - xstart1
        M1[2, t] = x1[2]
        M2[1, t] = x2[1] - xstart2
        M2[2, t] = x2[2]

        # do in two steps so the trajectories are non-equal
        n1 = update_polarity(cell1, dt, tissue.Dφ)
        n2 = update_polarity(cell2, dt, tissue.Dφ)

        x1 += dt * intrinsic_motion(x1, tissue) * n1
        x2 += dt * intrinsic_motion(x2, tissue) * n2

        cell1 = FCell(x1, n1, 0, 0, 0, 0)
        cell2 = FCell(x2, n2, 0, 0, 0, 0)

    end

    miny = minimum([minimum(M1[2, :]), minimum(M2[2, :])]) - 1
    maxy = maximum([maximum(M1[2, :]), maximum(M2[2, :])]) + 1

    minx = minimum([minimum(M1[1, :]), minimum(M2[1, :])]) - 1
    maxx = maximum([maximum(M1[1, :]), maximum(M2[1, :])]) + 1
    

    fig = Figure()

    ax1 = Axis(fig[1,1]; aspect=DataAspect(), backgroundcolor = :transparent)

    lines!(ax1, M1[1, :], M1[2, :])
    xlims!(minx, maxx)
    ylims!(miny, maxy)

    ax2 = Axis(fig[1,2]; aspect=DataAspect(), backgroundcolor = :transparent)

    lines!(ax2, M2[1, :], M2[2, :])
    xlims!(minx, maxx)
    ylims!(miny, maxy)

    hidedecorations!(ax1)
    hidedecorations!(ax2)

    save("figures/trace_of_poster_and_anterior_motility_seed$(seed).svg", fig)

    fig

end

# plot the advection velocity profile and annotate the slopes
begin

    fig = Figure(resolution=(500, 500))
    
    ax1 = Axis(
        fig[1,1], 
        aspect=AxisAspect(1.67),
        ylabel=L"v_{d}(x) \; \text{(μm⋅min^{-1})}",
        xlabel=L"(x - x_{a})/L")

    vlines!([0.7], linestyle=:dash, color=:black)

    X = tissue.xa:(tissue.xa+tissue.L)
    nx = [(x - tissue.xa) / tissue.L for x in X]
    lines!(ax1, nx, [-advection_velocity([x, 0, 0], tissue)[1] for x in X], color=Makie.wong_colors()[1])

    text!(ax1, L"x = x_{q}", position=(0.71, 0.1), fontsize=20)
    text!(ax1, L"\frac{dv_{d}}{dx} = \frac{v_{p}(1 - x_{q}) - v_{a}}{x_{q}\cdot L}", position=(0.1, 0.6), fontsize=20)
    text!(ax1, L"\frac{dv_{d}}{dx} = -\frac{v_{p}}{L}", position=(0.75, 0.9), fontsize=20)

    save("figures/20231204_celladvectionprofile.pdf", fig)

    fig
end
