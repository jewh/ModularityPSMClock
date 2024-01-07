# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Code for replicating plots

using CSV, DataFrames, Statistics, HypothesisTests, Distributions, CairoMakie, FileIO, Glob, Dictionaries

# load in code needed for some figures
include("cell_structs.jl")
include("cell_addition.jl") 
include("density_gradient.jl")
include("dynamic_tissue_geometry.jl")
include("simulation_analysis.jl")

# helper functions for plotting the same thing repeatedly
"""
Function that associates each simulation with a color,
depending on what the type of the simulation is.
"""
function color_of_box(simulation_type::AbstractString)

    if simulation_type == "tissue_simulation" || simulation_type == "tissue_simulation_with_delay"
        
        return 1

    else

        return 2

    end

end

"""
Function that defines my legend, which is going to be common to many different plots
"""
function make_my_legend!(fig::Figure)

    # make squares of correct colour
    elem_1 = PolyElement(
        color = (to_colormap(:Egypt)[color_of_box("tissue_simulation")], 0.7), 
        strokecolor = :black,
        strokewidth = 1
        )
    elem_2 = PolyElement(
        color = (to_colormap(:Egypt)[color_of_box("adding_zero_cells")], 0.7), 
        strokecolor = :black,
        strokewidth = 1
        )

    # set the labels and title Here
    title = "Phase of new cells"
    labels = ["Random", "Constant"] # order here must match order of coloured elems

    Legend(
        fig[2, 1:2], 
        [elem_1, elem_2], 
        labels, 
        title,
        orientation = :horizontal,
        nbanks = 1,
        tellheight = true,
        valign = :top
        )
end

"""Get median and quantiles"""
function Statistics.median(V::Vector{Vector{Float64}})

    same_length = true

    for v in V

        if length(v) != length(V[1])

            same_length = false

            break

        end

    end

    @assert same_length

    medians = Vector{Float64}(undef, length(V[1]))

    for i in eachindex(V[1])

        _v = Vector{Float64}(undef, length(V))

        for (j, v) in enumerate(V)

            _v[j] = v[i]

        end

        medians[i] = median(_v)

    end

    medians

end

function Statistics.quantile(V::Vector{Vector{Float64}}, p::Float64)

    same_length = true

    for v in V

        if length(v) != length(V[1])

            same_length = false

            break

        end

    end

    @assert same_length

    quantiles = Vector{Float64}(undef, length(V[1]))

    for i in eachindex(V[1])

        _v = Vector{Float64}(undef, length(V))

        for (j, v) in enumerate(V)

            _v[j] = v[i]

        end

        if any(isnan, _v)

            quantiles[i] = NaN

        else

            quantiles[i] = quantile(_v, p)

        end

    end

    quantiles

end

function get_first_NaN(V::Vector{Float64})

    for i in eachindex(V)

        if isnan(V[i])

            return i

        end

    end

end

function get_first_NaN(V::Vector{Vector{Float64}})

    ns = Vector{Int64}(undef, length(V))

    for i in eachindex(V)

        if get_first_NaN(V[i]) === nothing

            ns[i] = length(V[i])

        else

            ns[i] = get_first_NaN(V[i])

        end

    end

    minimum(ns)

end

function Base.parse(T::Type{Tuple{PhaseParams, TissueParams}}, str::String)
    # args for PhaseParams
    bottom = maximum(findfirst("PhaseParams(", str)) + 1
    top = minimum(findfirst("),", str[bottom:end])) - 2 + bottom
    pargs = [parse(Float64, x) for x in split(str[bottom:top], ",")]
    pargs[end] == floor(pargs[end]) ? pargs[end] = parse(Int64, pargs[end]) : pargs[end] = 0

    # args for TissueParams
    bottom = maximum(findfirst("TissueParams(", str)) + 1
    top = minimum(findfirst(")", str[bottom:end])) - 2 + bottom
    targs = [parse(Float64, x) for x in split(str[bottom:top], ",")]

    # length of TissueParams should be 41 args
    if length(targs) < 41

        num_left = 41 - length(targs)
        # append zeros to the end of it. This does not affect analysis, but is needed to prevent errors.
        append!(targs, [0.0 for _ = 1:num_left])

    end

    (PhaseParams(pargs...), TissueParams(targs...))
end

function Base.parse(T::Type{Vector{Float64}}, s::String)
    s = replace(s, "[" => "")
    s = replace(s, "]" => "")
    s = replace(s, " " => "")
    s = eachsplit(s, ",")
    parse.(Float64, s) 
end

# ugly fix
Base.parse(T::Type{Int64}, x::Float64) = Int64(floor(x))

function parse_line(s::String)

    if occursin("BoundsError", str)

        phase_nan = PhaseParams(NaN, NaN, NaN, NaN, NaN, NaN)

        tissue_nan = TissueParams(
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            0.,
            0.,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN
        )

        return NaN, [NaN], [NaN], (phase_nan, tissue_nan)

    else

        bottom = maximum(findfirst("(", str)) + 1

        vals = split(str[1:bottom-2], "\t")

        Ns = parse(Float64, vals[1])

        Z = parse(Vector{Float64}, string(vals[2]))

        freqs = parse(Vector{Float64}, string(vals[3]))

        phase, tissue = parse(Tuple{PhaseParams, TissueParams}, str[bottom:end-1])

        return Ns, Z, freqs, (phase, tissue)

    end

end

# ugly but handles the case where there's no parsing needed
parse_column(T::Type{Float64}, Vs::Vector{Float64}) = Vs

function parse_column(T::Type{Float64}, Vs::Vector{String})

    out = Vector{Float64}(undef, length(Vs))

    for i in eachindex(Vs)

        if occursin("Error", Vs[i])

            out[i] = NaN

        else

            out[i] = parse(Float64, Vs[i])

        end

    end

    out

end

function parse_column(T::Type{Vector{Float64}}, Vs::Vector{String})

    out = Vector{Vector{Float64}}(undef, length(Vs))

    for i in eachindex(Vs)

        if occursin("Error", Vs[i])

            out[i] = [NaN]

        else

            out[i] = parse(Vector{Float64}, Vs[i])

        end

    end

    out

end

function parse_column(T::Type{Tuple{PhaseParams, TissueParams}}, Vs::Vector{String})

    out = Vector{Tuple{PhaseParams, TissueParams}}(undef, length(Vs))

    for i in eachindex(Vs)

        if occursin("Error", Vs[i])

            phase_nan = PhaseParams(NaN, NaN, NaN, NaN, NaN, 0)

            tissue_nan = TissueParams(
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    0.,
                    0.,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN
                )

            out[i] = phase_nan, tissue_nan

        else

            out[i] = parse(Tuple{PhaseParams, TissueParams}, Vs[i])

        end

    end

    out

end

function file_to_vect(path; header=true)
    V = []
    open(path, "r") do io
        for (i, line) in enumerate(readlines(io))
            if header && (i == 1 || i == 2)
                continue
            end
            push!(V, parse_line(line))
        end
    end
    V
end

function Base.isnan(v::Vector{Float64})

    for i in eachindex(v)

        if isnan(v[i])

            return true

        end

    end

    false

end

Base.isnan(P::PhaseParams) = P.k == NaN
    
Base.isnan(T::TissueParams) = T.dc == NaN

Base.isnan(x::Tuple{PhaseParams, TissueParams}) = isnan(x[1]) || isnan(x[2])

function file_to_df(path; drop_nans=true)

    @assert path[end-3:end] == ".tsv"

    data = DataFrame(CSV.File(path, delim="\t", header=2, skipto=3))

    data[!, :Ns] = parse_column(Float64, data[!, :Ns])

    data[!, :Zs] = parse_column(Vector{Float64}, data[!, :Zs])

    data[!, :freqs] = parse_column(Vector{Float64}, data[!, :freqs])

    data[!, Symbol("(phase, tissue)")] = parse_column(Tuple{PhaseParams, TissueParams}, Vector{String}(data[!, Symbol("(phase, tissue)")]))

    if drop_nans

        data = filter(row->!any(isnan(x) for x in row), data)

    end

    data

end

# allow saving of .pdf plots
# CairoMakie.activate!(type = "pdf")

# define the theme constant to all plots
publication_theme = Theme(
    fontsize=20, 
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

# load in my nice data
path = raw"data\20221208_stochastic_simulations.tsv"
# path = raw"data\20230125_stochastic_delay_simulations.tsv"

data = DataFrame(CSV.File(path, delim="\t"))

data = data[completecases(data), :]

# plot the first figure - comparing frequency and synchrony for different morphogenetic conditions
begin

    addition_position_fig = Figure(resolution=(750, 550), backgroundcolor=:grey90)

    synch_ax = Axis(
        addition_position_fig[1, 1], 
        aspect = AxisAspect(1), 
        xlabel = "Position of new cells", 
        ylabel = L"r",
        xticks = (
            [1, 2, 3], ["Random", "DP", "DP + LV"]),
        title = "A",
        )

    ylims!(synch_ax, 0, 1)

    function filter_2(x)
        x in ["tissue_simulation", "adding_posterior_zero_cells_simulation", "adding_posterior_lateral_zero_cells_simulation"]
    end

    function position_of_x(simulation_type::AbstractString)

        if simulation_type == "tissue_simulation"

            return 1

        elseif simulation_type == "adding_posterior_zero_cells_simulation"

            return 2

        elseif simulation_type == "adding_posterior_lateral_zero_cells_simulation"

            return 3

        end

    end

    df = filter(Symbol("Simulation Type") => filter_2, data)

    CairoMakie.boxplot!(
        synch_ax, 
        position_of_x.(df[!, Symbol("Simulation Type")]), 
        df[!, Symbol("Synchrony (r)")]; 
        color = color_of_box.(df[!, Symbol("Simulation Type")]),  # this chooses which colour pair of the scheme to use.
        )

    freq_ax = Axis(
        addition_position_fig[1, 2], 
        aspect = AxisAspect(1), 
        xlabel = "Position of new cells", 
        ylabel = L"\frac{d\theta}{dt} \; \text{(min^{-1})}",
        xticks = (
            [1, 2, 3], ["Random", "DP", "DP + LV"]),
        title = "B",
        )

    CairoMakie.boxplot!(
        freq_ax, 
        position_of_x.(df[!, Symbol("Simulation Type")]), 
        df[!, Symbol("Mean frequency (/min)")]; 
        color = color_of_box.(df[!, Symbol("Simulation Type")]),  # this chooses which colour pair of the scheme to use.
        )

    make_my_legend!(addition_position_fig)

    rowsize!(addition_position_fig.layout, 2, Relative(0.2)) # to move legend up

    trim!(addition_position_fig.layout)

    # resize_to_layout!(addition_position_fig)

    save("figures/20221208_stochastic_simulations_cell_addition_fig.pdf", addition_position_fig)

    addition_position_fig

end

# plot the second figure comparing different effects of mitosis
begin

    mitosis_fig1 = Figure(resolution=(750, 750), backgroundcolor=:grey90)

    ax1 = Axis(
        mitosis_fig1[1, 1], 
        aspect = AxisAspect(1), 
        xlabel = "Dividing Cells", 
        ylabel = L"r",
        xticks = (
            [1, 2], ["Absent", "Present"]),
        title = "A",
        )

    ylims!(ax1, 0, 1)

    filter_2(x) = x == "mitosing_simulation" || x == "tissue_simulation"

    function position_of_x(simulation_type::AbstractString)

        if simulation_type == "tissue_simulation"

            return 1

        elseif simulation_type == "mitosing_simulation"

            return 2

        end

    end

    # control colour from here
    function colour(simulation_type::AbstractString)

        if simulation_type == "tissue_simulation"

            return Makie.wong_colors()[1]

        else

            return Makie.wong_colors()[3]

        end

    end

    df = filter(Symbol("Simulation Type") => filter_2, data)

    CairoMakie.boxplot!(
        ax1, 
        position_of_x.(df[!, Symbol("Simulation Type")]), 
        df[!, Symbol("Synchrony (r)")]; 
        color = colour.(df[!, Symbol("Simulation Type")])
        )

    ax2 = Axis(
        mitosis_fig1[1, 2], 
        aspect = AxisAspect(1), 
        xlabel = "Dividing Cells", 
        ylabel = L"\frac{d\theta}{dt} \; \text{(min^{-1})}",
        xticks = ([1, 2], ["Absent", "Present"]),
        title = "B",
        )

    CairoMakie.boxplot!(
        ax2, 
        position_of_x.(df[!, Symbol("Simulation Type")]), 
        df[!, Symbol("Mean frequency (/min)")]; 
        color = colour.(df[!, Symbol("Simulation Type")])
        )

    filter_3(x) = occursin("mitosing_simulation_x_div", x)

    x_div(x::AbstractString) = parse(Float64, x[26:end])

    df = filter(Symbol("Simulation Type") => filter_3, data)

    x = unique(x_div.(df[!, Symbol("Simulation Type")]))
    
    ax3 = Axis(
        mitosis_fig1[2, 1], 
        aspect = AxisAspect(1), 
        xlabel = L"x_{div} \; (\text{\mu\!\>m})",
        ylabel = L"r",
        xticks = (x, string.(x)),
        title = "C",
        )

    xlims!(ax3, -35, 385) # -35 so equal padding either side of figure
    ylims!(ax3, 0, 1)

    CairoMakie.boxplot!(
        ax3, 
        x_div.(df[!, Symbol("Simulation Type")]), 
        df[!, Symbol("Synchrony (r)")]; 
        color = colour.(df[!, Symbol("Simulation Type")]),  # this chooses which colour pair of the scheme to use.
        width=75.
        )

    ax4 = Axis(
        mitosis_fig1[2, 2], 
        aspect = AxisAspect(1), 
        xlabel = L"x_{div} \; \text{(\mu\!\>m)}",
        ylabel = L"\frac{d\theta}{dt} \; \text{(min^{-1})}",
        xticks = (x, string.(x)),
        title = "D",
        )

    xlims!(ax4, -35, 385) # -35 so equal padding either side of figure

    CairoMakie.boxplot!(
        ax4, 
        x_div.(df[!, Symbol("Simulation Type")]), 
        df[!, Symbol("Mean frequency (/min)")]; 
        color = colour.(df[!, Symbol("Simulation Type")]), # this chooses which colour pair of the scheme to use.
        width=75.
        )

    # now plot the long figure showing effect on frequency and synchrony
    mitosis_fig2 = Figure(resolution=(848, 600), backgroundcolor=:grey90)

    filter_4(x) = occursin("mitosing_simulation_TM", x)

    # round because the floating point error makes the numbers very long and unreadable
    get_TM(simulation_type::AbstractString) = round(parse(Float64, simulation_type[23:end]); digits=4)

    df = filter(Symbol("Simulation Type") => filter_4, data)

    x = unique(get_TM.(df[!, Symbol("Simulation Type")])) 

    ax5 = Axis(
        mitosis_fig2[1, 1:2], 
        # aspect = AxisAspect(1), 
        xlabel = L"T_{M} \; \text{(mins)}",
        ylabel = L"r",
        xticks = ([x for x in 0:10:50], string.([x for x in 0:10:50])),
        title = "E",
        )

    CairoMakie.boxplot!(
        ax5, 
        get_TM.(df[!, Symbol("Simulation Type")]), 
        df[!, Symbol("Synchrony (r)")]; 
        color = colour.(df[!, Symbol("Simulation Type")])
        )

    ax6 = Axis(
        mitosis_fig2[2, 1:2], 
        # aspect = AxisAspect(1), 
        xlabel = L"T_{M} \; \text{(mins)}",
        ylabel = L"\frac{d\theta}{dt} \; \text{(min^{-1})}",
        xticks = ([x for x in 0:10:50], string.([x for x in 0:10:50])),
        # yticks = 
        title = "F",
        )

    CairoMakie.boxplot!(
        ax6, 
        get_TM.(df[!, Symbol("Simulation Type")]), 
        df[!, Symbol("Mean frequency (/min)")]; 
        color = colour.(df[!, Symbol("Simulation Type")])
        )

    save("figures/20221208_stochastic_simulations_mitosis_fig1.pdf", mitosis_fig1)
    save("figures/20221208_stochastic_simulations_mitosis_fig2.pdf", mitosis_fig2)

    mitosis_fig1

end

# plot the effect of changing the density of the tissue
begin

    colors = [Makie.wong_colors()[1], Makie.wong_colors()[3], Makie.wong_colors()[4]]

    fig = Figure(resolution=(750, 400), backgroundcolor=:grey90)

    ax1 = Axis(
        fig[1, 1], 
        aspect = AxisAspect(1), 
        xlabel = L"ρ_{0} \; \text{(μm^{-3})}", 
        ylabel = L"r",
        xticks = ([1, 2, 3], ["0.00075", "0.0015", "0.003"]),
        title = "A",
        )

    ylims!(ax1, 0, 1)

    ax2 = Axis(
        fig[1, 2], 
        aspect = AxisAspect(1), 
        xlabel = L"ρ_{0} \; \text{(μm^{-3})}", 
        ylabel = L"\frac{d\theta}{dt} \; \text{(min^{-1})}",
        xticks = ([1, 2, 3], ["0.00075", "0.0015", "0.003"]),
        title = "B",
        )

    filter_rho_experiment(x) = occursin("_changing_density_rho", x) || x in ["tissue_simulation", "adding_posterior_zero_cells_simulation",  "adding_posterior_lateral_zero_cells_simulation"]

    function get_rho(x)
        if findfirst("_rho", x) === nothing
            return 0.0015
        else
            start = findfirst("_rho", x)[end] + 1
            return parse(Float64, x[start:end])
        end
    end

    df = filter(Symbol("Simulation Type") => filter_rho_experiment, data)

    random_addition(x) = occursin("tissue_simulation", x)
    dp_addition(x) = occursin("adding_posterior_zero_cells", x)
    dplv_addition(x) = occursin("adding_posterior_lateral_zero_cells", x) 

    dodge(x) = random_addition(x) ? 1 : dp_addition(x) ? 2 : dplv_addition(x) ? 3 : nothing
    x_position(x) = get_rho(x) == 0.00075 ? 1 : get_rho(x) == 0.0015 ? 2 : get_rho(x) == 0.003 ? 3 : 4 

    x = x_position.(df[!, Symbol("Simulation Type")])

    dodges = dodge.(df[!, Symbol("Simulation Type")])

    CairoMakie.boxplot!(
        ax1,
        x, 
        df[!, Symbol("Synchrony (r)")], 
        dodge = dodges, 
        color = [colors[i] for i in dodges]
    )

    CairoMakie.boxplot!(
        ax2,
        x, 
        df[!, Symbol("Mean frequency (/min)")], 
        dodge = dodges, 
        color = [colors[i] for i in dodges]
    )

    save("figures/20221208_stochastic_simulations_changing_density_fig.pdf", fig)

    fig 

end

# now figure for effect of cell density gradient

begin

    density_figure1 = Figure(
        resolution = (800, 400),
        backgroundcolor=:grey90
        )

    # plot density function and measured density

    tissue = TissueParams(11.0, 8.71, 1.67, 3.0, 0.7, 385.0, 25.0, 0.4, 1.0, 3.0, 0.026, 20.0, 1.0, 325.0, 110.0, 50.0, 25.0, 60.0, 0.0015, 2321, 2321, 172.5, 15.0, 100.0, 0.75, 1.1, 325.0, 0.0, 0.002123, 0.0, 27.6, 0.0, 253.6, 410.0, 195.0, 75.0, 0.0, 9.2, 0.0, 11.0, 536.4)

    ax1 =  Axis(
        density_figure1[1, 1], 
        aspect = AxisAspect(1), 
        xlabel = L"x",
        ylabel = L"d(x) \; \text{(\mu\!\>m)}",
        xticks = ([tissue.xa, tissue.Xc, tissue.xa + tissue.L], [L"x_{a}", L"X_{c}", L"x_{a} + L"]),
        yticks = ([11, 13, 15], ["11", "13", "15"]),
        title = "A",
        )

    x = tissue.xa:tissue.xa+tissue.L

    for dₘ in (0.65, 1.1, 2.2, 4.4)

        tissue.dm = dₘ

        y = [d(t, tissue) for t in x]

        CairoMakie.lines!(
            ax1,
            x,
            y,
            label = "$(tissue.dm)"
        )

    end

    xlims!(ax1, 0, tissue.xa + tissue.L)
    ylims!(ax1, 11, 16.5)

    # axislegend(L"d_{m}" ; position = :lt)

    density_figure1

    function df_to_vectFCell(df)

        X = Vector{FCell}(undef, length(unique(df[!, :TrackID])))
    
        for (i, tid) in enumerate(unique(df[!, :TrackID]))
    
            dtid = df[df[!, :TrackID] .== tid, :]
    
            x = SVector{3, Float64}(
                dtid[!, Symbol("Position X Reference Frame")][1],
                dtid[!, Symbol("Position Y Reference Frame")][1],
                dtid[!, Symbol("Position Z Reference Frame")][1],
            )
    
            n = SVector{3, Float64}([0., 0., 0.])
    
            X[i] = FCell(
                x, n, 0., 0., 0., tid
            )
    
        end
    
        X
    
    end

    """Get density for the left hand side of the tissue."""
    function left_hand_tissue_density(file::String)

        tissue = TissueParams(11.0, 8.71, 1.67, 3.0, 0.7, 385.0, 25.0, 0.4, 1.0, 3.0, 0.026, 20.0, 1.0, 325.0, 110.0, 50.0, 25.0, 60.0, 0.0015, 2321, 2321, 172.5, 15.0, 100.0, 0.75, 1.1, 325.0, 0.0, 0.002123, 0.0, 27.6, 0.0, 253.6, 410.0, 195.0, 75.0, 0.0, 9.2, 0.0, 11.0, 536.4)
    
        df = DataFrame(CSV.File(file))

        # subset the data to isolate the PSM
        df = df[df[!, :Time_Mins] .== 999., :]
        df = df[tissue.xa .<= df[!, Symbol("Position Y Reference Frame")] .<= tissue.Xc, :]
        df = df[df[!, Symbol("Position Y Reference Frame")] .< tissue.Yc, :]

        ldensity = Vector{Float64}(undef, length(tissue.xa+tissue.dc:tissue.dc:tissue.Xc))

        for (i, x) in enumerate(tissue.xa+tissue.dc:tissue.dc:tissue.Xc)

            n = 0

            for row in eachrow(df)

                if x - tissue.dc <= row[Symbol("Position X Reference Frame")] < x

                    n += 1

                end

            end

            ldensity[i] = n

        end

        ldensity / ( π * tissue.r^2 * tissue.dc)

    end

    ax2 = Axis(
        density_figure1[1,2], 
        aspect = AxisAspect(1), 
        xlabel = L"x\; \text{(\mu\!\>m)}",
        ylabel = L"\rho \; \text{(\mu\!\>m^{-3})}",
        xticks = ([x for x in 0:100:1000], string.([x for x in 0:100:1000])),
        title = "B",
        )
    
    xlims!(ax2, 0, 325)
    ylims!(ax2, 0.0005, 0.0025)

    x = [x for x = tissue.xa:tissue.dc:tissue.Xc-tissue.dc]

    # define which files to extract density from

    for val in [0.65, 1.1, 2.2, 4.4]

        files = Glob.glob("data/20221201_frequency_tracking_*_volume_simulation_seed*_dm$(val).csv")

        v = Vector{Vector{Float64}}(undef, length(files))

        for (i, f) in enumerate(files)

            v[i] = left_hand_tissue_density(f)

        end

        lower_q = quantile(v, 0.25)
        upper_q = quantile(v, 0.75)

        band!(ax2, x, lower_q, upper_q)

        CairoMakie.lines!(ax2, x, median(v))

    end

    Legend(density_figure1[1, 3], ax1, L"d_{m}", orientation = :vertical, tellwidth = true, tellheight = false)

    # define helper functions

    density_figure2 = Figure(
        resolution = (600, 600),
        backgroundcolor=:grey90
        )

    filter_rand(x) = occursin("increasing_extracellular_volume_simulation_dm", x)

    filter_dp(x) = occursin("increasing_extracellular_volume_adding_posterior_zero_cells_simulation_dm", x)

    filter_dplv(x) = occursin("increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation_dm", x)

    filter_5(x) = filter_rand(x) || filter_dp(x) || filter_dplv(x)

    dodge(sim_type::AbstractString) = filter_rand(sim_type) ? 1 : filter_dp(sim_type) ? 2 : filter_dplv(sim_type) ? 3 : nothing

    function get_dm(x::AbstractString)
        start = findfirst("_dm", x)[end] + 1
        parse(Float64, x[start:end])
    end

    x_position(x::AbstractString) = get_dm(x) == 0.65 ? 1 : get_dm(x) == 1.1 ? 2 : get_dm(x) == 2.2 ? 3 : 4 

    ax3 = Axis(
        density_figure2[1,1:2], 
        # aspect = AxisAspect(1), 
        xlabel = L"d_{m} \; \text{(\mu\!\>m)}",
        ylabel = L"r",
        xticks = ([1, 2, 3, 4], string.([0.65, 1.1, 2.2, 4.4])),
        title = "C",
        )

    ax4 = Axis(
        density_figure2[2,1:2], 
        # aspect = AxisAspect(1), 
        xlabel = L"d_{m} \; \text{(\mu\!\>m)}",
        ylabel = L"\frac{d\theta}{dt}  \; \text{(min^{-1})}",
        xticks = ([1, 2, 3, 4], string.([0.65, 1.1, 2.2, 4.4])),
        title = "D",
        )

    df = filter(Symbol("Simulation Type") => filter_5, data)

    x = x_position.((df[!, Symbol("Simulation Type")]))
    dodges = dodge.(df[!, Symbol("Simulation Type")])

    CairoMakie.boxplot!(
        ax3,
        x, 
        df[!, Symbol("Synchrony (r)")], 
        dodge = dodges, 
        color = dodges,
        colormap = (:Kandinsky, 0.7)
    )

    ylims!(ax3, 0, 1)

    CairoMakie.boxplot!(
        ax4,
        x, 
        df[!, Symbol("Mean frequency (/min)")], 
        dodge = dodges, 
        color = dodges,
        colormap = (:Kandinsky, 0.7)
    )

    # make a legend for the last two plots

    elem_1 = PolyElement(
        color = (to_colormap(:Kandinsky)[dodge("increasing_extracellular_volume_simulation_dm")], 0.7), 
        strokecolor = :black,
        strokewidth = 1
        )
    elem_2 = PolyElement(
        color = (to_colormap(:Kandinsky)[dodge("increasing_extracellular_volume_adding_posterior_zero_cells_simulation_dm")], 0.7), 
        strokecolor = :black,
        strokewidth = 1
        )
    elem_3 = PolyElement(
        color = (to_colormap(:Kandinsky)[dodge("increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation_dm")], 0.7), 
        strokecolor = :black,
        strokewidth = 1
        )

    # set the labels and title Here
    title = "Cell Ingression"
    labels = ["Random", "DP", "DP + LV"] # order here must match order of coloured elems

    Legend(
        density_figure2[1:2, 3], 
        [elem_1, elem_2, elem_3], 
        labels, 
        title,
        orientation = :vertical,
        # nbanks = 1,
        tellwidth = true,
        tellheight = false,
        )

    save("figures/20221208_stochastic_simulations_density_gradient_fig1.pdf", density_figure1)
    save("figures/20221208_stochastic_simulations_density_gradient_fig2.pdf", density_figure2)

    density_figure2

end

# now do figure for convergence-extension

begin

    include("dynamic_tissue_geometry.jl")

    shrinking_tissue_figure = Figure(
        resolution = (1000, 600),
        backgroundcolor=:grey90
        )

    ds = DataFrame(CSV.File(raw"schroteretal2008dataset.csv"))

    data_ax = Axis(
        shrinking_tissue_figure[1, 1], 
        # aspect = AxisAspect(1), 
        xlabel = L"t \text{(mins)}", 
        ylabel = "somites",
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        # title = "A",
        )

    ylims!(data_ax, 0, 34)
    
    scatter!(data_ax, ds[!, Symbol("Time (min)")], ds[!, :Somite], color=(:deeppink3, 0.7))

    y = [somite(t) for t = -50:800]

    lines!(data_ax, -50:800, y, color=(:black, 1.0))

    eqn = L"s(t) = 6 + \frac{t}{24.7} - 0.5001e^{0.0049(t - 300)}"

    text!(data_ax, eqn, position=(-2, 30), space = :data)

    # now make a sub grid layout with the three model functions

    fa = shrinking_tissue_figure[1, 2] = GridLayout(title="B")

    fax1 = Axis(
        fa[1,1],
        xlabel=L"t \; \text{(mins)}",
        ylabel=L"L(t)  \; \text{(\mu\!\>m)}", 
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        yticks = (
            [x for x = 0:100:350], string.([x for x = 0:100:350])
        )
        )

    ylims!(fax1, 0, 350)

    fax2 = Axis(
        fa[2,1],
        xlabel=L"t\; \text{(mins)}",
        ylabel=L"r(t)\; \text{(\mu\!\>m)}", 
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        yticks = (
            [x for x = 0:15:45], string.([x for x = 0:15:45])
        )       
        )

    ylims!(fax2, 0, 30)

    fax3 = Axis(
        fa[3,1],
        xlabel=L"t\; \text{(mins)}",
        ylabel=L"\rho(t)\; \text{(\mu\!\>m^{-3})}", 
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        yticks = (
            [x for x = 0:0.002:0.004], string.([x for x = 0:0.002:0.004])
            )
        )
    
    ylims!(fax3, 0, 0.005)
       
    ts = 0:772

    t_shrink = 253.6

    Lt = [Length(t, 14.1, 325, t_shrink) for t in ts]

    rt = [radius(t, 1.1625, 27.6, t_shrink) for t in ts]

    dt = [density(t, 0.000121, 0.002123, t_shrink) for t in ts]

    lines!(fax1, ts, Lt)
    lines!(fax2, ts, rt)
    lines!(fax3, ts, dt)

    vlines!(fax1, [t_shrink]; color=:black, linestyle=:dash)
    vlines!(fax2, [t_shrink]; color=:black, linestyle=:dash)
    vlines!(fax3, [t_shrink]; color=:black, linestyle=:dash)

    save("figures/shrinking_tissue_figure.pdf", shrinking_tissue_figure)

    df = DataFrame(CSV.File(raw"data\20230224_stochastic_shrinking.tsv", delim="\t"))

    function Base.parse(T::Type{Vector{Float64}}, s::String)
        s = replace(s, "[" => "")
        s = replace(s, "]" => "")
        s = replace(s, " " => "")
        s = eachsplit(s, ",")
        parse.(Float64, s) 
    end

    df[!, Symbol("Synchrony (r)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Synchrony (r)")]]

    df[!, Symbol("Mean frequency (/min)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Mean frequency (/min)")]]

    df[!, Symbol("Density (μm⁻³)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Density (μm⁻³)")]]

    df[!, Symbol("No. Cells")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("No. Cells")]]

    df[!, Symbol("xa (μm)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("xa (μm)")]]

    df[!, Symbol("Cell Additions")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Cell Additions")]]

    filter_6(x) = "tissue_simulation_shrinking_control" == x
    filter_7(x) = "tissue_simulation_shrinking_constantavdection" == x # typo in the file
    filter_8(x) = "tissue_simulation_shrinking" == x

    shrinking_figure = Figure(
        resolution = (1200, 500),
        backgroundcolor=:grey90
        )

    ax1 = Axis(
        shrinking_figure[1, 1], 
        # ax[1, 1],
        # aspect = AxisAspect(1), 
        xlabel = L"t \; \text{(mins)}", 
        ylabel = L"r",
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        title = "A",
        )

    ax2 = Axis(
        shrinking_figure[1, 2], 
        # ax[1, 2],
        # aspect = AxisAspect(1), 
        xlabel = L"t \; \text{(mins)}", 
        ylabel = L"\frac{d\theta}{dt} \; \text{(min^{-1})}",
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        title = "B",
        )

    vlines!(ax1, [t_shrink]; color=:black, linestyle=:dash)
    vlines!(ax2, [t_shrink]; color=:black, linestyle=:dash)

    # control colour here
    function color(type::AbstractString)

        if type == "tissue_simulation_shrinking_control"

            return Makie.wong_colors()[1]

        elseif type == "tissue_simulation_shrinking_constantavdection"

            return Makie.wong_colors()[2]

        end

    end

    nstop = 656 # beyond this index shrinking sim fails

    for f in [filter_6, filter_7]

        dd = filter(Symbol("Simulation Type") => f, df)

        type = dd[1, Symbol("Simulation Type")]

        nt = get_first_NaN(dd[!, Symbol("Synchrony (r)")]) - 1

        nt = nstop

        for row in eachrow(dd)

            row[Symbol("Synchrony (r)")] = row[Symbol("Synchrony (r)")][1:nt]
            row[Symbol("Mean frequency (/min)")] = row[Symbol(Symbol("Mean frequency (/min)"))][1:nt]

        end

        x = 0:nt-1

        med = median(dd[!, Symbol("Synchrony (r)")])
        lower_q = quantile(dd[!, Symbol("Synchrony (r)")], 0.25) 
        upper_q = quantile(dd[!, Symbol("Synchrony (r)")], 0.75) 


        CairoMakie.lines!(
            ax1,
            x,
            med,
            color=(color(type), 1.0)
        )

        band!(ax1, x, lower_q, upper_q, color=(color(type), 0.5))

        ylims!(ax1, 0, 1)

        med = median(dd[!, Symbol("Mean frequency (/min)")])
        lower_q = quantile(dd[!, Symbol("Mean frequency (/min)")], 0.25) 
        upper_q = quantile(dd[!, Symbol("Mean frequency (/min)")], 0.75) 

        CairoMakie.lines!(
            ax2,
            x,
            med,
            color=(color(type), 1.0)
        )

        band!(ax2, x, lower_q, upper_q, color=(color(type), 0.5))
    
    end

    xlims!(ax1, 0, 700)
    xlims!(ax2, 0, 700)

    # annotate
    text!(ax1, L"t = t_{shrink}", position=(270, 0.05))
    text!(ax2, L"t = t_{shrink}", position=(270, 0.12))

    # make legend
    elem_1 = PolyElement(
        color = (color("tissue_simulation_shrinking_control"), 0.7), 
        strokecolor = :black,
        strokewidth = 1
        )
    elem_2 = PolyElement(
        color = (color("tissue_simulation_shrinking_constantavdection"), 0.7), 
        strokecolor = :black,
        strokewidth = 1
        )

    # set the labels and title Here
    title = "Legend"
    labels = ["Constant PSM", "Compacting PSM"] # order here must match order of coloured elems

    Legend(
        shrinking_figure[2, 1:2], 
        [elem_1, elem_2], 
        labels, 
        title,
        orientation = :horizontal,
        nbanks = 1,
        tellheight = true,
        valign = :center
        )

    save("figures/20230224_shrinking.pdf", shrinking_figure)

    shrinking_figure

end


# plot across multiple time points to remove the NaNs - this is fine

begin
    # figure for shrinking tissue with different densities
    df = DataFrame(CSV.File(raw"data\20230331_stochastic_shrinking.tsv", delim="\t"))

    function Base.parse(T::Type{Vector{Float64}}, s::String)
        s = replace(s, "[" => "")
        s = replace(s, "]" => "")
        s = replace(s, " " => "")
        s = eachsplit(s, ",")
        parse.(Float64, s) 
    end

    df[!, Symbol("Synchrony (r)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Synchrony (r)")]]

    df[!, Symbol("Mean frequency (/min)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Mean frequency (/min)")]]

    df[!, Symbol("Density (μm⁻³)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Density (μm⁻³)")]]

    df[!, Symbol("No. Cells")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("No. Cells")]]

    df[!, Symbol("xa (μm)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("xa (μm)")]]

    df[!, Symbol("Cell Additions")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Cell Additions")]]

    filter_7(x) = "tissue_simulation_shrinking_constantavdection" == x
    filter_9(x) = "tissue_simulation_shrinking_constantavdection_d00.001123" == x # typo in the file
    filter_10(x) = "tissue_simulation_shrinking_constantavdection_d00.003123" == x

    function color(type::AbstractString)

        if type == "tissue_simulation_shrinking_constantavdection"

            return Makie.wong_colors()[2]

        elseif type == "tissue_simulation_shrinking_constantavdection_d00.001123"

            return Makie.wong_colors()[3]

        elseif type == "tissue_simulation_shrinking_constantavdection_d00.003123"

            return Makie.wong_colors()[4]

        end

    end

    fig = Figure(
        resolution = (800, 800),
        backgroundcolor=:grey90
        )

    # plot the different density profiles
    ax1 = Axis(
        fig[1, 2],
        xlabel = L"t \; \text{(mins)}", 
        ylabel = L"ρ \; \text{(μm⁻³)}",
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        title = "A",
        )

    md = 0.000121
    t_shrink = 253.6
    t_max = 690.0

    ctrl_density = [density(t, md, 0.002123, t_shrink) for t = 0:t_max]
    low_density = [density(t, md, 0.001123, t_shrink) for t = 0:t_max]
    high_density = [density(t, md, 0.003123, t_shrink) for t = 0:t_max]

    lines!(ax1, 0:t_max, ctrl_density, color = color("tissue_simulation_shrinking_constantavdection"))
    lines!(ax1, 0:t_max, low_density, color = color("tissue_simulation_shrinking_constantavdection_d00.001123"))
    lines!(ax1, 0:t_max, high_density, color = color("tissue_simulation_shrinking_constantavdection_d00.003123"))

    vlines!(ax1, [t_shrink]; color=:black, linestyle=:dash)

    xlims!(ax1, 0, 700)
    ylims!(ax1, 0, 0.005)

    ax2 = Axis(
        fig[2, 1], 
        # ax[1, 1],
        # aspect = AxisAspect(1), 
        xlabel = L"t \; \text{(mins)}", 
        ylabel = L"r",
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        title = "B",
        )

    ax3 = Axis(
        fig[2, 2], 
        # ax[1, 2],
        # aspect = AxisAspect(1), 
        xlabel = L"t \; \text{(mins)}", 
        ylabel = L"\frac{d\theta}{dt} \; \text{(min^{-1})}",
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        title = "C",
        )

    vlines!(ax2, [t_shrink]; color=:black, linestyle=:dash)
    vlines!(ax3, [t_shrink]; color=:black, linestyle=:dash)

    for f in [filter_7, filter_9, filter_10]

        dd = filter(Symbol("Simulation Type") => f, df)

        type = dd[1, Symbol("Simulation Type")]

        nt = get_first_NaN(dd[!, Symbol("Synchrony (r)")]) - 1

        # nt = 690

        for row in eachrow(dd)

            row[Symbol("Synchrony (r)")] = row[Symbol("Synchrony (r)")][1:nt]
            row[Symbol("Mean frequency (/min)")] = row[Symbol(Symbol("Mean frequency (/min)"))][1:nt]

        end

        x = 0:nt-1

        med = median(dd[!, Symbol("Synchrony (r)")])
        lower_q = quantile(dd[!, Symbol("Synchrony (r)")], 0.25) 
        upper_q = quantile(dd[!, Symbol("Synchrony (r)")], 0.75) 


        CairoMakie.lines!(
            ax2,
            x,
            med,
            color=(color(type), 1.0)
        )

        band!(ax2, x, lower_q, upper_q, color=(color(type), 0.5))

        ylims!(ax2, 0, 1)

        med = median(dd[!, Symbol("Mean frequency (/min)")])
        lower_q = quantile(dd[!, Symbol("Mean frequency (/min)")], 0.25) 
        upper_q = quantile(dd[!, Symbol("Mean frequency (/min)")], 0.75) 

        CairoMakie.lines!(
            ax3,
            x,
            med,
            color=(color(type), 1.0)
        )

        band!(ax3, x, lower_q, upper_q, color=(color(type), 0.5))
    
    end

    xlims!(ax2, 0, 700)
    xlims!(ax3, 0, 700)

    # annotate
    text!(ax1, L"t = t_{shrink}", position=(270, 0.0005))
    text!(ax2, L"t = t_{shrink}", position=(270, 0.05))
    text!(ax3, L"t = t_{shrink}", position=(270, 0.12))

    # make legend
    elem_1 = PolyElement(
        color = (color("tissue_simulation_shrinking_constantavdection"), 0.7), 
        strokecolor = :black,
        strokewidth = 1
        )
    elem_2 = PolyElement(
        color = (color("tissue_simulation_shrinking_constantavdection_d00.001123"), 0.7), 
        strokecolor = :black,
        strokewidth = 1
        )
    elem_3 = PolyElement(
        color = (color("tissue_simulation_shrinking_constantavdection_d00.003123"), 0.7), 
        strokecolor = :black,
        strokewidth = 1
        )

    # set the labels and title Here
    title = L"d_{0} \; \text{(μm⁻³)}"
    labels = ["0.001123", "0.002123", "0.003123"] # order here must match order of coloured elems

    Legend(
        fig[1, 1], 
        [elem_2, elem_1, elem_3], 
        labels, 
        title,
        orientation = :horizontal,
        nbanks = 3,
        tellheight = false,
        valign = :center
        )

    save("figures/20230331_comparing_starting_density_shrinking.pdf", fig)

   fig
end

# plot for the change in cell size

begin
    # figure for shrinking tissue with different densities
    df = DataFrame(CSV.File(raw"data\20230331_stochastic_shrinking.tsv", delim="\t"))

    function Base.parse(T::Type{Vector{Float64}}, s::String)
        s = replace(s, "[" => "")
        s = replace(s, "]" => "")
        s = replace(s, " " => "")
        s = eachsplit(s, ",")
        parse.(Float64, s) 
    end

    df[!, Symbol("Synchrony (r)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Synchrony (r)")]]

    df[!, Symbol("Mean frequency (/min)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Mean frequency (/min)")]]

    df[!, Symbol("Density (μm⁻³)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Density (μm⁻³)")]]

    df[!, Symbol("No. Cells")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("No. Cells")]]

    df[!, Symbol("xa (μm)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("xa (μm)")]]

    df[!, Symbol("Cell Additions")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Cell Additions")]]

    filter_7(x) = "tissue_simulation_shrinking_constantavdection" == x
    filter_12(x) = "tissue_simulation_shrinking_constantavdection_stepwisecellshrinking" == x # typo in the file

    fig = Figure(
        resolution = (800, 800),
        backgroundcolor=:grey90
        )

    ax1 = Axis(
        fig[2, 1], 
        # ax[1, 1],
        # aspect = AxisAspect(1), 
        xlabel = L"t \; \text{(mins)}", 
        ylabel = L"r",
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        title = "B",
        )

    ax2 = Axis(
        fig[2, 2], 
        # ax[1, 2],
        # aspect = AxisAspect(1), 
        xlabel = L"t \; \text{(mins)}", 
        ylabel = L"\frac{d\theta}{dt} \; \text{(min^{-1})}",
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        title = "C",
        )

    t_shrink = 253.6
    t_stop = 536.4

    vlines!(ax1, [t_shrink]; color=:black, linestyle=:dash)
    vlines!(ax2, [t_shrink]; color=:black, linestyle=:dash)

    vlines!(ax1, [t_stop]; color=:black, linestyle=:dot)
    vlines!(ax2, [t_stop]; color=:black, linestyle=:dot)

    # control colour here
    function color(type::AbstractString)

        if filter_7(type)

            return Makie.wong_colors()[2]

        elseif filter_12(type)

            return Makie.wong_colors()[5]

        end

    end

    for f in [filter_7, filter_12]

        dd = filter(Symbol("Simulation Type") => f, df)

        type = dd[1, Symbol("Simulation Type")]

        nt = get_first_NaN(dd[!, Symbol("Synchrony (r)")]) - 1

        # nt = nstop

        for row in eachrow(dd)

            row[Symbol("Synchrony (r)")] = row[Symbol("Synchrony (r)")][1:nt]
            row[Symbol("Mean frequency (/min)")] = row[Symbol(Symbol("Mean frequency (/min)"))][1:nt]

        end

        x = 0:nt-1

        med = median(dd[!, Symbol("Synchrony (r)")])
        lower_q = quantile(dd[!, Symbol("Synchrony (r)")], 0.25) 
        upper_q = quantile(dd[!, Symbol("Synchrony (r)")], 0.75) 


        CairoMakie.lines!(
            ax1,
            x,
            med,
            color=(color(type), 1.0)
        )

        band!(ax1, x, lower_q, upper_q, color=(color(type), 0.5))

        ylims!(ax1, 0, 1)

        med = median(dd[!, Symbol("Mean frequency (/min)")])
        lower_q = quantile(dd[!, Symbol("Mean frequency (/min)")], 0.25) 
        upper_q = quantile(dd[!, Symbol("Mean frequency (/min)")], 0.75) 

        CairoMakie.lines!(
            ax2,
            x,
            med,
            color=(color(type), 1.0)
        )

        band!(ax2, x, lower_q, upper_q, color=(color(type), 0.5))
    
    end

    xlims!(ax1, 0, 700)
    xlims!(ax2, 0, 700)

    # add axis showing the plot of cell size
    ax3 = Axis(
        fig[1, 2], 
        # ax[1, 2],
        # aspect = AxisAspect(1), 
        xlabel = L"t \; \text{(mins)}", 
        ylabel = L"d_{c} \; \text{(μm)}",
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        title = "A",
        )

    dc₀ = 9.2
    mdc = 0.2
    t_max = 690

    endless_shrinking = [cell_diameter(t, mdc, dc₀, t_shrink, t_max) for t = 0:t_max]
    stepwise_shrinking = [cell_diameter(t, mdc, dc₀, t_shrink, t_stop) for t = 0:t_max]

    vlines!(ax3, [t_shrink]; color=:black, linestyle=:dash)
    vlines!(ax3, [t_stop]; color=:black, linestyle=:dot)

    lines!(ax3, 0:t_max, endless_shrinking, label="continuous shrinking", color = color("tissue_simulation_shrinking_constantavdection"))
    lines!(ax3, 0:t_max, stepwise_shrinking, label="stepwise shrinking", color = color("tissue_simulation_shrinking_constantavdection_stepwisecellshrinking"))

    xlims!(ax3, 0, 700)
    ylims!(ax3, 0, 10)

    # annotate
    text!(ax1, L"t = t_{shrink}", position=(270, 0.05))
    text!(ax2, L"t = t_{shrink}", position=(270, 0.12))
    text!(ax3, L"t = t_{shrink}", position=(270, 1))

    text!(ax1, L"t = t_{stop}", position=(550, 0.05))
    text!(ax2, L"t = t_{stop}", position=(550, 0.12))
    text!(ax3, L"t = t_{stop}", position=(550, 1))

    # make legend
    elem_1 = PolyElement(
        color = (color("tissue_simulation_shrinking_constantavdection"), 0.7), 
        strokecolor = :black,
        strokewidth = 1
        )
    elem_2 = PolyElement(
        color = (color("tissue_simulation_shrinking_constantavdection_stepwisecellshrinking"), 0.7), 
        strokecolor = :black,
        strokewidth = 1
        )

    # set the labels and title Here
    title = "Legend"
    labels = ["continuous shrinking", "stepwise shrinking"] # order here must match order of coloured elems

    Legend(
        fig[3, 1:2], 
        [elem_1, elem_2], 
        labels, 
        title,
        orientation = :horizontal,
        nbanks = 1,
        tellheight = true,
        valign = :center
        )

    save("figures/20230331_shrinking_cells_stepwise_vs_cont.pdf", fig)

    fig
end


# figure for DDE simulations
begin
    
    data = DataFrame(CSV.File(raw"data\20230301_stochastic_delay_simulations.tsv", delim="\t"))

    dde_addition_position_fig = Figure(backgroundcolor=:grey90)

    synch_ax = Axis(
        dde_addition_position_fig[1, 1], 
        aspect = AxisAspect(1), 
        xlabel = "Position of new cells", 
        ylabel = L"r",
        xticks = (
            [1, 2, 3], ["Random", "DP", "DP + LV"]),
        title = "A",
        )

    ylims!(synch_ax, 0, 1)

    function filter_2(x)
        x in ["tissue_simulation_with_delay", "adding_posterior_zero_cells_simulation_with_delay", "adding_posterior_lateral_zero_cells_simulation_with_delay"]
    end

    function position_of_x(simulation_type::AbstractString)

        if simulation_type == "tissue_simulation_with_delay"

            return 1

        elseif simulation_type == "adding_posterior_zero_cells_simulation_with_delay"

            return 2

        elseif simulation_type == "adding_posterior_lateral_zero_cells_simulation_with_delay"

            return 3

        end

    end

    df = filter(Symbol("Simulation Type") => filter_2, data)

    CairoMakie.boxplot!(
        synch_ax, 
        position_of_x.(df[!, Symbol("Simulation Type")]), 
        df[!, Symbol("Synchrony (r)")]; 
        color = color_of_box.(df[!, Symbol("Simulation Type")]),  # this chooses which colour pair of the scheme to use.
        )

    freq_ax = Axis(
        dde_addition_position_fig[1, 2], 
        aspect = AxisAspect(1), 
        xlabel = "Position of new cells", 
        ylabel = L"\frac{d\theta}{dt} \; \text{(min^{-1})}",
        xticks = (
            [1, 2, 3], ["Random", "DP", "DP + LV"]),
        title = "B",
        )

    CairoMakie.boxplot!(
        freq_ax, 
        position_of_x.(df[!, Symbol("Simulation Type")]), 
        df[!, Symbol("Mean frequency (/min)")]; 
        color = color_of_box.(df[!, Symbol("Simulation Type")]),  # this chooses which colour pair of the scheme to use.
        )

    make_my_legend!(dde_addition_position_fig)

    rowsize!(dde_addition_position_fig.layout, 2, Relative(0.2)) # to move legend up

    trim!(dde_addition_position_fig.layout)

    # resize_to_layout!(addition_position_fig)

    dde_addition_position_fig

    save("figures/202230303_stochastic_simulations_cell_addition_fig_with_delay.pdf", dde_addition_position_fig)
end

# plot the frequency and synchrony for changing nτ
begin

    fname = "data/20231127_delaymitosis_changing_ntau_simulations.tsv"

    data = DataFrame(CSV.File(fname, delim="\t"))

    function nt(seed)
        if 5700 <= seed < 5800
            return 21
        elseif 6600 <= seed < 6700
            return 26
        elseif 6700 <= seed < 6800
            return 0
        end
    end

    data[!, :nt] = nt.(data[!, :Seed])

    fig = Figure()
    ax1 = Axis(fig[1,1], xlabel="τ (mins)", ylabel="r")
    ylims!(0, 1)
    boxplot!(ax1, data[!, :nt], data[!, Symbol("Synchrony (r)")])
    ax2 = Axis(fig[2,1], xlabel="τ (mins)", ylabel="dθ/dt")
    boxplot!(ax2, data[!, :nt], data[!, Symbol("Mean frequency (/min)")])
    save("figures/20231127_delaymitosis_changing_ntau_simulations.svg")
    fig

end

# plot the effect of length
begin

    colors = [Makie.wong_colors()[1], Makie.wong_colors()[3], Makie.wong_colors()[4]]

    path_rand = raw"data\20221212KappaMotilityParamSpaceTissueSimulation.tsv"
    df_rand = file_to_df(path_rand; drop_nans=true)

    path_dp = raw"data\20221212KappaMotilityParamSpaceDorsalPosteriorSimulation.tsv"
    df_dp = file_to_df(path_dp; drop_nans=true)

    path_dplv = raw"data\20221212KappaMotilityParamSpacePosteriorLateralSimulation.tsv"
    df_dplv = file_to_df(path_dplv; drop_nans=true)

    # isolate out length and plot a series of synchrony vs length, frequency vs PSM length. Then plot a parameter space
    # for κ versus L

    # isolate the anterior values
    df_rand[!, :ant_z] = [df_rand[!, :Zs][i][1] for i in eachindex(df_rand[!, :Zs])]
    df_dp[!, :ant_z] = [df_dp[!, :Zs][i][1] for i in eachindex(df_dp[!, :Zs])]
    df_dplv[!, :ant_z] = [df_dplv[!, :Zs][i][1] for i in eachindex(df_dplv[!, :Zs])]

    df_rand[!, :ant_f] = [df_rand[!, :freqs][i][1] for i in eachindex(df_rand[!, :freqs])]
    df_dp[!, :ant_f] = [df_dp[!, :freqs][i][1] for i in eachindex(df_dp[!, :freqs])]
    df_dplv[!, :ant_f] = [df_dplv[!, :freqs][i][1] for i in eachindex(df_dplv[!, :freqs])]

    # isolate lengths and kappas
    df_rand[!, :L] = [df_rand[!, Symbol("(phase, tissue)")][i][2].L for i in eachindex(df_rand[!, Symbol("(phase, tissue)")])]
    df_dp[!, :L] = [df_dp[!, Symbol("(phase, tissue)")][i][2].L for i in eachindex(df_dp[!, Symbol("(phase, tissue)")])]
    df_dplv[!, :L] = [df_dplv[!, Symbol("(phase, tissue)")][i][2].L for i in eachindex(df_dplv[!, Symbol("(phase, tissue)")])]

    df_rand[!, :kappa] = [df_rand[!, Symbol("(phase, tissue)")][i][1].k0 for i in eachindex(df_rand[!, Symbol("(phase, tissue)")])]
    df_dp[!, :kappa] = [df_dp[!, Symbol("(phase, tissue)")][i][1].k0 for i in eachindex(df_dp[!, Symbol("(phase, tissue)")])]
    df_dplv[!, :kappa] = [df_dplv[!, Symbol("(phase, tissue)")][i][1].k0 for i in eachindex(df_dplv[!, Symbol("(phase, tissue)")])]

    # combine into one DataFrame
    df_rand[!, :type] .= "Random"
    df_dp[!, :type] .= "DP"
    df_dplv[!, :type] .= "DP+LV"

    append!(df_rand, df_dp)
    append!(df_rand, df_dplv)

    # filter df by κ and remove changing vs simulations
    df = filter(:kappa => x -> x == 0.07, df_rand)
    df[!, :vs] = [df[!, Symbol("(phase, tissue)")][i][2].vs for i in eachindex(df[!, Symbol("(phase, tissue)")])]
    filter!(:vs => x -> x == 1.0, df)

    function dodge(x)
        if x == "Random"
            return 1
        elseif x == "DP"
            return 2
        elseif x == "DP+LV"
            return 3
        end
    end

    fig1 = Figure(size=(1200, 500))
    ax1 = Axis(fig1[1,1], xlabel="L (μm)", ylabel=L"r")
    ylims!(0, 1)
    boxplot!(ax1, 
    df[!, :L], 
    df[!, :ant_z], 
    dodge=dodge.(df[!, :type]), 
    dodgegap=50, 
    width=40, 
    color=[colors[i] for i in dodge.(df[!, :type])]
    )
    ax2 = Axis(fig1[2,1], xlabel="L (μm)", ylabel="dθ/dt (min⁻¹)")
    boxplot!(ax2, 
    df[!, :L], 
    df[!, :ant_f], 
    dodge=dodge.(df[!, :type]), 
    dodgegap=50, 
    width=40, 
    color=[colors[i] for i in dodge.(df[!, :type])],
    )


    # make squares of correct colour
    elem_1 = PolyElement(
        color = colors[1],
        strokecolor = :black,
        strokewidth = 1
        )

    elem_2 = PolyElement(
        color = colors[2],
        strokecolor = :black,
        strokewidth = 1
        )

    elem_3 = PolyElement(
        color = colors[3],
        strokecolor = :black,
        strokewidth = 1
        )

    # set the labels and title Here
    title = "Cell Addition"
    labels = ["Random", "DP", "DP+LV"] # order here must match order of coloured elems

    Legend(
        fig1[3, 1], 
        [elem_1, elem_2, elem_3], 
        labels, 
        title,
        orientation = :horizontal,
        nbanks = 1,
        tellheight = true,
        valign = :top
        )
    

    save("figures/20231201_length_figure.svg", fig1)
    fig1


    # now plot parameter space for κ vs L
    df_rand[!, :vs] = [df_rand[!, Symbol("(phase, tissue)")][i][2].vs for i in eachindex(df_rand[!, Symbol("(phase, tissue)")])]
    filter!(:vs => x -> x == 1, df_rand)

    n_κ = length(Set(df_rand[!, :kappa]))
    n_L = length(Set(df_rand[!, :L]))

    sorted_κ = sort!(collect(Set(df_rand[!, :kappa])))
    sorted_L = sort!(collect(Set(df_rand[!, :L])))

    randr_M = Matrix{Float64}(undef, n_κ, n_L)
    randf_M = Matrix{Float64}(undef, n_κ, n_L)

    dpr_M = Matrix{Float64}(undef, n_κ, n_L)
    dpf_M = Matrix{Float64}(undef, n_κ, n_L)

    dplvr_M = Matrix{Float64}(undef, n_κ, n_L)
    dplvf_M = Matrix{Float64}(undef, n_κ, n_L)

    for (i, L) in enumerate(sorted_L)

        s_df = filter(:L => x -> x == L, df_rand)

        for (j, κ) in enumerate(sorted_κ)

            k_df = filter(:kappa => x -> x == κ, s_df)

            # for each ingression scenario add to each matrix

            rd_df = filter(:type => x -> x == "Random", k_df)
            dp_df = filter(:type => x -> x == "DP", k_df)
            dplv_df = filter(:type => x -> x == "DP+LV", k_df)

            randr_M[j, i] = median(rd_df[!, :ant_z])
            randf_M[j, i] = median(rd_df[!, :ant_f])

            dpr_M[j, i] = median(dp_df[!, :ant_z])
            dpf_M[j, i] = median(dp_df[!, :ant_f])

            dplvr_M[j, i] = median(dplv_df[!, :ant_z])
            dplvf_M[j, i] = median(dplv_df[!, :ant_f])

        end

    end

    fig2 = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif")

    fig3 = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif")

    raxs = [Axis(fig2[1, i], aspect=AxisAspect(1)) for i = 1:3]

    faxs = [Axis(fig3[1, i], aspect=AxisAspect(1)) for i = 1:3]

    raxs[1].ylabel = "κ (min⁻¹)"
    faxs[1].ylabel = "κ (min⁻¹)"

    for j = 1:3

        raxs[j].xlabel = "L (μm)"
        faxs[j].xlabel = "L (μm)"
    
    end

    raxs[1].title = "Random"
    raxs[2].title = "DP"
    raxs[3].title = "DP + LV"

    faxs[1].title = "Random"
    faxs[2].title = "DP"
    faxs[3].title = "DP + LV"


    p1 = heatmap!(raxs[1],
        sorted_L,
        sorted_κ,
        randr_M',
        colorrange=(0, 1)
    )

    heatmap!(raxs[2],
        sorted_L,
        sorted_κ,
        dpr_M',
        colorrange=(0, 1)
    )

    heatmap!(raxs[3],
        sorted_L,
        sorted_κ,
        dplvr_M',
        colorrange=(0, 1)
    )

    cbar = Colorbar(fig2[1, 4], p1, label="r")

    rowsize!(fig2.layout, 1, Aspect(2, 1))

    # heatmap!(faxs[1],
    #     sorted_L,
    #     sorted_κ,
    #     randf_M',
    #     colorrange=(0, 1)
    # )

    # heatmap!(raxs[2],
    #     sorted_L,
    #     sorted_κ,
    #     dpr_M',
    #     colorrange=(0, 1)
    # )

    # heatmap!(raxs[3],
    #     sorted_L,
    #     sorted_κ,
    #     dplvr_M',
    #     colorrange=(0, 1)
    # )

    save("figures/20231201_lengthvskappa.svg", fig2)
    fig2

end

# define functions for counting the number of additions that have occured total,
# as well as the number of additions per unit length
begin

    function Base.parse(T::Type{Vector{Float64}}, s::AbstractString)
        s = replace(s, "[" => "")
        s = replace(s, "]" => "")
        # s = replace(s, " " => "")
        s = eachsplit(s, " ", keepempty=false)
        parse.(Float64, s) 
    end

    function Base.parse(T::Type{Matrix{Float64}}, s::AbstractString)
        
        s = eachsplit(s, ";")
        mapreduce(permutedims, vcat, parse.(Vector{Float64}, s)) 

    end

    # function that takes a matrix and counts the number of additions within a dc-width region of space and dt-long period of time
    function process_addition_matrix(M, x_min, x_max, dc, t_min, t_max, dt)

        X = [x for x in x_min:dc:x_max]
        T = [t for t in t_min:dt:t_max]

        out = zeros(length(X), length(T))

        for col in eachcol(M)

            x = col[1]
            t = col[4]

            if x <= maximum(X) && t <= maximum(T)

                # get the time and space coordinates
                i = findfirst(foo -> foo <= x <= foo + dc, X[1:end-1])
                j = findfirst(foo -> foo <= t <= foo + dt, T[1:end-1])

                # now update the number of added cells
                out[i, j] += 1

            end

        end

        out

    end

    function Statistics.median(V::Vector{Matrix{Float64}})

        # check all matrices are same size
        s = size(V[1])

        for M in V

            @assert size(M) == s

        end

        # create large array from V and then take median along each axis
        V = reduce(hcat, V)
        V = reshape(V, s[1], s[2], :)

        out = zeros(s)

        for i in axes(out, 1)
            for j in axes(out, 2)

                out[i, j] = Statistics.median(V[i, j, :])

            end
        end

        out

    end

    
    function Statistics.mean(V::Vector{Matrix{Float64}})

        # check all matrices are same size
        s = size(V[1])

        for M in V

            @assert size(M) == s

        end

        # create large array from V and then take median along each axis
        V = reduce(hcat, V)
        V = reshape(V, s[1], s[2], :)

        out = zeros(s)

        for i in axes(out, 1)
            for j in axes(out, 2)

                out[i, j] = Statistics.mean(V[i, j, :])

            end
        end

        out

    end

end


# make figure 
begin

    file = "data/20231206_changinglength_control_addition_tracking_simulations.tsv"

    data = DataFrame(CSV.File(file, delim="\t"))

    # Need to make a new column to avoid error
    data[:, :processed_added_cells] = [parse(Matrix{Float64}, row[:added_cells]) for row in eachrow(data)]

end

begin

    # first make a figure showing how the number of additions changes with the length of the tissue, across three different ingression scenarios

    simulation_types = [
        "tracking_addition_adding_posterior_zero_cells_simulation_Xc325.0",
        "tracking_addition_tissue_simulation_random_cells_Xc425.0",
        "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc225.0",
        "tracking_addition_adding_posterior_zero_cells_simulation_Xc225.0",
        "tracking_addition_tissue_simulation_random_cells_Xc325.0",
        "tracking_addition_adding_posterior_zero_cells_simulation_Xc425.0",
        "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc325.0",
        "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc425.0",
        "tracking_addition_tissue_simulation_random_cells_Xc225.0"
    ]

    df = data[[x in simulation_types for x in data[!, Symbol("Simulation Type")]], :]

    changing_length_fig = Figure(resolution=(750, 700))

    axs = [Axis(changing_length_fig[i, j], aspect=AxisAspect(1)) for i = 1:3, j = 1:3]

    axs[1, 1].title = "Random"
    axs[1, 2].title = "DP"
    axs[1, 3].title = "DP + LV"

    axs[3, 1].xlabel = "(x - xₐ)/L"
    axs[3, 2].xlabel = "(x - xₐ)/L"
    axs[3, 3].xlabel = "(x - xₐ)/L"

    axs[2, 1].ylabel="No. cells added (μm⁻¹⋅min⁻¹)"

    for ax in axs

        xlims!(ax, 0, 1)
        ylims!(ax, 0, 0.1)

    
    end

    # go through each sim type and position manually

    sim_type = "tracking_addition_tissue_simulation_random_cells_Xc225.0"

    ds = df[df[!, Symbol("Simulation Type")] .== sim_type, :]

    x_min = 25
    x_max = 310
    dc = 11
    t_min = 0
    t_max = 999
    dt = 1.0

    X = [(x - 25) / x_max for x in x_min:dc:x_max]
    T = t_min:dt:t_max

    V = [process_addition_matrix(M, x_min, x_max, dc, t_min, t_max, dt) for M in ds[!, :processed_added_cells]]

    M = Statistics.mean(V)

    y = [Statistics.mean(M[i, :]) for i in axes(M, 1)] / dc

    xd = 

    lines!(axs[1], X, y)

    sim_type = "tracking_addition_tissue_simulation_random_cells_Xc325.0"

    ds = df[df[!, Symbol("Simulation Type")] .== sim_type, :]

    x_min = 25
    x_max = 410
    dc = 11
    t_min = 0
    t_max = 999
    dt = 1.0

    X = [(x - 25) / x_max for x in x_min:dc:x_max]
    T = t_min:dt:t_max

    V = [process_addition_matrix(M, x_min, x_max, dc, t_min, t_max, dt) for M in ds[!, :processed_added_cells]]

    M = Statistics.mean(V)

    y = [Statistics.mean(M[i, :]) for i in axes(M, 1)] / dc

    lines!(axs[2], X, y)

    sim_type = "tracking_addition_tissue_simulation_random_cells_Xc425.0"

    ds = df[df[!, Symbol("Simulation Type")] .== sim_type, :]

    x_min = 25
    x_max = 510
    dc = 11
    t_min = 0
    t_max = 999
    dt = 1.0

    X = [(x - 25) / x_max for x in x_min:dc:x_max]
    T = t_min:dt:t_max

    V = [process_addition_matrix(M, x_min, x_max, dc, t_min, t_max, dt) for M in ds[!, :processed_added_cells]]

    M = Statistics.mean(V)

    y = [Statistics.mean(M[i, :]) for i in axes(M, 1)] / dc

    lines!(axs[3], X, y)

    sim_type = "tracking_addition_adding_posterior_zero_cells_simulation_Xc225.0"

    ds = df[df[!, Symbol("Simulation Type")] .== sim_type, :]

    x_min = 25
    x_max = 310
    dc = 11
    t_min = 0
    t_max = 999
    dt = 1.0

    X = [(x - 25) / x_max for x in x_min:dc:x_max]
    T = t_min:dt:t_max

    V = [process_addition_matrix(M, x_min, x_max, dc, t_min, t_max, dt) for M in ds[!, :processed_added_cells]]

    M = Statistics.mean(V)

    y = [Statistics.mean(M[i, :]) for i in axes(M, 1)] / dc

    lines!(axs[4], X, y)

    sim_type = "tracking_addition_adding_posterior_zero_cells_simulation_Xc325.0"

    ds = df[df[!, Symbol("Simulation Type")] .== sim_type, :]

    x_min = 25
    x_max = 410
    dc = 11
    t_min = 0
    t_max = 999
    dt = 1.0

    X = [(x - 25) / x_max for x in x_min:dc:x_max]
    T = t_min:dt:t_max

    V = [process_addition_matrix(M, x_min, x_max, dc, t_min, t_max, dt) for M in ds[!, :processed_added_cells]]

    M = Statistics.mean(V)

    y = [Statistics.mean(M[i, :]) for i in axes(M, 1)] / dc

    lines!(axs[5], X, y)

    sim_type = "tracking_addition_adding_posterior_zero_cells_simulation_Xc425.0"

    ds = df[df[!, Symbol("Simulation Type")] .== sim_type, :]

    x_min = 25
    x_max = 510
    dc = 11
    t_min = 0
    t_max = 999
    dt = 1.0

    X = [(x - 25) / x_max for x in x_min:dc:x_max]
    T = t_min:dt:t_max

    V = [process_addition_matrix(M, x_min, x_max, dc, t_min, t_max, dt) for M in ds[!, :processed_added_cells]]

    M = Statistics.mean(V)

    y = [Statistics.mean(M[i, :]) for i in axes(M, 1)] / dc

    lines!(axs[6], X, y)

    sim_type = "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc225.0"

    ds = df[df[!, Symbol("Simulation Type")] .== sim_type, :]

    x_min = 25
    x_max = 310
    dc = 11
    t_min = 0
    t_max = 999
    dt = 1.0

    X = [(x - 25) / x_max for x in x_min:dc:x_max]
    T = t_min:dt:t_max

    V = [process_addition_matrix(M, x_min, x_max, dc, t_min, t_max, dt) for M in ds[!, :processed_added_cells]]

    M = Statistics.mean(V)

    y = [Statistics.mean(M[i, :]) for i in axes(M, 1)] / dc

    lines!(axs[7], X, y)

    sim_type = "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc325.0"

    ds = df[df[!, Symbol("Simulation Type")] .== sim_type, :]

    x_min = 25
    x_max = 410
    dc = 11
    t_min = 0
    t_max = 999
    dt = 1.0

    X = [(x - 25) / x_max for x in x_min:dc:x_max]
    T = t_min:dt:t_max

    V = [process_addition_matrix(M, x_min, x_max, dc, t_min, t_max, dt) for M in ds[!, :processed_added_cells]]

    M = Statistics.mean(V)

    y = [Statistics.mean(M[i, :]) for i in axes(M, 1)] / dc

    lines!(axs[8], X, y)

    sim_type = "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc425.0"

    ds = df[df[!, Symbol("Simulation Type")] .== sim_type, :]

    x_min = 25
    x_max = 510
    dc = 11
    t_min = 0
    t_max = 999
    dt = 1.0

    X = [(x - 25) / x_max for x in x_min:dc:x_max]
    T = t_min:dt:t_max

    V = [process_addition_matrix(M, x_min, x_max, dc, t_min, t_max, dt) for M in ds[!, :processed_added_cells]]

    M = Statistics.mean(V)

    y = [Statistics.mean(M[i, :]) for i in axes(M, 1)] / dc

    lines!(axs[9], X, y)

    save("figures/20231215_number_of_additions_over_space.svg", changing_length_fig)
    save("figures/20231215_number_of_additions_over_space.pdf", changing_length_fig)

    changing_length_fig

end

# now a figure plotting the effect of displacing cell addition posteriorly 

begin

    colors = [Makie.wong_colors()[1], Makie.wong_colors()[3], Makie.wong_colors()[4]]

    simulation_types = [x for x in Set(data[!, Symbol("Simulation Type")]) if occursin("_xd", x)]

    df = data[[x in simulation_types for x in data[!, Symbol("Simulation Type")]], :]

    function type_sim(x)
        if occursin("tissue_simulation_random_cells", x)
            return "Random"
        elseif occursin("adding_posterior_zero_cells_simulation", x)
            return "DP"
        elseif occursin("adding_posterior_lateral_zero_cells_simulation", x)
            return "DP+LV"
        end
    end

    function L(x)
        start = findfirst("_Xc", x)[end]
        stop = findfirst("_", x[start+1:end])[1]
        parse(Float64, x[start+1:start + stop - 1]) - 25 + 85
    end 

    df[!, :type] = type_sim.(df[!, Symbol("Simulation Type")])
    df[!, :L] = L.(df[!, Symbol("Simulation Type")])

    function dodge(x)
        if x == "Random"
            return 1
        elseif x == "DP"
            return 2
        elseif x == "DP+LV"
            return 3
        end
    end

    fig1 = Figure(size=(1200, 500))
    ax1 = Axis(fig1[1,1], xlabel="L (μm)", ylabel=L"r")
    ylims!(0, 1)
    boxplot!(ax1, 
    df[!, :L], 
    df[!, Symbol("Synchrony (r)")], 
    dodge=dodge.(df[!, :type]), 
    dodgegap=50, 
    width=40, 
    color=[colors[i] for i in dodge.(df[!, :type])]
    )
    ax2 = Axis(fig1[2,1], xlabel="L (μm)", ylabel="dθ/dt (min⁻¹)")
    boxplot!(ax2, 
    df[!, :L], 
    df[!, Symbol("Mean frequency (/min)")], 
    dodge=dodge.(df[!, :type]), 
    dodgegap=50, 
    width=40, 
    color=[colors[i] for i in dodge.(df[!, :type])],
    )


    # make squares of correct colour
    elem_1 = PolyElement(
        color = colors[1],
        strokecolor = :black,
        strokewidth = 1
        )

    elem_2 = PolyElement(
        color = colors[2],
        strokecolor = :black,
        strokewidth = 1
        )

    elem_3 = PolyElement(
        color = colors[3],
        strokecolor = :black,
        strokewidth = 1
        )

    # set the labels and title Here
    title = "Cell Addition"
    labels = ["Random", "DP", "DP+LV"] # order here must match order of coloured elems

    Legend(
        fig1[3, 1], 
        [elem_1, elem_2, elem_3], 
        labels, 
        title,
        orientation = :horizontal,
        nbanks = 1,
        tellheight = true,
        valign = :top
        )
    

    save("figures/20231215_length_controllingxd_figure.svg", fig1)
    save("figures/20231215_length_controllingxd_figure.pdf", fig1)
    fig1

end

# and now a figure tracking the number of additions across each scenario of cell ingression and changing advection velocity

begin

    colors = [Makie.wong_colors()[1], Makie.wong_colors()[3], Makie.wong_colors()[4]]
    
    simulation_types = [x for x in Set(data[!, Symbol("Simulation Type")]) if (occursin("_Xc325.0", x) && !occursin("_Xc325.0_", x)) || occursin("_va", x)]

    df = data[[x in simulation_types for x in data[!, Symbol("Simulation Type")]], :]

    # calculate the number of additions in each simulation
    df[!, :n_added] = [size(M, 2) for M in df[!, :processed_added_cells]]

    function type_sim(x)
        if occursin("tissue_simulation_random_cells", x)
            return "Random"
        elseif occursin("adding_posterior_zero_cells_simulation", x)
            return "DP"
        elseif occursin("adding_posterior_lateral_zero_cells_simulation", x)
            return "DP+LV"
        end
    end

    df[!, :type] = type_sim.(df[!, Symbol("Simulation Type")])

    function va(x)
        if occursin("_va", x)
            return 3.0
        else
            return 1.67
        end
    end

    df[!, :va] = va.(df[!, Symbol("Simulation Type")])

    function dodge(x)
        if x == "Random"
            return 1
        elseif x == "DP"
            return 2
        elseif x == "DP+LV"
            return 3
        end
    end

    fig = Figure()

    ax = Axis(fig[1,1], xlabel="vₐ (μm·min⁻¹)", ylabel="Number of cells added per simulation")
    ax.xticks = ([1.67, 3.0], ["1.67", "3.0"])
    boxplot!(ax, 
    df[!, :va], 
    df[!, :n_added], 
    dodge=dodge.(df[!, :type]),
    # width=40, 
    color=[colors[i] for i in dodge.(df[!, :type])]
    )

    elem_1 = PolyElement(
        color = colors[1],
        strokecolor = :black,
        strokewidth = 1
        )

    elem_2 = PolyElement(
        color = colors[2],
        strokecolor = :black,
        strokewidth = 1
        )

    elem_3 = PolyElement(
        color = colors[3],
        strokecolor = :black,
        strokewidth = 1
        )

    # set the labels and title Here
    title = "Cell Addition"
    labels = ["Random", "DP", "DP+LV"] # order here must match order of coloured elems

    Legend(
        fig[2, 1], 
        [elem_1, elem_2, elem_3], 
        labels, 
        title,
        orientation = :horizontal,
        nbanks = 1,
        tellheight = true,
        valign = :top
        )

    save("figures/20231215_number_of_additions_total_changing_advection.svg", fig)
    save("figures/20231215_number_of_additions_total_changing_advection.pdf", fig)

    fig

end

begin

    file = "data/20231216_stochastic_shrinking.tsv"

    df = DataFrame(CSV.File(file, delim="\t"))

    function Base.parse(T::Type{Vector{Float64}}, s::String)
        s = replace(s, "[" => "")
        s = replace(s, "]" => "")
        s = replace(s, " " => "")
        s = eachsplit(s, ",")
        parse.(Float64, s) 
    end

    # now add the older data with increasing density into df
    file = "data/20230224_stochastic_shrinking.tsv"

    odf = DataFrame(CSV.File(file, delim="\t"))
    odf = odf[odf[!, Symbol("Simulation Type")] .== "tissue_simulation_shrinking_constantavdection", :]

    df = vcat(df, odf)

    df[!, Symbol("Synchrony (r)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Synchrony (r)")]]

    df[!, Symbol("Mean frequency (/min)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Mean frequency (/min)")]]

    df[!, Symbol("Density (μm⁻³)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Density (μm⁻³)")]]

    df[!, Symbol("No. Cells")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("No. Cells")]]

    df[!, Symbol("xa (μm)")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("xa (μm)")]]

    df[!, Symbol("Cell Additions")] = [parse(Vector{Float64}, s) for s in df[!, Symbol("Cell Additions")]]

    Set(df[!, Symbol("Simulation Type")])

    fig = Figure(
        resolution = (1200, 800),
        # backgroundcolor=:grey90
        )

    ax1 = Axis(
        fig[1, 1], 
        # ax[1, 1],
        # aspect = AxisAspect(1), 
        xlabel = L"t \; \text{(min)}", 
        ylabel = L"r",
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        title = "B",
        )

    ax2 = Axis(
        fig[1, 2], 
        # ax[1, 2],
        # aspect = AxisAspect(1), 
        xlabel = L"t \; \text{(min)}", 
        ylabel = L"\frac{d\theta}{dt} \; \text{(min^{-1})}",
        xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
        title = "C",
        )

    t_shrink = 253.6

    vlines!(ax1, [t_shrink]; color=:black, linestyle=:dash)
    vlines!(ax2, [t_shrink]; color=:black, linestyle=:dash)

    # control colour here
    function color(type::AbstractString)

        if type == "tissue_simulation_shrinking_constantavdection"

            return Makie.wong_colors()[2]

        elseif type == "tissue_simulation_shrinking_constantavdection_constantdensity"

            return Makie.wong_colors()[3]

        elseif type == "tissue_simulation_shrinking_constantavdection_constantdensity_constant_dc"

            return Makie.wong_colors()[6]

        end

    end

    for type in ["tissue_simulation_shrinking_constantavdection", "tissue_simulation_shrinking_constantavdection_constantdensity", "tissue_simulation_shrinking_constantavdection_constantdensity_constant_dc"]

        dd = filter(Symbol("Simulation Type") => x -> x == type, df)

        nt = get_first_NaN(dd[!, Symbol("Synchrony (r)")]) - 1

        # nt = nstop

        for row in eachrow(dd)

            row[Symbol("Synchrony (r)")] = row[Symbol("Synchrony (r)")][1:nt]
            row[Symbol("Mean frequency (/min)")] = row[Symbol(Symbol("Mean frequency (/min)"))][1:nt]

        end

        x = 0:nt-1

        med = median(dd[!, Symbol("Synchrony (r)")])
        lower_q = quantile(dd[!, Symbol("Synchrony (r)")], 0.25) 
        upper_q = quantile(dd[!, Symbol("Synchrony (r)")], 0.75) 


        CairoMakie.lines!(
            ax1,
            x,
            med,
            color=(color(type), 1.0)
        )

        band!(ax1, x, lower_q, upper_q, color=(color(type), 0.5))

        ylims!(ax1, 0, 1)

        med = median(dd[!, Symbol("Mean frequency (/min)")])
        lower_q = quantile(dd[!, Symbol("Mean frequency (/min)")], 0.25) 
        upper_q = quantile(dd[!, Symbol("Mean frequency (/min)")], 0.75) 

        CairoMakie.lines!(
            ax2,
            x,
            med,
            color=(color(type), 1.0)
        )

        band!(ax2, x, lower_q, upper_q, color=(color(type), 0.5))
    
    end

    xlims!(ax1, 0, 700)
    xlims!(ax2, 0, 700)

    # make legend
    elem_1 = PolyElement(
        color = color("tissue_simulation_shrinking_constantavdection"),
        strokecolor = :black,
        strokewidth = 1
        )

    elem_2 = PolyElement(
        color = color("tissue_simulation_shrinking_constantavdection_constantdensity"),
        strokecolor = :black,
        strokewidth = 1
        )

    elem_3 = PolyElement(
        color = color("tissue_simulation_shrinking_constantavdection_constantdensity_constant_dc"),
        strokecolor = :black,
        strokewidth = 1
        )
    # set the labels and title Here
    title = "Legend"
    labels = ["Increasing ρ₀ + decreasing dc", "Constant ρ₀ + decreasing dc", "Constant ρ₀ + constant dc"] # order here must match order of coloured elems

    Legend(
        fig[3, 1:2], 
        [elem_1, elem_2, elem_3], 
        labels, 
        title,
        orientation = :horizontal,
        nbanks = 3,
        tellheight = true,
        valign = :center
        )

    save("figures/20231216_increasing_vs_constant_density.pdf", fig)

    fig

end


