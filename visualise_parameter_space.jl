# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Plots outputs of systematic parameter sweeps as 2D heatmaps for synchrony or frequency

include("visualise_parameter_space_functions.jl")

begin

    for x_ind = 1:27

        visualise_advection_parameter_space(x_ind)

        # visualise_kappa_motility_parameter_space(x_ind)

        visualise_motility_parameter_space(x_ind)

    end

end

# begin
    
#     path = raw"data\20221130AdvectionParamSpaceDorsalPosteriorSimulation.tsv"

#     df = file_to_df(path; drop_nans=true)

#     result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

#     for i in eachindex(result)

#         result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

#     end

#     plots = visualise_advection_parameter_space(result; colorbar=false)

#     # now make a big column of plots for each variable

#     N_plot = plot([pl[1] for pl in plots]..., layout = grid(11, 1), size=(350, 3000))


#     cbar = Plots.scatter(
#         [0, 0], [0, 1], zcolor=[0, 1], clims=(0, 1), xlims=(1.05, 1.06), label="", c=:viridis, 
#         colorbartitle="r", framestyle=:none,
#         )

#     l = @layout [Plots.grid(11, 1) a{0.035w}]

#     r_plot = plot([pl[2] for pl in plots]..., layout = grid(11, 1), size=(300, 3000))

#     freq_plot = plot([pl[3] for pl in plots]..., layout = grid(11, 1), size=(350, 3000))

# end

begin
    
    path = raw"data\20221130AdvectionParamSpacePosteriorLateralSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    plots = visualise_advection_parameter_space(result; colorbar=false)

    # now make a big column of plots for each variable

    N_plot = plot([pl[1] for pl in plots]..., layout = grid(11, 1), size=(350, 3000))


    cbar = Plots.scatter(
        [0, 0], [0, 1], zcolor=[0, 1], clims=(0, 1), xlims=(1.05, 1.06), label="", c=:viridis, 
        colorbartitle="r", framestyle=:none,
        )

    l = @layout [Plots.grid(11, 1) a{0.035w}]

    r_plot = plot([pl[2] for pl in plots]..., layout = grid(11, 1), size=(300, 3000))

    freq_plot = plot([pl[3] for pl in plots]..., layout = grid(11, 1), size=(350, 3000))

end

# load in the default params to annotate on figures
begin
    tissue = TissueParams(
        11.0, # cell diameter in μm
        8.71, # coefficient for intercellular force
        1.67, # cell advection velocity in anterior 
        3.0, # cell advection speed in posterior
        0.7, # % x coordinate dividing the PSM into anterior and posterior
        300.0 + 60.0 + 25.0, # length of the tissue in μm
        0., # x-coordinate of the most anterior cell
        0.4, # length scale of the cell motility gradient. Inf64 to set motility constant across tissue
        1.0, # intrinsic cell motility in the posterior of the tissue. This is from Uriu et al. 2017, not 2021! (1.0 in latter)
        3.0, # coefficient for shape of motility gradient (needs to be 1.0 for no gradient)
        0.026, # phase noise for cell directionality (/min)
        20.0, # coefficient for boundary force
        1.0, # length-scale of boundary force,
        300.0, # x coord of centre of PSM torus 
        60.0 + 25., # y coord of centre of PSM torus 
        25., # z coord of centre of PSM torus 
        25.0, # radius of PSM cylinder
        60.0, # radius of PSM torus
        0.0015, # desired density (from C code)
        2 * 883 + 555, # number of cells in tissue
        2 * 883 + 555, # total number of cells that have been in the tissue
        187.5 - 15., # duration of growth phase (mins)
        15., # duration of M phase (mins)
        100., # anterior limit of cell addition
        0.07, # coupling to the segmentation clock
        )
    
    phase = PhaseParams(
        0.07,
        0 * 0.0013,
        0.2094,
        0.66, # 0.66 by default
        3.07,
        0.0011,
        )
end

# code for motility gradient shape param sweep figure
begin
    # load input
    path = raw"data\20221130AdvectionParamSpaceTissueSimulation.tsv"

    # could write a function to automatically read out the param space
    # easier just to specify here though
    xvs = 0.2:0.05:0.6
    hs = 1:0.5:5

    result = file_to_vect(path)

    Nfig, rfig, freqfig = VisualiseMotilityParameterSpace(result, xvs, hs; colorbar=false)

    title!(rfig, "Random")

    # add in cross marking default parameters 
    scatter!(rfig, [tissue.Xᵥ], [tissue.h], marker=:cross, label=false, color=:black)
    scatter!(freqfig, [tissue.Xᵥ], [tissue.h], marker=:cross, label=false, color=:black)
    
    # load input
    path = raw"data/20220428MotilityParamSpaceEarlyZFMorphogenesisSimulation.txt"

    # could write a function to automatically read out the param space
    # easier just to specify here though
    xvs = 0.2:0.05:0.6
    hs = 1:0.5:5

    result = file_to_vect(path)

    Nfigdp, rfigdp, freqfigdp = VisualiseMotilityParameterSpace(result, xvs, hs; colorbar=false)

    title!(rfigdp, "DP + LV")

    scatter!(rfigdp, [tissue.Xᵥ], [tissue.h], marker=:cross, label=false, color=:black)
    scatter!(freqfigdp, [tissue.Xᵥ], [tissue.h], marker=:cross, label=false, color=:black)

    path = raw"data\20220719MotilityParamSpaceDorsalPosteriorSimulation.txt"

    # could write a function to automatically read out the param space
    # easier just to specify here though
    xvs = 0.2:0.05:0.6
    hs = 1:0.5:5

    result = file_to_vect(path)

    Nfigp, rfigp, freqfigp = VisualiseMotilityParameterSpace(result, xvs, hs; colorbar=false)

    title!(rfigp, "Dorso-Posterior")

    scatter!(rfigp, [tissue.Xᵥ], [tissue.h], marker=:cross, label=false, color=:black)
    scatter!(freqfigp, [tissue.Xᵥ], [tissue.h], marker=:cross, label=false, color=:black)

    # create two colorbars
    cbarr = Plots.scatter(
            [0, 0], [0, 1], zcolor=[0, 0], clims=(0, 1), xlims=(1, 1.1), label="", c=:viridis, 
            colorbartitle="r", framestyle=:none, right_margin = 6Plots.mm
            )

    cbartitle = L"$\frac{d \theta}{dt} - \frac{1}{d_{c}} \int_{x_{a}}^{x_{a} + d_{c}} \omega(x) dx$"

    cbartitle = L"$\Delta\omega$"

    cbarf = Plots.scatter(
        [0, 0], [0, 1], zcolor=[0, 0], clims=(-0.02, 0.02), xlims=(1, 1.1), label="", c=:balance, 
        colorbartitle=cbartitle, framestyle=:none,
        thickness_scaling = 1, right_margin = 6Plots.mm
        )

    # combine into one figure
    l = @layout [grid(1, 3) a{0.035w}]

    tfig = plot(rfig, rfigp, rfigdp, cbarr, layout=l)

    l = @layout [grid(1, 3) a{0.035w}]

    bfig = plot(freqfig, freqfigp, freqfigdp, cbarf, layout=l)

    bigfig = plot(tfig, bfig, layout=grid(2, 1), size=(1100, 600), left_margin=5Plots.mm, bottom_margin=5Plots.mm)

    savefig(bigfig, "figures/202204280717MotilityParamSpaceControlvsDorsalPosterior.svg")
end


# code for kappa motility param sweep figure
begin
    # load input
    path = "data/20220617KappaMotilityParamSpaceControlSimulation.txt"

    # could write a function to automatically read out the param space
    # easier just to specify here though
    v0s = 0:0.5:5
    κs = 0:0.01:0.14

    result = file_to_vect(path)

    Nfig, rfig, freqfig = VisualiseKappaMotilityParameterSpace(result, κs, v0s; colorbar=false)

    title!(rfig, "Random")

    scatter!(rfig, [tissue.vs], [phase.k0], marker=:cross, label=false, color=:black)
    scatter!(freqfig, [tissue.vs], [phase.k0], marker=:cross, label=false, color=:black)

    # load input
    path = raw"data\20220719KappaMotilityParamSpacePosteriorDorsalSimulation.txt"

    # could write a function to automatically read out the param space
    # easier just to specify here though
    v0s = 0:0.5:5
    κs = 0:0.01:0.14

    result = file_to_vect(path)

    Nfigdp, rfigdp, freqfigdp = VisualiseKappaMotilityParameterSpace(result, κs, v0s; colorbar=false)

    title!(rfigdp, "Dorso-Posterior")

    scatter!(rfigdp, [tissue.vs], [phase.k0], marker=:cross, label=false, color=:black)
    scatter!(freqfigdp, [tissue.vs], [phase.k0], marker=:cross, label=false, color=:black)

    path = raw"data\20220717KappaMotilityParamSpacePosteriorLateralCellsSimulation.txt"

    # could write a function to automatically read out the param space
    # easier just to specify here though
    v0s = 0:0.5:5
    κs = 0:0.01:0.14

    result = file_to_vect(path)

    Nfigdplv, rfigdplv, freqfigdplv = VisualiseKappaMotilityParameterSpace(result, κs, v0s; colorbar=false)

    title!(rfigdplv, "DP + LV")

    scatter!(rfigdplv, [tissue.vs], [phase.k0], marker=:cross, label=false, color=:black)
    scatter!(freqfigdplv, [tissue.vs], [phase.k0], marker=:cross, label=false, color=:black)

    # create two colorbars
    cbarr = Plots.scatter(
            [0, 0], [0, 1], zcolor=[0, 0], clims=(0, 1), xlims=(1, 1.1), label="", c=:viridis, 
            colorbartitle="r", framestyle=:none, right_margin = 6Plots.mm
            )

    cbartitle = L"$\frac{d \theta}{dt} - \frac{1}{d_{c}} \int_{x_{a}}^{x_{a} + d_{c}} \omega(x) dx$"

    cbartitle = L"$\Delta\omega$"

    cbarf = Plots.scatter(
        [0, 0], [0, 1], zcolor=[0, 0], clims=(-0.02, 0.02), xlims=(1, 1.1), label="", c=:balance, 
        colorbartitle=cbartitle, framestyle=:none,
        thickness_scaling = 1, right_margin = 6Plots.mm
        )

    # combine into one figure
    l = @layout [grid(1, 3) a{0.035w}]

    tfig = plot(rfig, rfigdp, rfigdplv, cbarr, layout=l)

    l = @layout [grid(1, 3) a{0.035w}]

    bfig = plot(freqfig, freqfigdp, freqfigdplv, cbarf, layout=l)

    bigfig = plot(tfig, bfig, layout=grid(2, 1), size=(600, 1000))

    savefig(bigfig, "figures/20220617KappaMotilityControlvsDorsalPosterior.svg")
end


# l = @layout [1 2; 3 4]

# tissue.Xᵥ = 0.2
# tissue.h = 1
# p1 = Plots.plot(0:tissue.L, [intrinsic_motion(x, tissue) for x in 0:tissue.L], label=false)
# Plots.xlabel!("x (μm)")
# Plots.ylabel!("Motility (min⁻¹)")
# Plots.annotate!(75, 0.9, "Xᵥ = $(tissue.Xᵥ)\nh = $(tissue.h)", 10)

# tissue.Xᵥ = 0.2
# tissue.h = 5
# p2 = Plots.plot(0:tissue.L, [intrinsic_motion(x, tissue) for x in 0:tissue.L], label=false)
# Plots.xlabel!("x (μm)")
# Plots.ylabel!("Motility (min⁻¹)")
# Plots.annotate!(75, 0.9, "Xᵥ = $(tissue.Xᵥ)\nh = $(tissue.h)", 10)

# tissue.Xᵥ = 0.6
# tissue.h = 1
# p3 = Plots.plot(0:tissue.L, [intrinsic_motion(x, tissue) for x in 0:tissue.L], label=false)
# Plots.xlabel!("x (μm)")
# Plots.ylabel!("Motility (min⁻¹)")
# Plots.annotate!(75, 0.9, "Xᵥ = $(tissue.Xᵥ)\nh = $(tissue.h)", 10)

# tissue.Xᵥ = 0.6
# tissue.h = 5
# p4 = Plots.plot(0:tissue.L, [intrinsic_motion(x, tissue) for x in 0:tissue.L], label=false)
# Plots.xlabel!("x (μm)")
# Plots.ylabel!("Motility (min⁻¹)")
# Plots.annotate!(75, 0.9, "Xᵥ = $(tissue.Xᵥ)\nh = $(tissue.h)", 10)

# pl = Plots.plot(p3, p4, p1, p2)
# savefig(pl, "figures/motilityspace.svg")