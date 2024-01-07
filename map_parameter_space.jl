# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Code for performing systematic parameter sweeps, using a HPC cluster that uses Slurm

println("Version is Julia $(VERSION)")

# code for parallelised slurm job
using Distributed, ClusterManagers

available_workers = parse(Int, ENV["SLURM_NTASKS"])

addprocs_slurm(available_workers)

println("Number of cores: ", nprocs())
println("Number of workers: ", nworkers())

# load packages and define functions
@everywhere begin
    using Statistics
    include("experiments.jl") # code for simulating
    include("save_simulations.jl") # load code for saving data as .csv
end


# couple of functions for analysis
#  """
# Get slices of the PSM along the AP axis.

# Argument side specifies whether taking a domain from the left or right hand side of the psm.
# """
@everywhere function thin_domain(x, tracks::Matrix{Float64}, m, Δx, r, R; side="l")
    cells_in_domain = Int64[]
    for i = 1:size(tracks, 1)
        if side == "l"
            if x + (m-1) * Δx <= tracks[i, 1] <= x + m * Δx && 0 <= tracks[i, 2] <= 2 * r && 0 <= tracks[i, 3] <= 2 * r
                push!(cells_in_domain, i)
            end
        elseif side == "r"
            if x + (m-1) * Δx <= tracks[i, 1] <= x + m * Δx && 2 * R <= tracks[i, 2] <= 2 * (r + R) && 0 <= tracks[i, 3] <= 2 * r
                push!(cells_in_domain, i)
            end
        end
    end
    domain = Matrix{Float64}(undef, length(cells_in_domain), size(tracks, 2))
    for (j, i) in enumerate(cells_in_domain)
        domain[j, :] .= tracks[i, :]
    end
    domain
end

# DEFUNCT
@everywhere function N(tracks::Vector{Matrix{Float64}}, tissue::TissueParams)
    mean(thin_domain(tissue.xa, tracks[end], 1, tissue.dc, tissue.r, tissue.R, side="l")[:, 6]) / 2π
end

# @everywhere function Z(tracks::Vector{Matrix{Float64}}, tissue::TissueParams)
#     # return a vector with average Z at that point in space, across the last 100 mins
#     n = 0.

#     out = zeros(Float64, length(tissue.xa:tissue.dc:tissue.Xc-tissue.dc))

#     for t = length(tracks) - 99 : length(tracks)

#         n += 1.

#         for (i, x) in enumerate(tissue.xa:tissue.dc:tissue.Xc-tissue.dc)

#             out[i] += phase_order(thin_domain(x, tracks[end], 1, tissue.dc, tissue.r, tissue.R, side="l")[:, 6])

#         end

#     end
#     out / n
# end

# Gets the synchrony along the x axis of the left hand side of the PSM
@everywhere function Z(tracks::Vector{Matrix{Float64}}, tissue::TissueParams)

    out = zeros(Float64, length(tissue.xa:tissue.dc:tissue.Xc-tissue.dc))

    for (i, x) in enumerate(tissue.xa:tissue.dc:tissue.Xc-tissue.dc)

        out[i] = phase_order(thin_domain(x, tracks[end], 1, tissue.dc, tissue.r, tissue.R, side="l")[:, 6])

    end

    out

end

# @everywhere function freq(tracks::Vector{Matrix{Float64}}, tissue::TissueParams)
#     n = 0.

#     out = zeros(Float64, length(tissue.xa:tissue.dc:tissue.Xc-tissue.dc))

#     for t = length(tracks) - 99 : length(tracks)

#         n += 1.

#         for (i, x) in enumerate(tissue.xa:tissue.dc:tissue.Xc-tissue.dc)

#             out[i] += mean(thin_domain(x, tracks[end], 1, tissue.dc, tissue.r, tissue.R, side="l")[:, 7])

#         end

#     end
#     out / n
# end

# Gets the mean frequency along the x axis of the left hand side of the PSM
@everywhere function freq(tracks::Vector{Matrix{Float64}}, tissue::TissueParams)

    out = zeros(Float64, length(tissue.xa:tissue.dc:tissue.Xc-tissue.dc))

    for (i, x) in enumerate(tissue.xa:tissue.dc:tissue.Xc-tissue.dc)

        out[i] = mean(thin_domain(x, tracks[end], 1, tissue.dc, tissue.r, tissue.R, side="l")[:, 7])

    end

    out

end

# function Z(tracks::Vector{Matrix{Float64}}, tissue::TissueParams)
#     out = Vector{Float64}(undef, length(tissue.xa:tissue.dc:tissue.Xc-tissue.dc))
#     for (i, x) in enumerate(tissue.xa:tissue.dc:tissue.Xc-tissue.dc)
#         out[i] = phase_order(thin_domain(x, tracks[end], 1, tissue.dc, tissue.r, tissue.R, side="l")[:, 6])
#     end
#     out
# end

@everywhere function Statistics.mean(X::Vector{Vector{Float64}})
    x = zeros(length(X[1]))
    for i in eachindex(X)
        x += X[i]
    end
    x / length(X)
end

# Functions to call the different experiments

@everywhere function TissueSimulation(seed::Int64, phase::PhaseParams, tissue::TissueParams)
    try
        # control time parameters from here
        t_min, t_max, dt = 0., 999., 0.01
        # need to copy params objects as they are modified dynamically
        new_phase, new_tissue = copy(phase), copy(tissue)
        tracks = frequency_tracking_tissue_simulation(seed, t_min, t_max, dt, new_tissue, new_phase)
        # return
        return N(tracks, new_tissue), Z(tracks, new_tissue), freq(tracks, new_tissue), (phase, tissue)
    catch e
        println("Error in TissueSimulation\nparams = $phase, $tissue")
        rethrow(e)
    end
end

@everywhere function PosteriorLateralTissueSimulation(seed::Int64, phase::PhaseParams, tissue::TissueParams)
    try
        # control time parameters from here
        t_min, t_max, dt = 0., 999., 0.01
        # need to copy params objects as they are modified dynamically
        new_phase, new_tissue = copy(phase), copy(tissue)
        tracks = frequency_tracking_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, new_tissue, new_phase)
        # return
        return N(tracks, new_tissue), Z(tracks, new_tissue), freq(tracks, new_tissue), (phase, tissue)
    catch e
        println("Error in PosteriorLateralTissueSimulation\nparams = $phase, $tissue")
        rethrow(e)
    end
end

@everywhere function PosteriorTissueSimulation(seed::Int64, phase::PhaseParams, tissue::TissueParams)
    try
        # control time parameters from here
        t_min, t_max, dt = 0., 999., 0.01
        # need to copy params objects as they are modified dynamically
        new_phase, new_tissue = copy(phase), copy(tissue)
        tracks = frequency_tracking_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, new_tissue, new_phase)
        # return
        return N(tracks, new_tissue), Z(tracks, new_tissue), freq(tracks, new_tissue), (phase, tissue)
    catch e
        println("Error in PosteriorTissueSimulation\nparams = $phase, $tissue")
        rethrow(e)
    end
end

# now create a function to simulate and analyse
@everywhere function delay_tissue_simulation(seed::Int64, phase::PhaseParams, tissue::TissueParams)
    try
        GC.gc()
        # control time parameters from here
        t_min, t_max, dt = 0., 999., 0.01
        # need to copy params objects as they are modified dynamically
        new_phase, new_tissue = copy(phase), copy(tissue)
        tracks = tissue_simulation_with_delay(seed, t_min, t_max, dt, new_tissue, new_phase)
        # return
        return N(tracks, new_tissue), Z(tracks, new_tissue), freq(tracks, new_tissue), (phase, tissue)
    catch e
        println("Error in TissueSimulation\nparams = $phase, $tissue")
        rethrow(e)
    end
end

@everywhere function delay_posterior_lateral_simulation(seed::Int64, phase::PhaseParams, tissue::TissueParams)
    try
        GC.gc()
        # control time parameters from here
        t_min, t_max, dt = 0., 999., 0.01
        # need to copy params objects as they are modified dynamically
        new_phase, new_tissue = copy(phase), copy(tissue)
        tracks = adding_posterior_lateral_zero_cells_simulation_with_delay(seed, t_min, t_max, dt, new_tissue, new_phase)
        # return
        return N(tracks, new_tissue), Z(tracks, new_tissue), freq(tracks, new_tissue), (phase, tissue)
    catch e
        println("Error in PosteriorLateralTissueSimulation\nparams = $phase, $tissue")
        rethrow(e)
    end
end

@everywhere function delay_posterior_simulation(seed::Int64, phase::PhaseParams, tissue::TissueParams)
    try
        GC.gc()
        # control time parameters from here
        t_min, t_max, dt = 0., 999., 0.01
        # need to copy params objects as they are modified dynamically
        new_phase, new_tissue = copy(phase), copy(tissue)
        tracks = adding_posterior_zero_cells_simulation_with_delay(seed, t_min, t_max, dt, new_tissue, new_phase)
        # return
        return N(tracks, new_tissue), Z(tracks, new_tissue), freq(tracks, new_tissue), (phase, tissue)
    catch e
        println("Error in PosteriorTissueSimulation\nparams = $phase, $tissue")
        rethrow(e)
    end
end

@everywhere function mitosing_simulation(seed::Int64, phase::PhaseParams, tissue::TissueParams)
    try
        GC.gc()
        # control time parameters from here
        t_min, t_max, dt = 0., 999., 0.01
        # need to copy params objects as they are modified dynamically
        new_phase, new_tissue = copy(phase), copy(tissue)
        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, new_tissue, new_phase)
        # return
        return N(tracks, new_tissue), Z(tracks, new_tissue), freq(tracks, new_tissue), (phase, tissue)
    catch e
        println("Error in mitosing simulation\nparams = $phase, $tissue")
        rethrow(e)
    end
end

# now define parameter space
# NOTE - below code will need customising per simulation

begin
    param_space = Tuple{Int64, PhaseParams, TissueParams}[]

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

    # shift the tissue 25 from the origin to allow 'escaped' cells to be simulated on the lattice
    tissue.xa += 25.
    tissue.Xc += 25.
    tissue.Yc += 25.
    tissue.Zc += 25.
    tissue.Lx += 25.
    tissue.Ly += 25.
    tissue.Lz += 25.

    # xqs = 0:0.1:1
    # γas = 0:0.835:5.845
    # γps = 0:0.5:6

    # for xq in xqs

    #     for γₐ ∈ γas

    #         for γₚ ∈ γps

    #             for seed = 1:100

    #                 np, nt = copy(phase), copy(tissue)

    #                 nt.xq, nt.γₐ, nt.γₚ = xq, γₐ, γₚ
                    
    #                 push!(param_space, (seed, np, nt))
                
    #             end

    #         end

    #     end

    # end

    # κs = 0:0.01:0.14
    # vss = 0:0.5:5
    # # Xcs = 100:100:500
    # τs = 1:100:2101

    # for κ ∈ κs

    #     for vs ∈ vss

    #         for τ ∈ τs

    #             for seed = 1:100

    #                 np, nt = copy(phase), copy(tissue)

    #                 # np.k0, nt.vs, nt.Xc = κ, vs, Xc + nt.xa

    #                 # nt.L = nt.Xc + nt.R + nt.r - nt.xa

    #                 # nt.Lx = nt.L + 25. + tissue.xa

    #                 # nt.N_cells = 2 * floor(nt.ρ₀ * π * nt.r^2 * (nt.Xc - nt.xa)) + floor(nt.ρ₀ * π^2 * nt.r^2 * nt.R)
    #                 # nt.total_cells = nt.N_cells

    #                 np.k0, nt.vs, np.nτ = κ, vs, τ
                    
    #                 push!(param_space, (seed, np, nt))
                
    #             end

    #         end

    #     end

    # end

    # xvs = 0:0.1:1
    # hs = 0:0.5:5
    # vss = 0:0.5:5

#     for xv ∈ xvs

#         for h ∈ hs

#             # for vs ∈ vss

#                 for seed = 1:100

#                     np, nt = copy(phase), copy(tissue)

#                     # nt.Xᵥ, nt.vs, nt.h = xv, vs, h

#                     nt.Xᵥ, nt.h = xv, h
                    
#                     push!(param_space, (seed, np, nt))
                
#                 end

#             # end

#         end

#     end

    # xdivs = tissue.xa:tissue.dc:tissue.Xc
    # TGs = 100:8.75:373.75

    # for x_div ∈ xdivs
    #     for TG ∈ TGs
    #         for seed = 1:100

    #             np, nt = copy(phase), copy(tissue)

    #             nt.TG = TG
    #             nt.x_div = x_div

    #             push!(param_space, (seed, np, nt))

    #         end
    #     end
    # end


end

# println("simulating tissue simulation...")

# result = pmap(x->delay_tissue_simulation(x...), param_space; on_error=identity)

# println("saving tissue simulation...")

# # Now save result to .tsv
# open("data/20230616_motility_delay_tissue_simulation_changing_delay.tsv", "w") do io
#     write(io, "20230616 Cells Added Randomly, Changing motility profile, changing Xv and h")
#     write(io, "\n" * "Ns" * "\t" * "Zs" * "\t" * "freqs" * "\t" * "(phase, tissue)")
#     for line in result
#         if isa(line, Exception)
#             write(io, "\n" * string(line) * "\t" * string(line) * "\t" * string(line) * "\t" * string(line))
#         else
#             write(io, "\n" * string(line[1]) * "\t" * string(line[2]) * "\t" * string(line[3]) * "\t" * string(line[4]))
#         end
#     end
# end

# println("simulating posterior simulation...")

# result = pmap(x->delay_posterior_simulation(x...), param_space; on_error=identity)

# println("saving posterior simulation...")

# # Now save result to .tsv
# open("data/20230616_motility_delay_posterior_simulation_changing_delay.tsv", "w") do io
#     write(io, "20230616 Cells added posterior-dorsally, and randomly in the PSM, Changing motility profile, changing Xv and h")
#     write(io, "\n" * "Ns" * "\t" * "Zs" * "\t" * "freqs" * "\t" * "(phase, tissue)")
#     for line in result
#         if isa(line, Exception)
#             write(io, "\n" * string(line) * "\t" * string(line) * "\t" * string(line) * "\t" * string(line))
#         else
#             write(io, "\n" * string(line[1]) * "\t" * string(line[2]) * "\t" * string(line[3]) * "\t" * string(line[4]))
#         end
#     end
# end

# println("starting posterior lateral simulation...")

# result = pmap(x->delay_posterior_lateral_simulation(x...), param_space; on_error=identity)

# println("saving posterior lateral simulation...")

# # Now save result to .tsv
# open("data/20230616_motility_delay_posterior_lateral_changing_delay.tsv", "w") do io
#     write(io, "20230616 Cells added lateral-ventrally and posterior-dorsally, Changing motility profile, changing Xv and h")
#     write(io, "\n" * "Ns" * "\t" * "Zs" * "\t" * "freqs" * "\t" * "(phase, tissue)")
#     for line in result
#         if isa(line, Exception)
#             write(io, "\n" * string(line) * "\t" * string(line) * "\t" * string(line) * "\t" * string(line))
#         else
#             write(io, "\n" * string(line[1]) * "\t" * string(line[2]) * "\t" * string(line[3]) * "\t" * string(line[4]))
#         end
#     end
# end

# println("Simulating mitosing simulation, param space 1/2...")

# result = pmap(x->mitosing_simulation(x...), param_space; on_error=identity)

# println("Saving mitosing simulation, param space 1/2...")

# # Now save result to .tsv
# open("data/20230821_mitosing_simulation_changing_TG.tsv", "w") do io
#     write(io, "20230821 Cells added randomly, Changing TG and x_div")
#     write(io, "\n" * "Ns" * "\t" * "Zs" * "\t" * "freqs" * "\t" * "(phase, tissue)")
#     for line in result
#         if isa(line, Exception)
#             write(io, "\n" * string(line) * "\t" * string(line) * "\t" * string(line) * "\t" * string(line))
#         else
#             write(io, "\n" * string(line[1]) * "\t" * string(line[2]) * "\t" * string(line[3]) * "\t" * string(line[4]))
#         end
#     end
# end

begin
    
    param_space = Tuple{Int64, PhaseParams, TissueParams}[]

    # xdivs = tissue.xa:tissue.dc:tissue.Xc
    TMs = 0:4.48:44.8
    ωs = [2 * pi / TM for TM in 4.48:4.48:89.6]

    for TM ∈ TMs
        for ω ∈ ωs
            for seed = 1:100

                np, nt = copy(phase), copy(tissue)

                nt.TM = TM
                nt.TG = 187.5 - TM
                np.ωₒ = ω

                push!(param_space, (seed, np, nt))

            end
        end
    end

end

# Now simulate for each parameter set

println("Simulating mitosing simulation, param space 2/2...")

result = pmap(x->mitosing_simulation(x...), param_space; on_error=identity)

println("Saving mitosing simulation, param space 2/2...")

# Now save result to .tsv
open("data/20230829_mitosing_simulation_changing_TM_and_freq.tsv", "w") do io
    write(io, "20230829 Cells added randomly, Changing TM and ω₀")
    write(io, "\n" * "Ns" * "\t" * "Zs" * "\t" * "freqs" * "\t" * "(phase, tissue)")
    for line in result
        if isa(line, Exception)
            write(io, "\n" * string(line) * "\t" * string(line) * "\t" * string(line) * "\t" * string(line))
        else
            write(io, "\n" * string(line[1]) * "\t" * string(line[2]) * "\t" * string(line[3]) * "\t" * string(line[4]))
        end
    end
end

println("Completed.")
