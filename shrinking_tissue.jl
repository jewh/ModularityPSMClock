# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Code for performing compaction-extension simulations on a HPC cluster using Slurm

println("Using Julia $(VERSION).")

# code for parallelised slurm job
using Distributed, ClusterManagers

available_workers = parse(Int, ENV["SLURM_NTASKS"])

addprocs_slurm(available_workers)

@everywhere println("Loading code...")

@everywhere begin
    include("experiments.jl") # code for simulating
    include("save_simulations.jl") # load code for saving data as .csv
end

@everywhere function tissue_simulation(seed)

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

    # specify the time range
    # simulate until 690 mins as beyond that we observe defects in PSM structure - gaps in the tissue form
    t_min, t_max, dt = 0.0, 690., 0.01

    # shift the tissue 25 from the origin to allow 'escaped' cells to be simulated on the lattice
    tissue.xa += 25.
    tissue.Xc += 25.
    tissue.Yc += 25.
    tissue.Zc += 25.

    tissue.r = tissue.r₀
    tissue.ρ₀ = tissue.d₀
    tissue.L = tissue.L₀
    tissue.xa = xa(0, tissue)
    tissue.dc = tissue.dc₀

    tissue.Lx += 25.
    tissue.Ly = 2 * tissue.r + 2. * tissue.R + 25.
    tissue.Lz = 2 * tissue.r + 25.
    tissue.lattice_side = tissue.dc₀

    tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa))
    tissue.N_cells += floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
    tissue.total_cells = tissue.N_cells

    tissue.t_shrink = 253.6
    tissue.t_stop = t_max # time at which stop shrinking the cells in the PSM


    tissue.uₐ = 14.1 # switch off length shrinkage by setting to 0
    tissue.mᵣ = 1.1625 # switch off radius shrinkage by setting to 0
    tissue.md = 0.000121 # switch off density increase by setting to 0
    tissue.mdc = 0.2

    if seed < 4800
        # PSM at length and radius as measured by Thomson et al. 2021

        tissue.dc₀ = 11.0
        tissue.dc = 11.0
        tissue.uₐ = 0.0 # switch off length shrinkage by setting to 0
        tissue.mᵣ = 0.0 # switch off radius shrinkage by setting to 0
        tissue.md = 0.0 # switch off density increase by setting to 0
        tissue.mdc = 0.0

        sim_type = "tissue_simulation_shrinking_control"

        tracks, n_additions = frequency_tracking_dynamic_tissue_simulation_tracking_addition(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20230331_frequency_tracking_dynamic_tissue_simulation_shrinking_control_randomnewcells_seed$(seed).csv"
    
    elseif seed < 4900
        # PSM at length and radius as measured by Thomson et al. 2021
        
        tissue.uₐ = 14.1 # switch off length shrinkage by setting to 0
        tissue.mᵣ = 1.1625 # switch off radius shrinkage by setting to 0
        tissue.md = 0.000121 # switch off density increase by setting to 0
        tissue.mdc = 0.2

        sim_type = "tissue_simulation_shrinking"

        tracks, n_additions = frequency_tracking_dynamic_tissue_simulation_tracking_addition(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20230331_frequency_tracking_dynamic_tissue_simulation_shrinking_randomnewcells_seed$(seed).csv"

    elseif seed < 5300

        # PSM at length and radius as measured by Thomson et al. 2021

        tissue.uₐ = 14.1 # switch off length shrinkage by setting to 0
        tissue.mᵣ = 1.1625 # switch off radius shrinkage by setting to 0
        tissue.md = 0.000121 # switch off density increase by setting to 0
        tissue.mdc = 0.2

        tissue.γₐ = 1.67 #- tissue.uₐ

        sim_type = "tissue_simulation_shrinking_constantavdection"

        tracks, n_additions = frequency_tracking_dynamic_tissue_simulation_constantadvection_tracking_addition(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20230331_frequency_tracking_dynamic_tissue_simulation_constantadvection_shrinking_randomnewcells_seed$(seed).csv"

    elseif seed < 6000

        # PSM at length and radius as measured by Thomson et al. 2021

        tissue.uₐ = 14.1 # switch off length shrinkage by setting to 0
        tissue.mᵣ = 1.1625 # switch off radius shrinkage by setting to 0
        tissue.md = 0.000121 # switch off density increase by setting to 0
        tissue.mdc = 0.2

        tissue.γₐ = 1.67 #- tissue.uₐ

        # stop shrinking the cells at 26ss
        tissue.t_stop = 536.4

        sim_type = "tissue_simulation_shrinking_constantavdection_stepwisecellshrinking"

        tracks, n_additions = frequency_tracking_dynamic_tissue_simulation_constantadvection_tracking_addition(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20230331_frequency_tracking_dynamic_tissue_simulation_constantadvection_shrinking_randomnewcells_stepwisecellshrinking_seed$(seed).csv"

    elseif seed < 6100

        tissue.d₀ = 0.001123
        tissue.ρ₀ = tissue.d₀

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa))
        tissue.N_cells += floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        tissue.uₐ = 14.1 # switch off length shrinkage by setting to 0
        tissue.mᵣ = 1.1625 # switch off radius shrinkage by setting to 0
        tissue.md = 0.000121 # switch off density increase by setting to 0
        tissue.mdc = 0.2

        tissue.γₐ = 1.67 #- tissue.uₐ

        sim_type = "tissue_simulation_shrinking_constantavdection_d0$(tissue.d₀)"

        tracks, n_additions = frequency_tracking_dynamic_tissue_simulation_constantadvection_tracking_addition(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20230331_frequency_tracking_dynamic_tissue_simulation_constantadvection_shrinking_randomnewcells__d0$(tissue.d₀)_seed$(seed).csv"

    elseif seed < 6200

        tissue.d₀ = 0.003123
        tissue.ρ₀ = tissue.d₀

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa))
        tissue.N_cells += floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        tissue.uₐ = 14.1 # switch off length shrinkage by setting to 0
        tissue.mᵣ = 1.1625 # switch off radius shrinkage by setting to 0
        tissue.md = 0.000121 # switch off density increase by setting to 0
        tissue.mdc = 0.2

        tissue.γₐ = 1.67 #- tissue.uₐ

        sim_type = "tissue_simulation_shrinking_constantavdection_d0$(tissue.d₀)"

        tracks, n_additions = frequency_tracking_dynamic_tissue_simulation_constantadvection_tracking_addition(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20230331_frequency_tracking_dynamic_tissue_simulation_constantadvection_shrinking_randomnewcells__d0$(tissue.d₀)_seed$(seed).csv"

    elseif 9200 <= seed < 9300
        # keep density constant

        tissue.uₐ = 14.1 # switch off length shrinkage by setting to 0
        tissue.mᵣ = 1.1625 # switch off radius shrinkage by setting to 0
        tissue.md = 0.0 # switch off density increase by setting to 0
        tissue.mdc = 0.2 # changing cell diameter too

        tissue.γₐ = 1.67 #- tissue.uₐ

        sim_type = "tissue_simulation_shrinking_constantavdection_constantdensity"

        tracks, n_additions = frequency_tracking_dynamic_tissue_simulation_constantadvection_tracking_addition(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231216_frequency_tracking_dynamic_tissue_simulation_constantadvection_constantdensity_shrinking_randomnewcells_seed$(seed).csv"

    elseif 9300 <= seed < 9400

        
        tissue.uₐ = 14.1 # switch off length shrinkage by setting to 0
        tissue.mᵣ = 1.1625 # switch off radius shrinkage by setting to 0
        tissue.md = 0.0 # switch off density increase by setting to 0
        tissue.mdc = 0.0 # switch of decrease in cell diameter

        tissue.γₐ = 1.67 #- tissue.uₐ

        sim_type = "tissue_simulation_shrinking_constantavdection_constantdensity_constant_dc"

        tracks, n_additions = frequency_tracking_dynamic_tissue_simulation_constantadvection_tracking_addition(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231216_frequency_tracking_dynamic_tissue_simulation_constantadvection__constantdensity_constant_dc_shrinking_randomnewcells_seed$(seed).csv"

    end

    if seed % 100 < 10

        #save the simulation
        tracking_freq_arr_to_csv(
            tracks, # simulation
            t_min, 
            t_max, # doesn;t matter anymore 
            fname
            )

    end

    anterior = xas(tracks, tissue)

    rs = anterior_synchrony(tracks, 25., 60., tissue.dc₀, anterior)

    freq = anterior_frequency(tracks, 25., 60., tissue.dc₀, anterior, tissue)

    densities = anterior_density(tracks, 25., 60., tissue.dc₀, anterior)

    number_of_cells = cells_in_psm_in_time(tracks, anterior)

    return sim_type, seed, rs, freq, densities, number_of_cells, anterior, n_additions

end

seeds = Int64[] # 🌱

# for n = 4700:4899

#     push!(seeds, n)

# end

# for n = 5200:5299

#     push!(seeds, n)

# end

# for n = 5900:5999

#     push!(seeds, n)

# end

# for n = 6000:6199
    
#     push!(seeds, n)

# end

for n = 9200:9399

    push!(seeds, n)

end


@everywhere println("Simulating...")

# simulate in parallel

result = pmap(seed -> tissue_simulation(seed), seeds; on_error=identity)

@everywhere println("Simulation complete. Saving...")
# Save output in big .csv

output_fname = "data/20231216_stochastic_shrinking.tsv"

write_output_more_columns(result, output_fname)

# print to check for errors in output
println("\n\n\n")
println(result)
println("\n\n\n")
println("Done!")

