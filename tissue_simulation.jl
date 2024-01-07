# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Purpose of this file is to simulate, in parallel, the computational experiments contained in experiments.jl
# Each simulation is associated with a unique seed to ensure reproducibility. This association is specified
# by the large function tissue_simulation.

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

#Function executing a given simulation and analysis for a given seed.
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
        2101,
        )

    # specify the time range
    t_min, t_max, dt = 0.0, 999.0, 0.01

    # shift the tissue 25 from the origin to allow 'escaped' cells to be simulated on the lattice
    tissue.xa += 25.
    tissue.Xc += 25.
    tissue.Yc += 25.
    tissue.Zc += 25.
    tissue.Lx += 25.
    tissue.Ly += 25.
    tissue.Lz += 25.
    
    if seed < 100

        sim_type = "tissue_simulation"

        tracks = frequency_tracking_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_tissue_simulation_randomnewcells_seed$(seed).csv"

    elseif seed < 200

        sim_type = "adding_zero_cells_simulation"

        tracks = frequency_tracking_adding_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_adding_zero_cells_simulation_seed$(seed).csv"

    elseif seed < 300

        sim_type = "adding_posterior_zero_cells_simulation"

        tracks = frequency_tracking_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_adding_posterior_zero_cells_simulation_seed$(seed).csv"

    elseif seed < 400

        sim_type = "adding_posterior_lateral_zero_cells_simulation"

        tracks = frequency_tracking_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_adding_posterior_lateral_zero_cells_simulation_seed$(seed).csv"

    elseif seed < 500

        sim_type = "mitosing_simulation"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed).csv"

    elseif seed < 600

        tissue.TM = 0
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_TM$(tissue.TM).csv"

    elseif seed < 700

        tissue.TM = 4.48
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_TM$(tissue.TM).csv"

    elseif seed < 800

        tissue.TM = 2. * 4.48
        tissue.TG = 187.5 - tissue.TM
        
        sim_type = "mitosing_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_TM$(tissue.TM).csv"

    elseif seed < 900

        tissue.TM = 3. * 4.48
        tissue.TG = 187.5 - tissue.TM
        
        sim_type = "mitosing_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_TM$(tissue.TM).csv"

    elseif seed < 1000

        tissue.TM = 4. * 4.48
        tissue.TG = 187.5 - tissue.TM
        
        sim_type = "mitosing_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_TM$(tissue.TM).csv"

    elseif seed < 1100

        tissue.TM = 5. * 4.48
        tissue.TG = 187.5 - tissue.TM
        
        sim_type = "mitosing_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_TM$(tissue.TM).csv"

    elseif seed < 1200

        sim_type = "mitosing_no_coupling_simulation"

        tracks = frequency_tracking_mitosing_no_coupling_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_mitosing_no_coupling_simulation_seed$(seed).csv"

    elseif seed < 1300

        tissue.TM = 6. * 4.48
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_TM$(tissue.TM).csv"

    elseif seed < 1400

        tissue.TM = 7. * 4.48
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_TM$(tissue.TM).csv"

    elseif seed < 1500

        tissue.TM = 8. * 4.48
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_TM$(tissue.TM).csv"

    elseif seed < 1600

        tissue.TM = 9. * 4.48
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_TM$(tissue.TM).csv"

    elseif seed < 1700

        tissue.TM = 10. * 4.48
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_TM$(tissue.TM).csv"

    elseif seed < 1800

        tissue.TM = 0
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_no_coupling_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_no_coupling_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_mitosing_no_coupling_simulation_seed$(seed)_TM$(tissue.TM).csv"      
    
    elseif seed < 1900

        tissue.TM = 4.48
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_no_coupling_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_no_coupling_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_no_coupling_simulation_seed$(seed)_TM$(tissue.TM).csv"      
    
    elseif seed < 2000

        tissue.TM = 4.48 * 2
        tissue.TG = 187.5 - tissue.TM
 
        sim_type = "mitosing_no_coupling_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_no_coupling_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_no_coupling_simulation_seed$(seed)_TM$(tissue.TM).csv"      
    
    elseif seed < 2100

        tissue.TM = 4.48 * 3
        tissue.TG = 187.5 - tissue.TM
 
        sim_type = "mitosing_no_coupling_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_no_coupling_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_no_coupling_simulation_seed$(seed)_TM$(tissue.TM).csv"      
    
    elseif seed < 2200

        tissue.TM = 4.48 * 4
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_no_coupling_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_no_coupling_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_no_coupling_simulation_seed$(seed)_TM$(tissue.TM).csv"      
    
    elseif seed < 2300

        tissue.TM = 4.48 * 5
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_no_coupling_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_no_coupling_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_no_coupling_simulation_seed$(seed)_TM$(tissue.TM).csv"      
    
    elseif seed < 2400

        tissue.TM = 4.48 * 6
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_no_coupling_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_no_coupling_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_no_coupling_simulation_seed$(seed)_TM$(tissue.TM).csv"      
    
    elseif seed < 2500

        tissue.TM = 4.48 * 7
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_no_coupling_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_no_coupling_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_no_coupling_simulation_seed$(seed)_TM$(tissue.TM).csv"      
    
    elseif seed < 2600

        tissue.TM = 4.48 * 8
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_no_coupling_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_no_coupling_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_no_coupling_simulation_seed$(seed)_TM$(tissue.TM).csv"      
    
    elseif seed < 2700

        tissue.TM = 4.48 * 9
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_no_coupling_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_no_coupling_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_no_coupling_simulation_seed$(seed)_TM$(tissue.TM).csv"      
    
    elseif seed < 2800

        tissue.TM = 4.48 * 10
        tissue.TG = 187.5 - tissue.TM

        sim_type = "mitosing_no_coupling_simulation_TM$(tissue.TM)"

        tracks = frequency_tracking_mitosing_no_coupling_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_no_coupling_simulation_seed$(seed)_TM$(tissue.TM).csv"      
    
    elseif seed < 2900

        tissue.ρ₀ = 0.003

        tissue.N_cells = floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R) + 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa))
        tissue.total_cells = tissue.N_cells

        sim_type = "tissue_simulation_changing_density_rho$(tissue.ρ₀)"

        tracks = frequency_tracking_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_tissue_simulation_randomnewcells_seed$(seed)_rho$(tissue.ρ₀).csv"

    elseif seed < 3000

        tissue.ρ₀ = 0.00075

        tissue.N_cells = floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R) + 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa))
        tissue.total_cells = tissue.N_cells

        sim_type = "tissue_simulation_changing_density_rho$(tissue.ρ₀)"

        tracks = frequency_tracking_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_tissue_simulation_randomnewcells_seed$(seed)_rho$(tissue.ρ₀).csv"

    elseif seed < 3100

        tissue.ρ₀ = 0.003

        tissue.N_cells = floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R) + 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa))
        tissue.total_cells = tissue.N_cells

        sim_type = "adding_posterior_zero_cells_simulation_changing_density_rho$(tissue.ρ₀)"

        tracks = frequency_tracking_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_adding_posterior_zero_cells_simulation_seed$(seed)_rho$(tissue.ρ₀).csv"

    elseif seed < 3200

        tissue.ρ₀ = 0.00075

        tissue.N_cells = floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R) + 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa))
        tissue.total_cells = tissue.N_cells

        sim_type = "adding_posterior_zero_cells_simulation_changing_density_rho$(tissue.ρ₀)"

        tracks = frequency_tracking_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_adding_posterior_zero_cells_simulation_seed$(seed)_rho$(tissue.ρ₀).csv"

    elseif seed < 3300

        tissue.ρ₀ = 0.003

        tissue.N_cells = floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R) + 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa))
        tissue.total_cells = tissue.N_cells

        sim_type = "adding_posterior_lateral_zero_cells_simulation_changing_density_rho$(tissue.ρ₀)"

        tracks = frequency_tracking_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_rho$(tissue.ρ₀).csv"

    elseif seed < 3400

        tissue.ρ₀ = 0.00075

        tissue.N_cells = floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R) + 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa))
        tissue.total_cells = tissue.N_cells

        sim_type = "adding_posterior_lateral_zero_cells_simulation_changing_density_rho$(tissue.ρ₀)"

        tracks = frequency_tracking_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_rho$(tissue.ρ₀).csv"

    elseif seed < 3500

        tissue.dm = 1.1
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_simulation_seed$(seed)_dm$(tissue.dm).csv"
    
    elseif seed < 3600

        tissue.dm = 0.65
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_simulation_seed$(seed)_dm$(tissue.dm).csv"
    
    elseif seed < 3700

        tissue.dm = 1.1
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_simulation_seed$(seed)_dm$(tissue.dm).csv"

    elseif seed < 3800

        tissue.dm = 2.2
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_simulation_seed$(seed)_dm$(tissue.dm).csv"
    
    elseif seed < 3900

        tissue.dm = 4.4
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_simulation_seed$(seed)_dm$(tissue.dm).csv"

    elseif seed < 4000

        tissue.dm = 0.65
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_adding_posterior_zero_cells_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_adding_posterior_zero_cells_simulation_seed$(seed)_dm$(tissue.dm).csv"

    elseif seed < 4100

        tissue.dm = 1.1
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_adding_posterior_zero_cells_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_adding_posterior_zero_cells_simulation_seed$(seed)_dm$(tissue.dm).csv"

    elseif seed < 4200

        tissue.dm = 2.2
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_adding_posterior_zero_cells_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_adding_posterior_zero_cells_simulation_seed$(seed)_dm$(tissue.dm).csv"

    elseif seed < 4300

        tissue.dm = 4.4
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_adding_posterior_zero_cells_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_adding_posterior_zero_cells_simulation_seed$(seed)_dm$(tissue.dm).csv"

    elseif seed < 4400

        tissue.dm = 0.65
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_dm$(tissue.dm).csv"

    elseif seed < 4500

        tissue.dm = 1.1
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_dm$(tissue.dm).csv"

    elseif seed < 4600
        
        tissue.dm = 2.2
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_dm$(tissue.dm).csv"

    elseif seed < 4700
    
        tissue.dm = 4.4
        tissue.N_cells = 2 * N_column(tissue) + N_torus(tissue)
        tissue.total_cells = tissue.N_cells
        tissue.lattice_side = tissue.dc + tissue.dm

        sim_type = "increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation_dm$(tissue.dm)"

        tracks = frequency_tracking_increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_dm$(tissue.dm).csv"

    elseif seed < 4800

        t_max = 772.
        # PSM at length and radius as measured by Thomson et al. 2021
        tissue.t_shrink = 253.6
        tissue.xd = 10.

        tissue.uₐ = 0.0 # switch off length shrinkage by setting to 0
        tissue.mᵣ = 0.0 # switch off radius shrinkage by setting to 0
        tissue.md = 0.0 # switch off density increase by setting to 0
        tissue.mdc = 0.0 # switch off cell diameter increase by setting to 0

        tissue.L = tissue.L₀
        tissue.Lx = tissue.L + 25.
        tissue.r = tissue.r₀
        tissue.Xc = tissue.L - tissue.R - tissue.r + tissue.xa
        tissue.γₐ = 1.67 #- tissue.uₐ
        tissue.ρ₀ = tissue.d₀
        tissue.Ly = 2 * tissue.r + 2. * tissue.R + 25.
        tissue.Lz = 2 * tissue.r + 25.
        tissue.lattice_side = tissue.dc₀

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa))
        tissue.N_cells += floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tissue_simulation_shrinking_control"

        tracks = frequency_tracking_dynamic_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_dynamic_tissue_simulation_shrinking_control_randomnewcells_seed$(seed).csv"
    
    elseif seed < 4900

        t_max = 772.
        # PSM at length and radius as measured by Thomson et al. 2021
        tissue.t_shrink = 253.6
        tissue.xd = 10.

        tissue.uₐ = 14.1 # switch off length shrinkage by setting to 0
        tissue.mᵣ = 1.1625 # switch off radius shrinkage by setting to 0
        tissue.md = 0.000121 # switch off density increase by setting to 0
        tissue.mdc = 0.2 # switch off cell diameter increase by setting to 0

        tissue.L = tissue.L₀
        tissue.Lx = tissue.L + 25.
        tissue.r = tissue.r₀
        tissue.Xc = tissue.L - tissue.R - tissue.r + tissue.xa
        tissue.γₐ = 1.67 #- tissue.uₐ
        tissue.ρ₀ = tissue.d₀
        tissue.Ly = 2 * tissue.r + 2. * tissue.R + 25.
        tissue.Lz = 2 * tissue.r + 25.
        tissue.lattice_side = tissue.dc₀

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa))
        tissue.N_cells += floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tissue_simulation_shrinking"

        tracks = frequency_tracking_dynamic_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20221208_frequency_tracking_dynamic_tissue_simulation_shrinking_randomnewcells_seed$(seed).csv"
    
    elseif seed < 5000

        tissue.x_div = tissue.xa
        
        sim_type = "mitosing_simulation_x_div$(tissue.x_div)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_x_div$(tissue.x_div).csv"

    elseif seed < 5100

        tissue.x_div = (tissue.xa + tissue.Xc) / 2.
        
        sim_type = "mitosing_simulation_x_div$(tissue.x_div)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_x_div$(tissue.x_div).csv"

    elseif seed < 5200

        tissue.x_div = tissue.Xc
        
        sim_type = "mitosing_simulation_x_div$(tissue.x_div)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20221208_frequency_tracking_mitosing_simulation_seed$(seed)_x_div$(tissue.x_div).csv"

    elseif seed < 5300

        t_max = 772.
        # PSM at length and radius as measured by Thomson et al. 2021
        tissue.t_shrink = 253.6
        tissue.xd = 10.

        tissue.uₐ = 14.1 # switch off length shrinkage by setting to 0
        tissue.mᵣ = 1.1625 # switch off radius shrinkage by setting to 0
        tissue.md = 0.000121 # switch off density increase by setting to 0
        tissue.mdc = 0.2 # switch off cell diameter increase by setting to 0

        tissue.L = tissue.L₀
        tissue.Lx = tissue.L + 25.
        tissue.r = tissue.r₀
        tissue.Xc = tissue.L - tissue.R - tissue.r + tissue.xa
        tissue.γₐ = 1.67 #- tissue.uₐ
        tissue.ρ₀ = tissue.d₀
        tissue.Ly = 2 * tissue.r + 2. * tissue.R + 25.
        tissue.Lz = 2 * tissue.r + 25.
        tissue.lattice_side = tissue.dc₀

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa))
        tissue.N_cells += floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tissue_simulation_shrinking_constantavdection"

        tracks = frequency_tracking_dynamic_nt_simulation_constantadvection(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20230116_frequency_tracking_dynamic_tissue_simulation_constantadvection_shrinking_randomnewcells_seed$(seed).csv"

    elseif seed < 5400

        sim_type = "tissue_simulation_with_delay"

        tracks = tissue_simulation_with_delay(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20230330_tissue_simulation_with_delay_randomnewcells_seed$(seed)_steep_gradient.csv"

    elseif seed < 5500

        sim_type = "adding_zero_cells_simulation_with_delay"

        tracks = adding_zero_cells_simulation_with_delay(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231127_adding_zero_cells_simulation_with_delay_seed$(seed).csv"

    elseif seed < 5600

        sim_type = "adding_posterior_zero_cells_simulation_with_delay"

        tracks = adding_posterior_zero_cells_simulation_with_delay(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231127_tissue_simulation_adding_posterior_zero_cells_simulation_with_delay_seed$(seed).csv"

    elseif seed < 5700

        sim_type = "adding_posterior_lateral_zero_cells_simulation_with_delay"

        tracks = adding_posterior_lateral_zero_cells_simulation_with_delay(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231127_adding_posterior_lateral_zero_cells_simulation_with_delay_seed$(seed).csv"

    elseif seed < 5800

        sim_type = "mitosing_simulation_with_delay_nτ$(phase.nτ)"

        tracks = mitosing_simulation_with_delay(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231127_mitosing_simulation_with_delay_seed$(seed)_nτ$(phase.nτ).csv"

    elseif seed < 5900

        tissue.μb = 100.

        seed = 5800

        t_max = 500.

        sim_type = "mitosing_simulation_varying_spatially"

        tracks = frequency_tracking_mitosing_simulation_varying_spatially(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20230131_mitosing_simulation_varying_spatially_seed$(seed).csv"

    elseif seed < 6000
        # we do the shrinking here
    elseif seed < 6100
        # shrinking
    elseif seed < 6200
        # shrinking
    elseif seed < 6300

        tissue.TG = 162.5

        sim_type = "mitosing_simulation_TG$(tissue.TG)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20230505_frequency_tracking_mitosing_simulation_seed$(seed)_TG$(tissue.TG).csv"

    elseif seed < 6400

        tissue.TG = 182.5

        sim_type = "mitosing_simulation_TG$(tissue.TG)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20230505_frequency_tracking_mitosing_simulation_seed$(seed)_TG$(tissue.TG).csv"

    elseif seed < 6500

        tissue.TG = 152.5

        sim_type = "mitosing_simulation_TG$(tissue.TG)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20230505_frequency_tracking_mitosing_simulation_seed$(seed)_TG$(tissue.TG).csv"

    elseif seed < 6600

        tissue.TG = 192.5

        sim_type = "mitosing_simulation_TG$(tissue.TG)"

        tracks = frequency_tracking_mitosing_simulation(seed, t_min, t_max, dt, tissue, phase)
        
        fname = "data/20230505_frequency_tracking_mitosing_simulation_seed$(seed)_TG$(tissue.TG).csv"

    elseif seed < 6700

        # change nτ
        phase.nτ = 2601

        sim_type = "mitosing_simulation_with_delay_nτ$(phase.nτ)"

        tracks = mitosing_simulation_with_delay(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231127_mitosing_simulation_with_delay_seed$(seed)_nτ$(phase.nτ).csv"

    elseif seed < 6800

        # change nτ
        phase.nτ = 1

        sim_type = "mitosing_simulation_with_delay_nτ$(phase.nτ)"

        tracks = mitosing_simulation_with_delay(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231127_mitosing_simulation_with_delay_seed$(seed)_nτ$(phase.nτ).csv"

    elseif seed < 6900

        Xc = 200

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_tissue_simulation_random_cells_Xc$(tissue.Xc)"

        tracks, added_cells = addition_tracking_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_tissue_simulation_seed$(seed)_Xc$(tissue.Xc).csv"

    elseif seed < 7000

        Xc = 300

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_tissue_simulation_random_cells_Xc$(tissue.Xc)"

        tracks, added_cells = addition_tracking_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_tissue_simulation_seed$(seed)_Xc$(tissue.Xc).csv"

    elseif seed < 7100

        Xc = 400

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_tissue_simulation_random_cells_Xc$(tissue.Xc)"

        tracks, added_cells = addition_tracking_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_tissue_simulation_seed$(seed)_Xc$(tissue.Xc).csv"

    elseif seed < 7200

        Xc = 200

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_zero_cells_simulation_Xc$(tissue.Xc)"

        tracks, added_cells = addition_tracking_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_adding_posterior_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc).csv"

    elseif seed < 7300

        Xc = 300

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_zero_cells_simulation_Xc$(tissue.Xc)"

        tracks, added_cells = addition_tracking_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_adding_posterior_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc).csv"


    elseif seed < 7400

        Xc = 400

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_zero_cells_simulation_Xc$(tissue.Xc)"

        tracks, added_cells = addition_tracking_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_adding_posterior_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc).csv"


    elseif seed < 7500

        Xc = 200

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc$(tissue.Xc)"

        tracks, added_cells = tracking_addition_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_tracking_addition_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc).csv"

    elseif seed < 7600

        Xc = 300

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc$(tissue.Xc)"

        tracks, added_cells = tracking_addition_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_tracking_addition_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc).csv"


    elseif seed < 7700

        Xc = 400

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc$(tissue.Xc)"

        tracks, added_cells = tracking_addition_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_tracking_addition_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc).csv"

    elseif seed < 7800

        tissue.γₐ = 3.0

        sim_type = "tracking_addition_tissue_simulation_random_cells_va$(tissue.γₐ)"

        tracks, added_cells = addition_tracking_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_tissue_simulation_seed$(seed)_va$(tissue.γₐ).csv"

    elseif seed < 7900

        tissue.γₐ = 3.0

        sim_type = "tracking_addition_adding_posterior_zero_cells_simulation_va$(tissue.γₐ)"

        tracks, added_cells = addition_tracking_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_adding_posterior_zero_cells_simulation_seed$(seed)_n_va$(tissue.γₐ).csv"

    elseif seed < 8000

        tissue.γₐ = 3.0

        sim_type = "tracking_addition_adding_posterior_lateral_zero_cells_simulation_va$(tissue.γₐ)"

        tracks, added_cells = tracking_addition_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_tracking_addition_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_va$(tissue.γₐ).csv"

    # now do simulations for each different length of the PSM but keep the length of tissue that cells can be added into the save_simulations

    elseif seed < 8100

        Xc = 200

        tissue.Xc = Xc + tissue.xa

        # set position at which cells can be added into the tissue
        tissue.xd = Xc - 100.0 # note that the function add_cells! and equivalent account for varying xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_tissue_simulation_random_cells_Xc$(tissue.Xc)_xd$(tissue.xd)"

        tracks, added_cells = addition_tracking_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_tissue_simulation_seed$(seed)_Xc$(tissue.Xc)_xd$(tissue.xd).csv"

    elseif seed < 8200

        Xc = 300

        tissue.Xc = Xc + tissue.xa

        # set position at which cells can be added into the tissue
        tissue.xd = Xc - 100.0 # note that the function add_cells! and equivalent account for varying xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_tissue_simulation_random_cells_Xc$(tissue.Xc)_xd$(tissue.xd)"

        tracks, added_cells = addition_tracking_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_tissue_simulation_seed$(seed)_Xc$(tissue.Xc)_xd$(tissue.xd).csv"

    elseif seed < 8300

        Xc = 400

        tissue.Xc = Xc + tissue.xa

        # set position at which cells can be added into the tissue
        tissue.xd = Xc - 100.0 # note that the function add_cells! and equivalent account for varying xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_tissue_simulation_random_cells_Xc$(tissue.Xc)_xd$(tissue.xd)"

        tracks, added_cells = addition_tracking_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_tissue_simulation_seed$(seed)_Xc$(tissue.Xc)_xd$(tissue.xd).csv"

    elseif seed < 8400

        Xc = 500

        tissue.Xc = Xc + tissue.xa

        # set position at which cells can be added into the tissue
        tissue.xd = Xc - 100.0 # note that the function add_cells! and equivalent account for varying xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_tissue_simulation_random_cells_Xc$(tissue.Xc)_xd$(tissue.xd)"

        tracks, added_cells = addition_tracking_tissue_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_tissue_simulation_seed$(seed)_Xc$(tissue.Xc)_xd$(tissue.xd).csv"

    elseif seed < 8500

        Xc = 200

        tissue.xd = Xc - 100.0

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_zero_cells_simulation_Xc$(tissue.Xc)_xd$(tissue.xd)"

        tracks, added_cells = addition_tracking_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_adding_posterior_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc)_xd$(tissue.xd).csv"

    elseif seed < 8600

        Xc = 300

        tissue.xd = Xc - 100.0

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_zero_cells_simulation_Xc$(tissue.Xc)_xd$(tissue.xd)"

        tracks, added_cells = addition_tracking_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_adding_posterior_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc)_xd$(tissue.xd).csv"

    elseif seed < 8700

        Xc = 400

        tissue.xd = Xc - 100.0

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_zero_cells_simulation_Xc$(tissue.Xc)_xd$(tissue.xd)"

        tracks, added_cells = addition_tracking_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_adding_posterior_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc)_xd$(tissue.xd).csv"

    elseif seed < 8800

        Xc = 500

        tissue.xd = Xc - 100.0

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_zero_cells_simulation_Xc$(tissue.Xc)_xd$(tissue.xd)"

        tracks, added_cells = addition_tracking_adding_posterior_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_addition_tracking_adding_posterior_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc)_xd$(tissue.xd).csv"

    elseif seed < 8900

        Xc = 200

        tissue.xd = Xc - 100.0

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc$(tissue.Xc)_xd$(tissue.xd)"

        tracks, added_cells = tracking_addition_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_tracking_addition_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc)_xd$(tissue.xd).csv"

    elseif seed < 9000

        Xc = 300

        tissue.xd = Xc - 100.0

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc$(tissue.Xc)_xd$(tissue.xd)"

        tracks, added_cells = tracking_addition_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_tracking_addition_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc)_xd$(tissue.xd).csv"

    elseif seed < 9100

        Xc = 400

        tissue.xd = Xc - 100.0

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc$(tissue.Xc)_xd$(tissue.xd)"

        tracks, added_cells = tracking_addition_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_tracking_addition_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc)_xd$(tissue.xd).csv"


    elseif seed < 9200

        Xc = 500

        tissue.xd = Xc - 100.0

        tissue.Xc = Xc + tissue.xa

        tissue.L = tissue.Xc + tissue.R + tissue.r - tissue.xa

        tissue.Lx = tissue.L + 25. + tissue.xa

        tissue.N_cells = 2 * floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) + floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R)
        tissue.total_cells = tissue.N_cells

        sim_type = "tracking_addition_adding_posterior_lateral_zero_cells_simulation_Xc$(tissue.Xc)_xd$(tissue.xd)"

        tracks, added_cells = tracking_addition_adding_posterior_lateral_zero_cells_simulation(seed, t_min, t_max, dt, tissue, phase)

        fname = "data/20231206_tracking_addition_adding_posterior_lateral_zero_cells_simulation_seed$(seed)_Xc$(tissue.Xc)_xd$(tissue.xd).csv"


    elseif seed < 9300

        # shrinking sims

    elseif seed < 9400

    elseif seed < 9500

    elseif seed < 9600

    end

    # save the first 10 simulations
    if seed % 100 < 10

        #save the simulation
        tracking_freq_arr_to_csv(
            tracks, # simulation
            t_min, 
            t_max, 
            fname
            )

    end

    data = thin_domain(tissue.xa, tracks[end], 1, tissue.dc, tissue.r, tissue.R)
    r = phase_order(data[:, 6])
    # exclude the cells that are in M phase if we're in a mitosing simulation
    freq = mean_frequency(data, tissue; excluding_mitosing=occursin("mitosing", sim_type))

    if occursin("shrinking", sim_type)

        anterior = xas(tracks, tissue)

        rs = anterior_synchrony(tracks, 25., 60., tissue.dc, anterior)

        freq = anterior_frequency(tracks, 25., 60., tissue.dc, anterior, tissue)

        return sim_type, seed, rs, freq

    elseif occursin("tracking_addition", sim_type)

        return sim_type, seed, r, freq, added_cells

    end

    return sim_type, seed, r, freq

end


used_seeds = [x for x in 0:6799]
seeds = [x for x in 6800:9199 if !(x in used_seeds)]
# append!(seeds, [x for x in 6600:6799])

@everywhere println("Simulating...")

# simulate in parallel

result = pmap(seed -> tissue_simulation(seed), seeds)

@everywhere println("Simulation complete. Saving...")
# Save output in big .csv
result[end]
output_fname = "data/20231206_changinglength_control_addition_tracking_simulations.tsv"

# write_output(result, output_fname)
write_output_with_cell_addition(result, output_fname) # just for when running simulations with tracking addition

# print to check for errors in output
println("\n\n\n")
println(result)
println("\n\n\n")
println("Done!")
