# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# In this file there are functions corresponding to computational 'experiments' for different scenarios of PSM morphogenesis, that take parameter objects as input and a seed for reproducibility. This is where we numerically solve the equations for phase and cell position.

include("cell_movements.jl") # load code for cell movements
include("kuramoto.jl") # load code for phase oscillator dynamics
include("cell_addition.jl") # load code for adding cells
include("cell_cycle.jl") # load code for simulating effect of mitosis.
include("density_gradient.jl") # code specifying the density of the tissue
include("dynamic_tissue_geometry.jl") # code for dynamically modifying the tissue's geometry
include("delay.jl") # functions for solving a DDE kuramoto model of clock phase
include("utils.jl") # misc housekeeping functions


"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Uses FCell structs to track the frequency in real time.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_tissue_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)

        euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        add_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Uses FCell structs to track the frequency in real time. Adding zero cells.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_adding_zero_cells_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)
    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
        
    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        
        # # call garbage collection at random to reduce memory burden
        # thr=0.01
        # rand(Uniform(0,1)) < thr && GC.gc()

        euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        add_zero_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

# function frequency_tracking_static_tissue_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)
#     # Get the initial condition
#     cells = initial_FCells_horseshoe_PSM(seed, phase, tissue)

#     # Now 'relax' the tissue for 10 mins
#     for t = 0:dt:10
#         initial_relaxation!(cells, phase, tissue, dt, boundary_horseshoe_psm)
#     end

#     # Initialise output array
#     tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))
#     # Simulate
#     k = 1
#     for t = 1:length(t_min:dt:t_max)

#         new_cells = copy(cells)

#         euler_step!(cells, new_cells, dt, phase, tissue) # simulate gene expression on cell positions at previous time point
#         # save after every minute, and at the initial condition (t = 0)
#         if (t - 1) % (1 / dt) == 0
#             M = Matrix{Float64}(undef, tissue.total_cells, 8) 
#             for i = 1:tissue.total_cells
#                 if i <= tissue.N_cells # this value updates throughout the simulation
#                     M[i, 1:3] .= cells[i].x
#                     M[i, 4] = t_min + dt * (t - 1)
#                     M[i, 5] = cells[i].trackid 
#                     M[i, 6] = cells[i].θ
#                     M[i, 7] = cells[i].dθdt
#                     M[i, 8] = cells[i].τ
#                 elseif tissue.N_cells < i <= tissue.total_cells
#                     M[i, 1:3] .= old_cells[i - tissue.N_cells].x
#                     M[i, 4] = t_min + dt * (t - 1)
#                     M[i, 5] = old_cells[i - tissue.N_cells].trackid 
#                     M[i, 6] = old_cells[i - tissue.N_cells].θ
#                     M[i, 7] = old_cells[i - tissue.N_cells].dθdt
#                     M[i, 8] = old_cells[i - tissue.N_cells].τ
#                 end
#             end
#             tracks[k] = M
#             k += 1
#         end

#         cells .= new_cells
#     end
#     tracks # return subset
# end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021. Here, 
new cells are added to the dorsum of the toroid tailbud if the tissue density falls below requisite threshold,
and phase θ of these new cells is set to θ = 0. To maintain density, also add cells with θ ∈ [0, 2π) in the LHS and RHS PSM.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_adding_posterior_zero_cells_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)
    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
        
    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        
        # # call garbage collection at random to reduce memory burden
        # thr=0.01
        # rand(Uniform(0,1)) < thr && GC.gc()

        euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        add_posterior_zero_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

# """
# Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021. Here, 
# new cells are added to the dorsum of the toroid tailbud if the tissue density falls below requisite threshold,
# and phase θ of these new cells is set to θ = 0.

# seed - seed for reproducing the initial condition, and simulation
# t_min - start time (mins)
# t_max - end time (mins)
# dt - time step (mins)
# tissue - parameters defining tissue geometry and cell movements
# phase - parameters specifying the kuramoto model for phase oscillators

# Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
# corresponds to a time point of the simulation, and the elements of n x 8 matrix 
# M correspond to:

# M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)
# M[i, 4] => t (timepoint of the simulation)
# M[i, 5] => trackid (ID of the cell i)
# M[i, 6] => θ (segmentation clock phase)
# M[i, 7] => dθᵢ/dt (frequency of phase oscillations)
# M[i, 8] => τ (cell cycle phase)
# """
# function frequency_tracking_adding_posterior_zero_cells_simulation(seed::Real, t_min, t_max, dt, ϕ, tissue::TissueParams, phase::PhaseParams)
#     @assert 0 <= ϕ <= π/2
#     # Get the initial condition
#     cells = initial_FCells_horseshoe_PSM(seed, phase, tissue)

#     # Now 'relax' the tissue for 10 mins, so that 
#     for t = 0:dt:10
#         initial_relaxation!(cells, phase, tissue, dt, boundary_horseshoe_psm)
#     end

#     # Create a waste array for cells anterior to the PSM
#     old_cells = Vector{FCell}(undef, (2 * 883 + 555) * 150);

#     # Initialise output array
#     tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))
#     # Simulate
#     k = 1
#     for t = 1:length(t_min:dt:t_max)
#         euler_step!(cells, dt, phase, tissue) # simulate gene expression on cell positions at previous time point
#         cell_movements!(cells, phase, tissue, dt, boundary_horseshoe_psm) # move cells in the psm
#         anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
#         remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells
#         add_posterior_zero_cells!(cells, ϕ, phase, tissue) # add new cells to replenish psm
#         # save after every minute, and at the initial condition (t = 0)
#         if (t - 1) % (1 / dt) == 0
#             M = Matrix{Float64}(undef, tissue.total_cells, 8) 
#             for i = 1:tissue.total_cells
#                 if i <= tissue.N_cells # this value updates throughout the simulation
#                     M[i, 1:3] .= cells[i].x
#                     M[i, 4] = t_min + dt * (t - 1)
#                     M[i, 5] = cells[i].trackid 
#                     M[i, 6] = cells[i].θ
#                     M[i, 7] = cells[i].dθdt
#                     M[i, 8] = cells[i].τ
#                 elseif tissue.N_cells < i <= tissue.total_cells
#                     M[i, 1:3] .= old_cells[i - tissue.N_cells].x
#                     M[i, 4] = t_min + dt * (t - 1)
#                     M[i, 5] = old_cells[i - tissue.N_cells].trackid 
#                     M[i, 6] = old_cells[i - tissue.N_cells].θ
#                     M[i, 7] = old_cells[i - tissue.N_cells].dθdt
#                     M[i, 8] = old_cells[i - tissue.N_cells].τ
#                 end
#             end
#             tracks[k] = M
#             k += 1
#         end
#     end
#     tracks
# end


"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021. Here, 
new cells are added to the dorsum of the toroid tailbud if the tissue density falls below requisite threshold,
and phase θ of these new cells is set to θ ~ [0, 2π).\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_adding_posterior_random_cells_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)
    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
        
    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        
        # # call garbage collection at random to reduce memory burden
        # thr=0.01
        # rand(Uniform(0,1)) < thr && GC.gc()

        euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        add_posterior_random_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

# """
# Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021. Here, 
# new cells are added to the dorsum of the toroid tailbud if the tissue density falls below requisite threshold,
# and phase θ of these new cells is set to θ ~ [0, 2π).

# seed - seed for reproducing the initial condition, and simulation
# t_min - start time (mins)
# t_max - end time (mins)
# dt - time step (mins)
# tissue - parameters defining tissue geometry and cell movements
# phase - parameters specifying the kuramoto model for phase oscillators

# Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
# corresponds to a time point of the simulation, and the elements of n x 8 matrix 
# M correspond to:

# M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)
# M[i, 4] => t (timepoint of the simulation)
# M[i, 5] => trackid (ID of the cell i)
# M[i, 6] => θ (segmentation clock phase)
# M[i, 7] => dθᵢ/dt (frequency of phase oscillations)
# M[i, 8] => τ (cell cycle phase)
# """
# function frequency_tracking_adding_posterior_random_cells_simulation(seed::Real, t_min, t_max, dt, ϕ, tissue::TissueParams, phase::PhaseParams)
#     @assert 0 <= ϕ <= π/2
#     # Get the initial condition
#     cells = initial_FCells_horseshoe_PSM(seed, phase, tissue)

#     # Now 'relax' the tissue for 10 mins, so that 
#     for t = 0:dt:10
#         initial_relaxation!(cells, phase, tissue, dt, boundary_horseshoe_psm)
#     end

#     # Create a waste array for cells anterior to the PSM
#     old_cells = Vector{FCell}(undef, (2 * 883 + 555) * 150);

#     # Initialise output array
#     tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))
#     # Simulate
#     k = 1
#     for t = 1:length(t_min:dt:t_max)
#         euler_step!(cells, dt, phase, tissue) # simulate gene expression on cell positions at previous time point
#         cell_movements!(cells, phase, tissue, dt, boundary_horseshoe_psm) # move cells in the psm
#         anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
#         remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells
#         add_posterior_random_cells!(cells, ϕ, phase, tissue) # add new cells to replenish psm
#         # save after every minute, and at the initial condition (t = 0)
#         if (t - 1) % (1 / dt) == 0
#             M = Matrix{Float64}(undef, tissue.total_cells, 8) 
#             for i = 1:tissue.total_cells
#                 if i <= tissue.N_cells # this value updates throughout the simulation
#                     M[i, 1:3] .= cells[i].x
#                     M[i, 4] = t_min + dt * (t - 1)
#                     M[i, 5] = cells[i].trackid 
#                     M[i, 6] = cells[i].θ
#                     M[i, 7] = cells[i].dθdt
#                     M[i, 8] = cells[i].τ
#                 elseif tissue.N_cells < i <= tissue.total_cells
#                     M[i, 1:3] .= old_cells[i - tissue.N_cells].x
#                     M[i, 4] = t_min + dt * (t - 1)
#                     M[i, 5] = old_cells[i - tissue.N_cells].trackid 
#                     M[i, 6] = old_cells[i - tissue.N_cells].θ
#                     M[i, 7] = old_cells[i - tissue.N_cells].dθdt
#                     M[i, 8] = old_cells[i - tissue.N_cells].τ
#                 end
#             end
#             tracks[k] = M
#             k += 1
#         end
#     end
#     tracks
# end


"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021. Here, 
new cells are added to the dorsum of the toroid tailbud or ventral-lateral sides of the PSM if the tissue density falls below requisite threshold,
and phase θ of these new cells is set to θ = 0.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_adding_posterior_lateral_zero_cells_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)
    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
        
    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        
        # # call garbage collection at random to reduce memory burden
        # thr=0.01
        # rand(Uniform(0,1)) < thr && GC.gc()

        euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        add_posterior_lateral_zero_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Cells increase cell cycle phase τ linearly and divide when τ ≥ TG + TM. If τ ∈ [TG, TG+TM) then dθᵢ/dt = 0.
New cells are added at random positions within the PSM if the tissue density after division falls below requisite threshold,
and phase θ of these new cells is a random value within [0, 2π).\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_mitosing_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)
    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
        
    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        
        # # call garbage collection at random to reduce memory burden
        # thr=0.01
        # rand(Uniform(0,1)) < thr && GC.gc()

        euler_step_with_mitosis!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells 
        mitosis!(cells, dt, tissue)      
        add_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Cells increase cell cycle phase τ linearly and divide when τ ≥ TG + TM. If τ ∈ [TG, TG+TM) then dθᵢ/dt = 0, and no adjacent
cell may couple its clock to cell i. New cells are added at random positions within the PSM if the tissue density after division falls below requisite threshold,
and phase θ of these new cells is a random value within [0, 2π).\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_mitosing_no_coupling_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)
    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
        
    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        
        # # call garbage collection at random to reduce memory burden
        # thr=0.01
        # rand(Uniform(0,1)) < thr && GC.gc()

        euler_step_with_mitosis_no_coupling!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells  
        mitosis!(cells, dt, tissue)     
        add_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Cells increase cell cycle phase τ according to a function which synchronises the τ with θ, and divide when τ ≥ TG + TM. If τ ∈ [TG, TG+TM) then dθᵢ/dt = 0.
New cells are added at random positions within the PSM if the tissue density after division falls below requisite threshold,
and phase θ of these new cells is a random value within [0, 2π).\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_mitosing_simulation_coupled_to_phase(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)
    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
        
    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        
        # # call garbage collection at random to reduce memory burden
        # thr=0.01
        # rand(Uniform(0,1)) < thr && GC.gc()

        euler_step_with_mitosis!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells     
        mitosis_coupled_to_phase!(cells, dt, tissue)  
        add_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Extracellular volume increases towards the posterior of the tissue. Cells added at random.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_increasing_extracellular_volume_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    cells = initial_FCells_horseshoe_PSM_increasing_extracellular_volume(seed, phase, tissue)

    # define lattice

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    for t = 1:dt:10
        initial_relaxation_with_increasing_extracellular_volume!(cells, tissue, dt, boundary_horseshoe_psm, d, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
    end

    old_cells = Vector{FCell}(undef, tissue.total_cells * 10);

    # Initialise output array
    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))
    # Simulate
    k = 1
    for t = 1:length(t_min:dt:t_max)

        euler_step!(cells, dt, phase, tissue, d, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, d, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells
        cleanup_escaped_cells!(cells, tissue) # remove cells that have escaped the PSM
        add_cells_changing_diameter!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

        # save after every minute, and at the initial condition (t = 0)
        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8) 

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M
            k += 1

        end

    end

    tracks # return subset

end


"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Extracellular volume increases towards the posterior of the tissue. new cells are added to the dorsum of the toroid tailbud if the tissue density falls below requisite threshold,
and phase θ of these new cells is set to θ = 0. To maintain density, also add cells with θ ∈ [0, 2π) in the LHS and RHS PSM.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_increasing_extracellular_volume_adding_posterior_zero_cells_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    cells = initial_FCells_horseshoe_PSM_increasing_extracellular_volume(seed, phase, tissue)

    # define lattice

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    for t = 1:dt:10
        initial_relaxation_with_increasing_extracellular_volume!(cells, tissue, dt, boundary_horseshoe_psm, d, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
    end

    old_cells = Vector{FCell}(undef, tissue.total_cells * 10);

    # Initialise output array
    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))
    # Simulate
    k = 1
    for t = 1:length(t_min:dt:t_max)

        euler_step!(cells, dt, phase, tissue, d, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, d, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells
        add_posterior_zero_cells_changing_diameter!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

        # save after every minute, and at the initial condition (t = 0)
        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8) 

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M
            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Extracellular volume increases towards the posterior of the tissue. new cells are added to the dorsum of the toroid tailbud or ventral-lateral sides of the PSM if the tissue density falls below requisite threshold,
and phase θ of these new cells is set to θ = 0.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_increasing_extracellular_volume_adding_posterior_lateral_zero_cells_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    cells = initial_FCells_horseshoe_PSM_increasing_extracellular_volume(seed, phase, tissue)

    # define lattice

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    for t = 1:dt:10
        initial_relaxation_with_increasing_extracellular_volume!(cells, tissue, dt, boundary_horseshoe_psm, d, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
    end

    old_cells = Vector{FCell}(undef, tissue.total_cells * 10);

    # Initialise output array
    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))
    # Simulate
    k = 1
    for t = 1:length(t_min:dt:t_max)

        euler_step!(cells, dt, phase, tissue, d, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, d, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells
        add_posterior_lateral_zero_cells_changing_diameter!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

        # save after every minute, and at the initial condition (t = 0)
        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8) 

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M
            k += 1

        end

    end

    tracks # return subset

end


"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Uses FCell structs to track the frequency in real time.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_resynchronising_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)
    # Get the initial condition
    cells = random_initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
        
    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        
        # # call garbage collection at random to reduce memory burden
        # thr=0.01
        # rand(Uniform(0,1)) < thr && GC.gc()

        euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        add_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Tissue evolves dynamically according to the functions in dynamic_tissue_geometry.jl.
Uses FCell structs to track the frequency in real time.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_dynamic_tissue_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)

        euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        add_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # update geometry
        # length
        tissue.L = Length(t_min + dt * (t - 1), tissue.uₐ, tissue.L₀, tissue.t_shrink)
        
        tissue.xa = xa(t_min + dt * (t - 1), tissue)

        # advection speed
        tissue.γₐ = 1.67 - wavefront_velocity(t_min + dt * (t - 1), tissue.uₐ, tissue.t_shrink)

        # radius
        tissue.r = radius(t_min + dt * (t - 1), tissue.mᵣ, tissue.r₀, tissue.t_shrink)

        # density
        tissue.ρ₀ = density(t_min + dt * (t - 1), tissue.md, tissue.d₀, tissue.t_shrink)
        
        # cell diameter
        tissue.dc = cell_diameter(t_min + dt * (t - 1), tissue.mdc, tissue.dc₀, tissue.t_shrink, tissue.t_stop)

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

        # terminate when tissue short enough
        if tissue.xa + tissue.xd >= tissue.Xc

            return tracks[1:k-1]

        end

    end

    tracks# return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Tissue evolves dynamically according to the functions in dynamic_tissue_geometry.jl. Tracks the addition of cells in time.
Uses FCell structs to track the frequency in real time.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_dynamic_tissue_simulation_tracking_addition(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    cell_additions = Vector{Int64}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    # track the number of cells added in a 1 minute window
    old_total = tissue.total_cells

    for t = 1:length(t_min:dt:t_max)

        euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        add_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # update geometry
        # length
        tissue.L = Length(t_min + dt * (t - 1), tissue.uₐ, tissue.L₀, tissue.t_shrink)
        
        tissue.xa = xa(t_min + dt * (t - 1), tissue)

        # advection speed
        tissue.γₐ = 1.67 - wavefront_velocity(t_min + dt * (t - 1), tissue.uₐ, tissue.t_shrink)

        # radius
        tissue.r = radius(t_min + dt * (t - 1), tissue.mᵣ, tissue.r₀, tissue.t_shrink)

        # density
        tissue.ρ₀ = density(t_min + dt * (t - 1), tissue.md, tissue.d₀, tissue.t_shrink)
 
        # cell diameter
        tissue.dc = cell_diameter(t_min + dt * (t - 1), tissue.mdc, tissue.dc₀, tissue.t_shrink, tissue.t_stop)

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            cell_additions[k] = tissue.total_cells - old_total

            k += 1

            old_total = tissue.total_cells

        end

        # terminate when tissue short enough
        if tissue.xa + tissue.xd >= tissue.Xc

            return tracks[1:k-1], cell_additions[1:k-1]

        end

    end

    tracks, cell_additions# return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Tissue evolves dynamically according to the functions in dynamic_tissue_geometry.jl, except advection velocity is constant,
and does not change with the speed of the wavefront.
Uses FCell structs to track the frequency in real time.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_dynamic_tissue_simulation_constantadvection_tracking_addition(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

   # Get the initial condition
   cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

   # make a lattice object

   nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

   ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

   nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

   lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

   update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

   # Now 'relax' the tissue for 10 mins

   for t = 0:dt:9

       initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
       update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

   end

   # Create a waste array for cells anterior to the PSM
   old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

   # Initialise output array

   tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

   cell_additions = Vector{Int64}(undef, length(t_min:t_max))

   # Simulate

   k = 1

   # track the number of cells added in a 1 minute window
   old_total = tissue.total_cells

   for t = 1:length(t_min:dt:t_max)

       euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
       cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
       anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
       remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
       add_cells!(cells, phase, tissue) # add new cells to replenish psm
       update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

       # update geometry
       # length
       tissue.L = Length(t_min + dt * (t - 1), tissue.uₐ, tissue.L₀, tissue.t_shrink)
       
       tissue.xa = xa(t_min + dt * (t - 1), tissue)

       # radius
       tissue.r = radius(t_min + dt * (t - 1), tissue.mᵣ, tissue.r₀, tissue.t_shrink)

       # density
       tissue.ρ₀ = density(t_min + dt * (t - 1), tissue.md, tissue.d₀, tissue.t_shrink)

       # cell diameter
       tissue.dc = cell_diameter(t_min + dt * (t - 1), tissue.mdc, tissue.dc₀, tissue.t_shrink, tissue.t_stop)

       # save after every minute, and at the initial condition (t = 0)

       if (t - 1) % (1 / dt) == 0

           M = Matrix{Float64}(undef, tissue.total_cells, 8)

           for i = 1:tissue.total_cells

               if i <= tissue.N_cells # this value updates throughout the simulation

                   M[i, 1:3] .= cells[i].x
                   M[i, 4] = t_min + dt * (t - 1)
                   M[i, 5] = cells[i].trackid 
                   M[i, 6] = cells[i].θ
                   M[i, 7] = cells[i].dθdt
                   M[i, 8] = cells[i].τ

               elseif tissue.N_cells < i <= tissue.total_cells

                   M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                   M[i, 4] = t_min + dt * (t - 1)
                   M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                   M[i, 6] = old_cells[i - tissue.N_cells].θ
                   M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                   M[i, 8] = old_cells[i - tissue.N_cells].τ

               end

           end

           tracks[k] = M

           cell_additions[k] = tissue.total_cells - old_total

           k += 1

           old_total = tissue.total_cells

       end

       # terminate when tissue short enough, or too narrow to contain cells
       if tissue.xa + tissue.xd >= tissue.Xc + tissue.R + tissue.r || 2 * tissue.r < tissue.dc

           return tracks[1:k-1], cell_additions[1:k-1]

       end

   end

   tracks, cell_additions# return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Uses FCell structs to track the frequency in real time. Phase dynamics incorporate a delay in signalling between cells.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function tissue_simulation_with_delay(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    end

    # now create initial condition for the phase DDE

    history = Vector{Vector{FCell}}(undef, phase.nτ)

    size_of_tissue = Vector{Int64}(undef, phase.nτ)

    for i = 1:phase.nτ

        history[i] = cells

        size_of_tissue[i] = tissue.N_cells

    end

    # Create a waste array for cells anterior to the PSM

    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)

        # call garbage collection at random to reduce memory burden
        thr=0.005
        rand(Uniform(0,1)) < thr && GC.gc()

        euler_step_with_delay!(cells, history, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        new_cell = add_cells!(cells, phase, tissue) # add new cells to replenish psm

        # add history for the new cell, if it was added at all
        if new_cell !== nothing

            for i = 1:phase.nτ

                size_of_tissue[i] += 1

                history[i][size_of_tissue[i]] = new_cell

            end

        end

        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # update history, both of cell number and states
        update_state!(tissue.N_cells, size_of_tissue)
        update_state!(cells, history)

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Uses FCell structs to track the frequency in real time. Phase dynamics incorporate a delay in signalling between cells.
Cells added with initial phase = 0.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function adding_zero_cells_simulation_with_delay(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    end

    # now create initial condition for the phase DDE

    history = Vector{Vector{FCell}}(undef, phase.nτ)

    size_of_tissue = Vector{Int64}(undef, phase.nτ)

    for i = 1:phase.nτ

        history[i] = cells

        size_of_tissue[i] = tissue.N_cells

    end

    # Create a waste array for cells anterior to the PSM

    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        
        # call garbage collection at random to reduce memory burden
        thr=0.005
        rand(Uniform(0,1)) < thr && GC.gc()

        euler_step_with_delay!(cells, history, dt, phase, tissue, lattice)
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        new_cell = add_zero_cells!(cells, phase, tissue) # add new cells to replenish psm

        # add history for the new cell, if it was added at all
        if new_cell !== nothing

            for i = 1:phase.nτ

                size_of_tissue[i] += 1

                history[i][size_of_tissue[i]] = new_cell

            end

        end

        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # update history, both of cell number and states
        update_state!(tissue.N_cells, size_of_tissue)
        update_state!(cells, history)

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Uses FCell structs to track the frequency in real time. Phase dynamics incorporate a delay in signalling between cells.
Cells added in dorsal-posterior region with initial phase = 0 and random elsewhere in the tissue.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function adding_posterior_zero_cells_simulation_with_delay(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    end

    # now create initial condition for the phase DDE

    history = Vector{Vector{FCell}}(undef, phase.nτ)

    size_of_tissue = Vector{Int64}(undef, phase.nτ)

    for i = 1:phase.nτ

        history[i] = cells

        size_of_tissue[i] = tissue.N_cells

    end

    # Create a waste array for cells anterior to the PSM

    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        
        # call garbage collection at random to reduce memory burden
        thr=0.005
        rand(Uniform(0,1)) < thr && GC.gc()

        euler_step_with_delay!(cells, history, dt, phase, tissue, lattice)
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        new_cell = add_posterior_zero_cells!(cells, phase, tissue) # add new cells to replenish psm

        # add history for the new cell, if it was added at all
        if new_cell !== nothing

            for i = 1:phase.nτ

                size_of_tissue[i] += 1

                history[i][size_of_tissue[i]] = new_cell

            end

        end

        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # update history, both of cell number and states
        update_state!(tissue.N_cells, size_of_tissue)
        update_state!(cells, history)

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Uses FCell structs to track the frequency in real time. Phase dynamics incorporate a delay in signalling between cells.
Cells added on dorsal surface of posterior and ventral of two lateral regions with initial phase = 0.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function adding_posterior_lateral_zero_cells_simulation_with_delay(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    end

    # now create initial condition for the phase DDE

    history = Vector{Vector{FCell}}(undef, phase.nτ)

    size_of_tissue = Vector{Int64}(undef, phase.nτ)

    for i = 1:phase.nτ

        history[i] = cells

        size_of_tissue[i] = tissue.N_cells

    end

    # Create a waste array for cells anterior to the PSM

    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)

        # call garbage collection at random to reduce memory burden
        thr=0.005
        rand(Uniform(0,1)) < thr && GC.gc()

        euler_step_with_delay!(cells, history, dt, phase, tissue, lattice)
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        new_cell = add_posterior_lateral_zero_cells!(cells, phase, tissue) # add new cells to replenish psm

        # add history for the new cell, if it was added at all
        if new_cell !== nothing

            for i = 1:phase.nτ

                size_of_tissue[i] += 1

                history[i][size_of_tissue[i]] = new_cell

            end

        end

        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # update history, both of cell number and states
        update_state!(tissue.N_cells, size_of_tissue)
        update_state!(cells, history)

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end


"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Uses FCell structs to track the frequency in real time. Phase dynamics incorporate a delay in signalling between cells,
and arrest of the clock during M phase of the cell cycle.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function mitosing_simulation_with_delay(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    end
    
    # now create initial condition for the phase DDE

    history = Vector{Vector{FCell}}(undef, phase.nτ)

    size_of_tissue = Vector{Int64}(undef, phase.nτ)

    for i = 1:phase.nτ

        history[i] = cells

        size_of_tissue[i] = tissue.N_cells

    end

    # Create a waste array for cells anterior to the PSM

    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)

        euler_step_with_mitosis_and_delay!(cells, history, dt, phase, tissue, lattice)
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells 
        
        # Add cells via cell division and if they're added, create a history for them.
        mitosis_with_history!(cells, history, dt, phase, tissue, size_of_tissue)

        # add other cells if density too low
        new_cell = add_cells!(cells, phase, tissue) # add new cells to replenish psm

        # add history for the new cell, if it was added at all
        if new_cell !== nothing

            for i = 1:phase.nτ

                size_of_tissue[i] += 1

                history[i][size_of_tissue[i]] = new_cell

            end

        end

        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # update history, both of cell number and states
        update_state!(tissue.N_cells, size_of_tissue)
        update_state!(cells, history)

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Cells increase cell cycle phase τ linearly and divide when τ ≥ TG + TM. If τ ∈ [TG, TG+TM) then dθᵢ/dt = 0. Tg and TM increase towards
the anterior of the PSM. New cells are added at random positions within the PSM if the tissue density after division falls below requisite threshold,
and phase θ of these new cells is a random value within [0, 2π).\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_mitosing_simulation_varying_spatially(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)
    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue; length=tissue.N_cells + 20000) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
        
    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        # # call garbage collection at random to reduce memory burden
        # thr=0.01
        # rand(Uniform(0,1)) < thr && GC.gc()

        euler_step_with_mitosis_varying_spatially!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        cleanup_escaped_cells!(cells, tissue) # ugly fix to prevent the rapid proliferation of cells creating errors
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells 
        mitosis_varying_spatially!(cells, dt, tissue)
        add_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice
        
        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end


# DEFUNCT functions, those I will not use but keeping in case they are useful in the future

# """DEFUNCT"""
# function frequency_tracking_mitosing_accelerated_phase_simulation(seed::Real, t_min, t_max, dt, n::Real, tissue::TissueParams, phase::PhaseParams)
#     # Get the initial condition
#     cells = initial_FCells_horseshoe_PSM(seed, phase, tissue)

#     # Now 'relax' the tissue for 10 mins, so that 
#     for t = 0:dt:10
#         initial_relaxation!(cells, phase, tissue, dt, boundary_horseshoe_psm)
#     end

#     # Create a waste array for cells anterior to the PSM
#     old_cells = Array{FCell}(undef, (2 * 883 + 555) * 150);

#     # Initialise output array
#     tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))
#     # Simulate
#     k = 1
#     for t = 1:length(t_min:dt:t_max)
#         euler_step_with_accelerated_phase_mitosis!(cells, n, dt, phase, tissue) # simulate gene expression on cell positions at previous time point
#         cell_movements!(cells, phase, tissue, dt, boundary_horseshoe_psm) # move cells in the psm
#         anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
#         remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells
#         mitosis!(cells, dt, tissue) # remaining cells undergo cell cycling and divide when ending M phase
#         add_cells!(cells, phase, tissue) # add new cells to replenish psm
#         # save after every minute, and at the initial condition (t = 0)
#         if (t - 1) % (1 / dt) == 0
#             M = Matrix{Float64}(undef, tissue.total_cells, 8) 
#             for i = 1:tissue.total_cells
#                 if i <= tissue.N_cells # this value updates throughout the simulation
#                     M[i, 1:3] .= cells[i].x
#                     M[i, 4] = t_min + dt * (t - 1)
#                     M[i, 5] = cells[i].trackid 
#                     M[i, 6] = cells[i].θ
#                     M[i, 7] = cells[i].dθdt
#                     M[i, 8] = cells[i].τ
#                 elseif tissue.N_cells < i <= tissue.total_cells
#                     M[i, 1:3] .= old_cells[i - tissue.N_cells].x
#                     M[i, 4] = t_min + dt * (t - 1)
#                     M[i, 5] = old_cells[i - tissue.N_cells].trackid 
#                     M[i, 6] = old_cells[i - tissue.N_cells].θ
#                     M[i, 7] = old_cells[i - tissue.N_cells].dθdt
#                     M[i, 8] = old_cells[i - tissue.N_cells].τ
#                 end
#             end
#             tracks[k] = M
#             k += 1
#         end
#     end
#     tracks
# end

# """
# DEFUNCT

# Simulate gene expression dynamics and cell movements along a time course in a flattened and curved PSM.

# cells - initial condition of simulation.
# old_cells - 'waste' array in which to store cells that advect to somites.
# t_min - start time (mins)
# t_max - end time (mins)
# dt - time step (mins)
# boundary_force - force that defines the tissue geometry
# tissue - tissue parameters
# phase - phase parameters

# Returns a n_cells x 6 x n_times array.
# """
# function curved_tissue_simulation(
#     cells::Union{Vector{Cell}, Vector{ICell}}, 
#     old_cells::Union{Vector{Cell}, Vector{ICell}}, 
#     t_min, 
#     t_max, 
#     dt, 
#     tissue::TissueParams, 
#     phase::PhaseParams,
#     d::Function,
#     φ::Function,
#     invφ::Function)
#     # Initialise output
#     tracks = Array{Float64, 3}(undef, length(cells) + length(old_cells), 6, length(t_min:t_max))
#     # Simulate
#     k = 1
#     for t = 1:length(t_min:dt:t_max)
#         try
#         xcells = copy(cells)
#         xtissue = copy(tissue)
#         # print("\ncount cells before movement is $(zcount(cells, tissue))")
#         euler_step!(cells, dt, phase, tissue) # simulate gene expression on cell positions at previous time point
#         cell_movements!(cells, phase, tissue, dt, boundary_force_spatula, d, φ, invφ) # move cells in the psm
#         anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
#         # print("\ncount cells before removal is $(zcount(cells, tissue))")
#         remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells
#         # print("\ncount cells after removal is $(zcount(cells, tissue))")
#         add_cells!(cells, d, φ, tissue) # add new cells to replenish psm
#         if isnan(cells[1:tissue.N_cells])
#             print("\nNaN found at t = $t")
#             return xcells, xtissue
#             break
#         end
#         if t % (1 / dt) == 1 # save every minute, and the initial condition
#             for i = 1:(length(cells) + length(old_cells))
#                 if i <= tissue.N_cells # this value updates throughout the simulation
#                     tracks[i, 1:3, k] .= cells[i].x
#                     tracks[i, 4, k] = t_min + dt * (t - 1)
#                     tracks[i, 5, k] = cells[i].trackid 
#                     tracks[i, 6, k] = cells[i].θ
#                 elseif tissue.N_cells < i <= tissue.total_cells
#                     tracks[i, 1:3, k] .= old_cells[i - tissue.N_cells].x
#                     tracks[i, 4, k] = t_min + dt * (t - 1)
#                     tracks[i, 5, k] = old_cells[i - tissue.N_cells].trackid 
#                     tracks[i, 6, k] = old_cells[i - tissue.N_cells].θ
#                 end
#             end
#             k += 1
#         end
#     catch err
#         return tracks[1:tissue.total_cells, :, :]
#         rethrow(err)
#     end
#     end
#     tracks[1:tissue.total_cells, :, :] # return subset
# end

# """
# DEFUNCT

# Simulate gene expression dynamics and cell movements within growing horseshoe PSM.

# cells - initial condition of simulation.
# old_cells - 'waste' array in which to store cells that advect to somites.
# t_min - start time (mins)
# t_max - end time (mins)
# dt - time step (mins)
# tissue - tissue parameters
# phase - phase parameters
# t_gstart - time to start growth
# t_gend - time to end growth
# t_dstart - time to start decrease
# t_dend - time to end decrease
# g - rate of growth of length (μm/min)
# d - rate of shrink of length (μm/min)

# Returns a n_cells x 6 x n_times array.
# """
# function growing_psm_simulation(
#     cells::Union{Vector{Cell}, Vector{ICell}}, 
#     old_cells::Union{Vector{Cell}, Vector{ICell}}, 
#     t_min, t_max, dt, 
#     tissue::TissueParams, 
#     phase::PhaseParams,
#     t_gstart, 
#     t_gend, 
#     t_dstart, 
#     t_dend, 
#     g, 
#     d
#     )
#     # Initialise output
#     tracks = Array{Float64, 3}(undef, length(cells) + length(old_cells), 6, length(t_min:t_max))
#     # Simulate
#     k = 1
#     for t = 1:length(t_min:dt:t_max)
#         euler_step!(cells, dt, phase, tissue) # simulate gene expression on cell positions at previous time point
#         cell_movements!(cells, phase, tissue, dt, boundary_force) # move cells in the psm
#         anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
#         remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells
#         add_cells!(cells, tissue) # add new cells to replenish psm
#         psm_growth!(
#             t, 
#             x -> linear_growth_and_decrease(x, t_gstart, t_gend, t_dstart, t_dend, g, d, tissue.L), 
#             tissue)
#         if t % (1 / dt) == 1 # save every minute, and the initial condition
#             for i = 1:(length(cells) + length(old_cells))
#                 if i <= tissue.N_cells # this value updates throughout the simulation
#                     tracks[i, 1:3, k] .= cells[i].x
#                     tracks[i, 4, k] = t_min + dt * (t - 1)
#                     tracks[i, 5, k] = cells[i].trackid 
#                     tracks[i, 6, k] = cells[i].θ
#                 elseif length(cells) < i <= tissue.total_cells
#                     tracks[i, 1:3, k] .= old_cells[i - length(cells)].x
#                     tracks[i, 4, k] = t_min + dt * (t - 1)
#                     tracks[i, 5, k] = old_cells[i - length(cells)].trackid 
#                     tracks[i, 6, k] = old_cells[i - length(cells)].θ
#                 end
#             end
#             k += 1
#         end
#     end
#     tracks[1:tissue.total_cells, :, :] # return subset
# end

# """
# DEFUNCT

# Simulate gene expression dynamics and cell movements within shrinking horseshoe PSM.

# cells - initial condition of simulation.
# old_cells - 'waste' array in which to store cells that advect to somites.
# t_min - start time (mins)
# t_max - end time (mins)
# dt - time step (mins)
# tissue - tissue parameters
# phase - phase parameters
# t_start - time to start shrinking the psm
# t_end - time to end shrinking the psm
# shrink_rate - rate of psm radius shrinking (μm/min)
# start_r - initial value of radius

# Returns a n_cells x 6 x n_times array.
# """
# function shrinking_psm_simulation(
#     cells::Union{Vector{Cell}, Vector{ICell}}, 
#     old_cells::Union{Vector{Cell}, Vector{ICell}}, 
#     t_min, t_max, dt, 
#     tissue::TissueParams, 
#     phase::PhaseParams,
#     t_start,
#     t_end,
#     shrink_rate, 
#     start_r
#     )
#     # Initialise output
#     tracks = Array{Float64, 3}(undef, length(cells) + length(old_cells), 6, length(t_min:t_max))
#     # Simulate
#     k = 1
#     for t = 1:length(t_min:dt:t_max)
#         euler_step!(cells, dt, phase, tissue) # simulate gene expression on cell positions at previous time point
#         cell_movements!(cells, phase, tissue, dt, boundary_force) # move cells in the psm
#         anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
#         remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells
#         add_cells!(cells, tissue) # add new cells to replenish psm
#         shrink_psm!(t, t_start, t_end, shrink_rate, start_r, tissue)
#         if t % (1 / dt) == 1 # save every minute, and the initial condition
#             for i = 1:(length(cells) + length(old_cells))
#                 if i <= tissue.N_cells # this value updates throughout the simulation
#                     tracks[i, 1:3, k] .= cells[i].x
#                     tracks[i, 4, k] = t_min + dt * (t - 1)
#                     tracks[i, 5, k] = cells[i].trackid 
#                     tracks[i, 6, k] = cells[i].θ
#                 elseif length(cells) < i <= tissue.total_cells
#                     tracks[i, 1:3, k] .= old_cells[i - length(cells)].x
#                     tracks[i, 4, k] = t_min + dt * (t - 1)
#                     tracks[i, 5, k] = old_cells[i - length(cells)].trackid 
#                     tracks[i, 6, k] = old_cells[i - length(cells)].θ
#                 end
#             end
#             k += 1
#         end
#     end
#     tracks[1:tissue.total_cells, :, :] # return subset
# end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Uses FCell structs to track the frequency in real time. Here the advection profile is constant and non-zero across the PSM.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_tissue_simulation_constant_advection(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)

        euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements_without_advection!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(cells, tissue, phase, dt) # move cells
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        add_cells!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Extracellular volume increases towards the posterior of the tissue. Cells added at random. Advection is constant and non-zero.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function frequency_tracking_increasing_extracellular_volume_simulation_with_constant_advection(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    cells = initial_FCells_horseshoe_PSM_increasing_extracellular_volume(seed, phase, tissue; length = tissue.N_cells + 3000)

    # define lattice

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    for t = 1:dt:10
        initial_relaxation_with_increasing_extracellular_volume!(cells, tissue, dt, boundary_horseshoe_psm, d, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
    end

    old_cells = Vector{FCell}(undef, tissue.total_cells * 10);

    # Initialise output array
    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))
    # Simulate
    k = 1
    for t = 1:length(t_min:dt:t_max)

        euler_step!(cells, dt, phase, tissue, d, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements_with_constant_advection_with_increasing_extracellular_volume!(cells, tissue, dt, boundary_horseshoe_psm, d, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells
        cleanup_escaped_cells!(cells, tissue) # remove escaped cells
        add_cells_changing_diameter!(cells, phase, tissue) # add new cells to replenish psm
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

        # save after every minute, and at the initial condition (t = 0)
        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8) 

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M
            k += 1

        end

    end

    tracks # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021.
Uses FCell structs to track the frequency in real time. Tracks the time and position at which new cells are added.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function addition_tracking_tissue_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)

    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # initialise array for tracking the number of additions

    added_cells = zeros(4, length(t_min:dt:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)

        euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        added_cell = add_cells!(cells, phase, tissue) # add new cells to replenish psm
        if added_cell != nothing
            added_cells[1:3, t] .= added_cell.x
            added_cells[4, t] = t_min + (t - 1) * dt
        end
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    # as added_cells contains many zero-elements trim it for performance
    # columns with all elements = 0 cannot exist unless xd <= 0
    if tissue.xd > 0

        added_cells = added_cells[:, vec(mapslices(col -> any(col .!= 0), added_cells, dims = 1))]

    end

    tracks, added_cells # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021. Here, 
new cells are added to the dorsum of the toroid tailbud if the tissue density falls below requisite threshold,
and phase θ of these new cells is set to θ = 0. To maintain density, also add cells with θ ∈ [0, 2π) in the LHS and RHS PSM.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function addition_tracking_adding_posterior_zero_cells_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)
    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
        
    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # intialise array for recording added cells
    added_cells = zeros(4, length(t_min:dt:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        
        # # call garbage collection at random to reduce memory burden
        # thr=0.01
        # rand(Uniform(0,1)) < thr && GC.gc()

        euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        added_cell = add_posterior_zero_cells!(cells, phase, tissue) # add new cells to replenish psm
        if added_cell != nothing
            added_cells[1:3, t] .= added_cell.x
            added_cells[4, t] = t_min + (t - 1) * dt
        end
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    # as added_cells contains many zero-elements trim it for performance
    # columns with all elements = 0 cannot exist unless xd <= 0
    if tissue.xd > 0

        added_cells = added_cells[:, vec(mapslices(col -> any(col .!= 0), added_cells, dims = 1))]

    end

    tracks, added_cells # return subset

end

"""
Simulate gene expression dynamics and cell movements along a time course, using the tissue geometry of Uriu et al. 2021. Here, 
new cells are added to the dorsum of the toroid tailbud or ventral-lateral sides of the PSM if the tissue density falls below requisite threshold,
and phase θ of these new cells is set to θ = 0. Tracks the addition of new cells into the tissue.\n\n

seed - seed for reproducing the initial condition, and simulation\n
t_min - start time (mins)\n
t_max - end time (mins)\n
dt - time step (mins)\n
tissue - parameters defining tissue geometry and cell movements\n
phase - parameters specifying the kuramoto model for phase oscillators\n\n

Returns a Vector{Matrix{Float64}}, where each element M of the vector V 
corresponds to a time point of the simulation, and the elements of n x 8 matrix 
M correspond to:\n\n

M[i, 1:3] => xᵢ, yᵢ, zᵢ (position of cell i in space)\n
M[i, 4] => t (timepoint of the simulation)\n
M[i, 5] => trackid (ID of the cell i)\n
M[i, 6] => θ (segmentation clock phase)\n
M[i, 7] => dθᵢ/dt (frequency of phase oscillations)\n
M[i, 8] => τ (cell cycle phase)\n
"""
function tracking_addition_adding_posterior_lateral_zero_cells_simulation(seed::Real, t_min, t_max, dt, tissue::TissueParams, phase::PhaseParams)
    # Get the initial condition
    cells = initial_FCells_horseshoe_PSM(seed, phase, tissue) # 6 alloc

    # make a lattice object

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    lattice = Dictionary{Int64, Vector{Int64}}(1:nx*ny*nz, [Int64[] for i in 1:nx*ny*nz])

    update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)

    # Now 'relax' the tissue for 10 mins

    for t = 0:dt:9

        initial_relaxation!(cells, tissue, dt, boundary_horseshoe_psm, lattice)
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny)
        
    end

    # Create a waste array for cells anterior to the PSM
    old_cells = Vector{FCell}(undef, tissue.N_cells * 150);

    # Initialise output array

    tracks = Vector{Matrix{Float64}}(undef, length(t_min:t_max))

    # intialise array for recording added cells
    added_cells = zeros(4, length(t_min:dt:t_max))

    # Simulate

    k = 1

    for t = 1:length(t_min:dt:t_max)
        
        # # call garbage collection at random to reduce memory burden
        # thr=0.01
        # rand(Uniform(0,1)) < thr && GC.gc()

        euler_step!(cells, dt, phase, tissue, lattice) # simulate gene expression on cell positions at previous time point
        cell_movements!(cells, tissue, dt, boundary_horseshoe_psm, lattice) # move cells in the psm
        anterior_cell_movement!(old_cells, tissue, phase, dt) # move cells in the anterior
        remove_cells!(cells, old_cells, tissue) # remove psm cells that have gone into anterior and place in old_cells       
        added_cell = add_posterior_lateral_zero_cells!(cells, phase, tissue) # add new cells to replenish psm
        if added_cell != nothing
            added_cells[1:3, t] .= added_cell.x
            added_cells[4, t] = t_min + (t - 1) * dt
        end
        update_lattice!(cells[1:tissue.N_cells], lattice, tissue.lattice_side, nx, ny) # update the lattice

        # save after every minute, and at the initial condition (t = 0)

        if (t - 1) % (1 / dt) == 0

            M = Matrix{Float64}(undef, tissue.total_cells, 8)

            for i = 1:tissue.total_cells

                if i <= tissue.N_cells # this value updates throughout the simulation

                    M[i, 1:3] .= cells[i].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = cells[i].trackid 
                    M[i, 6] = cells[i].θ
                    M[i, 7] = cells[i].dθdt
                    M[i, 8] = cells[i].τ

                elseif tissue.N_cells < i <= tissue.total_cells

                    M[i, 1:3] .= old_cells[i - tissue.N_cells].x
                    M[i, 4] = t_min + dt * (t - 1)
                    M[i, 5] = old_cells[i - tissue.N_cells].trackid 
                    M[i, 6] = old_cells[i - tissue.N_cells].θ
                    M[i, 7] = old_cells[i - tissue.N_cells].dθdt
                    M[i, 8] = old_cells[i - tissue.N_cells].τ

                end

            end

            tracks[k] = M

            k += 1

        end

    end

    # as added_cells contains many zero-elements trim it for performance
    # columns with all elements = 0 cannot exist unless xd <= 0
    if tissue.xd > 0

        added_cells = added_cells[:, vec(mapslices(col -> any(col .!= 0), added_cells, dims = 1))]

    end

    tracks, added_cells # return subset

end
