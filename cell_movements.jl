# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# The file contains functions for movement of cells within the tissue in response to cell-cell repulsion, cell advection, random motility, and the tissue boundary.

include("cell_structs.jl")
include("cell_direction.jl")
include("tissue_boundary.jl")
include("box_splitting.jl")
include("density_gradient.jl")

"""
F(xᵢ, xⱼ, params)

Scalar magnitude of intercellular force between cells i and j.
"""
function F(xᵢ::AbstractVector{Float64}, xⱼ::AbstractVector{Float64}, params::TissueParams)
    d = norm(xᵢ - xⱼ) / params.dc
    if d <= 1
        return params.μ * (d - 1)
    else # not adjacent!
        return 0.0
    end
end

"""
F(xᵢ, xⱼ, diam, params)

Scalar magnitude of intercellular force between cells i and j, with increasing radius of repulsion d/2.
"""
@inline function F(xᵢ::AbstractVector{Float64}, xⱼ::AbstractVector{Float64}, diam::Function, tissue::TissueParams)

    d = 2 * norm(xᵢ - xⱼ) / ( diam(xᵢ, tissue) + diam(xⱼ, tissue) )

    if d <= 1

        return tissue.μ * (d - 1)

    else # not adjacent!

        return 0.0

    end

end

"""
intercellular_force(xᵢ, xⱼ, params)

Returns the force between two cells.
"""
intercellular_force(xᵢ::AbstractVector{Float64}, xⱼ::AbstractVector{Float64}, params::TissueParams) = (xⱼ - xᵢ) * (F(xᵢ, xⱼ, params) / norm(xⱼ - xᵢ))

"""
intercellular_force(xᵢ, xⱼ, diam, params)

Returns the force between two cells, where cell diameter controlled by function diam.
"""
intercellular_force(xᵢ::AbstractVector{Float64}, xⱼ::AbstractVector{Float64}, diam::Function, params::TissueParams) = (xⱼ - xᵢ) * (F(xᵢ, xⱼ, diam, params) / norm(xⱼ - xᵢ))

"""
vd(xᵢ, params)

Compute the scalar value of advection to the anterior.
"""
function vd(xᵢ::AbstractVector{Float64}, params::TissueParams)
    x = (xᵢ[1] - params.xa) / params.L
    if x <= params.xq
        return params.γₐ - (params.γₐ - (params.γₚ * (1 - params.xq))) * x / params.xq 
    elseif x > params.xq
        return - params.γₚ * x + params.γₚ
    end
end

"""
vd(xᵢ, invφ, params)

Compute the scalar value of advection to the anterior.
"""
function vd(xᵢ::AbstractVector{Float64}, invφ::Function, params::TissueParams)
    φ = atan((xᵢ[1] - params.Xc) / xᵢ[3])
    x = (invφ(φ) - params.xa) / params.L
    if x <= params.xq
        return params.γₐ - (params.γₐ - (params.γₚ * (1 - params.xq))) * x / params.xq 
    elseif x > params.xq
        return - params.γₚ * x + params.γₚ
    end
end

"""
advection_velocity(xᵢ, params)

Velocity vector for cell movements.
"""
advection_velocity(xᵢ::AbstractVector{Float64}, params::TissueParams) = SVector{3, Float64}(-vd(xᵢ, params), 0, 0)

"""
advection_velocity(xᵢ, invφ, params)

Velocity vector for cell movements with angle φ.
"""
function advection_velocity(xᵢ::AbstractVector{Float64}, invφ::Function, params::TissueParams)
    φ = atan((xᵢ[1] - params.Xc) / xᵢ[3])
    SVector{3, Float64}(-vd(xᵢ, invφ, params) * cos(φ), 0, vd(xᵢ, invφ, params) * sin(φ))
end

"""
intrinstic_motion(xᵢ, params)

Get the scalar intrinsic motion of each cell, modelled as a gradient in space.
"""
function intrinsic_motion(xᵢ, params::TissueParams)
    frac = (1 - (xᵢ[1] - params.xa) / params.L) / params.Xᵥ
    return params.vs / (1 + frac ^ params.h)
end

"""
Cell advection anterior to the PSM.
"""
function anterior_cell_movement!(oldX::Vector{FCell}, tissue_params::TissueParams, phase_params::PhaseParams, dt::Float64)
    new_oldX = Vector{FCell}(undef, length(oldX))
    for i = 1:tissue_params.total_cells - tissue_params.N_cells # total number of removed cells
        x = copy(oldX[i].x)
        x -= dt * tissue_params.γₐ * SVector{3, Float64}(1, 0, 0) # assuming no change in psm length and segment size 50 μm.
        new_oldX[i] = FCell(x, oldX[i].n, oldX[i].dθdt, oldX[i].θ, oldX[i].τ, oldX[i].trackid)
    end
    oldX .= new_oldX
end

"""
Slower (O(n²)) method

Move cells in the vector X according to one forward euler step of the numerical solution to the equation of motion used by Uriu et al. 2021.\n\n

Arguments:\n\n

X::Vector{FCell} => Cells in the tissue\n
tissue_params::TissueParams => TissueParams object containing tissue-wide parameters for the simulation\n
dt::Float64 => time-step for the euler scheme\n
boundary_force::Function => boundary force containing cells within the tissue\n

"""
function cell_movements!(X::Vector{FCell}, tissue_params::TissueParams, dt::Float64, boundary_force::Function)
    new_X = Vector{FCell}(undef, length(X))
    for i = 1:tissue_params.N_cells
        x = copy(X[i].x) # 0 allocations
        for j = 1:tissue_params.N_cells
            if i != j
                x += dt * intercellular_force(X[i].x, X[j].x, tissue_params) # 0 allocs
            end
        end
        x += dt * advection_velocity(X[i].x, tissue_params) # 0 allocs
        x += dt * intrinsic_motion(X[i].x, tissue_params) * X[i].n # 0 allocs
        x += dt * boundary_force(X[i].x, tissue_params) # 1 allocs
        # Now update the polarity vector
        n = update_polarity(X[i], dt, tissue_params.Dφ) # 0 allocs
        # Create new cell object
        new_X[i] = FCell(x, n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid)
    end
    # Now output
    X .= new_X
end

"""
Faster (O(nlogn)) method

Move cells in the vector X according to one forward euler step of the numerical solution to the equation of motion used by Uriu et al. 2021.\n\n

Arguments:\n\n

X::Vector{FCell} => Cells in the tissue\n
tissue_params::TissueParams => TissueParams object containing tissue-wide parameters for the simulation\n
dt::Float64 => time-step for the euler scheme\n
boundary_force::Function => boundary force containing cells within the tissue\n
lattice::Dictionary{Int64, Vector{Int64}} => object mapping regions of space onto the indices of the cells contained within them\n

"""
function cell_movements!(X::Vector{FCell}, tissue_params::TissueParams, dt::Float64, boundary_force::Function, lattice::Dictionary{Int64, Vector{Int64}})
    
    new_X = Vector{FCell}(undef, length(X))

    for i = 1:tissue_params.N_cells

        x = copy(X[i].x) # 0 allocations

        for j in get_candidate_neighbour_cells(X[i], lattice, tissue_params)

            if i != j

                x += dt * intercellular_force(X[i].x, X[j].x, tissue_params) # 0 allocs

            end

        end

        x += dt * advection_velocity(X[i].x, tissue_params) # 0 allocs

        x += dt * intrinsic_motion(X[i].x, tissue_params) * X[i].n # 0 allocs

        x += dt * boundary_force(X[i].x, tissue_params) # 1 allocs

        # Now update the polarity vector

        n = update_polarity(X[i], dt, tissue_params.Dφ) # 0 allocs

        # Create new cell object
        new_X[i] = FCell(x, n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid)

    end

    # Now output

    X .= new_X

end

"""
Faster (O(nlogn)) method\n\n

X::Vector{FCell} => cells within the tissue\n
tissue_params::TissueParams => parameters for tissue-wide properties\n
dt::Float64 => time-step for the euler scheme\n
boundary_force::Function => function defining the tissue boundaries\n
diam::Function => function describing the diameters of cells\n
lattice::Dictionary{Int64, Vector{Int64}} => object mapping regions of space to the cells within them\n

"""
function cell_movements!(X::Vector{FCell}, tissue_params::TissueParams, dt::Float64, boundary_force::Function, diam::Function, lattice::Dictionary{Int64, Vector{Int64}})
    
    new_X = Vector{FCell}(undef, length(X))

    for i = 1:tissue_params.N_cells

        x = copy(X[i].x) # 0 allocations

        for j in get_candidate_neighbour_cells(X[i], lattice, tissue_params)

            if i != j

                x += dt * intercellular_force(X[i].x, X[j].x, diam, tissue_params) # 0 allocs

            end

        end

        x += dt * advection_velocity(X[i].x, tissue_params) # 0 allocs

        x += dt * intrinsic_motion(X[i].x, tissue_params) * X[i].n # 0 allocs

        x += dt * boundary_force(X[i].x, tissue_params) # 1 allocs

        # Now update the polarity vector

        n = update_polarity(X[i], dt, tissue_params.Dφ) # 0 allocs

        # Create new cell object
        new_X[i] = FCell(x, n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid)

    end

    # Now output
    X .= new_X
end

"""
Slow (O(n²)) method

Move cells in the vector X according to one forward euler step of 
the numerical solution to the equation of motion used by Uriu et al. 2021, but without tissue advection, and confining cells
within the tissue at the anterior boundary.\n\n

Arguments:\n\n

X::Vector{FCell} => Cells in the tissue\n
tissue_params::TissueParams => TissueParams object containing tissue-wide parameters for the simulation\n
dt::Float64 => time-step for the euler scheme\n
boundary_force::Function => boundary force containing cells within the tissue\n

"""
function initial_relaxation!(X::Vector{FCell}, tissue_params::TissueParams, dt::Float64, boundary_force::Function)

    new_X = Vector{FCell}(undef, length(X))

    for i = 1:tissue_params.N_cells

        x = copy(X[i].x) # 0 allocations

        for j = 1:tissue_params.N_cells

            if i != j

                x += dt * intercellular_force(X[i].x, X[j].x, tissue_params) # 0 allocs

            end

        end

        x += dt * intrinsic_motion(X[i].x, tissue_params) * X[i].n # 0 allocs

        x += dt * boundary_force(X[i].x, tissue_params) # 0 allocs

        # Contain cells in the PSM tubes by a boundary force at the anterior.

        if X[i].x[1] < tissue_params.Xc

            x += dt * tissue_params.μb * exp(tissue_params.xa - X[i].x[1]) / tissue_params.rb * SVector{3, Float64}(1, 0, 0)

        end

        # Now update the polarity vector

        n = update_polarity(X[i], dt, tissue_params.Dφ) # 0 allocs

        # Create new cell object

        new_X[i] = FCell(x, n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid)

    end

    # Now output

    X .= new_X

end


"""
Faster (O(nlogn)) method

Move cells in the vector X according to one forward euler step of 
the numerical solution to the equation of motion used by Uriu et al. 2021, but without tissue advection, and confining cells
within the tissue at the anterior boundary.\n\n

Arguments:\n\n

X::Vector{FCell} => Cells in the tissue\n
tissue_params::TissueParams => TissueParams object containing tissue-wide parameters for the simulation\n
dt::Float64 => time-step for the euler scheme\n
boundary_force::Function => boundary force containing cells within the tissue\n
lattice::Dictionary{Int64, Vector{Int64}} => object mapping regions of space onto the indices of the cells contained within them\n

"""
function initial_relaxation!(X::Vector{FCell}, tissue_params::TissueParams, dt::Float64, boundary_force::Function, boxes::Dictionary{Int64, Vector{Int64}})

    new_X = Vector{FCell}(undef, length(X))

    for i = 1:tissue_params.N_cells

        x = copy(X[i].x) # 0 allocations

        for j in get_candidate_neighbour_cells(X[i], boxes, tissue_params) # 3 allocs

            if i != j # otherwise Infs

                x += dt * intercellular_force(X[i].x, X[j].x, tissue_params) # 0 allocs

            end

        end

        x += dt * intrinsic_motion(X[i].x, tissue_params) * X[i].n # 0 allocs

        x += dt * boundary_force(X[i].x, tissue_params) # 0 allocs

        # Contain cells in the PSM tubes by a boundary force at the anterior.
        if X[i].x[1] < tissue_params.Xc

            x += dt * tissue_params.μb * exp(tissue_params.xa - X[i].x[1]) / tissue_params.rb * SVector{3, Float64}(1, 0, 0)

        end

        # Now update the polarity vector

        n = update_polarity(X[i], dt, tissue_params.Dφ) # 0 allocs

        # Create new cell object
        new_X[i] = FCell(x, n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid)

    end

    # Now output
    X .= new_X

end

# DEFUNCT
function initial_relaxation_with_increasing_extracellular_volume!(X::Vector{FCell}, tissue_params::TissueParams, dt::Float64, boundary_force::Function, diam::Function, lattice::Dictionary{Int64, Vector{Int64}})
    
    new_X = Vector{FCell}(undef, length(X))

    for i = 1:tissue_params.N_cells

        x = copy(X[i].x) # 0 allocations

        for j in get_candidate_neighbour_cells(X[i], lattice, tissue_params)

            if i != j

                x += dt * intercellular_force(X[i].x, X[j].x, diam, tissue_params) # 0 allocs

            end

        end

        x += dt * intrinsic_motion(X[i].x, tissue_params) * X[i].n # 0 allocs

        x += dt * boundary_force(X[i].x, tissue_params) # 0 allocs

        # Contain cells in the PSM tubes by a boundary force at the anterior.
        if X[i].x[1] < tissue_params.Xc
            
            x += dt * tissue_params.μb * exp(tissue_params.xa - X[i].x[1]) / tissue_params.rb * SVector{3, Float64}(1, 0, 0)
        
        end

        # Now update the polarity vector
        n = update_polarity(X[i], dt, tissue_params.Dφ) # 0 allocs

        # Create new cell object
        new_X[i] = FCell(x, n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid)

    end

    # Now output
    X .= new_X

end

"""
Slow (O(n²)) method

Move cells in the vector X according to one forward euler step of 
the numerical solution to the equation of motion used by Uriu et al. 2021, but without tissue advection.\n\n

Arguments:\n\n

X::Vector{FCell} => Cells in the tissue\n
tissue_params::TissueParams => TissueParams object containing tissue-wide parameters for the simulation\n
dt::Float64 => time-step for the euler scheme\n
boundary_force::Function => boundary force containing cells within the tissue\n

"""
function cell_movements_without_advection!(X::Vector{FCell}, tissue_params::TissueParams, dt::Float64, boundary_force::Function)

    new_X = Vector{FCell}(undef, length(X))

    for i = 1:tissue_params.N_cells

        x = copy(X[i].x) # 0 allocations

        for j = 1:tissue_params.N_cells

            if i != j

                x += dt * intercellular_force(X[i].x, X[j].x, tissue_params) # 0 allocs

            end

        end

        x += dt * intrinsic_motion(X[i].x, tissue_params) * X[i].n # 0 allocs

        x += dt * boundary_force(X[i].x, tissue_params) # 0 allocs

        # Now update the polarity vector

        n = update_polarity(X[i], dt, tissue_params.Dφ) # 0 allocs

        # Create new cell object

        new_X[i] = FCell(x, n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid)

    end

    # Now output

    X .= new_X

end

# """DEFUNCT Second method, with the arguments for boundary_force, d and φ."""
# function cell_movements!(
#     X::Vector{ICell},
#     phase_params::PhaseParams,
#     tissue_params::TissueParams, 
#     dt::Float64, 
#     boundary_force::Function, 
#     d::Function, 
#     φ::Function,
#     invφ::Function)
#     new_X = Vector{ICell}(undef, length(X))
#     if isnan(X[1:tissue_params.N_cells])
#         throw(ArgumentError("Passed a vector containing NaNs."))
#     end
#     for i = 1:tissue_params.N_cells
#         if any(isnan, X[i].x)
#             throw(ArgumentError("Passed a NaN-containing position, xᵢ = $(X[i].x), for i = $i."))
#         end
#         x = copy(X[i].x)
#         for j = 1:tissue_params.N_cells
#             if i != j
#                 x += dt * intercellular_force(X[i].x, X[j].x, tissue_params)
#                 if any(isnan, x)
#                     throw(ArgumentError("intercellular force creates NaNs, for i = $i, xᵢ = $(X[i].x), j = $j, xⱼ = $(X[j].x)."))
#                 elseif any(isinf, x)
#                     throw(ArgumentError("intercellular force creates Infs, for i = $i, xᵢ = $(X[i].x), j = $j, xⱼ = $(X[j].x)."))
#                 end
#             end
#         end
#         x += dt * advection_velocity(X[i].x, invφ, tissue_params) # 0 allocs
#         if any(isnan, x)
#             throw(ArgumentError("advection velocity creates NaNs, for i = $i, xᵢ = $(X[i].x)."))
#         elseif any(isinf, x)
#             throw(ArgumentError("advection velocity creates Infs, for i = $i, xᵢ = $(X[i].x)."))
#         end
#         x += dt * intrinsic_motion(X[i].x, tissue_params) * X[i].n # 0 allocs
#         if any(isnan, x)
#             throw(ArgumentError("intrinsic motion creates NaNs, for i = $i, xᵢ = $(X[i].x)."))
#         elseif any(isinf, x)
#             throw(ArgumentError("intrinsic motion creates Infs, for i = $i, xᵢ = $(X[i].x)."))
#         end
#         try
#             x += dt * boundary_force(X[i].x, d, φ, invφ, tissue_params) # 2 allocs
#             if any(isnan, x)
#                 throw(ArgumentError("boundary force creates NaNs, for i = $i, xᵢ = $(X[i].x)."))
#             elseif any(isinf, x)
#                 throw(ArgumentError("boundary force creates Infs, for i = $i, xᵢ = $(X[i].x)."))
#             end
#         catch err
#             println("Error trying to multiply boundary_force(X[i].x, d, φ, tissue_params = $(boundary_force(X[i].x, d, φ, invφ, tissue_params))")
#             println("Error in cell i = $i")
#             rethrow(err)
#         end
#         # Now update the polarity vector
#         n = update_polarity(X[i], dt, tissue_params.Dφ) # 0 allocs
#         if any(isnan, n)
#             print("\nNaN at n = $(X[i].n), i = $i")
#         end
#         # And update the frequency ω
#         if new_X[i].x[1] >= tissue_params.xa # avoids NaNs for -ve (x - xa)
#             ω = frequency(new_X[i].x[1], tissue_params.xa, tissue_params.L, phase_params.ωₒ, phase_params.σ, phase_params.k) # 0 allocs
#         else
#             ω = X[i].ω
#         end
#         # Create new cell object
#         new_X[i] = ICell(x, n, ω, X[i].θ, X[i].τ, X[i].trackid)
#     end
#     # Now output
#     X .= new_X
# end

# """
# DEFUNCT

# Cell advection anterior to the PSM, at an angle φ.
# """
# function anterior_cell_movement!(oldX::Vector{ICell}, tissue_params::TissueParams, phase_params::PhaseParams, dt::Float64, φ::Float64)
#     new_oldX = Vector{ICell}(undef, length(oldX))
#     for i = 1:tissue_params.total_cells - tissue_params.N_cells # total number of removed cells
#         x = copy(oldX[i].x)
#         x -= dt * tissue_params.γₐ * SVector{3, Float64}(cos(φ), 0, -sin(φ)) # assuming no change in psm length and segment size 50 μm.
#         new_oldX[i] = ICell(x, oldX[i].n, oldX[i].ω, oldX[i].θ, oldX[i].τ, oldX[i].trackid)
#     end
#     oldX .= new_oldX
# end


"""
Faster (O(nlogn)) method

Move cells in the vector X according to one forward euler step of 
the numerical solution to the equation of motion used by Uriu et al. 2021, but without tissue advection.\n\n

Arguments:\n\n

X::Vector{FCell} => Cells in the tissue\n
tissue_params::TissueParams => TissueParams object containing tissue-wide parameters for the simulation\n
dt::Float64 => time-step for the euler scheme\n
boundary_force::Function => boundary force containing cells within the tissue\n
lattice::Dictionary{Int64, Vector{Int64}} => object mapping regions of space onto the indices of the cells contained within them\n

"""
function cell_movements_without_advection!(X::Vector{FCell}, tissue_params::TissueParams, dt::Float64, boundary_force::Function, boxes::Dictionary{Int64, Vector{Int64}})

    new_X = Vector{FCell}(undef, length(X))

    for i = 1:tissue_params.N_cells

        x = copy(X[i].x) # 0 allocations

        for j in get_candidate_neighbour_cells(X[i], boxes, tissue_params) # 3 allocs

            if i != j # otherwise Infs

                x += dt * intercellular_force(X[i].x, X[j].x, tissue_params) # 0 allocs

            end

        end

        x += dt * intrinsic_motion(X[i].x, tissue_params) * X[i].n # 0 allocs

        x += dt * boundary_force(X[i].x, tissue_params) # 0 allocs

        # Now update the polarity vector

        n = update_polarity(X[i], dt, tissue_params.Dφ) # 0 allocs

        # Create new cell object
        new_X[i] = FCell(x, n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid)

    end

    # Now output
    X .= new_X

end

function cell_movements_without_advection_with_increasing_extracellular_volume!(X::Vector{FCell}, tissue_params::TissueParams, dt::Float64, boundary_force::Function, diam::Function, lattice::Dictionary{Int64, Vector{Int64}})
    
    new_X = Vector{FCell}(undef, length(X))

    for i = 1:tissue_params.N_cells

        x = copy(X[i].x) # 0 allocations

        for j in get_candidate_neighbour_cells(X[i], lattice, tissue_params)

            if i != j

                x += dt * intercellular_force(X[i].x, X[j].x, diam, tissue_params) # 0 allocs

            end

        end

        x += dt * intrinsic_motion(X[i].x, tissue_params) * X[i].n # 0 allocs

        x += dt * boundary_force(X[i].x, tissue_params) # 0 allocs

        # Now update the polarity vector
        n = update_polarity(X[i], dt, tissue_params.Dφ) # 0 allocs

        # Create new cell object
        new_X[i] = FCell(x, n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid)

    end

    # Now output
    X .= new_X

end

function cell_movements_with_constant_advection_with_increasing_extracellular_volume!(X::Vector{FCell}, tissue_params::TissueParams, dt::Float64, boundary_force::Function, diam::Function, lattice::Dictionary{Int64, Vector{Int64}})
    
    new_X = Vector{FCell}(undef, length(X))

    for i = 1:tissue_params.N_cells

        x = copy(X[i].x) # 0 allocations

        for j in get_candidate_neighbour_cells(X[i], lattice, tissue_params)

            if i != j

                x += dt * intercellular_force(X[i].x, X[j].x, diam, tissue_params) # 0 allocs

            end

        end

        # advection
        x -= dt * tissue.γₐ * SVector{3, Float64}(1, 0, 0)

        x += dt * intrinsic_motion(X[i].x, tissue_params) * X[i].n # 0 allocs

        x += dt * boundary_force(X[i].x, tissue_params) # 0 allocs

        # Now update the polarity vector
        n = update_polarity(X[i], dt, tissue_params.Dφ) # 0 allocs

        # Create new cell object
        new_X[i] = FCell(x, n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid)

    end

    # Now output
    X .= new_X

end
