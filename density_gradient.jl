# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# DEFUNCT

# Functions which attempt to create a density gradient in the tissue. Unsuccessfully!

include("cell_structs.jl")

@inline function d(x, xa, Xc, dm, dc)

    if x < xa

        return 0

    elseif x > Xc

        return dm + dc

    else

        return dm * ( ( x - xa ) / (Xc - xa) ) + dc

    end

end

d(x::SVector{3, Float64}, dm, tissue::TissueParams) = d(x[1], tissue.xa, tissue.Xc, dm, tissue.dc)

d(x::SVector{3, Float64}, tissue::TissueParams) = d(x, tissue.dm, tissue)

d(x::Real, tissue::TissueParams) = d(x, tissue.xa, tissue.Xc, tissue.dm, tissue.dc)

sphere_packing(tissue::TissueParams) = tissue.ρ₀ * tissue.dc^3 / 6.

@inline function N_column(tissue::TissueParams)
    N = 3. * sphere_packing(tissue) * π * tissue.r^2 * (tissue.Xc - tissue.xa) / tissue.dm
    N *= (1 / d(tissue.xa, tissue)^2 - 1 / d(tissue.Xc, tissue)^2) * 1.5
    ceil(N)
end

# N_column(tissue::TissueParams) = 10

density_column(tissue::TissueParams) = N_column(tissue) / (π * tissue.r^2 * (tissue.Xc - tissue.xa))

N_torus(tissue::TissueParams) = ceil(6. * sphere_packing(tissue) * π^2 * tissue.r^2 * tissue.R / (tissue.dm + tissue.dc)^3 * 1.5)

# N_torus(tissue::TissueParams) = 10

density_torus(tissue::TissueParams) = 6. * sphere_packing(tissue) / (tissue.dm + tissue.dc)^3

function μ(x::Real, xa, Xc, μ)

    if x < xa

        return 0

    elseif x > Xc

        return 30 * μ

    else

        return 29 * μ * ( ( x - xa ) / (Xc - xa) ) + μ

    end

end

μ(x::SVector{3, Float64}, tissue::TissueParams) = μ(x[1], tissue.xa, tissue.Xc, tissue.μ)

"""
F(xᵢ, xⱼ, diam, params)

Scalar magnitude of intercellular force between cells i and j, with increasing radius of repulsion d/2.
"""
@inline function F_x(xᵢ::AbstractVector{Float64}, xⱼ::AbstractVector{Float64}, diam::Function, tissue::TissueParams)

    d = 2 * norm(xᵢ - xⱼ) / ( diam(xᵢ, tissue) + diam(xⱼ, tissue) )

    if d <= 1

        return μ(xᵢ, tissue) * (d - 1)

    else # not adjacent!

        return 0.0

    end

end

"""
intercellular_force(xᵢ, xⱼ, diam, params)

Returns the force between two cells, where cell diameter controlled by function diam, and force magnitude varies with space.
"""
intercellular_force(xᵢ::AbstractVector{Float64}, xⱼ::AbstractVector{Float64}, diam::Function, params::TissueParams) = (xⱼ - xᵢ) * (F_x(xᵢ, xⱼ, diam, params) / norm(xⱼ - xᵢ))

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
