# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023


# Functions to evolve the vector n for cell cell which determines its random directionality

include("cell_structs.jl")

using LinearAlgebra

"""
Update the cell polarity vector n based on the Langevin Equations in Uriu et al. 2021.
"""
function update_polarity(cell::FCell, dt, Dφ)
    # Don't test equality of vectors below, as that creates allocations. Instead do this ugly thing.
    if cell.n[1] == 0 && cell.n[2] == 0 && (cell.n[3] == 1 || cell.n[3] == -1)  # cell pointing directly north or south
        mx = SVector{3, Float64}(0, -1, 0)
        my = SVector{3, Float64}(-1, 0, 0) # clearest to see this when writing mx and my in spherical coordinates
    else
        ez = SVector{3, Float64}(0, 0, 1)
        mx = SVector{3, Float64}(cross(cell.n, ez)) / norm(cross(cell.n, ez))
        my = SVector{3, Float64}(cross(mx, cell.n)) / norm(cross(mx, cell.n))
    end
    return cell.n - dt * 2 * Dφ * cell.n + sqrt(2 * Dφ) * (mx * randn() * sqrt(dt) + my * randn() * sqrt(dt))
end

"""
Update the cell polarity vector n based on the Langevin Equations in Uriu et al. 2021.
"""
function update_polarity(n::SVector{3, Float64}, dt, Dφ)
    # Don't test equality of vectors below, as that creates allocations. Instead do this ugly thing.
    if n[1] == 0 && n[2] == 0 && (n[3] == 1 || n[3] == -1)  # cell pointing directly north or south
        mx = SVector{3, Float64}(0, -1, 0)
        my = SVector{3, Float64}(-1, 0, 0) # clearest to see this when writing mx and my in spherical coordinates
    else
        ez = SVector{3, Float64}(0, 0, 1)
        mx = SVector{3, Float64}(cross(n, ez)) / norm(cross(n, ez))
        my = SVector{3, Float64}(cross(mx, n)) / norm(cross(mx, n))
    end
    return n - dt * 2 * Dφ * n + sqrt(2 * Dφ) * (mx * randn() * sqrt(dt) + my * randn() * sqrt(dt))
end
