# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

using Distributions

# Code used in my model:
# Below is code actually called in my simulations. At the end of the file
# there is a body of functions which I do not use, but have kept in case
# they are useful in the future.

# _______________________________________________________________________
# Used Code
# _______________________________________________________________________

# Code for calculating cell density

"""
Get the density of cells in the left and right columns of the PSM, and the tailbud.

Implementation from Uriu et al. 2021.
"""
function cell_density(X::Vector{FCell}, tissue::TissueParams)
    N_L = 0
    N_R = 0
    N_TB = 0
    for i = 1:tissue.N_cells
        if X[i].x[1] <= tissue.Xc && X[i].x[2] <= tissue.Yc # left hand
            N_L += 1
        elseif X[i].x[1] <= tissue.Xc && X[i].x[2] > tissue.Yc # right hand
            N_R += 1
        elseif X[i].x[1] > tissue.Xc # tailbud
            N_TB += 1
        end
    end
    density_L = N_L / (π * tissue.r^2 * (tissue.Xc - tissue.xa)) # volume of cylinder
    density_R = N_R / (π * tissue.r^2 * (tissue.Xc - tissue.xa)) # volume of cylinder
    density_TB = N_TB / (π^2 * tissue.r^2 * tissue.R) # half the volume of a torus
    return SVector{3, Float64}(density_L, density_R, density_TB)
end

"""Get the difference between actual density in each tissue region and desired value."""
function density_deltas(X::Vector{FCell}, tissue::TissueParams)
    cell_density(X, tissue) .- tissue.ρ₀
end

"""
Get the index corresponding to the tissue region with lowest density (1 - left PSM cylinder,
2 - right PSM cylinder, 3 - half-torus).

If there is a tie between these regions, it is broken by choosing one of the tied values at random.
"""
function get_lowest_density_region(X::Vector{FCell}, tissue::TissueParams)

    dd = density_deltas(X, tissue)
    min_d = minimum(dd)

    n_equal = 0

    for d in dd

        if d == min_d

            n_equal += 1

        end

    end

    q = rand()

    for (j, d) in enumerate(dd)

        if d == min_d && q <= j / n_equal

            return j

        end

    end


end

"""
Get the index corresponding to the tissue region with lowest density (1 - left PSM cylinder,
2 - right PSM cylinder, 3 - half-torus).

If there is a tie between these regions, it is broken by choosing one of the tied values at random.
"""
function get_lowest_density_region(X::Vector{FCell}, TB_density, col_density, tissue::TissueParams)

    dd = density_deltas(X, tissue, TB_density, col_density)
    min_d = minimum(dd)

    n_equal = 0

    for d in dd

        if d == min_d

            n_equal += 1

        end

    end

    for (j, d) in enumerate(dd)

        if d == min_d && rand() <= j / n_equal # check if less than or equal to 1 / n, i.e. equal probability of choosing any j

            return j

        end

    end

end


# Code for adding cells to tissue
#_______________________________________________________________________

"""
Add cells to X if density below desired threshold, with θ ∈ [0, 2π].
"""
function add_cells!(X::Vector{FCell}, phase::PhaseParams, tissue::TissueParams)
    i = get_lowest_density_region(X, tissue)
    if density_deltas(X, tissue)[i] < 0 #&& tissue.N_cells < length(X) # latter to avoid indexing errors
        if i == 1 && tissue.xa + tissue.xd < tissue.Xc # left hand domain
            # generate random y, z in polar coords
            rᵢ = rand(Uniform(0, tissue.r))
            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc - tissue.R + rᵢ * cos(pᵢ), # left hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )


            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            # θ = rand() * 2 * π # ∈ [0, 2π)
            θ = rand() * 2π
            τ = rand() * (tissue.TG + tissue.TM)
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)


        elseif i == 2 && tissue.xa + tissue.xd < tissue.Xc # right hand domain
            # generate random y, z in polar coords
            rᵢ = rand(Uniform(0, tissue.r))
            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc + tissue.R + rᵢ * cos(pᵢ), # right hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )


            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            # θ = rand() * 2 * π # ∈ [0, 2π)
            θ = rand() * 2π
            τ = rand() * (tissue.TG + tissue.TM)
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        elseif i == 3 || (tissue.xa + tissue.xd >= tissue.Xc && density_deltas(X, tissue)[3] < 0)
            # generate random x, y, z in spherical coords
            pᵢ = rand(Uniform(-π/2, π/2))
            # qᵢ = rand() * 2 * π # ∈ [0, 2π)
            qᵢ = rand() * 2π
            rᵢ = rand(Uniform(0, tissue.r))
            position = SVector{3, Float64}(
                tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
                tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
                tissue.Zc + rᵢ * sin(qᵢ)
            )

            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            # θ = rand() * 2 * π # ∈ [0, 2π)
            θ = rand() * 2π
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        end

    end
end

"""
Add cells to X if density below desired threshold, with θ = 0.
"""
function add_zero_cells!(X::Vector{FCell}, phase::PhaseParams, tissue::TissueParams)

    i = get_lowest_density_region(X, tissue)

    if density_deltas(X, tissue)[i] < 0 #&& tissue.N_cells < length(X) # latter to avoid indexing errors
        if i == 1 && tissue.xa + tissue.xd < tissue.Xc # left hand domain
            # generate random y, z in polar coords
            rᵢ = rand(Uniform(0, tissue.r))
            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc - tissue.R + rᵢ * cos(pᵢ), # left hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )
            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase, and initial phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            θ = 0.
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        elseif i == 2 && tissue.xa + tissue.xd < tissue.Xc # right hand domain
            # generate random y, z in polar coords
            rᵢ = rand(Uniform(0, tissue.r))
            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc + tissue.R + rᵢ * cos(pᵢ), # right hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )
            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase, and initial phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            θ = 0.
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        elseif i == 3 || (tissue.xa + tissue.xd >= tissue.Xc && density_deltas(X, tissue)[3] < 0) # tailbud
            # generate random x, y, z in spherical coords
            pᵢ = rand(Uniform(-π/2, π/2))
            # qᵢ = rand() * 2 * π # ∈ [0, 2π)
            qᵢ = rand() * 2π
            rᵢ = rand(Uniform(0, tissue.r))
            position = SVector{3, Float64}(
                tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
                tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
                tissue.Zc + rᵢ * sin(qᵢ)
            )
            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase, and initial phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            θ = 0.
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        end
    end
end

"""
Add cells to X if density below desired threshold, with θ = 0.
Cells added ventrally, i.e. pᵢ ∈ [5π/4, 7π/4] in the PSM, and dorso-posteriorly in the tailbud.
"""
function add_posterior_lateral_zero_cells!(X::Vector{FCell}, phase::PhaseParams, tissue::TissueParams)
    i = get_lowest_density_region(X, tissue)
    if density_deltas(X, tissue)[i] < 0 #&& tissue.N_cells < length(X) # latter to avoid indexing errors
        if i == 1 && tissue.xa + tissue.xd < tissue.Xc # left hand domain
            # generate random y, z in polar coords
            rᵢ = tissue.r
            pᵢ = rand(Uniform(5π/4, 7π/4))
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc - tissue.R + rᵢ * cos(pᵢ), # left hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )
            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase, and initial phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            θ = 0.
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)
        elseif i == 2 && tissue.xa + tissue.xd < tissue.Xc # right hand domain
            # generate random y, z in polar coords
            rᵢ = tissue.r
            pᵢ = rand(Uniform(5π/4, 7π/4))
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc + tissue.R + rᵢ * cos(pᵢ), # right hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )
            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase, and initial phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            θ = 0.
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)
        elseif i == 3 || (tissue.xa + tissue.xd >= tissue.Xc && density_deltas(X, tissue)[3] < 0) # tailbud
            # generate random x, y, z in spherical coords
            pᵢ = rand(Uniform(-π/2, π/2))
            # qᵢ = rand() * 2 * π # ∈ [0, 2π)
            qᵢ = rand(Uniform(π/4, 3π/4))
            rᵢ = tissue.r
            position = SVector{3, Float64}(
                tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
                tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
                tissue.Zc + rᵢ * sin(qᵢ)
            )
            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase, and initial phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            θ = 0.
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)
        end
    end
end

"""
Add cells to X if density below desired threshold, with θ = 0, at dorsal-posterior point, 
and at random in the PSM tubes.
"""
function add_posterior_zero_cells!(X::Vector{FCell}, phase::PhaseParams, tissue::TissueParams)
    i = get_lowest_density_region(X, tissue)
    if density_deltas(X, tissue)[i] < 0 #&& tissue.N_cells < length(X) # latter to avoid indexing errors
        if i == 1 && tissue.xa + tissue.xd < tissue.Xc # left hand domain
            # generate random y, z in polar coords
            rᵢ = rand(Uniform(0, tissue.r))
            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc - tissue.R + rᵢ * cos(pᵢ), # left hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )            
            
            # assign initial phase

            θ = rand() * 2π

            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        elseif i == 2 && tissue.xa + tissue.xd < tissue.Xc # right hand domain
            # generate random y, z in polar coords
            rᵢ = rand(Uniform(0, tissue.r))
            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc + tissue.R + rᵢ * cos(pᵢ), # right hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )

            # assign initial phase

            θ = rand() * 2π

            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        elseif i == 3 || (tissue.xa + tissue.xd >= tissue.Xc && density_deltas(X, tissue)[3] < 0)
            # generate random x, y, z in spherical coords
            pᵢ = rand(Uniform(-π/2, π/2))
            # qᵢ = rand() * 2 * π # ∈ [0, 2π)
            qᵢ = rand(Uniform(π/4, 3π/4))
            rᵢ = tissue.r
            position = SVector{3, Float64}(
                tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
                tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
                tissue.Zc + rᵢ * sin(qᵢ)
            )

            # assign initial phase
            θ = 0.

            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        end
        
    end

end

"""Get the difference between density in each tissue region and desired value."""
@inline function density_deltas(X, tissue::TissueParams, TB_density, col_density)
    d = cell_density(X, tissue)
    SVector{3, Float64}(d[1] - col_density, d[2] - col_density, d[3] - TB_density)
end

# Code for removing cells at the anterior
# _______________________________________________________________________

"""
Remove cells if cells exceed anterior boundary.
"""
function remove_cells!(X::Vector{FCell}, oldX::Vector{FCell}, tissue::TissueParams)
    new_X = Vector{FCell}(undef, length(X))
    j, k = 0, 0
    for i = 1:tissue.N_cells # as not every element in X is a real cell!
        if X[i].x[1] >= tissue.xa
            k += 1
            new_X[k] = X[i]
        else
            j += 1
            n_removed = tissue.total_cells - tissue.N_cells + j
            oldX[n_removed] = X[i] # add to the 'waste' array
        end
    end
    X .= new_X
    tissue.N_cells -= j # update the number of cells in the tissue
    k, j
end

# Initial conditions
#_______________________________________________________________________

"""Get the initial condition for a straight horseshoe PSM."""
function initial_FCells_horseshoe_PSM(seed, phase::PhaseParams, tissue::TissueParams; length=tissue.N_cells + 1000)
    # specify the initial condition
    N_torus = floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R) # cells in torus
    N_column = floor(tissue.ρ₀ * π * tissue.r^2 * (tissue.Xc - tissue.xa)) # cells in each column
    cells = Array{FCell}(undef, length) # adding 1000 because 500 fails
    Random.seed!(seed)
    for i = 1:tissue.N_cells
        if i <= N_torus
            # generate random x, y, z in spherical coords
            pᵢ = rand(Uniform(-π/2, π/2))
            # qᵢ = rand() * 2 * π # ∈ [0, 2π)
            qᵢ = rand() * 2π
            rᵢ = rand(Uniform(0, tissue.r))
            position = SVector{3, Float64}(
                tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
                tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
                tissue.Zc + rᵢ * sin(qᵢ)
            )
        elseif N_torus < i <= N_torus + N_column && tissue.Xc > tissue.xa
            # generate random y, z in polar coords
            rᵢ = rand(Uniform(0, tissue.r))
            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa, tissue.Xc)),
                tissue.Yc + tissue.R + rᵢ * cos(pᵢ), # right hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )
        elseif N_torus + N_column < i <= N_torus + 2 * N_column && tissue.Xc > tissue.xa
            # generate random y, z in polar coords
            rᵢ = rand(Uniform(0, tissue.r))
            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa, tissue.Xc)),
                tissue.Yc - tissue.R + rᵢ * cos(pᵢ), # left hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )
        else
            break
        end
        # ϕ = rand() * 2 * π # ∈ [0, 2π)
        ϕᵢ = rand()* 2π
        φᵢ = rand(Uniform(0, π))
        polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
            sin(φᵢ) * cos(ϕᵢ),
            sin(φᵢ) * sin(ϕᵢ),
            cos(φᵢ),)
        # θ = rand() * 2π
        # θ = 0
        # θ = π * ω / phase.ωₒ
        θ = 3. * π / 2.
        # Specify unformly distributed cell cycle phase
        if position[1] < tissue.x_div + tissue.xa
            τ = 0.0
        else
            τ = rand() * (tissue.TG + tissue.TM)
        end
        cells[i] = FCell(position, polarity, 0, θ, τ, i)
    end
    cells
end

"""Get the initial condition for a straight horseshoe PSM."""
function random_initial_FCells_horseshoe_PSM(seed, phase::PhaseParams, tissue::TissueParams; length=tissue.N_cells + 1000)
    # specify the initial condition
    N_torus = floor(tissue.ρ₀ * π^2 * tissue.r^2 * tissue.R) # cells in torus
    N_column = floor(tissue.ρ₀ * π * tissue.r^2 * tissue.Xc) # cells in each column
    cells = Array{FCell}(undef, length) # adding 1000 because 500 fails
    Random.seed!(seed)
    for i = 1:tissue.N_cells
        if i <= N_torus
            # generate random x, y, z in spherical coords
            pᵢ = rand(Uniform(-π/2, π/2))
            # qᵢ = rand() * 2 * π # ∈ [0, 2π)
            qᵢ = rand() * 2π
            rᵢ = rand(Uniform(0, tissue.r))
            position = SVector{3, Float64}(
                tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
                tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
                tissue.Zc + rᵢ * sin(qᵢ)
            )
        elseif N_torus < i <= N_torus + N_column && tissue.Xc > tissue.xa
            # generate random y, z in polar coords
            rᵢ = rand(Uniform(0, tissue.r))
            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa, tissue.Xc)),
                tissue.Yc + tissue.R + rᵢ * cos(pᵢ), # right hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )
        elseif N_torus + N_column < i <= N_torus + 2 * N_column && tissue.Xc > tissue.xa
            # generate random y, z in polar coords
            rᵢ = rand(Uniform(0, tissue.r))
            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa, tissue.Xc)),
                tissue.Yc - tissue.R + rᵢ * cos(pᵢ), # left hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )
        else
            break
        end
        # ϕ = rand() * 2 * π # ∈ [0, 2π)
        ϕᵢ = rand()* 2π
        φᵢ = rand(Uniform(0, π))
        polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
            sin(φᵢ) * cos(ϕᵢ),
            sin(φᵢ) * sin(ϕᵢ),
            cos(φᵢ),)
        θ = rand() * 2π
        # θ = 0
        # θ = π * ω / phase.ωₒ
        # θ = 3. * π / 2.
        # Specify unformly distributed cell cycle phase
        if position[1] < tissue.x_div + tissue.xa
            τ = 0.0
        else
            τ = rand() * (tissue.TG + tissue.TM)
        end
        cells[i] = FCell(position, polarity, 0, θ, τ, i)
    end
    cells
end


"""Get the initial condition for a straight horseshoe PSM, with increasing extracellular volume."""
function initial_FCells_horseshoe_PSM_increasing_extracellular_volume(seed, phase::PhaseParams, tissue::TissueParams; length=tissue.N_cells + 1000)
    
    # specify the initial condition

    cells = Array{FCell}(undef, length) # adding 1000 because 500 fails

    Random.seed!(seed)

    for i = 1:tissue.N_cells

        if i <= N_torus(tissue)

            # generate random x, y, z in spherical coords

            pᵢ = rand(Uniform(-π/2, π/2))

            qᵢ = rand() * 2π

            rᵢ = rand(Uniform(0, tissue.r))

            position = SVector{3, Float64}(
                tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
                tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
                tissue.Zc + rᵢ * sin(qᵢ)
            )

        elseif N_torus(tissue) < i <= N_torus(tissue) + N_column(tissue) && tissue.Xc > tissue.xa

            # generate random y, z in polar coords

            rᵢ = rand(Uniform(0, tissue.r))

            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)

            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa, tissue.Xc)),
                tissue.R + tissue.Yc + rᵢ * cos(pᵢ), # right hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )

        elseif N_torus(tissue) + N_column(tissue) < i <= N_torus(tissue) + 2 * N_column(tissue) && tissue.Xc > tissue.xa

            # generate random y, z in polar coords

            rᵢ = rand(Uniform(0, tissue.r))

            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)

            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa, tissue.Xc)),
                tissue.Yc - tissue.R + rᵢ * cos(pᵢ), # left hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )

        else
            break
        end
        # ϕ = rand() * 2 * π # ∈ [0, 2π)
        ϕᵢ = rand() * 2π

        φᵢ = rand(Uniform(0, π))

        polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
            sin(φᵢ) * cos(ϕᵢ),
            sin(φᵢ) * sin(ϕᵢ),
            cos(φᵢ),)

        # θ = rand() * 2π
        # θ = 0
        # θ = π * ω / phase.ωₒ
        θ = 3. * π / 2.

        # Specify unformly distributed cell cycle phase
        if position[1] < tissue.x_div + tissue.xa
            τ = 0.0
        else
            τ = rand() * (tissue.TG + tissue.TM)
        end

        cells[i] = FCell(position, polarity, 0, θ, τ, i)

    end

    cells

end


# DEFUNCT CODE
#_______________________________________________________________________

function TB_density(X, ϕ, tissue::TissueParams)
    N_TB = 0
    for i = 1:tissue.N_cells
        if X[i].x[1] > tissue.Xc # && ϕ >= atan(X[i].x[2] / X[i].x[2]) >= -ϕ # tailbud
            N_TB += 1
        end
    end
    N_TB / (2π * ϕ * tissue.r^2 * tissue.R) # half the volume of a torus
end

function density_deltas(X, ϕ, tissue::TissueParams)
    TB_density(X, ϕ, tissue) .- tissue.ρ₀
end

"""Get a random SVector for a cell's position within the flattened toroid."""
function cell_position_toroid(d::Function, tissue::TissueParams)
    Rᵢ = rand(Uniform(tissue.R - tissue.r, tissue.R + tissue.r + d(tissue.Xc)))
    if Rᵢ <= tissue.R
        pᵢ = rand(Uniform(-π/2, π/2))
        qᵢ = rand(Uniform(π/2, 3π/2))
        rᵢ = rand(Uniform(0, tissue.r))
        return SVector{3, Float64}(
            tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
            tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
            tissue.Zc + rᵢ * sin(qᵢ)
        )
    elseif tissue.R + d(tissue.Xc) >= Rᵢ > tissue.R
        pᵢ = rand(Uniform(-π/2, π/2))
        return SVector{3, Float64}(
            tissue.Xc + Rᵢ * cos(pᵢ),
            tissue.Yc + Rᵢ * sin(pᵢ),
            rand(Uniform(tissue.Zc - tissue.r, tissue.Zc + tissue.r))
        )
    elseif Rᵢ > tissue.R + d(tissue.Xc)
        pᵢ = rand(Uniform(-π/2, π/2))
        qᵢ = rand(Uniform(-π/2, π/2))
        rᵢ = rand(Uniform(0, tissue.r))
        return SVector{3, Float64}(
            tissue.Xc + (tissue.R + d(tissue.Xc)) * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
            tissue.Yc + (tissue.R + d(tissue.Xc)) * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
            tissue.Zc + rᵢ * sin(qᵢ)
        )
    else
        throw(DomainError(Rᵢ))
    end
end

"""Get a random SVector for a cell's position within the right hand flattened psm tube."""
function cell_position_right_psm(d::Function, φ::Function, tissue::TissueParams)
    xᵢ = rand(Uniform(tissue.xa, tissue.Xc))
    yᵢ = rand(Uniform(2 * tissue.R + d(tissue.Xc), 2 * tissue.R + d(tissue.Xc) + d(xᵢ) + 2 * tissue.r))
    if yᵢ <= 2 * tissue.R + d(tissue.Xc) + tissue.r
        qᵢ = rand(Uniform(π/2, 3π/2))
        rᵢ = rand(Uniform(0, tissue.r))
        return SVector{3, Float64}(
            tissue.Xc - sin(φ(xᵢ)) * (tissue.Zc + rᵢ * sin(qᵢ)),
            rᵢ * cos(qᵢ) + 2 * tissue.R + d(tissue.Xc) + tissue.r,
            -cos(φ(xᵢ)) * (tissue.Zc + rᵢ * sin(qᵢ))
        )
    elseif 2 * tissue.R + d(tissue.Xc) + d(xᵢ) + tissue.r >= yᵢ > 2 * tissue.R + d(tissue.Xc) + tissue.r
        return SVector{3, Float64}(
            tissue.Xc - sin(φ(xᵢ)) * (tissue.Zc + rand(Uniform(-tissue.r, tissue.r))),
            yᵢ,
            -cos(φ(xᵢ)) * (tissue.Zc + rand(Uniform(-tissue.r, tissue.r)))
        )
    elseif yᵢ > 2 * tissue.R + d(tissue.Xc) + d(xᵢ) + tissue.r
        qᵢ = rand(Uniform(-π/2, π/2))
        rᵢ = rand(Uniform(0, tissue.r))
        return SVector{3, Float64}(
            tissue.Xc - sin(φ(xᵢ)) * (tissue.Zc + rᵢ * sin(qᵢ)),
            2 * tissue.R + d(tissue.Xc) + d(xᵢ) + tissue.r + rᵢ * cos(qᵢ),
            -cos(φ(xᵢ)) * (tissue.Zc + rᵢ * sin(qᵢ))
        )
    else
        throw(DomainError(yᵢ))
    end
end

"""Get a random SVector for a cell's position within the left hand flattened psm tube."""
function cell_position_left_psm(d::Function, φ::Function, tissue::TissueParams)
    xᵢ = rand(Uniform(tissue.xa, tissue.Xc))
    yᵢ = rand(Uniform(d(tissue.Xc) - d(xᵢ), 2 * tissue.r + d(tissue.Xc)))
    if yᵢ <= tissue.r + d(tissue.Xc) - d(xᵢ)
        qᵢ = rand(Uniform(π/2, 3π/2))
        rᵢ = rand(Uniform(0, tissue.r))
        return SVector{3, Float64}(
            tissue.Xc - sin(φ(xᵢ)) * (tissue.Zc + rᵢ * sin(qᵢ)),
            rᵢ * cos(qᵢ) + tissue.r + d(tissue.Xc) - d(xᵢ),
            -cos(φ(xᵢ)) * (tissue.Zc + rᵢ * sin(qᵢ))
        )
    elseif tissue.r + d(tissue.Xc) >= yᵢ > tissue.r + d(tissue.Xc) - d(xᵢ)
        return SVector{3, Float64}(
            tissue.Xc - sin(φ(xᵢ)) * (tissue.Zc + rand(Uniform(-tissue.r, tissue.r))),
            yᵢ,
            -cos(φ(xᵢ)) * (tissue.Zc + rand(Uniform(-tissue.r, tissue.r)))
        )
    elseif yᵢ > tissue.r + d(tissue.Xc)
        qᵢ = rand(Uniform(-π/2, π/2))
        rᵢ = rand(Uniform(0, tissue.r))
        return SVector{3, Float64}(
            tissue.Xc - sin(φ(xᵢ)) * (tissue.Zc + rᵢ * sin(qᵢ)),
            tissue.r + d(tissue.Xc) + rᵢ * cos(qᵢ),
            -cos(φ(xᵢ)) * (tissue.Zc + rᵢ * sin(qᵢ))
        )
    else
        throw(DomainError(yᵢ))
    end
end

"""
DEFUNCT

Add cells to X if density below desired threshold, with θ = 0, at dorsal-posterior point.
"""
function add_posterior_zero_cells!(X::Vector{FCell}, ϕ::Float64, phase::PhaseParams, tissue::TissueParams)
    if density_deltas(X, tissue)[3] < 0 #&& tissue.N_cells < length(X) # latter to avoid indexing errors
        # generate random x, y, z in spherical coords
        pᵢ = rand(Uniform(-ϕ, ϕ))
        # qᵢ = rand() * 2 * π # ∈ [0, 2π)
        qᵢ = rand(Uniform(π/4, 3π/4))
        rᵢ = tissue.r
        position = SVector{3, Float64}(
            tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
            tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
            tissue.Zc + rᵢ * sin(qᵢ)
        )
        # Now assign polarity
        # ϕ = rand() * 2 * π # ∈ [0, 2π)
        ϕ = rand()* 2π
        φ = rand(Uniform(0, π))
        polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
            sin(φ) * cos(ϕ),
            sin(φ) * sin(ϕ),
            cos(φ),
        )
        # assign intrinsic frequency, cell cycle phase, and initial phase
        ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
        θ = 0.
        τ = rand() * (tissue.TG + tissue.TM)
        # Get the trackid by updating the number of cells
        tissue.total_cells += 1
        # Now update the number of cells in the tissue
        tissue.N_cells += 1
        # Add a new cell
        X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)
    end
end

"""
DEFUNCT

Add cells to X if density below desired threshold, with θ = [0, 2π), at dorsal-posterior point.
"""
function add_posterior_random_cells!(X::Vector{FCell}, phase::PhaseParams, tissue::TissueParams)
    if density_deltas(X, tissue)[3] < 0 #&& tissue.N_cells < length(X) # latter to avoid indexing errors
        # generate random x, y, z in spherical coords
        pᵢ = rand(Uniform(-π/2, π/2))
        # qᵢ = rand() * 2 * π # ∈ [0, 2π)
        qᵢ = rand(Uniform(π/4, 3π/4))
        rᵢ = tissue.r
        position = SVector{3, Float64}(
            tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
            tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
            tissue.Zc + rᵢ * sin(qᵢ)
        )
        # Now assign polarity
        # ϕ = rand() * 2 * π # ∈ [0, 2π)
        ϕ = rand() * 2π
        φ = rand(Uniform(0, π))
        polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
            sin(φ) * cos(ϕ),
            sin(φ) * sin(ϕ),
            cos(φ),
        )
        # assign intrinsic frequency, cell cycle phase, and initial phase
        ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
        θ = rand() * 2π
        τ = rand() * (tissue.TG + tissue.TM)
        # Get the trackid by updating the number of cells
        tissue.total_cells += 1
        # Now update the number of cells in the tissue
        tissue.N_cells += 1
        # Add a new cell
        X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)
    end
end

"""
DEFUNCT

Add cells to X if density below desired threshold, with θ = [0, 2π), at dorsal-posterior point.
"""
function add_posterior_random_cells!(X::Vector{FCell}, ϕ::Float64, phase::PhaseParams, tissue::TissueParams)
    if density_deltas(X, tissue)[3] < 0 #&& tissue.N_cells < length(X) # latter to avoid indexing errors
        # generate random x, y, z in spherical coords
        pᵢ = rand(Uniform(-ϕ, ϕ))
        # qᵢ = rand() * 2 * π # ∈ [0, 2π)
        qᵢ = rand(Uniform(π/4, 3π/4))
        rᵢ = tissue.r
        position = SVector{3, Float64}(
            tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
            tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
            tissue.Zc + rᵢ * sin(qᵢ)
        )
        # Now assign polarity
        # ϕ = rand() * 2 * π # ∈ [0, 2π)
        ϕ = rand() * 2π
        φ = rand(Uniform(0, π))
        polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
            sin(φ) * cos(ϕ),
            sin(φ) * sin(ϕ),
            cos(φ),
        )
        # assign intrinsic frequency, cell cycle phase, and initial phase
        ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
        θ = rand() * 2π
        τ = rand() * (tissue.TG + tissue.TM)
        # Get the trackid by updating the number of cells
        tissue.total_cells += 1
        # Now update the number of cells in the tissue
        tissue.N_cells += 1
        # Add a new cell
        X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)
    end
end

"""
DEFUNCT

Get the initial condition for a curved PSM.
"""
function initial_cells_curved_PSM(seed, d, φ, phase::PhaseParams, tissue::TissueParams)
    # specify the initial condition
    N_torus = 555 # cells in torus
    N_column = 883 # cells in each column
    N_cells = N_torus + 2 * N_column
    cells = Array{ICell}(undef, N_cells + 500)
    Random.seed!(seed)
    for i = 1:N_cells
        if i <= N_torus
            position = cell_position_toroid(d, tissue)
        elseif N_torus < i <= N_torus + N_column
            position = cell_position_right_psm(d, φ, tissue)
        elseif N_torus + N_column < i <= N_torus + 2 * N_column
            position = cell_position_left_psm(d, φ, tissue)
        else
            break
        end
        # ϕ = rand() * 2 * π # ∈ [0, 2π)
        ϕᵢ = rand()* 2π
        φᵢ = rand(Uniform(0, π))
        polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
            sin(φᵢ) * cos(ϕᵢ),
            sin(φᵢ) * sin(ϕᵢ),
            cos(φᵢ),
        )
        ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
        # θ = rand() * 2π
        # θ = 0
        # θ = π * ω / phase.ωₒ
        θ = 3 * π / 2
        # Specify unformly distributed cell cycle phase
        τ = rand() * (tissue.TG + tissue.TM)
        cells[i] = ICell(position, polarity, ω, θ, τ, i)
    end
    cells
end

"""
Add cells to X if density below desired threshold (according to changing diameter), with θ ∈ [0, 2π].
"""
function add_cells_changing_diameter!(X::Vector{FCell}, phase::PhaseParams, tissue::TissueParams)

    i = get_lowest_density_region(X, density_torus(tissue), density_column(tissue), tissue)

    if density_deltas(X, tissue, density_torus(tissue), density_column(tissue))[i] < 0 #&& tissue.N_cells < length(X) # latter to avoid indexing errors
        
        if i == 1 && tissue.xd + tissue.xa < tissue.Xc # left hand domain

            # generate random y, z in polar coords

            rᵢ = rand(Uniform(0, tissue.r))

            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)

            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc - tissue.R + rᵢ * cos(pᵢ), # left hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )

            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π

            φ = rand(Uniform(0, π))

            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )

            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)

            # θ = rand() * 2 * π # ∈ [0, 2π)
            θ = rand() * 2π

            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end

            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        elseif i == 2 && tissue.xd + tissue.xa < tissue.Xc # right hand domain

            # generate random y, z in polar coords

            rᵢ = rand(Uniform(0, tissue.r))

            pᵢ = rand() * 2π
            # pᵢ = rand() * 2 * π # ∈ [0, 2π)

            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc + tissue.R + rᵢ * cos(pᵢ), # right hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )

            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π

            φ = rand(Uniform(0, π))

            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )

            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)

            # θ = rand() * 2 * π # ∈ [0, 2π)
            θ = rand() * 2π

            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end

            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        elseif i == 3 || (tissue.xa + tissue.xd >= tissue.Xc && density_deltas(X, tissue)[3] < 0) # tailbud

            # generate random x, y, z in spherical coords

            pᵢ = rand(Uniform(-π/2, π/2))

            # qᵢ = rand() * 2 * π # ∈ [0, 2π)
            qᵢ = rand() * 2π

            rᵢ = rand(Uniform(0, tissue.r))

            position = SVector{3, Float64}(
                tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
                tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
                tissue.Zc + rᵢ * sin(qᵢ)
            )

            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π

            φ = rand(Uniform(0, π))

            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )

            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)

            # θ = rand() * 2 * π # ∈ [0, 2π)
            θ = rand() * 2π

            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end

            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        end

    end

end

"""
Add cells to X if density below desired threshold, with θ = 0, at dorsal-posterior point, 
and at random in the PSM tubes.
"""
function add_posterior_zero_cells_changing_diameter!(X::Vector{FCell}, phase::PhaseParams, tissue::TissueParams)
    
    i = get_lowest_density_region(X, density_torus(tissue), density_column(tissue), tissue)

    if density_deltas(X, tissue, density_torus(tissue), density_column(tissue))[i] < 0
        if i == 1 && tissue.xa + tissue.xd < tissue.Xc # left hand domain
            # generate random y, z in polar coords
            rᵢ = rand(Uniform(0, tissue.r))
            pᵢ = rand() * 2π
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc - tissue.R + rᵢ * cos(pᵢ), # left hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )            
            
            # assign initial phase

            θ = rand() * 2π

            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        elseif i == 2 && tissue.xa + tissue.xd < tissue.Xc # right hand domain
            # generate random y, z in polar coords
            rᵢ = rand(Uniform(0, tissue.r)) # ∈ [0, r]
            pᵢ = rand() * 2π # ∈ [0, 2π)
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc + tissue.R + rᵢ * cos(pᵢ), # right hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )

            # assign initial phase

            θ = rand() * 2π

            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        elseif i == 3 || (tissue.xa + tissue.xd >= tissue.Xc && density_deltas(X, tissue)[3] < 0) # tailbud
            # generate random x, y, z in spherical coords
            pᵢ = rand(Uniform(-π/2, π/2))
            # qᵢ = rand() * 2 * π # ∈ [0, 2π)
            qᵢ = rand(Uniform(π/4, 3π/4))
            rᵢ = tissue.r
            position = SVector{3, Float64}(
                tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
                tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
                tissue.Zc + rᵢ * sin(qᵢ)
            )

            # assign initial phase
            θ = 0.

            # Now assign polarity
            # ϕ = rand() * 2 * π # ∈ [0, 2π)
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        end
        
    end

end

"""
Add cells to X if density below desired threshold, with θ = 0.
Cells added ventrally, i.e. pᵢ ∈ [5π/4, 7π/4] in the PSM, and dorso-posteriorly in the tailbud.
"""
function add_posterior_lateral_zero_cells_changing_diameter!(X::Vector{FCell}, phase::PhaseParams, tissue::TissueParams)

    i = get_lowest_density_region(X, density_torus(tissue), density_column(tissue), tissue)

    if density_deltas(X, tissue, density_torus(tissue), density_column(tissue))[i] < 0

        if i == 1 && tissue.xa + tissue.xd < tissue.Xc # left hand domain
            # generate random y, z in polar coords
            rᵢ = tissue.r
            pᵢ = rand(Uniform(5π/4, 7π/4))
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc - tissue.R + rᵢ * cos(pᵢ), # left hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )

            # Now assign polarity
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase, and initial phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            θ = 0.
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        elseif i == 2 && tissue.xa + tissue.xd < tissue.Xc # right hand domain
            # generate random y, z in polar coords
            rᵢ = tissue.r
            pᵢ = rand(Uniform(5π/4, 7π/4))
            position = SVector{3, Float64}(
                rand(Uniform(tissue.xa + tissue.xd, tissue.Xc)),
                tissue.Yc + tissue.R + rᵢ * cos(pᵢ), # right hand column
                tissue.Zc + rᵢ * sin(pᵢ)
            )

            # Now assign polarity
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase, and initial phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            θ = 0.
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        elseif i == 3 || (tissue.xa + tissue.xd >= tissue.Xc && density_deltas(X, tissue)[3] < 0) # tailbud
            # generate random x, y, z in spherical coords
            pᵢ = rand(Uniform(-π/2, π/2))
            qᵢ = rand(Uniform(π/4, 3π/4))
            rᵢ = tissue.r
            position = SVector{3, Float64}(
                tissue.Xc + tissue.R * cos(pᵢ) + rᵢ * cos(pᵢ) * cos(qᵢ),
                tissue.Yc + tissue.R * sin(pᵢ) + rᵢ * sin(pᵢ) * cos(qᵢ),
                tissue.Zc + rᵢ * sin(qᵢ)
            )

            # Now assign polarity
            ϕ = rand() * 2π
            φ = rand(Uniform(0, π))
            polarity = SVector{3, Float64}( # random cell orientation - note that this must be on the unit sphere!
                sin(φ) * cos(ϕ),
                sin(φ) * sin(ϕ),
                cos(φ),
            )
            # assign intrinsic frequency, cell cycle phase, and initial phase
            ω = frequency(position[1], tissue.xa, tissue.L, phase.ωₒ, phase.σ, phase.k)
            θ = 0.
            if position[1] < tissue.x_div + tissue.xa
                τ = 0.0
            else
                τ = rand() * (tissue.TG + tissue.TM)
            end
            # Get the trackid by updating the number of cells
            tissue.total_cells += 1
            # Now update the number of cells in the tissue
            tissue.N_cells += 1
            # Add a new cell
            X[tissue.N_cells] = FCell(position, polarity, ω, θ, τ, tissue.total_cells)

        end
        
    end
end

