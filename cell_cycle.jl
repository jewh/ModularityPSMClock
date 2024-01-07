# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# This file contains functions for implementing cell division within the tissue.

include("cell_structs.jl")
include("tissue_boundary.jl")

"""
p(τ, tissue)

Get the cell cycle phase for time τ.
"""
function p(τ, tissue::TissueParams)
    if 0 <= τ <= tissue.TG
        return "G"
    elseif tissue.TG <= τ <= tissue.TG + tissue.TM
        return "M"
    end
end

"""
Checks if the proposed position of a new cell that is rₘ away from a point x in the direction specified by the angles ϕ, φ is in the tissue specified by the parameters in tissue::TissueParams.

Use this to prevent cell division creating cells that lie outside the tissue.
"""
function is_in_tissue(rₘ::Float64, ϕ::Float64, φ::Float64, x::SVector{3, Float64}, tissue::TissueParams)

    if x[1] ≥ tissue.Xc

        rᵢ, pᵢ, qᵢ = spherical_coords_torus(x, tissue) # from tissue_boundary

        new_distance = (rᵢ * cos(pᵢ) * cos(qᵢ) + rₘ * cos(ϕ) * sin(φ))^2
        new_distance += (rᵢ * sin(pᵢ) * cos(qᵢ) + rₘ * sin(ϕ) * sin(φ))^2
        new_distance += (rᵢ * sin(qᵢ) + rₘ * cos(φ))^2

        return new_distance <= tissue.r^2

    elseif x[2] < tissue.Yc # lhs

        yᵢ = x[2]
        zᵢ = x[3]

        new_distance = (yᵢ - tissue.Yc + tissue.R + rₘ * sin(ϕ) * sin(φ)) ^ 2 
        new_distance += (zᵢ - tissue.Zc + rₘ * cos(φ)) ^ 2

        return new_distance <= tissue.r^2

    elseif x[2] > tissue.Yc # rhs

        yᵢ = x[2]
        zᵢ = x[3]

        new_distance = (yᵢ - tissue.Yc - tissue.R + rₘ * sin(ϕ) * sin(φ)) ^ 2 
        new_distance += (zᵢ - tissue.Zc + rₘ * cos(φ)) ^ 2

        return new_distance <= tissue.r^2

    end

end


"""

Goes through each cell in the vector cells::Vector{FCell} and checks if they've come to the end of their cell cycle, and adds spawns a new daughter cell adjacent to them if so. Also updates the cell cycle phase τ for each cell.

"""
function mitosis!(cells::Vector{FCell}, dt::Float64, tissue::TissueParams)

    new_cells = copy(cells)
    # define two variables here to avoid chance of infinite loop

    N_cells = tissue.N_cells

    total_cells = tissue.total_cells

    for i = 1:tissue.N_cells

        if cells[i].x[1] < tissue.x_div + tissue.xa

            new_cells[i] = FCell(
                cells[i].x,
                cells[i].n, 
                cells[i].dθdt, 
                cells[i].θ, 
                0.0, # set tau to 0
                cells[i].trackid
            )

        elseif cells[i].τ + dt > tissue.TG + tissue.TM

            # update tau
            τ = cells[i].τ + dt - tissue.TG - tissue.TM

            # cell undergoes mitosis, create a sister cell

            # define direction and distance from mother cell

            ϕ = 2 * π * rand() # ∈ [0, 2π)

            φ = rand(Uniform(0, π)) # ∈ [0, π]

            rₘ = 0.1 * tissue.dc

            # need to check that the new position of the sister cell is not outside the tissue

            if !is_in_tissue(rₘ, ϕ, φ, cells[i].x, tissue)

                if !is_in_tissue(rₘ, ϕ + π, π - φ, cells[i].x, tissue) && is_in_tissue(rₘ, ϕ + π/2, ϕ + π/2, cells[i].x, tissue)

                    ϕ = ϕ + π/2
                    φ = φ + π/2

                elseif !is_in_tissue(rₘ, ϕ + π, π - φ, cells[i].x, tissue) && is_in_tissue(rₘ, ϕ - π/2, φ - π/2, cells[i].x, tissue)

                    ϕ = ϕ - π/2
                    φ = φ - π/2

                else

                    ϕ = ϕ + π
                    φ = π - φ

                end

            end

            new_cells[N_cells + 1] = FCell( # going to try adding cells in the same place.
                cells[i].x + rₘ * SVector{3, Float64}(
                    cos(ϕ) * sin(φ),
                    sin(ϕ) * sin(φ),
                    cos(φ)
                )
                ,
                cells[i].n, 
                cells[i].dθdt, 
                cells[i].θ, 
                τ,
                total_cells + 1
            )

            # update the number of cells
            N_cells += 1
            total_cells += 1
            # re-assign 'parent' cell to update τ
            new_cells[i] = FCell(
                cells[i].x,
                cells[i].n, 
                cells[i].dθdt, 
                cells[i].θ, 
                τ,
                cells[i].trackid
            )

        else # just update cell cycle phase
            new_cells[i] = FCell(
                cells[i].x,
                cells[i].n, 
                cells[i].dθdt, 
                cells[i].θ, 
                cells[i].τ + dt,
                cells[i].trackid
            )
        end

    end

    # re-assign 
    cells .= new_cells
    tissue.total_cells = total_cells
    tissue.N_cells = N_cells

end

"""Method for delay simulations where one must create a history for the newly created cells."""
# needs to edit history
function mitosis_with_history!(cells::Vector{FCell}, history::Vector{Vector{FCell}}, dt::Float64, phase::PhaseParams, tissue::TissueParams, size_of_tissue::Vector{Int64})

    new_cells = copy(cells)

    N_cells = tissue.N_cells

    total_cells = tissue.total_cells

    for i = 1:tissue.N_cells

        if cells[i].x[1] < tissue.x_div + tissue.xa

            new_cells[i] = FCell(
                cells[i].x,
                cells[i].n, 
                cells[i].dθdt, 
                cells[i].θ, 
                0.0, # set tau to 0
                cells[i].trackid
            )

        elseif cells[i].τ + dt > tissue.TG + tissue.TM

            # update tau
            τ = cells[i].τ + dt - tissue.TG - tissue.TM

            # cell undergoes mitosis, create a sister cell

            # define direction and distance from mother cell

            ϕ = 2 * π * rand() # ∈ [0, 2π)

            φ = rand(Uniform(0, π)) # ∈ [0, π]

            rₘ = 0.1 * tissue.dc

            # need to check that the new position of the sister cell is not outside the tissue

            if !is_in_tissue(rₘ, ϕ, φ, cells[i].x, tissue)

                if !is_in_tissue(rₘ, ϕ + π, π - φ, cells[i].x, tissue) && is_in_tissue(rₘ, ϕ + π/2, ϕ + π/2, cells[i].x, tissue)

                    ϕ = ϕ + π/2
                    φ = φ + π/2

                elseif !is_in_tissue(rₘ, ϕ + π, π - φ, cells[i].x, tissue) && is_in_tissue(rₘ, ϕ - π/2, φ - π/2, cells[i].x, tissue)

                    ϕ = ϕ - π/2
                    φ = φ - π/2

                else

                    ϕ = ϕ + π
                    φ = π - φ

                end

            end

            new_cells[N_cells + 1] = FCell( # going to try adding cells in the same place.
                cells[i].x + rₘ * SVector{3, Float64}(
                    cos(ϕ) * sin(φ),
                    sin(ϕ) * sin(φ),
                    cos(φ)
                )
                ,
                cells[i].n, 
                cells[i].dθdt, 
                cells[i].θ, 
                τ,
                total_cells + 1
            )

            # re-assign 'parent' cell to update τ
            new_cells[i] = FCell(
                cells[i].x,
                cells[i].n, 
                cells[i].dθdt, 
                cells[i].θ, 
                τ,
                cells[i].trackid
            )

            # now need to create a history for this

            for j = 1:phase.nτ

                size_of_tissue[j] += 1

                history[j][size_of_tissue[j]] = new_cells[N_cells + 1]

            end

            # update the number of cells
            N_cells += 1
            total_cells += 1

        else # just update cell cycle phase
            new_cells[i] = FCell(
                cells[i].x,
                cells[i].n, 
                cells[i].dθdt, 
                cells[i].θ, 
                cells[i].τ + dt,
                cells[i].trackid
            )
        end

    end

    # re-assign 
    cells .= new_cells
    tissue.total_cells = total_cells
    tissue.N_cells = N_cells

end

dτdt(τ, θ, kc, TG, TM) = 1 + kc * sin((θ + π) / 6. - 2π * τ / (TG + TM))

function mitosis_coupled_to_phase!(cells::Vector{FCell}, dt::Float64, tissue::TissueParams)
    new_cells = copy(cells)
    # define two variables here to avoid chance of infinite loop
    N_cells = tissue.N_cells
    total_cells = tissue.total_cells
    for i = 1:tissue.N_cells
        if cells[i].τ + dt * dτdt(cells[i].τ, cells[i].θ, tissue.kc, tissue.TG, tissue.TM) > tissue.TG + tissue.TM
            # update tau
            τ = cells[i].τ + dt * dτdt(cells[i].τ, cells[i].θ, tissue.kc, tissue.TG, tissue.TM) - tissue.TG - tissue.TM
            # cell undergoes mitosis, create a sister cell
            ϕ = 2 * π * rand() # ∈ [0, 2π)
            φ = rand(Uniform(0, π)) # ∈ [0, π]
            # add sister cell to end of array
            new_cells[N_cells + 1] = FCell( # going to try adding cells in the same place.
                cells[i].x + 0.1 * tissue.dc * SVector{3, Float64}(
                    cos(ϕ) * sin(φ),
                    sin(ϕ) * sin(φ),
                    cos(φ)
                )
                ,
                cells[i].n, 
                cells[i].dθdt, 
                cells[i].θ, 
                τ,
                total_cells + 1
            )
            # update the number of cells
            N_cells += 1
            total_cells += 1
            # re-assign cell to update τ
            new_cells[i] = FCell(
                cells[i].x,
                cells[i].n, 
                cells[i].dθdt, 
                cells[i].θ, 
                τ,
                cells[i].trackid
            )
        else # just update cell cycle phase
            new_cells[i] = FCell(
                cells[i].x,
                cells[i].n, 
                cells[i].dθdt, 
                cells[i].θ, 
                cells[i].τ + dt * dτdt(cells[i].τ, cells[i].θ, tissue.kc, tissue.TG, tissue.TM),
                cells[i].trackid
            )
        end
    end
    # re-assign 
    cells .= new_cells
    tissue.total_cells = total_cells
    tissue.N_cells = N_cells
end

