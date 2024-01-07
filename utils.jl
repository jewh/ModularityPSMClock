# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Misc utility functions

include("cell_structs.jl")
include("tissue_boundary.jl")

"""
Remove cells if cells escape the tissue.\n

This function should only be used in really extreme scenarios, e.g. where the number of cells escaping the tissue is small relative to the total number,
    and are causing unresolvable numerical errors. Increasing the boundary force tissue.μb is always preferable but not always suitable, as it changes
    the movement dynamics of cells within the tissue. Warns the user if >1 cell is lost from the tissue in a timestep, indicating broader issues to the
    simulation.
"""
function cleanup_escaped_cells!(X::Vector{FCell}, tissue::TissueParams)
    new_X = Vector{FCell}(undef, length(X))
    j, k = 0, 0
    for i = 1:tissue.N_cells # as not every element in X is a real cell!
        if in_tissue(X[i].x, tissue)
            k += 1
            new_X[k] = X[i]
        else
            j += 1
        end
    end
    X .= new_X
    tissue.N_cells -= j # update the number of cells in the tissue

    if j > 1
        @warn "$j cells escaped tissue in one timestep. Check tissue for loss of density if this happens repeatedly."
    end

    k, j
end
