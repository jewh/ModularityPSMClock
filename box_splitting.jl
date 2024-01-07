# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# box_splitting.jl

# This splits the spatial domain into boxes for performance. We then use a dictionary object to efficiently assign each box to a list of cells which are in that box. This dictionary is read by functions such as euler_step! or cell_movements! where the function needs to know the cells which might be adjacent to a given cell. This is done by checking the distance between the given cell and those in the ≤27 boxes that include and border the box containing the given cell. 

using StaticArrays, LinearAlgebra, Random, Dictionaries

"""
Assigns each point in space a positive integer, corresponding the nth box dividing the spatial domain.

x::SVector{3, Float64} => point in space
d::Float64 => side of boxes
nx::Int64 => number of boxes in the x direction (floor(maximum(x) / dc) + 1)
ny::Int64 => number of boxes in the y direction (floor(maximum(y) / dc) + 1)

"""
@inline function box(x::SVector{3, Float64}, d::Float64, nx::Int64, ny::Int64)

    xs = floor(x[1] / d) + 1
    ys = floor(x[2] / d) + 1
    zs = floor(x[3] / d) + 1

    Int64((zs - 1) * nx * ny + (ys - 1) * nx + xs)

end

"""
Gets the list of integers corresponding to boxes which border a given box. Returned vector also includes the given box.

WARNING - if a box is near to the edge of the domain or in a corner, some of these values may be non-positive.

box::Int64 => given box
nx::Int64 => number of boxes in the x direction (floor(maximum(x) / dc) + 1)
ny::Int64 => number of boxes in the y direction (floor(maximum(y) / dc) + 1)

"""
@inline function neighbouring_boxes(box::Int64, nx::Int64, ny::Int64)

    neigbs = @SVector [
        box, 
        box + 1, 
        box - 1, 
        box + nx, 
        box - nx, 
        box + nx * ny, 
        box - nx * ny,
        box + 1 + nx * ny,
        box - 1 + nx * ny,
        box + 1 - nx * ny,
        box - 1 - nx * ny,
        box + nx + nx * ny,
        box - nx + nx * ny,
        box + nx - nx * ny,
        box - nx - nx * ny,
        box + 1 - nx,
        box + 1 + nx,
        box + 1 + nx + nx * ny,
        box + 1 + nx - nx * ny,
        box + 1 - nx + nx * ny,
        box + 1 - nx - nx * ny,
        box - 1 - nx,
        box - 1 + nx, 
        box - 1 + nx + nx * ny,
        box - 1 + nx - nx * ny,
        box - 1 - nx + nx * ny,
        box - 1 - nx - nx * ny
        ]

end

# Utility function
"""
Removes all elements from a dictionary
"""
@inline function clear_dict!(d::Dictionary{Int64, Vector{Int64}})

    for i in eachindex(d)

        resize!(d[i], 0)

    end

end


"""

Updates a dictionary d::Dictionary{Int64, Vector{Int64}} with the indices of cells that are within each box of the spatial domain. 

NOTE - the indices of cells refer to their position in the vector cells::Vector{FCell}.

cells::Vector{FCell} => cells to assign to regions
d::Dictionary{Int64, Vector{Int64}} => dictionary mapping domains to the indices of cells within it
dc::Real => side of boxes
nx::Int64 => number of boxes in the x direction (floor(maximum(x) / dc) + 1)
ny::Int64 => number of boxes in the y direction (floor(maximum(y) / dc) + 1)

"""
function assign_cells!(cells::Vector{FCell}, d::Dictionary{Int64, Vector{Int64}}, dc::Real, nx::Int64, ny::Int64)
    
    for (i, cell) in enumerate(cells)

        if any(isnan, cell.x)

            throw(ArgumentError("NaN value in cell at position $i, cell = $cell."))

        elseif box(cell.x, dc, nx, ny) <= 0

            throw(ArgumentError("Cell $i's position $(cell.x) outside lattice domain."))

        else

            box_index = box(cell.x, dc, nx, ny)

            if box_index > length(d)

                throw(ArgumentError("Cell $i's position $(cell.x) outside lattice domain."))

            end

            push!(d[box_index], i)

        end

    end

end

"""

Clears and re-initialises a dictionary object lattice::Dictionary{Int64, Vector{Int64}} with the indices of cells that are within each box of the spatial domain. 

X::Vector{FCell} => cells to assign to regions
lattice::Dictionary{Int64, Vector{Int64}} => dictionary mapping domains to the indices of cells within it
dc::Real => side of boxes

"""
@inline function update_lattice!(X::Vector{FCell}, lattice::Dictionary{Int64, Vector{Int64}}, dc::Real, nx::Int64, ny::Int64)
    clear_dict!(lattice)
    assign_cells!(X, lattice, dc, nx, ny)
end

"""

Get indices of cells in the ≤26 surrounding cubes adjacent to the cube that a cell is in, and the cells within that cube itself.\n\n

NOTE - the result of this function will contain the index for argument cell itself. This must be accounted for or excluded before use, or it may invoke numerical errors.

cell::FCell => cell of interest
d::Dictionary{Int64, Vector{Int64}} => dictionary object mapping each box dividing the spatial domain to a vector of positions in a vector that contains all the cells in the tissue
tissue::TissueParams => parameters object

"""
function get_candidate_neighbour_cells(cell::FCell, d::Dictionary{Int64, Vector{Int64}}, tissue::TissueParams)

    nx = Int64(floor(tissue.Lx / tissue.lattice_side) + 1)

    ny = Int64(floor(tissue.Ly / tissue.lattice_side) + 1)

    nz = Int64(floor(tissue.Lz / tissue.lattice_side) + 1)

    box_index = box(cell.x, tissue.lattice_side, nx, ny) # 0 allocs

    neighbours = neighbouring_boxes(box_index, nx, ny) # 0 allocs

    neighbour_cells = [d[i] for i in neighbours if 1 <= i <= nx * ny * nz] # 87 allocs?

    cands = Vector{Int64}(undef, sum(length, neighbour_cells)) # 1 alloc

    offset = 0

    for i in neighbours

        if 1 <= i <= nx * ny * nz

            for j in eachindex(d[i])

                cands[offset + j] = d[i][j]

            end

            offset += length(d[i])

        end

    end

    cands

end
