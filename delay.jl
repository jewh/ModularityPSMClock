# Author: James Hammond (github: jewh) james.hammond@merton.ox.ac.uk

# This file contains helper functions and functions specifically related to solving a 
# delay differential equation (DDE) model of the segmentation clock. Its structure is different from other files
# in this repository as functions here correspond to a wide range of functionalities which are divided across
# files otherwise.

include("cell_structs.jl")
include("kuramoto.jl")

"""
Given a cell at a specified index in a reference vector, return its position in a query vector. 

Return an error if not found.
"""
function find_old_index(index::Int64, reference::Vector{FCell}, query::Vector{FCell})

    for j ∈ eachindex(query)

        if query[j].trackid == reference[index].trackid

            return j

        end

    end

    throw(ArgumentError("Cell $(reference[index].trackid) not found in query vector."))

end

"""
Check whether or not a cell has at least nτ instances recorded in the history vector. Use this to decide whether
we can couple to a cell or not. Returns boolean true or false.

This is slow.
"""
function has_history(trackid::Real, history::Vector{Vector{FCell}}, phase::PhaseParams)

    for cell in history[phase.nτ]

        if cell.trackid == trackid

            return true

        end

    end

    return false

end

"""
Check which indices in the vector have cells that have a history.
"""
function get_cells_with_history(cells::Vector{FCell}, history::Vector{Vector{FCell}}, phase::PhaseParams, tissue::TissueParams)

    n = 0

    v = Vector{Int64}(undef, tissue.N_cells)

    for (ind, cell) in enumerate(cells[1:tissue.N_cells])

        if has_history(cell.trackid, history, phase)

            n += 1

            v[n] = ind

        end

    end

    v[1:n]

end

"""
Add a new element to the first position of state vector, move each element of the state vector along the vector by one position,
and remove the end element.
"""
function update_state!(new::Union{Int64, Vector{FCell}}, state::Union{Vector{Int64}, Vector{Vector{FCell}}})

    for i ∈ reverse(eachindex(state)[2:end])

        state[i] = state[i-1]

    end

    state[1] = copy(new)

end


"""
Faster (O(nlogn)) method\n\n

Calculate a step of the euler-maruyama numerical scheme to solve for θᵢ according to the kuramoto model,
incorporating a delay in signalling between cells.\n\n

cells::Vector{FCell} => cells within the tissue\n
history::Vector{Vector{FCell}} => previous states of the tissue\n
cells_with_history::Vector{Int64} => cells which have enough history to be interact with.\n
dt::Float64 => time-step for the euler-maruyama scheme\n
phase::PhaseParams => parameters for the phase model\n
tissue::TissueParams => parameters for tissue-wide properties\n
lattice::Dictionary{Int64, Vector{Int64}} => object mapping regions of space to the cells within them\n

"""
function euler_step_with_delay!(
    cells::Vector{FCell},
    history::Vector{Vector{FCell}},
    dt::Float64, 
    phase::PhaseParams, 
    tissue::TissueParams, 
    lattice::Dictionary{Int64, Vector{Int64}}
    )

    # # get which cells we can actually couple to
    # cells_w_history = get_cells_with_history(cells, history, phase, tissue)

    new_cells = Vector{FCell}(undef, length(cells))

    for i = 1:tissue.N_cells

        phase_coupling = 0.0

        num_neighbours = 0.0

        for j in get_candidate_neighbour_cells(cells[i], lattice, tissue)

            if j != i #&& j in cells_w_history

                # if adjacent (and thus coupled)

                if euclid_distance(cells[i].x, cells[j].x) <= tissue.dc #&& X[j].x[1] >= cell_params.xa

                    k = find_old_index(j, cells, history[phase.nτ])

                    phase_coupling += phase.k0 * sin(history[phase.nτ][k].θ - cells[i].θ)
                    num_neighbours += 1.0

                end

            end 

        end

        # As we're dealing with FCell objects, want to record the frequency

        stoch_freq = stochastic(phase.Dθ)

        if num_neighbours > 0.0

            det_freq = frequency(cells[i], phase, tissue) + phase_coupling / num_neighbours

        else

            det_freq = frequency(cells[i], phase, tissue)

        end

        θ = cells[i].θ + dt * det_freq + stoch_freq * sqrt(dt)

        new_cells[i] = FCell(cells[i].x, cells[i].n, det_freq + stoch_freq, θ, cells[i].τ, cells[i].trackid)

    end

    cells .= new_cells

end

"""
Faster (O(nlogn)) method\n\n

Calculate a step of the euler-maruyama numerical scheme to solve for θᵢ according to the kuramoto model,
incorporating a delay in signalling between cells, and mitosis.\n\n

cells::Vector{FCell} => cells within the tissue\n
history::Vector{Vector{FCell}} => previous states of the tissue\n
dt::Float64 => time-step for the euler-maruyama scheme\n
phase::PhaseParams => parameters for the phase model\n
tissue::TissueParams => parameters for tissue-wide properties\n
lattice::Dictionary{Int64, Vector{Int64}} => object mapping regions of space to the cells within them\n

"""
function euler_step_with_mitosis_and_delay!(
    cells::Vector{FCell},
    history::Vector{Vector{FCell}},
    dt::Float64, 
    phase::PhaseParams, 
    tissue::TissueParams, 
    lattice::Dictionary{Int64, Vector{Int64}}
    )

    # # get which cells we can actually couple to
    # cells_w_history = get_cells_with_history(cells, history, phase, tissue)

    new_cells = Vector{FCell}(undef, length(cells))

    for i = 1:tissue.N_cells

        if tissue.TG < cells[i].τ <= tissue.TG + tissue.TM

            new_cells[i] = FCell(cells[i].x, cells[i].n, cells[i].dθdt, cells[i].θ, cells[i].τ, cells[i].trackid) # keeps phase frozen

        # If not in M phase, update phase according to kuramoto model

        else

            phase_coupling = 0.0

            num_neighbours = 0.0

            for j in get_candidate_neighbour_cells(cells[i], lattice, tissue)

                if j != i #&& j in cells_w_history

                    # if adjacent (and thus coupled)

                    if euclid_distance(cells[i].x, cells[j].x) <= tissue.dc #&& X[j].x[1] >= cell_params.xa

                        k = find_old_index(j, cells, history[phase.nτ])

                        phase_coupling += phase.k0 * sin(history[phase.nτ][k].θ - cells[i].θ)
                        num_neighbours += 1.0

                    end

                end 

            end

            # As we're dealing with FCell objects, want to record the frequency

            stoch_freq = stochastic(phase.Dθ)

            if num_neighbours > 0.0

                det_freq = frequency(cells[i], phase, tissue) + phase_coupling / num_neighbours

            else

                det_freq = frequency(cells[i], phase, tissue)

            end

            θ = cells[i].θ + dt * det_freq + stoch_freq * sqrt(dt)

            new_cells[i] = FCell(cells[i].x, cells[i].n, det_freq + stoch_freq, θ, cells[i].τ, cells[i].trackid)

        end

    end

    cells .= new_cells

end


# function for mitosis with delay