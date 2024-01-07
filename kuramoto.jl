# Simulation of model of Uriu et al 2021, using the algorithm presented in that paper. 
using Random

# include("process_tracks.jl")
include("cell_structs.jl")
include("box_splitting.jl")
include("density_gradient.jl")

"""
translate phase

Translate phase to [0, 1]
"""
translate_phase(θ::Float64) = (1 + sin(θ)) / 2.0

"""
frequency 

Gets the frequency for a cell at a specific point in space, for a specified frequency gradient.
"""
frequency(x, xa, L, ω, σ, k) = ω * (σ + (1 - σ) * (1 - exp(- k * (x - xa) / L)) / (1 - exp(- k)))

frequency(x, pp::PhaseParams, tp::TissueParams) = frequency(x, tp.xa, tp.L, pp.ωₒ, pp.σ, pp.k)

frequency(x::FCell, pp::PhaseParams, tp::TissueParams) = frequency(x.x[1], tp.xa, tp.L, pp.ωₒ, pp.σ, pp.k)

"""
stochastic 

Returns a value for stochastic term of the euler_maruyama scheme.

Input: 

dt (float) - time step.
Dθ (float) - phase noise intensity for the model.

\nOutput: float
"""
stochastic(dt, Dθ) = sqrt(2 * Dθ) * randn() * sqrt(dt)

stochastic(Dθ) = sqrt(2 * Dθ) * randn()

euclid_distance((xᵢ, yᵢ, zᵢ), (xⱼ, yⱼ, zⱼ)) = sqrt((xⱼ - xᵢ)^2 + (yⱼ - yᵢ)^2 + (zⱼ - zᵢ)^2)

"""Calculate a step of the euler-maruyama numerical scheme to solve for θᵢ according to the kuramoto model."""
function euler_step!(X::Vector{FCell}, dt::Float64, params::PhaseParams, cell_params::TissueParams)
    # note that tracks is the time slice at a specific point in time.
    new_X = Vector{FCell}(undef, length(X))
    for i = 1:cell_params.N_cells
        phase_coupling = 0.0
        num_neighbours = 0.0
        for j = 1:cell_params.N_cells
            if j != i
                # if adjacent (and thus coupled)
                if euclid_distance(X[i].x, X[j].x) <= cell_params.dc #&& X[j].x[1] >= cell_params.xa
                    phase_coupling += params.k0 * sin(X[j].θ - X[i].θ)
                    num_neighbours += 1.0
                end
            end 
        end
        # As we're dealing with FCell objects, want to record the frequency
        stoch_freq = stochastic(params.Dθ)
        if num_neighbours > 0.0
            det_freq = frequency(X[i], params, cell_params) + phase_coupling / num_neighbours
        else
            det_freq = frequency(X[i], params, cell_params)
        end
        θ = X[i].θ + dt * det_freq + stoch_freq * sqrt(dt)
        new_X[i] = FCell(X[i].x, X[i].n, det_freq + stoch_freq, θ, X[i].τ, X[i].trackid)
    end
    X .= new_X
end

"""
Faster (O(nlogn)) method\n\n

Calculate a step of the euler-maruyama numerical scheme to solve for θᵢ according to the kuramoto model.\n\n

X::Vector{FCell} => cells within the tissue\n
dt::Float64 => time-step for the euler-maruyama scheme\n
params::PhaseParams => parameters for the phase model\n
cell_params::TissueParams => parameters for tissue-wide properties\n
lattice::Dictionary{Int64, Vector{Int64}} => object mapping regions of space to the cells within them\n

"""
function euler_step!(X::Vector{FCell}, dt::Float64, params::PhaseParams, cell_params::TissueParams, lattice::Dictionary{Int64, Vector{Int64}})

    new_X = Vector{FCell}(undef, length(X))

    for i = 1:cell_params.N_cells

        phase_coupling = 0.0

        num_neighbours = 0.0

        for j in get_candidate_neighbour_cells(X[i], lattice, cell_params)

            if j != i

                # if adjacent (and thus coupled)

                if euclid_distance(X[i].x, X[j].x) <= cell_params.dc #&& X[j].x[1] >= cell_params.xa

                    phase_coupling += params.k0 * sin(X[j].θ - X[i].θ)
                    num_neighbours += 1.0

                end

            end 

        end

        # As we're dealing with FCell objects, want to record the frequency

        stoch_freq = stochastic(params.Dθ)

        if num_neighbours > 0.0

            det_freq = frequency(X[i], params, cell_params) + phase_coupling / num_neighbours

        else

            det_freq = frequency(X[i], params, cell_params)

        end

        θ = X[i].θ + dt * det_freq + stoch_freq * sqrt(dt)

        new_X[i] = FCell(X[i].x, X[i].n, det_freq + stoch_freq, θ, X[i].τ, X[i].trackid)

    end

    X .= new_X

end


"""
Faster (O(nlogn)) method\n\n

Calculate a step of the euler-maruyama numerical scheme to solve for θᵢ according to the kuramoto model, arresting
θᵢ dynamics when a cell enters the M-phase of the cell cycle.\n\n

X::Vector{FCell} => cells within the tissue\n
dt::Float64 => time-step for the euler-maruyama scheme\n
params::PhaseParams => parameters for the phase model\n
cell_params::TissueParams => parameters for tissue-wide properties\n
lattice::Dictionary{Int64, Vector{Int64}} => object mapping regions of space to the cells within them\n

"""
function euler_step_with_mitosis!(X::Vector{FCell}, dt::Float64, params::PhaseParams, cell_params::TissueParams, lattice::Dictionary{Int64, Vector{Int64}})

    new_X = Vector{FCell}(undef, length(X))

    for i = 1:cell_params.N_cells

        # Hold the phase constant if the cell is in M phase

        if cell_params.TG < X[i].τ <= cell_params.TG + cell_params.TM

            new_X[i] = FCell(X[i].x, X[i].n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid) # keeps phase frozen

        # If not in M phase, update phase according to kuramoto model

        else

            phase_coupling = 0.0
            num_neighbours = 0.0

            for j in get_candidate_neighbour_cells(X[i], lattice, cell_params)

                if j != i

                    # if adjacent (and thus coupled)

                    if euclid_distance(X[i].x, X[j].x) <= cell_params.dc #&& X[j].x[1] >= cell_params.xa

                        phase_coupling += params.k0 * sin(X[j].θ - X[i].θ)
                        num_neighbours += 1.0

                    end

                end 

            end

            dWt = stochastic(params.Dθ)
            
            if num_neighbours > 0.0

                det_freq = frequency(X[i], params, cell_params) + phase_coupling / num_neighbours

            else

                det_freq = frequency(X[i], params, cell_params)

            end

            θ = X[i].θ + dt * det_freq + dWt * sqrt(dt)

            new_X[i] = FCell(X[i].x, X[i].n, det_freq + dWt, θ, X[i].τ, X[i].trackid)

        end

    end

    X .= new_X

end

"""
Faster (O(nlogn)) method\n\n

Calculate a step of the euler-maruyama numerical scheme to solve for θᵢ according to the kuramoto model, arresting
θᵢ dynamics when a cell enters the M-phase of the cell cycle, and disallowing nearby cells j ∈ {j s.t. |xᵢ - xⱼ| ≤ dc} 
to couple to cell i during this.\n\n

X::Vector{FCell} => cells within the tissue\n
dt::Float64 => time-step for the euler-maruyama scheme\n
params::PhaseParams => parameters for the phase model\n
cell_params::TissueParams => parameters for tissue-wide properties\n
lattice::Dictionary{Int64, Vector{Int64}} => object mapping regions of space to the cells within them\n

"""
function euler_step_with_mitosis_no_coupling!(X::Vector{FCell}, dt::Float64, params::PhaseParams, cell_params::TissueParams, lattice::Dictionary{Int64, Vector{Int64}})

    new_X = Vector{FCell}(undef, length(X))

    for i = 1:cell_params.N_cells

        if cell_params.TG < X[i].τ <= cell_params.TG + cell_params.TM

            new_X[i] = FCell(X[i].x, X[i].n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid) # keeps phase frozen

        # If not in M phase, update phase according to kuramoto model

        else

            phase_coupling = 0.0
            num_neighbours = 0.0

            for j in get_candidate_neighbour_cells(X[i], lattice, cell_params)

                # check cell not identical, and not in M phase

                if j != i && !(cell_params.TG < X[j].τ <= cell_params.TG + cell_params.TM)

                    # if adjacent (and thus coupled)

                    if euclid_distance(X[i].x, X[j].x) <= cell_params.dc #&& X[j].x[1] >= cell_params.xa

                        phase_coupling += params.k0 * sin(X[j].θ - X[i].θ)
                        num_neighbours += 1.0

                    end

                end 

            end

            dWt = stochastic(params.Dθ)
            
            if num_neighbours > 0.0

                det_freq = frequency(X[i], params, cell_params) + phase_coupling / num_neighbours

            else

                det_freq = frequency(X[i], params, cell_params)

            end

            θ = X[i].θ + dt * det_freq + dWt * sqrt(dt)

            new_X[i] = FCell(X[i].x, X[i].n, det_freq + dWt, θ, X[i].τ, X[i].trackid)

        end

    end

    X .= new_X

end

"""
Faster (O(nlogn)) method\n\n

Calculate a step of the euler-maruyama numerical scheme to solve for θᵢ according to the kuramoto model.\n\n

X::Vector{FCell} => cells within the tissue\n
dt::Float64 => time-step for the euler-maruyama scheme\n
params::PhaseParams => parameters for the phase model\n
cell_params::TissueParams => parameters for tissue-wide properties\n
diam::Function => function describing the diameters of cells\n
lattice::Dictionary{Int64, Vector{Int64}} => object mapping regions of space to the cells within them\n

"""
function euler_step!(X::Vector{FCell}, dt::Float64, params::PhaseParams, cell_params::TissueParams, diam::Function, lattice::Dictionary{Int64, Vector{Int64}})

    new_X = Vector{FCell}(undef, length(X))

    for i = 1:cell_params.N_cells

        phase_coupling = 0.0

        num_neighbours = 0.0

        for j in get_candidate_neighbour_cells(X[i], lattice, cell_params)

            if j != i

                # if adjacent (and thus coupled)

                if euclid_distance(X[i].x, X[j].x) <= cell_params.dc #&& X[j].x[1] >= cell_params.xa

                    phase_coupling += params.k0 * sin(X[j].θ - X[i].θ)
                    num_neighbours += 1.0

                end

            end 

        end

        # As we're dealing with FCell objects, want to record the frequency

        stoch_freq = stochastic(params.Dθ)

        if num_neighbours > 0.0

            det_freq = frequency(X[i], params, cell_params) + phase_coupling / num_neighbours

        else

            det_freq = frequency(X[i], params, cell_params)

        end

        θ = X[i].θ + dt * det_freq + stoch_freq * sqrt(dt)

        new_X[i] = FCell(X[i].x, X[i].n, det_freq + stoch_freq, θ, X[i].τ, X[i].trackid)

    end

    X .= new_X

end

# Old things
# _________________________________________________________________________

"""DEFUNCT Simulate Kuramoto model on vector of cells, holding phase constant advanced if a cell is in M phase."""
function euler_step_with_accelerated_phase_mitosis!(X::Vector{FCell}, n::Real, dt::Float64, params::PhaseParams, cell_params::TissueParams)
    # note that tracks is the time slice at a specific point in time.
    new_X = Vector{FCell}(undef, length(X))
    for i = 1:cell_params.N_cells
        # Hold the phase constant if the cell is in M phase, accelerating when it first enters M phase
        if cell_params.TG < X[i].τ <= cell_params.TG + dt
            # calculate accelerated phase
            θ = X[i].θ + frequency(X[i], params, cell_params) * cell_params.TM * n
            new_X[i] = FCell(X[i].x, X[i].n, X[i].dθdt, θ, X[i].τ, X[i].trackid)
        elseif cell_params.TG < X[i].τ <= cell_params.TG + cell_params.TM
            new_X[i] = FCell(X[i].x, X[i].n, X[i].dθdt, X[i].θ, X[i].τ, X[i].trackid) # keeps phase frozen
        # If not in M phase, update phase according to kuramoto model
        else
            phase_coupling = 0.0
            num_neighbours = 0.0
            for j = 1:cell_params.N_cells
                if j != i
                    # if adjacent (and thus coupled)
                    if euclid_distance(X[i].x, X[j].x) <= cell_params.dc #&& X[j].x[1] >= cell_params.xa
                        phase_coupling += params.k0 * sin(X[j].θ - X[i].θ)
                        num_neighbours += 1.0
                    end
                end 
            end

            dWt = stochastic(params.Dθ)
            
            if num_neighbours > 0.0
                det_freq = frequency(X[i], params, cell_params) + phase_coupling / num_neighbours
            else
                det_freq = frequency(X[i], params, cell_params)
            end
            θ = X[i].θ + dt * det_freq + dWt * sqrt(dt)
            new_X[i] = FCell(X[i].x, X[i].n, det_freq + dWt, θ, X[i].τ, X[i].trackid)
        end
    end
    X .= new_X
end

function phase_order(X::Vector{Float64})
    z = 0.0 + 0.0im
    for θ in X
        z += exp(θ * im)
    end
    abs(z / length(X))
end

phase_order(X::Vector{FCell}) = phase_order([x.θ for x in X])

"""
Calculate average phase for collection of oscillators.
"""
function avg_phase(X::Vector{Float64})
    z = 0.0 + 0.0im
    for θ in X
        z += exp(θ * im)
    end
    z /= length(X)
    # Now get phase
    x, y = z.re, z.im
    if x == 0
        if y > 0
            φ = π / 2
        elseif y < 0
            φ = 3 * π / 2
        end
    elseif y == 0
        if x > 0
            φ = 0
        elseif x < 0
            φ = π
        end
    end
    # Assign phase based on quadrants
    if x > 0 && y > 0 # first quadrant
        φ = atan(y / x)
    elseif x < 0 && y > 0 # second
        φ = π + atan(y / x) # atan(y / x) < 0 here
    elseif x < 0 && y < 0 # third
        φ = π + atan(y / x) # atan(y / x) > 0 here
    elseif x > 0 && y < 0 # fourth
        φ = 2 * π + atan(y / x) # atan(y / x) < 0 here
    end
    # Now return value, only if defined
    if @isdefined φ
        return φ
    end
end
