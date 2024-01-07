# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# This file contains functions for analysing the outputs of simulations.

using LinearAlgebra

# mean value of ω
"""Integral of ω(x) across the interval [x, x+dc], where dc diameter of a cell, divided by dc."""
function avgω(x, xa, L, σ, ωₒ, dc, k) 
    avgω = ωₒ * σ * dc
    avgω += ωₒ * (1 - σ) * dc / (1 - exp(-k))
    avgω += ωₒ * L * (1 - σ) / (k * (1 - exp(-k))) * (exp(-k * (x + dc - xa) / L) - exp(-k * (x - xa) / L))
    avgω / dc
end

"""Integral of ω(x) across the interval [x, x+dc], where dc diameter of a cell, divided by dc."""
avgω(x, P::PhaseParams, T::TissueParams) = avgω(x, T.xa, T.L, P.σ, P.ωₒ, T.dc, P.k)

"""
Function to get the mean frequency, excluding the cells that are in M phase.
"""
function mean_frequency(M::Matrix{Float64}, tissue::TissueParams; excluding_mitosing=false)

    if excluding_mitosing

        avg = 0.
        n = 0.

        for row in eachrow(M)

            if tissue.TG <= row[8] < tissue.TG + tissue.TM

                continue

            else

                n += 1.
                avg += row[7]

            end

        end

        if n == 0.

            return avg

        else

            return avg / n

        end

    else

        return mean(M[:, 7])

    end

end

"""
Get slices of the PSM along the AP axis.

Argument side specifies whether taking a domain from the left or right hand side of the psm.
"""
function thin_domain(x, tracks::Matrix{Float64}, m, Δx, r, R; side="l")
    cells_in_domain = Int64[]
    for i = 1:size(tracks, 1)
        if side == "l"
            if x + (m-1) * Δx <= tracks[i, 1] <= x + m * Δx && 0 <= tracks[i, 2] <= 2 * r && 0 <= tracks[i, 3] <= 2 * r
                push!(cells_in_domain, i)
            end
        elseif side == "r"
            if x + (m-1) * Δx <= tracks[i, 1] <= x + m * Δx && 2 * R <= tracks[i, 2] <= 2 * (r + R) && 0 <= tracks[i, 3] <= 2 * r
                push!(cells_in_domain, i)
            end
        end
    end
    domain = Matrix{Float64}(undef, length(cells_in_domain), size(tracks, 2))
    for (j, i) in enumerate(cells_in_domain)
        domain[j, :] .= tracks[i, :]
    end
    domain
end

"""
Get local synchrony around a point in space.
"""
function local_synchrony(x, tracks::Matrix{Float64}, Δx, r, R; side="l")
    synchrony = 0.0
    for m = 1:5
        Ωₘ = thin_domain(x, tracks, m, Δx, r, R, side=side)
        synchrony += phase_order(Ωₘ[:, 6])
    end
    synchrony / 5.0
end

"""
Get local phase around a point in space.
"""
function local_phase(x, tracks::Matrix{Float64}, Δx, r, R; side="l")
    phase = 0.0
    for m = 1:5
        Ωₘ = thin_domain(x, tracks, m, Δx, r, R, side=side)
        phase += avg_phase(Ωₘ[:, 6])
    end
    phase / 5.0
end

"""
Get synchrony along PSM at a specific time.
"""
function AP_synchrony(tracks::Matrix{Float64}, tissue::TissueParams, Δx::Float64, R::Float64, r::Float64)
    synchrony = Vector{Float64}(undef, length(tissue.xa:Δx:tissue.Xc))
    for (i, x) in enumerate(tissue.xa:Δx:tissue.Xc)
        synchrony[i] = local_synchrony(x, tracks, Δx, tissue.r, tissue.R)
    end
    synchrony
end

"""
Get frequency along PSM at a specific time.
"""
function AP_frequency(tracks::Matrix{Float64}, tissue::TissueParams, Δx::Float64)
    frequency = Vector{Float64}(undef, length(tissue.xa:Δx:tissue.Xc))
    for (i, x) in enumerate(tissue.xa:Δx:tissue.Xc)
        Ωₘ = thin_domain(x, tracks, m, Δx, tissue.r, tissue.R)
        frequency[i] = local_synchrony(x, tracks, Δx, tissue.r, tissue.R)
    end
    synchrony
end

function anterior_frequency(M::Matrix{Float64}, max_r::Float64, max_R::Float64, Δx::Float64, xa::Float64, tissue::TissueParams)

    mean_frequency(thin_domain(xa, M, 1, Δx, max_r, max_R), tissue)

end

function anterior_frequency(tracks::Vector{Matrix{Float64}}, max_r::Float64, max_R::Float64, Δx::Float64, xas::Vector{Float64}, tissue::TissueParams)

    @assert length(xas) == length(tracks)

    out = Vector{Float64}(undef, length(tracks))

    for (i, M) in enumerate(tracks)

        out[i] = anterior_frequency(M, max_r, max_R, Δx, xas[i], tissue)

    end

    out

end

function anterior_synchrony(tracks::Vector{Matrix{Float64}}, max_r::Float64, max_R::Float64, Δx::Float64, xas::Vector{Float64})

    @assert length(xas) == length(tracks)

    out = Vector{Float64}(undef, length(tracks))

    for (i, M) in enumerate(tracks)

        out[i] = local_synchrony(xas[i], M, Δx, max_r, max_R)

    end

    out

end

function anterior_density(tracks::Vector{Matrix{Float64}}, max_r::Float64, max_R::Float64, Δx::Float64, xas::Vector{Float64})

    @assert length(xas) == length(tracks)

    out = Vector{Float64}(undef, length(tracks))

    for (i, M) in enumerate(tracks)

        out[i] = size(thin_domain(xas[i], M, 1, Δx, max_r, max_R), 1) / (π * max_r ^ 2 * Δx)

    end

    out

end

in_psm(x::Float64, xa::Float64) = x >= xa

function in_psm(M::Matrix{Float64}, xa::Float64)

    out = 0

    for V in eachrow(M)

        if in_psm(V[1], xa)

            out += 1

        end

    end

    out

end

function cells_in_psm_in_time(tracks::Vector{Matrix{Float64}}, xas::Vector{Float64})

    @assert length(xas) == length(tracks)

    out = Vector{Float64}(undef, length(tracks))

    for (i, M) in enumerate(tracks)

        out[i] = in_psm(M, xas[i])

    end

    out

end

# function plot_AP_synchrony(tracks::Matrix{Float64}, tissue::TissueParams, Δx::Float64)
#     plot(
#         tissue.xa:Δx:tissue.Xc,
#         AP_synchrony(tracks, tissue, Δx)
#         )
#     ylims!(0, 1)
# end

# function plot_anterior_synchrony(tracks::Array{Float64, 3}, tissue::TissueParams, Δx::Float64, time_course)
#     zs = Vector{Float64}(undef, size(tracks, 3))
#     for t = 1:size(tracks, 3)
#         zs[t] = local_synchrony(tissue.xa, remove_zero_rows(tracks[:, :, t]), Δx, tissue.r, tissue.R)
#     end
#     Plots.plot(
#         time_course,
#         zs,
#         label = false
#         )
#     ylims!(0, 1)
#     xlabel!("Time (mins)")
#     ylabel!("r")
# end

# function plot_anterior_synchrony(tracks::Vector{Matrix{Float64}}, tissue::TissueParams, Δx::Float64, time_course)
#     zs = Vector{Float64}(undef, length(tracks))
#     for t = 1:length(tracks)
#         zs[t] = local_synchrony(tissue.xa, remove_zero_rows(tracks[t]), Δx, tissue.r, tissue.R)
#     end
#     Plots.plot(
#         time_course,
#         zs,
#         label = false
#         )
#     ylims!(0, 1)
#     xlabel!("Time (mins)")
#     ylabel!("r")
# end

# function animate_AP_synchrony(tracks::Array{Float64, 3}, tissue::TissueParams, Δx::Float64, fname::String)
#     anim = @animate for t = 1:size(tracks, 3)
#         Plots.plot(
#             tissue.xa:Δx:tissue.Xc,
#             AP_synchrony(remove_zero_rows(tracks[:, :, t]), tissue, Δx),
#             label=false,
#             legend=:outertopright
#         )
#         ylims!(0, 1)
#         xlabel!("x (μm)")
#         ylabel!("r")
#         annotate!(307, 1., Plots.text("t = $(t-1) mins", :left, 10))
#     end
#     gif(anim, fname)
# end

# function animate_AP_synchrony(tracks::Vector{Matrix{Float64}}, tissue::TissueParams, Δx::Float64, fname::String)
#     anim = @animate for t = 1:length(tracks)
#         Plots.plot(
#             tissue.xa:Δx:tissue.Xc,
#             AP_synchrony(remove_zero_rows(tracks[t]), tissue, Δx),
#             label=false,
#             legend=:outertopright
#         )
#         ylims!(0, 1)
#         xlabel!("x (μm)")
#         ylabel!("r")
#         annotate!(307, 1., Plots.text("t = $(t-1) mins", :left, 10))
#     end
#     gif(anim, fname)
# end

"""
Measure the average phase along the segmented paraxial mesoderm.
"""
function paraxial_phase(tracks::Matrix{Float64}, Δx, xmin, xmax, side::String, tissue::TissueParams)
    pphase = Matrix{Float64}(undef, 2, length(xmin:Δx:xmax))
    for (i, x) in enumerate(xmin:Δx:xmax)
        Ω = thin_domain(x, tracks, 1, tissue.dc, tissue.r, tissue.R, side=side)
        φ = avg_phase(Ω[:, 6])
        pphase[1, i] = x
        pphase[2, i] = φ
    end
    pphase
end

"""
Get segment boundaries, i.e. positions such that φ = 3π/2
"""
function segment_boundaries(tracks::Matrix{Float64}, Δx, xmin, xmax, side::String, tissue::TissueParams)
    phases = paraxial_phase(tracks, Δx, xmin, xmax, side, tissue)
    phases[2, :] = mod2pi.(phases[2, :])
    positions = Float64[]
    for i = 1:length(xmin:Δx:xmax)
        if abs(abs(phases[2, i]) - 3 * π / 2) <= 0.001
            push!(positions, phases[1, i])
        end
    end
    positions
end

function filter_for_phi(arr)
    subs = abs.(abs.(arr) .- 3 * π / 2)
    for φ in arr
        if abs(abs(φ) - 3 * π / 2) == minimum(subs)
            return φ
        end
    end
end

"""
Get the positions of segment boundaries for 
"""
function segment_boundaries(phases::Matrix{Float64}, tissue::TissueParams, ϵ::Float64)
    phases[2, :] = mod2pi.(phases[2, :])
    # Partition PSM into disjoint intervals
    filtered = phases[:,  abs.(abs.(phases[2, :]) .- 3 * π / 2) .<= ϵ]
    # For each interval, find closest value to 3π/2
    positions = Float64[]
    j = 1
    for i = 1:size(filtered, 2)-1
        if filtered[1, i+1] - filtered[1, i] > tissue.dc # find minimum of all phases back to j
            φ = filter_for_phi(filtered[2, j:i])
            x = minimum(phases[1, phases[2, :] .== φ])
            # save x value
            push!(positions, x)
            # update j
            j = i+1
            # check if at the penultimate point in set, in which output final segment
        end
        if i == size(filtered, 2) - 1
            φ = filter_for_phi(filtered[2, j:size(filtered, 2)])
            x = minimum(phases[1, phases[2, :] .== φ])
            # save x value
            push!(positions, x)

        end
    end
    positions
end

"""
Get the number of segments formed in the simulation on a given side of the PSM.
"""
function total_number_of_segments(tracks::Matrix{Float64}, Δx, xmin, xmax, side::String, tissue::TissueParams)
    length(segment_boundaries(tracks, Δx, xmin, xmax, side, tissue)) - 1
    # Note that, if we assume a segment only to be tissue bounded either side by a region of phase = 3π/2, AND that
    # the anterior most region of the tissue has phase = 3π/2, then the total number of segments = number boundaries - 1.
    # If the anterior-most region does not satisfy this, then the number of segments = number boundaries, if the anterior-most
    # segment is permitted to not be bounded by a boundary at its anterior end.
end

"""
Get the segment lengths in order from anterior to posterior, for a given side of the PSM.
"""
function segment_lengths(tracks::Matrix{Float64}, Δx, xmin, xmax, side::String, tissue::TissueParams)
    seg_bs = segment_boundaries(tracks, Δx, xmin, xmax, side, tissue)
    lengths = Vector{Float64}(undef, length(seg_bs)-1)
    for i = 1:length(seg_bs)-1
        lengths[i] = seg_bs[i+1] - seg_bs[i]
    end
    lengths
end

"""
Quantify the asymmetry in segment length between the left and right hand sides of the PSM.
"""
function left_right_asymmetry(LHS_segment_boundaries::Vector{Float64}, RHS_segment_boundaries::Vector{Float64})
    minr = minimum((length(LHS_segment_boundaries), length(RHS_segment_boundaries)))
    lra = 0.
    for i = 1:minr
        lra += abs(RHS_segment_boundaries[i] - LHS_segment_boundaries[i])
    end
    lra /= minr 
    # Now add the uncounted segments (if mismatch in segment number)
    if minr != length(LHS_segment_boundaries)
        for i = minr+1:length(LHS_segment_boundaries)
            lra += abs(LHS_segment_boundaries[i] - LHS_segment_boundaries[i-1])
        end
        lra /= length(LHS_segment_boundaries) - minr
    end
    if minr != length(RHS_segment_boundaries)
        for i = minr+1:length(RHS_segment_boundaries)
            lra += abs(RHS_segment_boundaries[i] - RHS_segment_boundaries[i-1])
        end
        lra /= length(RHS_segment_boundaries) - minr
    end
    # return output
    lra
end

"""
Measures the total time spent by each cell in the PSM.
"""
function cell_tracking(tracks, tissue::TissueParams)
    new_tracks = copy(tracks)
    for i = 1:size(tracks, 1)
        for t = 1:size(tracks, 3)
            if tracks[i, 1, t] < tissue.xa # if in segmented mesoderm
                new_tracks[i, 6, t] = tracks[i, 4, t]
                break
            end
        end
    end
    new_tracks
end

"""
Convert a list of matrices to a 3D array.

WARNING: If the matrices are of differing size then there will be undef values in the final array.
"""
function list_of_matrices_to_arr(list_of_matrices::Vector{Matrix{Float64}})
    dim1 = maximum([size(M, 1) for M in list_of_matrices])
    dim2 = maximum([size(M, 2) for M in list_of_matrices])
    dim3 = length(list_of_matrices)
    arr = Array{Float64, 3}(undef, dim1, dim2, dim3)
    for (k, M) in enumerate(list_of_matrices)
        for i = 1:size(M, 1)
            for j = 1:size(M, 2)
                arr[i, j, k] = M[i, j]
            end
        end
    end
    arr
end

function small_region_in_space_and_time(x, tracks::Array{Float64, 3}, tissue::TissueParams, m; side="l")
    mats = Vector{Matrix{Float64}}(undef, size(tracks, 3))
    for t = 1:size(tracks, 3)
        mats[t] = thin_domain(x, tracks[:, :, t], m, tissue.dc, tissue.r, tissue.R; side)
    end
    list_of_matrices_to_arr(mats)
end

"""
Get a 2 x t matrix of time in mins, and average phase value for a small region in space.

NOTE - this is not average phase per cell, as cells will move through a region in space over time. Rather this is the average phase for a
small domain along the AP axis, over time.
"""
function average_phase_in_time(x, tracks::Array{Float64, 3}, tissue::TissueParams, m; side="l")
    area = Vector{Matrix{Float64}}(undef, size(tracks, 3))
    for t = 1:size(tracks, 3)
        area[t] = remove_zero_rows(thin_domain(x, tracks[:, :, t], m, tissue.dc, tissue.r, tissue.R; side))
    end
    # Initialise output
    average_phase_in_time = Matrix{Float64}(undef, 2, length(area))
    for t = 1:length(area)
        average_phase_in_time[1, t] = area[t][1, 4]
        average_phase_in_time[2, t] = avg_phase(area[t][:, 6])
    end
    average_phase_in_time
end

"""Local phase in time."""
function local_phase_in_time(x, tracks::Array{Float64, 3}, tissue::TissueParams; side="l")
    locphase = Matrix{Float64}(undef, 2, size(tracks, 3))
    for t = 1:size(tracks, 3)
        locphase[1, t] = tracks[1, 4, t]
        locphase[2, t] = local_phase(x, tracks[:, :, t], tissue.dc, tissue.r, tissue.R; side)
    end
    locphase
end

"""
Get the linear curves that define the saw-edge graph of average phase.

Returns a vector of vectors.
"""
function get_phase_peaks(x, tracks::Array{Float64, 3}, tissue::TissueParams, m; side="l")
    avg_ps = average_phase_in_time(x, tracks, tissue, m; side)
    vect_of_vect = Vector{Vector{Float64}}(undef, size(avg_ps, 2))
    size_vect = 1
    size_vect_vect = 1
    vect = Vector{Float64}(undef, size(avg_ps, 2))
    vect[1] = avg_ps[2, 1]
    for t = 2:size(avg_ps, 2)
        size_vect += 1
        vect[size_vect] = avg_ps[2, t]
        if avg_ps[2, t] - avg_ps[2, t-1] < -π
            vect_of_vect[size_vect_vect] = vect[1:size_vect-1]
            size_vect_vect += 1
            vect = Vector{Float64}(undef, size(avg_ps, 2))
            vect[1] = avg_ps[2, t]
            size_vect = 1
            continue
        end
    end
    vect_of_vect[size_vect_vect] = vect[1:size_vect]
    vect_of_vect[1:size_vect_vect]
end

"""
Get the linear curves that define the saw-edge graph of average phase.

Returns a vector of vectors.
"""
function get_phase_peaks(x, tracks::Array{Float64, 3}, tissue::TissueParams; side="l")
    avg_ps = local_phase_in_time(x, tracks, tissue; side)
    vect_of_vect = Vector{Vector{Float64}}(undef, size(avg_ps, 2))
    size_vect = 1
    size_vect_vect = 1
    vect = Vector{Float64}(undef, size(avg_ps, 2))
    vect[1] = avg_ps[2, 1]
    for t = 2:size(avg_ps, 2)
        size_vect += 1
        vect[size_vect] = avg_ps[2, t]
        if avg_ps[2, t] - avg_ps[2, t-1] < -π
            vect_of_vect[size_vect_vect] = vect[1:size_vect-1]
            size_vect_vect += 1
            vect = Vector{Float64}(undef, size(avg_ps, 2))
            vect[1] = avg_ps[2, t]
            size_vect = 1
            continue
        end
    end
    vect_of_vect[size_vect_vect] = vect[1:size_vect]
    vect_of_vect[1:size_vect_vect]
end


"""
Perform an OLS regression for a vector of independent variable values x and observations y.

Returns vector [c, m], defining the parameters for the equation y = mx + c.
"""
function ordinary_least_squares(x::Vector{Float64}, y::Vector{Float64})
    @assert length(x) == length(y)
    X = Matrix{Float64}(undef, length(y), 2)
    X[:, 1] .= 1
    for (i, x) in enumerate(x)
        X[i, 2] = x
    end
    try
        inv(transpose(X) * X) * transpose(X) * y
    catch err
        println("Error for X = $X.")
        println("determinant = $(det(X)).")
        rethrow(err)
    end
end

"""
Infer the frequency per phase cycle via OLS regression.

Returns a Vector{Float64} object of length Nc, where Nc is the number of cycles for >1 mins.
"""
function infer_frequency(x, tracks::Array{Float64, 3}, tissue::TissueParams, m; side="l")
    phase_cycles = get_phase_peaks(x, tracks, tissue, m; side)
    frequencies = Vector{Float64}(undef, length(phase_cycles))
    nc = 0
    @inbounds for i = 1:length(frequencies)
        if length(phase_cycles[i]) > 1
            nc += 1
            frequencies[nc] = ordinary_least_squares(
                [Float64(x) for x in 0:length(phase_cycles[i])-1], # 'time' array - time defined from starting the cycle to its end
            phase_cycles[i]
            )[2] # [2] as we want just the magnitude.
        end
    end
    if nc == 0
        throw(ErrorException("No phase cycles of length > 1 in the input. Cannot infer frequency."))
    end
    frequencies[1:nc]
end

"""
Infer the frequency per phase cycle via OLS regression.

Returns a Vector{Float64} object of length Nc, where Nc is the number of cycles for >1 mins.
"""
function infer_frequency(x, tracks::Array{Float64, 3}, tissue::TissueParams; side="l")
    phase_cycles = get_phase_peaks(x, tracks, tissue; side)
    frequencies = Vector{Float64}(undef, length(phase_cycles))
    nc = 0
    @inbounds for i = 1:length(frequencies)
        if length(phase_cycles[i]) > 1
            nc += 1
            frequencies[nc] = ordinary_least_squares(
                [Float64(x) for x in 0:length(phase_cycles[i])-1], # 'time' array - time defined from starting the cycle to its end
            phase_cycles[i]
            )[2] # [2] as we want just the magnitude.
        end
    end
    if nc == 0
        throw(ErrorException("No phase cycles of length > 1 in the input. Cannot infer frequency."))
    end
    frequencies[1:nc]
end

# """Plot a line graph of frequency across cycles."""
# function plot_frequency(x, tracks::Array{Float64, 3}, tissue::TissueParams, m; side="l")
#     freqs = infer_frequency(x, tracks, tissue, m; side)
#     Plots.plot(
#         1:length(freqs),
#         freqs,
#         label=false
#     )
#     xlims!(1, length(freqs))
#     xlabel!("Cycle Number")
#     ylabel!("Frequency (min⁻¹)")
#     title!("Frequency per Cycle, at x = $x μm")
# end

# """Plot frequency along the AP axis."""
# function plot_frequency_AP(tracks::Array{Float64, 3}, tissue::TissueParams, m; side="l", label=false)
#     x = tissue.xa:tissue.Xc
#     freqs = Vector{Float64}(undef, length(x))
#     for i = 1:length(freqs)
#         try
#             inf_freqs = infer_frequency(x[i], tracks, tissue, m; side)
#             if length(inf_freqs) == 3
#                 freqs[i] = infer_frequency(x[i], tracks, tissue, m; side)[2]
#             elseif length(inf_freqs) >= 4
#                 freqs[i] = mean(infer_frequency(x[i], tracks, tissue, m; side)[2:end-1])
#             else # if you cannot calculate a value, interpolate it
#                 print("Length inf_freqs = $(length(inf_freqs))")
#                 if i > 1
#                     freqs[i] = freqs[i-1]
#                 else
#                     freqs[i] = 0.
#                 end
#             end
#         catch err
#             println("Error for i = $i")
#             println("Error for freqs[i] = $(freqs[i])")
#             rethrow(err)
#         end
#     end
#     Plots.plot(
#         x,
#         freqs,
#         label=label
#     )
# end

# function plot_frequency_AP!(tracks::Array{Float64, 3}, tissue::TissueParams, m; side="l", label=false)
#     x = tissue.xa:tissue.Xc
#     freqs = Vector{Float64}(undef, length(x))
#     for i = 1:length(freqs)
#         try
#             inf_freqs = infer_frequency(x[i], tracks, tissue, m; side)
#             if length(inf_freqs) == 3
#                 freqs[i] = infer_frequency(x[i], tracks, tissue, m; side)[2]
#             elseif length(inf_freqs) >= 4
#                 freqs[i] = mean(infer_frequency(x[i], tracks, tissue, m; side)[2:end-1])
#             else # if you cannot calculate a value, interpolate it
#                 if i > 1
#                     freqs[i] = freqs[i-1]
#                 else
#                     freqs[i] = 0.
#                 end
#             end
#         catch err
#             println("Error for i = $i")
#             println("Error for freqs[i] = $(freqs[i])")
#             rethrow(err)
#         end
#     end
#     Plots.plot!(
#         x,
#         freqs,
#         label=label
#     )
# end

# """Plot frequency along the AP axis."""
# function plot_frequency_AP(tracks::Array{Float64, 3}, tissue::TissueParams; side="l")
#     x = tissue.xa:tissue.Xc
#     freqs = Vector{Float64}(undef, length(x))
#     for i = 1:length(freqs)
#         try
#             inf_freqs = infer_frequency(x[i], tracks, tissue; side)
#             if length(inf_freqs) == 3
#                 freqs[i] = infer_frequency(x[i], tracks, tissue; side)[2]
#             elseif length(inf_freqs) >= 4
#                 freqs[i] = mean(infer_frequency(x[i], tracks, tissue; side)[2:end-1])
#             else # if you cannot calculate a value, interpolate it
#                 print("Length inf_freqs = $(length(inf_freqs))")
#                 if i > 1
#                     freqs[i] = freqs[i-1]
#                 else
#                     freqs[i] = 0.
#                 end
#             end
#         catch err
#             println("Error for i = $i")
#             println("Error for freqs[i] = $(freqs[i])")
#             rethrow(err)
#         end
#     end
#     Plots.plot(
#         x,
#         freqs,
#         label=false
#     )
# end

# function plot_frequency_AP!(tracks::Array{Float64, 3}, tissue::TissueParams; side="l")
#     x = tissue.xa:tissue.Xc
#     freqs = Vector{Float64}(undef, length(x))
#     for i = 1:length(freqs)
#         try
#             inf_freqs = infer_frequency(x[i], tracks, tissue; side)
#             if length(inf_freqs) == 3
#                 freqs[i] = infer_frequency(x[i], tracks, tissue; side)[2]
#             elseif length(inf_freqs) >= 4
#                 freqs[i] = mean(infer_frequency(x[i], tracks, tissue; side)[2:end-1])
#             else # if you cannot calculate a value, interpolate it
#                 print("Length inf_freqs = $(length(inf_freqs))")
#                 if i > 1
#                     freqs[i] = freqs[i-1]
#                 else
#                     freqs[i] = 0.
#                 end
#             end
#         catch err
#             println("Error for i = $i")
#             println("Error for freqs[i] = $(freqs[i])")
#             rethrow(err)
#         end
#     end
#     Plots.plot!(
#         x,
#         freqs,
#         label=false
#     )
# end

# """Plot the frequency term measured by dθᵢ/dt."""
# function plot_measured_frequency_AP(tracks::Array{Float64, 3}, tissue::TissueParams; side="l")
#     X = tissue.xa:tissue.Xc
#     freqs = Vector{Float64}(undef, length(X))
#     # Get nearby area, then average freq value over Time
#     for (i, x) in enumerate(X)
#         f = 0
#         for t = 1:size(tracks, 3)
#             M = remove_zero_rows(thin_domain(x, tracks[:, :, t], 1, tissue.dc, tissue.r, tissue.R, side=side))
#             f += mean(M[:, 7])
#         end
#         freqs[i] = f / Float64(size(tracks, 3))
#     end

#     # Now plot
#     Plots.plot(
#         X,
#         freqs,
#         label="dθᵢ/dt"
#     )
# end

# """Plot the frequency term measured by dθᵢ/dt."""
# function plot_measured_frequency_AP(tracks::Vector{Matrix{Float64}}, tissue::TissueParams; side="l")
#     X = tissue.xa:tissue.Xc
#     freqs = Vector{Float64}(undef, length(X))
#     # Get nearby area, then average freq value over Time
#     for (i, x) in enumerate(X)
#         f = 0
#         for t = 1:length(tracks)
#             M = remove_zero_rows(thin_domain(x, tracks[t], 1, tissue.dc, tissue.r, tissue.R, side=side))
#             f += mean(M[:, 7])
#         end
#         freqs[i] = f / Float64(length(tracks))
#     end

#     # Now plot
#     Plots.plot(
#         X,
#         freqs,
#         label="dθᵢ/dt"
#     )
# end

# """Plot the frequency term measured by dθᵢ/dt."""
# function plot_measured_frequency_AP(tracks::Matrix{Float64}, tissue::TissueParams; side="l", label="dθᵢ/dt")
#     X = tissue.xa:tissue.Xc
#     freqs = Vector{Float64}(undef, length(X))
#     # Get nearby area, then average freq value over Time
#     for (i, x) in enumerate(X)
#         f = 0
#         M = remove_zero_rows(thin_domain(x, tracks, 1, tissue.dc, tissue.r, tissue.R, side=side))
#         f += mean(M[:, 7])
#         freqs[i] = f / Float64(size(tracks, 3))
#     end

#     # Now plot
#     Plots.plot(
#         X,
#         freqs,
#         label=label,
#         #color=:deepskyblue2
#     )
# end

# """Plot the frequency term measured by dθᵢ/dt."""
# function plot_measured_frequency_AP!(tracks::Array{Float64, 3}, tissue::TissueParams; side="l")
#     X = tissue.xa:tissue.Xc
#     freqs = Vector{Float64}(undef, length(X))
#     # Get nearby area, then average freq value over Time
#     for (i, x) in enumerate(X)
#         f = 0
#         for t = 1:size(tracks, 3)
#             M = remove_zero_rows(thin_domain(x, tracks[:, :, t], 1, tissue.dc, tissue.r, tissue.R, side=side))
#             f += mean(M[:, 7])
#         end
#         freqs[i] = f / Float64(size(tracks, 3))
#     end

#     # Now plot
#     Plots.plot!(
#         X,
#         freqs,
#         label="dθᵢ/dt"
#     )
# end

# """Create an animation for the measured frequency."""
# function animate_frequency_AP(tracks::Array{Float64, 3}, phase::PhaseParams, tissue::TissueParams, fname::String; side="l", label="dθᵢ/dt")
#     anim = @animate for t = 1:size(tracks, 3)
#         plot_measured_frequency_AP(tracks[:, :, t], tissue; side, label)
#         plot!(tissue.xa:tissue.Xc, [avgω(x, phase, tissue) for x in tissue.xa:tissue.Xc], label="̅ω(x)", legend=:outertopright)
#         ylabel!("Frequency (min⁻¹)")
#         xlabel!("x (μm)")
#         ylims!(0.12, 0.22)
#         # add time stamp
#         annotate!(325, 0.195, Plots.text("t = $(t-1) mins", :left, 10))
#     end
#     gif(anim, fname)
# end

# """Create an animation for the measured frequency."""
# function animate_frequency_AP(tracks::Vector{Matrix{Float64}}, phase::PhaseParams, tissue::TissueParams, fname::String; side="l", label="dθᵢ/dt")
#     anim = @animate for t = 1:length(tracks)
#         plot_measured_frequency_AP(tracks[t][:, :], tissue; side, label)
#         plot!(tissue.xa:tissue.Xc, [avgω(x, phase, tissue) for x in tissue.xa:tissue.Xc], label="̅ω(x)", color=:sienna, legend=:outertopright)
#         ylabel!("Frequency (min⁻¹)")
#         xlabel!("x (μm)")
#         ylims!(0.12, 0.22)
#     end
#     gif(anim, fname)
# end

# function plot_mitotic_index(tracks::Vector{Matrix{Float64}}, tissue::TissueParams)
#     mi = Vector{Float64}(undef, length(tracks))
#     for t = 1:length(tracks)
#         n = 0.
#         for i = 1:size(tracks[t], 1)
#             if tissue.TG < tracks[t][i, 8] <= tissue.TG + tissue.TM
#                 n += 1
#             end
#         end
#         mi[t] = n / size(tracks[t], 1)
#     end

#     Plots.plot(1:length(tracks), mi)

# end

"""Calculate, for a given dataset, the difference in phase between a cell and its neighbours after mitosis."""
function get_phase_offset_after_mitosis(X::Vector{Matrix{Float64}}, dt::Float64, tissue::TissueParams)
    _max = maximum([size(M, 1) for M in X])
    vals = Vector{Float64}(undef, _max * length(X))
    k = 0
    for M in X
        for i = 1:size(M, 1)
            # check if cell has just divided or not
            if 0 <= M[i, 8] <= dt         
                phase_offset = 0.
                N = 0.
                for j = 1:size(M, 1)
                    if euclid_distance(SVector{3, Float64}(M[j, 1], M[j, 2], M[j, 3]), SVector{3, Float64}(M[i, 1], M[i, 2], M[i, 3])) <= tissue.dc && i != j
                        phase_offset += M[j, 6] - M[i, 6]
                        N += 1
                    end
                end
                if N != 0
                    push!(vals, phase_offset / N)
                else
                    push!(vals, 0.)
                end
                k += 1
            end
        end
    end
    vals[1:k]
end

"""Calculate, for a given dataset, the difference in phase between a cell and its neighbours after mitosis. Returns
the phase offset and frequency per cell."""
function get_phase_offset_after_mitosis(X::Vector{Matrix{Float64}}, dt::Float64, tissue::TissueParams)
    _max = maximum([size(M, 1) for M in X])
    vals = Vector{Float64}(undef, _max * length(X))
    k = 0
    for M in X
        for i = 1:size(M, 1)
            # check if cell has just divided or not
            if 0 <= M[i, 8] <= dt         
                phase_offset = 0.
                N = 0.
                for j = 1:size(M, 1)
                    if euclid_distance(SVector{3, Float64}(M[j, 1], M[j, 2], M[j, 3]), SVector{3, Float64}(M[i, 1], M[i, 2], M[i, 3])) <= tissue.dc && i != j
                        phase_offset += M[j, 6] - M[i, 6]
                        N += 1
                    end
                end
                if N != 0
                    push!(vals, phase_offset / N)
                else
                    push!(vals, 0.)
                end
                k += 1
            end
        end
    end
    vals[1:k]
end

"""
Function to convert dataframe to vector of matrices.

df_to_mat(df::DataFrame, time_course)

time_course is the iterable containing all the time points at which the data is saved.
"""
function df_to_mat(df::DataFrame, time_course)
    V = Vector{Matrix{Float64}}(undef, length(time_course))
    for (t, time) in enumerate(time_course)
        V[t] = Matrix{Float64}(df[df[!, :Time_Mins] .== time, :])
    end
    V
end

""""""

