# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Functions to read the large .tsv files that are outputs of systematic parameter sweeps, and plot as 2D heatmaps for synchrony or frequency

include("experiments.jl") # code for simulating
include("save_simulations.jl") # load code for saving data as .csv
include("simulation_analysis.jl") # functions for analysis
using LaTeXStrings, CSV, DataFrames, CairoMakie

# need to define how we parse strings into input
function Base.parse(T::Type{Tuple{PhaseParams, TissueParams}}, str::String)
    # args for PhaseParams
    bottom = maximum(findfirst("PhaseParams(", str)) + 1
    top = minimum(findfirst("),", str[bottom:end])) - 2 + bottom
    pargs = [parse(Float64, x) for x in split(str[bottom:top], ",")]
    pargs[end] == floor(pargs[end]) ? pargs[end] = parse(Int64, pargs[end]) : pargs[end] = 0

    # args for TissueParams
    bottom = maximum(findfirst("TissueParams(", str)) + 1
    top = minimum(findfirst(")", str[bottom:end])) - 2 + bottom
    targs = [parse(Float64, x) for x in split(str[bottom:top], ",")]

    # length of TissueParams should be 41 args
    if length(targs) < 41

        num_left = 41 - length(targs)
        # append zeros to the end of it. This does not affect analysis, but is needed to prevent errors.
        append!(targs, [0.0 for _ = 1:num_left])

    end

    (PhaseParams(pargs...), TissueParams(targs...))
end

function Base.parse(T::Type{Vector{Float64}}, s::String)
    s = replace(s, "[" => "")
    s = replace(s, "]" => "")
    s = replace(s, " " => "")
    s = eachsplit(s, ",")
    parse.(Float64, s) 
end

# ugly fix
Base.parse(T::Type{Int64}, x::Float64) = Int64(floor(x))

function parse_line(s::String)

    if occursin("BoundsError", str)

        phase_nan = PhaseParams(NaN, NaN, NaN, NaN, NaN, NaN)

        tissue_nan = TissueParams(
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            0.,
            0.,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN
        )

        return NaN, [NaN], [NaN], (phase_nan, tissue_nan)

    else

        bottom = maximum(findfirst("(", str)) + 1

        vals = split(str[1:bottom-2], "\t")

        Ns = parse(Float64, vals[1])

        Z = parse(Vector{Float64}, string(vals[2]))

        freqs = parse(Vector{Float64}, string(vals[3]))

        phase, tissue = parse(Tuple{PhaseParams, TissueParams}, str[bottom:end-1])

        return Ns, Z, freqs, (phase, tissue)

    end

end

# ugly but handles the case where there's no parsing needed
parse_column(T::Type{Float64}, Vs::Vector{Float64}) = Vs

function parse_column(T::Type{Float64}, Vs::Vector{String})

    out = Vector{Float64}(undef, length(Vs))

    for i in eachindex(Vs)

        if occursin("Error", Vs[i])

            out[i] = NaN

        else

            out[i] = parse(Float64, Vs[i])

        end

    end

    out

end

function parse_column(T::Type{Vector{Float64}}, Vs::Vector{String})

    out = Vector{Vector{Float64}}(undef, length(Vs))

    for i in eachindex(Vs)

        if occursin("Error", Vs[i])

            out[i] = [NaN]

        else

            out[i] = parse(Vector{Float64}, Vs[i])

        end

    end

    out

end

function parse_column(T::Type{Tuple{PhaseParams, TissueParams}}, Vs::Vector{String})

    out = Vector{Tuple{PhaseParams, TissueParams}}(undef, length(Vs))

    for i in eachindex(Vs)

        if occursin("Error", Vs[i])

            phase_nan = PhaseParams(NaN, NaN, NaN, NaN, NaN, 0)

            tissue_nan = TissueParams(
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    0.,
                    0.,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN
                )

            out[i] = phase_nan, tissue_nan

        else

            out[i] = parse(Tuple{PhaseParams, TissueParams}, Vs[i])

        end

    end

    out

end

function file_to_vect(path; header=true)
    V = []
    open(path, "r") do io
        for (i, line) in enumerate(readlines(io))
            if header && (i == 1 || i == 2)
                continue
            end
            push!(V, parse_line(line))
        end
    end
    V
end

function Base.isnan(v::Vector{Float64})

    for i in eachindex(v)

        if isnan(v[i])

            return true

        end

    end

    false

end

Base.isnan(P::PhaseParams) = P.k == NaN
    
Base.isnan(T::TissueParams) = T.dc == NaN

Base.isnan(x::Tuple{PhaseParams, TissueParams}) = isnan(x[1]) || isnan(x[2])

function file_to_df(path; drop_nans=true)

    @assert path[end-3:end] == ".tsv"

    data = DataFrame(CSV.File(path, delim="\t", header=2, skipto=3))

    data[!, :Ns] = parse_column(Float64, data[!, :Ns])

    data[!, :Zs] = parse_column(Vector{Float64}, data[!, :Zs])

    data[!, :freqs] = parse_column(Vector{Float64}, data[!, :freqs])

    data[!, Symbol("(phase, tissue)")] = parse_column(Tuple{PhaseParams, TissueParams}, Vector{String}(data[!, Symbol("(phase, tissue)")]))

    if drop_nans

        data = filter(row->!any(isnan(x) for x in row), data)

    end

    data

end










function get_range_γₐ(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = tissue.γₐ

    end

    x = sort!(collect(Set(vals)))
    @assert length(x) > 1
    step = x[2] - x[1]

    minimum(vals), step, maximum(vals)

end

function get_range_γₚ(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = tissue.γₚ

    end

    x = sort!(collect(Set(vals)))
    @assert length(x) > 1
    step = x[2] - x[1]

    minimum(vals), step, maximum(vals)

end

function get_range_xq(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = tissue.xq

    end

    x = sort!(collect(Set(vals)))
    @assert length(x) > 1
    step = x[2] - x[1]

    minimum(vals), step, maximum(vals)

end

function get_range_κ(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = phase.k0

    end

    x = sort!(collect(Set(vals)))
    @assert length(x) > 1
    step = x[2] - x[1]

    minimum(vals), step, maximum(vals)

end

function get_range_vs(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = tissue.vs

    end

    x = sort!(collect(Set(vals)))

    if length(x) > 1
        
        step = x[2] - x[1]

        return minimum(vals), step, maximum(vals)

    elseif length(x) == 1

        return minimum(vals), 1, maximum(vals)

    else

        return 0, 1, 0

    end

end

function get_range_L(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = tissue.L

    end

    x = sort!(collect(Set(vals)))
    @assert length(x) > 1
    step = x[2] - x[1]

    minimum(vals), step, maximum(vals)

end

function get_range_Xv(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = tissue.Xᵥ

    end

    x = sort!(collect(Set(vals)))
    @assert length(x) > 1
    step = x[2] - x[1]

    minimum(vals), step, maximum(vals)

end

function get_range_h(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = tissue.h

    end

    x = sort!(collect(Set(vals)))
    @assert length(x) > 1
    step = x[2] - x[1]

    minimum(vals), step, maximum(vals)

end

function get_range_nτ(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = phase.nτ

    end

    x = sort!(collect(Set(vals)))
    @assert length(x) > 1
    step = x[2] - x[1]

    minimum(vals), step, maximum(vals)

end

function get_range_TM(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = tissue.TM

    end

    x = sort!(collect(Set(vals)))
    @assert length(x) > 1
    step = x[2] - x[1]

    minimum(vals), step, maximum(vals)

end

function get_range_TG(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = tissue.TG

    end

    x = sort!(collect(Set(vals)))
    @assert length(x) > 1
    step = x[2] - x[1]

    minimum(vals), step, maximum(vals)

end

function get_range_xdiv(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = tissue.x_div

    end

    x = sort!(collect(Set(vals)))
    @assert length(x) > 1
    step = x[2] - x[1]

    minimum(vals), step, maximum(vals)

end

function get_range_ω₀(V)

    vals = Vector{Float64}(undef, length(V))

    for (i, v) in enumerate(V)

        Ns, Z, freqs, (phase, tissue) = v

        vals[i] = phase.ωₒ

    end

    x = sort!(collect(Set(vals)))
    @assert length(x) > 1
    # step = x[2] - x[1]

    # minimum(vals), step, maximum(vals)

    x

end



function matrix_advection_velocity_parameterspace(result; x_ind = 1)

    x = (x_ind - 1) * 11.

    minγₐ, stepγₐ, maxγₐ = get_range_γₐ(result)

    γₐs = minγₐ:stepγₐ:maxγₐ

    minγₚ, stepγₚ, maxγₚ = get_range_γₚ(result)

    γₚs = minγₚ:stepγₚ:maxγₚ

    minxq, stepxq, maxxq = get_range_xq(result)

    xqs = minxq:stepxq:maxxq

    Nout = Array{Float64, 4}(undef, length(xqs), length(γₐs), length(γₚs), 100)

    Zout = Array{Float64, 4}(undef, length(xqs), length(γₐs), length(γₚs), 100)

    freqout = Array{Float64, 4}(undef, length(xqs), length(γₐs), length(γₚs), 100)

    # make array to count instances of a given parameter tuple
    n_arr = Array{Int64, 3}(undef, length(xqs), length(γₐs), length(γₚs))

    n_arr[:, :, :] .= 0.

    for (N, Z, freq, (phase, tissue)) in result

        i = findfirst(x->x==tissue.xq, xqs)
        j = findfirst(x->x==tissue.γₐ, γₐs)
        k = findfirst(x->x==tissue.γₚ, γₚs)

        n_arr[i, j, k] += 1

        Nout[i, j, k, n_arr[i, j, k]] = N
        Zout[i, j, k, n_arr[i, j, k]] = Z[x_ind]
        freqout[i, j, k, n_arr[i, j, k]] = (freq[x_ind] - avgω(x, phase, tissue))

    end

    N = Array{Float64, 3}(undef, length(xqs), length(γₐs), length(γₚs))

    Z = Array{Float64, 3}(undef, length(xqs), length(γₐs), length(γₚs))

    freq = Array{Float64, 3}(undef, length(xqs), length(γₐs), length(γₚs))

    for i in eachindex(xqs)
        for j in eachindex(γₐs)
            for k in eachindex(γₚs)

                N[i, j, k] = mean(Nout[i, j, k, :])
                Z[i, j, k] = median(Zout[i, j, k, :])
                freq[i, j, k] = mean(freqout[i, j, k, :])

            end
        end
    end

    N, Z, freq
    
end

function matrix_kappa_motility_parameterspace(result; x_ind = 1)

    x = (x_ind - 1) * 11.

    minκ, stepκ, maxκ = get_range_κ(result)

    κs = minκ:stepκ:maxκ

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minL, stepL, maxL = get_range_L(result)

    Ls = minL:stepL:maxL

    Nout = Array{Float64, 4}(undef, length(Ls), length(κs), length(vss), 100)

    Zout = Array{Float64, 4}(undef, length(Ls), length(κs), length(vss), 100)

    freqout = Array{Float64, 4}(undef, length(Ls), length(κs), length(vss), 100)

    # make array to count instances of a given parameter tuple
    n_arr = Array{Int64, 3}(undef, length(Ls), length(κs), length(vss))

    n_arr[:, :, :] .= 0

    for (N, Z, freq, (phase, tissue)) in result

        i = findfirst(x->x==tissue.L, Ls)
        j = findfirst(x->x==phase.k0, κs)
        k = findfirst(x->x==tissue.vs, vss)

        n_arr[i, j, k] += 1

        Nout[i, j, k, n_arr[i, j, k]] = N
        Zout[i, j, k, n_arr[i, j, k]] = Z[x_ind]
        freqout[i, j, k, n_arr[i, j, k]] = (freq[x_ind] - avgω(x, phase, tissue))

    end

    N = Array{Float64, 3}(undef, length(Ls), length(κs), length(vss))

    Z = Array{Float64, 3}(undef, length(Ls), length(κs), length(vss))

    freq = Array{Float64, 3}(undef, length(Ls), length(κs), length(vss))

    for i in eachindex(Ls)
        for j in eachindex(κs)
            for k in eachindex(vss)

                N[i, j, k] = mean(Nout[i, j, k, :])
                Z[i, j, k] = median(Zout[i, j, k, :])
                freq[i, j, k] = mean(freqout[i, j, k, :])

            end
        end
    end

    N, Z, freq
    
end

function matrix_kappa_motility_delay_parameterspace(result; x_ind = 1)

    x = (x_ind - 1) * 11.

    minκ, stepκ, maxκ = get_range_κ(result)

    κs = minκ:stepκ:maxκ

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minnτ, stepnτ, maxnτ = get_range_nτ(result)

    nτs = minnτ:stepnτ:maxnτ

    Nout = Array{Float64, 4}(undef, length(nτs), length(κs), length(vss), 100)

    Zout = Array{Float64, 4}(undef, length(nτs), length(κs), length(vss), 100)

    freqout = Array{Float64, 4}(undef, length(nτs), length(κs), length(vss), 100)

    # make array to count instances of a given parameter tuple
    n_arr = Array{Int64, 3}(undef, length(nτs), length(κs), length(vss))

    n_arr[:, :, :] .= 0

    for (N, Z, freq, (phase, tissue)) in result

        i = findfirst(x->x==phase.nτ, nτs)
        j = findfirst(x->x==phase.k0, κs)
        k = findfirst(x->x==tissue.vs, vss)

        n_arr[i, j, k] += 1

        Nout[i, j, k, n_arr[i, j, k]] = N
        Zout[i, j, k, n_arr[i, j, k]] = Z[x_ind]
        freqout[i, j, k, n_arr[i, j, k]] = (freq[x_ind] - avgω(x, phase, tissue))

    end

    N = Array{Float64, 3}(undef, length(nτs), length(κs), length(vss))

    Z = Array{Float64, 3}(undef, length(nτs), length(κs), length(vss))

    freq = Array{Float64, 3}(undef, length(nτs), length(κs), length(vss))

    for i in eachindex(nτs)
        for j in eachindex(κs)
            for k in eachindex(vss)

                N[i, j, k] = mean(Nout[i, j, k, :])
                Z[i, j, k] = median(Zout[i, j, k, :])
                freq[i, j, k] = mean(freqout[i, j, k, :])

            end
        end
    end

    N, Z, freq
    
end


function matrix_motility_parameterspace(result; x_ind = 1)

    x = (x_ind - 1) * 11.

    minxv, stepxv, maxxv = get_range_Xv(result)

    xvs = minxv:stepxv:maxxv

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minh, steph, maxh = get_range_h(result)

    hs = minh:steph:maxh

    Nout = Array{Float64, 4}(undef, length(vss), length(xvs), length(hs), 100)

    Zout = Array{Float64, 4}(undef, length(vss), length(xvs), length(hs), 100)

    freqout = Array{Float64, 4}(undef, length(vss), length(xvs), length(hs), 100)

    # make array to count instances of a given parameter tuple
    n_arr = Array{Int64, 3}(undef, length(vss), length(xvs), length(hs))

    n_arr[:, :, :] .= 0

    for (N, Z, freq, (phase, tissue)) in result

        i = findfirst(x->x==tissue.vs, vss)
        j = findfirst(x->x==tissue.Xᵥ, xvs)
        k = findfirst(x->x==tissue.h, hs)

        n_arr[i, j, k] += 1

        Nout[i, j, k, n_arr[i, j, k]] = N
        Zout[i, j, k, n_arr[i, j, k]] = Z[x_ind]
        freqout[i, j, k, n_arr[i, j, k]] = (freq[x_ind] - avgω(x, phase, tissue))

    end

    N = Array{Float64, 3}(undef, length(vss), length(xvs), length(hs))

    Z = Array{Float64, 3}(undef, length(vss), length(xvs), length(hs))

    freq = Array{Float64, 3}(undef, length(vss), length(xvs), length(hs))

    for i in eachindex(vss)
        for j in eachindex(xvs)
            for k in eachindex(hs)

                N[i, j, k] = mean(Nout[i, j, k, :])
                Z[i, j, k] = median(Zout[i, j, k, :])
                freq[i, j, k] = mean(freqout[i, j, k, :])

            end
        end
    end

    N, Z, freq
    
end


function matrix_x_div_TM_parameterspace(result; x_ind = 1)

    x = (x_ind - 1) * 11.

    minTM, stepTM, maxTM = get_range_TM(result)

    TMs = minTM:stepTM:maxTM

    minxdiv, stepxdiv, maxxdiv = get_range_xdiv(result)

    xdivs = minxdiv:stepxdiv:maxxdiv

    Nout = Array{Float64, 3}(undef, length(TMs), length(xdivs), 100)

    Zout = Array{Float64, 3}(undef, length(TMs), length(xdivs), 100)

    freqout = Array{Float64, 3}(undef, length(TMs), length(xdivs), 100)

    # make array to count instances of a given parameter tuple
    n_arr = Array{Int64, 2}(undef, length(TMs), length(xdivs))

    n_arr[:, :, :] .= 0

    for (N, Z, freq, (phase, tissue)) in result

        i = findfirst(x->x==tissue.TM, TMs)
        j = findfirst(x->x==tissue.x_div, xdivs)

        n_arr[i, j] += 1

        Nout[i, j, n_arr[i, j]] = N
        Zout[i, j, n_arr[i, j]] = Z[x_ind]
        freqout[i, j, n_arr[i, j]] = (freq[x_ind] - avgω(x, phase, tissue))

    end

    N = Array{Float64, 2}(undef, length(TMs), length(xdivs))

    Z = Array{Float64, 2}(undef, length(TMs), length(xdivs))

    freq = Array{Float64, 2}(undef, length(TMs), length(xdivs))

    for i in eachindex(TMs)
        for j in eachindex(xdivs)

            N[i, j] = mean(Nout[i, j, :])
            Z[i, j] = median(Zout[i, j, :])
            freq[i, j] = mean(freqout[i, j, :])

        end
    end

    N, Z, freq
    
end

function IQR(x)

    q1 = quantile(x, 0.25)
    q2 = quantile(x, 0.75)

    q2 - q1

end

function IQRmatrix_x_div_TM_parameterspace(result; x_ind = 1)

    x = (x_ind - 1) * 11.

    minTM, stepTM, maxTM = get_range_TM(result)

    TMs = minTM:stepTM:maxTM

    minxdiv, stepxdiv, maxxdiv = get_range_xdiv(result)

    xdivs = minxdiv:stepxdiv:maxxdiv

    Nout = Array{Float64, 3}(undef, length(TMs), length(xdivs), 100)

    Zout = Array{Float64, 3}(undef, length(TMs), length(xdivs), 100)

    freqout = Array{Float64, 3}(undef, length(TMs), length(xdivs), 100)

    # make array to count instances of a given parameter tuple
    n_arr = Array{Int64, 2}(undef, length(TMs), length(xdivs))

    n_arr[:, :, :] .= 0

    for (N, Z, freq, (phase, tissue)) in result

        i = findfirst(x->x==tissue.TM, TMs)
        j = findfirst(x->x==tissue.x_div, xdivs)

        n_arr[i, j] += 1

        Nout[i, j, n_arr[i, j]] = N
        Zout[i, j, n_arr[i, j]] = Z[x_ind]
        freqout[i, j, n_arr[i, j]] = (freq[x_ind] - avgω(x, phase, tissue))

    end

    N = Array{Float64, 2}(undef, length(TMs), length(xdivs))

    Z = Array{Float64, 2}(undef, length(TMs), length(xdivs))

    freq = Array{Float64, 2}(undef, length(TMs), length(xdivs))

    for i in eachindex(TMs)
        for j in eachindex(xdivs)

            N[i, j] = IQR(Nout[i, j, :])
            Z[i, j] = IQR(Zout[i, j, :])
            freq[i, j] = IQR(freqout[i, j, :])

        end
    end

    N, Z, freq
    
end

function matrix_x_div_TG_parameterspace(result; x_ind = 1)

    x = (x_ind - 1) * 11.

    minTG, stepTG, maxTG = get_range_TG(result)

    TGs = minTG:stepTG:maxTG

    minxdiv, stepxdiv, maxxdiv = get_range_xdiv(result)

    xdivs = minxdiv:stepxdiv:maxxdiv

    Nout = Array{Float64, 3}(undef, length(TGs), length(xdivs), 100)

    Zout = Array{Float64, 3}(undef, length(TGs), length(xdivs), 100)

    freqout = Array{Float64, 3}(undef, length(TGs), length(xdivs), 100)

    # make array to count instances of a given parameter tuple
    n_arr = Array{Int64, 2}(undef, length(TGs), length(xdivs))

    n_arr[:, :, :] .= 0

    for (N, Z, freq, (phase, tissue)) in result

        i = findfirst(x->x==tissue.TG, TGs)
        j = findfirst(x->x==tissue.x_div, xdivs)

        n_arr[i, j] += 1

        Nout[i, j, n_arr[i, j]] = N
        Zout[i, j, n_arr[i, j]] = Z[x_ind]
        freqout[i, j, n_arr[i, j]] = (freq[x_ind] - avgω(x, phase, tissue))

    end

    N = Array{Float64, 2}(undef, length(TGs), length(xdivs))

    Z = Array{Float64, 2}(undef, length(TGs), length(xdivs))

    freq = Array{Float64, 2}(undef, length(TGs), length(xdivs))

    for i in eachindex(TGs)
        for j in eachindex(xdivs)

            N[i, j] = mean(Nout[i, j, :])
            Z[i, j] = median(Zout[i, j, :])
            freq[i, j] = mean(freqout[i, j, :])

        end
    end

    N, Z, freq
    
end

function matrix_ω₀_TM_parameterspace(result; x_ind = 1)

    x = (x_ind - 1) * 11.

    minTM, stepTM, maxTM = get_range_TM(result)

    TMs = minTM:stepTM:maxTM

    ω₀s = get_range_ω₀(result)

    Nout = Array{Float64, 3}(undef, length(TMs), length(ω₀s), 100)

    Zout = Array{Float64, 3}(undef, length(TMs), length(ω₀s), 100)

    freqout = Array{Float64, 3}(undef, length(TMs), length(ω₀s), 100)

    # make array to count instances of a given parameter tuple
    n_arr = Array{Int64, 2}(undef, length(TMs), length(ω₀s))

    n_arr[:, :, :] .= 0

    for (N, Z, freq, (phase, tissue)) in result

        i = findfirst(x->x==tissue.TM, TMs)
        j = findfirst(x->isapprox(x, phase.ωₒ), ω₀s)

        n_arr[i, j] += 1

        Nout[i, j, n_arr[i, j]] = N
        Zout[i, j, n_arr[i, j]] = Z[x_ind]
        freqout[i, j, n_arr[i, j]] = (freq[x_ind] - avgω(x, phase, tissue))

    end

    N = Array{Float64, 2}(undef, length(TMs), length(ω₀s))

    Z = Array{Float64, 2}(undef, length(TMs), length(ω₀s))

    freq = Array{Float64, 2}(undef, length(TMs), length(ω₀s))

    for i in eachindex(TMs)
        for j in eachindex(ω₀s)

            N[i, j] = mean(Nout[i, j, :])
            Z[i, j] = median(Zout[i, j, :])
            freq[i, j] = mean(freqout[i, j, :])

        end
    end

    N, Z, freq
    
end







# function to calculate synchrony plot for advection parameter space
function visualise_advection_parameter_space(x_ind)

    path = raw"data\20221130AdvectionParamSpaceTissueSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_advection_velocity_parameterspace(result; x_ind = x_ind)

    minγₐ, stepγₐ, maxγₐ = get_range_γₐ(result)

    γₐs = minγₐ:stepγₐ:maxγₐ

    minγₚ, stepγₚ, maxγₚ = get_range_γₚ(result)

    γₚs = minγₚ:stepγₚ:maxγₚ

    minxq, stepxq, maxxq = get_range_xq(result)

    xqs = minxq:stepxq:maxxq

    rfig = Figure(resolution=(1250, 3000), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)
    ffig = Figure(resolution=(1250, 3000), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    raxs = [Axis(rfig[i, j], aspect=DataAspect()) for i = 1:11, j = 1:3]
    faxs = [Axis(ffig[i, j], aspect=DataAspect()) for i = 1:11, j = 1:3]

    for j = 1:3

        for i = 1:11

            raxs[i, j].xlabel = "vₐ (μm·min⁻¹)"
            raxs[i, j].ylabel = "vₚ (μm·min⁻¹)"

            faxs[i, j].xlabel = "vₐ (μm·min⁻¹)"
            faxs[i, j].ylabel = "vₚ (μm·min⁻¹)"

        end
    
    end

    raxs[1, 1].title = "Random"
    raxs[1, 2].title = "DP"
    raxs[1, 3].title = "DP + LV"

    faxs[1, 1].title = "Random"
    faxs[1, 2].title = "DP"
    faxs[1, 3].title = "DP + LV"

    rhms = [
        CairoMakie.heatmap!(raxs[i, 1], γₐs, γₚs, resultZ[i, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect_ratio=:equal, 
        xlims=(γₐs[1], γₐs[end])
        ) for i = 1:11]

    fhms = [
        CairoMakie.heatmap!(faxs[i, 1], γₐs, γₚs, resultfreq[i, :, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1), 
        aspect_ratio=:equal, 
        xlims=(γₐs[1], γₐs[end])
        ) for i = 1:11]

    path = raw"data\20221130AdvectionParamSpaceDorsalPosteriorSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_advection_velocity_parameterspace(result; x_ind = x_ind)

    rhms = [
        CairoMakie.heatmap!(raxs[i, 2], γₐs, γₚs, resultZ[i, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect_ratio=:equal, 
        xlims=(γₐs[1], γₐs[end])
        ) for i = 1:11]

    fhms = [
        CairoMakie.heatmap!(faxs[i, 2], γₐs, γₚs, resultfreq[i, :, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1), 
        aspect_ratio=:equal, 
        xlims=(γₐs[1], γₐs[end])
        ) for i = 1:11]

    path = raw"data\20221130AdvectionParamSpacePosteriorLateralSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_advection_velocity_parameterspace(result; x_ind = x_ind)

    rhms = [
        CairoMakie.heatmap!(raxs[i, 3], γₐs, γₚs, resultZ[i, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect_ratio=:equal, 
        xlims=(γₐs[1], γₐs[end])
        ) for i = 1:11]

    fhms = [
        CairoMakie.heatmap!(faxs[i, 3], γₐs, γₚs, resultfreq[i, :, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1), 
        aspect_ratio=:equal, 
        xlims=(γₐs[1], γₐs[end])
        ) for i = 1:11]

    Colorbar(rfig[4:8, 4], rhms[1], label=L"r")
    Colorbar(ffig[4:8, 4], fhms[1], label=L"\Delta\frac{d\theta}{dt}")

    [Label(rfig[i, 1, Left()], "xq = $(xqs[i])", padding=(0, 0, 0, 20)) for i = 1:11]
    [Label(ffig[i, 1, Left()], "xq = $(xqs[i])", padding=(0, 0, 0, 20)) for i = 1:11]

    colgap!(rfig.layout, 1)
    rowgap!(rfig.layout, 10)

    colgap!(ffig.layout, 1)
    rowgap!(ffig.layout, 10)

    save("figures/advection_velocity/synchrony/20221130AdvectionParamSpaceTissueSimulation_synchrony_x$(x_ind).pdf", rfig)
    save("figures/advection_velocity/frequency/20221130AdvectionParamSpaceTissueSimulation_frequency_x$(x_ind).pdf", ffig)

end


# function to calculate synchrony plot for advection parameter space
function visualise_kappa_motility_parameter_space(x_ind)

    path = raw"data\20221212KappaMotilityParamSpaceTissueSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_parameterspace(result; x_ind = x_ind)

    minκ, stepκ, maxκ = get_range_κ(result)

    κs = minκ:stepκ:maxκ

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minL, stepL, maxL = get_range_L(result)

    Ls = minL:stepL:maxL

    fig = Figure(resolution=(1250, 1250), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    axs = [Axis(fig[i, j], aspect=AxisAspect(1)) for i = 1:5, j = 1:3]

    for j = 1:3

        for i = 1:5

            axs[i, j].xlabel = "κ (min⁻¹)"
            axs[i, j].ylabel = "vs (μm·min⁻¹)"

        end
    
    end

    axs[1, 1].title = "Random"
    axs[1, 2].title = "DP"
    axs[1, 3].title = "DP + LV"

    hms = [
        CairoMakie.heatmap!(axs[i, 1], κs, vss, resultZ[i, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end])
        ) for i = 1:5]

    path = raw"data\20221212KappaMotilityParamSpaceDorsalPosteriorSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_parameterspace(result; x_ind = x_ind)

    hms = [CairoMakie.heatmap!(axs[i, 2], κs, vss, resultZ[i, :, :], colormap=:viridis, colorrange=(0, 1), 
    aspect=:equal, 
    xlims=(κs[1], κs[end])) for i = 1:5]

    path = raw"data\20221212KappaMotilityParamSpacePosteriorLateralSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_parameterspace(result; x_ind = x_ind)

    hms = [CairoMakie.heatmap!(axs[i, 3], κs, vss, resultZ[i, :, :], colormap=:viridis, colorrange=(0, 1), 
    aspect=:equal, 
    xlims=(κs[1], κs[end])) for i = 1:5]

    kL = findfirst(x->x==385, Ls)

    for j = 1:3

        text!(axs[kL, j], "+", position = (0.07, 1.0), align = (:center, :center))

    end

    Colorbar(fig[2:4, 4], hms[1], label="r")

    [Label(fig[i, 1, Left()], "L = $(Ls[i]) μm", padding=(0, 0, 0, 20)) for i = 1:5]

    colgap!(fig.layout, 1)
    rowgap!(fig.layout, 10)

    fig

    save("figures/kappa_motility/synchrony/20221212KappaMotilityParamSpaceTissueSimulation_synchrony_x$(x_ind).pdf", fig)

end

# function to calculate synchrony plot for advection parameter space, but just the plots for L = 385
function visualise_kappa_motility_parameter_space()

    path = raw"data\20221212KappaMotilityParamSpaceTissueSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_parameterspace(result; x_ind = 1) # yielding errors

    minκ, stepκ, maxκ = get_range_κ(result)

    κs = minκ:stepκ:maxκ

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minL, stepL, maxL = get_range_L(result)

    Ls = minL:stepL:maxL

    Lind = findfirst(x->x==385., Ls)

    rfig = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    ffig = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    raxs = [Axis(rfig[1, i], aspect=AxisAspect(1)) for i = 1:3]

    faxs = [Axis(ffig[1, i], aspect=AxisAspect(1)) for i = 1:3]

    raxs[1].ylabel = "vs (μm·min⁻¹)"
    faxs[1].ylabel = "vs (μm·min⁻¹)"

    for j = 1:3

        raxs[j].xlabel = "κ (min⁻¹)"
        faxs[j].xlabel = "κ (min⁻¹)"
    
    end

    raxs[1].title = "Random"
    raxs[2].title = "DP"
    raxs[3].title = "DP + LV"

    faxs[1].title = "Random"
    faxs[2].title = "DP"
    faxs[3].title = "DP + LV"

    CairoMakie.heatmap!(raxs[1], κs, vss, resultZ[Lind, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end]),)

    CairoMakie.contour!(raxs[1], κs, vss, resultZ[Lind, :, :], color=:black, levels=0:0.2:1, labels=true, labelsize=20, labelcolor = :black)

    CairoMakie.heatmap!(faxs[1], κs, vss, resultfreq[Lind, :, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end]))

    path = raw"data\20221212KappaMotilityParamSpaceDorsalPosteriorSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_parameterspace(result; x_ind = 1)

    CairoMakie.heatmap!(raxs[2], κs, vss, resultZ[Lind, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end]),
        )

    CairoMakie.contour!(raxs[2], κs, vss, resultZ[Lind, :, :], color=:black, levels=0:0.2:1, labels=true, labelsize=20, labelcolor = :black)

    CairoMakie.heatmap!(faxs[2], κs, vss, resultfreq[Lind, :, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end]))

    path = raw"data\20221212KappaMotilityParamSpacePosteriorLateralSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_parameterspace(result; x_ind = 1)

    rhm = CairoMakie.heatmap!(raxs[3], κs, vss, resultZ[Lind, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end]))

    CairoMakie.contour!(raxs[3], κs, vss, resultZ[Lind, :, :], color=:black, levels=0:0.2:1, labels=true, labelsize=20, labelcolor = :black)

    fhm = CairoMakie.heatmap!(faxs[3], κs, vss, resultfreq[Lind, :, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end]))

    for j = 1:3

        text!(raxs[j], "+", position = (0.07, 1.0), align = (:center, :center), fontsize=30)
        text!(faxs[j], "+", position = (0.07, 1.0), align = (:center, :center), fontsize=30)

    end

    cbar = Colorbar(rfig[1, 4], rhm, label="r")
    cbar = Colorbar(ffig[1, 4], fhm, label="Δdθ/dt (min⁻¹)")
    # cbar.height = Relative(2/3)

    colgap!(rfig.layout, 20)
    rowgap!(rfig.layout, 10)

    rowsize!(rfig.layout, 1, Aspect(2, 1))

    colgap!(ffig.layout, 20)
    rowgap!(ffig.layout, 10)

    rowsize!(ffig.layout, 1, Aspect(2, 1))

    save("figures/kappa_motility/synchrony/20221212KappaMotilityParamSpaceTissueSimulation_synchrony_L385_with_contours.pdf", rfig)
    save("figures/kappa_motility/frequency/20221212KappaMotilityParamSpaceTissueSimulation_frequency_L385.pdf", ffig)

    ffig
end

# function to calculate synchrony plot for advection parameter space
function visualise_motility_parameter_space(x_ind)

    path = raw"data\20230113MotilityParamSpaceTissueSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_motility_parameterspace(result; x_ind = x_ind)

    minxv, stepxv, maxxv = get_range_Xv(result)

    xvs = minxv:stepxv:maxxv

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minh, steph, maxh = get_range_h(result)

    hs = minh:steph:maxh

    fig = Figure(resolution=(1500, 3000), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    axs = [Axis(fig[i, j], aspect=AxisAspect(1)) for i = 1:11, j = 1:3]

    for j = 1:3

        for i = 1:11

            axs[i, j].xlabel = "Xᵥ"
            axs[i, j].ylabel = "h"

        end
    
    end

    axs[1, 1].title = "Random"
    axs[1, 2].title = "DP"
    axs[1, 3].title = "DP + LV"

    hms = [
        CairoMakie.heatmap!(axs[i, 1], xvs, hs, resultZ[i, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(xvs[1], xvs[end])
        ) for i = 1:11]

    path = raw"data\20230113MotilityParamSpaceDorsalPosteriorSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_motility_parameterspace(result; x_ind = x_ind)

    hms = [
        CairoMakie.heatmap!(axs[i, 2], xvs, hs, resultZ[i, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(xvs[1], xvs[end])
        ) for i = 1:11]

    path = raw"data/20230113MotilityParamSpacePosteriorLateralSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_motility_parameterspace(result; x_ind = x_ind)

    hms = [
        CairoMakie.heatmap!(axs[i, 3], xvs, hs, resultZ[i, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(xvs[1], xvs[end])
        ) for i = 1:11]
    
    kvs = findfirst(x->x==1, vss)

    for j = 1:3

        text!(axs[kvs, j], "+", position = (0.4, 3.0), align = (:center, :center))

    end

    Colorbar(fig[4:8, 4], hms[1], label="r")

    [Label(fig[i, 1, Left()], "vs = $(vss[i]) μm⋅min⁻¹", padding=(0, 0, 0, 20)) for i = 1:11]

    colgap!(fig.layout, 1)
    rowgap!(fig.layout, 10)

    fig

    save("figures/motility/synchrony/20230113MotilityParamSpaceTissueSimulation_synchrony_x$(x_ind).pdf", fig)

end

# function to calculate synchrony plot for advection parameter space
function visualise_motility_parameter_space()

    path = raw"data\20230113MotilityParamSpaceTissueSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_motility_parameterspace(result; x_ind = 1)

    minxv, stepxv, maxxv = get_range_Xv(result)

    xvs = minxv:stepxv:maxxv

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minh, steph, maxh = get_range_h(result)

    hs = minh:steph:maxh

    vsind = findfirst(x->x==1., vss)

    rfig = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    ffig = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    raxs = [Axis(rfig[1, i], aspect=AxisAspect(1)) for i = 1:3]

    faxs = [Axis(ffig[1, i], aspect=AxisAspect(1)) for i = 1:3]

    raxs[1].ylabel = "h"
    faxs[1].ylabel = "h"

    for j = 1:3

        raxs[j].xlabel = "Xᵥ"
        faxs[j].xlabel = "Xᵥ"
    
    end

    raxs[1].title = "Random"
    raxs[2].title = "DP"
    raxs[3].title = "DP + LV"

    faxs[1].title = "Random"
    faxs[2].title = "DP"
    faxs[3].title = "DP + LV"

    CairoMakie.heatmap!(raxs[1], xvs, hs, resultZ[vsind, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(xvs[1], xvs[end])
        )

    CairoMakie.heatmap!(faxs[1], xvs, hs, resultfreq[vsind, :, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1), 
        aspect=:equal, 
        xlims=(xvs[1], xvs[end])
        )

    path = raw"data\20230113MotilityParamSpaceDorsalPosteriorSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_motility_parameterspace(result; x_ind = 1)

    CairoMakie.heatmap!(raxs[2], xvs, hs, resultZ[vsind, :, :], 
    colormap=:viridis, 
    colorrange=(0, 1), 
    aspect=:equal, 
    xlims=(xvs[1], xvs[end])
    )

    CairoMakie.heatmap!(faxs[2], xvs, hs, resultfreq[vsind, :, :], 
    colormap=:balance, 
    colorrange=(-0.1, 0.1), 
    aspect=:equal, 
    xlims=(xvs[1], xvs[end])
    )

    path = raw"data/20230113MotilityParamSpacePosteriorLateralSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_motility_parameterspace(result; x_ind = 1)

    rhm = CairoMakie.heatmap!(raxs[3], xvs, hs, resultZ[vsind, :, :], 
    colormap=:viridis, 
    colorrange=(0, 1), 
    aspect=:equal, 
    xlims=(xvs[1], xvs[end])
    )

    fhm = CairoMakie.heatmap!(faxs[3], xvs, hs, resultfreq[vsind, :, :], 
    colormap=:balance, 
    colorrange=(-0.1, 0.1), 
    aspect=:equal, 
    xlims=(xvs[1], xvs[end])
    )

    for j = 1:3

        text!(raxs[j], "+", position = (0.4, 3.0), align = (:center, :center))
        text!(faxs[j], "+", position = (0.4, 3.0), align = (:center, :center))

    end

    cbar = Colorbar(rfig[1, 4], rhm, label=L"r")
    cbar = Colorbar(ffig[1, 4], fhm, label=L"\Delta\frac{d\theta}{dt} \; \text{(min^{-1})}")
    # cbar.height = Relative(2/3)

    colgap!(rfig.layout, 20)
    rowgap!(rfig.layout, 10)

    rowsize!(rfig.layout, 1, Aspect(2, 1))

    colgap!(ffig.layout, 20)
    rowgap!(ffig.layout, 10)

    rowsize!(ffig.layout, 1, Aspect(2, 1))

    save("figures/motility/synchrony/20230113MotilityParamSpaceTissueSimulation_synchrony_vs1.pdf", rfig)
    save("figures/motility/frequency/20230113MotilityParamSpaceTissueSimulation_frequency_vs1.pdf", ffig)

end

# function to calculate synchrony plot for parameter space, but just the plots for L = 385
function visualise_kappa_motility_parameter_space_with_delay()

    path = raw"data\20230125_kappa_motility_delay_tissue_simulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_parameterspace(result; x_ind = 1) # yielding errors

    minκ, stepκ, maxκ = get_range_κ(result)

    κs = minκ:stepκ:maxκ

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minL, stepL, maxL = get_range_L(result)

    Ls = minL:stepL:maxL

    Lind = findfirst(x->x==385., Ls)

    rfig = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    ffig = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    raxs = [Axis(rfig[1, i], aspect=AxisAspect(1)) for i = 1:3]

    faxs = [Axis(ffig[1, i], aspect=AxisAspect(1)) for i = 1:3]

    raxs[1].ylabel = "vs (μm·min⁻¹)"
    faxs[1].ylabel = "vs (μm·min⁻¹)"

    for j = 1:3

        raxs[j].xlabel = "κ (min⁻¹)"
        faxs[j].xlabel = "κ (min⁻¹)"
    
    end

    raxs[1].title = "Random"
    raxs[2].title = "DP"
    raxs[3].title = "DP + LV"

    faxs[1].title = "Random"
    faxs[2].title = "DP"
    faxs[3].title = "DP + LV"

    CairoMakie.heatmap!(raxs[1], κs, vss, resultZ[Lind, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end]))

    CairoMakie.heatmap!(faxs[1], κs, vss, resultfreq[Lind, :, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end]))

    path = raw"data\20230125_kappa_motility_delay_posterior_simulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_parameterspace(result; x_ind = 1)

    println(maximum(resultfreq))

    CairoMakie.heatmap!(raxs[2], κs, vss, resultZ[Lind, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end]))

    CairoMakie.heatmap!(faxs[2], κs, vss, resultfreq[Lind, :, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end]))

    path = raw"data\20230125_kappa_motility_delay_posterior_lateral_simulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_parameterspace(result; x_ind = 1)

    rhm = CairoMakie.heatmap!(raxs[3], κs, vss, resultZ[Lind, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end]))

    fhm = CairoMakie.heatmap!(faxs[3], κs, vss, resultfreq[Lind, :, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end]))

    for j = 1:3

        text!(raxs[j], "+", position = (0.07, 1.0), align = (:center, :center))
        text!(faxs[j], "+", position = (0.07, 1.0), align = (:center, :center))

    end

    cbar = Colorbar(rfig[1, 4], rhm, label="r")
    cbar = Colorbar(ffig[1, 4], fhm, label="dθ/dt (min⁻¹)")
    # cbar.height = Relative(2/3)

    colgap!(rfig.layout, 20)
    rowgap!(rfig.layout, 10)

    rowsize!(rfig.layout, 1, Aspect(2, 1))

    colgap!(ffig.layout, 20)
    rowgap!(ffig.layout, 10)

    rowsize!(ffig.layout, 1, Aspect(2, 1))

    save("figures/kappa_motility_delay/synchrony/20230125_kappamotilitydelaysimulation_synchrony_L385.pdf", rfig)
    save("figures/kappa_motility_delay/frequency/20230125_kappamotilitydelaysimulation_frequency_L385.pdf", ffig)

    ffig
end



# function to plot a square of increasing delay from 0 to 21 mins
function visualise_kappa_motility_parameter_space_changing_delay(x_ind)

    path = raw"data\20230314_kappa_motility_delay_tissue_simulation_changing_delay.tsv"

    path_to_increasing_delay = raw"data\20230418_kappa_motility_delay_tissue_simulation_changing_delay.tsv"

    df = file_to_df(path; drop_nans=true)

    df_inc_del = file_to_df(path_to_increasing_delay; drop_nans=true)

    # remove rows in df_inc_del with nτ = 2101, as this is a duplicate of some of the data in the df dataframe
    function find_duplicate(row)
        phase, tissue = row
        phase.nτ !== 2101
    end

    filter!(Symbol("(phase, tissue)") => find_duplicate, df_inc_del)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df) + nrow(df_inc_del))

    for i = 1:nrow(df)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    for i = 1:nrow(df_inc_del)

        result[i + nrow(df)] = (df_inc_del[i, :Ns], df_inc_del[i, :Zs], df_inc_del[i, :freqs], df_inc_del[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_delay_parameterspace(result; x_ind = 1)

    minκ, stepκ, maxκ = get_range_κ(result)

    κs = minκ:stepκ:maxκ

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minnτ, stepnτ, maxnτ = get_range_nτ(result)

    nτs = minnτ:stepnτ:maxnτ

    fig = Figure(resolution=(1550, 1750), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    # need axis background colour same as figure so we can hide axes without data
    axs = [Axis(fig[i, j], aspect=AxisAspect(1), backgroundcolor = :grey90) for i = 1:7, j = 1:5]

    for j = 1:5

        for i = 1:7

            axs[i, j].xlabel = "κ (min⁻¹)"
            axs[i, j].ylabel = "vs (μm·min⁻¹)"

        end
    
    end

    # axs[1, 1].title = "Random"
    # axs[1, 2].title = "DP"
    # axs[1, 3].title = "DP + LV"

    k = 0

    for i = 1:7

        for j = 1:5

            k += 1

            if k <= 31
    
                # need to put this into global so that we can get a colorbar later
                # hideous thing to have to do, but quick to write
                global hm = CairoMakie.heatmap!(
                    axs[i, j], κs, vss, resultZ[k, :, :], 
                    colormap=:viridis, 
                    colorrange=(0, 1), 
                    aspect=:equal, 
                    xlims=(κs[1], κs[end])
                    )

                text!(axs[i, j], "+", position = (0.07, 1.0), align = (:center, :center))

                text!(axs[i, j], "τ = $(round(nτs[k] / 100.0; digits=0)) mins", position = (0.09, 4.5), align = (:center, :center))

            end

        end

    end

    # knτ = findfirst(x->x==2101, nτs)

    Colorbar(fig[2:5, 6], hm, label="r")

    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 10)

    # hide axes for the final plots without data
    i = 7
    for j = 2:5
        hidedecorations!(axs[i, j])
        hidespines!(axs[i, j])
    end

    save("figures/kappa_motility_delay/synchrony/20230314_20230418_kappa_motility_delay_tissue_simulation_changing_delay_x$(x_ind).pdf", fig)

    fig

end


function visualise_kappa_motility_parameter_space_changing_delay_side_by_side(x_ind)

    path = raw"data\20230418_kappa_motility_delay_tissue_simulation_changing_delay.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_delay_parameterspace(result; x_ind = 1)

    minκ, stepκ, maxκ = get_range_κ(result)

    κs = minκ:stepκ:maxκ

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minnτ, stepnτ, maxnτ = get_range_nτ(result)

    nτs = minnτ:stepnτ:maxnτ

    fig = Figure(resolution=(1050, 2550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    # need axis background colour same as figure so we can hide axes without data
    axs = [Axis(fig[i, j], aspect=AxisAspect(1), backgroundcolor = :grey90) for i = 1:10, j = 1:3]

    for j = 1:3

        for i = 1:10

            axs[i, j].xlabel = "κ (min⁻¹)"
            axs[i, j].ylabel = "vs (μm·min⁻¹)"

        end
    
    end

    axs[1, 1].title = "Random"
    axs[1, 2].title = "DP"
    axs[1, 3].title = "DP + LV"

    for i = 1:10

        global hm = CairoMakie.heatmap!(
            axs[i, 1], κs, vss, resultZ[i, :, :], 
            colormap=:viridis, 
            colorrange=(0, 1), 
            aspect=:equal, 
            xlims=(κs[1], κs[end])
            )

        text!(axs[i, 1], "+", position = (0.07, 1.0), align = (:center, :center))

        text!(axs[i, 1], "τ = $(round(nτs[i] / 100.0; digits=0)) mins", position = (0.09, 4.5), align = (:center, :center))

    end

    path = raw"data\20230418_kappa_motility_delay_posterior_simulation_changing_delay.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_delay_parameterspace(result; x_ind = 1)
    
    for i = 1:10

        global hm = CairoMakie.heatmap!(
            axs[i, 2], κs, vss, resultZ[i, :, :], 
            colormap=:viridis, 
            colorrange=(0, 1), 
            aspect=:equal, 
            xlims=(κs[1], κs[end])
            )

        text!(axs[i, 2], "+", position = (0.07, 1.0), align = (:center, :center))

        text!(axs[i, 2], "τ = $(round(nτs[i] / 100.0; digits=0)) mins", position = (0.09, 4.5), align = (:center, :center))

    end

    path = raw"data\20230418_kappa_motility_delay_posterior_lateral_changing_delay.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_delay_parameterspace(result; x_ind = 1)
    
    for i = 1:10

        global hm = CairoMakie.heatmap!(
            axs[i, 3], κs, vss, resultZ[i, :, :], 
            colormap=:viridis, 
            colorrange=(0, 1), 
            aspect=:equal, 
            xlims=(κs[1], κs[end])
            )

        text!(axs[i, 3], "+", position = (0.07, 1.0), align = (:center, :center))

        text!(axs[i, 3], "τ = $(round(nτs[i] / 100.0; digits=0)) mins", position = (0.09, 4.5), align = (:center, :center))

    end

    # knτ = findfirst(x->x==2101, nτs)

    Colorbar(fig[3:8, 4], hm, label="r")

    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 10)



    # hide axes for the final two plots without data
    # hidedecorations!(axs[6, 3])
    # hidedecorations!(axs[6, 4])
    # hidespines!(axs[6, 3])
    # hidespines!(axs[6, 4])

    save("figures/kappa_motility_delay/synchrony/20230418_kappa_motility_delay_tissue_simulation_increasing_delay_x$(x_ind).pdf", fig)

    fig

end

function visualise_kappa_motility_parameter_space_changing_delay_side_by_side(x_ind, nτ)

    path = raw"data\20230418_kappa_motility_delay_tissue_simulation_changing_delay.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_delay_parameterspace(result; x_ind = 1)

    minκ, stepκ, maxκ = get_range_κ(result)

    κs = minκ:stepκ:maxκ

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minnτ, stepnτ, maxnτ = get_range_nτ(result)

    nτs = minnτ:stepnτ:maxnτ

    fig = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    # need axis background colour same as figure so we can hide axes without data
    axs = [Axis(fig[1, j], aspect=AxisAspect(1), backgroundcolor = :grey90) for j = 1:3]

    axs[1].ylabel = "vs (μm·min⁻¹)"

    for j = 1:3

        axs[j].xlabel = "κ (min⁻¹)"
    
    end

    axs[1].title = "Random"
    axs[2].title = "DP"
    axs[3].title = "DP + LV"

    nτ_ind = findfirst(x->x==nτ, nτs)

    hm = CairoMakie.heatmap!(
        axs[1], κs, vss, resultZ[nτ_ind, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end])
        )

    text!(axs[1], "+", position = (0.07, 1.0), align = (:center, :center))

    text!(axs[1], "τ = $(round(nτs[nτ_ind] / 100.0; digits=0)) mins", position = (0.09, 4.5), align = (:center, :center))


    path = raw"data\20230418_kappa_motility_delay_posterior_simulation_changing_delay.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_delay_parameterspace(result; x_ind = 1)


    hm = CairoMakie.heatmap!(
        axs[2], κs, vss, resultZ[nτ_ind, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end])
        )

    text!(axs[2], "+", position = (0.07, 1.0), align = (:center, :center))

    text!(axs[2], "τ = $(round(nτs[nτ_ind] / 100.0; digits=0)) mins", position = (0.09, 4.5), align = (:center, :center))


    path = raw"data\20230418_kappa_motility_delay_posterior_lateral_changing_delay.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_kappa_motility_delay_parameterspace(result; x_ind = 1)
    
    hm = CairoMakie.heatmap!(
        axs[3], κs, vss, resultZ[nτ_ind, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(κs[1], κs[end])
        )

    text!(axs[3], "+", position = (0.07, 1.0), align = (:center, :center))

    text!(axs[3], "τ = $(round(nτs[nτ_ind] / 100.0; digits=0)) mins", position = (0.09, 4.5), align = (:center, :center))

    # knτ = findfirst(x->x==2101, nτs)

    Colorbar(fig[1, 4], hm, label="r")

    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 10)

    rowsize!(fig.layout, 1, Aspect(2, 1))

    save("figures/kappa_motility_delay/synchrony/20230418_kappa_motility_delay_tissue_simulation_increasing_delay_x$(x_ind)_nt$(nτ).pdf", fig)

    fig

end

function visualise_motility_parameter_space_with_delay()

    path = raw"data\20230616_motility_delay_posterior_simulation_changing_delay.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_motility_parameterspace(result; x_ind = 1)

    minxv, stepxv, maxxv = get_range_Xv(result)

    xvs = minxv:stepxv:maxxv

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minh, steph, maxh = get_range_h(result)

    hs = minh:steph:maxh

    vsind = findfirst(x->x==1., vss)

    rfig = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    ffig = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    raxs = [Axis(rfig[1, i], aspect=AxisAspect(1)) for i = 1:3]

    faxs = [Axis(ffig[1, i], aspect=AxisAspect(1)) for i = 1:3]

    raxs[1].ylabel = "h"
    faxs[1].ylabel = "h"

    for j = 1:3

        raxs[j].xlabel = "Xᵥ"
        faxs[j].xlabel = "Xᵥ"
    
    end

    raxs[1].title = "Random"
    raxs[2].title = "DP"
    raxs[3].title = "DP + LV"

    faxs[1].title = "Random"
    faxs[2].title = "DP"
    faxs[3].title = "DP + LV"

    CairoMakie.heatmap!(raxs[1], xvs, hs, resultZ[vsind, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(xvs[1], xvs[end])
        )

    CairoMakie.heatmap!(faxs[1], xvs, hs, resultfreq[vsind, :, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1), 
        aspect=:equal, 
        xlims=(xvs[1], xvs[end])
        )

    path = raw"data\20230616_motility_delay_posterior_simulation_changing_delay.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_motility_parameterspace(result; x_ind = 1)

    CairoMakie.heatmap!(raxs[2], xvs, hs, resultZ[vsind, :, :], 
    colormap=:viridis, 
    colorrange=(0, 1), 
    aspect=:equal, 
    xlims=(xvs[1], xvs[end])
    )

    CairoMakie.heatmap!(faxs[2], xvs, hs, resultfreq[vsind, :, :], 
    colormap=:balance, 
    colorrange=(-0.1, 0.1), 
    aspect=:equal, 
    xlims=(xvs[1], xvs[end])
    )

    path = raw"data\20230616_motility_delay_posterior_lateral_changing_delay.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_motility_parameterspace(result; x_ind = 1)

    rhm = CairoMakie.heatmap!(raxs[3], xvs, hs, resultZ[vsind, :, :], 
    colormap=:viridis, 
    colorrange=(0, 1), 
    aspect=:equal, 
    xlims=(xvs[1], xvs[end])
    )

    fhm = CairoMakie.heatmap!(faxs[3], xvs, hs, resultfreq[vsind, :, :], 
    colormap=:balance, 
    colorrange=(-0.1, 0.1), 
    aspect=:equal, 
    xlims=(xvs[1], xvs[end])
    )

    for j = 1:3

        text!(raxs[j], "+", position = (0.4, 3.0), align = (:center, :center))
        text!(faxs[j], "+", position = (0.4, 3.0), align = (:center, :center))

    end

    cbar = Colorbar(rfig[1, 4], rhm, label=L"r")
    cbar = Colorbar(ffig[1, 4], fhm, label=L"\Delta\frac{d\theta}{dt} \; \text{(min^{-1})}")
    # cbar.height = Relative(2/3)

    colgap!(rfig.layout, 20)
    rowgap!(rfig.layout, 10)

    rowsize!(rfig.layout, 1, Aspect(2, 1))

    colgap!(ffig.layout, 20)
    rowgap!(ffig.layout, 10)

    rowsize!(ffig.layout, 1, Aspect(2, 1))

    save("figures/motility/synchrony/20230616_motility_delay_changing_profile_synchrony_vs1.pdf", rfig)
    save("figures/motility/frequency/20230616_motility_delay_changing_profile_frequency_vs1.pdf", ffig)

end


function visualise_motility_parameter_space(x_ind)

    path = raw"data\20230113MotilityParamSpaceTissueSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_motility_parameterspace(result; x_ind = x_ind)

    minxv, stepxv, maxxv = get_range_Xv(result)

    xvs = minxv:stepxv:maxxv

    minvs, stepvs, maxvs = get_range_vs(result)

    vss = minvs:stepvs:maxvs

    minh, steph, maxh = get_range_h(result)

    hs = minh:steph:maxh

    fig = Figure(resolution=(1500, 3000), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    axs = [Axis(fig[i, j], aspect=AxisAspect(1)) for i = 1:11, j = 1:3]

    for j = 1:3

        for i = 1:11

            axs[i, j].xlabel = "Xᵥ"
            axs[i, j].ylabel = "h"

        end
    
    end

    axs[1, 1].title = "Random"
    axs[1, 2].title = "DP"
    axs[1, 3].title = "DP + LV"

    hms = [
        CairoMakie.heatmap!(axs[i, 1], xvs, hs, resultZ[i, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(xvs[1], xvs[end])
        ) for i = 1:11]

    path = raw"data\20230113MotilityParamSpaceDorsalPosteriorSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_motility_parameterspace(result; x_ind = x_ind)

    hms = [
        CairoMakie.heatmap!(axs[i, 2], xvs, hs, resultZ[i, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(xvs[1], xvs[end])
        ) for i = 1:11]

    path = raw"data/20230113MotilityParamSpacePosteriorLateralSimulation.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_motility_parameterspace(result; x_ind = x_ind)

    hms = [
        CairoMakie.heatmap!(axs[i, 3], xvs, hs, resultZ[i, :, :], 
        colormap=:viridis, 
        colorrange=(0, 1), 
        aspect=:equal, 
        xlims=(xvs[1], xvs[end])
        ) for i = 1:11]
    
    kvs = findfirst(x->x==1, vss)

    for j = 1:3

        text!(axs[kvs, j], "+", position = (0.4, 3.0), align = (:center, :center))

    end

    Colorbar(fig[4:8, 4], hms[1], label="r")

    [Label(fig[i, 1, Left()], "vs = $(vss[i]) μm⋅min⁻¹", padding=(0, 0, 0, 20)) for i = 1:11]

    colgap!(fig.layout, 1)
    rowgap!(fig.layout, 10)

    fig

    save("figures/motility/synchrony/20230113MotilityParamSpaceTissueSimulation_synchrony_x$(x_ind).pdf", fig)

end

# function to calculate synchrony plot for advection parameter space
function visualise_mitosis_TG_TM_xdiv(x_ind)

    path = raw"data\20230821_mitosing_simulation_changing_TG.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_x_div_T_parameterspace(result; x_ind = 1)

    minTG, stepTG, maxTG = get_range_TG(result)

    TGs = minTG:stepTG:maxTG

    minxdiv, stepxdiv, maxxdiv = get_range_xdiv(result)

    xdivs = minxdiv:stepxdiv:maxxdiv


    fig1 = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    aspect = length(TGs) / length(xdivs)

    ax1 = [Axis(fig1[1, i], aspect=AxisAspect(aspect)) for i = 1:2]

    ax1[1].ylabel = L"x_{div} \; (\mu m)"
    ax1[2].ylabel = L"x_{div} \; (\mu m)"

    ax1[1].xlabel = L"T_{G} \; (min)"
    ax1[2].xlabel = L"T_{G} \; (min)"

    rfig1 = CairoMakie.heatmap!(ax1[1], TGs, xdivs, resultZ[:, :], 
        colormap=:viridis, 
        colorrange=(0, 1)
        )

    ffig1 = CairoMakie.heatmap!(ax1[2], TGs, xdivs, resultfreq[:, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1)
        )

    # Colorbar(fig1[1, 3], rfig1; label=L"r")
    # Colorbar(fig1[1, 5], ffig1; label=L"\Delta\frac{d\theta}{dt}")

    save("figures/20230821_mitosing_simulation_changing_TG_and_xdiv.svg", fig1)

    path = raw"data\20230823_mitosing_simulation_changing_TM.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_x_div_TM_parameterspace(result; x_ind = 1)

    minTM, stepTM, maxTM = get_range_TM(result)

    TMs = minTM:stepTM:maxTM

    minxdiv, stepxdiv, maxxdiv = get_range_xdiv(result)

    xdivs = minxdiv:stepxdiv:maxxdiv

    fig2 = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    aspect = length(TMs) / length(xdivs)

    ax2 = [Axis(fig2[1, i], aspect=AxisAspect(aspect)) for i = 1:2]

    ax2[1].ylabel = L"x_{div} \; (\mu m)"
    ax2[2].ylabel = L"x_{div} \; (\mu m)"

    ax2[1].xlabel = L"T_{M} \; (min)"
    ax2[2].xlabel = L"T_{M} \; (min)"

    rfig2 = CairoMakie.heatmap!(ax2[1], TMs, xdivs, resultZ[:, :], 
    colormap=:viridis, 
    colorrange=(0, 1)
    )

    ffig2 = CairoMakie.heatmap!(ax2[2], TMs, xdivs, resultfreq[:, :], 
    colormap=:balance, 
    colorrange=(-0.1, 0.1)
    )

    Colorbar(fig2[1, 3], rfig2; label=L"r")
    Colorbar(fig2[1, 5], ffig2; label=L"\Delta\frac{d\theta}{dt}")

    # figure 3
    resultN, resultZ, resultfreq = IQRmatrix_x_div_TM_parameterspace(result; x_ind = 1)

    fig3 = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    ax3 = [Axis(fig3[1, i], aspect=AxisAspect(aspect)) for i = 1:2]

    ax3[1].ylabel = L"x_{div} \; (\mu m)"
    ax3[2].ylabel = L"x_{div} \; (\mu m)"

    ax3[1].xlabel = L"T_{M} \; (min)"
    ax3[2].xlabel = L"T_{M} \; (min)"

    rfig3 = CairoMakie.heatmap!(ax3[1], TMs, xdivs, resultZ[:, :], 
    colormap=:viridis, 
    colorrange=(0, 1)
    )

    ffig3 = CairoMakie.heatmap!(ax3[2], TMs, xdivs, resultfreq[:, :], 
    colormap=:balance,
    colorrange=(-0.1, 0.1)
    )

    Colorbar(fig3[1, 3], rfig3; label=L"IQR(r)")
    Colorbar(fig3[1, 5], ffig3; label=L"IQR(\Delta\frac{d\theta}{dt})")

    save("figures/20230823_mitosing_simulation_changing_TM_xdiv.svg", fig2)
    save("figures/20230823_mitosing_simulation_changing_TM_xdiv_IQR.svg", fig3)

end


function visualise_mitosis_ω₀_TM(x_ind)

    path = raw"data\20230829_mitosing_simulation_changing_TM_and_freq.tsv"

    df = file_to_df(path; drop_nans=true)

    result = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}, Tuple{PhaseParams, TissueParams}}}(undef, nrow(df))

    for i in eachindex(result)

        result[i] = (df[i, :Ns], df[i, :Zs], df[i, :freqs], df[i, Symbol("(phase, tissue)")])

    end

    resultN, resultZ, resultfreq = matrix_ω₀_TM_parameterspace(result; x_ind = 1)

    minTM, stepTM, maxTM = get_range_TM(result)

    TMs = minTM:stepTM:maxTM

    ω₀s = get_range_ω₀(result)

    periods = [2 * pi / f for f in ω₀s]

    fig1 = Figure(resolution=(1250, 550), fontsize=20, font="CMU Serif", backgroundcolor=:grey90)

    aspect = length(TMs) / length(ω₀s)

    ax1 = [Axis(fig1[1, i], aspect=AxisAspect(aspect)) for i = 1:2]

    ax1[1].ylabel = L"2 \pi / \omega_{0} \; (min)"
    ax1[2].ylabel = L"2 \pi / \omega_{0} \; (min)"

    ax1[1].xlabel = L"T_{M} \; (min)"
    ax1[2].xlabel = L"T_{M} \; (min)"

    rfig1 = CairoMakie.heatmap!(ax1[1], TMs, periods, resultZ[:, :], 
        colormap=:viridis, 
        colorrange=(0, 1)
        )

    # add y = 2x
    y1 = [2 * T for T in TMs]
    y2 = [T for T in TMs]
    lines!(ax1[1], TMs, y1, color=:black, linestyle=:dash)
    lines!(ax1[1], TMs, y2, color=:black, linestyle=:dash)

    xlims!(ax1[1], minimum(TMs), maximum(TMs))
    ylims!(ax1[1], minimum(periods), maximum(periods))

    ffig1 = CairoMakie.heatmap!(ax1[2], TMs, periods, resultfreq[:, :], 
        colormap=:balance, 
        colorrange=(-0.1, 0.1)
        )

    Colorbar(fig1[1, 3], rfig1; label=L"r")
    Colorbar(fig1[1, 5], ffig1; label=L"\Delta\frac{d\theta}{dt}")

    save("figures/20230829_mitosing_simulation_changing_TM_freq.svg", fig1)

end
