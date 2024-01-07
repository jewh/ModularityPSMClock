# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Functions needed to save simulation outputs as .csv files for rendering.

using CSV, DataFrames

include("simulation_analysis.jl")

function is_zero_arr(x::Vector{Float64})
    for i = 1:length(x)
        if x[i] != 0
            return false
        end
    end
    return true
end

function remove_zero_rows(M::Matrix{Float64})
    n = 0
    for i = 1:size(M, 1)
        if is_zero_arr(M[i, 1:end .!= 4]) || M[i, 5] == 0
            n += 1
        end
    end
    # now create array
    A = Matrix{Float64}(undef, size(M, 1) - n, size(M, 2))
    k = 1
    for i = 1:size(M, 1)
        if !is_zero_arr(M[i, :]) && M[i, 5] != 0
            A[k, :] .= M[i, :]
            k += 1
        end
    end
    # return a
    A
end

function arr_to_csv(tracks::Array{Float64, 3}, t_min::Float64, t_max::Float64, fname::String)
    out = DataFrame()
    for t = 1:length(t_min:t_max)
        # remove zero rows, where trackID == 0., no idea where they come from!
        df = DataFrame(remove_zero_rows(tracks[:, :, t]), :auto)
        rename!(df, [ # set the column namess
            Symbol("Position X Reference Frame"),
            Symbol("Position Y Reference Frame"),
            Symbol("Position Z Reference Frame"),
            Symbol("Time_Mins"),
            Symbol("TrackID"),
            Symbol("Phase_sim"),
            # Symbol("CellCyclePhase")
        ])
        # df[!, :Phase_sim] .= (translate_phase(θ) for θ in df.Phase_sim)
        append!(out, df)
    end
    CSV.write(fname, out)
end

function arr_to_csv(tracks::Vector{Matrix{Float64}}, t_min::Float64, t_max::Float64, fname::String)
    out = DataFrame()
    for t = 1:length(t_min:t_max)
        # remove zero rows, where trackID == 0., no idea where they come from!
        df = DataFrame(remove_zero_rows(tracks[:, :, t]), :auto)
        rename!(df, [ # set the column namess
            Symbol("Position X Reference Frame"),
            Symbol("Position Y Reference Frame"),
            Symbol("Position Z Reference Frame"),
            Symbol("Time_Mins"),
            Symbol("TrackID"),
            Symbol("Phase_sim"),
            # Symbol("CellCyclePhase")
        ])
        # df[!, :Phase_sim] .= (translate_phase(θ) for θ in df.Phase_sim)
        append!(out, df)
    end
    CSV.write(fname, out)
end

function tracking_freq_arr_to_csv(tracks::Array{Float64, 3}, t_min::Float64, t_max::Float64, fname::String)
    out = DataFrame()
    for t = 1:length(t_min:t_max)
        # remove zero rows, where trackID == 0., no idea where they come from!
        df = DataFrame(remove_zero_rows(tracks[:, :, t]), :auto)
        rename!(df, [ # set the column namess
            Symbol("Position X Reference Frame"),
            Symbol("Position Y Reference Frame"),
            Symbol("Position Z Reference Frame"),
            Symbol("Time_Mins"),
            Symbol("TrackID"),
            Symbol("Phase_sim"),
            Symbol("dθ/dt")
        ])
        # df[!, :Phase_sim] .= (translate_phase(θ) for θ in df.Phase_sim)
        append!(out, df)
    end
    CSV.write(fname, out)
end

function tracking_freq_arr_to_csv(tracks::Vector{Matrix{Float64}}, fname::String)
    out = DataFrame()
    for t = 1:length(tracks)
        # remove zero rows, where trackID == 0., no idea where they come from!
        df = DataFrame(tracks[t], :auto)
        rename!(df, [ # set the column namess
            Symbol("Position X Reference Frame"),
            Symbol("Position Y Reference Frame"),
            Symbol("Position Z Reference Frame"),
            Symbol("Time_Mins"),
            Symbol("TrackID"),
            Symbol("Phase_sim"),
            Symbol("dθ/dt"),
            Symbol("Cell Cycle Phase")
        ])
        # df[!, :Phase_sim] .= (translate_phase(θ) for θ in df.Phase_sim)
        append!(out, df)
    end
    CSV.write(fname, out)
end

function tracking_freq_arr_to_csv(tracks::Vector{Matrix{Float64}}, t_min::Real, t_max::Real, fname::String)
    out = DataFrame()
    for t = 1:length(tracks)
        # remove zero rows, where trackID == 0., no idea where they come from!
        df = DataFrame(tracks[t], :auto)
        if size(tracks[t], 2) == 8
            rename!(df, [ # set the column namess
                Symbol("Position X Reference Frame"),
                Symbol("Position Y Reference Frame"),
                Symbol("Position Z Reference Frame"),
                Symbol("Time_Mins"),
                Symbol("TrackID"),
                Symbol("Phase_sim"),
                Symbol("dθ/dt"),
                Symbol("Cell Cycle Phase")
            ])
        elseif size(tracks[t], 2) == 7
            rename!(df, [ # set the column namess
            Symbol("Position X Reference Frame"),
            Symbol("Position Y Reference Frame"),
            Symbol("Position Z Reference Frame"),
            Symbol("Time_Mins"),
            Symbol("TrackID"),
            Symbol("Phase_sim"),
            Symbol("dθ/dt"),
        ])
        elseif size(tracks[t], 2) == 6
            rename!(df, [ # set the column namess
            Symbol("Position X Reference Frame"),
            Symbol("Position Y Reference Frame"),
            Symbol("Position Z Reference Frame"),
            Symbol("Time_Mins"),
            Symbol("TrackID"),
            Symbol("Phase_sim"),
        ])
        end
        # df[!, :Phase_sim] .= (translate_phase(θ) for θ in df.Phase_sim)
        append!(out, df)
    end
    CSV.write(fname, out)
end

"""Convert the specified .csv file to a 3-dimensional array in the format expected for 'tracks' arrays."""
function csv_to_arr(path::String)
    df = DataFrame(CSV.File(path))
    # Check column names
    if "Phase" in names(df)
        pname = "Phase"
    elseif "Phase_sim" in names(df)
        pname = "Phase_sim"
    end
    if "Time" in names(df)
        tname = "Time"
    elseif "Time_Mins" in names(df)
        tname = "Time_Mins"
    elseif "Time_Mins_sim" in names(df)
        tname = "Time_Mins_sim"
    end
    if "CellCyclePhase" in names(df)
        cell_cycle = true
    else
        cell_cycle = false
    end
    # Get the number of cells
    n_cells = length(Set(df[!, "TrackID"]))
    # Get the number of time points
    n_times = length(Set(df[!, tname]))
    # Initialise output
    arr = Array{Float64, 3}(undef, n_cells, size(df, 2), n_times)

    for (i, cell) in enumerate(sort!(collect(Set(df[!, "TrackID"]))))
        for (j, time) in enumerate(sort!(collect(Set(df[!, tname]))))

            dti = df[(df[!, "TrackID"] .== cell) .& (df[!, tname] .== time), :]

            if size(dti, 1) > 1
                throw(ErrorException("More than one row in df with TrackID = $cell, $tname = $time"))
            end

            if cell_cycle
                arr[i, :, j] .= (
                    dti[!, "Position X Reference Frame"][1],
                    dti[!, "Position Y Reference Frame"][1],
                    dti[!, "Position Z Reference Frame"][1],
                    dti[!, pname][1],
                    dti[!, tname][1],
                    dti[!, "TrackID"][1],
                    dti[!, "CellCyclePhase"][1]
                )
            else
                arr[i, :, j] .= (
                    dti[!, "Position X Reference Frame"][1],
                    dti[!, "Position Y Reference Frame"][1],
                    dti[!, "Position Z Reference Frame"][1],
                    dti[!, pname][1],
                    dti[!, tname][1],
                    dti[!, "TrackID"][1]
                )
            end
        end
    end
    arr
end

# function to write to .csv

function write_output(v::AbstractArray, fname::AbstractString)

    open(fname, "w") do io
        write(io, "Simulation Type\tSeed\tSynchrony (r)\tMean frequency (/min)")
        for line in v
            if isa(line, Exception)
                write(io, "\n" * string(line) * "\t" *  string(line) * "\t" *  string(line) * "\t" *  string(line))
            elseif isa(line, Tuple{String, Int64, Any, Any})
                write(io, "\n" * string(line[1]) * "\t" *  string(line[2]) * "\t" *  string(line[3]) * "\t" *  string(line[4]))
            end
        end
    end

end

function write_output_more_columns(v::AbstractArray, fname::AbstractString)

    open(fname, "w") do io
        write(io, "Simulation Type\tSeed\tSynchrony (r)\tMean frequency (/min)\tDensity (μm⁻³)\tNo. Cells\txa (μm)\tCell Additions")
        for line in v
            if isa(line, Exception)
                write(io, "\n" * "\t" * string(line)^8)
            else
                write(io, "\n" * string(line[1]) * "\t" *  string(line[2]) * "\t" *  string(line[3]) * "\t" *  string(line[4]) * "\t" *  string(line[5]) * "\t" *  string(line[6]) * "\t" *  string(line[7]) * "\t" *  string(line[8]))
            end
        end
    end

end

"Write the model results to .tsv with the "
function write_output_with_cell_addition(v::AbstractArray, fname::AbstractString)

    open(fname, "w") do io
        write(io, "Simulation Type\tSeed\tSynchrony (r)\tMean frequency (/min)\tadded_cells")
        for line in v
            if isa(line, Exception)
                write(io, "\n" * "\t" * string(line)^8)
            else
                write(io, "\n" * string(line[1]) * "\t" *  string(line[2]) * "\t" *  string(line[3]) * "\t" *  string(line[4]) * "\t" * string(line[5]))
            end
        end
    end

end