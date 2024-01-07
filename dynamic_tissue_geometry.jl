# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Code for changing the geometry of the PSM dynamically according to parameters in Thomson et al. 2021

include("cell_structs.jl")

# Parameterised equation describing number of somites over time, from 16ss onwards in Zebrafish
somite(t) = 6. + (t / 24.7) - 0.5001 * exp(0.0049 * (t - 300.))

# reciprocal
time(s) = 24.7 * (s - 6) + 0.3381 * exp(0.1858 * s)

wavefront_velocity(t, uₐ) = uₐ * (1. / 24.7 - 0.00245049 * exp(0.0049 * (t - 300.)))

function wavefront_velocity(t, uₐ, t_shrink)

    if t < t_shrink

        return 0.0

    else

        return wavefront_velocity(t, uₐ)

    end

end

Length(t, uₐ, L₀) = L₀ - uₐ * (somite(t) - 16.)

function Length(t, uₐ, L₀, t_shrink)

    if t <= t_shrink

        return L₀

    else

        return Length(t, uₐ, L₀)

    end

end

radius(t, mᵣ, r₀) = r₀ - mᵣ * (somite(t) - 16.)

function radius(t, mᵣ, r₀, t_shrink)

    if t <= t_shrink

        return r₀

    else

        return radius(t, mᵣ, r₀)

    end

end

density(t, md, d₀) = d₀ + md * (somite(t) - 16.)

function density(t, md, d₀, t_shrink)

    if t <= t_shrink

        return density(t_shrink, md, d₀) 

    else

        return density(t, md, d₀)

    end

end

function xa(t::Real, tissue::TissueParams; start=25.)

    if t < tissue.t_shrink

        return start

    else

        return start + tissue.uₐ * (somite(t) - somite(tissue.t_shrink))

    end

end


function xa(t::Vector{<:Real}, tissue::TissueParams) 

    out = Vector{Float64}(undef, length(t)) 

    for i in eachindex(t)

        out[i] = xa(t[i], tissue)

    end

    out
end

xas(tracks::Vector{Matrix{Float64}}, tissue::TissueParams) = xa([i - 1 for i in 1:length(tracks)], tissue)

effective_volume_cylinder(t::Real, tissue::TissueParams; start=25.) = π * radius(t, tissue.mᵣ, tissue.r₀, tissue.t_shrink) ^ 2 * (tissue.Xc - xa(t, tissue; start=start))

effective_volume_torus(t::Real, tissue::TissueParams) = π^2 * radius(t, tissue.mᵣ, tissue.r₀, tissue.t_shrink) ^ 2 * tissue.R

effective_volume(time::Vector{Float64}, tissue::TissueParams; start=25.) = [2. * effective_volume_cylinder(t, tissue; start=start) + effective_volume_torus(t, tissue) for t in time]

cell_diameter(t, mdc, dc₀) = dc₀ - mdc * (somite(t) - 16.)

function cell_diameter(t, mdc, dc₀, t_shrink)

    if t < t_shrink

        return dc₀

    else

        return cell_diameter(t, mdc, dc₀)

    end

end

function cell_diameter(t, mdc, dc₀, t_shrink, t_stop)

    if t < t_shrink

        return dc₀

    elseif t < t_stop

        return cell_diameter(t, mdc, dc₀)

    else

        return cell_diameter(t_stop, mdc, dc₀)

    end

end
