# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Functions defining the boundary force for cells within a tissue.

include("cell_structs.jl")


"""
boundary_cylinder

Frictionless boundary force in one of two open-ended cylinders.
"""
function boundary_cylinder((xᵢ, yᵢ, zᵢ), params::TissueParams)
    # Calculate distance from periphery
    if yᵢ < params.Yc
        dy = yᵢ - params.Yc + params.R
    else
        dy = yᵢ - params.Yc - params.R
    end
    dz = zᵢ - params.Zc
    # Now check non-zero (to avoid NaNs)
    if dy == 0 && dz == 0
        return SVector{3, Float64}(0., 0., 0.)
    else
        # Now calculate polar coordinate
        qᵢ = atan(dz / dy)
        # check which quadrant cell is in and update accordingly
        # do nothing if dz > 0 && qᵢ > 0!
        if dz > 0 && qᵢ < 0
            qᵢ += pi
        elseif dz < 0 && qᵢ > 0
            qᵢ += pi
        elseif dz < 0 && qᵢ < 0
            qᵢ += 2 * pi
        end
        δy = abs(params.r * cos(qᵢ) - dy)
        δz = abs(params.r * sin(qᵢ) - dz)
        # Now return boundary force
        return SVector{3, Float64}(
            0.,
            - params.μb * exp(- δy / params.rb) * cos(qᵢ),
            - params.μb * exp(- δz / params.rb) * sin(qᵢ)
        )
    end
end

"""
spherical_coords_torus

Get the spherical coordinates rᵢ, pᵢ, qᵢ, for a point in the toroid tailbud domain.
"""
function spherical_coords_torus((xᵢ, yᵢ, zᵢ), params::TissueParams)
    dx = xᵢ - params.Xc
    dy = yᵢ - params.Yc
    dz = zᵢ - params.Zc
    # calculate spherical coordinates
    if dx > 0
        pᵢ = atan(dy / dx)
    elseif dy > 0
        pᵢ = pi / 2
    else
        pᵢ = - pi / 2
    end
    if abs(cos(pᵢ)) >= abs(sin(pᵢ))
        qᵢ = atan(dz * cos(pᵢ) / (dx - params.R * cos(pᵢ))) # ∈ [-π/2, π/2]
    else
        qᵢ = atan(dz * sin(pᵢ) / (dy - params.R * sin(pᵢ))) # ∈ [-π/2, π/2]
    end
    # Now assign qᵢ to [0, 2π]
    if dz > 0 && qᵢ < 0
        qᵢ += pi 
    elseif dz < 0 && qᵢ > 0
        qᵢ += pi
    elseif dz < 0 && qᵢ < 0
        qᵢ += 2 * pi
    end
    # Now calculate r 
    rᵢ = dz / sin(qᵢ)
    # Now return 
    rᵢ, pᵢ, qᵢ
end

"""
boundary_torus

Frictionless boundary force in half of open-ended torus.
"""
function boundary_torus((xᵢ, yᵢ, zᵢ), params::TissueParams)
    dx = xᵢ - params.Xc
    dy = yᵢ - params.Yc
    dz = zᵢ - params.Zc
    # calculate spherical coordinates
    rᵢ, pᵢ, qᵢ = spherical_coords_torus((xᵢ, yᵢ, zᵢ), params)
    # Now calculate δ
    δx = abs(params.R * cos(pᵢ) + params.r * cos(pᵢ) * cos(qᵢ) - dx)
    δy = abs(params.R * sin(pᵢ) + params.r * sin(pᵢ) * cos(qᵢ) - dy)
    δz = abs(params.r * sin(qᵢ) - dz)
    # Now return
    if params.R * cos(pᵢ) == dx && params.R * sin(pᵢ) == dy && dz == 0
        return SVector{3, Float64}(0., 0., 0.)
    else
        return SVector{3, Float64}(
            - params.μb * exp(- δx / params.rb) * cos(pᵢ) * cos(qᵢ),
            - params.μb * exp(- δy / params.rb) * sin(pᵢ) * cos(qᵢ),
            - params.μb * exp(- δz / params.rb) * sin(qᵢ)
        )
    end
end

"""
Get the boundary force for a horseshoe-shaped PSM, as in Uriu et al. 2021.
"""
function boundary_horseshoe_psm((xᵢ, yᵢ, zᵢ), params::TissueParams)
    if xᵢ <= params.Xc # <= because I imagine there's issues in the torus if x == Xc?
        return boundary_cylinder((xᵢ, yᵢ, zᵢ), params)
    else
        return boundary_torus((xᵢ, yᵢ, zᵢ), params)
    end
end

"""
Check if cell is within torus domain.
"""
function in_torus((xᵢ, yᵢ, zᵢ), params::TissueParams)
    # calculate spherical coordinates
    rᵢ, pᵢ, qᵢ = spherical_coords_torus((xᵢ, yᵢ, zᵢ), params)
    # return 
    rᵢ <= params.r
end

"""
Check if cell is within either cylinder domain.
"""
function in_cylinder((xᵢ, yᵢ, zᵢ), params::TissueParams)
    # Calculate distance from periphery
    if yᵢ < params.Yc
        dy = yᵢ - params.Yc + params.R
    else
        dy = yᵢ - params.Yc - params.R
    end
    dz = zᵢ - params.Zc
    # return
    dy^2 + dz^2 <= params.r^2
end

"""
Check if cell in tissue.
"""
in_tissue((xᵢ, yᵢ, zᵢ), params::TissueParams) = in_cylinder((xᵢ, yᵢ, zᵢ), params) || in_torus((xᵢ, yᵢ, zᵢ), params)


# DEFUNCT STUFF
#______________________________________________________________________________________

# Functions below get the distances from the periphery in x, y, z; for a cube with side 2L centred around (0,0,0).

"""
δ(x)

Distance of cell at position x along a line from the boundary along that line. 
"""
function cδ(x::Float64, L::Float64)
    if x >= 0
        return L - x
    elseif x < 0
        return L + x
    end
end

# Small function to infer the sign of the boundary force
function h(x)
    if x >= 0
        return -1
    else
        return 1
    end
end

"""Get the boundary force for a cube of side 2L centred around (0,0,0)."""
function boundary_cube((xᵢ, yᵢ, zᵢ), params)
    SVector{3, Float64}(
       h(xᵢ) * params.μb * exp(- cδ(xᵢ, params.L) / params.rb),
       h(yᵢ) * params.μb * exp(- cδ(yᵢ, params.L) / params.rb),
       h(zᵢ) * params.μb * exp(- cδ(zᵢ, params.L) / params.rb)
    )
end

"""Get the boundary force for a cube of side L centred around Xc, Yc, Zc."""
function boundary_cube((xᵢ, yᵢ, zᵢ), params)
    SVector{3, Float64}(
       h(xᵢ - params.Xc) * params.μb * exp(- cδ(xᵢ - params.Xc, params.L/2) / params.rb),
       h(yᵢ - params.Yc) * params.μb * exp(- cδ(yᵢ - params.Yc, params.L/2) / params.rb),
       h(zᵢ - params.Zc) * params.μb * exp(- cδ(zᵢ - params.Zc, params.L/2) / params.rb)
    )
end

# Functions below allow boundary force for a curved and spatulate PSM.

diff(x, f, dx) = (f(x + dx) - f(x)) / dx

# Create specific exception class for handling errors due to non-convergence
mutable struct ConvergenceError <: Exception
    x0::Real
    err_msg::String
end

Base.showerror(io::IO, e::ConvergenceError) = print(io, e.err_msg)

"""Find the root of a given function f, for an initial guess x0 and tolerance for convergence ϵ."""
function newton_rhapson(x0::Float64, f::Function, ϵ::Float64, dx::Float64)
    init_guess = x0
    xn = x0 - f(x0) / diff(x0, f, dx)
    n = 1
    while abs(xn - x0) > ϵ && n <= 999
        x0 = xn
        xn = xn - f(xn) / diff(xn, f, dx)
        n += 1
    end
    if n > 999 # did not converge, in this case
        throw(ConvergenceError(init_guess, "Newton-Rhapson scheme did not converge for x0 = $init_guess"))
    else
        return xn
    end
end

"""
The distance of a cell in the internal part of the PSM 'tube' from the periphery.

xᵢ - x-coordinate of the cell i
zᵢ - z-coordinate of the cell i
x - x-coordinate of intersection of normal to PSM centre with the central curve.
φ - rotation function value
xc - tissue parameter determining centre of the PSM along the z-axis.
zc - tissue parameter determining centre of the PSM along the z-axis.
"""
function r(xᵢ, zᵢ, φ, xc, zc)
    sqrt((xᵢ - xc + sin(φ) * zc)^2 + (zᵢ - cos(φ) * zc)^2)
end

function sign(xᵢ, zᵢ, xc, zc, φ)
    if φ <= π/4
        z = zᵢ / cos(φ)
    elseif π/2 >= φ > π/4
        z = (xc - xᵢ) / sin(φ)
    else
        throw(DomainError(φ, "Function sign() only accepts φ ∈ [0, π/2]!"))
    end
    if z >= zc
        return -1
    elseif z < zc
        return 1
    end
end

"""Return a vector for the boundary force within the internal part of the PSM 'tubes'."""
function internal_boundary_force((xᵢ, yᵢ, zᵢ), params::TissueParams)
    try
        φ = atan((params.Xc - xᵢ) / zᵢ)
        rᵢ = r(xᵢ, zᵢ, φ, params.Xc, params.Zc)
        SVector{3, Float64}(
            - params.μb * sign(xᵢ, zᵢ, params.Xc, params.Zc, φ) * exp(-(params.r - rᵢ) / params.rb) * sin(φ),
            0.,
            params.μb * sign(xᵢ, zᵢ, params.Xc, params.Zc, φ) * exp(-(params.r - rᵢ) / params.rb) * cos(φ)
        )
    catch err
        if isa(err, ConvergenceError)
            println("Convergence failure for xᵢ, yᵢ, zᵢ = $xᵢ, $yᵢ, $zᵢ")
            println(err.err_msg)
        elseif isa(err, DomainError)
            println("Failure for (xᵢ - params.Xc) / params.Zc = ($xᵢ - $(params.Xc)) / $(params.Zc)")
            println(err.err_msg)
        else
            rethrow()
        end
    end
end

function isapprox(x, y)
    for (a, b) in zip(x, y)
        if round(a, digits=5) != round(b, digits=5)
            return false
        end
    end
    true
end

function is_internal_boundary_force_working(params::TissueParams)
    xc = params.Xc - sin(π/4) * params.Zc
    zc = cos(π/4) * params.Zc
    if ~isapprox(internal_boundary_force((xc, 0, zc), params), SVector{3, Float64}(0, 0, 0))
        return false
    elseif ~isapprox(internal_boundary_force((tissue.Xc - sin(π/4) * (tissue.Zc + tissue.r), 0, cos(π/4) * (tissue.Zc + tissue.r)), params), SVector{3, Float64}(tissue.μb * sin(π/4), 0, - tissue.μb * cos(π/4)))
        return false
    elseif ~isapprox(internal_boundary_force((tissue.Xc - sin(π/4) * (tissue.Zc - tissue.r), 0, cos(π/4) * (tissue.Zc - tissue.r)), params), SVector{3, Float64}(- tissue.μb * sin(π/4), 0, tissue.μb * cos(π/4)))
        return false
    end
    true
end

# """Get the frictionless boundary force for a flattened tube curved in the x-z plane by a function φ."""
# function left_psm_boundary_force((xᵢ, yᵢ, zᵢ), φ::Function, invφ::Function, d::Function, params::TissueParams)
#     dm = d(params.Xc) # assumes that dm == max{d(x) | x ∈ [0, Xc]}
#     # get the angle of rotation φ and centre of tissue in z axis, Zc
#     try # Errors from the fact that boundary force is written incorrectly.
#         φᵢ = - atan((xᵢ - params.Xc) / zᵢ)
#         x = invφ(φᵢ)
#         xc = tissue.Xc - sin(φᵢ) * tissue.Zc
#         zc =  cos(φᵢ) * params.Zc
#         if yᵢ >= params.r + dm # then at the innermost edge
#             # get the polar coordinates rᵢ and qᵢ
#             rᵢ = sqrt((xᵢ - xc)^2 + (yᵢ - params.r - dm)^2 + (zᵢ - zc)^2)
#             if yᵢ - params.r - dm == 0 && zᵢ >= zc
#                 qᵢ = π/2
#             elseif yᵢ - params.r - dm == 0 && zᵢ < zc
#                 qᵢ = -π/2
#             else
#                 qᵢ = atan((zᵢ - zc) / (yᵢ - params.r - dm))
#             end
#             # return boundary force
#             return SVector{3, Float64}(
#             - params.μb * exp(-(params.r - rᵢ) / params.rb) * sin(qᵢ) * sin(φᵢ),
#             - params.μb * exp(-(params.r - rᵢ) / params.rb) * cos(qᵢ),
#             - params.μb * exp(-(params.r - rᵢ) / params.rb) * sin(qᵢ) * cos(φᵢ)
#         )
#         elseif params.r + dm > yᵢ >= params.r + dm - d(x) # in the interior
#             return internal_boundary_force((xᵢ, yᵢ, zᵢ), params)
#         elseif params.r + dm - d(x) > yᵢ # outermost edge
#             # get the polar coordinates rᵢ and qᵢ
#             rᵢ = sqrt((xᵢ - xc)^2 + (yᵢ - params.r - dm + d(x))^2 + (zᵢ - zc)^2)
#             if yᵢ - params.r - dm + d(x) == 0 && zᵢ >= zc
#                 qᵢ = π/2
#             elseif yᵢ - params.r - dm + d(x) == 0 && zᵢ < zc
#                 qᵢ = -π/2
#             else
#                 qᵢ = π + atan((zᵢ - zc) / (yᵢ - params.r - dm + d(x)))
#             end
#             # return boundary force
#             return SVector{3, Float64}(
#             - params.μb * exp(-(params.r - rᵢ) / params.rb) * sin(qᵢ) * sin(φᵢ),
#             - params.μb * exp(-(params.r - rᵢ) / params.rb) * cos(qᵢ),
#             - params.μb * exp(-(params.r - rᵢ) / params.rb) * sin(qᵢ) * cos(φᵢ)
#         )
#         end
#     catch err
#         if isa(err, ConvergenceError)
#             println("Convergence failure for xᵢ, yᵢ, zᵢ = $xᵢ, $yᵢ, $zᵢ")
#             println(err.err_msg)
#         elseif isa(err, DomainError)
#             println("Failure for (xᵢ - params.Xc) / params.Zc = ($xᵢ - $(params.Xc)) / $(params.Zc)")
#             println(err.err_msg)
#         else
#             rethrow()
#         end
#     end
# end

# """Get the frictionless boundary force for a flattened tube curved in the x-z plane by a function φ."""
# function right_psm_boundary_force((xᵢ, yᵢ, zᵢ), φ::Function, invφ::Function, d::Function, params::TissueParams)
#     dm = d(params.Xc) # assumes that dm == max{d(x) | x ∈ [0, Xc]}
#     # get the angle of rotation φ and centre of tissue in z axis, Zc
#     try
#         φᵢ = atan((xᵢ - params.Xc) / zᵢ)
#         x = invφ(φ)
#         xc = tissue.Xc - sin(φᵢ) * tissue.Zc
#         zc = cos(φᵢ) * params.Zc
#         if yᵢ <= 2 * params.R + params.r + dm # then at the innermost edge
#             # get the polar coordinates rᵢ and qᵢ
#             rᵢ = sqrt((xᵢ - xc)^2 + (yᵢ - 2 * params.R - params.r - dm)^2 + 2 * (zᵢ - zc)^2)
#             if yᵢ - 2 * params.R - params.r - dm == 0 && zᵢ >= zc
#                 qᵢ = π/2
#             elseif yᵢ - 2 * params.R - params.r - dm == 0 && zᵢ < zc
#                 qᵢ = -π/2
#             else
#                 qᵢ = atan((zᵢ - zc) / (yᵢ - 2 * params.R - params.r - dm))
#             end
#             # return boundary force
#             return SVector{3, Float64}(
#             - params.μb * exp(-(params.r - rᵢ) / params.rb) * sin(qᵢ) * sin(φᵢ),
#             - params.μb * exp(-(params.r - rᵢ) / params.rb) * cos(qᵢ),
#             - params.μb * exp(-(params.r - rᵢ) / params.rb) * sin(qᵢ) * cos(φᵢ)
#         )
#         elseif 2 * params.R + params.r + dm + d(x) >= yᵢ > 2 * params.R + params.r + dm # in the interior
#             return internal_boundary_force((xᵢ, yᵢ, zᵢ), params)
#         elseif yᵢ > 2 * params.R + params.r + dm + d(x) # outermost edge
#             # get the polar coordinates rᵢ and qᵢ
#             rᵢ = sqrt((xᵢ - xc)^2 + (yᵢ - 2 * params.R - params.r - dm - d(x))^2 + (zᵢ - zc)^2)
#             if yᵢ - 2 * params.R - params.r - dm - d(x) == 0 && zᵢ >= zc
#                 qᵢ = π/2
#             elseif yᵢ - 2 * params.R - params.r - dm - d(x) == 0 && zᵢ < zc
#                 qᵢ = -π/2
#             else
#                 qᵢ = π + atan((zᵢ - zc) / (yᵢ - 2 * params.R - params.r - dm - d(x)))
#             end
#             # return boundary force
#             return SVector{3, Float64}(
#             - params.μb * exp(-(params.r - rᵢ) / params.rb) * sin(qᵢ) * sin(φᵢ),
#             - params.μb * exp(-(params.r - rᵢ) / params.rb) * cos(qᵢ),
#             - params.μb * exp(-(params.r - rᵢ) / params.rb) * sin(qᵢ) * cos(φᵢ)
#         )
#         end
#     catch err
#         if isa(err, ConvergenceError)
#             println("Convergence failure for xᵢ, yᵢ, zᵢ = $xᵢ, $yᵢ, $zᵢ")
#             println(err.err_msg)
#         elseif isa(err, DomainError)
#             println("Failure for (xᵢ - params.Xc) / params.Zc = ($xᵢ - $(params.Xc)) / $(params.Zc)")
#             println(err.err_msg)
#         else
#             rethrow()
#         end
#     end
# end

"""Get the frictionless boundary force for a flattened tube curved in the x-z plane by a function φ."""
function left_psm_boundary_force((xᵢ, yᵢ, zᵢ), φ::Function, invφ::Function, d::Function, params::TissueParams)
    dm = d(params.Xc) # assumes that dm == max{d(x) | x ∈ [0, Xc]}
    # get the angle of rotation φ and centre of tissue in z axis, Zc
    try # Errors from the fact that boundary force is written incorrectly.
        φᵢ = atan((params.Xc - xᵢ) / zᵢ)
        x = - invφ(φᵢ)
        xc = tissue.Xc - sin(φᵢ) * tissue.Zc
        zc =  cos(φᵢ) * params.Zc
        if yᵢ >= params.r + dm # then at the innermost edge
            # get the polar coordinates rᵢ and qᵢ
            rᵢ = sqrt((xᵢ - xc)^2 + (yᵢ - params.r - dm)^2 + (zᵢ - zc)^2)
            if yᵢ - params.r - dm == 0 && zᵢ >= zc
                qᵢ = π/2
            elseif yᵢ - params.r - dm == 0 && zᵢ < zc
                qᵢ = -π/2
            else
                qᵢ = atan((zᵢ - zc) / (yᵢ - params.r - dm))
            end
            # return boundary force
            return SVector{3, Float64}(
            params.μb * exp(-(params.r - rᵢ) / params.rb) * (xc - xᵢ) / rᵢ,
            params.μb * exp(-(params.r - rᵢ) / params.rb) * (params.r + dm - yᵢ) / rᵢ,
            params.μb * exp(-(params.r - rᵢ) / params.rb) * (zc - zᵢ) / rᵢ
        )
        elseif params.r + dm > yᵢ >= params.r + dm - d(x) # in the interior
            return internal_boundary_force((xᵢ, yᵢ, zᵢ), params)
        elseif params.r + dm - d(x) > yᵢ # outermost edge
            # get the polar coordinates rᵢ and qᵢ
            rᵢ = sqrt((xᵢ - xc)^2 + (yᵢ - params.r - dm + d(x))^2 + (zᵢ - zc)^2)
            if yᵢ - params.r - dm + d(x) == 0 && zᵢ >= zc
                qᵢ = π/2
            elseif yᵢ - params.r - dm + d(x) == 0 && zᵢ < zc
                qᵢ = -π/2
            else
                qᵢ = π + atan((zᵢ - zc) / (yᵢ - params.r - dm + d(x)))
            end
            # return boundary force
            return SVector{3, Float64}(
            params.μb * exp(-(params.r - rᵢ) / params.rb) * (xc - xᵢ) / rᵢ,
            params.μb * exp(-(params.r - rᵢ) / params.rb) * (params.r + dm - d(x) - yᵢ) / rᵢ,
            params.μb * exp(-(params.r - rᵢ) / params.rb) * (zc - zᵢ) / rᵢ
        )
        end
    catch err
        if isa(err, ConvergenceError)
            println("Convergence failure for xᵢ, yᵢ, zᵢ = $xᵢ, $yᵢ, $zᵢ")
            println(err.err_msg)
        elseif isa(err, DomainError)
            println("Failure for (xᵢ - params.Xc) / params.Zc = ($xᵢ - $(params.Xc)) / $(params.Zc)")
            println(err.err_msg)
        else
            rethrow()
        end
    end
end

function is_external_boundary_force_working(yc, lhs::Bool, φ::Function, invφ::Function, d::Function, params::TissueParams)
    x = 150.
    xc = params.Xc - sin(π/4 + π/2) * params.Zc
    zc = - cos(π/4 + π/2) * params.Zc
    if lhs
        qᵢ = π/4 + π/2
    else
        qᵢ = - π/4 - π/2
    end
    (xᵢ, yᵢ, zᵢ) = (xc, yc, zc) .+ params.r .* (
        cos(π/4 + π/2),
        cos(qᵢ),
        sin(π/4 + π/2) * sin(qᵢ)
    )
    if ~isapprox(
        left_psm_boundary_force((xᵢ, yᵢ, zᵢ), φ, invφ, d, params),
        - params.μb .* (cos(π/4 + π/2), cos(qᵢ), sin(π/4 + π/2) * sin(qᵢ))
    )
    return false
    end
    true
end

function is_left_psm_boundary_force_working(φ::Function, invφ::Function, d::Function, params::TissueParams)
    xc = params.Xc - sin(π/4 + π/2) * params.Zc
    zc = cos(π/4 + π/2) * params.Zc
    
    if ~isapprox(internal_boundary_force((xc, params.r + d(tissue.Xc), zc), params), left_psm_boundary_force((xc, params.r + d(tissue.Xc), zc), φ, invφ, d, params))
        return false
    elseif ~isapprox(internal_boundary_force((xc, params.r + d(tissue.Xc) - d(x), zc), params), left_psm_boundary_force((xc, params.r + d(tissue.Xc) - d(x), zc), φ, invφ, d, params))
        return false
    else
    end
    true
end

function right_psm_boundary_force((xᵢ, yᵢ, zᵢ), φ::Function, invφ::Function, d::Function, params::TissueParams)
    dm = d(params.Xc) # assumes that dm == max{d(x) | x ∈ [0, Xc]}
    # get the angle of rotation φ and centre of tissue in z axis, Zc
    try
        φᵢ = atan((params.Xc - xᵢ) / zᵢ)
        x = -invφ(φᵢ)
        xc = tissue.Xc - sin(φᵢ) * tissue.Zc
        zc = cos(φᵢ) * params.Zc
        if yᵢ <= 2 * params.R + params.r + dm # then at the innermost edge
            # get the polar coordinates rᵢ and qᵢ
            rᵢ = sqrt((xᵢ - xc)^2 + (yᵢ - 2 * params.R - params.r - dm)^2 + 2 * (zᵢ - zc)^2)
            if yᵢ - 2 * params.R - params.r - dm == 0 && zᵢ >= zc
                qᵢ = π/2
            elseif yᵢ - 2 * params.R - params.r - dm == 0 && zᵢ < zc
                qᵢ = -π/2
            else
                qᵢ = atan((zᵢ - zc) / (yᵢ - 2 * params.R - params.r - dm))
            end
            # return boundary force
            return SVector{3, Float64}(
            params.μb * exp(-(params.r - rᵢ) / params.rb) * (xc - xᵢ) / rᵢ,
            params.μb * exp(-(params.r - rᵢ) / params.rb) * (2 * params.R + params.r + dm - yᵢ) / rᵢ,
            params.μb * exp(-(params.r - rᵢ) / params.rb) * (zc - zᵢ) / rᵢ
        )
        elseif 2 * params.R + params.r + dm + d(x) >= yᵢ > 2 * params.R + params.r + dm # in the interior
            return internal_boundary_force((xᵢ, yᵢ, zᵢ), params)
        elseif yᵢ > 2 * params.R + params.r + dm + d(x) # outermost edge
            # get the polar coordinates rᵢ and qᵢ
            rᵢ = sqrt((xᵢ - xc)^2 + (yᵢ - 2 * params.R - params.r - dm - d(x))^2 + (zᵢ - zc)^2)
            if yᵢ - 2 * params.R - params.r - dm - d(x) == 0 && zᵢ >= zc
                qᵢ = π/2
            elseif yᵢ - 2 * params.R - params.r - dm - d(x) == 0 && zᵢ < zc
                qᵢ = -π/2
            else
                qᵢ = π + atan((zᵢ - zc) / (yᵢ - 2 * params.R - params.r - dm - d(x)))
            end
            # return boundary force
            return SVector{3, Float64}(
            params.μb * exp(-(params.r - rᵢ) / params.rb) * (xc - xᵢ) / rᵢ,
            params.μb * exp(-(params.r - rᵢ) / params.rb) * (2 * params.R + params.r + dm + d(x) - yᵢ) / rᵢ,
            params.μb * exp(-(params.r - rᵢ) / params.rb) * (zc - zᵢ) / rᵢ
        )
        end
    catch err
        if isa(err, ConvergenceError)
            println("Convergence failure for xᵢ, yᵢ, zᵢ = $xᵢ, $yᵢ, $zᵢ")
            println(err.err_msg)
        elseif isa(err, DomainError)
            println("Failure for (xᵢ - params.Xc) / params.Zc = ($xᵢ - $(params.Xc)) / $(params.Zc)")
            println(err.err_msg)
        else
            rethrow()
        end
    end
end

"""Get the frictionless boundary force for a flattened half torus."""
function flat_torus_boundary_force((xᵢ, yᵢ, zᵢ), d::Function, params::TissueParams)
    dm = d(params.Xc)
    # get the toroidal coordinates, Rᵢ, pᵢ
    if yᵢ - params.Yc > 0 && xᵢ - params.Xc == 0
        pᵢ = π/2
    elseif yᵢ - params.Yc < 0 && xᵢ - params.Xc == 0
        pᵢ = -π/2
    else
        pᵢ = atan((yᵢ - params.Yc) / (xᵢ - params.Xc))
    end
    Rᵢ = sqrt((xᵢ - params.Xc)^2 + (yᵢ - params.Yc)^2)
    # Now use this to infer where in the tissue the cell is, and calculate the corresponding force.
    if Rᵢ <= params.R # interior periphery
        rᵢ = sqrt((xᵢ - params.R * cos(pᵢ))^2 + (yᵢ - params.R * sin(pᵢ))^2 + (zᵢ - params.Zc)^2)
        if rᵢ == 0
            qᵢ = 0.
        else
            qᵢ = asin((zᵢ - params.Zc) / rᵢ)
        end
        δx = abs((params.r - rᵢ) * cos(pᵢ) * cos(qᵢ))
        δy = abs((params.r - rᵢ) * sin(pᵢ) * cos(qᵢ))
        δz = abs((params.r - rᵢ) * sin(qᵢ))
        return SVector{3, Float64}(
            - params.μb * exp(- δx / params.rb) * cos(pᵢ) * cos(qᵢ),
            - params.μb * exp(- δy / params.rb) * sin(pᵢ) * cos(qᵢ),
            - params.μb * exp(- δz / params.rb) * sin(qᵢ)
        )
    elseif params.R + dm >= Rᵢ > params.R # interior
        if zᵢ >= params.Zc
            δz = abs(params.Zc + params.r - zᵢ)
            return SVector{3, Float64}(
                0,
                0,
                - params.μb * exp(- δz / params.rb)
            )
        else
            δz = abs(zᵢ - params.Zc + params.r)
            return SVector{3, Float64}(
                0,
                0,
                params.μb * exp(- δz / params.rb)
            )
        end
    elseif Rᵢ > params.R + dm
        rᵢ = sqrt((xᵢ - (params.R + dm) * cos(pᵢ))^2 + (yᵢ - (params.R + dm) * sin(pᵢ))^2 + (zᵢ - params.Zc)^2)
        if rᵢ == 0
            qᵢ = 0.
        else
            qᵢ = asin((zᵢ - params.Zc) / rᵢ)
        end
        δx = abs((params.r - rᵢ) * cos(pᵢ) * cos(qᵢ))
        δy = abs((params.r - rᵢ) * sin(pᵢ) * cos(qᵢ))
        δz = abs((params.r - rᵢ) * sin(qᵢ))
        return SVector{3, Float64}(
            - params.μb * exp(- δx / params.rb) * cos(pᵢ) * cos(qᵢ),
            - params.μb * exp(- δy / params.rb) * sin(pᵢ) * cos(qᵢ),
            - params.μb * exp(- δz / params.rb) * sin(qᵢ)
        )
    end
end

"""Get the boundary force in a spatulate, curved PSM."""
function boundary_force_spatula((xᵢ, yᵢ, zᵢ), d::Function, φ::Function, invφ::Function, params::TissueParams)
    if xᵢ >= params.Xc
        try
            return flat_torus_boundary_force((xᵢ, yᵢ, zᵢ), d, params)
        catch err
            println("Error in toroid domain.")
            rethrow()
        end
    elseif yᵢ < params.Yc
        try
        return left_psm_boundary_force((xᵢ, yᵢ, zᵢ), φ, invφ, d, params)
        catch err
            println("Error in left hand psm.")
            rethrow()
        end
    elseif yᵢ > params.Yc
        try
        return right_psm_boundary_force((xᵢ, yᵢ, zᵢ), φ, invφ, d, params)
        catch err
            println("Error in right hand psm.")
            rethrow()
        end
    end
end