# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# This file contains the objects containing cell attributes and the parameters describing the clock and tissue.

using StaticArrays

"""
FCell

Struct that records the frequency dθ/dt for a cell.

x::SVector{3, Float64} => position of cell in 3 dimensions
n::SVector{3, Float64} => direction of intrinsic cell movement
dθdt::Float64 => frequency of phase oscillations
θ::Float64 => clock phase
τ::Float64 => cell cycle phase (∈ [0, TG+TM))
trackid::Float64 => identity of the cell

"""
struct FCell # called 'frequency tracking cells'
    x::SVector{3, Float64}
    n::SVector{3, Float64}
    dθdt::Float64
    θ::Float64
    τ::Float64
    trackid::Float64
end

Base.copy(x::FCell) = FCell(x.x, x.n, x.dθdt, x.θ, x.τ, x.trackid)
Base.copy(X::Vector{FCell}) = copy.(X)

"""
TissueParams

Encodes parameters for the tissue.

dc: diameter of Cells (in μm)
μ: coefficient for interCellular force
γₐ: parameter for Cell advection speed in anterior
γₚ: parameter for Cell advection speed in posterior
xq: position defining anterior and posterior PSM
L: length of tissue
xa: x coordinate of most anterior Cell
Xᵥ: length scale of the Cell motility gradient
vs: max velocity in posterior of the tailbud
h: coefficient controlling shape of Cell motility gradient
Dφ: coefficient controlling noise in Cell directionality
μb: coefficient of the boundary force
rb: lengthscale of the boundary force
Xc: x coordinate demarcating the psm columns and torus
Yc: y coordinate demarcating the psm columns and torus
Zc: z coordinate demarcating the psm columns and torus
r: radius of psm columns and torus
R: radius from torus centre to midline
ρ₀: desired constant density of the tissue
N_cells: number of Cells in the tissue
total_cells: total Cells that have passed through the tissue
TG: duration of growth phase of cell cycle (mins)
TM: duration of M phase of cell cycle (mins)
xd: anterior limit of cell addition
kc: coupling strength of cell cycle to segmentation clock phase
L₀: initial length of PSM at t = 0
uₐ: rate of PSM length shrinkage
d₀: density at t=0
md: rate of density decrease
r₀: initial PSM radius
mᵣ: rate of decrease in PSM radius
t_shrink: time at which tissue begins to shrink
Lx: length of lattice domain in x
Ly: length of lattice domain in y
Lz: length of lattice domain in z
dm: controls distance of cell-cell repulsion, where one wishes this to vary
dc₀: initial cell diameter, if one wishes this to vary linearly over time
mdc: magnitude of linear cell diameter change, if one wishes this to change
lattice_side: controls the length of lattice cube domains
t_stop: time at which to stop shrinking cell diameter
"""
mutable struct TissueParams
    dc::Float64 # diameter of Cells (in μm)
    μ::Float64 # coefficient for interCellular force
    γₐ::Float64 # parameter for Cell advection speed in anterior
    γₚ::Float64 # parameter for Cell advection speed in posterior
    xq::Float64 # position defining anterior and posterior PSM
    L::Float64 # length of tissue
    xa::Float64 # x coordinate of most anterior Cell
    Xᵥ::Float64 # length scale of the Cell motility gradient
    vs::Float64 # max velocity in posterior of the tailbud
    h::Float64 # coefficient controlling shape of Cell motility gradient
    Dφ::Float64 # coefficient controlling noise in Cell directionality
    μb::Float64 # coefficient of the boundary force
    rb::Float64 # lengthscale of the boundary force
    Xc::Float64 # x coordinate demarcating the psm columns and torus
    Yc::Float64 # y coordinate demarcating the psm columns and torus
    Zc::Float64 # z coordinate demarcating the psm columns and torus
    r::Float64 # radius of psm columns and torus
    R::Float64 # radius from torus centre to midline
    ρ₀::Float64 # desired constant density of the tissue
    N_cells::Int64 # number of Cells in the tissue
    total_cells::Int64 # total Cells that have passed through the tissue
    TG::Float64 # duration of growth phase of cell cycle (mins)
    TM::Float64 # duration of M phase of cell cycle (mins)
    xd::Float64 # anterior limit of cell addition
    kc::Float64 # coupling strength of cell cycle to segmentation clock phase
    dm::Float64 # controls distance of repulsion for cells, where one wishes this to vary
    L₀::Float64 # initial length of PSM at t = 0
    uₐ::Float64 # rate of PSM length shrinkage
    d₀::Float64 # density at t=0
    md::Float64 # rate of density decrease
    r₀::Float64 # initial PSM radius
    mᵣ::Float64 # rate of decrease in PSM radius
    t_shrink::Float64 # time at which tissue begins to shrink
    Lx::Float64 # length of lattice domain in x
    Ly::Float64 # length of lattice domain in y
    Lz::Float64 # length of lattice domain in z
    x_div::Float64 # anterior limit of cell division
    dc₀::Float64 # initial cell diameter, if one wishes this to vary linearly over time
    mdc::Float64 # magnitude of linear cell diameter change, if one wishes this to change
    lattice_side::Float64 # controls the length of lattice cube domains
    t_stop::Float64 # time at which to stop shrinking cell diameter
end

Base.copy(tissue::TissueParams) = TissueParams(
    tissue.dc,
    tissue.μ,
    tissue.γₐ,
    tissue.γₚ,
    tissue.xq,
    tissue.L,
    tissue.xa,
    tissue.Xᵥ,
    tissue.vs,
    tissue.h,
    tissue.Dφ,
    tissue.μb,
    tissue.rb,
    tissue.Xc,
    tissue.Yc,
    tissue.Zc,
    tissue.r,
    tissue.R,
    tissue.ρ₀,
    tissue.N_cells,
    tissue.total_cells,
    tissue.TG, 
    tissue.TM,
    tissue.xd,
    tissue.kc,
    tissue.dm,
    tissue.L₀,
    tissue.uₐ,
    tissue.d₀,
    tissue.md,
    tissue.r₀,
    tissue.mᵣ,
    tissue.t_shrink,
    tissue.Lx,
    tissue.Ly,
    tissue.Lz,
    tissue.x_div,
    tissue.dc₀,
    tissue.mdc,
    tissue.lattice_side,
    tissue.t_stop
)

"""
PhaseParams

Encodes parameters for the kuramoto model.

k0: coupling constant for the kuramoto model
Dθ: phase noise intensity
ωₒ: intrinsic oscillator frequency in posterior of tissue
σ: fold change in frequency between anterior and posterior of tissue
k: shape of frequency gradient (dimensionless parameter)
nτ: delay of the kuramoto model, in discrete time steps of the numerical scheme
"""
mutable struct PhaseParams
    k0::Float64 # coupling constant for the kuramoto model
    Dθ::Float64 # phase noise intensity
    ωₒ::Float64 # intrinsic oscillator frequency in posterior of tissue
    σ::Float64 # fold change in frequency between anterior and posterior of tissue
    k::Float64 # shape of frequency gradient (dimensionless parameter)
    nτ::Int64 # delay of the kuramoto model, in discrete time steps of the numerical scheme
end

Base.copy(x::PhaseParams) = PhaseParams(
    x.k0,
    x.Dθ,
    x.ωₒ,
    x.σ,
    x.k,
    x.nτ,
)
