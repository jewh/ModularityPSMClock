# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# Code for calculating rate of somitogenesis over time based on data from Schroter et al 2008

using CSV, DataFrames, LsqFit, Plots

Plots.theme(
    :vibrant; 
    framestyle=:box
    )

df = DataFrame(CSV.File("data/schroteretal2008somitenumber.csv"))

begin
    # define model form

    model(somite_number, params) = 24.7 * (somite_number .- 6) + params[1] * exp.(params[2] * somite_number)

    # define guesses

    guess_params = [0.1, 0.1]

    # Perform fitting

    fit = curve_fit(model, df[!, :Somite], df[!, :Time], guess_params)


    # check results 

    fit.param

    estimate_covar(fit)

    estimate_errors(fit)

    correct = fit.param
end

# define model form
begin

    s(time, params) = time / 24.7 .+ 6 - params[1] * exp.(params[2] * (time .- 300))

    # define guesses

    guess_params = [0.1, 0.1]

    # Perform fitting

    fit = curve_fit(s, df[!, :Time], df[!, :Somite], guess_params)


    # check results 

    fit.param

    estimate_covar(fit)

    estimate_errors(fit)
end

# Plot to visually check sanity

plot(
    df[!, :Somite],
    [model(x, fit.param) for x in df[!, :Somite]],
    label="fitted model",
    legend=:outertopright
)
plot!(
    df[!, :Somite],
    [24.7* (x - 6)  for x in df[!, :Somite]],
    label="linear trendline",
    color=:black
)
scatter!(df[!, :Somite], df[!, :Time], label=false, color=:red)
ylims!((0, 800))
xlims!((5, 35))
ylabel!("Time (mins)")
xlabel!("Somite number")
savefig("figures/20221117_fitted_model_schroteretal2008.svg")

#----------------------------------
plot(
    df[!, :Time],
    s(df[!, :Time], fit.param),
    label="fitted model",
    legend=:outertopright
)

plot!(
    df[!, :Time],
    [x/24.7 + 6  for x in df[!, :Time]],
    label="constant rate",
    color=:black
)

scatter!(df[!, :Time], df[!, :Somite], label="data", color=:red)
xlims!((0, 800))
ylims!((5, 35))
ylabel!("Somite Number")
xlabel!("Time (mins)")
savefig("figures/20221117_fitted_model_schroteretal2008.svg")

# it works!

# Now plot to check PSM volume over time

# volume of PSM over time, according to Thomson et al. 2021
V(s) = 16.7e6 * exp(-0.159 * s)

s = 15:32
t = [model(s, fit.param) for s in s]
v = [V(s) for s in s]

plot(t, v, label=false)
xlabel!("t (mins)")
ylabel!("V (μm³)")

# Now compare with function of just decreasing radius of cylinder

L(s::Real) = -14.1 * s + 155

function model_V(somite_number, params)
    out = Vector{Float64}(undef, length(somite_number))
    for i in eachindex(somite_number)
        out[i] = L(somite_number[i]) * π * (params[1] - params[2] * somite_number[i])
    end
    out
end

guess_params = [25, 0.1]

fit_V = curve_fit(model_V, s, v, guess_params)

# check results 

fit_V.param

estimate_covar(fit_V)

estimate_errors(fit_V)

# visually check fit

mv = [model_V(ss, fit_V.param) for ss in s]
plot(s, v, label=false)
plot!(s, mv, label=false)

println("\007")
