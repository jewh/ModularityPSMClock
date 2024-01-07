# James Hammond (james.hammond@biology.ox.ac.uk) 6th Jan 2023

# plot the functions that feed into compaction-extension simulations

using CairoMakie

include("dynamic_tissue_geometry.jl")

begin

    colors = [Makie.wong_colors()[1], Makie.wong_colors()[2]]

    fig = Figure()

    ax1 = Axis(fig[1,1],
            xlabel=L"t \; \text{(mins)}",
            ylabel=L"L(t)  \; \text{(\mu\!\>m)}", 
            xticks = (
                [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
            yticks = (
                [x for x = 0:100:350], string.([x for x = 0:100:350])
            ))

    ylims!(ax1, 0, 350)

    ax2 = Axis(fig[1,2],xlabel=L"t\; \text{(mins)}",
            ylabel=L"r(t)\; \text{(\mu\!\>m)}", 
            xticks = (
                [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
            yticks = (
                [x for x = 0:15:45], string.([x for x = 0:15:45])
            ))

    ylims!(ax2, 0, 30)

    ax3 = Axis(fig[2,1],
            xlabel=L"t\; \text{(mins)}",
            ylabel=L"\rho(t)\; \text{(\mu\!\>m^{-3})}", 
            xticks = (
                [x for x = 0:200:1000], string.([x for x = 0:200:1000])),
            yticks = (
                [x for x = 0:0.002:0.004], string.([x for x = 0:0.002:0.004])
                ))

    ylims!(ax3, 0, 0.005)

    ax4 = Axis(fig[2,2],
            xlabel = L"t \; \text{(mins)}", 
            ylabel = L"d_{c} \; \text{(μm)}",
            xticks = (
            [x for x = 0:200:1000], string.([x for x = 0:200:1000])
            ))

    ylims!(ax4, 0, 11)

    ts = 0:772

    t_shrink = 253.6
    t_stop = 536.4 # time at which cell size stops shrinking

    # annotate time points
    vlines!(ax1, [t_shrink]; color=:black, linestyle=:dash)
    vlines!(ax2, [t_shrink]; color=:black, linestyle=:dash)
    vlines!(ax3, [t_shrink]; color=:black, linestyle=:dash)
    vlines!(ax4, [t_shrink]; color=:black, linestyle=:dash)
    vlines!(ax4, [t_stop]; color=:black, linestyle=:dash)

    # now plot for non-shrinking case
    Lt = [Length(t, 0.0, 325, t_shrink) for t in ts]
    rt = [radius(t, 0.0, 27.6, t_shrink) for t in ts]
    dt = [density(t, 0.0, 0.002123, t_shrink) for t in ts]
    stepwise_shrinking = [cell_diameter(t, 0.0, 9.2, t_shrink, t_stop) for t in ts]

    lines!(ax1, ts, Lt, color=colors[1])
    lines!(ax2, ts, rt, color=colors[1])
    lines!(ax3, ts, dt, color=colors[1])
    lines!(ax4, ts, stepwise_shrinking, color=colors[1])

    # plot for shrinking case
    Lt = [Length(t, 14.1, 325, t_shrink) for t in ts]
    rt = [radius(t, 1.1625, 27.6, t_shrink) for t in ts]
    dt = [density(t, 0.000121, 0.002123, t_shrink) for t in ts]
    stepwise_shrinking = [cell_diameter(t, 0.2, 9.2, t_shrink, maximum(ts)) for t in ts]

    lines!(ax1, ts, Lt, color=colors[2])
    lines!(ax2, ts, rt, color=colors[2])
    lines!(ax3, ts, dt, color=colors[2])
    lines!(ax4, ts, stepwise_shrinking, color=colors[2])

    # annotate
    text!(ax1, L"t = t_{shrink}", position=(270, 32))
    text!(ax2, L"t = t_{shrink}", position=(270, 2.72))
    text!(ax3, L"t = t_{shrink}", position=(270, 0.00045))
    text!(ax4, L"t = t_{shrink}", position=(270, 1))
    text!(ax4, L"t = t_{stop}", position=(550, 1))

    save("figures/20231204_shrinking_model_functions.svg", fig)

    fig
end
