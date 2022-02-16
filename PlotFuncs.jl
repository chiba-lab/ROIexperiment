function eventwindow!(ax, times, data, onset, offset; kwargs...)
    plot!(ax,times.-onset, data, kwargs...)
    vlines!(ax, 0; color=:blue, linewidth=2, linestyle=:dash)
    vlines!(ax, offset-onset; color=:red, linewidth=2, linestyle=:dash)
    ax.xticks = ([times[1]-onset, 0, offset-onset, times[end]-onset], ["$(times[1]-onset)", "Onset", "Offset", "$(times[end]-onset)"])
end

function errorbands!(ax, x, y, error; kwargs...)
    l=lines!(ax, x, y; kwargs...)
    band!(ax, x, y-error, y+error; color=(l.color, 0.4))
end

