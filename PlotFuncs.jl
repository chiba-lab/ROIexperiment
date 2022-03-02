function eventwindow!(ax, times, data, onset, offset; kwargs...)
    plot!(ax,times.-onset, data, kwargs...)
    vlines!(ax, 0; color=:blue, linewidth=2, linestyle=:dash)
    vlines!(ax, offset-onset; color=:red, linewidth=2, linestyle=:dash)
    ax.xticks = ([times[1]-onset, 0, offset-onset, times[end]-onset], ["$(times[1]-onset)", "Onset", "Offset", "$(times[end]-onset)"])
end

function errorbands!(ax, x, y, error; kwargs...)
    l=lines!(ax, x, y; kwargs...)
    b=band!(ax, x, y-error, y+error; color=(l.color, 0.4))
    return l, b
end



# function pgram(ax, x, y, error; kwargs...)

#     ax = Axis(ax, yscale = log10,
#         yminorticksvisible = true, yminorgridvisible = true,
#         yminorticks = IntervalsBetween(8))
#     errorbands!(ax, x, vec(m), vec(s); linewidth = 2)
# end