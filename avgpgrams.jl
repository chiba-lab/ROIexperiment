#* 1 filter and group data 
#* 2 run DSP
#* 3 averaging 
#* 4 plots
#* 5 Stats
##==============================================================================
include("./ROIontology.jl")

##==============================================================================
#Load data

ROIexp = load("./Data/ROIexp.jld2", "ROIexp");
lfp = ROIexp[:lfp_data]
events = ROIexp[:behavioral_events]
trials = ROIexp[:trials]
sessions = ROIexp[:sessions]

include("./Links.jl")

##==============================================================================
# #univariate analysis

using FourierAnalysis

const WD = filter(x -> window_data_event(window_data(x)) ∈ dom(event_trial), collect(dom(window_data)))
const WC = filter(x -> window_data_event(x) ∈ dom(event_trial), collect(codom(window_data)))

# WT = event_trial∘window_data_event
# using JLD2 
# @save "TempData/event_lfp_windows_1_1.jld2"  WC
##==============================================================================
#split events into pre-peri-post with 512 sample pre and post windows
#events shorter than 512 are expanded
wd_pre = GMap(x -> x[498:1010], WC)
wd_peri = GMap(x -> length(x) <= 2532 ? x[round(Int, length(x) / 2)-256:round(Int, length(x) / 2)+256] : x[1010:end-1010], WC)
wd_post = GMap(x -> x[end-1010:end-498], WC)

#Freq Dom
# pre_spec = GMap(x -> spectra(x, 256, 512), wd_pre.(WC))
peri_spec = GMap(x -> spectra(x, 256, 512), wd_peri.(WC))
# post_spec = GMap(x -> spectra(x, 256, 512), wd_post.(WC)) 
##==============================================================================
using CairoMakie
include("./PlotFuncs.jl")
using SplitApplyCombine

using Statistics

function fuckyouplot(e, tc, r)
    # dataregion = MOB
    event_name = e
    # agent_types = ["Rat", "Robot", "Object"]
    trial_condition = tc
    region = r

    idxs = findall(x -> window_data_event(x).behavior.name == event_name, WC) ∩ findall(x -> event_trial(window_data_event(x)).condition.name == trial_condition, WC) ∩ findall(x -> window_data_region(x) == region, WC)

    g = group(x -> agenttype(event_trial(window_data_event(x)).condition), WC[idxs])
    dct = Dict([(Object, "Object"), (Rat, "Rat"), (Robot, "Robot")])


    f = Figure()
    ax = Axis(f[1, 1], yscale = log10,
        yminorticksvisible = true, yminorgridvisible = true,
        yminorticks = IntervalsBetween(8))

    for i in collect(1:1:length(g))
        k = collect(keys(g))[i]
        spcs = g[k] |> x -> wd_peri.(x) |> x -> peri_spec.(x)
        fl = spcs[1].flabels
        spcs = map(x -> x.y, spcs) |> x -> hcat(x...)
        m = mean(spcs, dims = 2)
        s = std(spcs, dims = 2) / sqrt(size(spcs, 2))
        # errorbars!(ax, fl, vec(m), vec(s); linewidth = 2, label = k)
        lines!(fl, vec(m), label = dct[k]; linewidth = 2)
    end
    f[1, 2] = Legend(f, ax, "Agent")
    ax.xlabel = "Freq (Hz)"
    ax.ylabel = "Amp"
    ax.title = "Condition:$trial_condition, Event:$event_name, Region:$(region.name)"
    save("./Mean specs/$(trial_condition)_$(event_name)_$(region.name).png", f)
end

for e in ["Groom", "Rear", "Baseline", "Immobility"]
    for tc in ["Free Roam"]
        for r in [MOB, CA2, AMG]
            fuckyouplot(e, tc, r)
        end
    end
end

##==============================================================================
# f = Figure();
# ax = Axis(f[1, 1], yscale = log10,
#     yminorticksvisible = true, yminorgridvisible = true,
#     yminorticks = IntervalsBetween(8))
# errorbands!(ax, fl, vec(m), vec(s); linewidth = 2, label = "")
# ax.xlabel = "Freq (Hz)"
# ax.ylabel = "Amp"


# # ##==============================================================================
# # #average
# # using Statistics
# # pspec = hcat(event_subs_power.(event_subs_lfp_windows)...)
# # m = mean(mat, dims = 2)
# # s = std(mat, dims = 2) / sqrt(size(mat, 2))


# f = Figure();
# ax = Axis(f[1, 1], yscale = log10,
#     yminorticksvisible = true, yminorgridvisible = true,
#     yminorticks = IntervalsBetween(8))
# errorbands!(ax, 1:1:256, vec(m), vec(s); linewidth = 2, label = "")
# ax.ylabel = "Freq (Hz)"
# ax.xlabel = "Amp"
# f