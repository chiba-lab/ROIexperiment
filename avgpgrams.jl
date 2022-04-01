#* 1 filter and group data 
#* 2 run DSP
#* 3 averaging 
#* 4 plots
#* 5 Stats
##==============================================================================
include("./ROIontology.jl")


##==============================================================================
#Load data

ROIexp = load("./Data/ROIexp_FR.jld2", "ROIexp");
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
# wd_pre = GMap(x -> x[498:1010], WC)
wd_peri = GMap(x -> length(x) <= 2532 ? x[round(Int, length(x) / 2)-256:round(Int, length(x) / 2)+256] : x[1010:end-1010], WC)
# wd_post = GMap(x -> x[end-1010:end-498], WC)

#Freq Dom
# pre_spec = GMap(x -> spectra(x, 256, 512), wd_pre.(WC))
peri_spec = GMap(x -> spectra(x, 256, 512), wd_peri.(WC))
# post_spec = GMap(x -> spectra(x, 256, 512), wd_post.(WC)) 
##==============================================================================
using CairoMakie
include("./PlotFuncs.jl")
import SplitApplyCombine as sac

using DataFrames
import MixedModels: term
using MixedModels
using Statistics
using ProgressBars



function fuckyouplot(ax, e, tc, r)
    # dataregion = MOB
    event_name = e
    # agent_types = ["Rat", "Robot", "Object"]
    trial_condition = tc
    region = r

    idxs = findall(x -> window_data_event(x).behavior.name == event_name, WC)  ∩ findall(x -> window_data_region(x) == region, WC)

    g = sac.group(x -> agenttype(event_trial(window_data_event(x)).condition), WC[idxs])
    dct = Dict([(Object, "Object"), (Rat, "Rat"), (Robot, "Robot")])
    dcc = Dict([(Object, :red), (Rat, :green), (Robot, :blue), ("Object", :red), ("Rat", :green), ("Robot", :blue)])




    # f = Figure()
    # ax = Axis(f[1, 1], xscale = log10, yscale = log10,
    #     yminorticksvisible = true, yminorgridvisible = true,
    #     yminorticks = IntervalsBetween(8))
    av = ["Object", "Rat", "Robot"]
    # for k in av
    # k = collect(keys(g))[i]
    agent_type = map(x -> agenttype(event_trial(window_data_event(x)).condition), WC[idxs])
    ratname = map(x -> window_lfp(inv(window_data)(x)).rat.name, WC[idxs])
    spcs = WC[idxs] |> x -> wd_peri.(x) |> x -> peri_spec.(x)
    fr = FourierAnalysis.fres(256, 512) .* (1010.1 / 256)
    fl = collect(fr:fr:fr*256)

    fls = fl |> x -> round.(Int, x .* 100)
    # uid = findfirst(x -> x > 150, fl)
    # fl = fl[1:uid]

    spcs = map(x -> x.y ./ sum(x.y), spcs) |> x -> hcat(x...)

    fcols = ["f_$i" for i = fls]
    spcstbl = DataFrame()
    spcstbl.rat = ratname
    spcstbl.agent_type = agent_type
    for i in 1:length(fcols)
        spcstbl[!, Symbol(fcols[i])] = spcs[i, :]
    end

    function mefx(x)
        fm = (term(x) ~ term(1) + term(:agent_type) + (term(1) | term(:rat)))
        fm1 = fit(MixedModel, fm, spcstbl)
        # return coef(fm1), stderror(fm1)
        return fm1
    end

    # if length(unique(ratname)) == 1
    #     m = mean(spcs, dims = 2)
    #     s = std(spcs, dims = 2) / sqrt(size(spcs, 2))
    # else
    #     ms = [mefx(fc) for fc = Symbol.(fcols)]
    # end


    ms = [mefx(fc) for fc = Symbol.(fcols)]



    # errorbars!(ax, fl, vec(m), vec(s); linewidth = 2, label = k)
    for nm in av
    
        # m = map(x -> x[1][k], ms) |> vcat
        # s = map(x -> x[2][k], ms) |> vcat
        # nms=split()
        k = findfirst(x -> x == nm, map(x -> x[2], split.(fixefnames(ms[1])[2:end], ": ")))
    
        m = map(x -> x.beta[k], ms)
        s = map(x -> x.stderror[k], ms)
       
    
        lines!(ax, fl, max.(m, 0.0001), label = nm, color = dcc[nm]; linewidth = 2)
        # axislegend(ax)
        band!(ax, fl, max.(m .- s, 0.0001), m .+ s, color = (dcc[nm], 0.4))
        xlims!(ax, fl[1], 50)
        ylims!(ax, 0.001, 1)
        ax.xticks = ([5, 10, 25, 45], ["5", "10", "25", "45"])
    
    
    end
    # end
    # fl = spcs[1].flabels .* 1010.1 / 256
    # xlims!(fl[1], 120)
    # f[1, 2] = Legend(f, ax, "Agent")
    # ax.xlabel = "Freq (Hz)"
    # ax.ylabel = "Amp"
    ax.title = "$event_name, $(dcn[region])"
    # save("./TempData/MeanSpecs/$(trial_condition)_$(event_name)_$(region.name).png", f)
    # f
end

dcn = Dict([(AMG, "Amygdala"), (MOB, "MOB"), (CA2, "CA")])

f = Figure();
axs = [Axis(f[i, j], yscale = log10, xscale = log10, yminorticksvisible = true, yminorgridvisible = true) for i = 1:3, j = 1:3]
cnt = 1
for e in tqdm(["Grooming", "Rearing", "Immobility"])
    for tc in ["Free Roam"]
        for r in [MOB, CA2, AMG]
            fuckyouplot(axs[cnt], e, tc, r)
            if e == "Grooming"
                axs[cnt].ylabel = "Amp (log10)"
            end
            if r == AMG
                axs[cnt].xlabel = "Freq (Hz)"
            end

            #     # axs[cnt].title = "$(dcn[r])"
            # end
            # axs[cnt].title = "$(dcn[r]), $e"
            # if e == "Grooming"
            #     axs[cnt].ylabel = "Amp"
            # end
            # if r == AMG
            #     axs[cnt].xlabel = "Freq (Hz)"
            # end
            cnt = cnt + 1

        end
    end
end

# f
# xlims!(ax, fl[1], 150)
#
save("./TempData/MeanSpecs/allspecs_corrected_50.png", f)
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
##
# using CairoMakie
# include("./PlotFuncs.jl")


# using Statistics
# using DataFrames

# idxs = findall(x -> window_data_event(x).behavior.name == "Grooming", WC) ∩ findall(x -> event_trial(window_data_event(x)).condition.name == "Free Roam", WC) ∩ findall(x -> window_data_region(x) == MOB, WC)

# dta = WC[idxs]
# ratname = map(x -> window_lfp(inv(window_data)(x)).rat.name, WC[idxs])
# agent_type = map(x -> agenttype(event_trial(window_data_event(x)).condition), WC[idxs])
# spcs = dta |> x -> wd_peri.(x) |> x -> peri_spec.(x)
# fl = spcs[1].flabels .* (1010.1 / 256) |> x -> round.(Int, x .* 100)
# spcs = map(x -> x.y, spcs) |> x -> hcat(x...)

# fcols = ["f_$i" for i = fl]
# spcstbl = DataFrame()
# spcstbl.rat = ratname
# spcstbl.agent_type = agent_type
# for i in 1:length(fcols)
#     spcstbl[!, Symbol(fcols[i])] = spcs[i, :]
# end
# spcstbl

# import MixedModels as mm

# fm = (mm.term(:f_197) ~  mm.term(:agent_type) + (mm.term(1) | mm.term(:rat)))
# fm1 = fit(MixedModel, fm, spcstbl, contrasts = Dict(:agent_type => fulldummy()))
# coef(fm1)
# stderror(fm1)



# function mefx(x)
#     fm = (mm.term(x) ~ mm.term(0) + mm.term(:agent_type) + (mm.term(1) | mm.term(:rat)))
#     fm1 = fit(MixedModel, fm, spcstbl)
#     return coef(fm1), stderror(fm1)
# end



# ms = [mefx(fc) for fc = Symbol.(fcols)]

# m = map(x -> x[1][2], ms) |> vcat
# s = map(x -> x[2][2], ms) |> vcat

# s
# m
# f = Figure()
# ax = Axis(f[1, 1], yscale = log10, xscale = log10, yminorticksvisible = true, yminorgridvisible = true,
#     yminorticks = IntervalsBetween(8))



# spcs = dta |> x -> wd_peri.(x) |> x -> peri_spec.(x)
# fl = spcs[1].flabels * 1010.1 / 128
# spcs = map(x -> x.y, spcs) |> x -> hcat(x...)
# m = mean(spcs, dims = 2)
# s = std(spcs, dims = 2) / sqrt(size(spcs, 2))
# # errorbars!(ax, fl, vec(m), vec(s); linewidth = 2, label = k)
# lines!(fl, vec(m); linewidth = 2)

# f
# # f[1, 2] = Legend(f, ax, "Agent")
# # ax.xlabel = "Freq (Hz)"
# # ax.ylabel = "Amp"
# # ax.title = "Condition:$trial_condition, Event:$event_name, Region:$(region.name)"
# # save("./Mean specs/$(trial_condition)_$(event_name)_$(region.name).png", f)