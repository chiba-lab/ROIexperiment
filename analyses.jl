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
#univariate analysis

using FourierAnalysis

WD =  collect(dom(window_data))
WC =  collect(codom(window_data))
# using JLD2 
# @save "TempData/event_lfp_windows_1_1.jld2"  WC
##==============================================================================
#split events into pre-peri-post with 512 sample pre and post windows
#events shorter than 512 are expanded
wd_pre=GMap(x->x[498:1010], WC)
wd_peri=GMap(x->length(x)<=2532 ? x[round(Int, length(x)/2)-256:round(Int, length(x)/2)+256] : x[1010:end-1010] , WC)
wd_post=GMap(x->x[end-1010:end-498], WC)

#Freq Dom
pre_spec=GMap(x->spectra(x, 256, 512), wd_pre.(WC))
peri_spec=GMap(x->spectra(x, 256, 512), wd_peri.(WC))
post_spec=GMap(x->spectra(x, 256, 512), wd_post.(WC))
##==============================================================================


using DataFrames
spctbl = DataFrame()
for f in pre_spec[1].flabels
    pre_spec[i] = DataFrame(pre_spec[i], colnames = ["freq", "power"])
end

es=filter(x->agenttype(event_trial(x).condition) == Rat, events)

d=inverse(window_data_event)(es[1])

window_data_region.(d)




##==============================================================================
using Interpolations
using CairoMakie
using ProgressBars

dataregion = MOB
dataevent="Sniff"

CairoMakie.activate!(type="png")
idxs=findall(x->window_data_event(x).behavior.name == dataevent,WC) ∩ findall(x->window_data_region(x)==dataregion,WC)

f= Figure();
axs = [Axis(f[i, j]) for i = 1:3, j=1:3]

for n in tqdm(1:9)
    idx=rand(idxs)
    # test = 1:1:length(WC[idx])|>x->map(x->sin(16*2*pi*x/1010.1),x)
    S=TFamplitude(WC[idx], 256, 0, 1)
    itp=interpolate(rotl90(S.y), BSpline(Cubic(Line(OnGrid()))))
    x=1:1:size(S.y, 2)
    y=exp2.(range(0, log2(252.99), 100))
    
    heatmap!(axs[n],[itp(x1,y1) for x1 = x, y1 = y]; colormap = :jet)
    vlines!(axs[n], 1010; color=:orange, linewidth=2, linestyle=:dash)
    vlines!(axs[n],length(x)-1010;color=:red, linewidth=2, linestyle=:dash)
    axs[n].xticks = ([1010, length(x)-1010], ["start", "end"])
    axs[n].yticks=([21, 100],  ["8", "127"])
    e=window_data_event(WC[idx]).behavior.name
    reg=window_data_region(WC[idx]).name
    axs[n].title = "$e $reg"
end
f


##==============================================================================
# #Time Freq Dom
# using ProgressBars
# # lfp_analytic_dict=Dict{LFPRecording, TFAnalyticSignal}() 
# # p=Planner(plan_exhaustive, 10.0, 1024, Float64)
# # for l in tqdm(image(window_lfp, WD))
# #     lfp_analytic_dict[l] = TFanalyticsignal(l.lfp, 512, 2^10, 1; fmax=128)
# # end


# using ProgressBars
# WCshort=filter(x->length(x) <= 2^13, WC)
# wd_analytic_dict=Dict{Vector{Float64},TFAnalyticSignal}()


# for w in @tdqm(WCshort)
#     try
#         wd_analytic_dict[w] = iseven(length(w)) ? TFanalyticsignal(w, 256,  0, 1) : TFanalyticsignal(w[1:end-1], 256,  0, 1)    
#     catch
#         println(inverse(window_data)(w))
#     end
# end




# save("./Data/wd_analytic.jld2", "wdtf",  wd_analytic_dict)