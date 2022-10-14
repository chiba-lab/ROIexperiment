
include("./FamiliarityTable.jl")
using FourierAnalysis
##==============================================================
elist = unique(eventtbl.event)
event_receiver_name = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_name[i]) for i in 1:size(eventtbl, 1)]))
event_receiver_exp_num = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_exp_num[i]) for i in 1:size(eventtbl, 1)]))
event_trial_exp_num = GMap(Dict([(eventtbl.event[i] => eventtbl.trial_exp_num[i]) for i in 1:size(eventtbl, 1)]))
makeeven(x) = length(x) % 2 == 0 ? x : x[1:end-1]
function peritrim(x)
    if length(x) <= 2532
        return makeeven(x[round(Int, length(x) / 2)-256:round(Int, length(x) / 2)+256])
    elseif (length(x) - 2020) > 2^14
        return makeeven(x[1010:1010+2^14])
    else
        return makeeven(x[1010:end-1010])
    end
end

wd_peri = GMap(x -> peritrim(x), WC)


length.(wd_peri.(WC)) .% 2
##==============================================================
using ProgressBars
codict = Dict()

inef = inverse(window_data_event)
using BenchmarkTools
using UncertainData
lbe = [2, 6, 12, 40, 60, 70, 100, 150] .* (256 / 1010.1)
fr = FourierAnalysis.fres(256, 512) .* (1010.1 / 256)
fl = collect(fr:fr:fr*256)


for e in tqdm(unique(eventtbl.event))
    # e = ev[ei]
    # for e in unique(eventtbl.event[1:500])
    eWC = inef(e) |> x -> isa(x[1], Vector) ? x : [x]
    rWC = map(x -> window_data_region(x).name, eWC)
    cWC = collect(coherence(hcat(wd_peri.(eWC)...), 256, 512).y)
    # if length(rWC) == 4
    #     print(ei)
    # end
    for i in 1:length(rWC)
        for j in 1:length(rWC)
            codict[e, rWC[i], rWC[j]] = map(x -> x[i, j], cWC) #|> x -> bin(mean, 1:5:150, fl, x)
        end
    end
end


codict[unique(eventtbl.event)[483], "Amygdala", "Ca2"]
# event_region_coherence = GMap(codict)
# event_region_coherence = []
##==============================================================
##==============================================================
#event_region_coherence((unique(eventtbl.event)[483], "Insula", "Ca2"))


##==============================================================
rlist = ["MOB", "Amygdala", "Ca2", "Insula"]
blist = ["respiratory", "theta", "beta", "low_gamma", "high-gamma"]
bdict = Dict(zip(blist, [(2, 6), (6, 12), (50, 40), (70, 100)]))
cotbl = DataFrame()
cvec = collect(keys(codict))
cmat = hcat(map(x -> [x...], cvec)...)
cotbl.event = vec(cmat[1, :])
cotbl.trial_type = map(e -> event_trial(e).condition.name, cotbl.event)
cotbl.behavior_type = map(e -> e.behavior.name, cotbl.event)
agmap = GMap(Dict([(Object, "Object"), (Rat, "Rat"), (Robot, "Robot"), (Missing, "missing")]))
cotbl.agent_type = agmap.(map(e -> agenttype(event_trial(e).condition), cotbl.event))
cotbl.trial_exp_num = map(e -> event_trial_exp_num(e), cotbl.event)
cotbl.receiver_exp_num = map(e -> event_receiver_exp_num(e), cotbl.event)
cotbl.trial_fam = map(x -> isa(x, Int) ? (x <= 2 ? "Novel" : "Familiar") : "missing", cotbl.trial_exp_num)
cotbl.receiver_fam = map(x -> isa(x, Int) ? (x <= 2 ? "Novel" : "Familiar") : "missing", cotbl.receiver_exp_num)
cotbl.fam = map(x -> x[1] == "missing" ? (x[2] == "missing" ? "missing" : x[2]) : x[1], [(cotbl.trial_fam[i], cotbl.receiver_fam[i]) for i = 1:size(cotbl, 1)])
cotbl.region_1 = vec(cmat[2, :])
cotbl.region_2 = vec(cmat[3, :])
cotbl.region_pair = map(x -> x[1] * "-" * x[2], [(cotbl.region_1[i], cotbl.region_2[i]) for i = 1:size(cotbl, 1)])
cotbl.co_vec = map(x -> codict[x], cvec)
cotbl.freq = fill(fl, length(cotbl.co_vec))
cotbl.respiratory = map(x -> mean(x[2 .<= fl .< 6]), cotbl.co_vec)
cotbl.theta = map(x -> mean(x[6 .<= fl .< 12]), cotbl.co_vec)
cotbl.beta = map(x -> mean(x[12 .<= fl .< 40]), cotbl.co_vec)
cotbl.low_gamma = map(x -> mean(x[50 .<= fl .< 60]), cotbl.co_vec)
cotbl.high_gamma = map(x -> mean(x[70 .<= fl .< 100]), cotbl.co_vec)



##==============================================================
# using CSV, FileIO, Missings
# cmctbl_reduc = cmctbl[:, Not([:event, :co_vec, :freq, :trial_exp_num, :receiver_exp_num])]
# @subset!(cmctbl_reduc, :region_1 .∈ Ref(("MOB", "Ca2")), :region_2 .∈ Ref(("MOB", "Ca2")))
# CSV.write("ca-mob_coherence.csv", cmctbl_reduc)

##=============================================================

using JLD2
save("cotbl.jld2", "cotable", cotbl)
##==============================================================

rlist = ["MOB", "Amygdala", "Ca2", "Insula"]
sort!(rlist)
rrlist = [rlist[i], rlist[j]] for i in 1:length(rlist) for j in 1:length(rlist) if i < j]

l = 14

sbsco = @subset(cotbl, :agent_type .!= "missing", :fam .!= "missing", :region_1 .∈ Ref(rlist), :region_2 .∈ Ref(rlist))
@select!(sbsco, :event, :behavior_type, :agent_type, :fam, :region_1, :region_2, :respiratory, :theta, :beta, :low_gamma, :high_gamma)
giantstate = Dict()
giantstateaxis = Dict()
for e in tqdm(unique(eventtbl.event))
    for i in 2:length(rlist)
        for j in 1:i-1
            if haskey(codict, (e, rlist[i], rlist[j]))
                if haskey(giantstate, e)
                    giantstate[e] = vcat(giantstate[e]..., codict[e, rlist[i], rlist[j]])
                else
                    giantstate[e] = codict[e, rlist[i], rlist[j]]
                end

            else
                if haskey(giantstate, e)
                    giantstate[e] = vcat(giantstate[e]..., fill(0, l))
                else
                    giantstate[e] = fill(0, l)
                end
            end
        end
    end
end

giantstate[unique(eventtbl.event)[200]]
gsmat = hcat(map(x -> giantstate[x], elist)...)
##=================================
event_receiver_name = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_name[i]) for i in 1:size(eventtbl, 1)]))
event_receiver_exp_num = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_exp_num[i]) for i in 1:size(eventtbl, 1)]))
event_trial_exp_num = GMap(Dict([(eventtbl.event[i] => eventtbl.trial_exp_num[i]) for i in 1:size(eventtbl, 1)]))
##==============================================================================
clusttbl = DataFrame()
clusttbl.trial_type = map(e -> event_trial(e).condition.name, elist)
clusttbl.behavior_type = map(e -> e.behavior.name, elist)
agmap = GMap(Dict([(Object, "Object"), (Rat, "Rat"), (Robot, "Robot"), (Missing, "missing")]))
clusttbl.agent_type = agmap.(map(e -> agenttype(event_trial(e).condition), elist))
clusttbl.trial_exp_num = map(e -> event_trial_exp_num(e), elist)
clusttbl.receiver_exp_num = map(e -> event_receiver_exp_num(e), elist)
clusttbl.trial_fam = map(x -> isa(x, Int) ? (x <= 2 ? "Novel" : "Familiar") : "missing", clusttbl.trial_exp_num)
clusttbl.receiver_fam = map(x -> isa(x, Int) ? (x <= 2 ? "Novel" : "Familiar") : "missing", clusttbl.receiver_exp_num)
clusttbl.fam = map(x -> x[1] == "missing" ? (x[2] == "missing" ? "missing" : x[2]) : x[1], [(clusttbl.trial_fam[i], clusttbl.receiver_fam[i]) for i in 1:size(clusttbl, 1)])
clusttbl.coh_vec = map(e -> giantstate[e], elist)
sbs = @subset(clusttbl, :agent_type .!= "missing", :fam .!= "missing")

using StatsModels, StatsAPI
using MultivariateStats

# fm = @formula(coh_vec ~ 0 + behavior_type * agent_type * fam)

# DM = apply_schema(fm, schema(clusttbl))

# X = modelmatrix(DM, clusttbl)'

Y = hcat(sbs.coh_vec...)


##==============================================================


# C = MultivariateStats.fit(CCA, X', Y)
# using CairoMakie
# f = Figure();
# ax = Axis(f[1, 1])
# CairoMakie.heatmap!(ax, C)
# ax.xticks = ([collect(1:84)]..., axislbs)
# f
##==============================================================
# fm = @formula(coh_vec ~0+ behavior_type&agent_type&fam)
Xtbl = @select(sbs, :behavior_type, :agent_type, :fam)
ux = [collect(unique(Xtbl)[i, :]) for i in 1:size(unique(Xtbl), 1)]
ux_b = [unique(Xtbl.behavior_type)[i] for i in 1:4]


# DM = apply_schema(fm, schema(clusttbl))

# X = modelmatrix(DM, clusttbl)
# ux = unique(X; dims=1)
# ux = [ux[i, :] for i in 1:size(ux, 1)]
# XX = [X[i, :] for i in 1:size(X, 1)]

classdict = Dict(zip(ux, [1:1:length(ux)...]))
classdict_b = Dict(zip(ux_b, [1:1:length(ux_b)...]))
gc = GMap(classdict)

xc = [classdict[[Xtbl[i, :behavior_type], Xtbl[i, :agent_type], Xtbl[i, :fam]]] for i = 1:size(Xtbl, 1)]
xc_b = [classdict_b[Xtbl[i, :behavior_type]] for i = 1:size(Xtbl, 1)]


MC = MultivariateStats.fit(MulticlassLDA, 24, Y, xc_b)


classmeans(MC)
projgs = MultivariateStats.transform(MC, Y)

using GLMakie

f = Figure();
ax = LScene(f[1, 1])
p = meshscatter!(ax, projgs[1, :], projgs[2, :], projgs[3, :], markersize=0.005, color=xc_b)


# ax.xticks = ([collect(1:84)]..., axislbs)

f
##==============================================================
classmeans(MC)
withclass_scatter(MC.stats)
##==============================================================

using GLMakie

bc = [mean(collect(1:10:150)[i-1:i]) for i = 2:length(collect(1:10:150))]

axislbs = [rlist[i] * "_" * rlist[j] * "_" * string(bc[k]) for i in 1:length(rlist) for j in 1:i-1 for k in 1:length(bc)]


f = Figure();
ax = Axis(f[1, 1])
heatmap!(ax, collect(1:84), collect(1:84), cmat)
ax.xticks = ([collect(1:84)]..., axislbs)
f

using MultivariateStats

M = fit(PCA, gsmat)

projgs = MultivariateStats.transform(M, gsmat)

using GLMakie

styledict_nov = Dict([("Sniff", :red), ("Immobility", :orange), ("Groom", :blue), ("Rear", :green, 0.9)])
cvec = vcat(map(x -> styledict_nov[event_behavior(x).name], elist)...)


f = Figure();
ax = LScene(f[1, 1])
meshscatter!(ax, projgs[1, :], projgs[2, :], projgs[3, :], color=cvec)
# ax.xticks = ([collect(1:84)]..., axislbs)

f

f = Figure();
ax = Axis(f[1, 1])
heatmap!(ax, projection(M))
# ax.xticks = ([collect(1:84)]..., axislbs)

f



using Clustering, Distances, AbstractTrees, Phylo
D = Distances.pairwise(Euclidean(), gsmat, gsmat, dims=2)
projgs = classical_mds(D, 3)
tt = Tree(hc)
nu = Nonultrametric(100)

plot(rand(nu))


# print_tree(hc)
import Plots as plts, StatsPlots

plts.plot(
    plts.plot(hc, xticks=false),
    plts.heatmap(gsmat[:, hc.order], colorbar=false, xticks=(1:2983, ["$i" for i in hcl1.order])),
    layout=grid(2, 1, heights=[0.2, 0.8])
)

import Flux: onehot, onehotbatch

Y = onehotbatch(map(x -> event_behavior(x).name, elist), unique(eventtbl.behavior_type))
Y = Y .* ones(size(Y))

C = MultivariateStats.fit(CCA, gsmat, Y)
using Statistics
projgs = MultivariateStats.xtransform(C, gsmat)


using GLMakie, RecipesBase

styledict_nov = Dict([("Sniff", :red), ("Immobility", :orange), ("Groom", :blue), ("Rear", :green, 0.9)])
cvec = vcat(map(x -> styledict_nov[event_behavior(x).name], elist)...)


f = Figure();
ax = LScene(f[1, 1])
meshscatter!(ax, projgs[1, :], projgs[2, :], projgs[3, :], color=cvec)
# ax.xticks = ([collect(1:84)]..., axislbs)

f

fit(ICA, gsmat', 5)


using StatsBase
o = StatsBase.indexmap(hc.order)
nodepos = Dict(-i => (float(o[i]), 0.0) for i in hc.order)


Tree(o)
