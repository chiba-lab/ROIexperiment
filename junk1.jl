
rlist = ["MOB", "Amygdala", "Ca2", "Insula"]
l = 14

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
