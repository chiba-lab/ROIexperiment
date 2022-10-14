include("./FamiliarityTable.jl")
##==============================================================
elist = unique(eventtbl.event)
event_receiver_name = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_name[i]) for i in 1:size(eventtbl, 1)]))
event_receiver_exp_num = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_exp_num[i]) for i in 1:size(eventtbl, 1)]))
event_trial_exp_num = GMap(Dict([(eventtbl.event[i] => eventtbl.trial_exp_num[i]) for i in 1:size(eventtbl, 1)]))
##==============================================================================

agmap=GMap(Dict([(Object, "Object"), (Rat, "Rat"), (Robot, "Robot"), (Missing, "Missing")]))

spctbl = DataFrame()
spctbl.rat = map(x -> window_lfp(inv(window_data)(x)).rat.name, WC)
spctbl.region = map(x -> window_data_region(x).name, WC)
spctbl.trial_type = map(x -> event_trial(window_data_event(x)).condition.name, WC)
spctbl.behavior_type = map(x -> window_data_event(x).behavior.name, WC)
spctbl.agent_type = map(x -> agenttype(event_trial(window_data_event(x)).condition), WC)|>x->agmap.(x)
spctbl.trial_exp_num = map(x -> event_trial_exp_num(window_data_event(x)), WC)
spctbl.receiver_exp_num = map(x -> event_receiver_exp_num(window_data_event(x)), WC)
##==============================================================
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

using FourierAnalysis

peri_spec = GMap(x -> spectra(x, 256, 512), wd_peri.(WC))

##==========================================================
using MixedModels
specs = peri_spec.(wd_peri.(WC))
spctbl.spec = map(x -> x.y./sum(x.y), specs)
fr = FourierAnalysis.fres(256, 512) .* (1010.1 / 256)
spctbl.freq = fill(round.(collect(fr:fr:fr*256), digits=1), (length(WC)))
spctbl1 = flatten(spctbl, [:spec, :freq])
# spctbl1.id=collect(1:size(spctbl1, 1))
# fm0 = @formula(spec ~ (freq|region))
# M0 = fit(MixedModel, fm0, spctbl1, contrasts=Dict([:freq => DummyCoding()]))
#  :region => DummyCoding(), :trial_type => DummyCoding()
 #+zerocorr(0+fulldummy(freq) | agent_type&region)+zerocorr(0+fulldummy(freq) | receiver_fam&region))

import Missings
# @auto_hash_equals struct NullAgent <: Agent
#     name
# end
# mv=Missings.coalesce(spctbl1.agent_type, NullAgent)
# findall(x->isa(x,NullAgent), mv)
##==========================================================
using MixedModels

sbs = @subset(spctbl1, :freq .< 60, :behavior_type .== "Sniff", :agent_type .!= "Missing", isa.(:receiver_exp_num, Int))


sbs.receiver_fam = map(x -> x <= 2 ? "Novel" : "Familiar", sbs.receiver_exp_num)


# freq & region & agent_type + freq & region & receiver_fam
fm1 = @formula(spec ~0 + freq&region&agent_type&receiver_fam+ zerocorr(0 + fulldummy(freq) | region&rat)) #+zerocorr(0+fulldummy(freq) | agent_type&region)+zerocorr(0+fulldummy(freq) | receiver_fam&region))
# fm2 = @formula(spec ~0+freq & region + zerocorr(0 + freq&region | rat)) #+zerocorr(0+f.ulldummy(freq) | agent_type&region)+zerocorr(0+fulldummy(freq) | receiver_fam&region))
# fm2 = @formula(spec ~ 0 + freq & region + freq & region & receiver_fam & agent_type + zerocorr(0+fulldummy(freq) |region & agent_type) + zerocorr(0+fulldummy(freq) | region & receiver_fam) + zerocorr(0 + fulldummy(freq) | region & rat))
M = fit(MixedModel, fm1, sbs, contrasts=Dict([:freq => DummyCoding(), :region => DummyCoding(), :agent_type => DummyCoding(), :receiver_fam => DummyCoding()]))
# M2 = fit(MixedModel, fm2, sbs, contrasts=Dict([:freq => DummyCoding(), :region => DummyCoding(), :agent_type => DummyCoding(;base=Object), :receiver_fam => DummyCoding(;base="Familiar")]))
# M2= fit(MixedModel, fm2, sbs, contrasts=Dict([:freq => DummyCoding()]));


sbs_new=copy(sbs)

# sbs_new.rat = fill("EigenRat", size(sbs_new, 1))


sbs_new.spec_fit = fitted(M)

ss=@subset(sbs, isequal.(:agent_type,"Object"), isequal.(:region, "Ca2"))
gdf = groupby(sbs_new, [:region, :agent_type, :receiver_fam]) 
gdf["Ca2", "Object" ,"Novel"]
# agmap=GMap(Dict([(Object, "Object"), (Rat, "Rat"), (Robot, "Robot")]))
# sbs_new.agent_type = map(x -> agmap(x), sbs.agent_type)


##==========================================================
# df = DataFrame(coeftable(M))

# splitnames = split.(df.Name, " & ")
# df.id=collect(1:size(df, 1))
# tblvec=map(x ->DataFrame([y[1]=>y[2] for y = split.(x, ": ")]...), splitnames)
# dfsplit=vcat(tblvec...,cols = :union)
# dfsplit.id=collect(1:size(df, 1))
# dffull=outerjoin(dfsplit, df[:,2:end]; on = :id)
# dffull.freq=map(x->parse(Float64, x), dffull.freq)
# dffull.agent_type=collect(Missings.replace(dffull.agent_type, "Object"))
# dffull.receiver_fam=collect(Missings.replace(dffull.receiver_fam, "Familiar"))
# gdf = groupby(dffull, [:freq, :region, :agent_type, :receiver_fam])







    

# df.splitnames = split.(df.Name, "& ")
# df2=DataFrame()
# tblvec=map(x ->DataFrame([y[1]=>y[2] for y = split.(x, ": ")]...) , df.splitnames)
# vcat(tblvec...,cols = :union)

# [x->tbl[x].id=x for x=1:size(tbl, 1)]
# dfbase = @subset(df, length.(:splitnames) .== 2)
# dffam = @subset(df, length.(:splitnames) .== 3, occursin.("receiver_fam", :coefnames))
# dfagnt = @subset(df, length.(:splitnames) .== 3, occursin.("agent_type", :coefnames))
# dffull = @subset(df, length.(:splitnames) .== 4)

# dfbase.cnames = map(x -> [y[1] for y = split.(x, ":")], dfbase.splitnames)
# dfbase.cvals = map(x -> [y[2] for y = split.(x, ":")], dfbase.splitnames)
# dfbase.freq = map(x -> parse(Float64, x[1]), dfbase.cvals)
# dfbase.region = map(x -> x[2], dfbase.cvals)

# dffam.cnames = map(x -> [y[1] for y = split.(x, ":")], dffam.splitnames)
# dffam.cvals = map(x -> [y[2] for y = split.(x, ":")], dffam.splitnames)
# dffam.freq = map(x -> parse(Float64, x[1]), dffam.cvals)
# dffam.region = map(x -> x[2], dffam.cvals)
# dffam.familiarity = map(x -> x[3], dffam.cvals)


VarCorr(M)

##==========================================================
using CairoMakie, RollingFunctions

# dcc = Dict([(Object, :red), (Rat, :green), (Robot, :blue), ("Object", :red), ("Rat", :green), ("Robot", :blue)])


# ax = Axis(f[1,1])
# axs = [Axis(f[1,1][i, j], yscale=log10, xscale=log10, yminorticksvisible=true, yminorgridvisible=true) for i = 1:2, j = 1:3]
Reg=sort!(collect(unique(sbs_new.region)))
agmap=GMap(Dict([(Object, "Object"), (Rat, "Rat"), (Robot, "Robot")]))
Agnt=sort!(collect(unique(sbs_new.agent_type)))
Fam=sort!(collect(unique(sbs_new.receiver_fam)))
gdf = groupby(sbs_new, [:region, :agent_type, :receiver_fam]) #|> x -> @combine(x, mspec = mean(:spec)) |> x -> groupby(x, [:region, :agent_type, :receiver_fam])

# styledict_fam = Dict([("Object",(:green,.9)), ("Rat",(:darkred,.9)), ("Robot", (:darkblue,.9))])
styledict_nov = Dict([("Object",(:green,.9)), ("Rat",(:red,.9)), ("Robot", (:blue,.9))])




f = Figure();







# axs =vec([Axis(f[i, j], yscale=log10, xscale=log10, yminorgridvisible=true) for i = 1:1, j = 1:1])



for i in 1:1
    f = Figure()
    axs = [Axis(f[1, 1], yscale=log10, xscale=log10, yminorgridvisible=true)]



    for j in 1:1
        # for k in 1:length(Fam)
        pdat_fam = sort(gdf[("Ca2", "Object", "Familiar")], :freq)
        pdat_nov = sort(gdf[("Ca2", "Object", "Novel")], :freq)
        # pdat_nov = sort(gdf[(Reg[i], Agnt[j], "Novel")], :freq)
        # pdat_fam = sort(gdf[(Reg[i], Agnt[j], "Familiar")], :freq)
        # pdat_nov = sort(gdf[(Reg[i], Agnt[j], "Novel")], :freq)
    
        # s=runmean(collect(pdat.spec), 3)
        scatter!(axs[1], pdat_fam.freq, pdat_fam.spec, label=Agnt[j] * " " * "Familiar")#; linewidth=2, color=styledict_nov[Agnt[j]], transparency=true)
        scatter!(axs[1], pdat_nov.freq, pdat_nov.spec, label=Agnt[j] * " " * "Novel")#; linewidth=2, color=styledict_nov[Agnt[j]], linestyle=:dash, transparency=true)
    
        # cv = [(pdat_fam.mspec[s], pdat_nov.mspec[s]) for s in 1:size(pdat_fam, 1)]
        # band!(axs[1], pdat_nov.freq, minimum.(cv), maximum.(cv); color=(styledict_nov[Agnt[j]][1], 0.1))
    
        # if i == 1 && j == 3
        #     axislegend(ax; position=:lb)
        # end
    
        # ax.title = Agnt[j] * " " * Fam[i]
        # if i == 1 && j == 1
        #     ax.ylabel = "Spectral Power"
        #     ax.yticks = ([10^(-1), 10^(-2), 10^(-3)], ["-1", "-2", "-3"])
        # else
        axs[1].yticks = ([], [])
        axs[1].xticks = ([3, 10, 30, 60, 120], ["3", "10", "30", "60", "120"])
        axs[1].title = "Sniff Spectrum:" * "Ca2"
    
        # axs[i].aspect = AxisAspect(1)
    
        # end
    
    
    
    
        # end
    
    
    end
    axislegend(axs[1])
    # save("sniff_sample_plots_" * Reg[i] * ".png", f)
end
f
# f[1, 4] = Legend(f, ax, "Region");

 pdat_nov = sort(gdf[("Ca2", "Object", "Novel")], :freq)
f
##==========================================================
fams=[[],"Novel", "Familiar"]
agnts=[[],"Object", "Rat", "Robot"]

for i in 1:3
    for a in 1:4
        rsubs = @subset(dffam, :familiarity .== f, :agent_type .== a)
        sort!(rsubs, :freq)
        lines!(axs[0], rsubs.freq, rsubs.coef, label=f + " " + a; linewidth=2)
        # band!(axs[0], rsubs.freq,  max.(rsubs.coef .- rsubs.err, 0.000001), rsubs.coef .+rsubs.err)
    end




# axislegend(ax)
# band!(ax, fl, max.(m .- s, 0.0001), m .+ s, )
# xlims!(ax, fl[1], 50)
# ylims!(ax, 0.001, 1)
# ax.xticks = ([5, 10, 25, 45], ["5", "10", "25", "45"])


axs = [Axis(f[i, j], yscale=log10, xscale=log10, yminorticksvisible=true, yminorgridvisible=true) for i = 1:3, j = 1:3]


##==========================================================
sbs = @subset(spctbl1, :freq .< 150, :behavior_type .!= "Sniff", :agent_type .<: Agent, isa.(:receiver_exp_num, Int))

sbs.receiver_fam = map(x -> x <= 2 ? "Novel" : "Familiar", sbs.receiver_exp_num)
fm = @formula(spec ~ 0 + freq & region + freq & region & agent_type + freq & region & receiver_exp_num + freq & region & receiver_exp_num & agent_type + zerocorr(0 + freq & region | rat))
M = fit(MixedModel, fm, sbs, contrasts=Dict([:freq => DummyCoding()]))


##==========================================================
using ProgressBars
codict = Dict()

inef = inverse(window_data_event)
using BenchmarkTools


Threads.@threads for e in tqdm(unique(eventtbl.event))
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
            codict[e, rWC[i], rWC[j]] = map(x -> x[i, j], cWC)
        end
    end
end


codict[unique(eventtbl.event)[483], "Insula", "Ca2"]
event_region_coherence = GMap(codict)
