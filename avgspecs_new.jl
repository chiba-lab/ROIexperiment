include("./FamiliarityTable.jl")
##==============================================================
elist = unique(eventtbl.event)
event_receiver_name = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_name[i]) for i in 1:size(eventtbl, 1)]))
event_receiver_exp_num = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_exp_num[i]) for i in 1:size(eventtbl, 1)]))
event_trial_exp_num = GMap(Dict([(eventtbl.event[i] => eventtbl.trial_exp_num[i]) for i in 1:size(eventtbl, 1)]))
##==============================================================================
spctbl = DataFrame()
spctbl.event = map(x -> window_data_event(x), WC)
spctbl.rat = map(x -> window_lfp(inv(window_data)(x)).rat.name, WC)
spctbl.region = map(x -> window_data_region(x).name, WC)
spctbl.trial_type = map(x -> event_trial(window_data_event(x)).condition.name, WC)
spctbl.behavior_type = map(x -> window_data_event(x).behavior.name, WC)
spctbl.agent_type = map(x -> agenttype(event_trial(window_data_event(x)).condition), WC)
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
using MixedModels, StatsBase
specs = peri_spec.(wd_peri.(WC))
spctbl.spec = map(x -> x.y ./ sum(x.y), specs)
fr = FourierAnalysis.fres(256, 512) .* (1010.1 / 256)
spctbl.freq = fill(round.(collect(fr:fr:fr*256), digits=1), (length(WC)))
idx = findlast(x -> x < 150, spctbl.freq[1])
spctbl.freq = map(x -> x[1:idx], spctbl.freq)
spctbl.spec = map(x -> x[1:idx], spctbl.spec)
spctbl = @subset(spctbl, :behavior_type .== "Sniff", :agent_type .<: Agent, isa.(:receiver_exp_num, Int))
gdf = groupby(spctbl, :region)

zt_mob = fit(ZScoreTransform, hcat(collect(gdf[(region="MOB",)].spec)...), dims=2)
cov_mob = cov(StatsBase.transform(zt_mob, hcat(collect(gdf[(region="MOB",)].spec)...))')
zt_amg = fit(ZScoreTransform, hcat(collect(gdf[(region="Amygdala",)].spec)...), dims=2)
cov_amg = cov(StatsBase.transform(zt_mob, hcat(collect(gdf[(region="Amygdala",)].spec)...))')
zt_ca = fit(ZScoreTransform, hcat(collect(gdf[(region="Ca2",)].spec)...), dims=2)
cov_ca = cov(StatsBase.transform(zt_mob, hcat(collect(gdf[(region="Ca2",)].spec)...))')
zt_ins = fit(ZScoreTransform, hcat(collect(gdf[(region="Insula",)].spec)...), dims=2)
cov_ins = cov(StatsBase.transform(zt_mob, hcat(collect(gdf[(region="Insula",)].spec)...))')
function fun(s, r)

    if r == "MOB"
        return StatsBase.transform(zt_mob, s)
    elseif r == "Amygdala"
        return StatsBase.transform(zt_amg, s)
    elseif r == "Ca2"
        return StatsBase.transform(zt_ca, s)
    elseif r == "Insula"
        return StatsBase.transform(zt_ins, s)
    end
end

function rfun(s, r)
    # s = reshape(s, 256, 1)
    if r == "MOB"
        return StatsBase.reconstruct(zt_mob, s)
    elseif r == "Amygdala"
        return StatsBase.reconstruct(zt_amg, s)
    elseif r == "Ca2"
        return StatsBase.reconstruct(zt_ca, s)
    elseif r == "Insula"
        return StatsBase.reconstruct(zt_ins, s)
    end
end


function cfun(s, r)
    # s = reshape(s, 256, 1)
    if r == "MOB"
        return cov_mob * s
    elseif r == "Amygdala"
        return cov_amg * s
    elseif r == "Ca2"
        return cov_ca * s
    elseif r == "Insula"
        return cov_ins * s
    end
end



# spctbl.spec=fun.(spctbl.spec,spctbl.region)



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


sbs = @subset(spctbl1, :freq .< 150)

sbs.receiver_fam = map(x -> x <= 2 ? "Novel" : "Familiar", sbs.receiver_exp_num)



# freq & region & agent_type + freq & region & receiver_fam
# sbs_mob = @subset(sbs, :region .== "MOB")
# sbs_amg = @subset(sbs, :region .== "Amygdala")
# sbs_ca = @subset(sbs, :region .== "Ca2")
# sbs_ins = @subset(sbs, :region .== "Insula")
fm1 = @formula(spec ~ 0 + freq & region + freq & region & agent_type + freq & region & receiver_fam + freq & region & agent_type & receiver_fam + zerocorr(0 + fulldummy(freq) | region & rat)) #+zerocorr(0+fulldummy(freq) | agent_type&region)+zerocorr(0+fulldummy(freq) | receiver_fam&region))
# fm2 = @formula(spec ~0+freq & region + zerocorr(0 + freq&region | rat)) #+zerocorr(0+fulldummy(freq) | agent_type&region)+zerocorr(0+fulldummy(freq) | receiver_fam&region))
# fm2 = @formula(spec ~ 0 + freq & region + freq & region & receiver_fam & agent_type + zerocorr(0+fulldummy(freq) |region & agent_type) + zerocorr(0+fulldummy(freq) | region & receiver_fam) + zerocorr(0 + fulldummy(freq) | region & rat))
M = fit(MixedModel, fm1, sbs, contrasts=Dict([:freq => DummyCoding(), :agent_type => DummyCoding(; base=Object)]))#:receiver_fam => DummyCoding(;base="Familiar")]))

# M_mob = fit(MixedModel, fm1, sbs_mob, contrasts=Dict([:freq => DummyCoding(), :agent_type => DummyCoding(;base=Object)]))#:receiver_fam => DummyCoding(;base="Familiar")]))
# M_amg = fit(MixedModel, fm1, sbs_amg, contrasts=Dict([:freq => DummyCoding(), :agent_type => DummyCoding(; base=Object)])) # :receiver_fam => DummyCoding(; base="Familiar")]))
# M_ca = fit(MixedModel, fm1, sbs_ca, contrasts=Dict([:freq => DummyCoding(), :agent_type => DummyCoding(; base=Object)]))#:receiver_fam => DummyCoding(; base="Familiar")]))
# M_ins = fit(MixedModel, fm1, sbs_ins, contrasts=Dict([:freq => DummyCoding(), :agent_type => DummyCoding(; base=Object)]))# :receiver_fam => DummyCoding(; base="Familiar")]))
# M2 = fit(MixedModel, fm2, sbs, contrasts=Dict([:freq => DummyCoding(), :region => DummyCoding(), :agent_type => DummyCoding(;base=Object), :receiver_fam => DummyCoding(;base="Familiar")]))
# M2= fit(MixedModel, fm2, sbs, contrasts=Dict([:freq => DummyCoding()]));
# DataFrame(coeftable(M))
##==========================================================
# sbs_mob.spec=fitted(M_mob)
# sbs_amg.spec=fitted(M_amg)
# sbs_ca.spec=fitted(M_ca)
# sbs_ins.spec=fitted(M_ins)

# sbs_new=vcat(sbs_mob, sbs_amg, sbs_ca, sbs_ins)
sbs_new = copy(sbs)

# 
# sbs_new.spec = fitted(M)


# sbs_new.rat = fill("EigenRat", size(sbs_new, 1))



agmap = GMap(Dict([(Object, "Object"), (Rat, "Rat"), (Robot, "Robot")]))

sbs_new.agent_type = map(x -> agmap(x), sbs.agent_type)
##==========================================================
sbs_new_grouped = groupby(sbs_new, [:rat, :event, :region, :trial_type, :behavior_type, :agent_type, :trial_exp_num, :receiver_exp_num]) |> x -> @combine(x, :respiratory = mean(:spec[2 .<= collect(:freq) .< 6]), :theta = mean(:spec[6 .<= collect(:freq) .< 12]), :beta = mean(:spec[12 .<= collect(:freq) .< 40]), :low_gamma = mean(:spec[50 .<= collect(:freq) .< 60]), :high_gamma = mean(:spec[70 .<= collect(:freq) .< 100]))
sbs_new_grouped.trial_fam = map(x -> isa(x, Int) ? (x <= 2 ? "Novel" : "Familiar") : "missing", sbs_new_grouped.trial_exp_num)

sbs_new_grouped.receiver_fam = map(x -> isa(x, Int) ? (x <= 2 ? "Novel" : "Familiar") : "missing", sbs_new_grouped.receiver_exp_num)
sbs_new_grouped = @select(sbs_new_grouped, :rat, :region, :trial_type, :behavior_type, :agent_type, :trial_fam, :receiver_fam, :respiratory, :theta, :beta, :low_gamma, :high_gamma)
CSV.write("./TempData/FreqBandsTbl.csv", sbs_new_grouped)

# @subset(sbs_new_grouped, :agent_type .== "Robot", :region .== "MOB", :behavior_type .== "Sniff", :receiver_fam .== "Novel")
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

unique(sbs_new_grouped.low_gamma)
gdf = groupby(sbs_new, [:event, :region])
mean(gdf[2].spec[70 .<= gdf[2].freq .< 100])

# |> x -> @combine(x, :respiratory = mean(:spec[2 .<= collect(:freq) .< 6]), :theta = y = mean(:spec[6 .<= collect(:freq) .< 12]), :beta = mean(:spec[12 .<= collect(:freq) .< 40]), :low_gamma = mean(:spec[50 .<= collect(:freq) .< 60]), :high_gamma = mean(:spec[70 .<= collect(:freq) .< 100]))






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



##==========================================================
using CairoMakie, RollingFunctions

# sbs_new = @subset(spctbl1, :freq .< 150, :behavior_type .== "Sniff", :agent_type .<: Agent, isa.(:receiver_exp_num, Int))

# sbs_new.receiver_fam = map(x -> x <= 2 ? "Novel" : "Familiar", sbs.receiver_exp_num)
# agmap=GMap(Dict([(Object, "Object"), (Rat, "Rat"), (Robot, "Robot")]))
# sbs_new.agent_type = map(x -> agmap(x), sbs.agent_type)

# dcc = Dict([(Object, :red), (Rat, :green), (Robot, :blue), ("Object", :red), ("Rat", :green), ("Robot", :blue)])


# ax = Axis(f[1,1])
# axs = [Axis(f[1,1][i, j], yscale=log10, xscale=log10, yminorticksvisible=true, yminorgridvisible=true) for i = 1:2, j = 1:3]
Reg = sort!(collect(unique(sbs_new.region)))
# agmap=GMap(Dict([(Object, "Object"), (Rat, "Rat"), (Robot, "Robot")]))
Agnt = sort!(collect(unique(sbs_new.agent_type)))
Fam = sort!(collect(unique(sbs_new.receiver_fam)))

gdf = groupby(sbs_new, [:freq, :region, :agent_type, :receiver_fam]) |> x -> @combine(x, mspec = mean(:spec), len = length(:spec)) |> x -> groupby(x, [:region, :agent_type, :receiver_fam])

# styledict_fam = Dict([("Object",(:green,.9)), ("Rat",(:darkred,.9)), ("Robot", (:darkblue,.9))])
styledict_nov = Dict([(("Rat", "Familiar"), (:red, 0.9)), (("Rat", "Novel"), (:orange, 0.9)), (("Robot", "Familiar"), (:blue, 0.9)), (("Robot", "Novel"), (:green, 0.9))])














# axs =vec([Axis(f[i, j], yscale=log10, xscale=log10, yminorgridvisible=true) for i = 1:1, j = 1:1])

Agnt2 = sort(filter(x -> x != "Object", Agnt))

for i in 1:length(Reg)

    f = Figure(resolution=(1200, 600))
    axs = [Axis(f[1, 1], xscale=log10, yscale=log10, yminorgridvisible=true), Axis(f[1, 2], xscale=log10, yminorgridvisible=true)]
    baseln = (gdf[(Reg[i], "Object", "Novel")].mspec .+ gdf[(Reg[i], "Object", "Familiar")].mspec) / 2
    hlines!(axs[2], 0; linewidth=3, color=:grey3, transparency=true, linestyle=:dash)
    hlines!(axs[2], log10(1.64), linewidth=2, color=:grey1, transparency=true, linestyle=:dot, label="p<0.05")
    hlines!(axs[2], -log10(1.64), linewidth=2, color=:grey1, transparency=true, linestyle=:dot, label="p<0.05")



    for j in 1:length(Agnt2)
        # for k in 1:length(Fam)
        pdat_fam = sort(gdf[(Reg[i], Agnt2[j], "Familiar")], :freq)
        pdat_nov = sort(gdf[(Reg[i], Agnt2[j], "Novel")], :freq)

        # s=runmean(collect(pdat.spec), 3)
        lines!(axs[1], pdat_fam.freq, pdat_fam.mspec, label=Agnt2[j] * " " * "Familiar"; linewidth=3, color=styledict_nov[(Agnt2[j], "Familiar")], transparency=true)
        lines!(axs[1], pdat_nov.freq, pdat_nov.mspec, label=Agnt2[j] * " " * "Novel"; linewidth=3, color=styledict_nov[(Agnt2[j], "Novel")], linestyle=:dash, transparency=true)
        lines!(axs[2], pdat_fam.freq, cfun(fun(pdat_fam.mspec, Reg[i]), Reg[i]) .* sqrt.(pdat_fam.len) |> x -> sign.(x) .* log10.(abs.(x)), label=Agnt2[j] * " " * "Familiar"; linewidth=3, color=styledict_nov[(Agnt2[j], "Familiar")], transparency=true)
        lines!(axs[2], pdat_nov.freq, cfun(fun(pdat_nov.mspec, Reg[i]), Reg[i]) .* sqrt.(pdat_nov.len) |> x -> sign.(x) .* log10.(abs.(x)), label=Agnt2[j] * " " * "Novel"; linewidth=3, color=styledict_nov[(Agnt2[j], "Novel")], linestyle=:dash, transparency=true)




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
        # axs[1].yticks = ([], [])
        axs[1].xticks = ([3, 10, 30, 60, 120], ["3", "10", "30", "60", "120"])
        axs[1].title = "Sniff Spectrum: Absolute\n" * Reg[i]
        axs[2].xticks = ([3, 10, 30, 60, 120], ["3", "10", "30", "60", "120"])
        axs[2].title = "Sniff Spectrum: Relative\n" * Reg[i]

        xlims!(axs[1], (minimum(pdat_fam.freq), maximum(pdat_fam.freq)))

        xlims!(axs[2], (minimum(pdat_fam.freq), maximum(pdat_fam.freq)))

        axs[1].xlabel = "Frequency (Hz)"
        axs[1].ylabel = "Power"
        axs[2].xlabel = "Frequency (Hz)"
        axs[2].ylabel = "sgn(Z)*log(|Z|)"
        # axs[i].aspect = AxisAspect(1)

        # end




        # end


    end
    axislegend(axs[1])
    save("sniff_plots_" * Reg[i] * ".pdf", f)
end
##==========================================================
using PDFmerger, Glob
fn = glob("sniff_plots_*.pdf")
merge_pdfs(fn, "sniff_plots.pdf", cleanup=true)
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
