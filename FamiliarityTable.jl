include("./ROIontology.jl")


ROIexp = load("./Data/ROIexp.jld2", "ROIexp");
lfp = ROIexp[:lfp_data]
events = ROIexp[:behavioral_events]
trials = ROIexp[:trials]
sessions = ROIexp[:sessions]

include("./Links.jl")

# WD = filter(x -> window_data_event(window_data(x)) ∈ dom(event_trial), collect(dom(window_data)))
WC = filter(x -> window_data_event(x) ∈ dom(event_trial), collect(codom(window_data)))

##==============================================================================

# using FourierAnalysis


# wd_peri = GMap(x -> length(x) <= 2532 ? x[round(Int, length(x) / 2)-256:round(Int, length(x) / 2)+256] : x[1010:end-1010], WC)

# peri_spec = GMap(x -> spectra(x, 256, 512), wd_peri.(WC))
# # post_spec = GMap(x -> spectra(x, 256, 512), wd_post.(WC)) 

using DataFrames, DataFramesMeta

# te=[inverse(event_trial)(x) for x in collect(codom(event_trial))]

# tr = map(te) do x
#     if isa(x, BehavioralEvent)
#         filter(x->isa(x, Agent), unique(map(x->x.receiver, [x])))|>x->map(x->x.name, x)|>x->isempty(x) ? [missing] : x
#     else
#          filter(x->isa(x, Agent), unique(map(x->x.receiver, x)))|>x->map(x->x.name, x)|>x->isempty(x) ? [missing] : x
#     end
# end

# td=Dict(zip(collect(codom(event_trial)), tr))
# makeeven(x) = length(x) % 2 == 0 ? x : x[1:end-1]
# function peritrim(x)
#     if length(x) <= 2532
#         return makeeven(x[round(Int, length(x) / 2)-256:round(Int, length(x) / 2)+256])
#     elseif (length(x) - 2020) > 2^14
#         return makeeven(x[1010:1010+2^14])
#     else
#         return makeeven(x[1010:end-1010])
#     end
# end

# wd_peri = GMap(x -> peritrim(x), WC)

# length.(wd_peri.(WC)) .% 2

# using FourierAnalysis

# peri_spec = GMap(x -> spectra(x, 256, 512), wd_peri.(WC))



eventtbl = DataFrame()
eventtbl.event = map(x -> window_data_event(x), WC)
eventtbl.rat = map(x -> window_lfp(inv(window_data)(x)).rat.name, WC)
eventtbl.session = map(x -> window_lfp(inv(window_data)(x)).session.filename, WC)
eventtbl.date = map(x -> window_lfp(inv(window_data)(x)).session.date, WC)
eventtbl.trial_onset = map(x -> event_trial(window_data_event(x)).start_time, WC)
eventtbl.trial_type = map(x -> event_trial(window_data_event(x)).condition.name, WC)
eventtbl.trial_agents = map(x -> trial_agents(event_trial(window_data_event(x))), WC)
eventtbl.onset = map(x -> inv(window_data)(x).onset, WC)
eventtbl.region = map(x -> window_data_region(x).name, WC)
eventtbl.behavior_type = map(x -> window_data_event(x).behavior.name, WC)
eventtbl.agent_type = map(x -> agenttype(event_trial(window_data_event(x)).condition), WC)
eventtbl.receiver_attribute = map(x -> window_data_event(x).receiver, WC)
# eventtbl.lfp=wd_peri.(WC)
# unstack(eventtbl, :event, :region,  :lfp)
# eventtbl.lfp=wd_peri.(WC)
# specs=peri_spec.(wd_peri.(WC))
# eventtbl.spec=map(x->x.y, specs)
# fr = FourierAnalysis.fres(256, 512) .* (1010.1 / 256)
# eventtbl.freq = fill(round.(collect(fr:fr:fr*256),digits=1), (length(WC)))
# eventtbl.respiratory = map(x -> mean(x, (3, 12) .* 256 ./ 1010.1), specs)
# eventtbl.theta = map(x -> mean(x, (5, 10) .* 256 ./ 1010.1), specs)
# eventtbl.beta = map(x -> mean(x, (15, 35) .* 256 ./ 1010.1), specs)
# eventtbl.gamma_low = map(x -> mean(x, (50, 59) .* 256 ./ 1010.1), specs)
# eventtbl.gamma_high = map(x -> mean(x, (70, 100) .* 256 ./ 1010.1), specs)

sort!(eventtbl, [:date, :trial_onset, :onset])

agids = filter(!ismissing, reduce(union, eventtbl.trial_agents))


# ssd = Dict()
# # ssd[("missing","missing")]=missing
# for aa in unique(eventtbl.rat)
#     for ai in agids
#         sbs = @subset(eventtbl, :rat .== aa, ai .∈ :trial_agents)
#         ss = sort(unique(sbs.date))
#         for jj in 1:length(ss)
#             ssd[(aa, ai, ss[jj])] = jj
#         end
#     end
# end



# sen = [map(x -> haskey(ssd, (eventtbl.rat[i], x, eventtbl.date[i])) ? ssd[(eventtbl.rat[i], x, eventtbl.date[i])] : [], eventtbl.trial_agents[i]) for i in 1:size(eventtbl, 1)]

# sess_exp_num = map(x -> filter(!isempty, x), sen)

# eventtbl.sess_exp_num = sess_exp_num

ttd = Dict()
# ttd[("missing", "missing")] = missing
for aa in unique(eventtbl.rat)
    for ai in agids
        sbt = @subset(eventtbl, :rat .== aa, ai .∈ :trial_agents)
        sbt = @select(sbt, :date, :trial_onset)
        tt = unique(sbt)
        tt = sort(tt, [:date, :trial_onset])
        for ii in 1:size(tt, 1)
            ttd[(aa, ai, tt.trial_onset[ii])] = ii
        end
    end
end

ten = [map(x -> haskey(ttd, (eventtbl.rat[i], x, eventtbl.trial_onset[i])) ? ttd[(eventtbl.rat[i], x, eventtbl.trial_onset[i])] : [], eventtbl.trial_agents[i]) for i in 1:size(eventtbl, 1)]

trial_exp_num = map(x -> filter(!isempty, x), ten)

eventtbl.trial_exp_num = trial_exp_num

function matchreceiver(n, r, rn, b)
    if isequal(b, "Sniff")

        if isempty(r)
            return missing
        end
        if ismissing(n)
            return r[1]
        end
        if isequal(n, "Novel")
            if rn == [1, 1]
                return r[1]
            else
                return r[rn.==minimum(rn)]
            end
        elseif isequal(n, "Familiar")
            if rn == [1, 1]
                return r[2]
            else
                return r[rn.==maximum(rn)]
            end
        else
            return r[1]
        end
    else
        missing
    end

end

eventtbl.receiver_name = [matchreceiver(eventtbl.receiver_attribute[i], eventtbl.trial_agents[i], eventtbl.trial_exp_num[i], eventtbl.behavior_type[i]) for i in 1:size(eventtbl, 1)]

# event_receiver_name = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_name[i]) for i in 1:size(eventtbl, 1)]))
# event_receiver_name(eventtbl.event[64])

function receiverexpnum(n, r, rn)
    if ismissing(n)
        return missing
    elseif isempty(n)
        return missing
    else
        return rn[r.==n][1]
    end
end


eventtbl.receiver_exp_num = [receiverexpnum(eventtbl.receiver_name[i], eventtbl.trial_agents[i], eventtbl.trial_exp_num[i]) for i in 1:size(eventtbl, 1)]

tmpfun(x) = (!ismissing(x) && length(x) == 1) ? x[1] : x

eventtbl.receiver_exp_num = tmpfun.(eventtbl[:, :receiver_exp_num])



eventtbl.trial_exp_num = tmpfun.(eventtbl[:, :trial_exp_num])
eventtbl.receiver_name = tmpfun.(eventtbl[:, :receiver_name])

# elist=unique(eventtbl.event)
# event_receiver_name = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_name[i]) for i in 1:size(eventtbl, 1)]))
# event_receiver_exp_num = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_exp_num[i]) for i in 1:size(eventtbl, 1)]))
# event_trial_exp_num = GMap(Dict([(eventtbl.event[i] => eventtbl.trial_exp_num[i]) for i in 1:size(eventtbl, 1)]))
# ##==============================================================================

# using JLD2

# # savedict=Dict([(:WC,WC), (:famtbl,eventtbl), (:wd_peri,wd_peri), (:peri_spec,peri_spec),(:event_region_coherence,event_region_coherence), (:event_mob_respiratory_tf,event_mob_respiratory_tf), (:event_low_gamma_tf,event_low_gamma_tf), (:event_high_gamma_tf, event_high_gamma_tf)])


# @save "./TempData//EventFamiliarity.jld2"
# # CSV.write("./TempData/EventFamiliarity.csv", eventtbl)
# ##==========================================================
# makeeven(x) = length(x) % 2 == 0 ? x : x[1:end-1]
# function peritrim(x)
#     if length(x) <= 2532
#         return makeeven(x[round(Int, length(x) / 2)-256:round(Int, length(x) / 2)+256])
#     elseif (length(x) - 2020) > 2^14
#         return makeeven(x[1010:1010+2^14])
#     else
#         return makeeven(x[1010:end-1010])
#     end
# end

# wd_peri = GMap(x -> peritrim(x), WC)

# length.(wd_peri.(WC)) .% 2

# using FourierAnalysis

# peri_spec = GMap(x -> spectra(x, 256, 512), wd_peri.(WC))


# ##==========================================================
# using ProgressBars
# codict = Dict()

# inef = inverse(window_data_event)
# using BenchmarkTools


# Threads.@threads for e in tqdm(unique(eventtbl.event))
#     # e = ev[ei]
#     # for e in unique(eventtbl.event[1:500])
#     eWC = inef(e) |> x -> isa(x[1], Vector) ? x : [x]
#     rWC = map(x -> window_data_region(x).name, eWC)
#     cWC = collect(coherence(hcat(wd_peri.(eWC)...), 256, 512).y)
#     # if length(rWC) == 4
#     #     print(ei)
#     # end
#     for i in 1:length(rWC)

#         for j in 1:length(rWC)
#             codict[e, rWC[i], rWC[j]] = map(x -> x[i, j], cWC)
#         end
#     end
# end


# codict[unique(eventtbl.event)[483], "Insula", "Ca2"]
# event_region_coherence = GMap(codict)

# ##==========================================================
# emrtf = Dict()
# Threads.@threads for e in tqdm(unique(eventtbl.event))
#     # e = ev[ei]
#     # for e in unique(eventtbl.event[1:500])
#     eWC = inef(e) |> x -> isa(x[1], Vector) ? x : [x]
#     rWC = map(x -> window_data_region(x).name, eWC)
#     idx = findfirst(x -> x == "MOB", rWC)
#     if !isnothing(idx)
#         Y = TFanalyticsignal(wd_peri(eWC[idx[1]]), 256, 0, 0.5; fmax=4)
#         emrtf[e] = vec(mean(extract(Y, (3, 12) .* 256 ./ 1010.1, :), dims=1))
#     end

# end
# emrtf[eventtbl.event[1]]

# event_mob_respiratory_tf = GMap(emrtf)

# elgtf = Dict()
# ehgtf = Dict()
# Threads.@threads for e in tqdm(unique(eventtbl.event))
#     # e = ev[ei]
#     # for e in unique(eventtbl.event[1:500])
#     eWC = inef(e) |> x -> isa(x[1], Vector) ? x : [x]
#     rWC = map(x -> window_data_region(x).name, eWC)
#     for i in 1:length(rWC)
#         Y = TFanalyticsignal(wd_peri(eWC[i]), 256, 0, 0.5; fmax=30)
#         elgtf[e, rWC[i]] = vec(mean(extract(Y, (50, 59) .* 256 ./ 1010.1, :), dims=1))
#         ehgtf[e, rWC[i]] = vec(mean(extract(Y, (70, 100) .* 256 ./ 1010.1, :), dims=1))
#     end
# end
# elgtf[eventtbl.event[1], "Amygdala"]

# event_low_gamma_tf = GMap(elgtf)
# event_high_gamma_tf = GMap(ehgtf)

# ##==========================================================
# # using JLD2

# # # savedict=Dict([(:WC,WC), (:famtbl,eventtbl), (:wd_peri,wd_peri), (:peri_spec,peri_spec),(:event_region_coherence,event_region_coherence), (:event_mob_respiratory_tf,event_mob_respiratory_tf), (:event_low_gamma_tf,event_low_gamma_tf), (:event_high_gamma_tf, event_high_gamma_tf)])


# # @save "./TempData/bunchastuff.jld2"
# # # CSV.write("./TempData/EventFamiliarity.csv", eventtbl)

# ##==========================================================







