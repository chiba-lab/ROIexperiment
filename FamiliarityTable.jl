include("./ROIontology.jl")


ROIexp = load("./Data/ROIexp.jld2", "ROIexp");
lfp = ROIexp[:lfp_data]
events = ROIexp[:behavioral_events]
trials = ROIexp[:trials]
sessions = ROIexp[:sessions]

include("./Links.jl")

WD = filter(x -> window_data_event(window_data(x)) ∈ dom(event_trial), collect(dom(window_data)))
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





eventtbl = DataFrame()
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
# specs=peri_spec.(wd_peri.(WC))
# eventtbl.spec=map(x->x.y, specs)
# fr = FourierAnalysis.fres(256, 512) .* (1010.1 / 256)
# eventtbl.freq = fill(collect(fr:fr:fr*256), (length(WC)))
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
            return r[rn.==minimum(rn)]
        elseif isequal(n, "Familiar")
            return r[rn.==maximum(rn)]
        else
            return r[1]
        end
    else
        missing
    end

end

eventtbl.receiver_name = [matchreceiver(eventtbl.receiver_attribute[i], eventtbl.trial_agents[i], eventtbl.trial_exp_num[i],eventtbl.behavior_type[i]) for i in 1:size(eventtbl, 1)]


function receiverexpnum(n, r, rn)
    if ismissing(n)
        return missing
    elseif isempty(n)
        return missing
    else
        return rn[r.==n]
    end
end

eventtbl.receiver_exp_num = [receiverexpnum(eventtbl.receiver_name[i], eventtbl.trial_agents[i], eventtbl.trial_exp_num[i]) for i in 1:size(eventtbl, 1)]





using JLD2, CSV


save("./TempData/EventFamiliarity.jld2", "events", eventtbl)

CSV.write("./TempData/EventFamiliarity.csv", eventtbl)








