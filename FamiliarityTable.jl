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

using FourierAnalysis


wd_peri = GMap(x -> length(x) <= 2532 ? x[round(Int, length(x) / 2)-256:round(Int, length(x) / 2)+256] : x[1010:end-1010], WC)

peri_spec = GMap(x -> spectra(x, 256, 512), wd_peri.(WC))
# # post_spec = GMap(x -> spectra(x, 256, 512), wd_post.(WC)) 

using DataFrames, DataFramesMeta

te=[inverse(event_trial)(x) for x in collect(codom(event_trial))]

tr = map(te) do x
    if isa(x, BehavioralEvent)
        filter(x->isa(x, Agent), unique(map(x->x.receiver, [x])))|>x->map(x->x.name, x)|>x->isempty(x) ? [missing] : x
    else
         filter(x->isa(x, Agent), unique(map(x->x.receiver, x)))|>x->map(x->x.name, x)|>x->isempty(x) ? [missing] : x
    end
end

td=Dict(zip(collect(codom(event_trial)), tr))




eventtbl=DataFrame()
eventtbl.rat=map(x -> window_lfp(inv(window_data)(x)).rat.name, WC)
eventtbl.session=map(x -> window_lfp(inv(window_data)(x)).session.filename, WC)
eventtbl.date=map(x -> window_lfp(inv(window_data)(x)).session.date, WC)
eventtbl.trial_onset=map(x -> event_trial(window_data_event(x)).start_time, WC)
eventtbl.trial_type = map(x -> event_trial(window_data_event(x)).condition.name, WC)
eventtbl.trial_agents=map(x -> td[event_trial(window_data_event(x))], WC)
eventtbl.onset=map(x -> inv(window_data)(x).onset, WC)
eventtbl.region = map(x -> window_data_region(x).name, WC)
eventtbl.behavior_type = map(x -> window_data_event(x).behavior.name, WC)
eventtbl.agent_type = map(x -> agenttype(event_trial(window_data_event(x)).condition), WC)
eventtbl.agent_name = map(x -> window_data_event(x).receiver|>x-> isa(x, Agent) ? x.name : "missing", WC)
eventtbl.lfp=wd_peri.(WC)
specs=peri_spec.(wd_peri.(WC))
eventtbl.spec=map(x->x.y, specs)
fr = FourierAnalysis.fres(256, 512) .* (1010.1 / 256)
eventtbl.freq = fill(collect(fr:fr:fr*256), (length(WC)))
eventtbl.respiratory = map(x -> mean(x, (3, 12) .* 256 ./ 1010.1), specs)
eventtbl.theta = map(x -> mean(x, (5, 10) .* 256 ./ 1010.1), specs)
eventtbl.beta = map(x -> mean(x, (15, 35) .* 256 ./ 1010.1), specs)
eventtbl.gamma_low = map(x -> mean(x, (50, 59) .* 256 ./ 1010.1), specs)
eventtbl.gamma_high = map(x -> mean(x, (70, 100) .* 256 ./ 1010.1), specs)


sort!(eventtbl, [:date, :trial_onset, :onset])

agids=filter(!ismissing, reduce(union, eventtbl.trial_agents))


ssd=Dict()
ssd[("missing","missing")]=missing
for ai in agids[2:end]
    sbs=@subset(eventtbl, :agent_name .== ai)
    ss=sort(unique(sbs.date))
    for i in 1:length(ss)
        ssd[(ai,ss[i])]=i
    end
end

eventtbl = @transform(eventtbl, :sess_exp_num = :agent_name ∈ agids[2:end] ? ssd[(:agent_name, :date)] : "missing")

ttd=Dict()
ttd[("missing", "missing")]=missing
for ai in agids[2:end]
    sbt=@subset(eventtbl, :agent_name .== ai)
    tt=sort(unique(sbt.trial_onset))
    for i in 1:length(tt)
        ttd[(ai,tt[i])]=i
    end
end

eventtbl = @transform(eventtbl, :trial_exp_num = :agent_name ∈ agids[2:end] ? ttd[(:agent_name, :trial_onset)] : "missing")



eventtbl


using JLD2, CSV


save("./TempData/EventSpecTbl.jld2", "EventSpecs", eventtbl)

CSV.write("./TempData/EventSpecTbl_csv.csv", eventtbl)








