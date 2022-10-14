
##==============================================================
rlist = ["MOB", "Ca2"]
blist = ["respiratory", "theta", "low_gamma", "high-gamma"]
bdict = Dict(zip(blist, [(2, 6), (6, 12), (50, 40), (70, 100)]))

cmctbl = DataFrame()
cvec = collect(keys(codict))
cmat = hcat(map(x -> [x...], cvec)...)

cmctbl.event = vec(cmat[1, :])
cmctbl.trial_type = map(e -> event_trial(e).condition.name, cmctbl.event)
cmctbl.behavior_type = map(e -> e.behavior.name, cmctbl.event)
agmap = GMap(Dict([(Object, "Object"), (Rat, "Rat"), (Robot, "Robot"), (Missing, "missing")]))
cmctbl.agent_type = agmap.(map(e -> agenttype(event_trial(e).condition), cmctbl.event))
cmctbl.trial_exp_num = map(e -> event_trial_exp_num(e), cmctbl.event)
cmctbl.receiver_exp_num = map(e -> event_receiver_exp_num(e), cmctbl.event)
cmctbl.trial_fam = map(x -> isa(x, Int) ? (x <= 2 ? "Novel" : "Familiar") : "missing", cmctbl.trial_exp_num)
cmctbl.receiver_fam = map(x -> isa(x, Int) ? (x <= 2 ? "Novel" : "Familiar") : "missing", cmctbl.receiver_exp_num)
# cmctbl.fam = map(x -> x[1] == "missing" ? (x[2] == "missing" ? "missing" : x[2]) : x[1], [(cmctbl.trial_fam[i], cmctbl.receiver_fam[i]) for i in 1:size(cmctbl, 1)])

cmctbl.region_1 = vec(cmat[2, :])
cmctbl.region_2 = vec(cmat[3, :])
cmctbl.co_vec = map(x -> codict[x], cvec)
cmctbl.freq = fill(fl, length(cmctbl.co_vec))
cmctbl.respiratory = map(x -> mean(x[2 .<= fl .< 6]), cmctbl.co_vec)
cmctbl.theta = map(x -> mean(x[2 .<= fl .< 6]), cmctbl.co_vec)
cmctbl.low_gamma = map(x -> mean(x[40 .<= fl .< 50]), cmctbl.co_vec)
cmctbl.high_gamma = map(x -> mean(x[70 .<= fl .< 100]), cmctbl.co_vec)


cmctbl_reduc=cmctbl[:, Not([:event, :co_vec, :freq, :trial_exp_num, :receiver_exp_num])]
@subset!(cmctbl_reduc, :region_1 .∈ Ref(("MOB", "Ca2")), :region_2 .∈ Ref(("MOB", "Ca2")))
##==============================================================
using CSV, FileIO, Missings
CSV.write("ca-mob_coherence.csv", cmctbl_reduc)