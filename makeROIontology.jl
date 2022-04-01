
include("./ROIontology.jl")
using Glob
using DSP

files = glob("*.csv", "Data/CSV/lfp_data")
tbls = [Table(CSV.File(f)) for f in files]
r2r = Dict([(AMG, :amyg), (MOB, :mob), (CA2, :ca2), (INS, :insula)])
flt = digitalfilter(Bandstop(59,61; fs=1010.1), Butterworth(4))
# t2lfp(t, r) = LFPRecording(Session(t.filename[1], Date(t.date[1]), 0, t.time[end] - t.time[1]), 0, length(t.time) / 1010.1, Rat(t.rat[1]), r, 1010.1, collect(getproperty(t, r2r[r])))
t2lfp(t, r) = LFPRecording(Session(t.filename[1], Date(t.date[1]), 0, t.time[end] - t.time[1]), 0, t.time[end] - t.time[1], Rat(t.rat[1]), r, 1010.1, filtfilt(flt, collect(getproperty(t, r2r[r]))))

lfpdata = map(x -> t2lfp.(tbls, Ref(x)), [AMG, MOB, CA2, INS])
lfpdata = vcat(lfpdata...)
sessions = map(x -> x.session, lfpdata)
unique!(sessions)

# condD = Dict([("EE", EE), ("Empty", EE), ("Freeroam", FR), ("Habituation", HBT), ("OF", OF), ("Interaction", ITR), ("Object", Object), ("Rat", Rat), ("Rat ", Rat), ("Robot", Robot), (missing, Missing)])
 condD = Dict([("EE", EE), ("Empty",EMPTY), ("Freeroam", FR), ("Habituation", HBT), ("OF", OF), ("Interaction", ITR), ("Object", Object), ("Rat", Rat), ("Rat ", Rat), ("Robot", Robot), (missing, Missing)])

macro name2cond(C, T)
    :($C($T))
end
function sesfun(i)
    strs(x) = split(x, "_") |> x -> [x[1], parse.(Int, x[2:end])]
    si = strs(i)
    f(j) = x -> strs(x.filename) == j
    s = filter(f(si), sessions)
    if isempty(s)
        return Session(i, missing, missing, missing)
    else
        return s[1]
    end
end

# tt = Table(vcat(CSV.File.("Data/CSV/trial_info/behavioral_trials.csv"; header = false), CSV.File.("Data/CSV/trial_info/freeroam_behavioral_trials.csv"; header = false)))
tt= Table(CSV.File.("Data/CSV/trial_info/behavioral_trials.csv"; header = false))
ta = Table(CSV.File.("Data/CSV/trial_info/behavioral_trial_agents.csv"; header = true))
trials = map(x -> Trial(sesfun(x.Column2), round(x.Column3, digits = 1), round(x.Column4, digits = 1), (condD[x.Column5])(condD[x.Column6])), tt)


trials
function str2sec(str)
    a = split(str, ":")
    b = parse.(Ref(Float64), a)
    b[1] * 60 + b[2]
end
function str2date(str)
    a = split(str, "_")
    b = parse.(Ref(Float64), a[2:4])
    Date(b[3] + 2000, b[1], b[2])
end
function matchnovelty(x)
    if ismissing(x.Category3)
        return missing
    end
    if isequal(x.Category3, x.TrialType3)
        return x.AgentName
    elseif isequal(x.Category3, x.TrialType4)
        return x.AgentName2
    else
        return x.AgentName
    end
end
function getreceiver(x)
    y = x.Category3
    if isequal(y, "Empty")
        "Empty"
    elseif isequal(y, "Object")
        Object(matchnovelty(x))
    elseif isequal(y, "Robot")
        Robot(matchnovelty(x))
    elseif isequal(y, "Rat")
        Rat(matchnovelty(x))
    else
        missing
    end
end

et = Table(CSV.File("Data/CSV/event_info/behavioral_events_labeled_header.csv"))
efilt = et[et.EventType.∈Ref(["Groom", "Rear", "Baseline", "Immobility", "Sniff", "Top", "Approach", "Retreat"])]
efilt_1 = efilt[efilt.EventType.∈Ref(["Groom", "Rear", "Baseline", "Immobility"])]
efilt_2 = efilt[efilt.EventType.∉Ref(["Groom", "Rear", "Baseline", "Immobility"])]
b1 = map(x -> BehavioralEvent(sesfun(x.VideoName), round(str2sec(x.StartTime),digits=1), round(str2sec(x.EndTime), digits=1),Behavior(x.EventType), Rat(x.RatName), missing), efilt_1)
b2 = map(x -> BehavioralEvent(sesfun(x.VideoName), round(str2sec(x.StartTime),digits=1), round(str2sec(x.EndTime), digits=1),Behavior(x.EventType), Rat(x.RatName), getreceiver(x)), efilt_2)

# et = Table(CSV.File("Data/CSV/event_info/freeroam_behavioral_events_labeled_header.csv"))
# efilt = et[et.EventType.∈Ref(["Grooming", "Rearing", "Immobility"])]
# # efilt_1 = efilt[efilt.EventType.∈Ref(["Groom", "Rear", "Baseline", "Immobility"])]
# # efilt_2 = efilt[efilt.EventType.∉Ref(["Groom", "Rear", "Baseline", "Immobility"])]
# b3 = map(x -> BehavioralEvent(sesfun(x.VideoName), round(x.StartTime, digits = 1), round(x.EndTime, digits = 1), Behavior(x.EventType), Rat(x.RatName), missing, FR(condD[x.Category1])), efilt)
# # b4 = map(x -> BehavioralEvent(sesfun(x.VideoName), x.StartTime, x.EndTime, Behavior(x.EventType), Rat(x.RatName), getreceiver(x)), efilt_2)


behavioral_events = vcat(b1, b2)
# behavioral_events = b3
# # sessions = StructArray(sessions)
# # trials = StructArray(trials)
# # lfpdata = StructArray(lfpdata)
# # behavioral_events = StructArray(behavioral_events)

ROIexp = Dict([(:sessions, sessions), (:trials, trials), (:behavioral_events, behavioral_events), (:lfp_data, lfpdata)])


save("./Data/ROIexp.jld2", "ROIexp", ROIexp)




