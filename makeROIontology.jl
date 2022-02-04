include("ROIontology.jl")

files = glob("*.csv", "Data/CSV/lfp_data")
tbls = [Table(CSV.File(f)) for f in files]
r2r = Dict([(AMG, :amyg), (MOB, :mob), (CA2, :ca2), (INS, :insula)])
t2lfp(t, r) = LFPRecording(Session(t.filename[1], Date(t.date[1]), 0, t.time[end]), t.time[1], t.time[end], Rat("EG7"), r, 1010.1, collect(getproperty(t, r2r[r])))
lfpdata = map(x -> t2lfp.(tbls, Ref(x)), [AMG, MOB, CA2, INS])
lfpdata=vcat(lfpdata...)
sessions = map(x -> x.session, lfpdata)

condD = Dict([("EE", EE), ("Empty", EE), ("Freeroam", FR), ("Habituation", HBT), ("OF", OF), ("Interaction", ITR), ("Object", Object), ("Rat", Rat), ("Rat ", Rat), ("Robot", Robot), (missing, missing)])
macro name2cond(C, T)
    :($C($T))
end
function sesfun(i)
    f(j) = x -> x.filename == j
    filter(f(i), sessions)
end
tt = Table(CSV.File("Data/CSV/trial_info/behavioral_trials.csv"; header = false))
trials = map(x -> Trial(sesfun(x.Column2), x.Column3, x.Column4, (condD[x.Column5])(condD[x.Column6])), tt)

function str2sec(str)
    a = split(str, ":")
    b = parse.(Ref(Float64), a)
    b[1]u"minute" + b[2]u"s"
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
        "Empty Enclosure"
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
b1 = map(x -> BehavioralEvent(sesfun(x.VideoName), str2sec(x.StartTime), str2sec(x.EndTime), Behavior(x.EventType), Rat(x.RatName), Rat(x.RatName)), efilt_1)
b2 = map(x -> BehavioralEvent(sesfun(x.VideoName), str2sec(x.StartTime), str2sec(x.EndTime), x.EventType, Rat(x.RatName), getreceiver(x)), efilt_2)

behaviors= vcat(b1, b2)

ROIexp = Dict([(:sessions, StructArray(sessions)), (:trials, StructArray(trials)), (:behaviors, behaviors), (:lfp_recordings,StructArray(lfpdata))])

save("./Data/ROIexp.jld2", "ROIexp", ROIexp)




