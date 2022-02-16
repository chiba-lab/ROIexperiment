include("./ROIontology.jl")
using Mappings

using Chain

ROIexp = load("./Data/ROIexp.jld2", "ROIexp")
lfp = ROIexp[:lfp_data]
events = ROIexp[:behavioral_events]
trials = ROIexp[:trials]
session = ROIexp[:sessions]






sesh = AssociativeMap(Dict(zip(events, map(x -> x.session, events))))
dom(sesh)

lfpsesh = AssociativeMap(Dict(zip(lfp, map(x -> x.session, lfp))))

seshin = inverse(sesh)

dom(seshin) ∩ codom(lfpsesh) |> x -> (preimage(lfpsesh, x), image(seshin, x))
# inv(inv(sesh) ∘ lfpsesh)

# \



# setdiff(dom(seshin), codom(lfpsesh))

# inv(sesh)(session[1])



function lfp_recording(event::BehavioralEvent, lfp)
    function f1(l)
        in_interval(event, l)
    end
    el = filter(x -> f1(x), lfp)
    !isempty(el) ? el : missing
end



function get_event_trial(event, trials)
    f1(t) = in_interval(event, t)
    et = filter(x -> f1(x), trials)
    !isempty(et) ? et : missing
end

function removemissings(x::Array{T,N}) where {T,N}
    filter(!ismissing, x) |> x -> convert(Array{typeof(x[1]),N}, x) |> x -> filter(!isempty, x)
end

function unelvec(x)
    map(x -> isempty(x) ? x : x, x)
end

e = get_event_trial.(events, Ref(trials)) |> removemissings |> unelvec
count(!ismissing(x) for x in lfp_recording.(events, Ref(lfp)))

EventLFPTable = Table(events = events, trials = get_event_trial.(events, Ref(trials)) |> unelvec, lfp = get_event_lfp_recording.(events, Ref(lfp)))

@chain EventLFPTable begin
    filter(x -> isa(x.trials, Trial), _)
    filter(x -> isa(x.trials, Trial), _)
    # filter(x->x.trials.condition == ITR(Robot), _) 
    filter(x -> x.events.behavior.name == "Immobility", _)
end




amyglfp = filter(x -> x.region.name == "Amygdala", lfp)
a = @chain events begin
    filter(x -> x.behavior.name == "Immobility" % %, _)
    get_event_lfp_recording.(_, Ref(amyglfp))
    filter(!ismissing, _)
    map(x -> x[1], _)
    convert(Vector{LFPRecording}, _)
end



