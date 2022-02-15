#* 1 filter and group data 
#* 2 run DSP
#* 3 averaging 
#* 4 plots
#* 5 Stats
##==============================================================================
include("./ROIontology.jl")
using Associates
import Associates: unwrapsinglet

import Base.(⊆)
(⊆)(i::Tuple{Session, ClosedInterval}, j::Tuple{Session, ClosedInterval}) = i[1]==j[1] ? i[2] ⊆ j[2] : false

##==============================================================================
#Load data

ROIexp = load("./Data/ROIexp.jld2", "ROIexp");
lfp = ROIexp[:lfp_data]
events = ROIexp[:behavioral_events]
trials = ROIexp[:trials]
sessions = ROIexp[:sessions]


##==============================================================================
#Make Links

lfp_session=GMap(x->x.session, lfp)
lfp_region=GMap(x->x.region, lfp)
lfp_rat=GMap(x->x.rat, lfp)


trial_session=GMap(x->x.session, trials)
trial_time=GMap(x->time_interval(x), trials)
trial_coord=GMap(x->(trial_session(x), trial_time(x)), trials)
trial_condition=GMap(x->x.condition, trials)
trial_condition_name=GMap(x->x.condition.name, trials)
trial_condition_agenttype=GMap(x->agenttype(x.condition), trials)


[TrialCondition{Missing}] ∪ [TrialCondition{Rat}] 
event_session=GMap(x->x.session, events)
event_time=GMap(x->time_interval(x), events)
event_coord=GMap(x->(event_session(x), event_time(x)), events)
event_behavior=GMap(x->x.behavior, events)
event_actor=GMap(x->x.actor, events)
event_reciever=GMap(x->x.receiver, events)

event_trial_dict=Dict{BehavioralEvent, Trial}()
for e in events
    t=filter(x->event_coord(e)⊆trial_coord(x), trials)
    if !isempty(t)
        event_trial_dict[e]=t[1]
    end
end
event_trial=GMap(event_trial_dict)
event_lfp=inv(lfp_session)∘event_session
struct EventData
    event
    data
end









##==============================================================================
#Selecting Events
using MacroTools: postwalk, @capture, isexpr 

macro where(ex) 
    expr = _where(ex)
    return expr
end

function _where(ex)
    postwalk(x -> (isexpr(x) && x.head==:tuple) ? quote preimage(($x)[1],($x)[2]) end : x,  ex)
end


##==============================================================================
#Selecting Events
# @where [(propertymap, values) or Set]  set operation [(propertymap, values) or Set] ...
event_subs = @where (trial_condition_name∘event_trial, "Habituation") ∩ (event_behavior, Groom) ∩ dom(event_lfp)



##==============================================================================
#Selecting LFP

lfp_subs = @where image(event_lfp, event_subs) ∩ (lfp_region, AMG)

lfp_subs











