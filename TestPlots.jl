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

event_lfp_window=GMap(e->map(x->DataWindow(e.start_time, e.end_time, x), event_lfp(e)), collect(dom(event_lfp)))
lfp_window_event_dict=Dict{DataWindow, BehavioralEvent}()
for i in collect(dom(event_lfp_window))
    for j in event_lfp_window(i)
      if  j.onset..j.offset ⊆ start_time(j.data)..end_time(j.data) 
            lfp_window_event_dict[j]=i
      end
    end
end
lfp_window_event=GMap(lfp_window_event_dict)
event_lfp_window=inv(lfp_window_event)

window_data=GMap(x->x.data,collect(dom(lfp_window_event)))



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
event_subs = @where (trial_condition_name∘event_trial, "Habituation") ∩ (event_behavior, Groom) ∩ dom(event_lfp_window)



##==============================================================================
#Selecting LFP

event_subs_lfp_windows = @where (lfp_window_event, event_subs) ∩ (lfp_region∘window_data, AMG) ∩ (lfp_rat∘window_data, RRSD18)


##==============================================================================
#Extract Data
#get_data(data_winndow, pre, post)
event_subs_data = GMap(x->get_data(x,1,1), event_subs_lfp_windows)

##==============================================================================

#run analysis
event_subs_spectra = GMap(x->FA.spectra(event_subs_data(x), 256, 512), event_subs_lfp_windows)

event_subs_power = GMap(x->FA.extract(event_subs_spectra(x),:), event_subs_lfp_windows)


##==============================================================================
#average
using Statistics
mat=hcat(event_subs_power.(event_subs_lfp_windows)...)
m=mean(mat, dims=2)
s=std(mat, dims=2)/sqrt(size(mat,2))



f = Figure();
ax=Axis(f[1, 1], yscale = log10,
        yminorticksvisible = true, yminorgridvisible = true,
        yminorticks = IntervalsBetween(8))
errorbands!(ax, 1:1:256, vec(m), vec(s); linewidth=2)
f


f = Figure();
ax=Axis(f[1, 1])
t=collect(1:10)
y=sin.(t)
eventwindow!(ax, 1:10, sin, 3, 7)
f