# include("./ROIontology.jl")


# ROIexp = load("./Data/ROIexp_FR.jld2", "ROIexp");
# lfp = ROIexp[:lfp_data]
# events = ROIexp[:behavioral_events]
# trials = ROIexp[:trials]
# sessions = ROIexp[:sessions]

lfp_session=GMap(x->x.session, lfp)
lfp_region=GMap(x->x.region, lfp)
lfp_rat=GMap(x->x.rat, lfp)


trial_session=GMap(x->x.session, trials)
trial_time=GMap(x->time_interval(x), trials)
trial_coord=GMap(x->(trial_session(x), trial_time(x)), trials)
trial_condition=GMap(x->x.condition, trials)
trial_condition_name=GMap(x->x.condition.name, trials)
trial_condition_agenttype=GMap(x->agenttype(x.condition), trials)

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
      if  (j.onset..j.offset ⊆ start_time(j.data)..end_time(j.data))
            lfp_window_event_dict[j]=i
      end
    end
end
lfp_window_event_dict
lfp_window_event=GMap(lfp_window_event_dict)
event_lfp_window=inv(lfp_window_event)

window_lfp=GMap(x->x.data,collect(dom(lfp_window_event)))

win=collect(dom(lfp_window_event))
d=get_data.(win,1.0,1.0)
# d=get_int_post_start.(win,1.0)
idx=Base.Iterators.findall(x->!any(isnan,x), d)
window_data=GMap(Dict(zip(win[idx], d[idx])), Dict(zip(d[idx], win[idx])))
win=collect(dom(window_data))
function isartifact(x::DataWindow)
    if x.data.region==MOB
        thresh=.6
    else
        thresh=.4
    end
    # if any(window_data(x).≥thresh)
    if sum((window_data(x) .≥ thresh)) > (length(window_data(x)) * .01)
        return true
    else
        return false
    end
end
win=filter(!isartifact, win)




window_data=GMap(Dict(zip(win, window_data.(win))), Dict(zip(window_data.(win), win)))
window_data_event=GMap(x->inverse(event_lfp_window)(inverse(window_data)(x)), collect(codom(window_data)))
window_data_region=GMap(x->inverse(window_data)(x).data.region, collect(codom(window_data)))
