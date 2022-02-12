#* 1 filter and group data 
#* 2 run DSP
#* 3 averaging 
#* 4 plots
#* 5 Stats
##==============================================================================
include("./ROIontology.jl")
using Associates
import Base.hash


##==============================================================================
# data  handeling 
#* link objs
#* load get_evet_lfp
#* filter 
#* select
#* remove artifacts
#* group

ROIexp = load("./Data/ROIexp.jld2", "ROIexp");
lfp = ROIexp[:lfp_data]
events = ROIexp[:behavioral_events]
trials = ROIexp[:trials]
session = ROIexp[:sessions]


lfp_hashmap = GMap(hash, lfp)
dom(lfp_hashmap)


# for (fun, field) = zip(Symbol.("get_".*String.(fnames)), fieldnames(LFPRecording))
#     eval(quote
#         $fun(k::UInt) = getfield($(inverse(lfp_hashmap))(k), $field)
#     end)
# end

events[1]



lfpsessions=GMap(x->x.session, lfp)




event_behavior=GMap(e->e.behavior, events)

dom(event_behavior)
inverse(event_behavior)



event_agent=GMap(e->e.actor, events)

for k in keys(inverse(event_agent).forward) 
   println(length(inverse(event_agent).forward[k]))
end

lfpsessino=GMap(x->x.session, lfp)

using FourierAnalysis
GMap(x->spectra(x.lfp, 256, 512), lfp)


session_event= inverse(GMap(x->x.session, events))

collect(session_event.forward)[1]

"(lfp<-session)(sessino<-event)(event<-rat)" 


Spectra








te=AssociativeMap(e -> filter(t -> in_interval(e, t), trials), events)
image(te, events[1])

links=Dict()
using ProfileView
@profview elfp = AssociativeMap(e -> filter(l -> in_interval(e, l), lfp), events)
getlfp(e::BehavioralEvent) = links["event->lfp"](e)

@btime el =filter(l -> in_interval(events[1500], l), lfp)
elfp = ans
[f for f in elfp(events[1500])]



links["lfp->rat"] = AssociativeMap(l -> l.rat, lfp)
links["lfp->region"] = AssociativeMap(l -> l.region, lfp)
rat(l::LFPRecording)=links["lfp->rat"](l)
region(l::LFPRecording)=links["lfp->region"](l)

links["event->lfp"].inv_amap
l2s = map(l -> l.session, lfp)
e2s=map(e -> e.session, events)

E=AssociativeMap(Dict(hash.(events).=>e2s))
L=AssociativeMap(Dict(hash.(lfp).=>l2s))

E ∘ L





links=Dict()
using JET
@btime AssociativeMap(D)
@ProfView links["event->behavior"] = AssociativeMap(e -> e.behavior, events)
links["event->conditon"] = AssociativeMap(t -> t.condition, trials) 
links["event->agent"] = AssociativeMap(e -> e.session, events);


inverse(links["event->behavior"])

links["event->conditon"](events[300])

links["trial->session"] = AssociativeMap(t -> filter(s -> in_interval(t, s), session), trials);
links["lfp->region"] = AssociativeMap(l -> l.region, lfp);
links["condition->agent"] = AssociativeMap(agenttype, [x.condition for x in trials]);
links["trial->condition"] = AssociativeMap(t -> t.condition, trials);
links["lfp-rat"] = AssociativeMap(l -> l.rat, lfp)



links["trial-condition"]
links["event-trial"]

using JLD2
count(!isempty(links["event-lfp"](e)) for e in events)
save("./Data/ROIexp.jld2", "ROIexp", ROIexp, "links", links)


trial_type
agent
condition
rat
region




