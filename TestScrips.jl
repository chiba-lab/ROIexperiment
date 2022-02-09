include("./ROIontology.jl")
using Chain

ROIexp = load("./Data/ROIexp.jld2", "ROIexp")
lfp = ROIexp[:lfp_data]
events = ROIexp[:behavioral_events]


function get_event_lfp_recording(event, lfp)
    f1(l) = in_interval(event, l)
    el = filter(x -> f1(x), lfp) 
    !isempty(el) ? el : missing
end

el = get_event_lfp_recording(events[500], lfp)

Tree(el)



events[2]
