struct LFP <: AbstractDataSetObject
    id
    file
    date
    rat 
    region 
    data::ContinuousTimeSeries
end

function LFP( file, date, rat ,region, data)
    LFP(abs(rand(Int32)), file, date, rat ,region, data)
end 

function id(l::LFP)
    return(l.id)
end

function data(l::LFP)
    return(l.data)
end
=>
function source(l::LFP)
    return file
end

function meta(l::LFP)
    return ( file=l.file, date=l.date, rat=l.rat, region=l.region, start_time=start_time(l.data), end_time=end_time(l.data))
end

function start_time(l::LFP)
    return(start_time(l.data))
end

function end_time(l::LFP)
    return(end_time(l.data))
end





# function getinterval(d::LFP, i::AbstractDataInterval, pre=0, post=0) 
#     if !∈(d,i)
#         return missing 
#     else
#         return LFP(d.file, d.date, d.rat, d.region, getinterval(d.data, i.start_time - pre, i.end_time + post))
#     end 
# end


struct NeuroPowerSpec <: AbstractDataSetObject
    id 
    timeseries
    data::PowerSpectrum
end

function NeuroPowerSpec(l::LFP)
    return NeuroPowerSpec(abs(rand(Int32)), data(l), compute_power_spectrum(data(l)))
end

function id(n::NeuroPowerSpec)
    return(n.id)
end

function data(n::NeuroPowerSpec)
    return(n.data)
end

function source(n::NeuroPowerSpec)
    return(n.timeseries)
end



struct Session <: AbstractDataInterval
    id 
    file 
    date 
    rat
    start_time
    end_time 
end


function Session(file, date, rat, start_time, end_time)
    return Session(abs(rand(Int32)), file, date, rat, start_time, end_time)
end

function data(s::Session)
    return nothing 
end

function source(s::Session)
    return nothing
end

function meta(s::Session)
    return (file=s.file, date=s.date, rat=s.rat, start_time=s.start_time, end_time=s.end_time)
end

function start_time(s::Session)
    return(s.start_time)
end

function end_time(s::Session)
    return(s.end_time)
end
    




struct Trial <: AbstractDataInterval
    id 
    file 
    date 
    rat
    trial_type
    agent 
    start_time
    end_time
end

function Trial(file, date, rat, trial_type, agent, start_time, end_time)
    return Trial(abs(rand(Int32)), file, date, rat, trial_type, agent, start_time, end_time)
end

function id(t::Trial)
    return(t.id)
end

function data(t::Trial)
    return nothing 
end

function source(t::Trial)
    return nothing
end

function meta(t::Trial)
    return (file=t.file, date=t.date, rat=t.rat, trial_type=t.trial_type, agent=t.agent, start_time=t.start_time, end_time=t.end_time)
end

function start_time(t::Trial)
    return(t.start_time)
end

function end_time(t::Trial)
    return(t.end_time)
end




struct Event<: AbstractDataInterval
    id 
    file 
    date 
    rat
    trial_type
    agent 
    behavior 
    start_time
    end_time
end

function Event(file, date, rat,  trial_type, agent, behavior, start_time, end_time)
    return Event(abs(rand(Int32)), file, date, rat,  trial_type, agent, behavior, start_time, end_time)
end

function id(e::Event)
    return(e.id)
end

function data(e::Event)
    return nothing 
end    

function source(e::Event)
    return nothing
end

function meta(e::Event)
    return (file=e.file, date=e.date, rat=e.rat, trial_type=e.trial_type, agent=e.agent, behavior=e.behavior, start_time=e.start_time, end_time=e.end_time)
end

function start_time(e::Event)
    return(e.start_time)
end

function end_time(e::Event)
    return(e.end_time)
end

struct EventTimeSeries <: AbstractDataSetObject
    id
    event
    timeseries
end

function EventTimeSeries(e, ts)
    return EventTimeSeries(abs(rand(Int32)), e, ts)
end

function id(et::EventTimeSeries)
    return(et.id)
end

function data(et::EventTimeSeries)
    return(et.timeseries.data(start_time(et.event), end_time(et.event)))
end

function data(et::EventTimeSeries, pre, post)
    return(et.timeseries.data(start_time(et.event)-pre, end_time(et.event)+post))
end

function source(et::EventTimeSeries)
    return nothing 
end

function meta(et::EventTimeSeries)
    e=et.event
    return (file=e.file, date=e.date, rat=e.rat, trial_type=e.trial_type, agent=e.agent, behavior=e.behavior, start_time=start_time(e), end_time=end_time(e))
end

function isartifact(et::EventTimeSeries)
    if et.timeseries.region=="mob" 
        thresh=.75
    else
        thresh=.6
    end
    if any(values(data(et)).≥thresh)
        return true
    else
        return false
    end
end



