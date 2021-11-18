abstract type AbstractTimeSeries <:AbstractArrayData end

function times(ts::AbstractTimeSeries) end

# function Base.values(ts::AbstractTimeSeries, t...; get_nearest=false)
#     if get_nearest 
#         return  values(ts)[findnearest.(times(ts), Ref(t))]
#     else
#         return  [values(ts)[times(ts).==i] for i in t]
#     end
# end


function axes(ts::AbstractTimeSeries)
    return (time = times(ts))
end

function start_time(ts::AbstractTimeSeries) end

function end_time(ts::AbstractTimeSeries) end

function duration(ts::AbstractTimeSeries) 
    return end_time(ts) - start_time(ts)
end

function nearesttime(ts::AbstractTimeSeries, t)
    return nearest(times(ts), t)
end

# function getinterval(ts::AbstractTimeSeries, t₁, t₂) end



struct TimeSeries <: AbstractTimeSeries 
    values 
    times 
end

function Base.values(ts::TimeSeries) 
    return ts.values 
end

function times(ts::TimeSeries) 
    return ts.times
end

function start_time(ts::TimeSeries) 
    return ts.times[1]
end

function end_time(ts::TimeSeries) 
    return ts.times[end]
end 

function getintervalidx(ts::TimeSeries, t₁, t₂) 
    return findall(t₁.≤ ts.times.≤t₂)
end

function (ts::TimeSeries)(axs)
    TimeSeries(values(ts, axs[:times]), axs[:times])
end

function (ts::TimeSeries)(t₁::Number, t₂::Number)
    idx = getintervalidx(ts, t₁, t₂)
    ts(times=idx)
end

function (ts::TimeSeries)(newdata)
    TimeSeries(newdata, ts.times)
end

function (ts::TimeSeries)(newdata::Array, newtimes::Vector)
    TimeSeries(newdata, newtimes)
end




struct ContinuousTimeSeries <: AbstractTimeSeries 
    values 
    fs
    start_time
    end_time
end

function Base.values(ts::ContinuousTimeSeries) 
    return ts.values
end

function times(ts::ContinuousTimeSeries) 
    return collect(ts.start_time:(1.0/(ts.fs)):ts.end_time)
end

function start_time(ts::ContinuousTimeSeries) 
    return ts.start_time
end

function end_time(ts::ContinuousTimeSeries) 
    return ts.end_time
end

function getintervalidx(ts::ContinuousTimeSeries, t₁, t₂) 
    int_start =  findnearest(times(ts), t₁)
    int_end = findnearest(times(ts), t₂)
    return int_start, int_end
end

# function (ts::ContinuousTimeSeries)(axs...)
#     ContinuousTimeSeries(values(ts)[axs[:times]], ts.fs, axs[:times][1][1], axs[:times][1][end])
# end

function (ts::ContinuousTimeSeries)(t₁::Number, t₂::Number)
    ids, ide = getintervalidx(ts, t₁, t₂)
    ContinuousTimeSeries(values(ts)[ids:ide], ts.fs, times(ts)[ids], times(ts)[ide])
end

function (ts::ContinuousTimeSeries)(newdata)
    ContinuousTimeSeries(newdata, ts.fs, ts.start_time, ts.end_time)
end

function (ts::ContinuousTimeSeries)(newdata::Array, newtimes)
    ContinuousTimeSeries(newdata, ts.fs, newtimes[1], newtimes[end])
end

function normalize(ts::ContinuousTimeSeries)
    ContinuousTimeSeries(StatsBase.normalize(values(ts)), ts.fs, ts.start_time, ts.end_time)
end

