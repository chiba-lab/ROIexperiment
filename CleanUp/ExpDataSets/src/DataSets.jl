import SplitApplyCombine as SAC

abstract type AbstractDataSet end

function make_metatable(d)
    g=SAC.group(x->split(string(typeof(x)),'.')[end], values(d))
    return Dict(keys(g).=>map(DataFrame, values(g)))
end


struct DataSet <: AbstractDataSet 
    data
    metatable 
    DataSet(data)=new(data, make_metatable(data))   
end

abstract type AbstractDataSubset <: AbstractDataSetObject end 

function source(s::AbstractDataSubset) 
    return data(s)
end

function getEventTimeSeries(ds::DataSet, eids; astable=true, pre=0, post=0)
    lfpids=ds.metatable["LFP"].id
    lfp=[ds.data[l] for l in lfpids]
    E=[ds.data[i] for i in eids]
    # return   [filter(x->x ∋ ds.data[e], lfp) for e in events]
    etv=[]
    # ett=[]
    for e in E
        for l in lfp
            if e ∈ l
                et=EventTimeSeries(e, l)
                if !isartifact(et)
                    push!(etv, et)
                end
            end
        end
    end
    if astable
        return leftjoin(DataFrame(event_id = map(x->x.event.id, etv), lfp_id= map(x->x.timeseries.id, etv), lfp = data.(etv, Ref(pre), Ref(post)), region = map(x->x.timeseries.region, etv)), ds.metatable["Event"], on=("event_id"=> "id"))
        # return DataFrame(meta.(etv))|>x->@transform(x, lfp=data.(etv))
    else
        return et
    end
end

# abstract type AbstractGroupedDataSet end 

# function groupkeys(g::AbstractGroupedDataSet) end

# function groupvalues(g::AbstractGroupedDataSet) end

# function data(g::AbstractGroupedDataSet) end

# function data(g::AbstractGroupedDataSet, (k,v)) end

# function metasubset(ds::DataSet, tname, args...)
#     return DataFrames.select(ds.metatable[tname], args...)
# end

# function metagroupby(ds::DataSet, tname, args...)
#     return DataFrames.groupby(ds.metatable[tname], args...)
# end

# function metacombine(ds::DataSet, tname, args...)
#     return DataFrames.combine(ds.metatable[tname], args...)
# end


# function events(ds::DataSet)
#     return ds.metatable["Events"]
# end

# function trials(ds::DataSet)
#     return ds.metatable["Trials"]
# end

# function lfp(ds::DataSet)
#     return ds.metatable["LFP"]
# end
