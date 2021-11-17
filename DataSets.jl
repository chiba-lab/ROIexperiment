abstract type AbstractDataSet end

import SplitApplyCombine as SAC

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

function getEventTimeSeries(ds::DataSet, eids; astable=true)
    lfpids=ds.metatable["LFP"].id
    lfp=[ds.data[l] for l in lfpids]
    E=[ds.data[i] for i in eids]
    # return   [filter(x->x ∋ ds.data[e], lfp) for e in events]
    et=[]
    ett=[]
    for e in E
        for l in lfp
            if l.file == e.file
               push!(et, l)
               push!(ett, e)
            end
        end
    end
    if astable
        return leftjoin(DataFrame(id = map(x->x.id, ett), lfp = data.(et), region = map(x->x.region, et)), ds.metatable["Event"], on="id")
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
