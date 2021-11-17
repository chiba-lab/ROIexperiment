module ROI

   import Base.values
   import Base.:∈
   import Base.:∋

   import Lazy.@forward

    using PyCall,  Dates, DataFrames, DataFramesMeta, StatsBase

    const filt = PyNULL()
    const spc = PyNULL()

    function __init__()
        copy!(filt, pyimport_conda("neurodsp.filt", "neurodsp"))
        copy!(spc, pyimport_conda("neurodsp.spectral", "neurodsp"))
    end

    export AbstractArrayData,  DataSet, amp, filter_signal , id, source, AbstractDataInterval, Event, axes, findnearest,  spc, AbstractDataSet, LFP, bandpower,freqs, meta, start_time, AbstractDataSetObject, NeuroPowerSpec ,compute_power_spectrum,  fstr2date, minstr2seconds,  times, AbstractDataSubset, PowerSpectrum, data, gendid, nearest, values, AbstractDataType, Session,  duration, genid, nearesttime, AbstractGroupedDataSet, Spectrum, end_time, getband,  phase, AbstractSpectrum, TimeSeries,  getintervalidx , power, AbstractTimeSeries, Trial, fbands, groupkeys, ContinuousTimeSeries, filt, groupvalues, ∈,  ∋, EventTimeSeries,  getEventTimeSeries, binps


    include("utils.jl")
    include("DataSetObjects.jl")
    include("DataSets.jl")
    include("DataTypes.jl")
    include("TimeSeries.jl")
    include("Spectra.jl")
    include("neurodsp.jl")
    include("DataIntervals.jl")
    include("ROIDataObjs.jl")

end
