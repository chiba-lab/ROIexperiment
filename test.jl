include("ROI.jl")
using .ROI, Dates, CSV, FileIO, JLD2, Chain, DataFrames, DataFramesMeta

ds = load("Data/freeroam_dataset.jld2", "freeroam_dataset");

E=@subset(ds.metatable["Event"], :behavior .== "Grooming")

et = getEventTimeSeries(ds, E.id)

@transform(et, @byrow lfp = binps(compute_power_spectrum(:lfp)))
