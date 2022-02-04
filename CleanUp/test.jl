using ExpDataSets, Dates, CSV, FileIO, JLD2, Chain, DataFrames, DataFramesMeta, StatsBase, Measurements, Plots 

ds = load("Data/freeroam_dataset.jld2", "freeroam_dataset");

# E=@subset(ds.metatable["Event"], :behavior .== "null")

et = getEventTimeSeries(ds, ds.metatable["Event"].id)


save("Data/freeroam_dataset.jld2", "freeroam_dataset", ds, "event_time_series", et)



function meanandstderr(v)
    measurement.(mean(v), (StatsBase.std(v))/sqrt(length(v)))
end

function meanps(df, by)
    ps=transform(df,  :lfp => ByRow.(x->values(binps(compute_power_spectrum(x))))=>:pspec)
    g=groupby(ps, by)
    @combine(g, mps=[meanandstderr(:pspec)])
end

function plotps(mps, args...)
    plot(Measurements.value.(mps), ribbon=Measurements.uncertainty.(mps), fillalpha=.25, args...)
end

function plotps!(mps, args...)
    plot!(Measurements.value.(mps), ribbon=Measurements.uncertainty.(mps), fillalpha=.25, args...)
end



