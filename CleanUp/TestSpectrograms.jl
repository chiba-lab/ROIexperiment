using ExpDataSets, Dates, CSV, FileIO, JLD2, Chain, DataFrames, DataFramesMeta, StatsBase, Measurements, CairoMakie, DSP, SignalAnalysis
include("other_utils.jl")

et = load("Data/behavior_dataset_1pre_1post.jld2", "event_time_series");


function mapreduce_ets(fun, ets, tscol, funname, window = (t, s, e) -> t)
    out = @chain ets begin
        transform(tscol, :start_time, :end_time => ByRow.(window) => tscol)
        transform(tscol => ByRow.(x -> normalize(x)) => tscol)
        # @subset((:end_time .- :start_time) .≥ 1.0)
        transform(tscol => ByRow.(x -> normalize(x)) => tscol)
        transform(tscol => ByRow.(fun) => funname)
    end
    # @select(:rat, :trial_type, :agent, :behavior, :region, :pspec)
end


s = spectrogram(values(test.lfp[1]); fs = 1010.1)
import SignalAnalysis as SA

s = SA.spectrogram(values(test.lfp[1]); fs = 1010.1)

using Plots

using Dates










test = mapreduce_ets(values, et, :lfp, :lfp_values)

length(values(test.lfp[1]))

window = (t, s, e) -> t(s, s + 0.001)

window

# function msemreduce(et, col, subs, gcols)
#     resultname="mean_"*String(col)
#     resultvalsname = rstr*"_vals"
#     resultsemname = rstr*"_err"
#     out =  @chain et begin
#     # @subset()
#     # groupby(gcols...)
#     @combine($resultname = [ meanandstderr($col)])
#     @transform($resultname => ByRow.(x -> Measurements.value.(x)) => $resultvals_name, :mps => ByRow.(x -> Measurements.uncertainty.(x)) => $resultsemname)
#     end
# end


# test=msemreduce(test, :lfp, [:region .== "mob"], :trial_type)







oi = @chain begin
    @subset(et, :agent .== "Object", :behavior .== "Immobility")
    select(:lfp, :file, :start_time, :end_time)
    @transform(@byrow :duration = duration(:lfp))
    minimum(_.duration)
end

ps = @chain et begin
    @subset((:end_time .- :start_time) .≥ 1.0)
    transform(:lfp => ByRow.(x -> normalize(x)) => :lfp)
    transform(:lfp => ByRow.(x -> binps(compute_power_spectrum(x), 60, 120)) => :pspec)
    @select(:rat, :trial_type, :agent, :behavior, :region, :pspec)
end

null_region_mps = @chain ps begin
    @subset(:behavior .== "null")
    transform(:pspec => ByRow.(x -> values(x)) => :pspec)
    groupby([:region])
    @combine(:mps = [meanandstderr(:pspec)])
    transform(:mps => ByRow.(x -> Measurements.value.(x)) => :mpsvals, :mps => ByRow.(x -> Measurements.uncertainty.(x)) => :err)
    select(:region, :mpsvals => :mps, :err)
end

behavior_mps = @chain ps begin
    @subset(:behavior .!= "null")
    transform(:pspec => ByRow.(x -> values(x)) => :pspec)
    groupby([:behavior, :region])
    @combine(:mps = [meanandstderr(:pspec)])
    transform(:mps => ByRow.(x -> Measurements.value.(x)) => :mpsvals, :mps => ByRow.(x -> Measurements.uncertainty.(x)) => :err)
    select(:behavior, :region, :mpsvals => :mps, :err)
end

freqs = 2:2:120;
b = "Immobility";
r = "Amyg";
s = @subset(behavior_mps, :behavior .== b, :region .== lowercase(r))
labels = r
specs = s.mps
errs = s.err

gr()
p = plot(
    freqs,
    specs,
    ribbon = errs[1],
    fillalpha = 0.2,
    label = b,
    xaxis = :log,
    yaxis = :log,
    ylabel = "Power Spectral Density",
    xlabel = "Frequency (Hz)",
    title = r,
    linewidth = 2.5,
    grid = false,
)
ns = @subset(null_region_mps, :region .== lowercase(r))
plot!(p, freqs, ns.mps, ribbon = ns.err[1], fillcolor = :black, fillalpha = 0.15, linewidth = 2, xlim = (2, 120), linestyle = :dash, linecolor = :black, label = "Baseline")


(respiratory = [3, 12], theta = [5, 10], beta = [15, 35], gamma_low = [50, 59], gamma_high = [70, 100])





(event_time_series, DSPfun, window, subset, (group)) -> (array(dsp), indicies)