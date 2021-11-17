# Make Data Set From CSV's 

# include("utils.jl")
# include("EventDataSets.jl")
# include("DataTypes.jl")
# include("neurodsp.jl")
include("ROI.jl")

using .ROI, DataFrames, Dates, CSV, FileIO, Glob, JLD2, DataFramesMeta, Chain, Random, StatsBase

#= ################################
 !Make Trials Info Table  
################################### =#
# fnames=["Data/CSV/freeroam_behavioral_trials.csv","Data/CSV/freeroam_approach_retreat_trials.csv"]
files=glob("freeroam*", "Data/CSV/trial_info") 
df = DataFrame.(CSV.File.(files;header=[:rat, :file_id, :start_time, :end_time, :trial_type,  :agent_type,:agent_1_novelty,:agent_2_novelty], missingstring=""))|>x->reduce(vcat, x)
trial_info = @chain df begin 
    @subset(:rat .∈ Ref(["EG7" "RRSD18" "RRSD28" "RRSD17"]))
    transform(:file_id => ByRow(fstr2date)=> :date)
    transform([:file_id, :start_time, :end_time] =>ByRow((x,y,z)->genid(x,y,z))=>:trial_id)
    select(:trial_id, :file_id, :rat, :date,:start_time, :end_time, :trial_type, :agent_type)
end


# #= ################################
#  Make Event Info Table 
# ################################### =#
# fnames=["Data/CSV/freeroam_behavioral_events.csv","Data/CSV/freeroam_approach_retreat_events.csv"]
files=glob("freeroam*", "Data/CSV/event_info") 
df = DataFrame.(CSV.File.(files; header=[:rat, :file_id, :start_time, :end_time, :behavior_type, :subtype_1, :subtype_2, :subtype_3, :trial_type, :agent_type,:agent_1_novelty, :agent_2_novelty, :agent_1, :agent_2], missingstring=""))|>x->reduce(vcat, x)
event_info = @chain df begin 
    transform(:file_id =>ByRow(fstr2date)=>:date)
    transform([:file_id, :start_time, :end_time] =>ByRow((x,y,z)->genid(x,y,z))=>:event_id)
    @select(:event_id, :file_id, :rat, :date, :start_time, :end_time,:behavior_type,:trial_type, :agent_type)
    DataFrames.rename(:start_time => :onset, :end_time => :offset)
    unique
end

event_info_null = event_info[sample(1:size(event_info,1),500),:]
event_info_null = @eachrow event_info_null begin
    shift=randn()*rand([1,-1])*10+20
    :onset = :onset + shift
    :offset = :offset + shift
    :event_id = genid(:file_id,:onset,:offset)
    :behavior_type = "null"
end

event_info=vcat(event_info, event_info_null)


#= ################################
 Make LFP Data Table 
################################### =#
files=glob("*.csv", "Data/CSV/lfp_data") 
files=files ∩ ("Data/CSV/lfp_data/".*event_info.file_id.*".csv")
df = CSV.File(files;dateformat="yyyy-mm-dd")|>DataFrame
df
lfp_data = @chain df begin 
    # select(Not(:insula))
    groupby([:rat, :date, :filename])
    # @combine(@astable :lfp = ContinuousTimeSeries(hcat(:amyg, :mob, :ca2), 1010.1, ["amyg","mob","ca2"]))
    @combine($AsTable=(amyg = ContinuousTimeSeries(:amyg, 1010.1, first(:time), last(:time)), mob = ContinuousTimeSeries(:mob, 1010.1, first(:time), last(:time)), ca2 = ContinuousTimeSeries(:ca2, 1010.1, first(:time), last(:time))))
    rename(Dict(:filename => :file_id))
    # # transform(:file_id => ByRow.(f -> f[begin:end - 4])=> :file_id)
    @transform(@byrow :amyg = filter_signal(:amyg, "bandstop", (59, 61), "fir", 60, nothing, false))
    @transform(@byrow :mob = filter_signal(:mob, "bandstop", (59, 61), "fir", 60, nothing, false))
    @transform(@byrow :ca2 = filter_signal(:ca2, "bandstop", (59, 61), "fir", 60, nothing, false))
end


#= ################################
make event timeseries 
################################### =#

# edata(e) = isempty(lfp_data[lfp_data.file_id.== e.file_id, :lfp])  ? missing : lfp_data[lfp_data.file_id.== e.file_id, :lfp][1]
# eets(d,start,stop)=ismissing(d) ? missing : EventTimeSeries(d, start, stop) 
# event_info.lfp=[eets(edata(e), e.onset, e.offset) for e ∈ eachrow(event_info)]
# event_info = @subset(event_info, :lfp .!= ismissing)


m(e) = isempty(lfp_data[lfp_data.file_id.== e.file_id, :mob])  ? 0 : 1
miss=[m(e) for e ∈ eachrow(event_info)]
event_info[findall(miss.==1),:]

#=################################
#   Merge Event and Trial Info
############################### =#

function whichtrial(trials, event)
    @subset(trials,:file_id.==event.file_id, :start_time .≤ event.onset, event.offset .≤ :end_time) |>x->isempty(x) ? missing : x[1,:trial_id]
end

function whichtrials(trials, events)
    transform(events, AsTable(:)=>ByRow(x->whichtrial(trials,x))=>:trial_id)
end


event_info = @chain begin 
    whichtrials(trial_info, event_info)
    select(:event_id, :trial_id, :) 
end

@subset(event_info, :trial_id .== ismissing)


# ###################################

#= ###############################
#   Make New Data Set
################################### =#

trials=Dict{Int32, Any}()
@eachrow trial_info begin 
    t=Trial(:file_id, :date, :rat,:trial_type, :agent_type, :start_time, :end_time)
    push!(trials, t.id=>t)
end

events=Dict{Int32, Any}()
@eachrow event_info begin 
    # file, date, rat,  trial_type, agent, behavior, start_time, end_time
    e=Event(:file_id, :date, :rat, :trial_type, :agent_type, :behavior_type, :onset, :offset)
    push!(events, e.id=>e)
end

lfps=Dict{Int32, Any}()
@eachrow lfp_data begin 
    la=LFP(:file_id, :date, :rat, "amyg", :amyg)
    lb=LFP(:file_id, :date, :rat, "mob", :mob)
    lc=LFP(:file_id, :date, :rat, "ca2", :ca2)
    push!(lfps, la.id=>la)
    push!(lfps, lb.id=>lb)
    push!(lfps, lc.id=>lc)
end

alldata=merge(trials, events, lfps)

ds=DataSet(alldata)

save("Data/freeroam_dataset.jld2", "freeroam_dataset", ds)

