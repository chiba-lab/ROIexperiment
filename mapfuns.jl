



function map_ets(fun, funname, window)
    @chain et begin
    # @transform!(:lfp => ByRow.(x -> normalize(x)) => :lfp)
    @subset((:end_time .- :start_time) .≥ 1.0)
    transform(:lfp => ByRow.(x -> normalize(x)) => :lfp)
    transform(:lfp => ByRow.(fun) => funname)
    # @select(:rat, :trial_type, :agent, :behavior, :region, :pspec)
end
end



function meanreduce(et, col, subs, grp)
    resame=Symbol("mean_"*String(col))
    @chain et begin
    @subset(subs...)
    groupby([grp])
    @combine(resname = [meanandstderr(col)])
    transform(:mps => ByRow.(x -> Measurements.value.(x)) => :mpsvals, :mps => ByRow.(x -> Measurements.uncertainty.(x)) => :err)
end
end 

behavior_mps = @chain ps begin
    @subset(:behavior .!= "null")
    transform(:pspec => ByRow.(x -> values(x)) => :pspec)
    groupby([:behavior, :region])
    @combine(:mps = [meanandstderr(:pspec)])
    transform(:mps => ByRow.(x -> Measurements.value.(x)) => :mpsvals, :mps => ByRow.(x -> Measurements.uncertainty.(x)) => :err)
    select(:behavior, :region, :mpsvals => :mps, :err)
end