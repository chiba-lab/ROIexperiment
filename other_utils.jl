function meanandstderr(v)
    measurement.(mean(v), (StatsBase.std(v))/sqrt(length(v)))
end

# function meanps(df, by)
#     ps=transform(df,  :lfp => ByRow.(x->values(binps(compute_power_spectrum(x))))=>:pspec)
#     g=groupby(ps, by)
#     @combine(g, mps=[meanandstderr(:pspec)])
# end

# function plotps(mps, args...)
#     plot(Measurements.value.(mps), ribbon=Measurements.uncertainty.(mps), fillalpha=.25, args...)
# end

# function plotps!(mps, args...)
#     plot!(Measurements.value.(mps), ribbon=Measurements.uncertainty.(mps), fillalpha=.25, args...)
# end

function agent_ps_plot(bmps, nmps, b, r)
    freqs=1:1:100 
    s=@subset(bmps, :behavior.== b, :region .== lowercase(r))
    s=@orderby(s, :agent)
    ns=@subset(nmps, :region.== lowercase(r))
    specs = s.mps
    errs =  s.err
    l =reshape(s.agent,1,length(s.agent))#|>x->map(y->"Agent: "*y, x)
    title=b*" "*r
    p=plot(freqs, specs, ribbon=errs, fillalpha=.25, labels=l, xaxis=:log, yaxis=:log, legend=false, linewidth=2)
    plot!(p, freqs, ns.mps, ribbon=ns.err, fillalpha=.25, labels="Baseline", fillcolor=:black, linestyle=:dash, linecolor=:black)
    return p
end



