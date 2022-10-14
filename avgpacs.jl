include("./FamiliarityTable.jl")
using FourierAnalysis
##==============================================================
elist = unique(eventtbl.event)
# event_receiver_name = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_name[i]) for i in 1:size(eventtbl, 1)]))
# event_receiver_exp_num = GMap(Dict([(eventtbl.event[i] => eventtbl.receiver_exp_num[i]) for i in 1:size(eventtbl, 1)]))
# event_trial_exp_num = GMap(Dict([(eventtbl.event[i] => eventtbl.trial_exp_num[i]) for i in 1:size(eventtbl, 1)]))
makeeven(x) = length(x) % 2 == 0 ? x : x[1:end-1]
function peritrim(x)
    if length(x) <= 2532
        return makeeven(x[round(Int, length(x) / 2)-256:round(Int, length(x) / 2)+256])
    elseif (length(x) - 2020) > 2^14
        return makeeven(x[1010:1010+2^14])
    else
        return makeeven(x[1010:end-1010])
    end
end

wd_peri = GMap(x -> peritrim(x), WC)

length.(wd_peri.(WC)) .% 2
##==============================================================
using ProgressBars
inef = inverse(window_data_event)
using BenchmarkTools
ecrtf = Dict()
Threads.@threads for e in tqdm(unique(eventtbl.event))
    # e = ev[ei]
    # for e in unique(eventtbl.event[1:500])
    eWC = inef(e) |> x -> isa(x[1], Vector) ? x : [x]
    rWC = map(x -> window_data_region(x).name, eWC)
    idx = findfirst(x -> x == "Ca2", rWC)
    if !isnothing(idx)
        Y = TFanalyticsignal(wd_peri(eWC[idx[1]]), 256, 0, 0.5; fmax=4)
        ecrtf[e] = vec(mean(extract(Y, (3, 12) .* 256 ./ 1010.1, :), dims=1))
    end

end
ecrtf[eventtbl.event[1]]



##==============================================================


# Threads.@threads for e in tqdm(unique(eventtbl.event))
#     # e = ev[ei]
#     # for e in unique(eventtbl.event[1:500])
#     eWC = inef(e) |> x -> isa(x[1], Vector) ? x : [x]
#     rWC = map(x -> window_data_region(x).name, eWC)

#     cWC = collect(coherence(hcat(wd_peri.(eWC)...), 256, 512).y)
#     # if length(rWC) == 4
#     #     print(ei)
#     # end
#     for i in 1:length(rWC)

#         for j in 1:length(rWC)
#             codict[e, rWC[i], rWC[j]] = map(x -> x[i, j], cWC)
#         end
#     end
# end


# codict[unique(eventtbl.event)[483], "Insula", "Ca2"]
# emrtf[eventtbl.event[1]]
##==============================================================
using NaNStatistics
elgtf = Dict()
ehgtf = Dict()
eWC = inef(elist[1]) |> x -> isa(x[1], Vector) ? x : [x]
frat = 256 / 1010.1


Threads.@threads for e in tqdm(elist)
    # e = ev[ei]
    # for e in unique(eventtbl.event[1:500])
    eWC = inef(e) |> x -> isa(x[1], Vector) ? x : [x]
    rWC = map(x -> window_data_region(x).name, eWC)
    lgWC = TFanalyticsignal.(wd_peri.(eWC), 256, 0, 0.5; fmin=48 * frat, fmax=62 * frat) |> Y -> [vec(mean(extract(y, (50, 59) .* frat, :), dims=1)) for y in Y]
    hgWC = TFanalyticsignal.(wd_peri.(eWC), 256; fmin=60 * frat, fmax=110 * frat) |> Y -> [vec(mean(extract(y, (70, 100) .* frat, :), dims=1)) for y in Y]
    # if length(rWC) == 4
    #     print(ei)
    # end
    for i in 1:length(rWC)

        elgtf[e, rWC[i]] = lgWC[i]
        ehgtf[e, rWC[i]] = hgWC[i]
    end
    # for i in 1:length(rWC)
    #     Y = TFanalyticsignal(wd_peri(eWC[i]), 256, 0, 0.5; fmin=45, fmax=59)
    #     elgtf[e, rWC[i]] = vec(mean(extract(Y, (50, 59) .* 256 ./ 1010.1, :), dims=1))
    #     ehgtf[e, rWC[i]] = vec(mean(extract(Y, (70, 100) .* 256 ./ 1010.1, :), dims=1))
    # end
end
# elgtf = Dict()
# ehgtf = Dict()
# Threads.@threads for e in tqdm(elist)
#     # e = ev[ei]
#     # for e in unique(eventtbl.event[1:500])
#     eWC = inef(e) |> x -> isa(x[1], Vector) ? x : [x]
#     rWC = map(x -> window_data_region(x).name, eWC)
#     for i in 1:length(rWC)
#         Y = TFanalyticsignal(wd_peri(eWC[i]), 256, 0, 0.5; fmax=30)
#         elgtf[e, rWC[i]] = vec(mean(extract(Y, (50, 59) .* 256 ./ 1010.1, :), dims=1))
#         ehgtf[e, rWC[i]] = vec(mean(extract(Y, (70, 100) .* 256 ./ 1010.1, :), dims=1))
#     end
# end
elgtf[eventtbl.event[1], "Amygdala"]

##==============================================================
event_ca2_respiratory_tf = GMap(ecrtf)
event_low_gamma_tf = GMap(elgtf)
event_high_gamma_tf = GMap(ehgtf)
##==============================================================================
using UncertainData
epaclg = Dict()
epachg = Dict()

pacelist = collect(dom(event_ca2_respiratory_tf))

Threads.@threads for e in tqdm(pacelist)
    # e = ev[ei]
    # for e in unique(eventtbl.event[1:500])
    ca2p = angle.(event_ca2_respiratory_tf(e))
    eWC = inef(e) |> x -> isa(x[1], Vector) ? x : [x]
    rWC = map(x -> window_data_region(x).name, eWC)
    for i in 1:length(rWC)
        epaclg[e, rWC[i]] = bin(mean, -π:2*π/180.0:π, ca2p, abs.(elgtf[e, rWC[i]]))
        epachg[e, rWC[i]] = bin(mean, -π:2*π/180.0:π, ca2p, abs.(ehgtf[e, rWC[i]]))
    end
end


epaclg[pacelist[1], "Amygdala"]
event_low_gamma_pac = GMap(epaclg)
event_high_gamma_pac = GMap(epachg)
##==============================================================================
function paccomvec(r, p)
    mean(r .* exp.(im .* p))
end
##==============================================================================
phs = collect(-π:2*π/180.0:π)
pactbl = DataFrame()
pacdom = [x for x = dom(event_low_gamma_pac)] ∩ [x for x = dom(event_high_gamma_pac)]
pactbl.event = map(x -> x[1], pacdom)
pactbl.region = map(x -> x[2], pacdom)
pactbl.low_gamma = map(x -> event_low_gamma_pac((x[1], x[2])), pacdom)
pactbl.high_gamma = map(x -> event_high_gamma_pac((x[1], x[2])), pacdom)
pactbl = pactbl[findall(x -> !any(isnan.(x)), pactbl.high_gamma)∩findall(x -> !any(isnan.(x)), pactbl.low_gamma), :]
pactbl.low_gamma_com_vec = map(x -> paccomvec(x, phs[1:end-1]), pactbl.low_gamma)
# pactbl.low_gamma_com_mag = map(x -> abs(x), pactbl.low_gamma_com_vec)
# pactbl.low_gamma_com_theta = map(x -> angle(x), pactbl.low_gamma_com_vec)
pactbl.high_gamma_com_vec = map(x -> paccomvec(x, phs[1:end-1]), pactbl.high_gamma)
# pactbl.high_gamma_com_mag = map(x -> abs(x), pactbl.high_gamma_com_vec)
# pactbl.high_gamma_com_theta = map(x -> angle(x), pactbl.high_gamma_com_vec)
pactbl = leftjoin(pactbl, eventtbl, on=[:event, :region])
agmap = GMap(Dict([(Object, "Object"), (Rat, "Rat"), (Robot, "Robot"), (Missing, "missing")]))
pactbl.agent_type = agmap.(pactbl.agent_type)
using JLD2
save("PACtbl_ca2.jld2", "pactbl", pactbl)
##==============================================================================




