include("ROIontology.jl")



using JLD2, DataFrames, DataFramesMeta, Statistics, Random

pactbl = load("PACtbl_ca2.jld2", "pactbl")
describe(pactbl)

##==============================================================================
function paccomvec(r, p)
    mean(r .* exp.(im .* p))
end

function shuffle(x)
    x[randperm(length(x))]
end

function nullpaccomvec(r, p)
    paccomvec(shuffle(r), p)
end

##==============================================================================

# sbs = @subset(pactbl, isa.(:trial_exp_num, Int))

# sbs.trial_fam = map(x -> x <= 2 ? "Novel" : "Familiar", sbs.trial_exp_num)
θ = [-π:2*π/180.0:π...]
pactbl.null_low_gamma_com_vec = map(x -> nullpaccomvec(x, θ[1:end-1]), pactbl.low_gamma)
pactbl.null_high_gamma_com_vec = map(x -> nullpaccomvec(x, θ[1:end-1]), pactbl.high_gamma)
pactbl.null_high_gamma_com_vec[1]

##==============================================================================
using HypothesisTests

#Sig of PAC
sbs = copy(pactbl)

sbs = @subset(sbs, :agent_type .!= "missing")
sbs.low_gamma_r_diff = abs.(sbs.low_gamma_com_vec) .- abs.(sbs.null_low_gamma_com_vec)
sbs.high_gamma_r_diff = abs.(sbs.high_gamma_com_vec) .- abs.(sbs.null_high_gamma_com_vec)
gdf = groupby(sbs, [:behavior_type, :agent_type, :region])
gdf = @combine(gdf, :low_gamma_pac_pval = pvalue(OneSampleTTest(:low_gamma_r_diff)), :high_gamma_pac_pval = pvalue(OneSampleTTest(:high_gamma_r_diff)))

using CSV, FileIO

CSV.write("./TempData/PACsig_ca2.csv", gdf)
##==============================================================================
PacHypoTbl = @select(pactbl, :rat, :region, :behavior_type, :agent_type, :trial_exp_num, :receiver_exp_num, :low_gamma, :high_gamma)
PacHypoTbl.trial_fam = map(x -> isa(x, Int) ? (x <= 2 ? "Novel" : "Familiar") : "missing", PacHypoTbl.trial_exp_num)
PacHypoTbl.receiver_fam = map(x -> isa(x, Int) ? (x <= 2 ? "Novel" : "Familiar") : "missing", PacHypoTbl.receiver_exp_num)
PacHypoTbl = stack(PacHypoTbl, [:low_gamma, :high_gamma], variable_name=:fband, value_name=:pac_vec)
PacHypoTbl.complex_avg = map(x -> paccomvec(x, θ[1:end-1]), PacHypoTbl.pac_vec)
PacHypoTbl.null_complex_avg = map(x -> nullpaccomvec(x, θ[1:end-1]), PacHypoTbl.pac_vec)
# PacHypoTbl.fband = map(x -> x == "low_gamma_com_vec" ? "Low Gamma" : "High Gamma", PacHypoTbl.fband)
PacHypoTbl.mag = map(x -> abs(x), PacHypoTbl.complex_avg)
PacHypoTbl.phase = map(x -> angle(x), PacHypoTbl.complex_avg)
PacHypoTbl.null_mag = map(x -> abs(x), PacHypoTbl.null_complex_avg)
PacHypoTbl.null_phase = map(x -> angle(x), PacHypoTbl.null_complex_avg)
PacHypoTbl[:, Not([:trial_exp_num, :pac_vec])]
CSV.write("./TempData/PACHypoTbl_ca2.csv", PacHypoTbl)
