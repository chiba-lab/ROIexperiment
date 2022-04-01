
##==============================================================================
include("./ROIontology.jl")

##==============================================================================
#Load data

ROIexp = load("./Data/ROIexp_FR.jld2", "ROIexp");
lfp = ROIexp[:lfp_data]
events = ROIexp[:behavioral_events]
trials = ROIexp[:trials]
sessions = ROIexp[:sessions]

include("./Links.jl")

##==============================================================================
# #univariate analysis

using FourierAnalysis

const WD = filter(x -> window_data_event(window_data(x)) ∈ dom(event_trial), collect(dom(window_data)))
const WC = filter(x -> window_data_event(x) ∈ dom(event_trial), collect(codom(window_data)))

# WT = event_trial∘window_data_event
# using JLD2 
# @save "TempData/event_lfp_windows_1_1.jld2"  WC
##==============================================================================
#split events into pre-peri-post with 512 sample pre and post windows
#events shorter than 512 are expanded
# wd_pre = GMap(x -> x[498:1010], WC)
wd_peri = GMap(x -> length(x) <= 2532 ? x[round(Int, length(x) / 2)-256:round(Int, length(x) / 2)+256] : x[1010:end-1010], WC)
# wd_post = GMap(x -> x[end-1010:end-498], WC)

# #Freq Dom
# # pre_spec = GMap(x -> spectra(x, 256, 512), wd_pre.(WC))
peri_spec = GMap(x -> spectra(x, 256, 512), wd_peri.(WC))
# # post_spec = GMap(x -> spectra(x, 256, 512), wd_post.(WC)) 

using DataFrames, DataFramesMeta

event_spectra = DataFrame()

event_spectra.rat = map(x -> window_lfp(inv(window_data)(x)).rat.name, WC)
event_spectra.region = map(x -> window_data_region(x).name, WC)
event_spectra.behavior_type = map(x -> window_data_event(x).behavior.name, WC)
event_spectra.trial_type = map(x -> event_trial(window_data_event(x)).condition.name, WC)
event_spectra.agent_type = map(x -> agenttype(event_trial(window_data_event(x)).condition), WC)
specs = map(x -> spectra(x, 256, 512), WC)
event_spectra.respiratory = map(x -> mean(x, (3, 12)), specs)
event_spectra.theta = map(x -> mean(x, (5, 10)), specs)
event_spectra.beta = map(x -> mean(x, (15, 35)), specs)
event_spectra.gamma_low = map(x -> mean(x, (50, 59)), specs)
event_spectra.gamma_high = map(x -> mean(x, (70, 100)), specs)

esnull = copy(event_spectra)
esnull.behavior_type = fill("null", length(esnull.behavior_type))

event_spectra = vcat(event_spectra, esnull)


##==============================================================================

using MixedModels

hyposub = @subset(event_spectra, :behavior_type .∈ Ref(["null", "Groom"]), :region .== "MOB")
fm = @formula(respiratory ~ 1 + behavior_type * 1 + (1 | rat))
fm1 = fit(MixedModel, fm, hyposub, contrasts = Dict(:behavior_type => DummyCoding(; base = "null")))