
using FourierAnalysis, ProgressBars, JLD2
# @load "./TempData/event_lfp_windows_1_1.jld2" WC



# WCshort=filter(x->length(x) <= 2^13, WC)
# wd_analytic_dict=Dict{Vector{Float64},TFAnalyticSignal}()
# cnt=1
# for w in tqdm(WCshort)
#     wd_analytic_dict[w] = iseven(length(w)) ? TFanalyticsignal(w, 256,  0, 1) : TFanalyticsignal(w[1:end-1], 256,  0, 1)  
#     if cnt%500==0||cnt==length(WCshort)  
#         jldopen("TempData/event_tf_widows_1_1_e13.jld2", "a+") do file 
#             file["wd_analytic_dict_$cnt"]=wd_analytic_dict
#         end
#         wd_analytic_dict=Dict{Vector{Float64},TFAnalyticSignal}()
#     end
#     cnt=cnt+1
# end



data = load("./TempData/event_tf_widows_1_1_e13.jld2", "wd_analytic_dict_500")
tfs = values(data)
tfs=collect(tfs)
TFphase(tfs[1]).flabels





