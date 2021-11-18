# include("utils.jl")
using DataFrames, DataFramesMeta, CSV, Chain, Missings

# function reformatter()
#     df=CSV.File("Data/CSV/event_info/behavioral_events.csv"; header=false)|>DataFrame
#     if ~isa(df[1,3], String)
#         print("already run")
#         return 
#     else
#         @chain begin 
#             CSV.File("Data/CSV/event_info/behavioral_events.csv"; header=false)
#             DataFrame
#             @transform(@byrow :Column3 = minstr2seconds(:Column3))
#             @transform(@byrow :Column4 = minstr2seconds(:Column4))
#             select(Not([5,6]))
#             CSV.write("Data/CSV/event_info/behavioral_events.csv", _, writeheader=false)
#         end
        
#         @chain begin 
#             CSV.File("Data/CSV/trial_info/behavioral_trials.csv"; header=false)
#             DataFrame
#             @transform(@byrow :Column3 = minstr2seconds(:Column3))
#             @transform(@byrow :Column4 = minstr2seconds(:Column4))
#             select(Not(5))
#             CSV.write("Data/CSV/trial_info/behavioral_trials.csv", _, writeheader=false)
#         end
        
#         @chain begin 
#             CSV.File("Data/CSV/event_info/freeroam_behavioral_events.csv"; header=false)
#             DataFrame
#             select(Not([5,6]))
#             CSV.write("Data/CSV/event_info/freeroam_behavioral_events.csv", _, writeheader=false)
#         end
        
#         @chain begin 
#             CSV.File("Data/CSV/trial_info/freeroam_behavioral_trials.csv"; header=false)
#             DataFrame
#             select(Not(5))
#             CSV.write("Data/CSV/trial_info/reeroam_behavioral_trials.csv", _, writeheader=false)
#         end
        
#         @chain begin 
#             CSV.File("Data/CSV/event_info/freeroam_approach_retreat_events.csv"; header=false)
#             DataFrame
#             select(Not([5,6]))
#             CSV.write("Data/CSV/event_info/freeroam_approach_retreat_events.csv", _, writeheader=false)
#         end
        
#         @chain begin 
#             CSV.File("Data/CSV/trial_info/freeroam_approach_retreat_trials.csv"; header=false)
#             DataFrame
#             @transform(@byrow :Column3 = Float64(:Column3))
#             @transform(@byrow :Column4 = Float64(:Column4))
#             select(Not(5))
#             CSV.write("Data/CSV/trial_info/freeroam_approach_retreat_trials.csv", _, writeheader=false)
#         end
#     end
# end

# reformatter()

@chain begin 
    CSV.File("Data/CSV/event_info/freeroam_approach_retreat_events.csv"; header=false)
    DataFrame
    select(Not([5,6]))
    # @subset(:Column1 .!= ismissing)
    CSV.write("Data/CSV/event_info/freeroam_approach_retreat_events.csv", _, writeheader=false)
end
