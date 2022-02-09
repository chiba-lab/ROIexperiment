include("./ROIontology.jl")
using Chain 

ROIexp = load("./Data/ROIexp.jld2", "ROIexp")
lfp = ROIexp[:lfp_data]
events  = ROIexp[:behavioral_events]


function get_event_lfp_recording(event, lfp)
    f1(l) = in_interval(event, l) 
    filter(x -> f1(x), lfp)
end

get_event_lfp_recording(events[2], lfp)



time_interval(events[1])
typeof(events[1]) <: Observation
ObservationInterval
BehavioralEvent <: Observation
events
# value(lfp[2])
# dims(lfp[2])

# Table(lfp)

# import DataFrames as df

# LFPdf = df.DataFrame(lfp)

# using DataFramesMeta

# df.groupby(LFPdf, :region)

# transform!(LFPdf, :lfp => ByRow.(x -> spectra(x, 256, 512)) => :spectra)



# LFPdf.spectra[1]

# using Plots
# spec = LFPdf.spectra[1:5]