include("./ROIontology.jl")

ROIexp=load("./Data/ROIexp.jld2", "ROIexp")


lfp=ROIexp[:lfp_recordings]
lfp.region
lfp.region