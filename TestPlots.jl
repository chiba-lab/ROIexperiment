include("./ROIontology.jl")

ROIexp = load("./Data/ROIexp.jld2", "ROIexp");
lfp = ROIexp[:lfp_data]
events = ROIexp[:behavioral_events]
trials = ROIexp[:trials]
sessions = ROIexp[:sessions]

include("./Links.jl")
##==============================================================================
const WC = collect(codom(window_data))
##=============================================================================


using Interpolations
using CairoMakie

dataregion = MOB
dataevent="Sniff"

CairoMakie.activate!(type="png")
idxs=findall(x->window_data_event(x).behavior.name == dataevent,WC) ∩ findall(x->window_data_region(x)==dataregion,WC)
f= Figure();
axs = [Axis(f[i, j]) for i = 1:3, j=1:3]
for n in 1:9
    idx=rand(idxs)
    test = 1:1:length(WC[idx])|>x->map(x->sin(16*2*pi*x/1010.1),x)
    S=TFamplitude(WC[idx], 256, 0, 1)
    itp=interpolate(rotl90(S.y), BSpline(Cubic(Line(OnGrid()))))
    x=1:1:size(S.y, 2)
    y=exp2.(range(0, log2(252.99), 100))
    
    heatmap!(axs[n],[itp(x1,y1) for x1 = x, y1 = y]; colormap = :jet)
    vlines!(axs[n], 1010; color=:orange, linewidth=2, linestyle=:dash)
    vlines!(axs[n],length(x)-1010;color=:red, linewidth=2, linestyle=:dash)
    axs[n].xticks = ([1010, length(x)-1010], ["start", "end"])
    axs[n].yticks=([21, 100],  ["8", "127"])
    e=window_data_event(WC[idx]).behavior.name
    region=window_data_region(WC[idx]).name
    axs[n].title = "$e $region"
end
f


# heatmap(1:.5:3357, .5:.5:60, rotl90(S.y); colormap = :jet)
x=inverse(window_data_event)(events[1])
S=TFanalyticsignal(x[1][1:1000], 256, 0, 1)
R=TFanalyticsignal(x[2][1:1000], 256, 0, 1)
using ComplexPhasePortrait

fig=Figure();
ax, h = heatmap(fig[1,1],angle.(rotl90(conj(R.y).*(S.y)))[:,1:60]; colormap = :roma, colorrange = (-π,π), backlight = 1f0, highclip = :black, interpolate=true)
Colorbar(fig[1, 2], h, ticks = ([-π, -π / 2, 0, π / 2, π], String.([-π, -π / 2, 0, π / 2, π])))
fig
##==============================================================================
# itp=interpolate(rotl90(S.y), BSpline(Cubic(Line(OnGrid()))))
# x=1:1:size(S.y, 2)
# y=exp2.(range(0, log2(252.99), 300))
# (253-1)/2+1

# f=Figure();
# ax=Axis(f[1, 1])
# heatmap!(ax,[itp(x1,y1) for x1 = x, y1 = y]; colormap = :jet)
# hidedecorations!(ax)
# f
#diverging_rainbow_bgymr_45_85_c67_n256, 
