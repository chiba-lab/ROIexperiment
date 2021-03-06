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
    S=FA.TFamplitude(WC[idx], 256, 0, 1)
    itp=interpolate(rotl90(S.y), BSpline(Cubic(Line(OnGrid()))))
    x=1:1:size(S.y, 2)
    y=exp2.(range(0, log2(239), 200))
    
    heatmap!(axs[n],[itp(x1,y1) for x1 = x, y1 = y]; colormap = :jet)
    vlines!(axs[n], 1010; color=:orange, linewidth=2, linestyle=:dash)
    vlines!(axs[n],length(x)-1010;color=:red, linewidth=2, linestyle=:dash)
    axs[n].xticks = ([1010, length(x)-1010], ["start", "end"])
    axs[n].yticks = ([1,99, 152, 200], ["1", "8", "32", "120"])
    e=window_data_event(WC[idx]).behavior.name
    region=window_data_region(WC[idx]).name
    axs[n].title = "$e $region"
end
f

##==============================================================================
h = heatmap(f[j-1,i], collect((1:1:l)./1010.1), 1:.1:128, itp[1:1:l, 1:.1:128] ; colormap = :romaO, colorrange = (-π,π), backlight = 1f0)
##==============================================================================
# heatmap(1:.5:3357, .5:.5:60, rotl90(S.y); colormap = :jet)
using CairoMakie
using Interpolations
using ProgressBars
import DSP as dsp
CairoMakie.activate!(type="png")

x=inverse(window_data_event)(events[1])

f= Figure();
hs=Dict()
count=1
# axs = [Axis(f[i, j]) for i = 1:4, j=1:4]
for i in tqdm(1:4)
        for j in i:4
                e=window_data_event(x[1]).behavior.name
                r1=window_data_region(x[i]).name
                r2=window_data_region(x[j]).name
                # x1=dsp.resample.(x, 1/10)

                S=FA.TFanalyticsignal(x[i], 256, 0, 1)
                R=FA.TFanalyticsignal(x[j], 256, 0, 1)
        
                l= size(S.y, 2)
                if i==j
                        ph = angle.(rotl90(S.y))
                else
                        ph = angle.(rotl90(conj(R.y).*(S.y)))
                end
                # amp = abs.(rotl90(conj(R.y).*(S.y)))[:, 1:128]

                itp=interpolate(ph, BSpline(Cubic(Line(OnGrid()))))

                ax, hs[count] = heatmap(f[j,i], collect((1:1:l)./1010.1), 1:.1:253, itp[1:1:l, 1:.1:253] ; colormap = :romaO, colorrange = (-π,π), backlight = 1f0)
                ax.yticks = ([1,253], ["1","127"])
                if i == 1 
                        ax.ylabel = "Freq (Hz)"
                end
                if j == 4
                        ax.xlabel = "Time (s)"
                end
                count+=1
             
                ax.title = "$e, $r1 × $r2"
                vlines!(ax, 1; color=:green, linewidth=2, linestyle=:dash)
                vlines!(ax, l/1010.1-1;color=:red, linewidth=2, linestyle=:dash)
                # Colorbar(h, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))
                # amp = abs.(rotl90(conj(R.y).*(S.y)))[:, 1:128]
                # ax2, h2 = heatmap(f[j-1,i+1], collect((1:1:l)./1010.1), 1:.1:128, amp[1:1:l, 1:.1:128] ; colormap = :jet)


        end
end
# Colorbar(f[3,4], h, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))
Colorbar(f[:,5], hs[1], ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))

f
##==============================================================================
x=inverse(window_data_event)(events[11])

f= Figure();
# axs = [Axis(f[i, j]) for i = 1:4, j=1:4]
for i in tqdm(1:4)
        for j in i:4
                e=window_data_event(x[1]).behavior.name
                r1=window_data_region(x[i]).name
                r2=window_data_region(x[j]).name
                # x1=dsp.resample.(x, 1/10)
                test = 1:1:length(x[1])|>x->map(x->sin(16*2*pi*x/1010.1),x)

                S=FA.TFanalyticsignal(x[i], 256, 0, 1)
                R=FA.TFanalyticsignal(x[j], 256, 0, 1)
        
                l= size(S.y, 2)
                # ph = angle.(rotl90(conj(R.y).*(S.y)))[:, 1:128]
                amp = abs.(rotl90(conj(R.y).*(S.y)))

                itp=interpolate(amp, BSpline(Cubic(Line(OnGrid()))))
                yind=exp2.(range(0, log2(252.99), 100))
                xind=1:1:size(S.y, 2)


                ax, h = heatmap(f[j,i], [itp(x1,y1) for x1 = xind, y1 = yind]; colormap =cgrad(:jet, scale=linear))                # ax.yticks = ([1,64,128], ["1","36", "64"])
                ax.yticks = ([], [])
                if i == 1 
                        ax.ylabel = "Freq (Hz)"
                        ax.yticks = ([1,37, 100], ["1", "8", "127"])
                end
                if j == 4
                        ax.xlabel = "Time (s)"
                end
                if i == 1 && j == 1
                         ax.title = "$e\n $r1 × $r2"
                else
                        ax.title = "$r1 × $r2"
                end

             
                vlines!(ax, 1010.1; color=:green, linewidth=2, linestyle=:dash)
                vlines!(ax, l-1010.1;color=:blue, linewidth=2, linestyle=:dash)
                # Colorbar(h, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))
             
                ax.xticks = ([], [])

                # ax2, h2 = heatmap(f[j-1,i+1], collect((1:1:l)./1010.1), 1:.1:128, amp[1:1:l, 1:.1:128] ; colormap = :jet)

        end
end
# Colorbar(f[3,4], h, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))
# Colorbar(f[:,4], h, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))

f
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
##==============================================================================
using CairoMakie
let
        x = y = -2:0.005:2
        f(z) = 1 / (z * (z^2 - z - 1 - 3im))
        fvals = [f(u + 1im * v) for u in x, v in y]
        fvalues = abs.(fvals)
        fargs = angle.(fvals)
        polya(x, y) = Point2f(real(f(x + 1im * y)), -imag(f(x + 1im * y)))

        fig = Figure(resolution = (900, 400))
        axs = [Axis(fig[1, i], aspect = 1) for i in 1:2]
        cmap = :roma
        streamplot!(axs[1], polya, -2 .. 2, -2 .. 2, colormap = (:black, :red),
                gridsize = (40, 40), arrow_size = 6, linewidth = 1)
        pltobj2 = heatmap!(axs[2], x, y, fargs, colorrange = (-π, π), colormap = cmap)
        streamplot!(axs[2], polya, -2 .. 2, -2 .. 2, colormap = (:black, :black),
                gridsize = (40, 40), arrow_size = 6, linewidth = 1)
        Colorbar(fig[1, 3], pltobj2, ticks = ([-π, -π / 2, 0, π / 2, π],
                [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))
        limits!(axs[1], -2, 2, -2, 2)
        limits!(axs[2], -2, 2, -2, 2)
        colsize!(fig.layout, 1, Aspect(1, 1.0))
        colsize!(fig.layout, 2, Aspect(1, 1.0))
        display(fig)
end

##==============================================================================
using CairoMakie
f = Figure(fontsize = 18);

Axis(f[1, 1],
    title = L"\frac{x + y}{\sin(k^2)}",
    xlabel = L"\sum_a^b{xy}",
    ylabel = L"\sqrt{\frac{a}{b}}"
)

f
##==============================================================================
dataregion = CA2
dataregion2 = MOB
dataevent="Sniff"

idxs_r1=findall(x->window_data_event(x).behavior.name == dataevent,WC) #∩ findall(x->window_data_region(x)==dataregion2,WC) ∩ findall(x->window_data_region(x)==dataregion,WC)

wd_pre=GMap(x->x[498:1010+511], WC)
pre_tf_r1=GMap(x->FA.TFanalyticsignal(x, 256, 0, 1), wd_pre.(WC[idxs]))


pre_tf_amp = map(pre_tf.(wd_pre.(WC[idxs]))) do x
        abs.(x.y)
end


pre_tf_amp[1]
using Statistics
m = mean(pre_tf_amp)

itp=interpolate(rotl90(m), BSpline(Cubic(Line(OnGrid()))))
yind=exp2.(range(0, log2(239), 200))
xind=1:1:size(m,2)

fl(x)=x/2+.5

f=Figure();
ax=Axis(f[1, 1])
heatmap!(ax, [itp(x1,y1) for x1 = xind, y1 = yind]; colormap = :jet)
vlines!(ax, 512, color=:black, linewidth=4, linestyle=:dash)
ax.ylabel = "Freq (Hz)"
ax.xlabel = "Time (s)"
ax.yticks = ([1,99, 152, 200], ["1", "8", "32", "120"])
ax.xticks = ([12,512,1012], ["-.5", "onset", ".5"])

e=dataevent
region=dataregion.name
ax.title = "$e $region"
f


