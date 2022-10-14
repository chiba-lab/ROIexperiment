include("ROIontology.jl")



using JLD2, DataFrames, DataFramesMeta, Statistics

pactbl = load("PACtbl_ca2.jld2", "pactbl")
describe(pactbl)

sbs = copy(pactbl)
# sbs.trial_fam = map(x -> x <= 2 ? "Novel" : "Familiar", sbs.trial_exp_num)





gdf = groupby(sbs, [:behavior_type, :region])
gdf = @combine(gdf, :lgpac = [mean(:low_gamma)], :hgpac = [mean(:high_gamma)]) |> x -> sort(x, :region) |> x -> groupby(x, [:behavior_type])

##==============================================================================
using Plots, UncertainData, Plots, CircularArrays, StatsBase


θ = CircularArray([-π:2*π/180.0:π...])
function circmovmean(x, w=1, it=1)
    cx = CircularArray(x)
    for i in 1:length(cx)*it
        cx[i] = mean(cx[i-w:i+w])
    end
    return cx
end


noneg(x, m) = x < 0 ? m : x
lg=gdf[1][:,:lgpac]

Plots.plot(θ[1:end], lg[1], proj=:polar)




using FileIO
# p2=plot
fnames = []
for g in gdf
    p1 = Plots.plot()
    p2 = Plots.plot()
    # p1 = Plots.plot!(p1, θ[1:end-1], fill(1 / 180, 180), proj=:polar, linewidth=2, color="black")
    for i in 1:size(g, 1)
        reg = g[i, :region]
        hg = g[i, :hgpac] |>x->standardize(ZScoreTransform, x)|>x->noneg.(x,mean(x))|>x -> circmovmean(x, 4, 4)
        # hg = g[i, :hgpac] |> x -> x ./ sum(x) #|> x -> noneg.(x) |> x -> circmovmean(x, 4, 4)
        p1 = Plots.plot!(p1, θ[1:end-1], hg, proj=:polar, label=reg, linewidth=2, fill=(0,.1))


        # p2 = plot!(p2, θ[1:end-1], hg, proj=:polar, label=reg, linewidth=2, fill=(0, 0.1), title="High Gamma")
    end

    ttl = "High Gamma\n" * g[1, :behavior_type] 
    p = Plots.plot!(p1, title=ttl)
    fname = "PACplots/" * g[1, :behavior_type]* "_PAC: Ca2 Theta - High Gamma.pdf"
    savefig(fname)

    fnames = [fname, fnames...]
end

fnames
using PDFmerger, Glob

# fn = ["PACplots/PACplots_high_gamma_nov - fam.pdf", "PACplots/PACplots_high_gamma.pdf", "PACplots/PACplots_low_gamma.pdf"]

# merge_pdfs(fn, "PACplots_low_gamma_ca2.pdf", cleanup=true)




# pol2cart(r, θ) = [r .* cos(θ) r .* sin(θ)]
# wXY = map(x -> vcat(map(pol2cart, vcat(x..., x[1]...), θ)...) |> x -> (vec(x[:, 1]), vec(x[:, 2])), gdf[1].lgpac)


# using StatsBase

# f = Figure()
# a = Axis(f[1, 1])

# # lines!(a, R * cos.(θ), R * sin.(θ), color="black", linewidth=0.5)
# noneg(x) = x < 0 ? 0 : x
# Rfit = map(x -> collect(noneg.(x)), [circmovmean(standardize(ZScoreTransform, r), 4, 4) for r in gdf[1].lgpac])
# [(a, r .* cos.(θ), r .* sin.(θ)) for r in Rfit]


# f
# # lines!(a, R*cos.(θ), R*sin.(θ), color=:black)
# a.aspect = AxisAspect(1)
# hidedecorations!(a)
# f
# mx = maximum(maximum.(filter.(!isnan, Rfit)))
# xlims!(a, (-1.1 * mx, 1.1 * mx))
# ylims!(a, (-1.1 * mx, 1.1 * mx))

# f
# lines!(a, 1.05 * mx .* cos.(θ[1:end+1]), 1.05 * mx .* sin.(θ[1:end+1]), color=:black)
# f