using Random, Statistics, StatsBase

function paccomvec(r, p)
    mean(r .* exp.(im .* p))
end

function shuffle(x)
    x[randperm(length(x))]
end

function nullpaccomvec(r, p)
    paccomvec(shuffle(r), p)
end


θ= [0:2*π/180.0:2π...] 
A= .5*sin.(θ).+1 .*(ones(1, length(θ)).*vcat(fill(-1, 90), fill(1, 91)))

v=paccomvec(A, θ)
vn=nullpaccomvec(A, θ)


using Plots

p1 = scatter( θ, A, proj=:polar, title="Simulated Data")
scatter!([angle(v)], [abs(v)],  proj=:polar, label="Mean Vector")
plot!(θ, fill(mean(A),length(θ)) , proj=:polar, label="Null")
p2 = scatter(θ, shuffle(A), proj=:polar, title="Shuffled Data")
scatter!([angle(vn)], [abs(vn)], proj=:polar, label="Mean Vector")
p=plot(p1, p2)



using FileIO
save("./TempData/data_shuffling_explainer.pdf", p)