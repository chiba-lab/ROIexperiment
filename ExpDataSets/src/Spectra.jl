using StatsBase
const  fbands = (delta=[1,4], respiratory = [3, 12], theta=[5,10], alpha=[8, 12], beta = [15, 35], gamma_low = [50, 59], gamma_high = [70, 100])

abstract type AbstractSpectrum <:AbstractArrayData  end

function Base.values(spec::AbstractSpectrum)
    return spec.values
end

function axes(spec::AbstractSpectrum)
    return (freqs = spec.freqs)
end

function freqs(spec::AbstractSpectrum)
    return spec.freqs
end

function getband(spec::AbstractSpectrum, f₁, f₂) end
  
function getband(spec::AbstractSpectrum, b::Symbol)
    return getband(spec, fbands[b][1], fbands[b][2])
end


struct Spectrum <:AbstractSpectrum 
    values
    freqs 
end

function phase(spec::Spectrum)
   angle.(spec.values)
end

function amp(spec::Spectrum)
   abs.(spec.values)
end

function getband(spec::Spectrum, f₁, f₂)
    idx_range = findall(f₁.≤ spec.freqs .≤f₂)
    return Spectrum(spec.values[idx_range], spec.freqs[idx_range])
end


function bandpower(spec::Spectrum, f₁, f₂)
    return mean(power(getband(spec, f₁, f₂)))
end


function bandpower(spec::Spectrum, b::Symbol)
    return mean(power(getband(spec, b)))
end



struct PowerSpectrum <:AbstractSpectrum 
    values
    freqs 
end

function getband(spec::PowerSpectrum, f₁, f₂)
    idx_range = findall(f₁.≤ spec.freqs .≤f₂)
    return PowerSpectrum(spec.values[idx_range], spec.freqs[idx_range])
end

function power(spec::Spectrum)
   return PowerSpectrum(amp.(spec)^2, spec.freqs)
end

function bandpower(spec::PowerSpectrum, f₁, f₂)
    return StatsBase.mean(values(getband(spec, f₁, f₂)))
end


function bandpower(spec::PowerSpectrum, b::Symbol)
    return  bandpower(spec, fbands[b][1], fbands[b][2])
end

function binps(spec::PowerSpectrum, n=120, fmax=n)
    leftedge=collect(1:fmax/n:fmax)
    rightedge=collect(0:fmax/n:(fmax-1))
    pb=[bandpower(spec, rightedge[i], leftedge[i]) for i in 1:n]
    return PowerSpectrum(pb, leftedge)
end

using RollingFunctions

function smoothps(spec::PowerSpectrum, window)
    return PowerSpectrum(RollingFunctions.rollmean(spec.values, window),RollingFunctions.rollmean(spec.freqs,window))
end

# function  meanps(x::PowerSpectrum,y::PowerSpectrum)
#     if x.freqs == y.freqs
#         return PowerSpectrum(measurements((values(x).+values(y))/2), freqs(x))
#     else
#         return nothing 
#     end
# end 



