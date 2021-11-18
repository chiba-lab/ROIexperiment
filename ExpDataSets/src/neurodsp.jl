# * Filters
# const filt = pyimport("neurodsp.filt")

function filter_signal(ts::Array, args...)
    return filt.filter_signal(ts, args...)
end

function filter_signal(ts::ContinuousTimeSeries, args...)
    return ContinuousTimeSeries(filter_signal(collect(values(ts)), ts.fs, args...), ts.fs, ts.start_time, ts.end_time)
end


# # * Time Frequency Analysis 
# tf = pyimport("neurodsp.timefrequency")
 
# function robust_hilbert(sig::Array, args...)
#     return tf.robust_hilbert(sig, args...)
# end

# function robust_hilbert(sig::KeyedArray, args...)
#     return similar(sig) .= robust_hilbert(Array(sig), args...)
# end

# function phase_by_time(sig::Array, args...)
#     return tf.phase_by_time(sig, args...)
# end

# function phase_by_time(sig::KeyedArray, args...)
#     return tf.phase_by_time(sig, args...)
# end

# function phase_by_time(sig::KeyedArray, args...)
#     return similar(sig) .= phase_by_time(Array(sig), args...)
# end

# function freq_by_time(sig::Array, args...)
#     return tf.freq_by_time(sig, args...)
# end

# function freq_by_time(sig::KeyedArray, args...)
#     return similar(sig) .= freq_by_time(Array(sig), args...)
# end

# function amp_by_time(sig::Array, args...)
#     return tf.amp_by_time(sig, args...)
# end

# function amp_by_time(sig::KeyedArray, args...)
#     return similar(sig) .= amp_by_time(Array(sig), args...)
# end

# # * wavelet 

# function compute_wavelet_transform(sig::Array, fs, freqs, args...)
#     return tf.compute_wavelet_transform(sig, fs, freqs, args...)
# end

# function compute_wavelet_transform(sig::KeyedArray, fs, freqs, args...)
#     return KeyedArray(compute_wavelet_transform(Array(sig), fs, freqs, args...); time=keys(sig), freq = freqs)
# end

# # * Spectral Analysis

# const spc = pyimport("neurodsp.spectral")

function compute_power_spectrum(ts::Array, args...)
    return spc.compute_spectrum(ts, args...)
end

function compute_power_spectrum(ts::ContinuousTimeSeries, args...)
    f,p=compute_power_spectrum(collect(values(ts)), ts.fs, args...)
    return PowerSpectrum(p, f)
end

# function compute_spectrum(sig::EventTimeSeries, args...)
#     specfreq=mapslices(x->all(isnan.(x)) ? NaN : compute_spectrum(x, sig.data.fs, args...), sig.data.data; dims=1)
#     return PowerSpectrum(nancat(map(x->x[2], specfreq)...; dims=2), specfreq[1][1], sig.data.channels)
# end


# function bandpower(spec::PowerSpectrum, bands; norm=false)
#     if norm
#         totalpower=sum(spec.data, dims=1)
#         bspec=[sum(spec.data[findall(b[1] .< spec.freqs .< b[2]),:],dims=1)./totalpower for b ∈ bands]|>x->vcat(x...)
#     else
#         bspec=[mean(spec.data[findall(b[1] .< spec.freqs .< b[2]),:],dims=1) for b ∈ bands]|>x->vcat(x...)
#     end
#     BandPowerSpectrum(bspec, bands, spec.channels)
# end



# function compute_scv(sig::Array, args...)
#     return spc.compute_scv(sig, args...)
# end 

# function compute_scv(sig::KeyedArray, args...)
#     (freq, cv) = compute_scv(Array(sig), args...)
#     return KeyedArray(cv; freq=freq, scv=cv)
# end

# function compute_scv_rs(sig::Array, args...)
#     return spc.compute_scv(sig, args...)
# end 

# function compute_scv_rs(sig::KeyedArray, args...)
#     (freq, t_ind, cv_rs) = compute_scv_rs(Array(sig), args...)
#     if isempty(t_ind)
#         return KeyedArray(cv_rs; freq=freq, sample=1:size(cv_rs,2))
#     else
#         return KeyedArray(cv_rs; freq=freq, time=t_ind.+keys(sig)[1])
#     end
# end

# function compute_spectral_hist(sig::Array, args...)
#     return spc.compute_spectral_hist(sig, args...)
# end

# function compute_spectral_hist(sig::KeyedArray, args...)
#     (freq, bins, hist) = compute_spectral_hist(Array(sig), args...)
#     return KeyedArray(hist; bin=bins, freq=freq)
# end 




