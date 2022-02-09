using FourierAnalysis

###################################################################

# (1)

# Check that correct amplitude spectra is returned at all discrete
# Fourier Frequency (using a rectangular taper).
# Sinusoids are created at all frequencies with amplitude 10 at the
# first frequency and then incrementing by 10 units along frequencies.
# NB when t is even, correct amplitude for the last frequency is obtained
# only if the sinusoidal has a particular phase.

sr, t, wl= 16, 32, 16
V=Matrix{Float64}(undef, t, wl)
for i=1:wl V[:, i]=sinusoidal(10*i, b2f(i, sr, t), sr, t, π/6) end

using FFTW
# using FFTW.jl only
P=plan_rfft(V, 1)*(2/t);
P
Σ=abs.(P*V)
using Plots
bar(Σ[brange(t, DC=true), :], labels="")

# using FourierAnalysis.jl
Σ2=spectra(V, sr, t; tapering=rectangular, func=√, DC=true)
using Plots
bar(Σ2.y[brange(t, DC=true), :], labels="")

#############################################################################

### Check amplitude spectra on long signals obtained by welch methods
# one sinusoidal is at an exact discrete Fourier Frequency and the other not
# Rectangular window
sr, t, f, a = 128, 128, 10, 0.5
v=sinusoidal(a, f, sr, t*16)+sinusoidal(a, f*3.5+0.5, sr, t*16)+randn(t*16);
Σ=spectra(v, sr, t; tapering=rectangular, func=√)
bar(Σ.y, labels="rectangular")

# harris4 window (default)
Σ2=spectra(v, sr, t; func=√)
bar!(Σ2.y, labels="harris4")

#smooth spectra
Σ3=smooth(blackmanSmoother, Σ2)
bar!(Σ3.y, labels="smoothed")

#############################################################################

function generateSomeData(sr::Int, t::Int; noise::Real=1.)
    # four sinusoids of length t samples and sr sampling rate
    # peak amplitude: 0.7, 0.6, 0.5, 0.4
    # frequency:        5,   7,  13,  27
    # phase:            0, π/4, π/2,   π
    v1=sinusoidal(0.7, 5,  sr, t, 0)
    v2=sinusoidal(0.6, 7,  sr, t, π/4)
    v3=sinusoidal(0.5, 13, sr, t, π/2)
    v4=sinusoidal(0.4, 27, sr, t, π)
    return hcat(v1, v2, v3, v4) + (randn(t, 4)*noise)
end

sr, wl, t = 128, 512, 8192
X=generateSomeData(sr, t)
# multivariate data matrix 8192x4

# compute spectra
S=spectra(X, sr, wl)

# check the spectrum of first series
S.y[:, 1]

# gather some plot attributes to get nice plots
using Plots.Measures
spectraArgs=(fmax = 32,
             left_margin = 2mm,
             bottom_margin = 2mm,
             xtickfont = font(11, "Times"),
             ytickfont = font(11, "Times"))
plot(S; spectraArgs...)
plot(S; xspace=2, spectraArgs...)

# use a custom simple taperig window
S=spectra(X, sr, wl; tapering=riesz)

# use Slepian's multi-tapering
S=spectra(X, sr, wl; tapering=slepians(sr, wl, 1.5))

# construct with smoothing
S=spectra(X, sr, wl; tapering=slepians(sr, wl, 1.5), smoothing=hannSmoother)

# compute Amplitude Spectra instead
S=spectra(X, sr, wl; tapering=slepians(sr, wl, 1.5), func=√)

# plot Aplitude spectra
plot(S; ytitle="Amplitude", spectraArgs...)

# smooth the spectra a-posteriori
S=smooth(blackmanSmoother, S)

# extract spectra in range (8Hz to 12Hz)
e=extract(S, (8, 12))

# extract spectra in range (8Hz to 12Hz) for series 1 and 2
e=extract(S, (8, 12))[:, 1:2]

# extract the spectra at 10Hz only (1 bin per series)
e=extract(S, 10)

# average spectra in the 8Hz-12Hz range
bar(mean(S, (8.0, 12.0)))

# average across series of the average spectra in the 8Hz-12Hz range
mean(mean(S, (8.0, 12.0)))

# average spectrum across all frequencies for each series
bar(mean(S, :))

# average spectra in equally-spaced 2Hz-band-pass regions for all series
plot(bands(S, 2))

# average spectra in equally-spaced 2Hz-band-pass regions for series 1 and 2
plot(bands(S, 2)[:, 1:2])

# (2)

# generate 3 multivariate data matrices 8192x4
X=[generateSomeData(sr, t) for i=1:3]

# Now the call to the spectra function will generate 3 Spectra objects
S=spectra(X, sr, wl)
plot(S[1]; spectraArgs...)
plot(S[2]; spectraArgs...)
plot(S[3]; spectraArgs...)

# when you want to compute the spectra of many data matrices you may want
# to do it using a fast FFTW plan (wait up to 10s for computing the plan)
plan=Planner(plan_exhaustive, 10.0, wl)
S=spectra(X, sr, wl; planner=plan)

# how faster is this?
using BenchmarkTools
@benchmark(spectra(X, sr, wl))
@benchmark(spectra(X, sr, wl; planner=plan))


# average spectra in range (8Hz-12Hz) for all series of all objects
M=mean(S, (8, 12))

# plot the average spectrum across all series for the three objects
# using Julia standard mean function
plot(mean(S[1].y[:, i] for i=1:size(S[1].y, 2)))
plot!(mean(S[2].y[:, i] for i=1:size(S[2].y, 2)))
plot!(mean(S[3].y[:, i] for i=1:size(S[3].y, 2)))

# average spectra in range (4Hz-32.5Hz) across objects for the 4 series
plot(mean(mean(S, (4, 32.5))))

# extract spectra in range (8Hz-12Hz) for all series and all subjects
extract(S, (8, 12))

# if you enter en illegal range, nothing will be done and you will get
# an error in the REPL explaining what is wrong in your argument
extract(S, (0, 12))
mean(S, (1, 128))

# extract 4Hz-band-pass average spectra for all series and all objects
bands(S, 4)

# Apply smoothing in the spectra computations
S=spectra(X, sr, t; smoothing=blackmanSmoother)
plot(S[1]; spectraArgs...)
plot(S[2]; spectraArgs...)
plot(S[3]; spectraArgs...)

# plot spectra in in 1Hz band-pass regions for all series in S[1]
plot(bands(S[1], 1))

# use slepian multi-tapering
S=spectra(X, sr, wl; tapering=slepians(sr, wl, 1.))
plot(S[1]; spectraArgs...)

# average spectra across objects
plot(mean(s.y for s ∈ S))

