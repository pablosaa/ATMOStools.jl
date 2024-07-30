# Part of ATMOStools.jl
# Contains function related to climatology
#
# See LICENSE.TXT
# ================================================================================

module CLIMA

using Dates
using JSON3
using FFTW
using CSV, DataFrames


"""
Function to convert date from fraction of a year to DateTime format.
USAGE:
```julia-repl
julia> Datum = PartialYear2DateTime(FracYear)
julia> Datum = PartialYear2DateTime(FracYear; period=Millisecond)
```
WHERE:
* ```FracYear::AbstractFloat``` is the year in fractional format e.g. 2004.4321
* ```period::Type(<:Period``` is the precision used for the conversion, optional variable, default Minute
* ```Datum::DateTime``` the converted date in julia format DateTime

"""
function PartialYear2DateTime(jahr::AbstractFloat; period::Type{<:Period}=Minute)
    year0, δ = divrem(jahr,1)
    Δ = reduce(-, DateTime.(year0 .+ [1, 0])) |> period |> Dates.value
    Δ *= δ
    return DateTime(year0) + (period∘round)(Δ)
end
# ----

#=
Reading JSON data file for ENSO index and converting it to DataFrame
=#
"""
Function to read ENSO and PDO indexes:

```julia-repl
julia> DF = load_climate_index("/data/enso_index.json")
julia> DF = load_climate_index("/data/enso_index.json"; Tlim=[Date(2012,11), today()])
``` 
Output ```DF::DataFrame``` with the column names :date, :idx

Optional arguments are:
* ```Tlim::Tuple{Date, Date}``` contains 2 dates to limit the dataset (default read all data),
* ```iFiels::Symbol``` to indicate which field of the JSON file to read (default :items),

Reading Arctic Oscilation index from:
```html
* "https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.ao.cdas.z1000.19500101_current.csv"
```
Reading ENSO and PDO indexes from:
```html
* "https://sealevel.jpl.nasa.gov/data/vital-signs/pacific-decadal-oscillation"
```
 
"""
function load_climate_index(filen::String; Tlim::Tuple{Date, Date}=(), iField=:items)

    fileext = split(filen, '.')[end] |> lowercase
    
    dfindex = if fileext=="json"
        @info "Reading a JSON file"
        tmp = JSON3.read(filen)[iField] |> DataFrame
        select!(tmp, [:x, :y] .=> ByRow(f->parse(Float32,f)) .=> [:date, :idx] )
        transform!(tmp, :date => ByRow(PartialYear2DateTime), renamecols=false)
        tmp
       
    elseif fileext=="csv"
        @info "Reading a CSV file"
        http_response = HTTP.get(filen)
        tmp = CSV.read(http_response.body, header=1, ntasks=1, types=Dict(:ao_index_cdas=>Float32), DataFrame)
        transform!(tmp, [:year, :month, :day] => ByRow(DateTime) => :date)
        filter!(d->!ismissing(d.ao_index_cdas), tmp)
        transform!(tmp, :ao_index_cdas => f->collect(skipmissing(f)), renamecols=false)
        select!(tmp, [:date, :ao_index_cdas] .=> [:date, :idx])
        tmp
    else
        @warn "Data file name not supported, needs to be CSV or JSON."
        none
    end
    # if optional argument given, select the time frame required:
    if !isempty(Tlim)
        filter!(d->Tlim[1]≤ d.date ≤Tlim[2] , dfindex)
    end 
    return dfindex
end
# ----

#=
Script to estimate the 3 most powered FFT frequencies from time series:
=#
"""
Function to perform a FFT to timeseries vector given the period of data sampling.

```julia-repl
julia> DF = FourierFrequencies(Ts, Yt)
julia> DF = FourierFrequencies(Ts, Yt; P=Week, dB₀=20, fullout=true)
``` 
Input:
* ```Ts::Vector{DateTime}``` vector of time series,
* ```Yt::AbstractArray``` vector of sampled data corresponding to ```Ts```,

Optional arguments are:
* ```P::Type{<:Period}``` period of sampling data (default ```Day```),
* ```dB₀::Float``` minimum threshold to consider in dB (default ```nothing```),
* ```fullout::Bool``` if true, output include 1st element of FFT (default ```false```)


Output ```DF::DataFrame``` with the column names :νₖ, :YdB, :Yfft

"""
function FourierFrequencies(T::Vector{DateTime}, yt::AbstractArray; dB₀=nothing, P::Type{<:Period}=Day, fullout=false)
    N = length(T)
    ΔT = extrema(T) |> t->t[2]-t[1]  # [Millisecoonds]
    Ft = typeof(ΔT)
    ΔT /= (Ft∘Dates.toms)(P(1))  # Ft(P(1))

    # sampling rate [#/P]
    fₛ = (N-1)/ΔT

    # estimating the half of FFT vector:
    N₂ = round(Int32, N/2)

    # Calculating FFT
    yfft = fft(yt)

    # Convert output to decibels:
    ydB = @. 10log10(abs(yfft))
    
    # calculating the frequency bins:
    νₖ = (0:N-1) |> k-> k/N*fₛ
    
    idx_out = let idx=2:N₂
        if !isnothing(dB₀)
            idx = findall(>(dB₀), ydB[idx])
        end
        fullout && (idx .+= 1)  # add 1 to count for the zero frequency
    end
    
    return νₖ[idx_out], yfft[idx_out], ydB[idx_out]
    
end
# ----/

# ==========
"""
===================================================================
 Function to apply running window average to minimize seasonality:


"""
function ave_window(y::Vector{<:AbstractFloat}; w::Number=6, 𝐹ave::Function=mean, sp=1)
    n = length(y)
    y_idx = range(sp, step=sp, stop=n)

    # defining weights length centered at data point:
    δw = round(Int8, w/2)

    # defining output vector:
    y_ave = Vector{eltype(y)}(under, length(y_idx))
    
    for (j, i) in enumerate(y_idx)  #eachindex(y)
	i0 = range(i-δw, i+δw)
	i1 = min.(n, max.(1, i0))
	ω = ones(length(i0))
	ω[i0.<1] .= 0.5
	ω[i0.>n] .= 0.5
	inan = findall(!isnan, y[i1])
        # calculating the window average accoring to function 𝐹ave with weights ω:
	y_ave[j] = 𝐹ave(y[i1][inan], weights(ω[inan]))
    end
    
    return y_ave
end
# ----/

end # end of module CLIMA
# ================================================================================
