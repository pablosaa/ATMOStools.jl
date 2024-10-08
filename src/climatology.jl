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
using HTTP
using StatsBase

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
* ```iFiels Symbol or String``` indicate which field of the JSON file to read (default :items),
* For ```.dat``` files with multiple climate indexes, ```iFields``` can change the default to other column name,
    e.g. ENSO SST data can be ```iFields="NINO3"``` to rename the Nino3's anomaly to ```:idx```.

Reading Arctic Oscilation index from:
```html
* "https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.ao.cdas.z1000.19500101_current.csv"
```
Reading ENSO and PDO indexes from JSON files provided by:
```html
* "https://sealevel.jpl.nasa.gov/data/vital-signs/pacific-decadal-oscillation"
```
 
Reading NOAA data for ENSO and PDO indexes as .dat files obtained from:
```html
* "https://www.ncei.noaa.gov/access/monitoring/products/"
```
* * e.g. for ENSO 3-month seasonal Oceanic Nino Index (ONI), or monthly Nino3.4 SST and Anomalies.
* * for PDO monthly index. 
"""
function load_climate_index(datain::String; Tlim::Tuple{T, T}=(DateTime(1000,1,1),now()), iField=nothing) where T<:DateTime

    # Checking if the file is present in local dir or can be downloaded:
    fileext = split(datain, '.')[end] |> lowercase
    
    filen = if isfile(datain)
        datain
    else
        @info "Trying to download the file"
        http_response = HTTP.get(datain)
        http_response.body
    end
    
    
    dfindex = if fileext=="json"
        @info "Reading a JSON file"

        isnothing(iField) && (iField = :items)
        tmp = JSON3.read(filen)[iField] |> DataFrame
        select!(tmp, [:x, :y] .=> ByRow(f->parse(Float32,f)) .=> [:date, :idx] )
        transform!(tmp, :date => ByRow(PartialYear2DateTime), renamecols=false)
        tmp
       
    elseif fileext=="csv"
        @info "Reading a CSV file"

        isnothing(iField) && (iField = :ao_index_cdas)
        
        tmp = CSV.read(filen, header=1, ntasks=1, types=Dict(iField=>Float32), DataFrame)
        transform!(tmp, [:year, :month, :day] => ByRow(DateTime) => :date)
        filter!(d->!ismissing(d.ao_index_cdas), tmp)
        transform!(tmp, iField => f->collect(skipmissing(f)), renamecols=false)
        select!(tmp, [:date, iField] .=> [:date, :idx])
        tmp
    elseif fileext=="dat"
        @info "Reading a DAT file"
        
        tmp = open(filen, "r") do IO
            kopf = readline(IO)
            # checking what kind of file is:
            dat = if contains(kopf, "PDO Index")
                dat = CSV.read(IO, DataFrame; missingstring="99.99", delim=' ', ignorerepeated=true)
            
                # stacking the 12 months into a column:
                dat = sort(stack(dat, variable_name=:month, value_name=:idx), :Year)
                # converting Year and Month columns into DateTime format with column name :date
                select(dat, [:Year, :month] => ByRow((y,m)->DateTime("$(y)-"*m, "yyyy-u")) => :date, :idx)

            elseif contains(kopf, "NINO")
                isnothing(iField) && (iField = "NINO3.4")
                header_names = split(kopf)
                iidx = findfirst(==(iField), header_names) + 1
                header_names[iidx] = "idx"

                dat = CSV.read(IO, delim=' ', ignorerepeated=true, DataFrame, header=Symbol.(header_names) )
                select(dat, [:YR, :MON] => ByRow((y,m)->DateTime(y,m)) => :date, names(dat)[3:end]...)

            elseif contains(kopf, "SEAS")
                isnothing(iField) && (iField = "ANOM")
                header_names = split(kopf)
                iidx = findfirst(==(iField), header_names)
                header_names[iidx] = "idx"

                dat = CSV.read(IO, delim=' ', ignorerepeated=true, DataFrame, header=Symbol.(header_names) )
                # defining the seasons as numeric from 1 to 12:
                monaten = Dict(ss=>i for (i,ss) in enumerate(["DJF","JFM","FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ"]))
                transform(dat, [:SEAS, :YR] => ByRow((x,y)->DateTime(y,monaten[x]))=>:date, renamecols=false)
            else
                @warn "DAT file not supported!"
                none
            end

        end
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
julia> DF = FourierFrequencies(Yt; fₛ=5)
julia> DF = FourierFrequencies(Ts, Yt)
julia> DF = FourierFrequencies(Ts, Yt; P=Week, tₕ=3.1, fullout=true)
``` 
Input:
* ```Ts::Vector{DateTime}``` DateTime corresponding for time series samples,
* ```Yt::AbstractArray``` sampled data (signal) corresponding to ```Ts```,

Optional arguments are:
* ```P::Type{<:Period}``` period type for output frequency νₖ vector (default ```Dates.Day```),
* ```tₕ::Real``` threshold to limit noise in phase calculation (default ```nothing``` which uses ```maximum(Yₖ)/10```),
* ```fullout::Bool``` if true, outputs the imaginary power spectrum including 1st element DC value (default ```false```)
* ```fₛ::Real``` sampling rate (default 1, if ```Ts``` is given, ```fₛ``` is estimated based on ```Ts```)
* ```norm_unit::Bool``` if ```true``` and ```Ts``` is given, then the power stectrum is normalized to the ```Ts``` resolution, otherwise is normalized to ```P``` units.

Output ```DF::DataFrame``` with the column names:
* ```νₖ::Vector{Float64}``` Frequency in units ```1/P``` (default ```P=day```),
* ```Pₖ::Vector{Float64}``` Period as ```νₖ⁻¹```,
* ```Yₖ::Vector{Float64}``` or ```Vector{Complex{Float64}}``` amplitude of FFT,
* ```ϕₖ::Vector{Float64}``` phase of FFT in units ```deg```,
* ```YdBₖ::Vector{Float64}``` power spectrum in units ```[dB P]``` (default ```P=day```)

"""
function FourierFrequencies(yt::AbstractArray; fₛ=1, tₕ=nothing, fullout=false)

    # Avoiding missing data:
    yt = (collect∘skipmissing)(yt)
    
    # Number of samples in signal:
    N = length(yt)

    # estimating the half of FFT vector:
    N₂ = round(Int32, N/2)
    iseven(N) && (N₂ += 1) 

    # Calculating |FFT|² and the phase ϕₖ of the signal:
    𝑌 = fft(yt) |> Y->Y[1:N₂]
    yfft, ϕ = let y= @. abs2(𝑌)   # abs2(Y) == abs(Y)²
        y ./= (N*fₛ)
        y[2:end-1] .*= 2
        
        tₕ = ifelse(isnothing(tₕ), maximum(y)/10, tₕ)

        phi = @. ifelse(y<tₕ, missing, atand(imag(𝑌)/real(𝑌)) )
        y, phi
    end
    
    fullout && (yfft = 𝑌)  # output as FFT complex vector

    # calculating the discrete frequency for k-bin:
    νₖ = (0:N-1) |> k-> k/N*fₛ
    
    # Convert output to decibels:
    ydB = @. 10log10(yfft)

    df = DataFrame(νₖ=νₖ[1:N₂], Yₖ=yfft, ϕₖ=ϕ, YdBₖ=ydB)
    
    return ifelse(fullout, df, df[2:end, :])
end
# or: 
function FourierFrequencies(T::Vector{DateTime}, yt::AbstractArray; P::Type{<:Period}=Day, tₕ=nothing, fullout=false, norm_unit=true)

    N = length(T)
    ΔT = extrema(T) |> t->t[2]-t[1]  # [Millisecoonds]
    Ft = typeof(ΔT)

    # converting the ΔT from sampling units (e.g. days) to period given by P (e.g. Year) :
    ΔT /= (Ft∘Dates.toms)(P(1))  # type of ΔT is Float64 and the unit is [P] (e.g Year)

    # sampling rate [#/P]
    fₛ = (N-1)/ΔT

    # Calculating |FFT|²
    df = FourierFrequencies(yt; fₛ=fₛ, tₕ=tₕ, fullout=fullout)

    # adding the period as νₖ⁻¹
    insertcols!(df, 2, :Pₖ => inv.(df.νₖ))
        
    # if option given, convert normalization factor 1/fₛ from P to ΔT units:
    if norm_unit
        # the median resolutin of time vector T [milliseconds]
        δt=diff(T) .|> Dates.toms |> median
        δt /= Dates.toms(P(1))
        df.Yₖ ./= δt
        df.YdBₖ .-= 10log10(δt);
    end

    # returning df
    return df
    
end
# ----/

# ==========
"""
===================================================================
 Function to apply moving window average to minimize seasonality:
```julia-repl
julia> Ysm = ave_window(Yin; w=5, 𝐹ave=median)
julia> Ysm = ave_window(Yin; w=[0.1, 0.4, 1, 0.6, 0.2], 𝐹ave=median)
```
Where:
* ```Yin::Vector``` data to be averaged,
* ```w::Number``` the width of the window (default=6) to moving average, or
* ```w::Vector``` the weigths used for the moving average applied to the middle point,
* ```𝐹ave::Function``` to use for the evarage (default ```mean```),
* ```sp::Number``` step to consider for the vector ```Yin```, (default=1 i.e. every element),
* ```edges::Bool``` in case ```w::Number``` wheter or not include the edges not covered by a full window w (default ```true```).
Output:
* ```Ysm::Vector``` the same length as ```Yin``` but smoothed by moving average ```𝐹ave``` function.
"""
function ave_window(y::Vector{<:AbstractFloat}; w=6, 𝐹ave::Function=mean, sp=1, edges=true)
    n = length(y)
    y_idx = range(sp, step=sp, stop=n)

    # defining weights length centered at data point:
    if typeof(w)<:Union{Number, Vector} && ndims(w)<2
        nw = length(w)
        iseven(nw) && @error("when w::Vector, the length(w) needs to be odd.")
    else
        @error("window parameter w must be a scalar or Vector!")
    end

    δw = ifelse(nw==1, w, nw) |> x-> floor(Int8, x/2)

    # defining output vector:
    TT = eltype(y) 
    y_ave = TT.(fill(NaN, length(y_idx)) )
    
    for (j, i) in enumerate(y_idx)
        i0 = if nw==1
            range(i-δw, i+δw)
        else
            range(i-δw, length(nw))
        end

        n0 = length(i0)
        i1 = filter(k->1≤k≤n, i0)
        # if edges=false, then skip indexes ≤ δw or ≥ n-δw:
        (!edges && n0!=length(i1)) && continue
        # definig weights:)
        ω = if nw==1
            @. 1 - abs(i1 - j)/n0
        else
            ij = findall(x->0 <x≤n, i0)
            w[ij]
        end

	# selecting only values of y that are not NaNs:
        inan = findall(!isnan, y[i1])
        # calculating the window average accoring to function 𝐹ave with weights ω:
	y_ave[j] = 𝐹ave(y[i1][inan], weights(ω[inan]))
    end
    
    return y_ave
end
# ----/

end # end of module CLIMA
# ================================================================================
