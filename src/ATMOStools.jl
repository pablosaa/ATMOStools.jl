module ATMOStools

using Printf
using Statistics

"""
ATMOStools a set of functions and constants useful for atmospheric physics and meteorology.

"""
ATMOStools

include("constants.jl")
include("boundary_layer.jl")
include("thermodynamic.jl")
include("radiation.jl")

# For module CLIMA:
Base.include(ATMOStools, "climatology.jl")

#= ******************************************************************

=#
"""
Function to convert Wind U and V component into  wind speed, and direction
USAGE:
```julia-repl
julia> WS, WD = wind_UV2speeddir(U, V)
```
WHERE:
U, V : Wind components U and V [m/s]
WS   : Wind speed [m/s]
WD   : Wind direction [deg]

Input variables can be Float, Vector{Float}, or Matrix{Float}

 (c) 2019, Pablo Saavedra G.
 Geophysical Institute, University of Bergen
"""
function wind_speed(U::T, V::T) where T<:AbstractFloat
    return sqrt(U^2 + V^2)
end
function wind_direction(U::T, V::T) where T<:AbstractFloat
    # according to http://colaweb.gmu.edu/dev/clim301/lectures/wind/wind-uv
    return mod(270 - atand(V, U), 360)
end
function wind_UV2speeddir(U::T, V::T) where T<:AbstractFloat

    WS = wind_speed(U, V)
    WD = wind_direction(U, V)
    return WS, WD
end
# -- OR
function wind_UV2speeddir(U::T, V::T) where T<:Matrix{AbstractFloat}
    
    WS = wind_speed.(U, V)
    WD = wind_direction.(U, V)
    return WS, WD
end
# -- OR
function wind_UV2speeddir(U::T, V::T) where T<:Vector{AbstractFloat}

    WS = wind_speed.(U, V)
    WD = wind_direction.(U, V)
    return WS, WD
end

#= ******************************************************************
 Calculating Integrated Water Vapour Transport
=#
"""
Function to compute Integrated Water Vapour Transport

IVT = calculate_IVT(Pa::Vector, Qv::Vector, WS::Vector; Pmax=300)
IVT = calculate_IVT(Pa::Matrix, Qv::Matrix, WS::Matrix)
IVT, IVT_U, IVT_V = calculate_IVT(Pa::T, Qv::T, (U=Wu::T, V=Wv::T)) where T<:Vector
IVT = calculate_IVT(rs::Dict)

INPUT:
* -> Pa : Pressure [hPa]
* -> WS : Wind speed [m/s] or (U,V) components [m/s]
* -> QV : specific humidity [g/g]
* -> Pmax: Optional parameter for integration limit (Default 300hPa)
or

* -> rs::Dict : Dictionary with Radiosonde or Model data

OUTPUT:
* <- IVT    : total Integrated Water Vapour Transport [kg/m/s]
* <- IVT_vec: IVT by wind component separated (meridional, zonal)
"""
function calculate_IVT(Pa::T, Qv::T, WS::T; Pmax=300) where T<:Vector
    ΔP = diff(Pa, dims=1) #Pa[2:end,:] - Pa[1:end-1,:]
    ΔP *= 1f2 # [hPa -> Pa]
    ϕ = @. Qv*WS
    # finding z_top index for integration limit:
    zₜ = min(findlast(≥(Pmax), Pa), size(ΔP, 1) )

    IVT = sum(ϕ[1:zₜ].* ΔP[1:zₜ])
    IVT /= -g₀
    
    return IVT
end
function calculate_IVT(Pa::T, Qv::T, WS::NamedTuple{(:U, :V), <:Tuple{T,T}}; Pmax=300) where T<:Vector
    
    IVT_u = calculate_IVT(Pa, Qv, WS.U, Pmax=Pmax)
    IVT_v = calculate_IVT(Pa, Qv, WS.V, Pmax=Pmax)

    # calculating total IVT:
    IVT = @. sqrt(IVT_u^2 + IVT_v^2)
    return IVT, IVT_u, IVT_v
end
function calculate_IVT(Pa::T, Qv::T, WS::T; Pmax=300) where T<:Matrix

    IVT = [calculate_IVT(Pa[:,i], Qv[:,i], WS[:,i], Pmax=Pmax) for i ∈ axes(Pa,2)]
    #IVT = reduce(hcat, IVT)
    return IVT
end
function calculate_IVT(rs::Dict; Pmax=300)
    !mapreduce(x->haskey(rs, x), &, [:Pa, :qv]) && error("Input missing keys :Pa or :qv")
    
    Wspd = if haskey(rs, :WSPD)
        rs[:WSPD]
    elseif  haskey(rs, :U) && haskey(rs, :V)
        (U=rs[:U], V=rs[:V])
    else
        error("Input Dictionary needs wind data as :WS or :U and :V variables")
    end

    IVT = calculate_IVT(10rs[:Pa], rs[:qv], Wspd, Pmax=Pmax)
    return IVT
end
# -----/


#= ********************************************************************
    Function to calculate deltaIVT/deltaPa
=#
"""
Function to return the derivative of IVT
given that
IVT = -1/g₀∫Qv*Ws*dP

> δIVT = calculate_∇WVT(P, WS, Qv, H)
returns:
  δIVT = -∇f/g₀

where:
* ∇f = ∂Φᵥ/∂z
* ∂Φᵥ = (Qv*Ws)dP

output:
 δIVT [g m⁻² s-¹]

Note:
P pressure in [hPa], WS in [m s⁻¹], Qv [g/g] and H [km]

"""
function calculate_∇WVT(Pa::T, WS::T, Qv::T, H::Vector) where T<:Vector
    ΔP = 1f2diff(Pa)  # [Pa]
    
    ∇Φ_qv = let Φᵥ = Qv.*WS
    	Δz = 1f3diff(H) # [m]
	∂zΦᵥ = (diff(Φᵥ)./Δz).*ΔP .+ Φᵥ[1:end-1,:].*(ΔP./Δz)
        -∂zΦᵥ/g₀
    end
	
    ∇Φ_qv *= 1f3   # [g m⁻² s-¹]
    return ∇Φ_qv
end
function calculate_∇WVT(Pa::Matrix, WS::Matrix, Qv::Matrix, H::Vector)
    WVT=[calculate_∇WVT(Pa[:,i], WS[:,i], Qv[:,i], H) for i ∈ axes(Pa,2)]
    WVT = reduce(hcat, WVT)
    return WVT
end
function calculate_∇WVT(rs::Dict)
    return calculate_∇WVT(10rs[:Pa], rs[:WSPD], rs[:qv], rs[:height])
end
# ----/

#= *********************************************************************
   Function to estimate the altitude of the maximum water vapour flux
   profile
=#
"""
Function to estimate the altitude of the maximum water vapour flux.
The estimation is based on the gradient of water vapour transport and is
by default considered only the part of the atmospheric profile below 600hPa (top).

USAGE:
> H = estimate_WVT_peak_altitude(rs::Dict)
> H, idx_max, idx_top = estimate_WVT_peak_altitude(rs::Dict, Hmax=8f3, get_index=true)
> H, maxWVT, snrWVT = estimate_WVT_peak_altitude(rs::Dict, WVT=∇WVT, get_maxsnr=true)
> H, maxWVT, snrWVT, idx_max, idx_top = estimate_WVT_peak_altitude(rs::Dict, Hmax=8f3, get_index=true, get_maxsnr=true)

WHERE:
* rs::Dict variable obtained from rs = ATMOStools.getSondeData(sonde_file),
Optional inputs:
* WVT::Matrix the vertical gradient of water vapour transport (default compute from rs)
* Hmax::Real or Hmax::Vector, top height [m], default height Hmax at 600 hPa,
* get_index::Bool flag to output the indexes for max WVT and Hmax (default false).
* get_maxsnr::Bool flag to output the values for max WVT and SNR=maxWVT/median(WVT) (default false).

OUTPUT:
* H::Vector with altitudes of maximum WVT flux below Hmax (same units as rs[:height])
* maxWVT::Vector containing the value of maximum WVT
* snrWVT::Vector containing the ration maxWVT/median(WVT), the median is computed discarting maxWVT 
* idx_max::Vector indexes of rs[:height] to produce H
* idx_top::Vector indexes of rs corresponding to Hmax::Vector


"""
function estimate_WVT_peak_altitude(rs::Dict; WVT=nothing, Hmax=nothing, Hmin=nothing, get_index=false, get_maxsnr=false)
    
    # checking input variables:

    function get_Hxx(Hi, Hdefault)
	Hout = if isnothing(Hi)
	    fill(Hdefault, length(rs[:time]))
	elseif typeof(Hi)<:Real
	    fill(Hi, length(rs[:time]))
	elseif typeof(Hi)<:Vector
	    length(Hi) != length(rs[:time]) && @warn "When Vector, Hmin/Hmax needsto be same length as rs[:time]"
	    Hi
	else
	    @warn "When given Hmin or Hmax, needs to be scalar or Vector"
	end
    end
	
    Hmin = get_Hxx(Hmin, 0)
    
    Hmax = get_Hxx(Hmax, 8000)

    ii_Pmax = [findlast(≤(H_in), rs[:height]) for H_in ∈ Hmax]

    ii_Pmin = [findfirst(>(H_in), rs[:height]) for H_in ∈ Hmin]
    
    # if not given as input, calculate the gradient of water vapour transport:
    ∇ₕWVT = isnothing(WVT) ? calculate_∇WVT(rs) : WVT

    # getting the lowest WVT below Hmax as starting point to search for maximum:
    tt = [let imin = ii_Pmin[j]
	      tmp = ∇ₕWVT[imin+1:imax, j] .-∇ₕWVT[imin:imax-1, j]
	      findlast(≤(-0.1), tmp) |> x-> (isnothing(x) || x<imin) ? imax : x
	  end	for (j, imax) in enumerate(ii_Pmax)]
    
    # estimate the index at which qv*WS get maximum below Pmax:    
    dummy = [let profx = fx[i0:k]
		 profx[isnan.(profx)] .= -9999
		 maxii = try
		     isnothing(k) ? 0 : findmax(profx)
		 catch
		     @error "$(i0) ; $k ; $(profx)"
		 end
                 snrii = filter(<(maxii[1]), profx) |> X-> isempty(X) ? NaN32 : median(X)
                 (maxii..., maxii[1]/snrii )
             end
	     for (fx, i0, k) ∈ zip(eachcol(∇ₕWVT), ii_Pmin, tt)]

    ii_wvt_max = getindex.(dummy, 2)
    maxWVT = getindex.(dummy, 1)
    snrWVT = getindex.(dummy, 3)
                 
    ii_wvt_max .+= ii_Pmin .-1

    H_wvt = [k>0 ? rs[:height][k] : NaN for k ∈ ii_wvt_max]

    ii_Pmax .+= ii_Pmin .-1
    
    ###ii_wvt_max = [isnothing(k) ? 0 : filter(!isnan, fx[i0:k]) |> argmax for (fx, k) ∈ zip(eachcol(∇ₕWVT), ii_Pmax)]
    ###ii_wvt_max .+= i0

    if get_index && get_maxsnr
        return H_wvt, maxWVT, snrWVT, ii_wvt_max, ii_Pmax
    elseif get_index && !get_maxsnr
        return H_wvt, ii_wvt_max, ii_Pmax
    elseif !get_index && get_maxsnr
        return H_wvt, maxWVT, snrWVT
    else
        return H_wvt
    end
        
end
#function estimate_WVT_peak_altitude(rs::Dict; WVT=nothing, Pmax= 600, get_index=true)

    ##∇ₕWVT = let δH = diff(1f3rs[:height])
    ##    ATMOStools.get_δIVT(rs)./δH
    ##end

    ##(old stuff) ii_IWV50 = let ρ_wv = ATMOStools.MassRatio2MassVolume.(rs[:qv], T_K, P_hPa)
    ##    δH = diff([0; rs[:height]])
    ##    [filter(!isnan, Z_wv.*δH) |> WV->cumsum(WV)./sum(WV) |> x->findfirst(≥(0.5), x ) for Z_wv ∈ eachcol(ρ_wv)]
    ##end

    #(old stuff) ii_wvt_max = [isnothing(k) ? 0 : filter(!isnan, fx[findall(Pa.≥ max(Pmax, Pa[k]))]) |> argmax for (fx,Pa,k) ∈ zip(eachcol(∇ₕWVT), eachcol(P_hPa), ii_IWV50)];

    # estimate the index at which qv*WS get maximum below Pmax:
    # if optional input variable "skipsurface=true" then only considers profiles from second layer:

# ----/


#= **********************************************************************
   Calculating gradient of atmospheric profile
=#
"""
Function to calculate the first gradient of a given variable respect of height
δT = ∇ₕT(T::AbstractMatrix, h::Vector)

INPUT:
* -> h : profile altitudes []
* -> T : Profile variable to calculate gradient

OUTPUT:
* <- δT : Gradient of variabe T with respect of variable h

"""
function ∇ₕT(T::AbstractMatrix, h::Vector)
    δH = diff(h)
    δT = diff(T, dims=1)

    return δT./δH
end
# ----/


#= *********************************************************************
   Capturing the inversion layers in a profile variable
=#
"""
Algorithm to estiamte the profile indexes for bottom and top of inversion
layers within the atmosphere.
A maximum of 4 inversion layers is detected by the algorithm.

BASIC USAGE:
> idx_bottom, idx_top = estimate_inversion_layers(T::AbstractMatrix, H::Vector)

will return a matrix of indices containing the positon of up to 4 inversion layers
in the input profile T, based on the altitude H in km

Optional parameters are:
* mxhg = 7, which indicates the maximum altitute to consider in km (default 7 km).
* δH = 0.08, sub-inversion-layers separated by less than δH km will be merged.
* δT = 0.5, negative/positive lapse-rate below δT are ignored and considered on layer.
* mxly = 4, maximum inversion layers to consider (default 4 layers)

"""
function estimate_inversion_layers(T::AbstractMatrix, H::Vector; mxhg=15, δH=0.08, δT=0.5, mxly=4)

    m, nt = size(T)
    mxidx = findfirst(x-> x ≥ mxhg, H)

    out_idx_bot = Matrix{Int32}(undef, mxly, nt) .= 0
    out_idx_top = Matrix{Int32}(undef, mxly, nt) .= 0
    
    # calculating the gradient respect to height h
    δTz = ∇ₕT(T, H)

    for tline ∈ (1:nt)

        all(isnan.(δTz[:, tline])) && continue
        
        idx_bot = Vector{Int32}()
        idx_top = Vector{Int32}()

        # defining dummy expression:
        ex =:(δTz < 0)
        
        for i ∈ (1:mxidx)

	    ex.args[2] = δTz[i, tline]
		
	    if !isnan(δTz[i, tline]) && !eval(ex)
            
	        if ex.args[1] == :<
		    push!(idx_bot, i) #i-1
		    ex.args[1] = :≥
                        
	        else
		    push!(idx_top, i)
		    ex.args[1] = :<
                end
	    end
        end

        # Checking if it failed to find top layer:
        ex.args[1] == :≥ && pop!(idx_bot)
        
        # Merging layers closer than 60m form top to bottom
        δH_tb = H[idx_bot[2:end]] .- H[idx_top[1:end-1]]
        idx_out = δH_tb |> x-> (x .≤ δH) |> findall
        deleteat!(idx_bot, idx_out.+1)
        deleteat!(idx_top, idx_out)

        # Merging consecutive layers with difference less than threshold:
        δTinv = T[idx_top[1:end-1], tline] .- T[idx_bot[2:end], tline]
	idx_inv = findall(x-> x ≤ δT, δTinv)	
	deleteat!(idx_bot, idx_inv.+1)
	deleteat!(idx_top, idx_inv)
        
        # Dismissing layers with T inversion strength < 0.5 °C
        δTinv = T[idx_top, tline] .- T[idx_bot, tline]
        idx_inv = findall(x-> x ≤ δT, δTinv)
        deleteat!(idx_bot, idx_inv)
        deleteat!(idx_top, idx_inv)

        # Returning the allowed layers:
        ninv = min(length(idx_bot), mxly)
        out_idx_bot[1:ninv, tline] = idx_bot[1:ninv]
        out_idx_top[1:ninv, tline] = idx_top[1:ninv]
    end
    
    return out_idx_bot, out_idx_top
end
# ----/

#= **********************************************************************
   extracting the inversion variables from the specified profiles
=#
"""
function to extract from atmospheric profile variables the bottom and
 the thickness of the inversion give the profile indeces obtained by
the ATMOStools.jl function estimate_inversion_layers()

Default atmospheric variables are: T, Pa and height
USAGE:

out::Dict = get_inversion_variables(idx_bottom, idx_top, RS::Dict)

or for only one variable, e.g. temperature:

out::Dict = get_inversion_variables(idx_bottom, idx_top, RS::Dict, vars=(:T, ))

or for T, height and Qv:

out::Dict = get_inversion_variables(idx_bottom, idx_top, RS::Dict, vars=(:T, :height, :qv))

note that the Symbols in vars=() need to match with the variables in RS:Dict().

The output is a dictionary with the same variables in capital and added the thickness
with a sufix Δ, e.g. for temperature the output Dict will have two related variables:
:T and :ΔT for temperature bottom and strength, respectively.

"""
function get_inversion_variables(idx_bot, idx_top, rs::Dict; addvars=())

    # Default variables to use:
    vars = (:T, :height, Symbol.(addvars)...)
    nv, nt = size(idx_bot)

            
    # creating a set of variables to map the inversion variables:
    # e.g. :Pa will be mapped to :PA for bottom and top inversion and
    # to :ΔPA for strength
    varset = Dict( x =>
                   let tmp = uppercase(String(x))
                   (Symbol(tmp), Symbol(:Δ, tmp))
                   end
                   for x in vars if haskey(rs,x))

    # Creating output dictionary:
    out = Dict()
    
    # filling output variable with inversion estimations:	
    foreach(vars) do T
	out[varset[T][1]] = Matrix{Float32}(undef, nv, nt) .= NaN
	out[varset[T][2]] = Matrix{Float32}(undef, nv, nt) .= NaN
		
	for i ∈ (1:nv)
	    
	    for j ∈ findall(idx_bot[i,:] .> 0)
		K₀ = idx_bot[i, j]
		K₁ = idx_top[i, j]
		if typeof(rs[T]) <: Vector
                    # variables with bottom value: 
		    out[varset[T][1]][i, j] = rs[T][K₀]

                    # variables with strength value:
		    out[varset[T][2]][i, j] = rs[T][K₁] - rs[T][K₀]
		else
		    # variables with bottom value: 
		    out[varset[T][1]][i, j] = rs[T][K₀, j]

                    # variables with strength value:
		    out[varset[T][2]][i, j] = rs[T][K₁, j] - rs[T][K₀, j]
		end
	    end
	end
    end
    
    return out
   
end
# ----/

# ********************************************************************
# ********************************************************************
#  Old version for profile IVT:
# ********************************************************************
function getIVT(Pa::Matrix, U::Matrix, V::Matrix, QV::Matrix)
    
    ΔP = Pa[2:end,:] - Pa[1:end-1,:]
    ΔP *= 1f2 # [Pa]
    tmp_u = QV .* U
    tmp_v = QV .* V
    z_top = 235  # index of the top layer to consider ≈ 306 hPa.
    # calculating IVT components u and v:
    IVT_u = map(1:z_top) do z
        sum(tmp_u[z:z_top,:] .* ΔP[z:z_top,:], dims=1)
    end
    IVT_v = map(1:z_top) do z
        sum(tmp_v[z:z_top,:] .* ΔP[z:z_top,:], dims=1)
    end
    
    IVT_u = vcat(IVT_u...)
    IVT_v = vcat(IVT_v...)
    IVT_vec = cat(IVT_u, IVT_v, dims=3)
    IVT_vec /= g₀
    # calculating total IVT:
    IVT = sqrt.(IVT_u.^2 + IVT_v.^2)
    IVT /= g₀
    
    return IVT, IVT_vec
end
function getIVT(rs::Dict)
    IVT, IVT_vec = getIVT(10f0*rs[:Pa], rs[:U], rs[:V], rs[:qv])
    return IVT, IVT_vec
end
# ----/

#= ******************************************************************
   Get Maximum IVT index and values
=#
"""
Function to obtain the altitude of the maximum IVT from profile
idx, IVTmax = getMaxIVT_θ(IVT::Matrix{Float64})

INPUT:
* IVT::Matrix{Float64} -> 2D array with cumulative total IVT
OUTPUT:
* idx -> 1D Vector with the index of the maximum from IVT Matrix
* IVTmax -> 1D Vector with maximum IVT
"""
function getMaxIVT_θ(IVT::Matrix)
    
    ii = argmax(IVT, dims=1)
    # converting the CartesianIndexes ii from 2-D to 1-D :
    ii = vcat(ii...)

    # retrieving the max of IVT into a vector:
    IVTmax = IVT[ii][:]
    
    return ii, IVTmax
end
# ---


# --
end #module

