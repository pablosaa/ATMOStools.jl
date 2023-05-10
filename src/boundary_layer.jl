# Part of ATMOStools package
# *********************************************
# +++++++ Boundary Layer funcitons ++++++++++++
# *********************************************

# ********************************************************************
# Calculating Richardson Number
#
"""
Function to estimate the Richardson Number
N², Ri = Ri_N(height, WSPD, QV, θ, T)

INPUT: 
* -> height: profile altitudes [m]
* -> U,V   : wind speed componetns U,V [m/s]
* -> QV    : specific humidity [kg/kg]
* -> θ     : potential temperature [K]
* -> T     : ambient temperature [°C]
OR
* -> rs::Dict : Dictionary with Radiosonde or Model data
OPTIONAL PARAMETERS:
* Tₛ -> a vector of surface temperature [K] (default first layer or T)
* Hₛ -> reference altitude [m] (default 0)
* WSPDₛ -> wind speed at the reference altitude [m/s] (default 0)

 OUTPUT:
* <- N² : Brunt-Väisälä frequency [s⁻²]
* <- Ri : The bulk Richardson Number [-]

"""
function Ri_N(height::Vector, U::Matrix, V::Matrix, QV::Matrix, θ::Matrix, T::Matrix;
              Tₛ::Vector=[], Hₛ = 0, Uₛ = 0, Vₛ=0, flag_mixr=false)

    # calculating virtual potential temperatures:
    θᵥ = VirtualTemperature(θ, QV, flag_mixr=flag_mixr)

    if isempty(Tₛ)
        TVₛ_K = T[1,:] .+ 273.15       
        θₛ = θᵥ[1,:]
        Uₛ = 0.0 #U[1,:]
        Vₛ = 0.0 #V[1,:]
    else
        @assert length(Tₛ)==size(T, 2) error("Surface temperature Tₛ length does not match dim 2 of T")
        
        θₛ = VirtualTemperature(θ(Tₛ, P₀), QV[1,:], flag_mixr=flag_mixr)
    end
    
    Δθᵥ = θᵥ .- θₛ'   # θᵥ[1,:]'  # K
    ΔZ = height .- Hₛ  # height[1]  # m
    ΔU = U .- Uₛ' # [1,:]'  # m/s
    ΔV = V .- Vₛ'
    
    # calculating the Brunt-Väisälä frequency:
    N² = @. Δθᵥ/ΔZ/θᵥ
    N² *= g₀
    
    # calculating Richardson-number: N²/wind shear gradient:
    Ri = @. N²/((ΔU/ΔZ)^2 + (ΔV/ΔZ)^2)
    
    return N², Ri
end
# When function called from RS, calculate θ and mixing ratio from radiosonde data:
function Ri_N(rs::Dict)
    θ_data = θ(rs[:T].+273.15, 10rs[:Pa])
    mixr = qx_to_mixr(rs[:qv])
    # by default rs[:height] is read from ARM data in km
    N², Ri = Ri_N(1f3*rs[:height], rs[:U], rs[:V], mixr, θ_data, rs[:T], flag_mixr=true)

    return N², Ri
end
# ----/


"""
Function to estimate the Planetary boundary layer height from bulk Ri

USAGE:
> PBLH = estimatei\\_Ri\\_PBLH(Ri::Matrix, height::Vector)
> PBLH = estimate\\_Ri\\_PBLH(Ri::Matrix, height::Vector; ξ\\_ri=0.11)
> PBLH = estimate\\_Ri\\_PBLH(rs::Dict; ξ\\_ri=0.08)

WHERE:
* Ri::Matrix is the bulk Richardson Number (height, time)
* height::Vector the altitude corresponding to Ri 1st dimension (height)
* rs::Dict Radiosonde data readed by ARMtools
Optional variables:
* ξ\\_ri::Real Threshold of Ri to consider BLH (default 0.25)
* Hmax::Real The maximum height if ξ\\_ri not found (default height[end-1])
* i0::Integer lowest index of height::Vector to start the estimation.

OUTPUT:
* PBLH:Vector Boundary layer height, same units of height or rs[:height]
"""
function estimate_Ri_PBLH(Ri::AbstractMatrix, height::Vector; ξ_ri=0.25, i0=1, Hmax=nothing)
    Hmax = isnothing(Hmax) ? height[end-1] : Hmax
    𝐻(p) = isnothing(p) ? Hmax : height[p+i0-1]
    PBLH = [findfirst(≥(ξ_ri), x) |> 𝐻 for x in eachcol(Ri[i0:end, :])]
    return PBLH
end
# or
function estimate_Ri_PBLH(rs::Dict; ξ_ri=0.25, i0=1, Hmax=nothing)
    N², Ri = ATMOStools.Ri_N(rs);
    return estimate_Ri_PBLH(Ri, rs[:height]; ξ_ri = ξ_ri, i0=i0, Hmax=Hmax)
end
# ----/

# ***************************************************************************
"""
Function to find the height for the mixing layer below cloud base.

USAGE:

CMLH = cloud_decoupling_height(height, CBH, θᵥ)
CMLH = cloud_decoupling_height(height, CBH, θᵥ; θ_thr=0.05)
TCMLH = cloud_decoupling_height(height, CTH, θᵥ; topmixlayer=true)

WHERE:
* height::Vector with the profile altitudes,
* CBH::Vector with the information of the cloud base height to consider,
* θᵥ::Matrix with the potential temperature (or virtual θ) to use for the estimation
OPTIONAL:
* θ_thr=0.05 threshold for the cummulative variance to consider (default 0.025 K²).
* topmixlayer=true then the mixing layer is estimaged above the given cloud top CTH
OUTPUT:
* CMLH::Vector cloud mixing layer height a.g.l., same units as height.
* TCMLH::Vector (if optional parameter topmixlayer=true) cloud top mixing layer height above the cloud

"""
function cloud_decoupling_height(rs_height::Vector, CBH::Vector, θᵥ::Matrix; θ_thr=0.025, topmixlayer=false)
    # aux function to estimate cumulative sum:
    Σₖ(x::Vector) = cumsum(x)./(1:length(x))
    decop_hgt = fill(NaN32, length(CBH))
	
    for (k, θp) in enumerate(eachcol(θᵥ))
	ii_below = if !topmixlayer
            findall(rs_height .≤ CBH[k]) |> reverse
        else
            findall(rs_height .≥ CBH[k])
        end
        
    	var_θv = (Σₖ(θp[ii_below].^2) .- (Σₖ(θp[ii_below])).^2)
    	all(isnan.(var_θv)) && continue
        
    	ii_decop = findlast(var_θv .≤ θ_thr) |> x-> ii_below[x]
        
    	decop_hgt[k] = if !topmixlayer
            ii_decop ≤ 1 ? 0 : rs_height[ii_decop]
        else
            ii_decop ≤ length(rs_height) ? rs_height[ii_decop] : rs_height[end]
        end
    end
    
    return decop_hgt
end
function cloud_decoupling_height(rs_height::Vector, CBH::Matrix, θᵥ::Matrix; θ_thr=0.025, topmixlayer=false)
    ntime, nlayer = size(CBH)
    decop_hgt = Matrix{eltype(rs_height)}(undef, ntime, nlayer)
    for i ∈ (1:nlayer)
        decop_hgt[:, i] = cloud_decoupling_height(rs_height, CBH[:,i], θᵥ, θ_thr=θ_thr, topmixlayer=topmixlayer)
    end
    
    return decop_hgt
end
# ----/


# **************************************************************************
# Sub-routine to extract Wind Direction and Speed based on maximum water
# vapour transport.
function Collect_WindDir_using_WVT(INV::Dict, rs::Dict,
                              CBH::Vector, Decop_H::Vector,
                              H_wvt::Vector; ΔH = 0.15, R_lim=50f0)

    # ΔH is a +/- altitude from H_wvt in km
    
    # defining output variables:
    wind_dir = Vector{Vector{Float64}}(undef, 0)
    wind_range = Vector{Vector{Float64}}(undef, 0)
    wind_spd = Vector{Vector{Float64}}(undef, 0)
    
    for tidx = eachindex(H_wvt)

        # temporal variables for inversion:
	Inv_bot = INV[:HEIGHT][:, tidx]
	Inv_top = INV[:HEIGHT][:, tidx] .+ INV[:ΔHEIGHT][:, tidx]

        # ii_wdinr are the Prodile indexes where H_wvt fits criterium:
        # * decopling < Hwvt indexes [decopling : Hwvt or CBH] if Hwvt < CBH
        # * decopling ≥ Hwvt indexes [Inv_bot : Hwvt] or Hwvt ± ΔH if Hwvt ∉ INV
	ii_wdir = let HWVT = H_wvt[tidx] #δH = 0.1
	    i1, i2 = if !isnan(Decop_H[tidx]) && (Decop_H[tidx] < HWVT)
		(Decop_H[tidx], max(CBH[tidx], HWVT) )
	    else
		findfirst(Inv_bot .< HWVT .< Inv_top) |> x->!isnothing(x) ? (Inv_bot[x], HWVT) : HWVT .+ ΔH*[-1, 1]
	    end

	    # retrieving all Profile indexes within [i1, i2] above the first one
            # to avoid effect of low atmosphere e.g. strong surface inversion:
	    findall(i1 .≤ rs[:height] .≤ i2) |> x->filter!(>(1), x)
	end

	#feeding output variables:
        # * Wind direction, Wind speed, Wind range normalized to 50km
	push!(wind_dir, rs[:WDIR][ii_wdir, tidx])
	push!(wind_spd, rs[:WSPD][ii_wdir, tidx])
	dummy = let ws= rs[:WSPD][ii_wdir, tidx]
	    R_lim*ws./maximum(ws)
	end
	push!(wind_range, dummy)
    end

    return wind_dir, wind_spd, wind_range
end
# end of file
