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
* -> WSPD  : wind speed [m/s]
* -> QV    : specific humidity [kg/m³]
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
function Ri_N(height::Vector, WSPD::Matrix, QV::Matrix, θ::Matrix, T::Matrix;
              Tₛ::Vector=[], Hₛ = 0, WSPDₛ = 0)

    # calculating virtual potential temperatures:
    θᵥ = VirtualTemperature(θ, QV)
    VT_K = VirtualTemperature(T .+ 273.15, QV)

    if isempty(Tₛ)
        TVₛ_K = T[1,:] .+ 273.15       
        θₛ = θᵥ[1,:]
        
    else
        @assert length(Tₛ)==size(T, 2) error("Surface temperature Tₛ length does not match dim 2 of T")
        
        θₛ = VirtualTemperature(θ(Tₛ, P₀), QV[1,:])
        TVₛ_K = VirtualTemperature(Tₛ, QV[1,:])
    end
    
    Δθᵥ = θᵥ .- θₛ'   # θᵥ[1,:]'  # K
    ΔZ = height .- Hₛ  # height[1]  # m
    ΔU = WSPD .- WSPDₛ' # [1,:]'  # m/s
    Tᵥ = similar(T)
    Tᵥ[2:end, :] = 0.5(VT_K[1:end-1,:] .+ VT_K[2:end,:])
    Tᵥ[1, :] = 0.5(VT_K[1,:] .+ TVₛ_K)

    # calculating the Brunt-Väisälä frequency:
    N² = (Δθᵥ.*ΔZ)./Tᵥ
    N² *= g₀
    
    # calculating Richardson-number: N²/wind shear gradient:
    Ri = N²./(ΔU).^2
    
    return N², Ri
end

function Ri_N(rs::Dict)
    N², Ri = Ri_N(1f3*rs[:height], rs[:WSPD], rs[:qv], rs[:θ], rs[:T])

    return N², Ri
end
# ----/

function estimate_Ri_PBLH(Ri::AbstractMatrix, height::Vector; thr_ri=0.25)
    PBLH = [findfirst(>(0.25), x) for x in eachcol(Ri)] |> x->height[x] ;
    return PBLH
end
# or
function estimate_Ri_PBLH(rs::Dict; thr_ri=0.25)
    N², Ri = ATMOS.Ri_N(rs);
    return estimate_Ri_PBLH(Ri, rs[:height]; thr_ri = thr_ri)
end
# ----/


# **************************************************************************
# Sub-routine to extract Wind Direction and Speed based on maximum water
# vapour transport.
function Extimate_windDir_WVT(INV::Dict, rs::Dict,
                              CBH::Vector, Decop_H::Vector,
                              H_wvt::Vector; ΔH = 0.15)

    # ΔH is a +/- altitude from H_wvt in km
    
    # defining output variables:
    wind_dir = Vector{Vector{Float64}}(undef, 0)
    wind_range = Vector{Vector{Float64}}(undef, 0)
    wind_spd = Vector{Vector{Float64}}(undef, 0)
    
    for tidx = eachindex(H_wvt)

        # temperal variables for inversion:
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

    return wind_dir, wind_spd
end
# end of file
