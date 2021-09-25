module ATMOS

using Printf

"""
ATMOStools a set of functions and constants useful for atmospheric physics and meteorology.

"""
ATMOS

include("thermodynamic.jl")
# ******************************************************************
# Calculating Integrated Water Vapour Transport
#
"""
Function to compute Integrated Water Vapour Transport

IVT, IVT_vec = getIVT(Pa::Matrix, U::Matrix, V::Matrix, QV::Matrix)
IVT, IVT_vec = getIVT(rs::Dict)

INPUT:
* -> Pa : Pressure [hPa]
* -> U,V: Wind speed component u and v [m/s]
* -> QV : specific humidity [g/g]

or

* -> rs::Dict : Dictionary with Radiosonde or Model data

OUTPUT:
* <- IVT    : total Integrated Water Vapour Transport [kg/m/s]
* <- IVT_vec: IVT by wind component separated (meridional, zonal)
"""
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
    IVT_vec /= gravity_0
    # calculating total IVT:
    IVT = sqrt.(IVT_u.^2 + IVT_v.^2)
    IVT /= gravity_0
    
    return IVT, IVT_vec
end
function getIVT(rs::Dict)
    IVT, IVT_vec = getIVT(10f0*rs[:Pa], rs[:U], rs[:V], rs[:qv])
    return IVT, IVT_vec
end
# ----/

# ******************************************************************
# Get Maximum IVT index and values
#
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

# ********************************************************************
# Calculating Richardson Number
#
"""
Function to estimate the Richardson Number
N², Ri = Ri_N(height, WSPD, QV, θ)

INPUT: 
* -> height: profile altitudes [m]
* -> WSPD  : wind speed [m/s]
* -> QV    : specific humidity [kg/m³]
* -> θ     : potential temperature [K]
OR
* -> rs::Dict : Dictionary with Radiosonde or Model data
 OUTPUT:
* <- N2 : Vaisaala frequency [Hz]
* <- Ri : The bulk Richardson Number [-]

"""
function Ri_N(height::Vector, WSPD::Matrix, QV::Matrix, θ::Matrix)
    θv = VirtualTemperature(θ, QV)
    Δθv = θv[2:end,:] .- θv[1,:]'  # K
    ΔZ = height[2:end] .- height[1]  # m
    ΔU = WSPD[2:end,:] .- WSPD[1,:]'  # m/s
    N² = (Δθv./ΔZ)./θv[1:end-1,:]
    N² *= gravity_0
    # wind shear gradient:
    Ri = N²./(ΔU./ΔZ).^2
    return N², Ri
end
function Ri_N(rs::Dict)
    #N², Ri = Ri_N(1f3*rs[:height], rs[:WSPD], rs[:qv], rs[:θ])
    height = 1f3*rs[:height]
    WSPD = rs[:WSPD]
    QV = rs[:qv]
    θ = rs[:θ]
    θv = VirtualTemperature(θ, QV)
    Δθv = θv .- θv[1,:]'  # K
    ΔZ = height #[2:end] .- height[1:end-1]  # m
    ΔU = WSPD #[2:end, :] .- WSPD[1:end-1, :]  # m/s
    N² = (Δθv.*ΔZ)./θv[1,:]'
    N² *= gravity_0
    # wind shear gradient:
    Ri = N²./(ΔU).^2
    return N², Ri
end
# ----/

# **********************************************************************
# Calculating gradient of atmospheric profile
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


# *********************************************************************
# Capturing the inversion layers in a profile variable
function Indices_Inversion_Layers(T::AbstractMatrix, H::Vector; mxhg=7, ΔH=0.06, ΔT=0.5)

    m, nt = size(T)
    mxidx = findfirst(x-> x ≥ mxhg, H)

    Max_inv_layer = 4
    out_idx_bot = Matrix{Int32}(undef, Max_inv_layer, nt) = 0
    out_idx_top = Matrix{Int32}(undef, Max_inv_layer, nt) = 0
    
    # calculating the gradient respect to height h
    δTz = ∇ₕT(T, H)

    # defining dummy expression:
    ex =:(δTz < 0)
    for tline ∈ (1:nt)
        idx_bot = Vector{Int32}()
        idx_top = Vector{Int32}()

        for i ∈ (1:mxidx)
	    ex.args[2] = δTz[i, tline]
		
	    if !eval(ex)
            
	        if ex.args[1] == :<
		    push!(idx_bot, i) #i-1
		    ex.args[1] = :≥
	                else
		    push!(idx_top, i)
		    ex.args[1] = :<
                        end
	    end
        end
		
        # Merging layers closer than 60m form top to bottom
        # rs[:height][idx_bot[2:end]] .- rs[:height][idx_top[1:end-1]] 
        δH_tb = H[idx_bot[2:end]] .- H[idx_top[1:end-1]]
        idx_out = δH_tb |> x-> (x .≤ ΔH) |> findall
        deleteat!(idx_bot, idx_out.+1)
        deleteat!(idx_top, idx_out)
	
        # Dismissing layers with T inversion strength < 0.5 °C
        # rs[:T][idx_top, tline] .- rs[:T][idx_bot, tline]
        δTinv = T[idx_top, tline] .- T[idx_bot, tline]
        idx_inv = findall(x-> x ≤ ΔT, δTinv)
        deleteat!(idx_bot, idx_inv)
        deleteat!(idx_top, idx_inv)

        ninv = min(length(idx_bot), Max_inv_layer)
        out_idx_bot[1:ninv, tline] = idx_bot
        out_idx_top[1:ninv, tline] = idx_top
    end
    
    return out_idx_bot, out_idx_top
end
# ----/

# --
end #module
