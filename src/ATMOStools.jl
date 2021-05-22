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
    N2 = (Δθv./ΔZ)./θv[1:end-1,:]
    N2 *= gravity_0
    # wind shear gradient:
    Ri = N2./(ΔU./ΔZ).^2
    return N2, Ri
end
function Ri_N(rs::Dict)
    N2, Ri = Ri_N(1f3*rs[:height], rs[:WSPD], rs[:qv], rs[:θ])
    return N2, Ri
end


end #module
