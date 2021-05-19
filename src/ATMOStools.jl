module ATMOStools

using Printf

"""
ATMOStools a set of functions and constants useful for atmospheric physics and meteorology.

"""
ATMOStools

include("thermodynamic.jl")
# ******************************************************************
# Calculating Integrated Water Vapour Transport
#
"""
Function to compute Integrated Water Vapour Transport
IVT, IVT_vec = getIVT(rs::Dict)

INPUT:
* rs::Dict -> Dictionary with Radiosonde or Model data

OUTPUT:
* IVT -> total Integrated Water Vapour Transport [kg/m/s]
* IVT_vec -> IVT by wind component separated (meridional, zonal)
"""
function getIVT(rs)
    ΔP = rs[:Pa][2:end,:] - rs[:Pa][1:end-1,:]
    ΔP *= 1e3 # [Pa]
    tmp_u = rs[:QV] .* rs[:U]
    tmp_v = rs[:QV] .* rs[:V]
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
    IVT_vec /= 9.81
    # calculating total IVT:
    IVT = sqrt.(IVT_u.^2 + IVT_v.^2)
    IVT /= 9.81
    
    return IVT, IVT_vec
end
# ---

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
function getMaxIVT_θ(IVT::Array{Float64,2})
    
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
N², Ri = Ri_N(rs::Dict)

INPUT: 
* rs::Dict() -> Dictionary with Radiosonde or Model data
 OUTPUT:
* N2 -> Matrix with the Vaisaala frequency
* Ri -> Matrix with the bulk Richardson Number

"""
function Ri_N(rs)
    θv = rs[:θ].*(1.0 .+ 0.6*rs[:QV])  # K
    Δθv = θv[2:end,:] .- θv[1,:]'  # K
    ΔZ = 1e3*(rs[:height][2:end] .- rs[:height][1])  # m
    ΔU = rs[:U][2:end,:] .- rs[:U][1,:]'  # m/s
    N2 = (Δθv./ΔZ)./θv[1:end-1,:]
    N2 *= 9.81

    tmp = (ΔU./ΔZ).^2
    Ri = N2./tmp
    return N2, Ri
end


end #module
