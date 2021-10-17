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
        idx_bot = Vector{Int32}()
        idx_top = Vector{Int32}()

        # defining dummy expression:
        ex =:(δTz < 0)

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

# **********************************************************************
# extracting the inversion variables from the specified profiles
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

# --
end #module


##v, nt = size(idx_bot)
##
##  # Variable type
##  VarType = Matrix{Float32}(undef, nv, nt)
##  
##  # creating a set of variables to map the inversion variables:
##  # e.g. :Pa will be mapped to :PA for bottom and top inversion and
##  # to :ΔPA for strength
##  varset = Dict( x =>
##                 let tmp = uppercase(String(x))
##                 (Symbol(tmp), Symbol(:Δ, tmp))
##                 end
##                 for x in vars)
##  
##  # converting output keys to uppercase symbols
##  out = let tmp0 = @. uppercase(String(vars)) |> Symbol
##      # for inversion base variables:
##      tmp1 = Dict(x => VarType for x ∈ tmp0)
##      
##      # for inversin top - bottom differences variables:
##      tmp2 = Dict(Symbol(:Δ, x) => VarType for x in tmp0)
##
##      # convining both set of variables:
##      merge!(tmp1, tmp2)
##      
##      # Initializing all variables with NaN:
##      foreach(x-> tmp1[x] .= NaN, keys(tmp1))
##
##      #returning variable:
##      tmp1
##  end
##  
##  # filling output variable with inversion estimations:
##  foreach(vars) do X
##      for i ∈ (1:nt)
##          for j ∈ (1:nv)
##              i_top = idx_top[j,i]
##              i_bot = idx_bot[j,i]
##
##              i_bot<1 && continue
##
##              # for difference variables:
##              out[varset[X][2]][j, i] = ndims(rs[X])==1 ? rs[X][i_top] .- rs[X][i_bot] : rs[X][i_top, i] .- rs[X][i_bot, i]
##
##              # for bottom variables:
##              out[varset[X][1]][j, i] = ndims(rs[X])==1 ? rs[X][i_bot] : rs[X][i_bot, i]
##
##          end
##      end
##          
##  end
##
##  return out
