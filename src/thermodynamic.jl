# This file is part of the ATMOStool.jl module and
# contains atmospheric thermodynamic functions
# 
# See LICENCE

### ********************************************************************
### Defining global Atmospheric constatnts:
###
### Boltzman gas constant R (dry air)
##const Rd = 287.058   # [J/kg/K]
##
### Rv gas constant water vapour:
##const Rv = 461.5  # [J/kg/K]
##
### Specific Heat capacity (isobaric)
##const c_p = 1e3  # [J/kg/K]
##
### Ratio R/c_p
##const κ = Rd/c_p  # ~ 0.286
##
### Reference meteorological standard Pressure P0 [hPa]
##const P₀ = 1013.25
##
### Meteorological standard Temperature reference level [K]
##const T₀ = 288.15
##
### Lapse-rate reference level [K/m]
##const Lp₀ = -0.0065
##
### Gravity acceleration g [m/s²]
##const g₀ = 9.81
##
##const ϵ = Rd/Rv   # ~ 0.622
##
### M = molar mass of Earth's air: 0.0289644 kg/mol
##const M_air =  0.0289644
##
### R* = universal gas constant: 8.3144598 J/(mol·K)
##const R = 8.3144598
    
# **********************************************************************
# Function potential temperature
"""
    θ = θ(T, P)

INPUT:
* T -> Temperature [K]
* P -> Pressure [hPa]

OUTPUT:
* θ -> Potential temperature [K]
"""
function θ(T, P)
    return T*(P₀/P)^κ
end
function θ(T::Array, P::Array)::Array
    return θ.(T, P)
end
#function θ(T::Vector, P::Vector)::Vector
#    return θ.(T, P)
#end

# ----/
# ********************************************************************
# Function Virtual Temperature
"""
 ____________________________________________________________________
 Subroutine to calculate the Virtual Temperature given T, MIXR

 --------------------------------------------------------------------
 Function convert Temperature to Virtual Temp.
 * -> T    : Temperature [K]
 * -> MIXR : Vapour Mixing ratio [kg/kg]
 * <- Tv   : Virtual Temperature [K]

 ---
 (c) 2020, Pablo Saavedra G.
 Geophysical Institute, University of Bergen
 See LICENSE
 ---
"""
function VirtualTemperature(T::Real, HU::Real; flag_mixr::Bool=true)

    # if qv give, then convert it to mixing ratio:
    mᵣ = flag_mixr ? HU : qx_to_mixr(HU)
    
    TV = T*(1.0 + mᵣ/ϵ )/(1.0 + mᵣ)

    return TV
end
# or
function VirtualTemperature(T::Array, HU::Array; flag_mixr::Bool=true)
    return VirtualTemperature.(T, HU, flag_mixr=flag_mixr)
end
# ----/

# *********************************************************************
# Function Virtual Potential Temperature
"""
Virtual Potential Temperature

θ_v = Theta_virtual(T, P, qv)

Input:
* T : Temperature [K]
* P : Pressure [hPa]
* qv: water vapour mixing ratio [g/g]
Output:
* θ_v: virtual potential temperature [K]
"""
function Theta_virtual(T, P, qv)
    # Function to compute the virtual Potential Temperature:
    
    T_virtual = VirtualTemperature(T, qv, flag_mixr=false)
    return θ(T_virtual, P)  
end
# or
function Theta_virtual(T::Array, P::Array, qv::Array)::Array
    return Theta_virtual.(T, P, qv)
end
# ----/

"""
! -------------------------------------------------------------------
! Function to convert Specific [kg/kg] to mixing ratio [kg/kg].
! Input, Q_x  : Specific content in [kg/kg]
! Output, mixr : Mixing ratio in [kg/kg]
! ---
"""
function qx_to_mixr(Q_x)
  
  return Q_x/(1.0 - Q_x)
end #function qx_to_mixr
function qx_to_mixr(Q_x::Array)::Array
    return qx_to_mixr.(Q_x)
end
# ----/

"""
! -------------------------------------------------------------------
! Function to convert mixing ratio [kg/kg] to Specific [kg/kg].
! Input,  mixr : Mixing ratio in [kg/kg]
! Output, Q_x  : Specific content in [kg/kg]
! ---
"""
function mixr_to_qx(mixr)
  
    return mixr/(1.0 + mixr)
end #function qx_to_mixr
function mixr_to_qx(mixr::Array)::Array
    return mixr_to_qx.(mixr)
end
# ----/


"""
-------------------------------------------------------------
 Function to calculate Partial pressure water vapour

```julia-repl
julia> eₛ = PWS(T)
```
WHERE:
```julia-repl
* T::AbstractFloat Temperature [K]
* eₛ::AbstractFloat Partial water vapour pressure [hPa]
```
---
"""
function PWS(T)

    COEFF = 2.16679 # [g K J^-1]
    Tc = 647.096  # critical temperature [K]
    Pc = 220640   # critical pressure [hPa]
    B  = 0.6219907 # constant for air [kg/kg]
    CC = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
    EE = (1.0, 1.5, 3.0, 3.5, 4.0, 7.5)

    η = 1.0 - T/Tc
    A = 0.0

    map(1:6) do i
        A = A + CC[i]*(η^EE[i])
    end
    
    eₛ = Pc*exp(A*Tc/T)
    return eₛ
end
# ----/

"""
 ELEMENTAL FUNCTION Convert mixing ratio to Relative Humidity.
 HUMIDITY CONVERSION FORMULAS
 Calculation formulas for humidity (B210973EN-F)
 By (c) VAISALA 2003
 https://www.hatchability.com/Vaisala.pdf

 ---
 (c) 2019 Pablo Saavedra G.
 Geophysical Institute, University of Bergen
 See LICENSE

 USAGE:
```julia-repl
julia> RH = mixt_to_rh(mᵣ, P, T)
```
WHERE:
* mᵣ::AbstractFloat Vapour mixing ratio [kg/kg]
* P::AbstractFloat Pressure [hPa]
* T::AbstractFloat Temperature [K]
* RH::AbstractFloat Relative Humidity [%]

The input variables mᵣ, P, and T can also be <:Matrix
"""
function mixr_to_rh(mᵣ, P, T)
    B  = 0.6219907 # constant for air [kg/kg]

    #RH = 100*PW/PWS
    RH = 100*mᵣ*P/(mᵣ + B)/PWS(T)
    return RH
end ##function mixr_to_rh
function mixr_to_rh(mᵣ::Matrix, P::Matrix, T::Matrix)::Matrix
    return mixr_to_rh.(mᵣ, P, T)
end
# ----/

"""
 ELEMENTAL FUNCTION Specific Humidity to Relative Humidity
USAGE:
```julia-repl
julia> RH = qv_to_rh(Qᵥ, P, T)
```
WHERE:
* Qᵥ : Specific Humidity [kg/kg]
* P  : Pressure [hPa]
* T  : Temperature [K]
* RH : Relative Humidity [%]
 ---
"""
function qv_to_rh(Qᵥ, P, T)
  
  # Converting Specific Humidity to Mixing Ratio:
  mᵣ = qx_to_mixr(Qᵥ)
  
  RH = mixr_to_rh(mᵣ, P, T)

  return RH
end #function qv_to_rh
# ----/

"""
 --------------------------------------------------------------------
 Function to convert kg/kg to kg/m^3

 USAGE:
```julia-repl
 julia> RHOx = MassRatio2MassVolume(Qx::Matrix, T::Matrix, P::Matrix)
 julia> RHOx = MassRatio2MassVolume(Qx::Number, T::Number, P::Number)
```
 where:
```julia-repl
 * Qx  : Specific quantity e.g. humidity [kg/kg]
 * T   : Temperature [K]
 * P   : Pressure [hPa]
 * RHOx: Output Specific quantity [kg/m^3]
```
"""
function MassRatio2MassVolume(Q_x, T, P)

  MIXR = qx_to_mixr(Q_x)
  Tv = VirtualTemperature(T, MIXR)
  RHO_air = (1f2P)/Tv/Rd  # [kg/m^3]
  RHO_x   = Q_x*RHO_air      # [kg/m^3]
  
  return RHO_x
end #function MassRatio2MassVolume
function MassRatio2MassVolume(Q_x::Array, T::Array, P::Array)
    
    return MassRatio2MassVolume.(Q_x, T, P)
end
# ----/

"""
! ----------------------------------------------------------------------
! Function to convert dew point to relative humidity
! -> T  : Temperature [K]
! -> Td : Dew point []
! <- RH : Relative Humidity [%]
! ---
"""
function dewpoint_to_rh(T, Td)

    E_0 = 6.11    # [hPa]
    T_0 = 273.15  # [K]
    L_Rv = 5423.0 # [K]

    # according to an approximation of the Clausius-Clapeyron Equation (CCE):
    CCE(Tx) = E_0*exp(L_Rv*(1/T_0 - 1/Tx))
    
    RH = 100.0*CCE(Td)/CCE(T)
    return RH
end
function dewpoint_to_rh(T::Array, Td::Array)::Array
    return dewpoint_to_rh.(T, Td)
end
# ----/

"""
! -----------------------------------------------------------------------
! Function to convert from relative humidity to mixing ratio
! -> RH  : relative humidity [%]
! -> T   : temperature [K]
! -> P   : pressure [hPa]
! <- mixr: mixing ratio [g/g]
"""
function rh_to_mixr(RH, T, P)

    B  = 0.6219907 # constant for air [kg/kg]

    Es = PWS(T)
    mixr = RH*B*Es/(100.0*P - RH*Es)
    return mixr
end
function rh_to_mixr(RH::Array, T::Array, P::Array)::Array
    return rh_to_mixr.(RH, T, P)
end
# ----/


# ***************************************************************
"""
Barometric Formula for Pressure
See [Barometric Formula](https://en.wikipedia.org/wiki/Barometric_formula)

Pa = barometric_formula(H)

where:
H -> altitude in [m]
Pa <- Pressure level in [hPa]

"""
function barometric_formula(H::Number)

    faktor = -g₀*M_air/R/Lp₀
    
    P_h = (T₀ + Lp₀*(H))/T₀
    P_h ^= faktor
    P_h *= P₀

    return P_h
end


# ********************************************************************
"""
--------------------------------------------------------------------
Function to estimate the Pressure level from given Altitude.

USAGE:
> altitude2pressure(H_in::Real)
> altitude2pressure(H_in::Real; H_ref=[1,2,3,4], Pa_ref=[1000,900,850,820])

INPUT:
* H_in : Altitude for which pressure is calculated [m]
* H_ref: (Optional) Vector with reference altitudes [m]
* Pa_ref: (Optional) Vector same size as H_ref with reference Pressure

OUTPUT:
* Pa: Pressure corresponding to altitude H_in, same units as Pa_ref but only
if reference Altitude (H_ref) and Pressure (Pa_ref) levels are given.
Otherwise the Barometric Formula is used for the input altitude H_in [in meters]
 and return pressure in units of hPa. 
See [Barometric Formula](https://en.wikipedia.org/wiki/Barometric_formula)


"""
function altitude2pressure(H_in::Real; H_ref=nothing, Pa_ref=nothing)
    
    if !(typeof(H_ref) <: AbstractVector) || !(typeof(Pa_ref) <: AbstractVector)
	@warn "Using barometric formula as default"

    elseif (extrema(H_ref) |>  x-> (x[1] ≤ H_in ≤ x[2]))
        
	i₁ = findlast(H_ref .≤ H_in)
        
	i₂ = findfirst(H_ref[i₁:end] .≥ H_in)
        
	# interpolating Pa according the given referencs H_ref and Pa_ref:
        return i₁ == i₂ ? Pa_ref[i₁] : (Pa_ref[i₂] - Pa_ref[i₁])/(H_ref[i₂] - H_ref[i₁])*(H_in - H_ref[1])+Pa_ref[1]
    else
        @error "$(H_in) lays outside the bounds of H_ref!"
    end
	
    # otherwise returning Pa from barometric formula 
    return ATMOStools.barometric_formula(H_in) 
end
# ----/

# ***************************************************************************
# end of file.
