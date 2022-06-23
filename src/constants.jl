# ********************************************************************
# Defining global Atmospheric constatnts:
#
# Boltzman gas constant R (dry air)
const Rd = 287.058   # [J/kg/K]

# Rv gas constant water vapour:
const Rv = 461.5  # [J/kg/K]

# Specific Heat capacity (isobaric)
const c_p = 1e3  # [J/kg/K]

# Ratio R/c_p
const κ = Rd/c_p  # ~ 0.286

# Reference meteorological standard Pressure P0 [hPa]
const P₀ = 1013.25

# Meteorological standard Temperature reference level [K]
const T₀ = 288.15

# Lapse-rate reference level [K/m]
const Lp₀ = -0.0065

# Gravity acceleration g [m/s²]
const g₀ = 9.81

const ϵ = Rd/Rv   # ~ 0.622

# M = molar mass of Earth's air: 0.0289644 kg/mol
const M_air =  0.0289644

# R* = universal gas constant: 8.3144598 J/(mol·K)
const R = 8.3144598

# Stephan-Boltzman constant: W m⁻² K⁻⁴
const σ = 5.670374419E−8
