# ATMOStools.jl
## Atmospheric Physics and Meteorological Tools

This is a Julia package containing a set to tools as functions and constants related to atmospheric physics and meteorology.

## Installation
```
julia> ]
(v1.0) pkg> add git@github.com:pablosaa/ATMOStools.jl.git
```

## Usage
```
julia> using ATMOStools, const ATM=ATMOStools
julia> ?
help?> ATM.calculate_∇WVT
```
it will show the help for that function as:
```
Function to return the derivative of IVT
given that
IVT = -1/g₀∫Qv*Ws*dP

> δIVT = calculate_∇WVT(P, WS, Qv, H)
returns:
  δIVT = -∇f/g₀

where:
* ∇f = ∂Φᵥ/∂z
* ∂Φᵥ = (Qv*Ws)dP

Note:
P pressure in [hPa], WS in [m s⁻¹], Qv [g/g] and H [km]
julia>
```

DISCLAMER: The suit of functions are being under constant development.

--<br>
(c) 2020, Pablo Saavedra Garfias<br>
[contact: pablo.saavedra [at] uni-leipzig.de](mailto:pablo.saavedra@uni-leipzig.de)<br>
University of Leipzig<br>
Faculty of Physics and Geosciences<br>

See LICENSE
