
#   This file is part of UnitSystems.jl. It is licensed under the MIT license
#   UnitSystems Copyright (C) 2020 Michael Reed

"""
$(convertext(:charge,"sqrt(planck(U,S)*current(U,S)/voltage(U,S))"))

Electric `charge` quantization (C, A⋅s), unit conversion factor.

```Julia
julia> charge(EMU,Metric) # C⋅abC⁻¹
$(charge(EMU,Metric))

julia> charge(EMU,ESU) # statC⋅abC⁻¹
$(charge(EMU,ESU))

julia> charge(ESU,Metric) # C⋅statC⁻¹
$(charge(ESU,Metric))

julia> charge(Metric,SI2019) # C⋅C⁻¹
$(charge(Metric,SI2019))

julia> charge(Hartree,SI2019) # C⋅𝘦⁻¹
$(charge(Hartree,SI2019))
```
"""
@pure charge(U::UnitSystem,S::UnitSystem) = unit(sqrt((planckreduced(S)*permeability(U)*lightspeed(U)*rationalization(U)*lorentz(U)^2)/(planckreduced(U)*permeability(S)*lightspeed(S)*rationalization(S)*lorentz(S)^2)))

"""
$(convertext(:current,"charge(U,S)/time(U,S)"))

Flow of electric `charge` per `time` or `current` (A, C⋅s⁻¹), unit conversion factor.

```Julia
julia> current(EMU,Metric) # A⋅Bi⁻¹
$(current(EMU,Metric))

julia> current(EMU,ESU) # statA⋅Bi⁻¹
$(current(EMU,ESU))

julia> current(ESU,Metric) # A⋅statA⁻¹
$(current(ESU,Metric))

julia> current(Metric,SI2019) # A⋅A⁻¹
$(current(Metric,SI2019))
```
"""
@pure current(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)/time(U,S))

"""
$(convertext(:voltage,"energy(U,S)/charge(U,S)"))

Electric potential difference or `energy` per `charge` (V, J⋅C⁻¹), unit conversion factor.

```Julia
julia> voltage(EMU,Metric) # V⋅abV⁻¹
$(voltage(EMU,Metric))

julia> voltage(EMU,ESU) # statV⋅abV⁻¹
$(voltage(EMU,ESU))

julia> voltage(ESU,Metric) # V⋅statV⁻¹
$(voltage(ESU,Metric))

julia> voltage(Metric,SI2019) # V⋅V⁻¹
$(voltage(Metric,SI2019))
```
"""
@pure voltage(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/charge(U,S))

"""
$(convertext(:capacitance,"charge(U,S)/voltage(U,S)"))

Electrical `capactiance` or `charge` per `voltage` (F, C⋅V⁻¹), unit conversion factor.

```Julia
julia> capacitance(EMU,Metric) # F⋅abF⁻¹
$(capacitance(EMU,Metric))

julia> capacitance(ESU,Metric) # F⋅statF⁻¹
$(capacitance(ESU,Metric))

julia> capactiance(Metric,SI2019) # F⋅F⁻¹
$(capacitance(Metric,SI2019))
```
"""
@pure capacitance(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)/voltage(U,S))

"""
$(convertext(:impedance,"voltage(U,S)/current(U,S)"))

Electrical `impedance` or `voltage` per `current` (Ω, S⁻¹, V⋅A⁻¹), unit conversion factor.

```Julia
julia> impedance(EMU,Metric) # Ω⋅abΩ⁻¹
$(impedance(EMU,Metric))

julia> impedance(ESU,Metric) # Ω⋅statΩ⁻¹
$(impedance(ESU,Metric))

julia> impedance(Metric,SI2019) # Ω⋅Ω⁻¹
$(impedance(Metric,SI2019))
```
"""
@pure impedance(U::UnitSystem,S::UnitSystem) = unit(voltage(U,S)/current(U,S))

"""
$(convertext(:conductance,"voltage(U,S)/current(U,S)"))

Electrical `conductance` or `current` per `voltage` (S, Ω⁻¹, A⋅V⁻¹), unit conversion factor.

```Julia
julia> conductance(EMU,Metric) # S⋅abS⁻¹
$(conductance(EMU,Metric))

julia> conductance(ESU,Metric) # S⋅statS⁻¹
$(conductance(ESU,Metric))

julia> conductance(Metric,SI2019) # S⋅S⁻¹
$(conductance(Metric,SI2019))
```
"""
@pure conductance(U::UnitSystem,S::UnitSystem) = unit(current(U,S)/voltage(U,S))

"""
$(convertext(:magneticflux,"energy(U,S)/lorentz(U,S)/current(U,S)"))

Surface `magneticflux` or `energy` per `current` (Wb, J⋅A⁻¹, V⋅s), unit conversion factor.

```Julia
julia> magneticflux(EMU,Metric) # Wb⋅Mx⁻¹
$(magneticflux(EMU,Metric))

julia> magneticflux(ESU,Metric) # Wb⋅statWb⁻¹
$(magneticflux(ESU,Metric))

julia> magneticflux(Metric,SI2019) # Wb⋅Wb⁻¹
$(magneticflux(Metric,SI2019))
```
"""
@pure magneticflux(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/lorentz(U,S)/current(U,S))

"""
$(convertext(:magneticinduction,"mass(U,S)/lorentz(U,S)/current(U,S)/time(U,S)^2"))

Magnetic `magneticinduction` or `magneticflux` per `area` (T, Wb⋅m⁻²), unit conversion factor.

```Julia
julia> magneticinduction(EMU,Metric) # T⋅G⁻¹
$(magneticinduction(EMU,Metric))

julia> magneticinduction(EMU,ESU) # statT⋅G⁻¹
$(magneticinduction(EMU,ESU))

julia> magneticinduction(Metric,SI2019) # T⋅T⁻¹
$(magneticinduction(Metric,SI2019))
```
"""
@pure magneticinduction(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/lorentz(U,S)/current(U,S)/time(U,S)^2)

"""
$(convertext(:inductance,"mass(U,S)*area(U,S)/charge(U,S)"))

Electro-`magneticflux` per `current` or `inductance` (H, Ω⋅s, Wb⋅A⁻¹), unit conversion factor.

```Julia
julia> inductance(EMU,Metric) # H⋅abH⁻¹
$(inductance(EMU,Metric))

julia> inductance(ESU,Metric) # H⋅statH⁻¹
$(inductance(ESU,Metric))

julia> inductance(Metric,SI2019) # H⋅H⁻¹
$(inductance(Metric,SI2019))
```
"""
@pure inductance(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)*area(U,S)/charge(U,S)^2)

# electromagnetics

"""
$(convertext(:electricinduction,"charge(U,S)/area(U,S)"))

Electric field displacement or `electricinduction` (C⋅m⁻²), unit conversion factor.

```Julia
julia> electricinduction(EMU,Metric) # C⋅cm²⋅abC⁻¹⋅m⁻²
$(electricinduction(EMU,Metric))

julia> electricinduction(ESU,Metric) # C⋅cm²⋅statC⁻¹⋅m⁼²
$(electricinduction(ESU,Metric))

julia> electricinduction(Metric,SI2019) # C⋅C⁻¹
$(electricinduction(Metric,SI2019))
```
"""
@pure electricinduction(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)/area(U,S))

"""
$(convertext(:chargedensity,"charge(U,S)/volume(U,S)"))

Volume `chargedensity` or `charge` per `volume` (C⋅m⁻³), unit conversion factor.

```Julia
julia> chargedensity(EMU,Metric) # C⋅cm³⋅abC⁻¹⋅m⁻³
$(chargedensity(EMU,Metric))

julia> chargedensity(ESU,Metric) # C⋅cm³⋅statC⁻¹⋅m⁼³
$(chargedensity(ESU,Metric))

julia> chargedensity(Metric,SI2019) # C⋅C⁻¹
$(chargedensity(Metric,SI2019))
```
"""
@pure chargedensity(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)/volume(U,S))

"""
$(convertext(:currentdensity,"current(U,S)/area(U,S)"))

Cross-section `currentdensity` or `current` per `area` (A⋅m⁻²), unit conversion factor.

```Julia
julia> currentdensity(EMU,Metric) # A⋅cm²⋅Bi⁻¹⋅m⁻²
$(currentdensity(EMU,Metric))

julia> currentdensity(ESU,Metric) # A⋅cm²⋅statA⁻¹⋅m⁼²
$(currentdensity(ESU,Metric))

julia> currentdensity(Metric,SI2019) # A⋅A⁻¹
$(currentdensity(Metric,SI2019))
```
"""
@pure currentdensity(U::UnitSystem,S::UnitSystem) = unit(current(U,S)/area(U,S))

"""
$(convertext(:conductivity,"conductance(U,S)/length(U,S)"))

Reciprocal `resistivity` or electrical `conductivity` (S⋅m⁻¹), unit conversion factor.

```Julia
julia> conductivity(EMU,Metric) # S⋅cm⋅abS⁻¹⋅m⁻¹
$(conductivity(EMU,Metric))

julia> conductivity(ESU,Metric) # S⋅cm⋅statS⁻¹⋅m⁼¹
$(conductivity(ESU,Metric))

julia> conductivity(Metric,SI2019) # S⋅S⁻¹
$(conductivity(Metric,SI2019))
```
"""
@pure conductivity(U::UnitSystem,S::UnitSystem) = unit(conductance(U,S)/length(U,S))

"""
$(convertext(:permittivity,"capacitance(U,S)*rationalization(U,S)/length(U,S)"))

Absolute `permittivity` or `capacitance` per `length` (F⋅m⁻¹), unit conversion factor.

```Julia
julia> permittivity(EMU,Metric) # F⋅cm⋅abF⁻¹⋅m⁻¹
$(permittivity(EMU,Metric))

julia> permittivity(ESU,Metric) # F⋅cm⋅statF⁻¹⋅m⁼¹
$(permittivity(ESU,Metric))

julia> permittivity(Metric,SI2019) # F⋅F⁻¹
$(permittivity(Metric,SI2019))
```
"""
@pure permittivity(U::UnitSystem,S::UnitSystem) = unit(capacitance(U,S)*rationalization(U,S)/length(U,S))

@doc """
$(convertext(:permeability,"permeability(S)/permeability(U)"))

Magnetic `permeability` or `inductance` per `length` (H⋅m⁻¹), unit conversion factor.

```Julia
julia> permeability(EMU,Metric) # H⋅cm⋅abH⁻¹⋅m⁻¹
$(permeability(EMU,Metric))

julia> permeability(ESU,Metric) # H⋅cm⋅statH⁻¹⋅m⁼¹
$(permeability(ESU,Metric))

julia> permeability(Metric,SI2019) # H⋅H⁻¹
$(permeability(Metric,SI2019))
```
""" permeability(::UnitSystem,::UnitSystem)

"""
    electricfield # V⋅m⁻¹
"""
@pure electricfield(U::UnitSystem,S::UnitSystem) = unit(voltage(U,S)/length(U,S))

"""
$(convertext(:magneticfield,"current(U,S)*rationalization(U,S)*lorentz(U,S)/length(U,S)"))

Magnetization or `magneticfield` or `current` per `length` (A⋅m⁻¹), unit conversion factor.

```Julia
julia> magneticfield(EMU,Metric) # A⋅m⁻¹⋅Oe⁻¹
$(magneticfield(EMU,Metric))

julia> magneticfield(ESU,Metric) # A⋅cm⋅m⁻¹⋅statA⁻¹
$(magneticfield(ESU,Metric))

julia> magneticfield(Metric,SI2019) # A⋅A⁻¹
$(magneticfield(Metric,SI2019))
```
"""
@pure magneticfield(U::UnitSystem,S::UnitSystem) = unit(current(U,S)*rationalization(U,S)*lorentz(U,S)/length(U,S))

"""
$(convertext(:exposure,"charge(U,S)/mass(U,S)"))

Ionizing radiation `exposure` or `charge` per `mass` (C⋅kg⁻¹), unit conversion factor.

```Julia
julia> exposure(EMU,Metric) # C⋅g⋅abC⁻¹⋅kg
$(exposure(EMU,Metric))

julia> exposure(EMU,ESU) # statC⋅abC⁻¹
$(exposure(EMU,ESU))

julia> expsure(ESU,Metric) # C⋅g⋅statC⁻¹⋅kg
$(exposure(ESU,Metric))

julia> exposure(Metric,SI2019) # C⋅C⁻¹
$(exposure(Metric,SI2019))
```
"""
@pure exposure(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)/mass(U,S))

"""
$(convertext(:resistivity,"impedance(U,S)*length(U,S)"))

Electrical `resistivity` or `impedance` by `length` (Ω⋅m), unit conversion factor.

```Julia
julia> impedance(EMU,Metric) # Ω⋅m⋅abΩ⁻¹⋅cm⁻¹
$(impedance(EMU,Metric))

julia> impedance(ESU,Metric) # Ω⋅m⋅statΩ⁻¹⋅cm⁻¹
$(impedance(ESU,Metric))

julia> impedance(Metric,SI2019) # Ω⋅Ω⁻¹
$(impedance(Metric,SI2019))
```
"""
@pure resistivity(U::UnitSystem,S::UnitSystem) = unit(impedance(U,S)*length(U,S))

"""
    linearchargedensity # C⋅m⁻¹
$(convertext(:linearchargedensity,"charge(U,S)/length(U,S)"))

Amount of `linearchargedensity` or `charge` per `length` (C⋅m⁻¹), unit conversion factor.

```Julia
julia> linearchargedensity(EMU,Metric) # C⋅cm⋅abC⁻¹⋅m⁻¹
$(linearchargedensity(EMU,Metric))

julia> linearchargedensity(ESU,Metric) # C⋅cm⋅statC⁻¹⋅m⁼¹
$(linearchargedensity(ESU,Metric))

julia> linearchargedensity(Metric,SI2019) # C⋅C⁻¹
$(linearchargedensity(Metric,SI2019))
```
"""
@pure linearchargedensity(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)/length(U,S))

"""
$(convertext(:magneticdipolemoment,"current(U,S)*lorentz(U,S)/area(U,S)"))

Magnetic dipole moment or `magneticdipolemoment` (J⋅T⁻¹, A⋅m²), unit conversion factor.

```Julia
julia> magneticdipolemoment(EMU,Metric) # J⋅G⋅T⁻¹⋅erg⁻¹
$(magneticdipolemoment(EMU,Metric))

julia> magneticdipolemoment(ESU,Metric) # J⋅statT⋅T⁻¹⋅erg⁼¹
$(magneticdipolemoment(ESU,Metric))

julia> magneticdipolemoment(Metric,SI2019) # A⋅A⁻¹⋅
$(magneticdipolemoment(Metric,SI2019))
```
"""
@pure magneticdipolemoment(U::UnitSystem,S::UnitSystem) = unit(current(U,S)*lorentz(U,S)*area(U,S))

"""
$(convertext(:mobility,"charge(U,S)*time(U,S)/mass(U,S)"))

Electron `mobility` in solid state physics (m²⋅V⁻¹⋅s⁻¹, A⋅s²⋅kg⁻¹), unit conversion factor.

```Julia
julia> mobility(EMU,Metric) # C⋅g⋅abC⁻¹⋅kg
$(mobility(EMU,Metric))

julia> mobility(EMU,ESU) # statC⋅abC⁻¹
$(mobility(EMU,ESU))

julia> mobility(ESU,Metric) # C⋅g⋅statC⁻¹⋅kg
$(mobility(ESU,Metric))

julia> mobility(Metric,SI2019) # C⋅C⁻¹
$(mobility(Metric,SI2019))
```
"""
@pure mobility(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)*time(U,S)/mass(U,S))

"""
$(convertext(:reluctance,"1/inductance(U,S)"))

Magnetic `reluctance` or magnetic resistance (H⁻¹), unit conversion factor.

```Julia
julia> reluctance(EMU,Metric) # abH⋅H⁻¹
$(reluctance(EMU,Metric))

julia> reluctance(ESU,Metric) # statH⋅H⁻¹
$(reluctance(ESU,Metric))

julia> reluctance(Metric,SI2019) # H⋅H⁻¹
$(reluctance(Metric,SI2019))
```
"""
@pure reluctance(U::UnitSystem,S::UnitSystem) = inductance(S,U)

"""
$(convertext(:vectorpotential,"magneticflux(U,S)/length(U,S)"))

Magnetic `vectorpotential` (Wb⋅m⁻¹), unit conversion factor.

```julia
julia> vectorpotential(EMU,Metric) # Wb⋅cm⋅Mx⁻¹⋅m⁻¹
$(vectorpotential(EMU,Metric))

julia> vectorpotential(ESU,Metric) # Wb⋅cm⋅statWb⁻¹⋅m⁻¹
$(vectorpotential(ESU,Metric))

julia> vectorpotential(Metric,SI2019) # Wb⋅Wb⁻¹
$(vectorpotential(Metric,SI2019))
```
"""
@pure vectorpotential(U::UnitSystem,S::UnitSystem) = unit(magneticflux(U,S)/length(U,S))

"""
$(convertext(:magneticmoment,"magneticflux(U,S)*length(U,S)"))

Amount of `magneticmoment` or `magneticflux` by `length` (Wb⋅m), unit conversion factor.

```Julia
julia> magneticmoment(EMU,Metric) # Wb⋅m⋅Mx⁻¹⋅cm⁻¹
$(magneticmoment(EMU,Metric))

julia> magneticmoment(ESU,Metric) # Wb⋅m⋅statWb⁻¹⋅cm⁻¹
$(magneticmoment(ESU,Metric))

julia> magneticmoment(Metric,SI2019) # Wb⋅Wb⁻¹
$(magneticmoment(Metric,SI2019))
```
"""
@pure magneticmoment(U::UnitSystem,S::UnitSystem) = unit(magneticflux(U,S)*length(U,S))

"""
$(convertext(:rigidity,"magneticinduction(U,S)*length(U,S)"))

Electromagnetic `rigidity` (T⋅m), unit conversion factor.

```Julia
julia> rigidity(EMU,Metric) # T⋅m⋅G⁻¹⋅cm⁻¹
$(rigidity(EMU,Metric))

julia> rigidity(EMU,ESU) # statT⋅G⁻¹
$(rigidity(EMU,ESU))

julia> rigidity(Metric,SI2019) # T⋅T⁻¹
$(rigidity(Metric,SI2019))
```
"""
@pure rigidity(U::UnitSystem,S::UnitSystem) = unit(magneticinduction(U,S)*length(U,S))

"""
$(convertext(:susceptibility,"1/permeability(U,S)"))

Magnetic `susceptibility` or `length` per `inductance` (m⋅H⁻¹), unit conversion factor.

```Julia
julia> susceptibility(EMU,Metric) # H⋅cm⋅abH⁻¹⋅m⁻¹
$(susceptibility(EMU,Metric))

julia> susceptibility(ESU,Metric) # H⋅cm⋅statH⁻¹⋅m⁼¹
$(susceptibility(ESU,Metric))

julia> susceptibility(Metric,SI2019) # H⋅H⁻¹
$(susceptibility(Metric,SI2019))
```
"""
@pure susceptibility(U::UnitSystem,S::UnitSystem) = permeability(S,U)

# WARNING unchecked: rigidity, magneticmoment, vectorpotential, reluctance, mobility, electricflux, linearchargedensity, exposure, currentdensity, chargedensity, conductivity

# CGS extra

"""
$(convertext(:electricflux,"voltage(U,S)*length(U,S)"))

Amount of `electricflux` or `voltage` by `length` (V⋅m) # wikipedia CGS page has error ?

```Julia
julia> electricflux(EMU,Metric) # V⋅m⋅abV⁻¹⋅cm⁻¹
$(electricflux(EMU,Metric))

julia> electricflux(EMU,ESU) # statV⋅abV⁻¹
$(electricflux(EMU,ESU))

julia> electricflux(ESU,Metric) # V⋅m⋅statV⁻¹⋅cm⁻¹
$(electricflux(ESU,Metric))

julia> electricflux(Metric,SI2019) # V⋅V⁻¹
$(electricflux(Metric,SI2019))
```
"""
@pure electricflux(U::UnitSystem,S::UnitSystem) = unit(voltage(U,S)*length(U,S))

"""
$(convertext(:electricdipolemoment,"charge(U,S)*length(U,S)"))

Electric dipole moment or `electricdipolemoment` (C⋅m), unit conversion factor.

```Julia
julia> electricdipolemoment(EMU,Metric) # C⋅m⋅abC⁻¹⋅cm⁻¹
$(electricdipolemoment(EMU,Metric))

julia> electricdipolemoment(ESU,Metric) # C⋅m⋅statC⁻¹⋅cm⁼¹
$(electricdipolemoment(ESU,Metric))

julia> electricdipolemoment(Metric,Gauss)/1e-18 # D⋅C⁻¹⋅m⁻¹
$(electricdipolemoment(Metric,Gauss)/1e-18)

julia> electricdipolemoment(Metric,SI2019) # C⋅C⁻¹⋅
$(electricdipolemoment(Metric,SI2019))
```
"""
@pure electricdipolemoment(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)*length(U,S))
