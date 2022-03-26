
#   This file is part of UnitSystems.jl
#   It is licensed under the MIT license
#   UnitSystems Copyright (C) 2020 Michael Reed
#       _           _                         _
#      | |         | |                       | |
#   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _
#  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` |
# | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| |
#  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_|
#
#   https://github.com/chakravala
#   https://crucialflow.com

@doc """
$(convertext(:charge,"sqrt(action(U,S)*current(U,S)/electricpotential(U,S))"))

Electric `charge` quantization (C, A⋅s), unit conversion factor.

```Julia
julia> charge(EMU,Metric) # C⋅abC⁻¹
$(charge(EMU,Metric))

julia> charge(EMU,ESU) # stC⋅abC⁻¹
$(charge(EMU,ESU))

julia> charge(ESU,Metric) # C⋅stC⁻¹
$(charge(ESU,Metric))

julia> charge(Metric,SI2019) # C⋅C⁻¹
$(charge(Metric,SI2019))

julia> charge(Hartree,SI2019) # C⋅𝘦⁻¹
$(charge(Hartree,SI2019))
```
""" charge(U::UnitSystem,S::UnitSystem), Q

@doc """
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
""" current

@doc """
$(convertext(:electricpotential,"energy(U,S)/charge(U,S)"))

Voltage or `electricpotential` or `energy` per `charge` (V, J⋅C⁻¹), unit conversion factor.

```Julia
julia> electricpotential(EMU,Metric) # V⋅abV⁻¹
$(electricpotential(EMU,Metric))

julia> electricpotential(EMU,ESU) # statV⋅abV⁻¹
$(electricpotential(EMU,ESU))

julia> electricpotential(ESU,Metric) # V⋅statV⁻¹
$(electricpotential(ESU,Metric))

julia> electricpotential(Metric,SI2019) # V⋅V⁻¹
$(electricpotential(Metric,SI2019))
```
""" electricpotential

@doc """
$(convertext(:capacitance,"charge(U,S)/electricpotential(U,S)"))

Electrical `capactiance` or `charge` per `electricpotential` (F, C⋅V⁻¹), unit conversion factor.

```Julia
julia> capacitance(EMU,Metric) # F⋅abF⁻¹
$(capacitance(EMU,Metric))

julia> capacitance(ESU,Metric) # F⋅cm⁻¹
$(capacitance(ESU,Metric))

julia> capactiance(Metric,SI2019) # F⋅F⁻¹
$(capacitance(Metric,SI2019))
```
""" capacitance

@doc """
$(convertext(:resistance,"electricpotential(U,S)/current(U,S)"))

Electrical `resistance` or `electricpotential` per `current` (Ω, S⁻¹, V⋅A⁻¹), unit conversion factor.

```Julia
julia> resistance(EMU,Metric) # Ω⋅abΩ⁻¹
$(resistance(EMU,Metric))

julia> resistance(ESU,Metric) # Ω⋅statΩ⁻¹
$(resistance(ESU,Metric))

julia> resistance(Metric,SI2019) # Ω⋅Ω⁻¹
$(resistance(Metric,SI2019))
```
""" resistance

@doc """
$(convertext(:conductance,"current(U,S)/electricpotential(U,S)"))

Electrical `conductance` or `current` per `electricpotential` (S, Ω⁻¹, A⋅V⁻¹), unit conversion factor.

```Julia
julia> conductance(EMU,Metric) # S⋅abS⁻¹
$(conductance(EMU,Metric))

julia> conductance(ESU,Metric) # S⋅statS⁻¹
$(conductance(ESU,Metric))

julia> conductance(Metric,SI2019) # S⋅S⁻¹
$(conductance(Metric,SI2019))
```
""" conductance(U::UnitSystem,S::UnitSystem)

@doc """
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
""" magneticflux(U::UnitSystem,S::UnitSystem)

@doc """
$(convertext(:magneticfluxdensity,"magneticflux(U,S)/area(U,S)"))

Magnetic induction or surface `magneticfluxdensity` (T, Wb⋅m⁻²), unit conversion factor.

```Julia
julia> magneticfluxdensity(EMU,Metric) # T⋅G⁻¹
$(magneticfluxdensity(EMU,Metric))

julia> magneticfluxdensity(EMU,ESU) # statT⋅G⁻¹
$(magneticfluxdensity(EMU,ESU))

julia> magneticfluxdensity(Metric,SI2019) # T⋅T⁻¹
$(magneticfluxdensity(Metric,SI2019))
```
""" magneticfluxdensity

@doc """
$(convertext(:inductance,"magneticflux(U,S)/current(U,S)"))

Electro-`magneticflux` per `current` or `inductance` (H, Ω⋅s, Wb⋅A⁻¹), unit conversion factor.

```Julia
julia> inductance(EMU,Metric) # H⋅abH⁻¹
$(inductance(EMU,Metric))

julia> inductance(ESU,Metric) # H⋅statH⁻¹
$(inductance(ESU,Metric))

julia> inductance(Metric,SI2019) # H⋅H⁻¹
$(inductance(Metric,SI2019))
```
""" inductance

# electromagnetics

@doc """
$(convertext(:electricfluxdensity,"charge(U,S)*rationalization(U,S)/area(U,S)"))

Electric field displacement or surface `electricfluxdensity` (C⋅m⁻²), unit conversion factor.

```Julia
julia> electricfluxdensity(EMU,Metric) # C⋅cm²⋅abC⁻¹⋅m⁻²
$(electricfluxdensity(EMU,Metric))

julia> electricfluxdensity(ESU,Metric) # C⋅cm²⋅statC⁻¹⋅m⁼²
$(electricfluxdensity(ESU,Metric))

julia> electricfluxdensity(Metric,SI2019) # C⋅C⁻¹
$(electricfluxdensity(Metric,SI2019))
```
""" electricfluxdensity

@doc """
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
""" chargedensity

@doc """
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
""" currentdensity

@doc """
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
""" conductivity

@doc """
$(convertext(:permittivity,"capacitance(U,S)*rationalization(U,S)/length(U,S)"))

Absolute `permittivity` or `capacitance` per `length` (F⋅m⁻¹), unit conversion factor.

```Julia
julia> permittivity(EMU,Metric) # F⋅cm⋅abF⁻¹⋅m⁻¹
$(permittivity(EMU,Metric))

julia> permittivity(ESU,Metric) # F⋅m⁼¹
$(permittivity(ESU,Metric))

julia> permittivity(Metric,SI2019) # F⋅F⁻¹
$(permittivity(Metric,SI2019))
```
""" permittivity(U::UnitSystem,S::UnitSystem)

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

@doc """
$(convertext(:electricfield,"electricpotential(U,S)/length(U,S)"))

The `electricpotential` per `length` or `electricfield` (V⋅m⁻¹), unit conversion factor.

```Julia
julia> electricfield(EMU,Metric) # V⋅cm⋅abV⁻¹⋅m⁻¹
$(electricfield(EMU,Metric))

julia> electricfield(EMU,ESU) # statV⋅abV⁻¹
$(electricfield(EMU,ESU))

julia> electricfield(ESU,Metric) # V⋅cm⋅statV⁻¹⋅m⁻¹
$(electricfield(ESU,Metric))

julia> electricfield(Metric,SI2019) # V⋅V⁻¹
$(electricfield(Metric,SI2019))
```
""" electricfield

@doc """
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
""" magneticfield

@doc """
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
""" exposure

@doc """
$(convertext(:resistivity,"resistance(U,S)*length(U,S)"))

Electrical `resistivity` or `resistance` by `length` (Ω⋅m), unit conversion factor.

```Julia
julia> resistance(EMU,Metric) # Ω⋅m⋅abΩ⁻¹⋅cm⁻¹
$(resistance(EMU,Metric))

julia> resistance(ESU,Metric) # Ω⋅m⋅statΩ⁻¹⋅cm⁻¹
$(resistance(ESU,Metric))

julia> resistance(Metric,SI2019) # Ω⋅Ω⁻¹
$(resistance(Metric,SI2019))
```
""" resistivity

@doc """
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
""" linearchargedensity

@doc """
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
""" magneticdipolemoment

@doc """
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
""" mobility

@doc """
$(convertext(:reluctance,"rationalization(U,S)*lorentz(U,S)^2/inductance(U,S)"))

Magnetic `reluctance` or magnetic resistance (H⁻¹, Gb⋅Mx⁻¹), unit conversion factor.

```Julia
julia> reluctance(EMU,Metric) # abH⋅H⁻¹
$(reluctance(EMU,Metric))

julia> reluctance(ESU,Metric) # statH⋅H⁻¹
$(reluctance(ESU,Metric))

julia> reluctance(Metric,SI2019) # H⋅H⁻¹
$(reluctance(Metric,SI2019))
```
""" reluctance # reciprocal: permeance -- different concept from inductace but same units

@doc """
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
""" vectorpotential

@doc """
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
""" magneticmoment

@doc """
$(convertext(:rigidity,"magneticfluxdensity(U,S)*length(U,S)"))

Electromagnetic `rigidity` (T⋅m), unit conversion factor.

```Julia
julia> rigidity(EMU,Metric) # T⋅m⋅G⁻¹⋅cm⁻¹
$(rigidity(EMU,Metric))

julia> rigidity(EMU,ESU) # statT⋅G⁻¹
$(rigidity(EMU,ESU))

julia> rigidity(Metric,SI2019) # T⋅T⁻¹
$(rigidity(Metric,SI2019))
```
""" rigidity

@doc """
$(convertext(:susceptibility,"1/rationalization(U,S)"))

Magnetic/electric volume `susceptibility` (dimensionless), unit conversion factor.

```Julia
julia> susceptibility(EMU,Metric)
$(susceptibility(EMU,Metric))

julia> susceptibility(ESU,Metric)
$(susceptibility(ESU,Metric))

julia> susceptibility(Metric,SI2019)
$(susceptibility(Metric,SI2019))
```
""" susceptibility # magneticdipolemoment(U,S)/magneticfield(U,S)/volume(U,S)

# WARNING unchecked: rigidity, magneticmoment, vectorpotential, mobility, linearchargedensity, exposure

# CGS extra: polarizability, permeance, magnetic-current? + density, magneticresistance

@doc """
$(convertext(:electricflux,"electricpotential(U,S)*length(U,S)"))

Amount of `electricflux` or `electricpotential` by `length` (V⋅m), unit conversion factor.

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
""" electricflux

@doc """
$(convertext(:electricdipolemoment,"charge(U,S)*length(U,S)"))

Electric dipole moment or `electricdipolemoment` (C⋅m), unit conversion factor.

```Julia
julia> electricdipolemoment(EMU,Metric) # C⋅m⋅abC⁻¹⋅cm⁻¹
$(electricdipolemoment(EMU,Metric))

julia> electricdipolemoment(ESU,Metric) # C⋅m⋅statC⁻¹⋅cm⁼¹
$(electricdipolemoment(ESU,Metric))

julia> electricdipolemoment(Metric,SI2019) # C⋅C⁻¹
$(electricdipolemoment(Metric,SI2019))
```
""" electricdipolemoment
#julia> electricdipolemoment(Metric,Gauss)/1e-18 # D⋅C⁻¹⋅m⁻¹
#$(electricdipolemoment(Metric,Gauss)/1e-18)

@doc """
$(convertext(:magneticpotential,"magneticflux(U,S)*reluctance(U,S)"))

Magnetomotive force or `magneticpotential` (A, Gb), unit conversion factor.

```Julia
julia> magneticpotential(EMU,Metric) # A⋅Gb⁻¹
$(magneticpotential(EMU,Metric))

julia> magneticpotential(Metric,SI2019) # A⋅A⁻¹
$(magneticpotential(Metric,SI2019))
```
""" magneticpotential

@doc """
$(convertext(:polestrength,"magneticdipolemoment(U,S)/length(U,S)"))

Magnetic `polestrength` is analogous to `charge` (A⋅m⁻¹), unit conversion factor.

```Julia
julia> polestrength(EMU,Metric) # A⋅m⁻¹⋅pole⁻¹
$(polestrength(EMU,Metric))

julia> polestrength(Metric,SI2019) # A⋅A⁻¹⋅
$(polestrength(Metric,SI2019))
```
""" polestrength

@doc """
$(convertext(:permeance,"1/reluctance(U,S)"))

Magnetic `permeance` or magnetic conductance (H, Mx⋅Gb⁻¹), unit conversion factor.

```Julia
julia> permeance(EMU,Metric) # abH⋅H⁻¹
$(permeance(EMU,Metric))

julia> permeance(ESU,Metric) # statH⋅H⁻¹
$(permeance(ESU,Metric))

julia> permeance(Metric,SI2019) # H⋅H⁻¹
$(permeance(Metric,SI2019))
```
""" permeance

@doc """
$(convertext(:specificsusceptibility,"susceptibility(U,S)/density(U,S)"))

Magnetic/electric mass specific `susceptibility` (m³⋅kg⁻¹), unit conversion factor.

```Julia
julia> specificsusceptibility(EMU,Metric) # m³⋅g⋅kg⁻¹⋅cm⁻³
$(specificsusceptibility(EMU,Metric))

julia> specificsusceptibility(ESU,Metric) # m³⋅g⋅kg⁻¹⋅cm⁻³
$(specificsusceptibility(ESU,Metric))

julia> specificsusceptibility(Metric,SI2019) # m³⋅kg⋅kg⁻¹⋅m⁻³
$(specificsusceptibility(Metric,SI2019))
```
""" specificsusceptibility

@doc """
$(convertext(:magnetizability,"magneticmoment(U,S)/magneticfluxdensity(U,S)"))

Quantity of `magneticmoment` per `magneticfluxdensity` (m⁻¹), unit conversion factor.

```Julia
julia> magnetizability(EMU,Metric) # cm⋅m⁻¹
$(magnetizability(EMU,Metric))

julia> magnetizability(ESU,Metric) # cm⋅m⁻¹
$(magnetizability(ESU,Metric))

julia> magnetizability(Metric,SI2019) # m⋅m⁻¹
$(magnetizability(Metric,SI2019))
```
""" magnetizability

@doc """
$(convertext(:electricpolarizability,"electricdipolemoment(U,S)/electricfield(U,S)"))

Polarizability or `electricdipolemoment` per `electricfield` (C⋅m²⋅V⁻¹), unit conversion factor.

```Julia
julia> electricpolarizability(EMU,Metric) # C⋅m²⋅abV⋅abC⁻¹⋅cm⁻²⋅V⁻¹
$(electricpolarizability(EMU,Metric))

julia> electricpolarizability(ESU,Metric) # C⋅m²⋅statV⋅statC⁻¹⋅cm⁼²⋅V⁻¹
$(electricpolarizability(ESU,Metric))

julia> electricpolarizability(Metric,Gauss) # D⋅cm²⋅V⁻¹⋅C⁻¹⋅m⁻²⋅abV⁻¹
$(electricpolarizability(Metric,Gauss))

julia> electricpolarizability(Metric,SI2019) # C⋅V⋅C⁻¹⋅V⁻¹
$(electricpolarizability(Metric,SI2019))
```
""" electricpolarizability

@doc """
$(convertext(:magneticpolarizability,"magneticdipolemoment(U,S)/magneticfield(U,S)"))

Polarizability or `magneticdipolemoment` per `magneticfield` (m³), unit conversion factor.

```Julia
julia> electricpolarizability(EMU,Metric) # m³⋅cm⁻³
$(electricpolarizability(EMU,Metric))

julia> electricpolarizability(ESU,Metric) # m³⋅cm⁼³
$(electricpolarizability(ESU,Metric))

julia> electricpolarizability(Metric,Gauss) # cm³⋅m⁻³
$(electricpolarizability(Metric,Gauss))

julia> electricpolarizability(Metric,SI2019)
$(electricpolarizability(Metric,SI2019))
```
""" magneticpolarizability

@doc """
$(convertext(:magnetization,"magneticmoment(U,S)/volume(U,S)"))

Amount of `magneticmoment` per `volume` (Wb⋅m⁻²), unit conversion factor.

```Julia
julia> magnetization(EMU,Metric) # Wb⋅cm²⋅Mx⁻¹⋅m⁻²
$(magnetization(EMU,Metric))

julia> magnetization(ESU,Metric) # Wb⋅cm²⋅statWb⁻¹⋅m⁻²
$(magnetization(ESU,Metric))

julia> magnetization(Metric,SI2019) # Wb⋅Wb⁻¹
$(magnetization(Metric,SI2019))
```
""" magnetization

# specificmagnetization, mass magnetization = 1?
# magneticfluxdensity(Metric,EMU)/density(Metric,EMU)
@doc """
$(convertext(:specificmagnetization,"magneticmoment(U,S)/mass(U,S)"))

Amount of `magneticmoment` per `mass` (Wb⋅m⋅kg⁻¹), unit conversion factor.

```Julia
julia> specificmagnetization(EMU,Metric) # Wb⋅m⋅g⋅Mx⁻¹⋅cm⁻¹⋅kg⁻¹
$(specificmagnetization(EMU,Metric))

julia> specificmagnetization(ESU,Metric) # Wb⋅m⋅g⋅statWb⁻¹⋅cm⁻¹⋅kg⁻¹
$(specificmagnetization(ESU,Metric))

julia> specificmagnetization(Metric,SI2019) # Wb⋅Wb⁻¹
$(specificmagnetization(Metric,SI2019))
```
""" specificmagnetization

@doc """
$(convertext(:demagnetizingfactor,"1/susceptibility(U,S)"))

Quantitiy of `demagnetizingfactor` (dimensionless), unit conversion factor.

```Julia
julia> demagnetizingfactor(EMU,Metric)
$(demagnetizingfactor(EMU,Metric))

julia> demagnetizingfactor(ESU,Metric)
$(demagnetizingfactor(ESU,Metric))

julia> demagnetizingfactor(Metric,SI2019)
$(demagnetizingfactor(Metric,SI2019))
```
""" demagnetizingfactor

# Gyrator-capacitor model alternative:
#@pure magneticfluxrate(U::UnitSystem,S::UnitSystem) = unit(magneticflux(U,S)/time(U,S))
#@pure magneticfluxratedensity(U::UnitSystem,S::UnitSystem) = unit(magneticfluxrate(U,S)/area(U,S))
#@pure magneticresistance(U::UnitSystem,S::UnitSystem) = unit(magneticpotential(U,S)/magneticfluxrate(U,S)) # not reluctance
