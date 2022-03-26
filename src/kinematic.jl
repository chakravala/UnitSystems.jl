
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

"""
    kilograms(m::Real) = $(slug)m

Converts mass `m` from slugs to kilogram (kg).
"""
@pure kilograms(m::Number,U::UnitSystem=English) = mass(m,Metric,U)

"""
    slugs(m::Real) = $(1/slug)m

Converts mass `m` from kilograms to slugs (slug).
"""
@pure slugs(m::Number,U::UnitSystem=Metric) = mass(m,English,U)

"""
    feet(d) = $(1/ft)d

Converts distance `d` from meters to feet (ft).
"""
@pure feet(d,U::UnitSystem=Metric) = length(d,English,U)

"""
    meters(d) = $(ft)d

Converts distance `d` from feet to meters (m).
"""
@pure meters(d,U::UnitSystem=English) = length(d,Metric,U)

# special exact cases

# IAU to SI
@pure length(U::UnitSystem{kB,ħ,𝘤} where {kB,ħ},S::UnitSystem{kB,ħ,DAY*𝘤/au} where {kB,ħ}) = length(U,S,1/au)
@pure time(U::UnitSystem{kB,ħ,𝘤} where {kB,ħ},S::UnitSystem{kB,ħ,DAY*𝘤/au} where {kB,ħ}) = time(U,S,1/DAY)

# SI to IAU
@pure length(U::UnitSystem{kB,ħ,DAY*𝘤/au} where {kB,ħ},S::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = length(U,S,au)
@pure time(U::UnitSystem{kB,ħ,DAY*𝘤/au} where {kB,ħ},S::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = time(U,S,DAY)

# IAU to CGS
@pure length(::UnitSystem{kB,ħ,100𝘤} where {kB,ħ},S::UnitSystem{kB,ħ,DAY*𝘤/au} where {kB,ħ}) = length(U,S,1/au)
@pure time(::UnitSystem{kB,ħ,100𝘤} where {kB,ħ},S::UnitSystem{kB,ħ,DAY*𝘤/au} where {kB,ħ}) = time(U,S,1/DAY)

# CGS to IAU
@pure length(U::UnitSystem{kB,ħ,DAY*𝘤/au} where {kB,ħ},S::UnitSystem{kB,ħ,100𝘤} where {kB,ħ}) = length(U,S,au)
@pure time(U::UnitSystem{kB,ħ,DAY*𝘤/au} where {kB,ħ},S::UnitSystem{kB,ħ,100𝘤} where {kB,ħ}) = time(U,S,DAY)

# IAU to English
@pure length(U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},S::UnitSystem{kB,ħ,DAY*𝘤/au} where {kB,ħ}) = length(U,S,ft/au)
@pure time(U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},S::UnitSystem{kB,ħ,DAY*𝘤/au} where {kB,ħ}) = time(U,S,1/DAY)

# English to IAU
@pure length(U::UnitSystem{kB,ħ,DAY*𝘤/au} where {kB,ħ},S::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = length(U,S,au/ft)
@pure time(U::UnitSystem{kB,ħ,DAY*𝘤/au} where {kB,ħ},S::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = time(U,S,DAY)

# base unit conversions

@pure length(U::UnitSystem,S::UnitSystem,l=1) = isquantity(U,S) ? evaldim(length)(U,S) : unit((turn(S)/turn(U))*(planckreduced(S)*electronmass(U)*lightspeed(U)*gravity(S))/(planckreduced(U)*electronmass(S)*lightspeed(S)*gravity(U)),l)
@doc """
$(convertext(:length,"planck(U,S)/mass(U,S)/speed(U,S)"))

Extent of one-dimensional shape or `length` (m), unit conversion factor.

```Julia
julia> length(CGS,Metric) # m⋅cm⁻¹
$(length(CGS,Metric))

julia> length(IAU,Metric) # m⋅au⁻¹
$(length(IAU,Metric))

julia> length(English,Metric) # m⋅ft⁻¹
$(length(English,Metric))

julia> length(Survey,English) # ft⋅ftUS⁻¹
$(length(Survey,English))

julia> length(PlanckGauss,Metric) # m⋅ℓP⁻¹
$(length(PlanckGauss,Metric))
```
""" length, ft, ftUS, L

@pure time(U::UnitSystem,S::UnitSystem,t=1) = isquantity(U,S) ? evaldim(time)(U,S) : unit(length(U,S)/lightspeed(U,S),1)
@doc """
$(convertext(:time,"length(U,S)/speed(U,S)"))

Dimension along which events are ordered or `time` (s), unit conversion factor.

```Julia
julia> time(IAU,Metric) # s⋅day⁻¹
$(time(IAU,Metric))

julia> time(PlanckGauss,Metric) # s⋅tP⁻¹
$(time(PlanckGauss,Metric))
```
""" time(::UnitSystem,::UnitSystem), T

# spacetime

@pure area(U::UnitSystem,S::UnitSystem) = unit(length(U,S)^2)
@pure volume(U::UnitSystem,S::UnitSystem) = unit(length(U,S)^3)
@pure wavenumber(U::UnitSystem,S::UnitSystem) = unit(length(S,U))
@pure angularwavenumber(U::UnitSystem,S::UnitSystem) = unit(angle(U,S)*length(S,U))
@pure fuelefficiency(U::UnitSystem,S::UnitSystem) = area(S,U)
@pure frequency(U::UnitSystem,S::UnitSystem) = time(S,U)
@pure angularfrequency(U::UnitSystem,S::UnitSystem) = unit(angle(U,S)*time(S,U))
@pure frequencydrift(U::UnitSystem,S::UnitSystem) = unit(time(S,U)^2)
@pure speed(U::UnitSystem,S::UnitSystem) = lightspeed(U,S)
@pure acceleration(U::UnitSystem,S::UnitSystem) = unit(speed(U,S)/time(U,S))
@pure jerk(U::UnitSystem,S::UnitSystem) = unit(speed(U,S)/time(U,S)^2)
@pure snap(U::UnitSystem,S::UnitSystem) = unit(speed(U,S)/time(U,S)^3)
@pure crackle(U::UnitSystem,S::UnitSystem) = unit(speed(U,S)/time(U,S)^4)
@pure pop(U::UnitSystem,S::UnitSystem) = unit(speed(U,S)/time(U,S)^5)
@pure volumeflow(U::UnitSystem,S::UnitSystem) = unit(area(U,S)*speed(U,S))
@pure specificenergy(U::UnitSystem,S::UnitSystem) = unit(speed(U,S)^2/gravity(U,S))

# kinematic

@pure inertia(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/gravity(U,S))
@pure energy(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)*specificenergy(U,S))
@pure power(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/time(U,S))
@pure force(U::UnitSystem,S::UnitSystem) = unit(inertia(U,S)*acceleration(U,S))
@pure gforce(U::UnitSystem,S::UnitSystem) = unit(acceleration(U,S)/gravity(U,S))
@pure pressure(U::UnitSystem,S::UnitSystem) = unit(force(U,S)/area(U,S))

# mechanical

@pure impulse(U::UnitSystem,S::UnitSystem) = unit(force(U,S)*time(U,S))
@pure momentum(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)*speed(U,S))
@pure angularmomentum(U::UnitSystem,S::UnitSystem) = unit(momentum(U,S)*length(U,S)*angle(U,S))
@pure yank(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)*jerk(U,S))
@pure areadensity(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/area(U,S))
@pure density(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/volume(U,S))
@pure specificweight(U::UnitSystem,S::UnitSystem) = unit(pressure(U,S)/speed(U,S)^2)
@pure specificvolume(U::UnitSystem,S::UnitSystem) = unit(volume(U,S)/mass(U,S))
@pure action(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)*time(U,S))
@pure stiffness(U::UnitSystem,S::UnitSystem) = unit(force(U,S)/length(U,S))
@pure intensity(U::UnitSystem,S::UnitSystem) = unit(power(U,S)/area(U,S))
@pure diffusivity(U::UnitSystem,S::UnitSystem) = unit(speed(U,S)*length(U,S))
@pure viscosity(U::UnitSystem,S::UnitSystem) = unit(force(U,S)/speed(U,S)/length(U,S))
@pure lineardensity(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/length(U,S))
@pure massflow(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/time(U,S))
@pure spectralflux(U::UnitSystem,S::UnitSystem) = unit(power(U,S)/length(U,S))
@pure powerdensity(U::UnitSystem,S::UnitSystem) = unit(power(U,S)/volume(U,S))
@pure compressibility(U::UnitSystem,S::UnitSystem) = pressure(S,U)
@pure fluence(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/area(U,S))
@pure rotationalinertia(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)*area(U,S))

# acoustic

@pure soundexposure(U::UnitSystem,S::UnitSystem) = unit(time(U,S)*pressure(U,S)^2)
@pure specificimpedance(U::UnitSystem,S::UnitSystem) = unit(pressure(U,S)/speed(U,S))
@pure impedance(U::UnitSystem,S::UnitSystem) = unit(specificimpedance(U,S)/area(U,S))
@pure admittance(U::UnitSystem,S::UnitSystem) = unit(area(U,S)/specificimpedance(U,S))
@pure compliance(U::UnitSystem,S::UnitSystem) = unit(time(U,S)^2/mass(U,S))
@pure inertance(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/length(U,S)^4)

include("kinematicdocs.jl")
