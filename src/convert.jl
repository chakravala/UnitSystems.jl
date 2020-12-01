
#   This file is part of UnitSystems.jl. It is licensed under the MIT license
#   UnitSystems Copyright (C) 2020 Michael Reed

@doc """
    kilograms(m::Real) = $(slug)m

Converts mass `m` from slugs to kilogram (kg).
""" kilograms, slug
@pure kilograms(m::Real,U::UnitSystem=English) = mass(m,Metric,U)

"""
    slugs(m::Real) = $(1/slug)m

Converts mass `m` from kilograms to slugs (slug).
"""
@pure slugs(m::Real,U::UnitSystem=Metric) = mass(m,English,U)

"""
    feet(d) = $(1/ft)d

Converts distance `d` from meters to feet (ft).
"""
@pure feet(d,U::UnitSystem=Metric) = length(d,English,U)

@doc """
    meters(d) = $(ft)d

Converts distance `d` from feet to meters (m).
""" meters, ft
@pure meters(d,U::UnitSystem=English) = length(d,Metric,U)

@doc """
    rankine*T == (9/5)*T

Converts temperature `T` from Kelvin to degrees Rankine (°R).
""" rankine

@doc """
    kelvin*T == (5/9)*T

Converts temperature `T` from degrees Rankine to Kelvin (K).
""" kelvin

"""
    moles(N::Real,U::UnitSystem=Metric) = N/avogadro(U)

Converts the number of molecules `N` to number of moles (mol).
"""
@pure moles(N::Real,U::UnitSystem=Metric) = N/avogadro(U)

"""
    molecules(n::Real,U::UnitSystem=Metric) = n*avogadro(U)

Converts the number of moles `n` to number of molecules (dimensionless).
"""
@pure molecules(n::Real,U::UnitSystem=Metric) = n*avogadro(U)

# CGS to SI
@pure length(::UnitSystem{kB,ħ,𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,100𝘤} where {kB,ħ}) = 100
@pure time(::UnitSystem{kB,ħ,𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,100𝘤} where {kB,ħ}) = 1
@pure temperature(::UnitSystem{kB,ħ,𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,100𝘤} where {kB,ħ}) = 1

# SI to CGS
@pure length(::UnitSystem{kB,ħ,100𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = 0.01
@pure time(::UnitSystem{kB,ħ,100𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = 1
@pure temperature(::UnitSystem{kB,ħ,100𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = 1

# IAU to SI
@pure length(::UnitSystem{kB,ħ,𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ}) = 1/au
@pure time(::UnitSystem{kB,ħ,𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ}) = day
@pure temperature(::UnitSystem{kB,ħ,𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ}) = 1

# SI to IAU
@pure length(::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = au
@pure time(::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = 1/day
@pure temperature(::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = 1

# SI to English
@pure length(::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = ft
@pure time(::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = 1
@pure temperature(::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = rankine

# English to SI
@pure length(::UnitSystem{kB,ħ,𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = 1/ft
@pure time(::UnitSystem{kB,ħ,𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = 1
@pure temperature(::UnitSystem{kB,ħ,𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = kelvin

# CGS to English
@pure length(::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},U::UnitSystem{kB,ħ,100𝘤} where {kB,ħ}) = 100ft
@pure time(::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},U::UnitSystem{kB,ħ,100𝘤} where {kB,ħ}) = 1
@pure temperature(::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},U::UnitSystem{kB,ħ,100𝘤} where {kB,ħ}) = rankine

# English to CGS
@pure length(::UnitSystem{kB,ħ,100𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = 0.01ft
@pure time(::UnitSystem{kB,ħ,100𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = 1
@pure temperature(::UnitSystem{kB,ħ,100𝘤} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = kelvin

# IAU to English
@pure length(::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},U::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ}) = ft/au
@pure time(::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},U::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ}) = day
@pure temperature(::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},U::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ}) = rankine

# English to IAU
@pure length(::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = au*ft
@pure time(::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = 1/day
@pure temperature(::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ},U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = kelvin

# to PlanckGauss
@pure length(::typeof(PlanckGauss),U::UnitSystem) = sqrt(planckreduced(U)*newton(U)/lightspeed(U)^3)
@pure time(::typeof(PlanckGauss),U::UnitSystem) = sqrt(planckreduced(U)*newton(U)/lightspeed(U)^5)
@pure temperature(::typeof(PlanckGauss),U::UnitSystem) = sqrt(planckreduced(U)*lightspeed(U)^5/newton(U)/boltzmann(U)^2)

# to Stoney
@pure length(::typeof(Stoney),U::UnitSystem) = sqrt(newton(U)*coulomb(U)*charge(U)^2/lightspeed(U)^4)
@pure time(::typeof(Stoney),U::UnitSystem) = sqrt(planckreduced(U)*newton(U)^2/lightspeed(U)^6)

# to Hartree
@pure length(::typeof(Hartree),U::UnitSystem) = 4π*permeability(U)*planckreduced(U)^2/electronmass(U)/charge(U)^2
@pure time(::typeof(Hartree),U::UnitSystem) = (4π*permeability(U)/charge(U)^2)^2*planckreduced(U)/electronmass(U)

# to NaturalGauss
@pure length(::typeof(NaturalGauss),U::UnitSystem) = planckreduced(U)/electronmass(U)/lightspeed(U)
@pure time(::typeof(NaturalGauss),U::UnitSystem) = planckreduced(U)/electronmass(U)/lightspeed(U)^2

# to QCDx
@pure length(::UnitSystem{1,1,1,μ₀,1/μₚₑ} where μ₀,U::UnitSystem) = planckreduced(U)/protonmass(U)/lightspeed(U)
@pure time(::UnitSystem{1,1,1,μ₀,1/μₚₑ} where μ₀,U::UnitSystem) = planckreduced(U)/protonmass(U)/lightspeed(U)^2
