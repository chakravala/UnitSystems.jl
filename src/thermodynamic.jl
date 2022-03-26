
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
    T/K == (9/5)*T

Converts temperature `T` from Kelvin to degrees Rankine (°R).
""" K

@doc """
    T/°R == (5/9)*T

Converts temperature `T` from degrees Rankine to Kelvin (K).
""" °R

"""
    moles(N::Real,U::UnitSystem=Metric) = N/avogadro(U)

Converts the number of molecules `N` to number of moles (mol).
"""
@pure moles(N::Number,U::UnitSystem=Metric) = N/avogadro(U)

"""
    molecules(n::Real,U::UnitSystem=Metric) = n*avogadro(U)

Converts the number of moles `n` to number of molecules (dimensionless).
"""
@pure molecules(n::Number,U::UnitSystem=Metric) = n*avogadro(U)

# thermodynamics

@pure temperature(U::UnitSystem,S::UnitSystem) = unit((boltzmann(U)*electronmass(S)*lightspeed(S)^2*gravity(U))/(boltzmann(S)*electronmass(U)*lightspeed(U)^2*gravity(S)))
@pure entropy(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/temperature(U,S))
@pure specificentropy(U::UnitSystem,S::UnitSystem) = unit(specificenergy(U,S)/temperature(U,S))
@pure volumeheatcapacity(U::UnitSystem,S::UnitSystem) = unit(entropy(U,S)/volume(U,S))
@pure thermalconductivity(U::UnitSystem,S::UnitSystem) = unit(force(U,S)/time(U,S)/temperature(U,S))
@pure thermalconductance(U::UnitSystem,S::UnitSystem) = unit(thermalconductivity(U,S)*length(U,S))
@pure thermalresistance(U::UnitSystem,S::UnitSystem) = thermalconductance(S,U)
@pure thermalexpansion(U::UnitSystem,S::UnitSystem) = temperature(S,U)
@pure lapserate(U::UnitSystem,S::UnitSystem) = unit(temperature(U,S)/length(U,S))

# molar

@pure molality(U::UnitSystem,S::UnitSystem) = molarmass(S,U)
@pure mole(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)*molality(U,S))
@pure molarity(U::UnitSystem,S::UnitSystem) = unit(mole(U,S)/volume(U,S))
@pure molarvolume(U::UnitSystem,S::UnitSystem) = unit(volume(U,S)/mole(U,S))
@pure molarentropy(U::UnitSystem,S::UnitSystem) = unit(entropy(U,S)/mole(U,S))
@pure molarenergy(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/mole(U,S))
@pure molarconductivity(U::UnitSystem,S::UnitSystem) = unit(conductivity(U,S)*area(U,S)/mole(U,S))
@pure molarsusceptibility(U::UnitSystem,S::UnitSystem) = unit(specificsusceptibility(U,S)*molarmass(U,S))
@pure catalysis(U::UnitSystem,S::UnitSystem) = unit(mole(U,S)/time(U,S))
@pure specificity(U::UnitSystem,S::UnitSystem) = unit(volume(U,S)/mole(U,S)/time(U,S))

# photometrics

@pure luminousflux(U::UnitSystem,S::UnitSystem) = unit(frequency(U,S)^2*(luminousefficacy(S)*planckreduced(S))/(luminousefficacy(U)*planckreduced(U)))
@pure luminance(U::UnitSystem,S::UnitSystem) = unit(luminousflux(U,S)/area(U,S))
@pure luminousenergy(U::UnitSystem,S::UnitSystem) = unit(frequency(U,S)*(luminousefficacy(S)*planckreduced(S))/(luminousefficacy(U)*planckreduced(U)))
@pure luminousexposure(U::UnitSystem,S::UnitSystem) = unit(luminance(U,S)*time(U,S))

include("thermodynamicdocs.jl")
