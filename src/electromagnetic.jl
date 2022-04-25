
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

@pure charge(U::UnitSystem,S::UnitSystem) = unit(sqrt((turn(S)/turn(U))*(planckreduced(S)*vacuumpermeability(U)*lightspeed(U)*rationalization(U)*lorentz(U)^2)/(planckreduced(U)*vacuumpermeability(S)*lightspeed(S)*rationalization(S)*lorentz(S)^2)))

@pure current(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)/time(U,S))

@pure electricpotential(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/charge(U,S))
const voltage = electricpotential
@pure capacitance(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)/electricpotential(U,S))
@pure resistance(U::UnitSystem,S::UnitSystem) = unit(electricpotential(U,S)/current(U,S))
@pure conductance(U::UnitSystem,S::UnitSystem) = unit(current(U,S)/electricpotential(U,S))
@pure magneticflux(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/lorentz(U,S)/current(U,S))
@pure magneticfluxdensity(U::UnitSystem,S::UnitSystem) = unit(magneticflux(U,S)/area(U,S))
@pure inductance(U::UnitSystem,S::UnitSystem) = unit(magneticflux(U,S)/current(U,S)*lorentz(U,S))

# electromagnetics

@pure electricfluxdensity(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)*rationalization(U,S)/area(U,S))
@pure chargedensity(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)/volume(U,S))
@pure currentdensity(U::UnitSystem,S::UnitSystem) = unit(current(U,S)/area(U,S))
@pure conductivity(U::UnitSystem,S::UnitSystem) = unit(conductance(U,S)/length(U,S))
@pure permittivity(U::UnitSystem,S::UnitSystem) = unit(capacitance(U,S)*rationalization(U,S)/length(U,S))
@pure electricfield(U::UnitSystem,S::UnitSystem) = unit(electricpotential(U,S)/length(U,S))
@pure magneticfield(U::UnitSystem,S::UnitSystem) = unit(current(U,S)*rationalization(U,S)*lorentz(U,S)/length(U,S))
@pure exposure(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)/mass(U,S))
@pure resistivity(U::UnitSystem,S::UnitSystem) = unit(resistance(U,S)*length(U,S))
@pure linearchargedensity(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)/length(U,S))
@pure magneticdipolemoment(U::UnitSystem,S::UnitSystem) = unit(current(U,S)*lorentz(U,S)*area(U,S)/gravity(U,S)/angle(U,S))
@pure mobility(U::UnitSystem,S::UnitSystem) = unit(length(U,S)*speed(U,S)*electricpotential(U,S))
@pure reluctance(U::UnitSystem,S::UnitSystem) = unit(rationalization(U,S)*lorentz(U,S)^2/inductance(U,S))
@pure vectorpotential(U::UnitSystem,S::UnitSystem) = unit(magneticflux(U,S)/length(U,S))
@pure magneticmoment(U::UnitSystem,S::UnitSystem) = unit(magneticflux(U,S)*length(U,S))
#@pure rigidity(U::UnitSystem,S::UnitSystem) = unit(magneticfluxdensity(U,S)*length(U,S))
@pure susceptibility(U::UnitSystem,S::UnitSystem) = unit(rationalization(S,U))

# WARNING unchecked: rigidity, magneticmoment, vectorpotential linearchargedensity

# CGS extra: polarizability, permeance, magnetic-current? + density, magneticresistance

@pure electricflux(U::UnitSystem,S::UnitSystem) = unit(electricpotential(U,S)*length(U,S))
@pure electricdipolemoment(U::UnitSystem,S::UnitSystem) = unit(charge(U,S)*length(U,S))
@pure magneticpotential(U::UnitSystem,S::UnitSystem) = unit(magneticflux(U,S)*reluctance(U,S))
@pure polestrength(U::UnitSystem,S::UnitSystem) = unit(magneticdipolemoment(U,S)/length(U,S))
@pure permeance(U::UnitSystem,S::UnitSystem) = reluctance(S,U)
@pure specificsusceptibility(U::UnitSystem,S::UnitSystem) = unit(magneticdipolemoment(U,S)/magneticfield(U,S)/mass(U,S))
#@pure magnetizability(U::UnitSystem,S::UnitSystem) = unit(magneticmoment(U,S)/magneticfluxdensity(U,S))
@pure electricpolarizability(U::UnitSystem,S::UnitSystem) = unit(electricdipolemoment(U,S)/electricfield(U,S))
@pure magneticpolarizability(U::UnitSystem,S::UnitSystem) = unit(magneticdipolemoment(U,S)/magneticfield(U,S))
#@pure magnetization(U::UnitSystem,S::UnitSystem) = unit(magneticmoment(U,S)/volume(U,S))

# specificmagnetization, mass magnetization = 1?
# magneticfluxdensity(Metric,EMU)/density(Metric,EMU)
@pure specificmagnetization(U::UnitSystem,S::UnitSystem) = unit(magneticmoment(S,U)/mass(S,U))
@pure demagnetizingfactor(U::UnitSystem,S::UnitSystem) = unit(rationalization(U,S))

# Gyrator-capacitor model alternative:
#@pure magneticfluxrate(U::UnitSystem,S::UnitSystem) = unit(magneticflux(U,S)/time(U,S))
#@pure magneticfluxratedensity(U::UnitSystem,S::UnitSystem) = unit(magneticfluxrate(U,S)/area(U,S))
#@pure magneticresistance(U::UnitSystem,S::UnitSystem) = unit(magneticpotential(U,S)/magneticfluxrate(U,S)) # not reluctance

include("electromagneticdocs.jl")
