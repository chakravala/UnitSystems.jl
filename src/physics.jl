
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

@pure avogadro(U::UnitSystem,C::Coupling=universe(U)) = molarmass(U,C)*electronunit(C)/electronmass(U,C)
@pure dalton(U::UnitSystem,C::Coupling=universe(U)) = electronmass(U,C)/electronunit(C)
@pure protonmass(U::UnitSystem,C::Coupling=universe(U)) = protonelectron(C)*electronmass(U,C)
@pure gaussgravitation(U::UnitSystem,C::Coupling=universe(U)) = sqrt(normal(gravitation(IAU)))*radian(U)/day(U)
@pure einstein(U::UnitSystem,C::Coupling=universe(U)) = two(U)^2*tau(U)*gravitation(U,C)/lightspeed(U,C)^4
#@pure einstein2(U::UnitSystem,C::Coupling=universe(U)) = two(U)^2*tau(U)*gravitation(U,C)/lightspeed(U,C)^2
@pure molargas(U::UnitSystem,C::Coupling=universe(U)) = boltzmann(U,C)*avogadro(U,C)
@pure stefan(U::UnitSystem,C::Coupling=universe(U)) = tau(U)^5/two(U)^4*boltzmann(U,C)^4/(three(U)*five(U)*planck(U,C)^3*lightspeed(U,C)^2)
@pure radiationdensity(U::UnitSystem,C::Coupling=universe(U)) = two(U)^2*stefan(U,C)/lightspeed(U,C)
@pure vacuumpermittivity(U::UnitSystem,C::Coupling=universe(U)) = inv(vacuumpermeability(U,C)*(lightspeed(U,C)*lorentz(U))^2)
@pure electrostatic(U::UnitSystem,C::Coupling=universe(U)) = rationalization(U)/(two(U)*tau(U))/vacuumpermittivity(U,C)
@pure biotsavart(U::UnitSystem,C::Coupling=universe(U)) = vacuumpermeability(U,C)*lorentz(U)*(rationalization(U)/(two(U)*tau(U)))
@pure magnetostatic(U::UnitSystem) = lorentz(U)*biotsavart(U)
@pure vacuumimpedance(U::UnitSystem,C::Coupling=universe(U)) = vacuumpermeability(U,C)*lightspeed(U,C)*rationalization(U)*lorentz(U)^2
@pure faraday(U::UnitSystem,C::Coupling=universe(U)) = elementarycharge(U,C)*avogadro(U,C)
@pure josephson(U::UnitSystem,C::Coupling=universe(U)) = two(U)*elementarycharge(U,C)*lorentz(U)/planck(U,C)
@pure magneticfluxquantum(U::UnitSystem,C::Coupling=universe(U)) = inv(josephson(U,C))
@pure klitzing(U::UnitSystem,C::Coupling=universe(U)) = planck(U,C)/elementarycharge(U,C)^2
@pure conductancequantum(U::UnitSystem,C::Coupling=universe(U)) = two(U)*elementarycharge(U,C)^2/planck(U,C)
@pure hartree(U::UnitSystem,C::Coupling=universe(U)) = electronmass(U,C)/gravity(U)*(lightspeed(U,C)*finestructure(C))^2
@pure rydberg(U::UnitSystem,C::Coupling=universe(U)) = hartree(U,C)/(two(U)*planck(U,C))/lightspeed(U,C)
@pure bohr(U::UnitSystem,C::Coupling=universe(U)) = planckreduced(U,C)*gravity(U)/electronmass(U,C)/lightspeed(U,C)/finestructure(C)
#@pure bohrreduced(U::UnitSystem,C::Coupling=universe(U)) = bohr(U,C)*(one(U)+inv(protonelectron(C)))
@pure electronradius(U::UnitSystem,C::Coupling=universe(U)) = finestructure(C)*planckreduced(U,C)*gravity(U)/electronmass(U,C)/lightspeed(U,C)
@pure magneton(U::UnitSystem,C::Coupling=universe(U)) = elementarycharge(U,C)*planckreduced(U,C)*lorentz(U)/(two(U)*electronmass(U,C))

@pure hyperfine(U::UnitSystem) = frequency(ΔνCs,U,Metric)
@pure hubble(U::UnitSystem) = time(one(U),Hubble,U)
@pure cosmological(U::UnitSystem,C::Coupling=universe(U)) = three(U)*darkenergydensity(C)*(hubble(U)/lightspeed(U,C))^2
@pure loschmidt(U::UnitSystem,P=atmosphere(U),T=T₀*temperature(SI2019,U)) = P/T/boltzmann(U)
@pure amagat(U::UnitSystem) = loschmidt(U)/avogadro(U)
@pure wienwavelength(U::UnitSystem) = planck(U)*lightspeed(U)/boltzmann(U)/Constant(4.965114231744276303)
@pure wienfrequency(U::UnitSystem) = Constant(2.821439372122078893)*boltzmann(U)/planck(U)
@pure eddington(U::UnitSystem) = mass(one(U),U,Cosmological)
@pure solarmass(U::UnitSystem) = mass(GM☉/G,U,Metric)
@pure earthmass(U::UnitSystem) = mass(GME/G,U,Metric)
@pure jupitermass(U::UnitSystem) = mass(GMJ/G,U,Metric)
@pure lunarmass(U::UnitSystem) = earthmass(U)/μE☾
@pure mechanicalheat(U::UnitSystem) = molargas(U)*normal(calorie(Metric)/molargas(Metric))
@pure gaussianyear(U::UnitSystem) = turn(U)/gaussgravitation(U)
@pure siderealyear(U::UnitSystem) = gaussianyear(U)/normal(sqrt(solarmass(IAU)+earthmass(IAU)+lunarmass(IAU)))
@pure gaussianmonth(U::UnitSystem) = tau(U)*sqrt(LD^3/GME)*time(Metric,U)
@pure siderealmonth(U::UnitSystem) = gaussianmonth(U)/normal(sqrt(earthmass(IAUE)+lunarmass(IAUE)))
@pure synodicmonth(U::UnitSystem) = inv(inv(siderealmonth(U))-inv(siderealyear(U)))
@pure jovianyear(U::UnitSystem) = day(U)*sqrt(normal(jupiterdistance(U)^3/solarmass(U)/gravitation(U)))*turn(U)/radian(U)/normal(sqrt(solarmass(IAU)+jupitermass(IAU)))

include("derived.jl")

include("physicsdocs.jl")
