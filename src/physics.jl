
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
@pure atomicmass(U::UnitSystem,C::Coupling=universe(U)) = electronmass(U,C)/electronunit(C)
@pure protonmass(U::UnitSystem,C::Coupling=universe(U)) = protonelectron(C)*electronmass(U,C)
@pure einstein(U::UnitSystem,C::Coupling=universe(U)) = two(U)*sphere(U)*newton(U,C)/lightspeed(U,C)^4
@pure einstein2(U::UnitSystem,C::Coupling=universe(U)) = two(U)*sphere(U)*newton(U,C)/lightspeed(U,C)^2
@pure universalgas(U::UnitSystem,C::Coupling=universe(U)) = boltzmann(U,C)*avogadro(U,C)
@pure stefan(U::UnitSystem,C::Coupling=universe(U)) = normal(turn(U))^4/two(U)^5*sphere(U)*boltzmann(U,C)^4/(three(U)*five(U)*planck(U,C)^3*lightspeed(U,C)^2)
@pure radiationdensity(U::UnitSystem,C::Coupling=universe(U)) = two(U)^2*stefan(U,C)/lightspeed(U,C)
@pure vacuumpermittivity(U::UnitSystem,C::Coupling=universe(U)) = inv(vacuumpermeability(U,C)*(lightspeed(U,C)*lorentz(U))^2)
@pure coulomb(U::UnitSystem,C::Coupling=universe(U)) = rationalization(U)/sphere(U)/vacuumpermittivity(U,C)
@pure biotsavart(U::UnitSystem,C::Coupling=universe(U)) = vacuumpermeability(U,C)*lorentz(U)*(rationalization(U)/sphere(U))
@pure ampere(U::UnitSystem) = lorentz(U)*biotsavart(U)
@pure vacuumimpedance(U::UnitSystem,C::Coupling=universe(U)) = vacuumpermeability(U,C)*lightspeed(U,C)*rationalization(U)*lorentz(U)^2
@pure faraday(U::UnitSystem,C::Coupling=universe(U)) = elementarycharge(U,C)*avogadro(U,C)
@pure josephson(U::UnitSystem,C::Coupling=universe(U)) = two(U)*elementarycharge(U,C)*lorentz(U)/planck(U,C)
@pure magneticfluxquantum(U::UnitSystem,C::Coupling=universe(U)) = inv(josephson(U,C))
@pure klitzing(U::UnitSystem,C::Coupling=universe(U)) = planck(U,C)/elementarycharge(U,C)^2
@pure conductancequantum(U::UnitSystem,C::Coupling=universe(U)) = two(U)*elementarycharge(U,C)^2/planck(U,C)
@pure hartree(U::UnitSystem,C::Coupling=universe(U)) = electronmass(U,C)*(lightspeed(U,C)*finestructure(C))^2
@pure rydberg(U::UnitSystem,C::Coupling=universe(U)) = hartree(U,C)/(two(U)*planck(U,C))/lightspeed(U,C)
@pure bohr(U::UnitSystem,C::Coupling=universe(U)) = planckreduced(U,C)/electronmass(U,C)/lightspeed(U,C)/finestructure(C)
@pure bohrreduced(U::UnitSystem,C::Coupling=universe(U)) = bohr(U,C)*(one(U)+inv(protonelectron(C)))
@pure electronradius(U::UnitSystem,C::Coupling=universe(U)) = finestructure(C)*planckreduced(U,C)/electronmass(U,C)/lightspeed(U,C)
@pure magneton(U::UnitSystem,C::Coupling=universe(U)) = elementarycharge(U,C)*planckreduced(U,C)*lorentz(U)/(two(U)*electronmass(U,C))

@pure hyperfine(U::UnitSystem) = frequency(ΔνCs,U,Metric)
@pure hubble(U::UnitSystem) = time(one(U),Hubble,U)
@pure cosmological(U::UnitSystem,C::Coupling=universe(U)) = three(U)*darkenergydensity(C)*(hubble(U)/lightspeed(U,C))^2
@pure standardgravity(U::UnitSystem) = acceleration(g₀,U,Metric)
@pure standardpressure(U::UnitSystem) = pressure(atm,U,Metric)
@pure standardtemperature(U::UnitSystem) = temperature(Tₛ,U,Metric)
@pure solarmass(U::UnitSystem) = mass(GM☉/G,U,Metric)
@pure earthmass(U::UnitSystem) = mass(GME/G,U,Metric)
@pure jupitermass(U::UnitSystem) = mass(GMJ/G,U,Metric)
@pure lunarmass(U::UnitSystem) = earthmass(U)/μE☾
@pure astronomicalunit(U::UnitSystem) = length(𝟏,U,IAU)
@pure lunardistance(U::UnitSystem) = length(LD,U,Metric)
@pure mile(U::UnitSystem) = length(two(U)^5*three(U)*five(U)*eleven(U),U,English)
@pure clarkemile(U::UnitSystem) = length(two(U)^6*five(U)*nineteen(U),U,English)
@pure nauticalmile(U::UnitSystem) = length(two(U)^4*five(U)^5/three(U)^3,U,Metric)
@pure meancalorie(U::UnitSystem) = energy(two(U)^2*five(U)*three(U)^2/fourtythree(U),U,InternationalMean)
@pure kilocalorie(U::UnitSystem) = energy(two(U)^5*five(U)^4*three(U)^2/fourtythree(U),U,International)
@pure calorie(U::UnitSystem) = kilocalorie(U)/(two(U)*five(U))^3

@pure thermalunit(U::UnitSystem) = mass(temperature(kilocalorie(U),Metric,English),Metric,English)
@pure tonsrefrigeration(U::UnitSystem) = frequency(two(U)*five(U)/three(U),U,Metric)*thermalunit(U)
@pure boilerhorsepower(U::UnitSystem) = frequency(Constant(1339)/(two(U)^4*three(U)^2),U,Metric)*thermalunit(U)
# thermalconductivity_water(British) ≈ 0.5778
@pure thermalconductivity_water(U::UnitSystem) = thermalconductivity((two(U)^2*three(U)*five(U))^2/thermalunit(U),U,Metric)

@pure gallon(U::UnitSystem) = volume(seven(U)*eleven(U)/two(U)^2,U,English)
@pure litre(U::UnitSystem) = volume(inv((two(U)*five(U))^3),U,Metric)
@pure horsepower(U::UnitSystem) = power(two(U)*five(U)^2*eleven(U),U,British)
@pure horsepowerwatt(U::UnitSystem) = power(two(U)^4*three(U)^3/five(U)*normal(twopi(U)),U,British)
@pure horsepowermetric(U::UnitSystem) = power(three(U)*five(U)^2,U,GravitationalMetric)
@pure electricalhorsepower(U::UnitSystem) = power(Constant(746),U,Metric)
@pure inchmercury(U::UnitSystem) = pressure(inHg,U,Metric)
@pure torr(U::UnitSystem) = pressure(atm/(two(U)^3*five(U)*nineteen(U)),U,Metric)

@pure second(U::UnitSystem) = time(one(U),U,Metric)
@pure minute(U::UnitSystem) = two(U)^2*three(U)*five(U)*second(U)
@pure hour(U::UnitSystem) = two(U)^2*three(U)*five(U)*minute(U)
@pure day(U::UnitSystem) = two(U)^3*three(U)*hour(U)
@pure year(U::UnitSystem) = aⱼ*day(U)
@pure gaussianyear(U::UnitSystem) = (τ/k)*day(U)
@pure siderealyear(U::UnitSystem) = (τ/k/√(solarmass(IAU)+earthmass(IAU)+lunarmass(IAU)))*day(U)
@pure lightyear(U::UnitSystem) = year(U)*lightspeed(U)
@pure parsec(U::UnitSystem) = astronomicalunit(U)*two(U)^7*three(U)^4*five(U)^3/turn(U)

include("physicsdocs.jl")

#=@pure kilogram(U::UnitSystem) = mass(Metric,U)
@pure slug(U::UnitSystem) = mass(English,U)

@pure meter(U::UnitSystem) = length(Metric,U)
@pure foot(U::UnitSystem) = length(English,U)=#

#rankine, kelvin, moles/molecules
#add gravitional units of weight??
