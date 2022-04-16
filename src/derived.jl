
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

@pure deka(U::UnitSystem) = two(U)*five(U)
@pure hecto(U::UnitSystem) = deka(U)^2
@pure kilo(U::UnitSystem) = deka(U)^3
@pure mega(U::UnitSystem) = (Constant(1.0)*kilo(U))^2
@pure giga(U::UnitSystem) = (Constant(1.0)*kilo(U))^3
@pure tera(U::UnitSystem) = (Constant(1.0)*kilo(U))^4
@pure peta(U::UnitSystem) = (Constant(1.0)*kilo(U))^5
@pure exa(U::UnitSystem) = (Constant(1.0)*kilo(U))^6
@pure zetta(U::UnitSystem) = (Constant(1.0)*kilo(U))^7
@pure yotta(U::UnitSystem) = (Constant(1.0)*kilo(U))^8
@pure deci(U::UnitSystem) = inv(deka(U))
@pure centi(U::UnitSystem) = inv(hecto(U))
@pure milli(U::UnitSystem) = inv(kilo(U))
@pure micro(U::UnitSystem) = inv(mega(U))
@pure nano(U::UnitSystem) = inv(giga(U))
@pure pico(U::UnitSystem) = inv(tera(U))
@pure femto(U::UnitSystem) = inv(peta(U))
@pure atto(U::UnitSystem) = inv(exa(U))
@pure zepto(U::UnitSystem) = inv(zetta(U))
@pure yocto(U::UnitSystem) = inv(yotta(U))
@pure kibi(U::UnitSystem) = two(U)^10
@pure mebi(U::UnitSystem) = two(U)^20
@pure gibi(U::UnitSystem) = two(U)^30
@pure tebi(U::UnitSystem) = two(U)^40
@pure pebi(U::UnitSystem) = two(U)^50
@pure exbi(U::UnitSystem) = two(U)^60
@pure zebi(U::UnitSystem) = (Constant(1.0)*two(U))^70
@pure yobi(U::UnitSystem) = (Constant(1.0)*two(U))^80

# angle

@pure radian(U::UnitSystem) = angle(one(U),U,Metric)
@pure steradian(U::UnitSystem) = solidangle(one(U),U,Metric)
@pure degree(U::UnitSystem) = angle(turn(U)/two(U)^3/three(U)^2/five(U),U,Metric)
@pure gradian(U::UnitSystem) = angle(turn(U)/two(U)^4/five(U)^2,U,Metric)
@pure arcminute(U::UnitSystem) = degree(U)/two(U)^2/three(U)/five(U)
@pure arcsecond(U::UnitSystem) = arcminute(U)/two(U)^2/three(U)/five(U)

@pure rpm(U::UnitSystem) = one(U)/minute(U)
#@pure rpd(U::UnitSystem) = turn(U)/day(U)

# time

@pure second(U::UnitSystem) = time(one(U),U,Metric)
@pure minute(U::UnitSystem) = two(U)^2*three(U)*five(U)*second(U)
@pure hour(U::UnitSystem) = two(U)^2*three(U)*five(U)*minute(U)
@pure day(U::UnitSystem) = two(U)^3*three(U)*hour(U)
@pure year(U::UnitSystem) = aⱼ*day(U)
@pure radarmile(U::UnitSystem) = two(U)*nauticalmile(U)/lightspeed(U)

# length

@pure meter(U::UnitSystem) = length(one(U),U,Metric)
@pure earthmeter(U::UnitSystem) = length(one(U),U,Meridian)
#@pure navigationmeter(U::UnitSystem) = greatcircle(U)/two(U)^9/five(U)^7
@pure angstrom(U::UnitSystem) = hecto(U)*pico(U)*meter(U)
@pure foot(U::UnitSystem) = length(one(U),U,English)
@pure inch(U::UnitSystem) = length(one(U),U,IPS)
#@pure rackunit(U::UnitSystem) = foot(U)*seven(U)/two(U)^4/three(U)
@pure yard(U::UnitSystem) = three(U)*foot(U)
@pure surveyfoot(U::UnitSystem) = length(one(U),U,Survey)
@pure statutemile(U::UnitSystem) = length(two(U)^5*three(U)*five(U)*eleven(U),U,Survey)
@pure earthradius(U::UnitSystem) = sqrt(earthmass(U)*gravitation(U)/gforce(U))
#@pure navigationradius(U::UnitSystem) = length(sqrt(earthmass(Metric)*gravitation(Metric)/(gforce(Metric) - gravitation(Metric)*(lunarmass(Metric)/lunardistance(Metric)^2+solarmass(Metric)/astronomicalunit(Metric)^2))),U,Metric)
@pure greatcircle(U::UnitSystem) = normal(turn(U))*earthradius(U)
#@pure greatcircle(U::UnitSystem) = normal(turn(U))*navigationradius(U)
@pure nauticalmile(U::UnitSystem) = length(one(U),U,Nautical)
#@pure navigationmile(U::UnitSystem) = greatcircle(U)/two(U)^5/three(U)^3/five(U)^2
@pure astronomicalunit(U::UnitSystem) = length(𝟏,U,IAU)
@pure lunardistance(U::UnitSystem) = length(LD,U,Metric)
@pure mile(U::UnitSystem) = length(two(U)^5*three(U)*five(U)*eleven(U),U,English)
@pure admiraltymile(U::UnitSystem) = length(two(U)^6*five(U)*nineteen(U),U,English)
@pure meridianmile(U::UnitSystem) = length(two(U)^4*five(U)^5/three(U)^3,U,Metric)
@pure lightyear(U::UnitSystem) = year(U)*lightspeed(U)
@pure parsec(U::UnitSystem) = astronomicalunit(U)*two(U)^7*three(U)^4*five(U)^3/turn(U)

# area

@pure barn(U::UnitSystem) = area((two(U)*five(U))^-28,U,Metric)
@pure hectare(U::UnitSystem) = area(hecto(U)*hecto(U),U,Metric)
@pure acre(U::UnitSystem) = area(two(U)^-7/five(U),U,MPH)
@pure surveyacre(U::UnitSystem) = area(two(U)^3*three(U)^2*five(U)*eleven(U)^2,U,Survey)
#@pure township(U::UnitSystem) = two(U)^9*three(U)^2*five(U)*surveyacre(U)
#@pure footballfield(U::UnitSystem) = area(two(U)^8*three(U)^2*five(U)^2,U,English)

# volume

@pure gallon(U::UnitSystem) = volume(three(U)*seven(U)*eleven(U),U,IPS)
@pure liter(U::UnitSystem) = volume(inv((two(U)*five(U))^3),U,Metric)
@pure quart(U::UnitSystem) = gallon(U)/two(U)^2
@pure pint(U::UnitSystem) = quart(U)/two(U)
@pure cup(U::UnitSystem) = pint(U)/two(U)
@pure fluidounce(U::UnitSystem) = cup(U)/two(U)^3
@pure teaspoon(U::UnitSystem) = five(U)*milli(U)*liter(U)
@pure tablespoon(U::UnitSystem) = three(U)*teaspoon(U)
#@pure oilbarrel(U::UnitSystem) = two(U)*three(U)*seven(U)*gallon(U)

# mass

@pure grain(U::UnitSystem) = milli(U)*pound(U)/seven(U)
@pure gram(U::UnitSystem) = mass(one(U),U,Gauss)
@pure earthgram(U::UnitSystem) = mass(milli(U),U,Meridian)
@pure kilogram(U::UnitSystem) = mass(one(U),U,Metric)
@pure tonne(U::UnitSystem) = mass(kilo(U),U,Metric)
@pure ton(U::UnitSystem) = mass(two(U)*kilo(U),U,English)
@pure pound(U::UnitSystem) = mass(one(U),U,English)
@pure ounce(U::UnitSystem) = mass(two(U)^-4,U,English)
@pure slug(U::UnitSystem) = mass(one(U),U,British)
@pure slinch(U::UnitSystem) = mass(one(U),U,IPS)
@pure hyl(U::UnitSystem) = mass(one(U),U,GravitationalMetric)

# force

@pure dyne(U::UnitSystem) = force(one(U),U,Gauss)
@pure newton(U::UnitSystem) = force(one(U),U,Metric)
@pure poundal(U::UnitSystem) = force(one(U),U,FPS)
@pure poundforce(U::UnitSystem) = force(one(U),U,English)
@pure kilopond(U::UnitSystem) = force(one(U),U,MetricEngineering)

# pressure

@pure pascal(U::UnitSystem) = pressure(one(U),U,Metric)
@pure bar(U::UnitSystem) = pressure(hecto(U)*kilo(U),U,Metric)
@pure barye(U::UnitSystem) = pressure(one(U),U,Gauss)
@pure psi(U::UnitSystem) = pressure(one(U),U,IPS)
@pure technicalatmosphere(U::UnitSystem) = kilopond(U)/(centi(U)*meter(U))^2
@pure atmosphere(U::UnitSystem) = pressure(atm,U,Metric)
@pure inchmercury(U::UnitSystem) = pressure(inHg,U,Metric)
@pure torr(U::UnitSystem) = pressure(atm/(two(U)^3*five(U)*nineteen(U)),U,Metric)

# energy

@pure electronvolt(U::UnitSystem) = elementarycharge(U)*electricpotential(one(U),U,SI2019)
@pure erg(U::UnitSystem) = energy(one(U),U,Gauss)
@pure joule(U::UnitSystem) = energy(one(U),U,Metric)
@pure footpound(U::UnitSystem) = poundforce(U)*foot(U)
@pure meancalorie(U::UnitSystem) = energy(two(U)^2*five(U)*three(U)^2/fourtythree(U),U,InternationalMean)
@pure kilocalorie(U::UnitSystem) = energy(two(U)^5*five(U)^4*three(U)^2/fourtythree(U),U,International)
@pure calorie(U::UnitSystem) = kilocalorie(U)/(two(U)*five(U))^3
@pure earthcalorie(U::UnitSystem) = molaramount(temperature(calorie(U),Metric,Meridian),Metric,Meridian)
@pure thermalunit(U::UnitSystem) = mass(temperature(kilocalorie(U),Metric,English),Metric,English)
@pure tontnt(U::UnitSystem) = giga(U)*calorie(U)
@pure gasgallon(U::UnitSystem) = two(U)*three(U)*nineteen(U)*kilo(U)*thermalunit(U)

# power

@pure watt(U::UnitSystem) = power(one(U),U,Metric)
@pure tonsrefrigeration(U::UnitSystem) = frequency(two(U)*five(U)/three(U),U,Metric)*thermalunit(U)
@pure boilerhorsepower(U::UnitSystem) = frequency(Constant(1339)/(two(U)^4*three(U)^2),U,Metric)*thermalunit(U)
# thermalconductivity_water(British) ≈ 0.5778
@pure thermalconductivity_water(U::UnitSystem) = thermalconductivity((two(U)^2*three(U)*five(U))^2/thermalunit(U),U,Metric)
@pure horsepower(U::UnitSystem) = power(two(U)*five(U)^2*eleven(U),U,British)
@pure horsepowerwatt(U::UnitSystem) = power(two(U)^4*three(U)^3/five(U)*normal(twopi(U)),U,British)
@pure horsepowermetric(U::UnitSystem) = power(three(U)*five(U)^2,U,GravitationalMetric)
@pure electricalhorsepower(U::UnitSystem) = power(Constant(746),U,Metric)

# electromagnetic

@pure coulomb(U::UnitSystem) = charge(one(U),U,Metric)
@pure ampere(U::UnitSystem) = current(one(U),U,Metric)
@pure volt(U::UnitSystem) = electricpotential(one(U),U,Metric)
@pure henry(U::UnitSystem) = inductance(one(U),U,Metric)
@pure ohm(U::UnitSystem) = resistance(one(U),U,Metric)
@pure siemens(U::UnitSystem) = conductance(one(U),U,Metric)
@pure farad(U::UnitSystem) = capacitance(one(U),U,Metric)
@pure weber(U::UnitSystem) = magneticflux(one(U),U,Metric)
@pure tesla(U::UnitSystem) = magneticfluxdensity(one(U),U,Metric)
@pure abcoulomb(U::UnitSystem) = charge(one(U),U,EMU)
@pure abampere(U::UnitSystem) = current(one(U),U,EMU)
@pure abvolt(U::UnitSystem) = electricpotential(one(U),U,EMU)
@pure abhenry(U::UnitSystem) = inductance(one(U),U,EMU)
@pure abohm(U::UnitSystem) = resistance(one(U),U,EMU)
@pure abmho(U::UnitSystem) = conductance(one(U),U,EMU)
@pure abfarad(U::UnitSystem) = capacitance(one(U),U,EMU)
@pure maxwell(U::UnitSystem) = magneticflux(one(U),U,EMU)
@pure gauss(U::UnitSystem) = magneticfluxdensity(one(U),U,EMU)
@pure oersted(U::UnitSystem) = magneticfield(one(U),U,EMU)
@pure gilbert(U::UnitSystem) = abampere(U)/two(U)/turn(U)
@pure statcoulomb(U::UnitSystem) = charge(one(U),U,ESU)
@pure statampere(U::UnitSystem) = current(one(U),U,ESU)
@pure statvolt(U::UnitSystem) = electricpotential(one(U),U,ESU)
@pure stathenry(U::UnitSystem) = inductance(one(U),U,ESU)
@pure statohm(U::UnitSystem) = resistance(one(U),U,ESU)
@pure statmho(U::UnitSystem) = conductance(one(U),U,ESU)
@pure statfarad(U::UnitSystem) = capacitance(one(U),U,ESU)
@pure statweber(U::UnitSystem) = magneticflux(one(U),U,ESU)
@pure stattesla(U::UnitSystem) = magneticfluxdensity(one(U),U,ESU)
@pure earthcoulomb(U::UnitSystem) = charge(one(U),U,Meridian)

# temperature

@pure sealevel(U::UnitSystem) = temperature(T₀+𝟑*𝟓,U,Metric)
@pure kelvin(U::UnitSystem) = temperature(one(U),U,Metric)
@pure rankine(U::UnitSystem) = temperature(one(U),U,English)
#@pure delisle(U::UnitSystem) = temperature(two(U)/three(U),U,Metric)
#@pure reaumur(U::UnitSystem) = temperature(five(U)/two(U)^2,U,Metric)

# mole

@pure mole(U::UnitSystem) = molaramount(one(U),U,Metric)
@pure earthmole(U::UnitSystem) = molaramount(one(U),U,Meridian)
@pure poundmole(U::UnitSystem) = molaramount(one(U),U,English)
@pure slugmole(U::UnitSystem) = molaramount(one(U),U,British)
@pure slinchmole(U::UnitSystem) = molaramount(one(U),U,IPS)

# photometric

@pure lumen(U::UnitSystem) = luminousflux(one(U),U,Metric)
@pure candela(U::UnitSystem) = luminousintensity(one(U),U,Metric)
@pure lux(U::UnitSystem) = illuminance(one(U),U,Metric)
@pure footcandle(U::UnitSystem) = illuminance(one(U),U,English)
@pure phot(U::UnitSystem) = illuminance(one(U),U,Gauss)
@pure nit(U::UnitSystem) = luminance(one(U),U,Metric)
@pure apostilb(U::UnitSystem) = luminance(two(U)/turn(U),U,Metric)
@pure stilb(U::UnitSystem) = luminance(one(U),U,Gauss)
@pure lambert(U::UnitSystem) = luminance(two(U)/turn(U),U,Gauss)
@pure footlambert(U::UnitSystem) = luminance(two(U)/turn(U),U,English)
@pure bril(U::UnitSystem) = centi(U)*nano(U)*lambert(U)

#const neper = Metric(one(U),log(𝟙))
#const bel = Metric(one(U),log10(𝟙))
#const decibel = Metric(one(U),dB(𝟙))
@pure hertz(U::UnitSystem) = one(U)/second(U)
@pure kayser(U::UnitSystem) = wavenumber(one(U),U,Gauss)
@pure diopter(U::UnitSystem) = wavenumber(one(U),U,Metric)
@pure bubnoff(U::UnitSystem) = meter(U)/year(U)
@pure gforce(U::UnitSystem) = specificforce(one(U),U,English)
@pure galileo(U::UnitSystem) = specificforce(one(U),U,Gauss)
@pure eotvos(U::UnitSystem) = specificforce(nano(U),U,Gauss)/length(one(U),U,Gauss)
@pure darcy(U::UnitSystem) = area(milli(U)/normal(atmosphere(Metric)),U,Gauss)
@pure poise(U::UnitSystem) = viscosity(one(U),U,Gauss)
@pure reyn(U::UnitSystem) = viscosity(one(U),U,IPS)
@pure stokes(U::UnitSystem) = diffusivity(one(U),U,Gauss)
@pure rayl(U::UnitSystem) = specificimpedance(one(U),U,Metric)
@pure katal(U::UnitSystem) = catalysis(one(U),U,Metric)
@pure mpge(U::UnitSystem) = mile(U)/gasgallon(U)
@pure langley(U::UnitSystem) = calorie(U)/(centi(U)*meter(U))^2
@pure jansky(U::UnitSystem) = fluence((Constant(1.0)*deci(U))^26,U,Metric)
@pure solarflux(U::UnitSystem) = hecto(U)^2*jansky(U)
@pure curie(U::UnitSystem) = Constant(37)*giga(U)*hertz(U)
@pure sievert(U::UnitSystem) = energy(one(U),U,Metric)/mass(U,Metric)
@pure rem(U::UnitSystem) = centi(U)*sievert(U)
@pure roentgen(U::UnitSystem) = chargedensity(one(U),U,ESU)/density(Constant(1.293),U,Metric)

