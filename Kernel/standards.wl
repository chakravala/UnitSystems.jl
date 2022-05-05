
(*   This file is part of UnitSystems                   *)
(*   It is licensed under the MIT license               *)
(*   UnitSystems Copyright (C) 2021 Michael Reed        *)
(*       _           _                         _        *)
(*      | |         | |                       | |       *)
(*   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _  *)
(*  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` | *)
(* | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| | *)
(*  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_| *)
(*                                                      *)
(*   https://github.com/chakravala                      *)
(*   https://crucialflow.com                            *)

Pre = <|
"deka" -> 10,
"hecto" -> 10^2,
"kilo" -> 10^3,
"mega" -> 10^6,
"giga" -> 10^9,
"tera" -> 10^12,
"peta" -> 10^15,
"exa" -> 10^18,
"zetta" -> 10^21,
"yotta" -> 10^24,
"deci" -> 10^-1,
"centi" -> 10^-2,
"milli" -> 10^-3,
"micro" -> 10^-6,
"nano" -> 10^-9,
"pico" -> 10^-12,
"femto" -> 10^-15,
"atto" -> 10^-18,
"zepto" -> 10^-21,
"yocto" -> 10^-24,
"byte" -> 2^3,
"kibi" -> 2^10,
"mebi" -> 2^20,
"gibi" -> 2^30,
"tebi" -> 2^40,
"pebi" -> 2^50,
"exbi" -> 2^60,
"zebi" -> 2^70,
"yobi" -> 2^80|>

(*hyperfine,cosmological,lunarmass,gaussianyear,sidereralyear,atmosphere,loschmidt,amagat,wienwavelength,wienfrequency,mechanicalheat,lunardistance,year,lightyear,inchmercury,torr,earthcalorie,thermalunit,gasgallon,tonsrefrigeration,boilerhorsepower,mpge,celsius,sealevel,boiling,darcy,bubnoff,roentgen,electronvolt*)

StandardUnits = <|
"HubbleParameter" -> Measure[1,d1/dT,"Hubble"],
"SolarMass" -> Measure[1,dM,"IAU"],
"EarthMass" -> Measure[1,dM,"IAUE"],
"JupiterMass" -> Measure[1,dM,"IAUJ"],
"GForce" -> Measure[1,USQ[SpecificForce],"English"],
"Radian" -> Measure[1,dA,"MetricEngineering"],
"Steradian" -> Measure[1,dA^2,"MetricEngineering"],
"Degree" -> Measure[Pi/180,dA,"MetricEngineering"],
"Gradian" -> Measure[Pi/200,dA,"MetricEngineering"],
"ArcMinute" -> Measure[Pi/180/60,dA,"MetricEngineering"],
"ArcSecond" -> Measure[Pi/180/60^2,dA,"MetricEngineering"],
"Meter" -> Measure[1,dL,"Metric"],
"Angstrom" -> Measure[Pre["hecto"] Pre["pico"],dL,"Metric"],
"Inch" -> Measure[1,dL,"IPS"],
"Foot" -> Measure[1,dL,"English"],
"SurveyFoot" -> Measure[1,dL,"Survey"],
"Yard" -> Measure[3,dL,"English"],
"Mile" -> Measure[2^5 3 5 11, dL, "English"],
"StatuteMile" -> Measure[2^5 3 5 11, dL, "Survey"],
"EarthRadius" -> Sqrt[Measure[1, dM, "IAUE"]["Metric"] GravitationalConstant["Metric"]/Measure[1, USQ[SpecificForce], "English"]["Metric"]],
"GreatCircle" -> 2 Pi Sqrt[Measure[1, dM, "IAUE"]["Metric"] GravitationalConstant["Metric"]/Measure[1, USQ[SpecificForce], "English"]["Metric"]],
"EarthMeter" -> Measure[1,dL,"Meridian"],
"NauticalMile" -> Measure[1,dL,"Nautical"],
"AdmiraltyMile" -> Measure[2^6 5 19,dL,"English"],
"MeridianMile" -> Measure[2^4 5^5/3^3,dL,"Metric"],
"AstronomicalUnit" -> Measure[1,dL,"IAU"],
"ParSec" -> Measure[2^6 3^4 5^3/Pi,dL,"IAU"],
"Second" -> Measure[1,dT,"Metric"],
"Minute" -> Measure[60,dT,"Metric"],
"Hour" -> Measure[60^2,dT,"Metric"],
"Day" -> Measure[1,dT,"IAU"],
"RadarMile" -> Measure[2,dL,"Nautical"]/SpeedOfLight["Nautical"],
"ips" -> Measure[1,dL/dT,"IPS"],
"fps" -> Measure[1,dL/dT,"British"],
"fpm" -> Measure[1,dL/dT,"British"],
"ms" -> Measure[1,dL/dT,"Metric"],
"kmh" -> Measure[1,dL/dT,"KKH"],
"mph" -> Measure[1,dL/dT,"MPH"],
"knot" -> Measure[1,dL/dT,"Nautical"],
"Barn" -> Measure[10^-28,dL^2,"Metric"],
"Hectare" -> Measure[Pre["hecto"] Pre["hecto"],dL^2,"Metric"],
"Acre" -> Measure[2^-7/5,dL^2,"MPH"],
"SurveyAcre" -> Measure[2^3 3^2 5 11^2,dL^2,"Survey"],
"Gallon" -> Measure[3 7 11,dL^3,"IPS"],
"Liter" -> Measure[10^-3,dL^3,"Metric"],
"Quart" -> Measure[3 7 11/4,dL^3,"IPS"],
"Pint" -> Measure[3 7 11/8,dL^3,"IPS"],
"Cup" -> Measure[3 7 11/16,dL^3,"IPS"],
"FluidOunce" -> Measure[3 7 11/2^7,dL^3,"IPS"],
"TeaSpoon" -> Measure[5/10^3,dL^3,"Metric"],
"TableSpoon" -> Measure[15/10^3,dL^3,"Metric"],
"Gram" -> Measure[10^-3,dM,"Metric"],
"EarthGram" -> Measure[10^-3,dM,"Meridian"],
"Kilogram" -> Measure[1,dM,"Metric"],
"Tonne" -> Measure[10^3,dM,"Metric"],
"Ton" -> Measure[2 10^3,dM,"English"],
"Pound" -> Measure[1,dM,"English"],
"Ounce" -> Measure[2^-4,dM,"English"],
"Grain" -> Measure[10^-3/7,dM,"English"],
"Slug" -> Measure[1,dM,"British"],
"Slinch" -> Measure[1,dM,"IPS"],
"Hyl" -> Measure[1,dM,"GravitationalMetric"],
"Dyne" -> Measure[1,dF,"Gauss"],
"Newton" -> Measure[1,dF,"Metric"],
"Poundal" -> Measure[1,dF,"FPS"],
"Kilopond" -> Measure[1,dF,"MetricEngineering"],
"PoundForce" -> Measure[1,dF,"English"],
"Psi" -> Measure[1,USQ[Pressure],"IPS"],
"Bar" -> Measure[10^5,USQ[Pressure],"Metric"],
"Barye" -> Measure[1,USQ[Pressure],"Gauss"],
"Pascal" -> Measure[1,USQ[Pressure],"Metric"],
"TechnicalAtmosphere" -> Measure[10^-4,USQ[Pressure],"MetricEngineering"],
"Erg" -> Measure[1,USQ[Energy],"Gauss"],
"Joule" -> Measure[1,USQ[Energy],"Metric"],
"FootPound" -> Measure[1,USQ[Energy],"English"],
"MeanCalorie" -> Measure[2^2 5 3^2/43,USQ[Energy],"InternationalMean"],
"KiloCalorie" -> Measure[2^5 5^4 3^2/43,USQ[Energy],"International"],
"Calorie" -> Measure[2^2 5 3^2/43,USQ[Energy],"International"],
"TonTNT" -> Measure[Pre["giga"] 2^2 5 3^2/43,USQ[Energy],"International"],
"Watt" -> Measure[1,USQ[Power],"Metric"],
"HorsePower" -> Measure[2 5^2 11,USQ[Power],"British"],
"HorsePowerWatt" -> Measure[Pi 2^5 3^3/5,USQ[Power],"British"],
"HorsePowerMetric" -> Measure[75,USQ[Power],"GravitationalMetric"],
"ElectricalHorsePower" -> Measure[746,USQ[Power],"Metric"],
"Coulomb" -> Measure[1,dQ,"Metric"],
"Ampere" -> Measure[1,USQ[ElectricCurrent],"Metric"],
"Volt" -> Measure[1,USQ[ElectricPotential],"Metric"],
"Henry" -> Measure[1,USQ[MagneticInductance],"Metric"],
"Ohm" -> Measure[1,USQ[ElectricResistance],"Metric"],
"Siemens" -> Measure[1,USQ[ElectricConductance],"Metric"],
"Farad" -> Measure[1,USQ[ElectricCapacitance],"Metric"],
"Weber" -> Measure[1,USQ[MagneticFlux],"Metric"],
"Tesla" -> Measure[1,USQ[MagneticFluxDensity],"Metric"],
"StatCoulomb" -> Measure[1,dQ,"ESU"],
"StatAmpere" -> Measure[1,USQ[ElectricCurrent],"ESU"],
"StatVolt" -> Measure[1,USQ[ElectricPotential],"ESU"],
"StatHenry" -> Measure[1,USQ[MagneticInductance],"ESU"],
"StatOhm" -> Measure[1,USQ[ElectricResistance],"ESU"],
"StatSiemens" -> Measure[1,USQ[ElectricConductance],"ESU"],
"StatFarad" -> Measure[1,USQ[ElectricCapacitance],"ESU"],
"StatWeber" -> Measure[1,USQ[MagneticFlux],"ESU"],
"StatTesla" -> Measure[1,USQ[MagneticFluxDensity],"ESU"],
"AbCoulomb" -> Measure[1,dQ,"EMU"],
"AbAmpere" -> Measure[1,USQ[ElectricCurrent],"EMU"],
"AbVolt" -> Measure[1,USQ[ElectricPotential],"EMU"],
"AbHenry" -> Measure[1,USQ[MagneticInductance],"EMU"],
"AbOhm" -> Measure[1,USQ[ElectricResistance],"EMU"],
"AbSiemens" -> Measure[1,USQ[ElectricConductance],"EMU"],
"AbFarad" -> Measure[1,USQ[ElectricCapacitance],"EMU"],
"Maxwell" -> Measure[1,USQ[MagneticFlux],"EMU"],
"Gauss" -> Measure[1,USQ[MagneticFluxDensity],"EMU"],
"Oersted" -> Measure[1,USQ[MagneticFieldStrength],"EMU"],
"Gilberg" -> Measure[1/4/Pi,dQ/dA,"EMU"],
"EarthCoulomb" -> Measure[1,dQ,"Meridian"],
"Kelvin" -> Measure[1,d0,"Metric"],
"Rankine" -> Measure[1,d0,"English"],
"Fahrenheit" -> Measure[45967/100,d0,"English"],
"Mole" -> Measure[1,dN,"Metric"],
"EarthMole" -> Measure[1,dN,"Meridian"],
"PoundMole" -> Measure[1,dN,"English"],
"SlugMole" -> Measure[1,dN,"British"],
"SlinchMole" -> Measure[1,dN,"IPS"],
"Lumen" -> Measure[1,USQ[LuminousFlux],"Metric"],
"Candela" -> Measure[1,USQ[LuminousIntensity],"Metric"],
"Lux" -> Measure[1,USQ[Illuminance],"Metric"],
"Phot" -> Measure[1,USQ[Illuminance],"Gauss"],
"FootCandle" -> Measure[1,USQ[Illuminance],"English"],
"Nit" -> Measure[1,USQ[Luminance],"Metric"],
"Apostilb" -> Measure[1/Pi,USQ[Luminance],"Metric"],
"Stilb" -> Measure[1,USQ[Luminance],"Gauss"],
"Lambert" -> Measure[1/Pi,USQ[Luminance],"Gauss"],
"FootLambert" -> Measure[1/Pi,USQ[Luminance],"English"],
"Bril" -> Measure[Pre["centi"] Pre["nano"]/Pi,USQ[Luminance],"Gauss"],
"Hertz" -> Measure[1,d1/dT,"Metric"],
"APM" -> Measure[1/60,d1/dT,"Metric"],
"RPM" -> Measure[Pi/30,d1/dT,"Metric"],
"Galileo" -> Measure[1,USQ[SpecificForce],"Gauss"],
"Eotvos" -> Measure[Pre["nano"],USQ[SpecificForce]/dL,"Gauss"],
"Poise" -> Measure[1,USQ[DynamicViscosity],"Gauss"],
"Reyn" -> Measure[1,USQ[DynamicViscosity],"IPS"],
"Diopter" -> Measure[1,d1/dL,"Metric"],
"Kayser" -> Measure[1,d1/dL,"Gauss"],
"Stokes" -> Measure[1,USQ[KinematicViscosity],"Gauss"],
"Katal" -> Measure[1,USQ[Catalysis],"Metric"],
"Curie" -> Measure[37 Pre["giga"],d1/dT,"Metric"],
"Sievert" -> Measure[1,USQ[Energy]/dM,"Metric"],
"REM" -> Measure[1/100,USQ[Energy]/dM,"Metric"],
"Rayl" -> Measure[1,USQ[SpecificAcousticImpedance],"Metric"],
"Langley" -> Measure[2^2 5 3^2/43,USQ[Energy],"International"]/Measure[10^-4,dL^2,"Metric"],
"Jansky" -> Measure[10^-26,USQ[EnergyPerArea],"Metric"],
"SolarFlux" -> Measure[10^-22,USQ[EnergyPerArea],"Metric"]|>



(*Map[(
transformBoxes[transform[#,USQ[Frequency]],#] := "Hz";
transformBoxes[transform[#,USQ[FrequencyDrift]],#] := "Hz"/"s";
transformBoxes[transform[#,USQ[Illuminance]],#] := "lx";
transformBoxes[transform[#,USQ[LuminousExposure]],#] := "lx" "s";
transformBoxes[#,USQ[Luminance]] := "nt";
) &, {"MetricEngineering","SI2019Engineering","GravitationalMetric","GravitationalSI2019"}]

Map[(
transformBoxes[transform[#,USQ[SpecificForce]],#] := "g0";
) &, {"MetricEngineering","SI2019Engineering","MeridianEngineering","English","Survey"}]

Map[(
transformBoxes[transform[#,USQ[Frequency]],#] := "Hz";
transformBoxes[transform[#,USQ[FrequencyDrift]],#] := "Hz"/"s";
transformBoxes[transform[#,USQ[Force]],#] := "N";
transformBoxes[transform[#,d1/USQ[Force]],#] := 1/"N";
transformBoxes[transform[#,USQ[Pressure]],#] := "Pa";
transformBoxes[transform[#,USQ[Compressibility]],#] := 1/"Pa";
transformBoxes[transform[#,USQ[Energy]],#] := "J";
transformBoxes[transform[#,d1/USQ[Energy]],#] := 1/"J";
transformBoxes[transform[#,USQ[Power]],#] := "W";
transformBoxes[transform[#,d1/USQ[Power]],#] := 1/"W";

transformBoxes[transform[#,USQ[ElectricPotential]],#] := "V";
transformBoxes[transform[#,d1/USQ[ElectricPotential]],#] := 1/"V";
transformBoxes[transform[#,USQ[ElectricCapacitance]],#] := "F";
transformBoxes[transform[#,d1/USQ[ElectricCapacitance]],#] := 1/"F";
transformBoxes[transform[#,USQ[ElectricResistance]],#] := "\[Omega]";
transformBoxes[transform[#,USQ[ElectricConductance]],#] := "S";
transformBoxes[transform[#,USQ[MagneticFlux]],#] := "Wb";
transformBoxes[transform[#,d1/USQ[MagneticFlux]],#] := "Hz"/"V";
transformBoxes[transform[#,USQ[MagneticFluxDensity]],#] := "T";
transformBoxes[transform[#,d1/USQ[MagneticFluxDensity]],#] := 1/"T";
transformBoxes[transform[#,USQ[MagneticPermeance]],#] := "H";
transformBoxes[transform[#,USQ[MagneticReluctance]],#] := 1/"H";

transformBoxes[transform[#,USQ[Catalysis]],#] := "kat";
transformBoxes[transform[#,USQ[MolarEnergy]],#] := "J"/"mol";
transformBoxes[transform[#,USQ[MolarEntropy]],#] := "J"/"K"/"mol";

transformBoxes[transform[#,USQ[LuminousEfficacyOfRadiation]],#] := "lm"/"W";
transformBoxes[transform[#,d1/USQ[LuminousEfficacyOfRadiation]],#] := "W"/"lm";
transformBoxes[transform[#,USQ[LuminousIntensity]],#] := "cd";
transformBoxes[transform[#,USQ[Illuminance]],#] := "lx";
transformBoxes[transform[#,USQ[LuminousExposure]],#] := "lx" "s";
transformBoxes[#,USQ[Luminance]] := "nt";

transformBoxes[transform[#,USQ[AngularMomentum]],#] := "J" "s";
transformBoxes[transform[#,USQ[Action] USQ[Speed]],#] := "J" "m";
transformBoxes[transform[#,USQ[Impulse]],#] := "N" "s";
transformBoxes[transform[#,USQ[ForceOnsetRate]],#] := "N"/"s";
transformBoxes[transform[#,USQ[EnergyPerArea]],#] := "N"/"m";
transformBoxes[transform[#,USQ[Compliance]],#] := "m"/"N";

transformBoxes[transform[#,USQ[DynamicViscosity]],#] := "Pa" "s";
transformBoxes[transform[#,USQ[Irradiance]],#] := "W"/"m"^2;
transformBoxes[transform[#,USQ[PowerDensity]],#] := "W"/"m"^3;
transformBoxes[transform[#,USQ[Irradiance]/d0^4],#] := "W"/"m"^2/"K"^4;
transformBoxes[transform[#,USQ[Pressure]/d0^4],#] := "J"/"m"^3/"K"^4;
transformBoxes[transform[#,d1/dT/d0],#] := "Hz"/"K";
transformBoxes[transform[#,USQ[Entropy]/dQ],#] := "V"/"K";
transformBoxes[transform[#,USQ[Entropy]],#] := "J"/"K";
transformBoxes[transform[#,USQ[SpecificEntropy]],#] := "J"/"K"/"kg";
transformBoxes[transform[#,USQ[SpecificEnergy]],#] := "J"/"kg";
transformBoxes[transform[#,USQ[ThermalConductivity]],#] := "W"/"m"/"K";
transformBoxes[transform[#,USQ[ThermalConductance]],#] := "W"/"K";
transformBoxes[transform[#,USQ[ThermalResistance]],#] := "K"/"W";
transformBoxes[transform[#,USQ[ThermalResistivity]],#] := "K" "m"/"W";
transformBoxes[transform[#,USQ[MolarConductivity]],#] := "S" "m"^2/"mol";

transformBoxes[transform[#,USQ[ElectricPotential]/dM],#] := "V"/"kg";
transformBoxes[transform[#,USQ[ElectricPotential] dL],#] := "V" "m";
transformBoxes[transform[#,USQ[ElectricFieldStrength]],#] := "V"/"m";
transformBoxes[transform[#,USQ[ElectricPermittivity]],#] := "F"/"m";
transformBoxes[transform[#,d1/USQ[ElectricPermittivity]],#] := "m"/"F";
transformBoxes[transform[#,USQ[MagneticPermeability]],#] := "H"/"m";
transformBoxes[transform[#,d1/USQ[MagneticPermeability]],#] := "m"/"H";
transformBoxes[transform[#,USQ[Resistivity]],#] := "\[Omega]" "m";
transformBoxes[transform[#,USQ[Conductivity]],#] := "S"/"m";
transformBoxes[transform[#,USQ[MagneticDipoleMoment]],#] := "J"/"T";
transformBoxes[transform[#,USQ[MagneticVectorPotential]],#] := "Wb"/"m";
transformBoxes[transform[#,USQ[MagneticMoment]],#] := "Wb" "m";
transformBoxes[transform[#,USQ[ElectricalMobility]],#] := "m"^2/"s"/"V";
) &, {"Metric","SI2019","CODATA","Conventional","International","InternationalMean"}]

transformBoxes[transform["Meridian",USQ[Frequency]],"Meridian"] := "Hz";
transformBoxes[transform["Meridian",USQ[FrequencyDrift]],"Meridian"] := "Hz"/"s";
transformBoxes[transform["Meridian",USQ[Force]],"Meridian"] := "eN";
transformBoxes[transform["Meridian",d1/USQ[Force]],"Meridian"] := 1/"eN";
transformBoxes[transform["Meridian",USQ[Pressure]],"Meridian"] := "ePa";
transformBoxes[transform["Meridian",USQ[Compressibility]],"Meridian"] := 1/"ePa";
transformBoxes[transform["Meridian",USQ[Energy]],"Meridian"] := "eJ";
transformBoxes[transform["Meridian",d1/USQ[Energy]],"Meridian"] := 1/"eJ";
transformBoxes[transform["Meridian",USQ[Power]],"Meridian"] := "eW";
transformBoxes[transform["Meridian",d1/USQ[Power]],"Meridian"] := 1/"eW";

transformBoxes[transform["Meridian",USQ[ElectricPotential]],"Meridian"] := "eV";
transformBoxes[transform["Meridian",d1/USQ[ElectricPotential]],"Meridian"] := 1/"eV";
transformBoxes[transform["Meridian",USQ[ElectricCapacitance]],"Meridian"] := "eF";
transformBoxes[transform["Meridian",d1/USQ[ElectricCapacitance]],"Meridian"] := 1/"eF";
transformBoxes[transform["Meridian",USQ[ElectricResistance]],"Meridian"] := "e\[Omega]";
transformBoxes[transform["Meridian",USQ[ElectricConductance]],"Meridian"] := "eS";
transformBoxes[transform["Meridian",USQ[MagneticFlux]],"Meridian"] := "eWb";
transformBoxes[transform["Meridian",d1/USQ[MagneticFlux]],"Meridian"] := "Hz"/"eV";
transformBoxes[transform["Meridian",USQ[MagneticFluxDensity]],"Meridian"] := "eT";
transformBoxes[transform["Meridian",d1/USQ[MagneticFluxDensity]],"Meridian"] := 1/"eT";
transformBoxes[transform["Meridian",USQ[MagneticPermeance]],"Meridian"] := "eH";
transformBoxes[transform["Meridian",USQ[MagneticReluctance]],"Meridian"] := 1/"eH";

transformBoxes[transform["Meridian",USQ[Catalysis]],"Meridian"] := "ekat";
transformBoxes[transform["Meridian",USQ[MolarEnergy]],"Meridian"] := "eJ"/"eg-mol";
transformBoxes[transform["Meridian",USQ[MolarEntropy]],"Meridian"] := "eJ"/"K"/"eg-mol";

transformBoxes[transform["Meridian",USQ[LuminousEfficacyOfRadiation]],"Meridian"] := "lm"/"eW";
transformBoxes[transform["Meridian",d1/USQ[LuminousEfficacyOfRadiation]],"Meridian"] := "eW"/"lm";
transformBoxes[transform["Meridian",USQ[LuminousIntensity]],"Meridian"] := "cd";
transformBoxes[transform["Meridian",USQ[Illuminance]],"Meridian"] := "lx";
transformBoxes[transform["Meridian",USQ[LuminousExposure]],"Meridian"] := "lx" "s";
transformBoxes[USQ[Luminance],"Meridian"] := "nt";

transformBoxes[transform["Meridian",USQ[AngularMomentum]],"Meridian"] := "eJ" "s";
transformBoxes[transform["Meridian",USQ[Action] USQ[Speed]],"Meridian"] := "eJ" "em";
transformBoxes[transform["Meridian",USQ[Impulse]],"Meridian"] := "eN" "s";
transformBoxes[transform["Meridian",USQ[ForceOnsetRate]],"Meridian"] := "eN"/"s";
transformBoxes[transform["Meridian",USQ[EnergyPerArea]],"Meridian"] := "eN"/"em";
transformBoxes[transform["Meridian",USQ[Compliance]],"Meridian"] := "em"/"eN";

transformBoxes[transform["Meridian",USQ[DynamicViscosity]],"Meridian"] := "ePa" "s";
transformBoxes[transform["Meridian",USQ[Irradiance]],"Meridian"] := "eW"/"em"^2;
transformBoxes[transform["Meridian",USQ[PowerDensity]],"Meridian"] := "eW"/"em"^3;
transformBoxes[transform["Meridian",USQ[Irradiance]/d0^4],"Meridian"] := "eW"/"em"^2/"K"^4;
transformBoxes[transform["Meridian",USQ[Pressure]/d0^4],"Meridian"] := "eJ"/"em"^3/"K"^4;
transformBoxes[transform["Meridian",d1/dT/d0],"Meridian"] := "Hz"/"K";
transformBoxes[transform["Meridian",USQ[Entropy]/dQ],"Meridian"] := "eV"/"K";
transformBoxes[transform["Meridian",USQ[Entropy]],"Meridian"] := "eJ"/"K";
transformBoxes[transform["Meridian",USQ[SpecificEntropy]],"Meridian"] := "eJ"/"K"/"keg";
transformBoxes[transform["Meridian",USQ[SpecificEnergy]],"Meridian"] := "eJ"/"keg";
transformBoxes[transform["Meridian",USQ[ThermalConductivity]],"Meridian"] := "eW"/"em"/"K";
transformBoxes[transform["Meridian",USQ[ThermalConductance]],"Meridian"] := "eW"/"K";
transformBoxes[transform["Meridian",USQ[ThermalResistance]],"Meridian"] := "K"/"eW";
transformBoxes[transform["Meridian",USQ[ThermalResistivity]],"Meridian"] := "K" "em"/"eW";
transformBoxes[transform["Meridian",USQ[MolarConductivity]],"Meridian"] := "eS" "em"^2/"eg-mol";

transformBoxes[transform["Meridian",USQ[ElectricPotential]/dM],"Meridian"] := "eV"/"keg";
transformBoxes[transform["Meridian",USQ[ElectricPotential] dL],"Meridian"] := "eV" "em";
transformBoxes[transform["Meridian",USQ[ElectricFieldStrength]],"Meridian"] := "eV"/"em";
transformBoxes[transform["Meridian",USQ[ElectricPermittivity]],"Meridian"] := "eF"/"em";
transformBoxes[transform["Meridian",d1/USQ[ElectricPermittivity]],"Meridian"] := "em"/"eF";
transformBoxes[transform["Meridian",USQ[MagneticPermeability]],"Meridian"] := "eH"/"em";
transformBoxes[transform["Meridian",d1/USQ[MagneticPermeability]],"Meridian"] := "em"/"eH";
transformBoxes[transform["Meridian",USQ[Resistivity]],"Meridian"] := "e\[Omega]" "em";
transformBoxes[transform["Meridian",USQ[Conductivity]],"Meridian"] := "eS"/"em";
transformBoxes[transform["Meridian",USQ[MagneticDipoleMoment]],"Meridian"] := "eJ"/"eT";
transformBoxes[transform["Meridian",USQ[MagneticVectorPotential]],"Meridian"] := "eWb"/"em";
transformBoxes[transform["Meridian",USQ[MagneticMoment]],"Meridian"] := "eWb" "em";
transformBoxes[transform["Meridian",USQ[ElectricalMobility]],"Meridian"] := "em"^2/"s"/"eV";

Map[(
transformBoxes[transform[#,USQ[Frequency]],#] := "Hz";
transformBoxes[transform[#,USQ[Force]],#] := "dyn";
transformBoxes[transform[#,d1/USQ[Force]],#] := 1/"dyn";
transformBoxes[transform[#,USQ[SpecificForce]],#] := "gal";
transformBoxes[transform[#,USQ[SpecificForce]/dL],#] := "gal"/"cm";
transformBoxes[transform[#,USQ[Pressure]],#] := "Ba";
transformBoxes[transform[#,USQ[Compressibility]],#] := 1/"Ba";
transformBoxes[transform[#,USQ[Energy]],#] := "erg";
transformBoxes[transform[#,d1/USQ[Energy]],#] := 1/"erg";
transformBoxes[transform[#,USQ[Power]],#] := "erg"/"s";
transformBoxes[transform[#,d1/USQ[Power]],#] := "s"/"erg";

transformBoxes[transform[#,USQ[Catalysis]],#] := "kat";
transformBoxes[transform[#,USQ[MolarEnergy]],#] := "erg"/"mol";
transformBoxes[transform[#,USQ[MolarEntropy]],#] := "erg"/"K"/"mol";

transformBoxes[transform[#,USQ[LuminousEfficacyOfRadiation]],#] := "lm" "s"/"erg";
transformBoxes[transform[#,d1/USQ[LuminousEfficacyOfRadiation]],#] := "erg"/"s"/"lm";
transformBoxes[transform[#,USQ[LuminousIntensity]],#] := "cd";
transformBoxes[transform[#,USQ[Illuminance]],#] := "ph";
transformBoxes[#,USQ[Luminance]] := "sb";

transformBoxes[transform[#,USQ[AngularMomentum]],#] := "erg" "s";
transformBoxes[transform[#,USQ[Action] USQ[Speed]],#] := "erg" "cm";
transformBoxes[transform[#,USQ[Impulse]],#] := "dyn" "s";
transformBoxes[transform[#,USQ[ForceOnsetRate]],#] := "dyn"/"s";
transformBoxes[transform[#,USQ[EnergyPerArea]],#] := "dyn"/"cm";
transformBoxes[transform[#,USQ[Compliance]],#] := "cm"/"dyn";

transformBoxes[transform[#,USQ[DynamicViscosity]],#] := "P";
transformBoxes[transform[#,USQ[KinematicViscosity]],#] := "St";
transformBoxes[transform[#,USQ[Irradiance]],#] := "erg"/"s"/"cm"^2;
transformBoxes[transform[#,USQ[PowerDensity]],#] := "erg"/"s"/"cm"^3;
transformBoxes[transform[#,USQ[Irradiance]/d0^4],#] := "erg"/"s"/"cm"^2/"K"^4;
transformBoxes[transform[#,USQ[Pressure]/d0^4],#] := "Ba"/"K"^4;
transformBoxes[transform[#,d1/dT/d0],#] := "Hz"/"K";
transformBoxes[transform[#,USQ[Entropy]],#] := "erg"/"K";
transformBoxes[transform[#,USQ[SpecificEntropy]],#] := "erg"/"K"/"g";
transformBoxes[transform[#,USQ[SpecificEnergy]],#] := "erg"/"g";
transformBoxes[transform[#,USQ[ThermalConductivity]],#] := "erg"/"s"/"cm"/"K";
transformBoxes[transform[#,USQ[ThermalConductance]],#] := "erg"/"s"/"K";
transformBoxes[transform[#,USQ[ThermalResistance]],#] := "K" "s"/"erg";
transformBoxes[transform[#,USQ[ThermalResistivity]],#] := "K" "cm" "s"/"erg";
) &, {"Gauss", "EMU", "ESU", "LorentzHeaviside"}]

transformBoxes[transform["EMU",USQ[ElectricCurrent]],"EMU"] := "Bi";
transformBoxes[transform["EMU",USQ[MagneticFlux]],"EMU"] := "Mx";
transformBoxes[transform["EMU",USQ[MagneticFluxDensity]],"EMU"] := "G";
transformBoxes[transform["EMU",USQ[MagneticFieldStrength]],"EMU"] := "Oe";
transformBoxes[transform["EMU",USQ[MagneticReluctance]],"EMU"] := "Bi"/"Mx";
transformBoxes[transform["EMU",USQ[MagneticDipoleMoment]],"EMU"] := "erg"/"G";
transformBoxes[transform["EMU",USQ[MagneticVectorPotential]],"EMU"] := "Mx"/"cm";
transformBoxes[transform["EMU",USQ[MagneticMoment]],"EMU"] := "Mx" "cm";
transformBoxes[transform["EMU",USQ[MagneticPoleStrength]],"EMU"] := "pole";

transformBoxes[transform["Gauss",USQ[ElectricCharge]],"Gauss"] := "Fr";
transformBoxes[transform["Gauss",USQ[MagneticFlux]],"Gauss"] := "Mx";
transformBoxes[transform["Gauss",USQ[MagneticFluxDensity]],"Gauss"] := "G";
transformBoxes[transform["Gauss",USQ[MagneticFieldStrength]],"Gauss"] := "Oe";
transformBoxes[transform["Gauss",USQ[MagneticReluctance]],"Gauss"] := "Fr"/"s"/"Mx";
transformBoxes[transform["Gauss",USQ[MagneticDipoleMoment]],"Gauss"] := "erg"/"G";
transformBoxes[transform["Gauss",USQ[MagneticVectorPotential]],"Gauss"] := "Mx"/"cm";
transformBoxes[transform["Gauss",USQ[MagneticMoment]],"Gauss"] := "Mx" "cm";

transformBoxes[transform["MTS",USQ[Force]],"MTS"] := "sn";
transformBoxes[transform["MTS",d1/USQ[Force]],"MTS"] := 1/"sn";
transformBoxes[transform["MTS",USQ[Pressure]],"MTS"] := "p";
transformBoxes[transform["MTS",USQ[Compressibility]],"MTS"] := 1/"pz";

transformBoxes[transform["GravitationalMetric",USQ[Mass]],"GravitationalMetric"] := "hyl";
transformBoxes[transform["GravitationalSI2019",USQ[Mass]],"GravitationalSI2019"] := "hyl";
transformBoxes[transform["GravitationalMeridian",USQ[Mass]],"GravitationalMeridian"] := "ehyl";
transformBoxes[transform["British",USQ[Mass]],"British"] := "slug";
transformBoxes[transform["IPS",USQ[Mass]],"IPS"] := "slinch";
transformBoxes[transform["FPS",USQ[Force]],"FPS"] := "pdl";
transformBoxes[transform["FPS",USQ[Pressure]],"FPS"] := "pdl"/"ft"^2;
(*transformBoxes[transform["British",USQ[Density]],"British"] := "slug"/"ft"^3;
transformBoxes[transform["IPS",USQ[Density]],"IPS"] := "slinch"/"in"^3;
transformBoxes[transform["GravitationalMetric",USQ[Density]],"GravitationalMetric"] := "hyl"/"m"^3;
transformBoxes[transform["GravitationalSI2019",USQ[Density]],"GravitationalSI2019"] := "hyl"/"m"^3;
transformBoxes[transform["GravitationalMeridian",USQ[Density]],"GravitationalMeridian"] := "ehyl"/"m"^3;*)

Map[(
transformBoxes[transform[#,USQ[Frequency]],#] := "Hz";
transformBoxes[transform[#,USQ[FrequencyDrift]],#] := "Hz"/"s";
transformBoxes[transform[#,d1/dT/d0],#] := "Hz"/"\[Degree]R";
) &, {"FPS", "IPS", "British", "English", "Survey"}]

Map[(
transformBoxes[transform[#,USQ[LuminousIntensity]],#] := "cd";
transformBoxes[transform[#,USQ[Illuminance]],#] := "fc";
) &, {"FPS", "British", "English", "Survey"}]
*)
