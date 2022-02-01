
(* This file is part of UnitSystems. It is licensed under the MIT license *)
(* UnitSystems Copyright (C) 2021 Michael Reed *)

UnitData = <|
"g0" -> 980665/10^5,
"ft" -> 3048/10000,
"ftUS" -> 1200/3937,
"lb" -> 45359237/10^8,
"rankine" -> 5/9,
"\[CapitalDelta]\[Nu]Cs" -> 9192631770,
"Kcd" -> 683 555016/555000,
"mP" -> Around[2.176434/10^8, 24/10^14],
"NA" -> 602214076 10^15,
"kB" -> 1380649/10^29,
"h" -> 662607015/10^42,
"c" -> 299792458,
"e" -> 1602176634/10^28,
"\[Mu]eu" -> 1/Around[1822.888486209, 53/10^9],
"\[Mu]pu" -> Around[1.007276466621, 53/10^12],
"\[Alpha]" -> 1/Around[137.035999084, 21/10^9],
"R\[Infinity]" -> Around[10973731.5681601, 21/10^6],
"RK1990" -> 25812807/1000,
"RK2014" -> Around[25812.8074555, 59/10^7],
"KJ1990" -> 4835979 10^8,
"KJ2014" -> Around[4.835978525 10^14, 3 10^6],
"GMsun" -> Around[1.32712442099 10^20, 9 10^9],
"day" -> 60^2 24,
"au" -> 149597870700,
"fur" -> 201168/1000,
"H0" -> Around[67.66, 0.42],
"\[CapitalOmega]\[CapitalLambda]" -> Around[0.6889, 0.0056]|>

AbstractUnitData = <|
"slug" -> "lb" "g0"/"ft",
"slugUS" -> "lb" "g0"/"ftUS",
"Ru" -> "NA" "kB",
"me" -> 2 "R\[Infinity]" "h"/"\[Alpha]"^2/"c",
"\[Mu]0" -> 2 "\[Alpha]" "h"/"c"/"e"^2,
"\[HBar]" -> "h"/(2 Pi),
"\[Mu]pe" -> "\[Mu]pu"/"\[Mu]eu",
"pc" -> "au" 3*60^3/Pi|>

AppendTo[AbstractUnitData, "\[Alpha]G" -> (AbstractUnitData["me"]/"mP")^2]
AppendTo[AbstractUnitData, "G" -> "c" AbstractUnitData["\[HBar]"]/"mP"^2]
AppendTo[AbstractUnitData, "Mu" -> "NA" AbstractUnitData["me"]/"\[Mu]eu"]
AppendTo[AbstractUnitData, "th" -> 1000 AbstractUnitData["pc"]/"H0"]
AppendTo[AbstractUnitData, "\[CapitalLambda]" -> 3 (1/AbstractUnitData["th"]/"c")^2 "\[CapitalOmega]\[CapitalLambda]"]
AppendTo[AbstractUnitData, "lc" -> PowerExpand[2 Sqrt[2 Pi/AbstractUnitData["\[CapitalLambda]"]]]]
AppendTo[AbstractUnitData, "mc" ->  PowerExpand["c"^2/(2 Sqrt[2 Pi AbstractUnitData["\[CapitalLambda]"]] AbstractUnitData["G"])]]
AppendTo[AbstractUnitData, "\[Rho]\[CapitalLambda]" -> "c"^4 AbstractUnitData["\[CapitalLambda]"]/(8 Pi AbstractUnitData["G"])]
AppendTo[AbstractUnitData, "lcq" -> PowerExpand[("c" AbstractUnitData["\[HBar]"]/AbstractUnitData["\[Rho]\[CapitalLambda]"])^(1/4)]]
AppendTo[AbstractUnitData, "mcq" -> PowerExpand[(AbstractUnitData["\[HBar]"]^3 AbstractUnitData["\[Rho]\[CapitalLambda]"]/"c"^5)^(1/4)]]
AppendTo[AbstractUnitData, "ecq" -> PowerExpand[("c"^3 AbstractUnitData["\[HBar]"]^3 AbstractUnitData["\[Rho]\[CapitalLambda]"])^(1/4)]]
AppendTo[AbstractUnitData, "tcq" -> PowerExpand[AbstractUnitData["lcq"] Sqrt[AbstractUnitData["mcq"]/AbstractUnitData["ecq"]]]]

EvalUnitData[x] := AbstractUnitData[x]/.Normal[UnitData]

Coupling["AbstractCoupling"] = Coupling["\[Alpha]G", "\[Alpha]", "\[Mu]eu", "\[Mu]pu", "\[CapitalOmega]\[CapitalLambda]"]
Coupling["AbstractUniverse"] = Coupling[AbstractUnitData["\[Alpha]G"], "\[Alpha]", "\[Mu]eu", "\[Mu]pu", "\[CapitalOmega]\[CapitalLambda]"]
Coupling["StandardModel"] = Coupling["AbstractUniverse"] /. Normal[UnitData]

MetricSystem[mu_, perm_] := MetricSystem[mu, perm, AbstractUnitData["Ru"]];
MetricSystem[mu_, perm_, ru_] := MetricSystem[mu, perm, ru, AbstractUnitData["\[HBar]"], AbstractUnitData["me"]];
MetricSystem[mu_, perm_, ru_, hbar_, mass_] := UnitSystem[ru mass/mu/"\[Mu]eu", hbar, "c", perm, mass];

MetricSystem["SI2019"] := MetricSystem[AbstractUnitData["Mu"], AbstractUnitData["\[Mu]0"]]
MetricSystem["Metric"] := MetricSystem[1/1000, 4 Pi/10^7]
MetricSystem["SI1976"] := MetricSystem[1/1000, 4 Pi/10^7, 831432/10^5]

ConventionalSystem[klitz_, joseph_] := Module[{hbar = 2/klitz/joseph^2/Pi},
	MetricSystem[1/1000, 2 "\[Alpha]" klitz/"c", AbstractUnitData["Ru"], hbar, 4 Pi "R\[Infinity]" hbar/"\[Alpha]"^2/"c"]]

ConventionalSystem["CODATA"] := ConventionalSystem["RK2014", "KJ2014"]
ConventionalSystem["Conventional"] := ConventionalSystem["RK1990", "KJ1990"]

GaussSystem[u_String] := ($Failed /; False)
GaussSystem[perm_, u__] := GaussSystem["AbstractMetric", perm, u]
GaussSystem[u_String, perm_, ratio_] := GaussSystem[u, perm, ratio, 1]
GaussSystem[u_String, perm_, ratio_, lorentz_] := GaussSystem[u, perm, ratio, lorentz, 1/100, 1/1000]
GaussSystem[u_String, perm_, ratio_, lorentz_, length_, mass_] :=
	GaussSystem[u, perm, ratio, lorentz, length, mass, PowerExpand[mass length^2]]
GaussSystem[u_String, perm_, ratio_, lorentz_, length_, mass_, energy_] :=
	UnitSystem[Sequence@@GaussSystem[u, perm, ratio, 1, length, mass, energy], lorentz]
GaussSystem[u_String, perm_, ratio_, 1, length_, mass_, energy_] :=
	UnitSystem[Sequence@@EntropySystem[u, 1, length, mass, 1, perm, energy], ratio]

GaussSystem["EMU"] := GaussSystem[1, 4 Pi]
GaussSystem["ESU"] := GaussSystem[(100 "c")^-2, 4 Pi]
GaussSystem["Gauss"] := GaussSystem[1, 4 Pi, 1/100/"c"]
GaussSystem["LorentzHeaviside"] := GaussSystem[1, 1, 1/100/"c"]
GaussSystem["Thomson"] := GaussSystem[1, 4 Pi, 1/2]
GaussSystem["Kennelly"] := GaussSystem[10^-7, 4 Pi, 1, 1, 1]

EnergySystem[u_, time_, length_, mass_] := EntropySystem[u, time, length, mass, 1]
EnergySystem[u_, time_, length_, mass_, energy_] :=
	EntropySystem[u, time, length, mass, 1, PowerExpand[MagneticConstant[u] time^2/energy], energy]
EntropySystem[u_, time_, length_, mass_, temp_] := Module[{energy = PowerExpand[mass length^2/time^2]},
	EntropySystem[u, time, length, mass, temp, PowerExpand[MagneticConstant[u] time^2/energy]]]
EntropySystem[u_, time_, length_, mass_, temp_, perm_] :=
	EntropySystem[u, time, length, mass, temp, perm, PowerExpand[mass length^2/time^2]]
EntropySystem[u_, time_, length_, mass_, temp_, perm_, energy_] :=
	UnitSystem[BoltzmannConstant[u] temp/energy, ReducedPlanckConstant[u]/time/energy, time SpeedOfLight[u]/length, perm, ElectronMass[u]/mass]

EntropySystem["MTS"] := EnergySystem["AbstractMetric", 1, 1, 1000]
EntropySystem["IAU"] := EnergySystem["AbstractMetric", "day", "au", "GMsun"/AbstractUnitData["G"]]
EntropySystem["Hubble"] := EnergySystem["AbstractMetric", AbstractUnitData["th"], "c" AbstractUnitData["th"], 1]
EntropySystem["Cosmological"] := EnergySystem["AbstractMetric", AbstractUnitData["lc"]/"c", AbstractUnitData["lc"], AbstractUnitData["mc"]]
EntropySystem["CosmologicalQuantum"] := EnergySystem["AbstractMetric", AbstractUnitData["tcq"], AbstractUnitData["lcq"], AbstractUnitData["mcq"], AbstractUnitData["ecq"]]
EntropySystem["EMU2019"] := EnergySystem["AbstractSI2019", 1, 1/100, 1/1000]
EntropySystem["ESU2019"] := EntropySystem["AbstractSI2019", 1, 1/100, 1/1000, 1, 1000 AbstractUnitData["\[Mu]0"]/"c"^2]
EntropySystem["Mixed"] := EntropySystem["AbstractMetric", 1, 1, 1, 1, AbstractUnitData["\[Mu]0"]]
EntropySystem["English"] := EntropySystem["AbstractSI2019", 1, "ft", AbstractUnitData["slug"], "rankine", 4 Pi]
EntropySystem["EnglishUS"] := EntropySystem["AbstractMetric", 1, "ftUS", AbstractUnitData["slugUS"], "rankine", 4 Pi]
EntropySystem["FFF"] := EntropySystem["AbstractMetric", 14 "day", "fur", 90 "lb", "rankine", 0]

DimensionSystemQ[_] := False
DimensionSystemQ[UnitSystem[_, _, _, _, "M", ___]] := True
DimensionSystem[mu_, lore_] := UnitSystem[Sequence@@DimensionSystem[mu], 1, lore]
DimensionSystem[mu_] := UnitSystem["M" "L"^2/"T"^2/"\[CapitalTheta]", "M" "L"^2/"T", "L"/"T", mu, "M"]

DimensionSystem[u_String] := If[AbstractUnitSystemQ[u], DimensionSystem[StringDelete[u, "Abstract"]], DimensionSystem["M" "L"/"T"^2/"I"^2]]
DimensionSystem["ISQES"] := DimensionSystem["T"/"L"^2/"I"^2]
DimensionSystem["EMU"] := DimensionSystem[1]
DimensionSystem["ESU"] := DimensionSystem["T"^2/"L"^2]
DimensionSystem["Gauss"] := DimensionSystem[1, "T"/"L"]
DimensionSystem["Stoney"] := DimensionSystem["ISQEM"] //. {"\[CapitalTheta]" -> "M" "L"^2/"T"^2, "T" -> "L"}
DimensionSystem["Electronic"] := DimensionSystem["Stoney"] /. "M" -> 1
DimensionSystem["PlanckGauss"] := DimensionSystem["Stoney"] /. "L" -> 1/"M"
DimensionSystem["Planck"] := DimensionSystem["PlanckGauss"] /. "I" -> "M"
DimensionSystem["Natural"] := DimensionSystem["Planck"] /. "M" -> 1
DimensionSystem["NaturalGauss"] := DimensionSystem["PlanckGauss"] /. "M" -> 1
DimensionSystem["Rydberg"] := DimensionSystem["ISQEM"] //. {"\[CapitalTheta]" -> "M" "L"^2/"T"^2, "T" -> "L"^2 "M"}
DimensionSystem["Hartree"] := DimensionSystem["Rydberg"] /. "M" -> 1
DimensionSystem["Hubble"] := DimensionSystem["ISQEM"] /. "T" -> "L"
DimensionSystem["CosmologicalQuantum"] := DimensionSystem["Hubble"] /. "L" -> 1/"M"
DimensionSystem["ESU2019"] := DimensionSystem["ISQES"]
DimensionSystem["LorentzHeaviside"] := DimensionSystem["Gauss"]
DimensionSystem["HLU"] := DimensionSystem["Gauss"]
DimensionSystem["CGS"] := DimensionSystem["Gauss"]
DimensionSystem["CGSm"] := DimensionSystem["EMU"]
DimensionSystem["CGSe"] := DimensionSystem["ESU"]
DimensionSystem["Thomson"] := DimensionSystem["EMU"]
DimensionSystem["Kennelly"] := DimensionSystem["EMU"]
DimensionSystem["Schrodinger"] := DimensionSystem["Rydberg"]
DimensionSystem["QCD"] := DimensionSystem["Planck"]
DimensionSystem["QCDGauss"] := DimensionSystem["PlanckGauss"]
DimensionSystem["QCDoriginal"] := DimensionSystem["PlanckGauss"]
DimensionSystem["Cosmological"] := DimensionSystem["Hubble"]

AbstractUnitSystem = <|
"AbstractUnits" -> UnitSystem["kB","\[HBar]","c","\[Mu]0","me","\[Lambda]","\[Alpha]L"],
"AbstractUnits1" -> UnitSystem["kB1", "\[HBar]1", "c1", "\[Mu]01", "me1", "\[Lambda]1", "\[Alpha]L1"],
"AbstractUnits2" -> UnitSystem["kB2", "\[HBar]2", "c2", "\[Mu]02", "me2", "\[Lambda]2", "\[Alpha]L2"]|>

AppendTo[AbstractUnitSystem, Map[(StringJoin["Dimension",#] -> DimensionSystem[#]) &,
{"ISQEM", "ISQES", "EMU", "ESU", "Gauss", "Stoney", "Electronic", "PlanckGauss", "Planck", "Natural", "NaturalGauss", "Rydberg", "Hartree", "Hubble", "CosmologicalQuantum"}]]

AppendSystems[f_, l_] := AppendTo[AbstractUnitSystem, Map[(# -> f[#]) &, l]]
AppendSystems[MetricSystem, {"Metric", "SI2019", "SI1976"}]
AppendSystems[ConventionalSystem, {"CODATA", "Conventional"}]
AppendSystems[GaussSystem, {"EMU", "ESU", "Gauss", "LorentzHeaviside", "Thomson", "Kennelly"}]
AppendSystems[EntropySystem, {"MTS", "IAU", "Hubble", "Cosmological", "CosmologicalQuantum", "EMU2019", "ESU2019", "Mixed", "English", "EnglishUS", "FFF"}]

AppendTo[AbstractUnitSystem, {
"Planck" -> UnitSystem[1, 1, 1, 1, PowerExpand[Sqrt[4 Pi AbstractUnitData["\[Alpha]G"]]]],
"PlanckGauss" -> UnitSystem[1, 1, 1, 4 Pi, PowerExpand[Sqrt[AbstractUnitData["\[Alpha]G"]]]],
"Stoney" -> UnitSystem[1, 1/"\[Alpha]", 1, 4 Pi, PowerExpand[Sqrt[AbstractUnitData["\[Alpha]G"]/"\[Alpha]"]]],
"Hartree" -> UnitSystem[1,1,1/"\[Alpha]",4 Pi "\[Alpha]"^2,1],
"Rydberg" -> UnitSystem[1,1,2/"\[Alpha]",Pi "\[Alpha]"^2,1/2],
"Schrodinger" -> UnitSystem[1, 1, 1/"\[Alpha]", 4 Pi "\[Alpha]"^2, PowerExpand[Sqrt[AbstractUnitData["\[Alpha]G"]/"\[Alpha]"]]],
"Electronic" -> UnitSystem[1, 1/"\[Alpha]", 1, 4 Pi, 1],
"Natural" -> UnitSystem[1, 1, 1, 1, 1, 1, 1, "1"],
"NaturalGauss" -> UnitSystem[1, 1, 1, 4 Pi, 1, 1, 1, "1"],
"QCD" -> UnitSystem[1, 1, 1, 1, 1/AbstractUnitData["\[Mu]pe"]],
"QCDGauss" -> UnitSystem[1, 1, 1, 4 Pi, 1/AbstractUnitData["\[Mu]pe"]],
"QCDoriginal" -> UnitSystem[1, 1, 1, 4 Pi "\[Alpha]", 1/AbstractUnitData["\[Mu]pe"]],
"SI" :> AbstractUnitSystem["SI2019"],
"MKS" :> AbstractUnitSystem["Metric"],
"CGS" :> AbstractUnitSystem["Gauss"],
"CGS2019" :> AbstractUnitSystem["EMU2019"],
"CGSm" :> AbstractUnitSystem["EMU"],
"CGSe" :> AbstractUnitSystem["ESU"],
"HLU" :> AbstractUnitSystem["LorentzHeaviside"],
"ISQ" :> AbstractUnitSystem["DimensionISQ"],
"DimensionISQ" :> AbstractUnitSystem["DimensionISQEM"]}]

Map[(UnitSystem[#] := AbstractUnitSystem[#]) &,
{"AbstractUnits","AbstractUnits1","AbstractUnits2","ISQ","SI","MKS","CGS","CGS2019","CGSm","CGSe","HLU"}]
UnitSystem["Natural"] = UnitSystem[1, 1, 1, 1, 1]
UnitSystem["NaturalGauss"] = UnitSystem[1, 1, 1, 4 Pi, 1]
Map[(UnitSystem[#] = AbstractUnitSystem[#] /. Normal[UnitData]) &,
{"Gauss","LorentzHeaviside","Thomson","Kennelly","ESU","ESU2019","EMU","EMU2019","MTS","Mixed","Metric","SI1976","SI2019","CODATA","Conventional","English","EnglishUS","IAU","FFF","Planck","PlanckGauss","Stoney","Hartree","Rydberg","Schrodinger","Electronic","QCD","QCDGauss","QCDoriginal","Hubble","Cosmological","CosmologicalQuantum"}]
