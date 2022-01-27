
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

AbstractUnitData = <||>
AppendTo[AbstractUnitData, "slug" -> "lb" "g0"/"ft"]
AppendTo[AbstractUnitData, "slugUS" -> "lb" "g0"/"ftUS"]
AppendTo[AbstractUnitData, "Ru" -> "NA" "kB"]
AppendTo[AbstractUnitData, "me" -> 2 "R\[Infinity]" "h"/"\[Alpha]"^2/"c"]
AppendTo[AbstractUnitData, "\[Mu]0" -> 2 "\[Alpha]" "h"/"c"/"e"^2]
AppendTo[AbstractUnitData, "\[HBar]" -> "h"/(2 Pi)]
AppendTo[AbstractUnitData, "\[Mu]pe" -> "\[Mu]pu"/"\[Mu]eu"]
AppendTo[AbstractUnitData, "\[Alpha]G" -> (AbstractUnitData["me"]/"mP")^2]
AppendTo[AbstractUnitData, "G" -> "c" AbstractUnitData["\[HBar]"]/"mP"^2]
AppendTo[AbstractUnitData, "Mu" -> "NA" AbstractUnitData["me"]/"\[Mu]eu"]
AppendTo[AbstractUnitData, "pc" -> "au" 3*60^3/Pi]
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

ConventionalSystem[klitz_, joseph_] := ConventionalSystem[klitz, Nothing, 2/klitz/joseph^2/Pi]
ConventionalSystem[klitz_, _, hbar_] := ConventionalSystem[klitz, Nothing, hbar, 4 Pi "R\[Infinity]" hbar/"\[Alpha]"^2/"c"]
ConventionalSystem[klitz_, _, hbar_, mass_] := MetricSystem[1/1000, 2 "\[Alpha]" klitz/"c", AbstractUnitData["Ru"], hbar, mass]

ConventionalSystem["CODATA"] := ConventionalSystem["RK2014", "KJ2014"]
ConventionalSystem["Conventional"] := ConventionalSystem["RK1990", "KJ1990"]

GaussSystem[u_String] := ($Failed /; False)
GaussSystem[perm_, u__] := GaussSystem["AbstractMetric", perm, u]
GaussSystem[u_String, perm_, ratio_] := GaussSystem[u, perm, ratio, 1]
GaussSystem[u_String, perm_, ratio_, lorentz_] := GaussSystem[u, perm, ratio, lorentz, 1/100, 1/1000]
GaussSystem[u_String, perm_, ratio_, lorentz_, length_, mass_] :=
	GaussSystem[u, perm, ratio, lorentz, length, mass, PowerExpand[mass length^2]]
GaussSystem[u_String, perm_, ratio_, 1, length_, mass_, energy_] :=
	UnitSystem[BoltzmannConstant[u]/energy, ReducedPlanckConstant[u]/energy, SpeedOfLight[u]/length, perm, ElectronMass[u]/mass, ratio]
GaussSystem[u_String, perm_, ratio_, lorentz_, length_, mass_, energy_] :=
	UnitSystem[BoltzmannConstant[u]/energy, ReducedPlanckConstant[u]/energy, SpeedOfLight[u]/length, perm, ElectronMass[u]/mass, ratio, lorentz]

GaussSystem["Gauss"] := GaussSystem[1, 4 Pi, 1/100/"c"]
GaussSystem["LorentzHeaviside"] := GaussSystem[1, 1, 1/100/"c"]
GaussSystem["Thomson"] := GaussSystem[1, 4 Pi, 1/2]
GaussSystem["Kennelly"] := GaussSystem[10^-7, 4 Pi, 1, 1, 1]
GaussSystem["ESU"] := GaussSystem[(100 "c")^-2, 4 Pi]
GaussSystem["EMU"] := GaussSystem[1, 4 Pi]

EnergySystem[u_, time_, length_, mass_] := EntropySystem[u, time, length, mass, 1]
EnergySystem[u_, time_, length_, mass_, energy_] :=
	EntropySystem[u, time, length, mass, 1, PowerExpand[MagneticConstant[u] time^2/energy], energy]
EntropySystem[u_, time_, length_, mass_, temp_] := Module[{energy = PowerExpand[mass length^2/time^2]},
	EntropySystem[u, time, length, mass, temp, PowerExpand[MagneticConstant[u] time^2/energy]]]
EntropySystem[u_, time_, length_, mass_, temp_, perm_] :=
	EntropySystem[u, time, length, mass, temp, perm, PowerExpand[mass length^2/time^2]]
EntropySystem[u_, time_, length_, mass_, temp_, perm_, energy_] :=
	UnitSystem[BoltzmannConstant[u] temp/energy, ReducedPlanckConstant[u]/time/energy, time SpeedOfLight[u]/length, perm, ElectronMass[u]/mass]

EnergySystem["MTS"] := EnergySystem["AbstractMetric", 1, 1, 1000]
EnergySystem["IAU"] := EnergySystem["AbstractMetric", "day", "au", "GMsun"/AbstractUnitData["G"]]
EnergySystem["Hubble"] := EnergySystem["AbstractMetric", AbstractUnitData["th"], "c" AbstractUnitData["th"], 1]
EnergySystem["Cosmological"] := EnergySystem["AbstractMetric", AbstractUnitData["lc"]/"c", AbstractUnitData["lc"], AbstractUnitData["mc"]]
EnergySystem["CosmologicalQuantum"] := EnergySystem["AbstractMetric", AbstractUnitData["tcq"], AbstractUnitData["lcq"], AbstractUnitData["mcq"], AbstractUnitData["ecq"]]
EnergySystem["EMU2019"] := EnergySystem["AbstractSI2019", 1, 1/100, 1/1000]
EntropySystem["ESU2019"] := EntropySystem["AbstractSI2019", 1, 1/100, 1/1000, 1, 1000 AbstractUnitData["\[Mu]0"]/"c"^2]
EntropySystem["Mixed"] := EntropySystem["AbstractMetric", 1, 1, 1, 1, AbstractUnitData["\[Mu]0"]]
EntropySystem["English"] := EntropySystem["AbstractSI2019", 1, "ft", AbstractUnitData["slug"], "rankine", 4 Pi]
EntropySystem["EnglishUS"] := EntropySystem["AbstractMetric", 1, "ftUS", AbstractUnitData["slugUS"], "rankine", 4 Pi]
EntropySystem["FFF"] := EntropySystem["AbstractMetric", 14 "day", "fur", 90 "lb", "rankine", 0]

AbstractUnitSystem = <|
"AbstractUnits" -> UnitSystem["kB","\[HBar]","c","\[Mu]0","me","\[Lambda]","\[Alpha]L"],
"AbstractUnits1" -> UnitSystem["kB1", "\[HBar]1", "c1", "\[Mu]01", "me1", "\[Lambda]1", "\[Alpha]L1"],
"AbstractUnits2" -> UnitSystem["kB2", "\[HBar]2", "c2", "\[Mu]02", "me2", "\[Lambda]2", "\[Alpha]L2"],
"SI2019" -> MetricSystem["SI2019"],
"Metric" -> MetricSystem["Metric"],
"SI1976" -> MetricSystem["SI1976"],
"CODATA" -> ConventionalSystem["CODATA"],
"Conventional" -> ConventionalSystem["Conventional"]|>

AppendTo[AbstractUnitSystem, {
"Gauss" -> GaussSystem["Gauss"],
"LorentzHeaviside" -> GaussSystem["LorentzHeaviside"],
"Thomson" -> GaussSystem["Thomson"],
"Kennelly" -> GaussSystem["Kennelly"],
"ESU" -> GaussSystem["ESU"],
"EMU" -> GaussSystem["EMU"],
"MTS" -> EnergySystem["MTS"],
"IAU" -> EnergySystem["IAU"],
"Hubble" -> EnergySystem["Hubble"],
"Cosmological" -> EnergySystem["Cosmological"],
"CosmologicalQuantum" -> EnergySystem["CosmologicalQuantum"],
"EMU2019" -> EnergySystem["EMU2019"],
"ESU2019" -> EntropySystem["ESU2019"],
"Mixed" -> EntropySystem["Mixed"],
"English" -> EntropySystem["English"],
"EnglishUS" -> EntropySystem["EnglishUS"],
"FFF" -> EntropySystem["FFF"],
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
"HLU" :> AbstractUnitSystem["LorentzHeaviside"]}]

UnitSystem["AbstractUnits"] := AbstractUnitSystem["AbstractUnits"]
UnitSystem["AbstractUnits1"] := AbstractUnitSystem["AbstractUnits1"]
UnitSystem["AbstractUnits2"] := AbstractUnitSystem["AbstractUnits2"]
UnitSystem["SI"] := UnitSystem["SI2019"]
UnitSystem["MKS"] := UnitSystem["Metric"]
UnitSystem["CGS"] := UnitSystem["Gauss"]
UnitSystem["CGS2019"] := UnitSystem["EMU2019"]
UnitSystem["CGSm"] := UnitSystem["EMU"]
UnitSystem["CGSe"] := UnitSystem["ESU"]
UnitSystem["HLU"] := UnitSystem["LorentzHeaviside"]
UnitSystem["Natural"] = UnitSystem[1, 1, 1, 1, 1]
UnitSystem["NaturalGauss"] = UnitSystem[1, 1, 1, 4 Pi, 1]
Map[(UnitSystem[#] = AbstractUnitSystem[#] /. Normal[UnitData]) &,
{"Gauss","LorentzHeaviside","Thomson","Kennelly","ESU","ESU2019","EMU","EMU2019","MTS","Mixed","Metric","SI1976","SI2019","CODATA","Conventional","English","EnglishUS","IAU","FFF","Planck","PlanckGauss","Stoney","Hartree","Rydberg","Schrodinger","Electronic","QCD","QCDGauss","QCDoriginal","Hubble","Cosmological","CosmologicalQuantum"}]
