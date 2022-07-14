
(*   This file is part of UnitSystems                   *)
(*   It is licensed under the MIT license               *)
(*   UnitSystems Copyright (C) 2021 Michael Reed        *)
(*       _           _                         _        *)
(*      | |         | |                       | |       *)
(*   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _  *)
(*  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` | *)
(* | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| | *)
(*  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_| *)

UnitData = <|
"g0" -> 980665/10^5,
"ft" -> 3048/10000,
"ftUS" -> 1200/3937,
"lb" -> 45359237/10^8,
"\[Degree]R" -> 5/9,
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
"KJ1990" -> 4835979 10^8,
"RK2014" -> Around[25812.8074555, 59/10^7],
"KJ2014" -> Around[4.835978525 10^14, 3 10^6],
"Ru2014" -> Around[8.3144598,48/10^7],
"GME" -> Around[3.986004418 10^14, 8 10^5],
"GMJ" -> Aronund[1.26686534 10^17, 9 10^9],
"day" -> 60^2 24,
"aj" -> 36525/100,
"au" -> Around[149597870700, 3],
"fur" -> 201168/1000,
"kG" -> 354818761/10^5,
"H0" -> Around[67.66, 0.42],
"\[CapitalOmega]\[CapitalLambda]" -> Around[0.6889, 0.0056],
"\[CapitalOmega]it" -> 1000495/10^6,
"Vit" -> 100033/10^5,
"ks" -> 5778/10^4,
"T0" -> 27315/100,
"inHg" -> 1000/3386389|>

AbstractUnitData = <|
"slug" -> "lb" "g0"/"ft",
"slugUS" -> "lb" "g0"/"ftUS",
"Ru" -> "NA" "kB",
"me" -> 2 "R\[Infinity]" "h"/"\[Alpha]"^2/"c",
"\[Mu]0" -> 2 "\[Alpha]" "h"/"c"/"e"^2,
"\[HBar]" -> "h"/(2 Pi),
"\[Mu]pe" -> "\[Mu]pu"/"\[Mu]eu",
"k" -> "kG" Pi/(2^6 3^4 5^3),
"em" -> Sqrt["GME"/"g0"] Pi/2^8/5^7,
"nm" -> Sqrt["GME"/"g0"] Pi/2^4/3^3/5^2,
"pc" -> "au" 3*60^3/Pi|>

AppendTo[AbstractUnitData, "GMsun" -> "au"^3 AbstractUnitData["k"]^2/"day"^2]
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
AppendTo[AbstractUnitData, "tcq" -> Module[{ecq = PowerExpand[("c"^3 AbstractUnitData["\[HBar]"]^3 AbstractUnitData["\[Rho]\[CapitalLambda]"])^(1/4)]},
	PowerExpand[AbstractUnitData["lcq"] Sqrt[AbstractUnitData["mcq"]/ecq]]]]

EvalUnitData[x] := AbstractUnitData[x]/.Normal[UnitData]

Coupling["AbstractCoupling"] = Coupling["\[Alpha]G", "\[Alpha]", "\[Mu]eu", "\[Mu]pu", "\[CapitalOmega]\[CapitalLambda]"]
Coupling["AbstractUniverse"] = Coupling[AbstractUnitData["\[Alpha]G"], "\[Alpha]", "\[Mu]eu", "\[Mu]pu", "\[CapitalOmega]\[CapitalLambda]"]
Coupling["StandardModel"] = Coupling["AbstractUniverse"] /. Normal[UnitData]

MM[Measure[v_,_,_]] := v
MM[v_] := v
swap[u_,s_String] := SWAP[u,StringJoin["Abstract",s]]
SWAP[Measure[v_, d_, _],s_String] := Measure[v,d,s]
SWAP[UnitSystem[kB_,h_,c_,mu_,me_,Mu_,Kcd_,a_,r_,l_,g_,t_,u_],s_] := UnitSystem[
	SWAP[kB,s],SWAP[h,s],SWAP[c,s],SWAP[mu,s],SWAP[me,s],SWAP[Mu,s],SWAP[Kcd,s],SWAP[a,s],SWAP[r,s],SWAP[l,s],SWAP[g,s],SWAP[t,s],u]
MeasureSystem[n_,kB_,h_,c_,mu_,me_,Mu_:1,Kcd_:1,a_:1,r_:1,l_:1,g_:1,t_:(2 Pi),u_:Coupling["AbstractUniverse"]] = UnitSystem[
	Measure[MeasureMagnitude[kB], dF*dL/d0, n],
	Measure[MeasureMagnitude[h], dF*dL*dT/dA, n],
	Measure[MeasureMagnitude[c], dL/dT, n],
	Measure[MeasureMagnitude[mu], dF*dT^2/dQ^2/dR*dC^2, n],
	Measure[MeasureMagnitude[me], dM, n],
	Measure[MeasureMagnitude[Mu], dM/dN, n],
	Measure[MeasureMagnitude[Kcd], dT*dJ/dF/dL, n],
	Measure[MeasureMagnitude[a], d1, n],
	Measure[MeasureMagnitude[r], dR, n],
	Measure[MeasureMagnitude[l], d1/dC, n],
	Measure[MeasureMagnitude[g], dM*dL/dF/dT^2, n],
	Measure[MeasureMagnitude[t], dA, n],u]

MetricSystem[name_, mu_, perm_, ru_:AbstractUnitData["Ru"], g0_:1] := MetricSystem[name, mu, perm, ru, g0, "h", AbstractUnitData["me"]];
MetricSystem[name_, mu_, perm_, ru_, g0_, h_, mass_] := MeasureSystem[name, ru mass/mu/"\[Mu]eu"/g0, h/2/Pi/g0, "c", perm, mass, mu, "Kcd"*(AbstractUnitData["me"]/mass)^2*(h/"h"), 1, 1, 1, g0, 2 Pi, Coupling["AbstractUniverse"]];
ConventionalSystem[name_, klitz_, joseph_, Ru_:AbstractUnitData["Ru"], g0_:1] := Module[{h = 4/klitz/joseph^2},
	MetricSystem[name, 1/1000, 2 "\[Alpha]" klitz/"c", Ru, g0, h, 2 "R\[Infinity]" h/"\[Alpha]"^2/"c"]]
RankineSystem[u_, s_, l_, m_, g0_:1] := EntropySystem[u,s,1,l,m,"\[Degree]R",MM@MagneticConstant[u]/m/l/g0, 1000 MM@MolarMassConstant[u],g0]

MetricSystem["SI2019"] := MetricSystem["AbstractSI2019",AbstractUnitData["Mu"], AbstractUnitData["\[Mu]0"]]
MetricSystem["Metric"] := MetricSystem["AbstractMetric",1/1000, 4 Pi/10^7]
MetricSystem["SI2019Engineering"] := MetricSystem["AbstractSI2019Engineering",AbstractUnitData["Mu"], AbstractUnitData["\[Mu]0"]/"g0", AbstractUnitData["Ru"], "g0"]
MetricSystem["MetricEngineering"] := MetricSystem["AbstractMetricEngineering",1/1000, 4 Pi/10^7/"g0", AbstractUnitData["Ru"], "g0"]
MetricSystem["SI1976"] := MetricSystem["AbstractSI1976",1/1000, 4 Pi/10^7, 831432/10^5]
MetricSystem["CODATA"] := ConventionalSystem["AbstractCODATA", "RK2014", "KJ2014", "Ru2014"]
MetricSystem["Conventional"] := ConventionalSystem["AbstractConventional", "RK1990", "KJ1990"]

RankineSystem["IPS"] := RankineSystem["AbstractMetric","IPS","ft"/12,"lb" "g0" 12/"ft"]
RankineSystem["British"] := RankineSystem["AbstractMetric","British","ft",AbstractUnitData["slug"]]
(*RankineSystem["British2019"] := RankineSystem["AbstractSI2019","British2019","ft",AbstractUnitData["slug"]]*)
RankineSystem["Survey"] := RankineSystem["AbstractMetric","Survey","ftUS","lb","g0"/"ftUS"]
(*RankineSystem["Survey2019"] := RankineSystem["AbstractSI2019","Survey2019","ftUS","lb","g0"/"ftUS"]*)
RankineSystem["English"] := RankineSystem["AbstractMetric","English","ft","lb","g0"/"ft"]
(*RankineSystem["English2019"] := RankineSystem["AbstractSI2019","English2019","ft","lb","g0"/"ft"]*)
RankineSystem["FPS"] := RankineSystem["AbstractMetric","FPS","ft","lb"]
(*RankineSystem["FPS2019"] := RankineSystem["AbstractSI2019","FPS2019","ft","lb"]*)

AstronomicalSystem[u_, s_] := EntropySystem[u,s,1,1,1/MM@GravitationalConstant[u]]
ElectricSystem[u_,s_,r_,v_] := EntropySystem[u,s,1,1,v^2/r,1,MM@MagneticConstant[u]/r]

GaussSystem[u_String, s_, perm_, ratio_] := GaussSystem[u, s, perm, ratio, 1]
GaussSystem[u_String, s_, perm_, ratio_, lorentz_] := GaussSystem[u, s, perm, ratio, lorentz, 1/100, 1/1000]
GaussSystem[u_String, s_, perm_, ratio_, lorentz_, length_, mass_] :=
	GaussSystem[u, s, perm, ratio, lorentz, length, mass, MM@GravityConstant[u]]
GaussSystem[u_String, s_, perm_, ratio_, lorentz_, length_, mass_, g0_] :=
	EntropySystem[u, s, 1, length, mass, 1, perm, MM@MolarMassConstant[u]/mass, g0, PowerExpand[mass length^2], ratio, lorentz]

GaussSystem["EMU"] := GaussSystem["AbstractMetric", "EMU", 1, 4 Pi]
GaussSystem["ESU"] := GaussSystem["AbstractMetric", "ESU", (100 "c")^-2, 4 Pi]
GaussSystem["Gauss"] := GaussSystem["AbstractMetric", "Gauss", 1, 4 Pi, 1/100/"c"]
GaussSystem["LorentzHeaviside"] := GaussSystem["AbstractMetric", "LorentzHeaviside", 1, 1, 1/100/"c"]
(*GaussSystem["Thomson"] := GaussSystem["AbstractMetric", "Thomson", 1, 4 Pi, 1/2]*)
GaussSystem["Kennelly"] := GaussSystem["AbstractMetric", "Kennelly", 10^-7, 4 Pi, 1, 1, 1]

EntropySystem[u_, s_, time_, length_, mass_, temp_:1] :=
	EntropySystem[u, s, time, length, mass, temp, PowerExpand[MM@MagneticConstant[u]/mass/length],MM@MolarMassConstant[u]/mass, MM@GravityConstant[u], PowerExpand[mass length^2/time^2]]
EntropySystem[u_, s_, time_, length_, mass_, temp_, perm_] :=
	EntropySystem[u, s, time, length, mass, temp, perm, MM@MolarMassConstant[u]]
EntropySystem[u_, s_, time_, length_, mass_, temp_, perm_, mol_] :=
	EntropySystem[u, s, time, length, mass, temp, perm, mol, MM@GravityConstant[u]]
EntropySystem[u_, s_, time_, length_, mass_, temp_, perm_, mol_, g0_] :=
	EntropySystem[u, s, time, length, mass, temp, perm, mol, g0, PowerExpand[mass length^2/time^2]]
EntropySystem[u_, s_, time_, length_, mass_, temp_, perm_, mol_, g0_, energy_, rat_:1, lor_:1] :=
	MeasureSystem[StringJoin["Abstract",s], MM@BoltzmannConstant[u] temp/energy/g0, MM@ReducedPlanckConstant[u]/time/energy/g0, time MM@SpeedOfLight[u]/length, perm, MM@ElectronMass[u]/mass, mol, MM@MonochromaticRadiation540THzLuminousEfficacy[u]*energy/time g0, MM@AngleConstant[u], rat, lor, g0, MM@Turn[u], Universe[u]]

(*EntropySystem["Astronomical"] := AstronomicalSystem["AbstractMetric","Astronomical"]*)
EntropySystem["International"] := ElectricSystem["AbstractMetric","International","\[CapitalOmega]it","Vit"]
EntropySystem["InternationalMean"] := ElectricSystem["AbstractMetric","InternationalMean",100049/10^5, 100034/10^5]

EntropySystem["Nautical"] := EntropySystem["AbstractMetric", "Nautical", 60^2, AbstractUnitData["nm"], AbstractUnitData["em"]^3, 1, 27 Pi/2^9/5^12,1/1000]
EntropySystem["Meridian"] := EntropySystem["AbstractMetric", "Meridian", 1, AbstractUnitData["em"], AbstractUnitData["em"]^3, 1, Pi/2^5/5^7,1/1000]
EntropySystem["MeridianEngineering"] := EntropySystem["AbstractMetric", "MeridianEngineering", 1, AbstractUnitData["em"], AbstractUnitData["em"]^3, 1, Pi/2^5/5^7 AbstractUnitData["em"]/"g0",1/1000, "g0"/AbstractUnitData["em"]]
EntropySystem["MPH"] := EntropySystem["AbstractFPS", "MPH", 60^2, 5280, 1]
EntropySystem["KKH"] := EntropySystem["AbstractMetric", "KKH", 60^2, 1000, 1]
EntropySystem["MTS"] := EntropySystem["AbstractMetric", "MTS", 1, 1, 1000]
EntropySystem["IAU"] := EntropySystem["AbstractMetric", "IAU", "day", "au", AbstractUnitData["GMsun"]/AbstractUnitData["G"]]
EntropySystem["IAUE"] := EntropySystem["AbstractMetric", "IAUE", "day", "au", "GME"/AbstractUnitData["G"]]
EntropySystem["IAUJ"] := EntropySystem["AbstractMetric", "IAUJ", "day", "au", "GMJ"/AbstractUnitData["G"]]
EntropySystem["Hubble"] := EntropySystem["AbstractMetric", "Hubble", AbstractUnitData["th"], "c" AbstractUnitData["th"], 1]
EntropySystem["Cosmological"] := EntropySystem["AbstractMetric", "Cosmological", AbstractUnitData["lc"]/"c", AbstractUnitData["lc"], AbstractUnitData["mc"]]
EntropySystem["CosmologicalQuantum"] := EntropySystem["AbstractMetric", "CosmologicalQuantum", AbstractUnitData["tcq"], AbstractUnitData["lcq"], AbstractUnitData["mcq"]]
(*EntropySystem["EMU2019"] := EntropySystem["AbstractSI2019", "EMU2019", 1, 1/100, 1/1000]
EntropySystem["ESU2019"] := EntropySystem["AbstractSI2019", "ESU2019", 1, 1/100, 1/1000, 1, 1000 AbstractUnitData["\[Mu]0"]/"c"^2]
EntropySystem["Mixed"] := EntropySystem["AbstractMetric", "Mixed", 1, 1, 1, 1, AbstractUnitData["\[Mu]0"]]*)
EntropySystem["GravitationalMetric"] := EntropySystem["AbstractMetric", "GravitationalMetric", 1, 1, "g0"]
EntropySystem["GravitationalSI2019"] := EntropySystem["AbstractSI2019", "GravitationalSI2019", 1, 1, "g0"]
EntropySystem["GravitationalMeridian"] := EntropySystem["AbstractMetric", "GravitationalMeridian", 1, AbstractUnitData["em"], "g0" AbstractUnitData["em"]^2, 1, Pi/2^5/5^7 AbstractUnitData["em"]/"g0",1/1000]
EntropySystem["FFF"] := EntropySystem["AbstractMetric", "FFF", 14 "day", "fur", 90 "lb", "\[Degree]R", 0, 1]

(*DimensionSystem[] := UnitSystem[dF dL/d0, dF dL/dT/dA, dL/dT, dF dT^2/dQ^2/dA^2/dR dC^2, dM, dM/dN, dJ dT/dL/dF, d1, dR*dA^2, d1/dC, dM dL/dT^2/dF, dA, Coupling["AbstractUniverse"]]*)
DimensionSystem[] := UnitSystem["F" "L"/"\[CapitalTheta]", "F" "L" "T"/"A", "L"/"T", "F" "T"^2/"Q"^2/"A"^2/"\[CapitalLambda]" "C"^2, "M", "M"/"N", "J" "T"/"L"/"F", 1, "\[CapitalLambda]"*"A"^2, 1/"C", "M" "L"/"T"^2/"F", 2 Pi "A", Coupling["AbstractUniverse"]]

DimensionSystemQ[_] = False
DimensionSystemQ[UnitSystem[_, _, _, _, "M", ___]] := True

DimensionRules[u_String,x_String] := DimensionRules[u,DimensionSystem[x]]

DimensionRules["MetricEngineering", x_] := x /. {"C" -> 1, "\[CapitalLambda]" -> 1}
DimensionRules["GravitationalMetric", x_] := DimensionRules["MetricEngineering", x] /. {"M" -> "F" "T"^2/"L", "A" -> 1}
DimensionRules["ISQEM", x_] := DimensionRules["MetricEngineering", x] /. {"F" -> "M" "L"/"T"^2, "A" -> 1}
DimensionRules["EMU", x_] := x /. "Q" -> Sqrt["L" "M"]
DimensionRules["ESU", x_] := x /. "Q" -> Sqrt["L"^3 "M"]/"T"
DimensionRules["Gauss", x_] := x /. "Q" -> Sqrt["L"^3 "M"]/"T"
DimensionRules["Stoney",x_] := x //. {"\[CapitalTheta]" -> "M", "N" -> "M", "T" -> "L"}
DimensionRules["Electronic",x_] := DimensionRules["Stoney", x] /. "M" -> 1
DimensionRules["PlanckGauss",x_] := DimensionRules["Stoney", x] /. "L" -> 1/"M"
DimensionRules["Planck",x_] := DimensionRules["PlanckGauss", x] /. "Q" -> "M"
(*DimensionRules["QCDoriginal",x_] := DimensionRules["PlanckGauss", x] /. "Q" -> 1*)
DimensionRules["NaturalGauss",x_] := DimensionRules["PlanckGauss", x] /. "M" -> 1
DimensionRules["Rydberg",x_] := x //. {"\[CapitalTheta]" -> "M" "L"^2/"T"^2, "T" -> "L"^2 "M"}
DimensionRules["Hartree",x_] := DimensionRules["Rydberg", x] /. "M" -> 1
DimensionRules["Hubble",x_] := x /. "T" -> "L"
DimensionRules["CosmologicalQuantum",x_] := DimensionRules["Hubble", x] /. "L" -> 1/"M"

DimensionSystem[u_String] := If[AbstractUnitSystemQ[u], DimensionSystem[StringDelete[u, "Abstract"]], DimensionRules["ISQEM",DimensionSystem[]]]
(*DimensionSystem["Rationalized"] := DimensionSystem["L" "M"/("I" "T" "C")^2/"\[CapitalLambda]", "\[CapitalLambda]", "C"]*)
(*DimensionSystem["ISQEM"] := DimensionSystem["M" "L"/"T"^2/"I"^2]*)
(*DimensionSystem["ISQES"] := DimensionSystem["T"/"L"^2/"I"^2]*)
DimensionSystem["MetricEngineering"] := DimensionRules["MetricEngineering",DimensionSystem[]]
DimensionSystem["GravitationalMetric"] := DimensionRules["GravitationalMetric",DimensionSystem[]]
DimensionSystem["EMU"] := DimensionRules["EMU","ISQEM"]
DimensionSystem["ESU"] := DimensionRules["ESU","ISQEM"]
DimensionSystem["Gauss"] := DimensionRules["Gauss","ISQEM"]
DimensionSystem["Stoney"] := DimensionRules["Stoney","ISQEM"]
DimensionSystem["Electronic"] := DimensionRules["Electronic","Stoney"]
DimensionSystem["PlanckGauss"] := DimensionRules["PlanckGauss","Stoney"]
DimensionSystem["Planck"] := DimensionRules["Planck","PlanckGauss"]
(*DimensionSystem["QCDoriginal"] := DimensionRules["QCDoriginal","PlanckGauss"]*)
DimensionSystem["QCDoriginal"] := DimensionSystem["PlanckGauss"]
DimensionSystem["Natural"] := AbstractUnitSystem["Natural"]
DimensionSystem["NaturalGauss"] := DimensionRules["NaturalGauss","PlanckGauss"]
DimensionSystem["Rydberg"] := DimensionRules["Rydberg","ISQEM"]
DimensionSystem["Hartree"] := DimensionRules["Hartree","Rydberg"]
DimensionSystem["Hubble"] := DimensionRules["Hubble","ISQEM"]
DimensionSystem["CosmologicalQuantum"] := DimensionRules["CosmologicalQuantum","Hubble"]
(*DimensionSystem["ESU2019"] := DimensionSystem["ISQES"]*)
DimensionSystem["LorentzHeaviside"] := Dimension["Gauss"]
DimensionSystem["HLU"] := DimensionSystem["Gauss"]
DimensionSystem["CGS"] := DimensionSystem["Gauss"]
DimensionSystem["CGSm"] := DimensionSystem["EMU"]
DimensionSystem["CGSe"] := DimensionSystem["ESU"]
(*DimensionSystem["Thomson"] := DimensionSystem["EMU"]*)
DimensionSystem["Kennelly"] := DimensionSystem["EMU"]
DimensionSystem["Schrodinger"] := DimensionSystem["Rydberg"]
DimensionSystem["QCD"] := DimensionSystem["Planck"]
DimensionSystem["QCDGauss"] := DimensionSystem["PlanckGauss"]
DimensionSystem["Cosmological"] := DimensionSystem["Hubble"]
DimensionSystem["SI2019Engineering"] := DimensionSystem["MetricEngineering"]
DimensionSystem["MeridianEngineering"] := DimensionSystem["MetricEngineering"]
DimensionSystem["GravitationalSI2019"] := DimensionSystem["GravitationalMetric"]
DimensionSystem["GravitationalMeridian"] := DimensionSystem["GravitationalMetric"]
DimensionSystem["IPS"] := DimensionSystem["GravitationalMetric"]
DimensionSystem["British"] := DimensionSystem["GravitationalMetric"]
(*DimensionSystem["British2019"] := DimensionSystem["GravitationalMetric"]*)
DimensionSystem["English"] := DimensionSystem["MetricEngineering"]
(*DimensionSystem["English2019"] := DimensionSystem["MetricEngineering"]*)
DimensionSystem["Survey"] := DimensionSystem["MetricEngineering"]
(*DimensionSystem["Survey2019"] := DimensionSystem["MetricEngineering"]*)

DimensionSystemQ[DimensionSystem["Hartree"]] = True
DimensionSystemQ[DimensionSystem["Electronic"]] = True
DimensionSystemQ[DimensionSystem["NaturalGauss"]] = True
DimensionSystemQ[DimensionSystem["Natural"]] = True

AbstractUnitSystem = <|
"AbstractUnits" -> UnitSystem["kB","\[HBar]","c","\[Mu]0","me","\[Lambda]","\[Alpha]L"],
"AbstractUnits1" -> UnitSystem["kB1", "\[HBar]1", "c1", "\[Mu]01", "me1", "\[Lambda]1", "\[Alpha]L1"],
"AbstractUnits2" -> UnitSystem["kB2", "\[HBar]2", "c2", "\[Mu]02", "me2", "\[Lambda]2", "\[Alpha]L2"],
"Unified" -> MeasureSystem["AbstractUnified",dF*dL/d0,dF*dL*dT/dA,dL/dT,dF*dT^2/dQ^2/dR*dC^2,dM,dM/dN,dT*dJ/dF/dL,dA,dR,d1/dC,dM*dL/dT^2/dF]|>

AppendTo[AbstractUnitSystem, Map[(StringJoin["Dimension",#] -> DimensionSystem[#]) &,
{"ISQEM", "ISQES", "EMU", "ESU", "Gauss", "Stoney", "Electronic", "PlanckGauss", "Planck", "Natural", "NaturalGauss", "Rydberg", "Hartree", "Hubble", "CosmologicalQuantum"}]]

AppendSystems[f_, l_] := AppendTo[AbstractUnitSystem, Map[(# -> f[#]) &, l]]
AppendSystems[MetricSystem, {"Metric", "SI2019", "MetricEngineering", "SI2019Engineering", "SI1976","CODATA","Conventional"}]
AppendSystems[RankineSystem, {"FPS", "IPS", "British", "English", "Survey"}]
AppendSystems[GaussSystem, {"EMU", "ESU", "Gauss", "LorentzHeaviside", "Kennelly"}]
AppendSystems[EntropySystem, {"International", "InternationalMean", "Nautical", "Meridian", "MeridianEngineering", "MPH", "KKH", "MTS", "IAU", "IAUE", "IAUJ", "Hubble", "Cosmological", "CosmologicalQuantum", "GravitationalMetric", "GravitationalSI2019", "GravitationalMeridian", "FFF"}]

AppendTo[AbstractUnitSystem, {
"Planck" -> MeasureSystem["Planck", 1, 1, 1, 1, PowerExpand[Sqrt[4 Pi AbstractUnitData["\[Alpha]G"]]]],
"PlanckGauss" -> MeasureSystem["PlanckGauss", 1, 1, 1, 4 Pi, PowerExpand[Sqrt[AbstractUnitData["\[Alpha]G"]]]],
"Stoney" -> MeasureSystem["Stoney", 1, 1/"\[Alpha]", 1, 4 Pi, PowerExpand[Sqrt[AbstractUnitData["\[Alpha]G"]/"\[Alpha]"]]],
"Hartree" -> MeasureSystem["Hartree", 1,1,1/"\[Alpha]",4 Pi "\[Alpha]"^2,1],
"Rydberg" -> MeasureSystem["Rydberg", 1,1,2/"\[Alpha]",Pi "\[Alpha]"^2,1/2],
"Schrodinger" -> SchrodingerSystem["Schrodinger", 1, 1, 1/"\[Alpha]", 4 Pi "\[Alpha]"^2, PowerExpand[Sqrt[AbstractUnitData["\[Alpha]G"]/"\[Alpha]"]]],
"Electronic" -> MeasureSystem["Electronic", 1, 1/"\[Alpha]", 1, 4 Pi, 1],
"Natural" -> MeasureSystem["Natural", 1, 1, 1, 1, 1],
"NaturalGauss" -> MeasureSystem["NaturalGauss", 1, 1, 1, 4 Pi, 1],
"QCD" -> MeasureSystem["QCD", 1, 1, 1, 1, 1/AbstractUnitData["\[Mu]pe"]],
"QCDGauss" -> MeasureSystem["QCDGauss", 1, 1, 1, 4 Pi, 1/AbstractUnitData["\[Mu]pe"]],
"QCDoriginal" -> MeasureSystem["QCDoriginal", 1, 1, 1, 4 Pi "\[Alpha]", 1/AbstractUnitData["\[Mu]pe"]],
"SI" :> AbstractUnitSystem["SI2019"],
"MKS" :> AbstractUnitSystem["Metric"],
"SIE" :> AbstractUnitSystem["SI2019Engineering"],
"ME" :> AbstractUnitSystem["MetricEngineering"],
"GSI2019" :> AbstractUnitSystem["GravitationalSI2019"],
"GSI" :> AbstractUnitSystem["GravitationalSI2019"],
"GM" :> AbstractUnitSystem["GravitationalMetric"],
"CGS" :> AbstractUnitSystem["Gauss"],
"CGS2019" :> AbstractUnitSystem["EMU2019"],
"CGSm" :> AbstractUnitSystem["EMU"],
"CGSe" :> AbstractUnitSystem["ESU"],
"HLU" :> AbstractUnitSystem["LorentzHeaviside"],
"EnglishEngineering" :> AbstractUnitSystem["English"],
"BritishGravitational" :> AbstractUnitSystem["British"],
"BG" :> AbstractUnitSystem["British"],
"EnglishUS" :> AbstractUnitSystem["Survey"],
"AbsoluteEnglish" :> AbstractUnitSystem["FPS"],
"AE" :> AbstractUnitSystem["FPS"],
"EE" :> AbstractUnitSystem["English"],
"ISQ" :> AbstractUnitSystem["DimensionISQ"],
"DimensionISQ" :> AbstractUnitSystem["DimensionISQEM"]}]

Map[(UnitSystem[#] := AbstractUnitSystem[#]) &,
{"AbstractUnits","AbstractUnits1","AbstractUnits2","Unified","ISQ","SI","MKS","SIE","ME","GSI2019","GSI","GM","CGS","CGSm","CGSe","HLU","EnglishEngineering","BritishGravitational","BG","EnglishUS","AbsoluteEnglish","AE","EE"}]
Map[(UnitSystem[#] = SWAP[AbstractUnitSystem[#] /. Normal[UnitData],#]) &,
{"Unified","Gauss","LorentzHeaviside","Kennelly","ESU","EMU","Nautical","Meridian","MeridianEngineering","MPH","KKH","MTS","Metric","SI2019","MetricEngineering","SI2019Engineering","GravitationalMetric","GravitationalSI2019","GravitationalMeridian","SI1976","CODATA","Conventional","British","Survey","English","FPS","IPS","IAU","IAUE","IAUJ","FFF","Planck","PlanckGauss","Stoney","Hartree","Rydberg","Schrodinger","Electronic","QCD","QCDGauss","QCDoriginal","International","InternationalMean","Hubble","Cosmological","CosmologicalQuantum","Natural","NaturalGauss"}]
