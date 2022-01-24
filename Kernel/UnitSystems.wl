(* ::Package:: *)
(* This file is part of UnitSystems. It is licensed under the MIT license *)
(* UnitSystems Copyright (C) 2021 Michael Reed *)

Unprotect[UnitSystem];
ProtectedList = {Length, Area, Volume, Power, Entropy};
UnitSystemsList = {"Metric", "SI2019", "CODATA", "Conventional", "MTS", "English",
    "EnglishUS", "IAU", "SI1976", "Mixed", "ESU2019", "EMU2019", "EMU", "ESU", "Gauss",
   "LorentzHeaviside", "Thomson", "Kennelly", "Planck", "PlanckGauss", "Stoney",
   "Hartree", "Rydberg", "Schrodinger", "Electronic", "Natural", "NaturalGauss",
   "QCD", "QCDGauss", "QCDoriginal"};
ConstantsList = {Hyperfine, LightSpeed, Planck, PlanckReduced, 
   ElectronMass, MolarMass, Boltzmann, Permeability, Rationalization, 
   Lorentz, LuminousEfficacy};
PhysicsList = {AtomicMass, ProtonMass, PlanckMass, Newton, Einstein, 
   Hartree, Rydberg, Bohr, BohrReduced, ElectronRadius, Avogadro, 
   Universal, Stefan, RadiationDensity, Permittivity, Coulomb, Ampere,
    BiotSavart, Charge, Faraday, Impedance, Conductance, Klitzing, 
   Josephson, MagneticFlux, Magneton};
KinematicList = {Time, Length, Area, Volume, WaveNumber, 
   FuelEfficiency, Frequency, FrequencyDrift, Speed, Acceleration, 
   Jerk, Snap, VolumeFlow};
MechanicalList = {Mass, MassFlow, LinearDensity, AreaDensity, Density,
    SpecificVolume, Force, Stiffness, Pressure, Compressibility, 
   Viscosity, Diffusivity, RotationalInertia, Momentum, 
   AngularMomentum, Yank, Energy, SpecificEnergy, Action, Fluence, 
   Power, PowerDensity, Intensity, SpectralFlux, SoundExposure, 
   Impedance, SpecificImpedance, Admittance, Compliance, Inertance};
ElectromagneticList = {Charge, ChargeDensity, LinearChargeDensity, 
   Exposure, Mobility, Capacitance, Inductance, Reluctance, Permeance,
    Permittivity, Permeability, Susceptibility, 
   SpecificSusceptibility, DemagnetizingFactor, VectorPotential, 
   ElectricPotential, MagneticPotential, ElectricField, MagneticField,
    ElectricFlux, MagneticFlux, ElectricFluxDensity, 
   MagneticFluxDensity, ElectricDipoleMoment, MagneticDipoleMoment, 
   ElectricPolarizability, MagneticPolarizability, MagneticMoment, 
   Magnetizability, Magnetization, SpecificMagnetization, Rigidity, 
   PoleStrength};
ThermodynamicList = {Temperature, Entropy, SpecificEntropy, 
   VolumeHeatCapacity, ThermalConductivity, ThermalConductance, 
   ThermalResistance, ThermalExpansion, LapseRate};
MolarList = {MolarMass, Molality, Mole, Molarity, MolarVolume, 
   MolarEntropy, MolarEnergy, MolarConductivity, MolarSusceptibility, 
   Catalysis, Specificity};
PhotometricList = {LuminousFlux, Luminance, LuminousEnergy, 
   LuminousExposure, LuminousEfficacy};
MechanicsList = Join[KinematicList, MechanicalList];
ConvertList = 
  Join[MechanicsList, ElectromagneticList, ThermodynamicList, 
   MolarList, PhotometricList];

OneQ[1] := True
OneQ[-1] := False
OneQ[0] := False
OneQ[_] := False
OneQ[x_?NumberQ] := x == 1
OneQ[True] := True
OneQ[False] := False

measure[x_] := x;
GravityCoupling[Coupling[\[Alpha]G_, ___]] := measure[\[Alpha]G];
FineStructure[Coupling[_, \[Alpha]_, ___]] := measure[\[Alpha]];
ElectronUnit[Coupling[_, _, \[Mu]eu_, ___]] := measure[\[Mu]eu];
ProtonUnit[Coupling[_, _, _, \[Mu]pu_, ___]] := measure[\[Mu]pu];
ProtonElectron[c_Coupling] := ProtonUnit[c]/ElectronUnit[c];

(*UnitSystem[kB_, hbar_, c_, mu0_, me_, lambda_] := 
  UnitSystem[kB, hbar, c, mu0, me, lambda, 1];
UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_] := 
  UnitSystem[kB, \[HBar], c, \[Mu]0, me, 1];*)
Boltzmann[UnitSystem[kB_, ___]] := kB;
PlanckReduced[UnitSystem[_, \[HBar]_, ___]] := \[HBar];
LightSpeed[UnitSystem[_, _, c_, ___]] := c;
Permeability[UnitSystem[_, _, _, \[Mu]0_, ___]] := \[Mu]0;
ElectronMass[UnitSystem[_, _, _, _, me_, ___]] := me;
Rationalization[UnitSystem[_, _, _, _, _, \[Lambda]_, ___]] := \[Lambda];
Lorentz[UnitSystem[_, _, _, _, _, _, \[Alpha]L_, ___]] := \[Alpha]L;

Lorentz[_UnitSystem] := 1
Rationalization[_UnitSystem] := 1
RationalizedQ[u_UnitSystem] := Rationalization[u] != 4 Pi

(*Universe[_] = StandardModel*)
Universe[UnitSystem[_?NumericQ, ___]] := StandardModel
Universe[UnitSystem[_Around, ___]] := StandardModel
Unit[x_, y_ : 1] := x;
Mass[u_UnitSystem, s_UnitSystem] := ElectronMass[u, s];
ElectronMass[h_?NumberQ] := UnitData["R\[Infinity]"] 2 h/UnitData["\[Alpha]"]^2/UnitData["c"];
ElectronMass[h_?NumberQ, u_Coupling] := (UnitData["R\[Infinity]"] 2 h)/(FineStructure[u]^2 UnitData["c"]);
PlanckMass[u_UnitSystem, c_Coupling] := ElectronMass[u, c]/Sqrt[GravityCoupling[c]];
Planck[u_UnitSystem, c_Coupling] := 2 Pi PlanckReduced[u];
Newton[u_UnitSystem, c_Coupling] := (
  LightSpeed[u, c] PlanckReduced[u, c])/PlanckMass[u, c]^2;
Charge[u_UnitSystem, c_Coupling] := Sqrt[2 Planck[u]/(Permeability[u]/FineStructure[u])/
  (LightSpeed[u] Rationalization[u] Lorentz[u]^2)];

Map[(#[u_UnitSystem] := #[u, Universe[u]]) &, {PlanckMass, Planck, Newton, Charge}];
Map[(#[u_UnitSystem] := #[Universe[u]]) &, {GravityCoupling, 
   FineStructure, ElectronUnit, ProtonUnit, ProtonElectron}];
Map[(#[u_UnitSystem, c_Coupling] := #[u]) &, {Boltzmann, 
   PlanckReduced, LightSpeed, Permeability, ElectronMass, MolarMass}];
Map[(#[u_UnitSystem, s_UnitSystem] := Unit[#[s]/#[u]]) &, ConstantsList];
Map[(If[! MemberQ[ProtectedList, #],
     #[v_?NumberQ, u_UnitSystem] := #[v, u, Metric];
     #[v_?NumberQ, u_UnitSystem, s_UnitSystem] := 
      Module[{n = #[u, s]}, If[OneQ[n], v, v/n]];
     #[v_?NumberQ, 
       u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___], 
       s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := v;
     #[u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___], 
       s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := 1;
     If[!MemberQ[Join[ProtectedList, ConstantsList, {Permittivity, Charge, MagneticFlux, 
       Impedance, Conductance}], #], #[u_UnitSystem] := #[u, Metric], Nothing];,
     UnitSystem /: #[v_?NumberQ, u_UnitSystem] := #[v, u, Metric];
     UnitSystem /: #[v_?NumberQ, u_UnitSystem, s_UnitSystem] := 
      Module[{n = #[u, s]}, If[OneQ[n], v, v/n]];
     UnitSystem /: #[v_?NumberQ, 
       u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___], 
       s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := v;
     UnitSystem /: #[
       u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___], 
       s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := 1;
     UnitSystem /: #[u_UnitSystem] := #[u, Metric];]
    ) &, ConvertList];

UnitData = <||>

AppendTo[UnitData, "g0" -> 9.80665]
AppendTo[UnitData, "ft" -> 3048/10000]
AppendTo[UnitData, "ftUS" -> 1200/3937]
AppendTo[UnitData, "lb" -> 0.45359237]
AppendTo[UnitData, "rankine" -> 5/9]
(*AppendTo[UnitData, "kelvin" -> 9/5]
AppendTo[UnitData, "atm" -> 101325]

AppendTo[UnitData, "kcalth" -> 4184]
AppendTo[UnitData, "kcal4" -> 4204]
AppendTo[UnitData, "kcal10" -> 4185+1/2]
AppendTo[UnitData, "kcal20" -> 4182]
AppendTo[UnitData, "kcalm" -> 4190]
AppendTo[UnitData, "kcalit" -> 4186+8/10]
AppendTo[UnitData, "cal4" -> UnitData[["kcal4"]]/1000]
AppendTo[UnitData, "cal10" -> UnitData[["kcal10"]]/1000]
AppendTo[UnitData, "cal20" -> UnitData[["kcal20"]]/1000]
AppendTo[UnitData, "calm" -> UnitData[["kcalm"]]/1000]
AppendTo[UnitData, "calit" -> UnitData[["kcalit"]]/1000]
AppendTo[UnitData, "calth" -> UnitData[["kcalth"]]/1000]
AppendTo[UnitData, "kcal" -> UnitData[["kcalth"]]]
AppendTo[UnitData, "cal" -> UnitData[["kcal"]]/1000]*)
(*AppendTo[UnitData, "calth" -> thermal calorie*)

AppendTo[UnitData, "\[CapitalDelta]\[Nu]Cs" -> 9192631770.0]
AppendTo[UnitData, "Kcd" -> 683 555.016/555]
AppendTo[UnitData, "mP" -> Around[2.176434 10^-8, 2.4 10^-13]]
AppendTo[UnitData, "NA" -> 6.02214076 10^23]
AppendTo[UnitData, "kB" -> 1.380649 10^-23]
AppendTo[UnitData, "h" -> 6.62607015 10^-34]
AppendTo[UnitData, "c" -> 299792458.]
AppendTo[UnitData, "e" -> 1.602176634 10^-19]
AppendTo[UnitData, "\[Mu]eu" -> 1/Around[1822.888486209, 5.3 10^-8]]
AppendTo[UnitData, "\[Mu]pu" -> Around[1.007276466621, 5.3 10^-11]]
AppendTo[UnitData, "\[Alpha]" -> 1/Around[137.035999084, 2.1 10^-8]]
AppendTo[UnitData, "R\[Infinity]" -> Around[10973731.5681601, 2.1 10^-5]]
AppendTo[UnitData, "RK1990" -> 25812.807]
AppendTo[UnitData, "RK2014" -> Around[25812.8074555, 5.9 10^-6]]
AppendTo[UnitData, "KJ1990" -> 4.835979 10^14]
AppendTo[UnitData, "KJ2014" -> Around[4.835978525 10^14, 3 10^6]]
AppendTo[UnitData, "GMsun" -> Around[1.32712442099 10^20, 9 10^9]]
AppendTo[UnitData, "day" -> 60^2 24]
AppendTo[UnitData, "au" -> 149597870.7 10^3]

AbstractUnitData["lbm"] = "g0"/"ft"
AbstractUnitData["lbmUS"] = "g0"/"ftUS"
AbstractUnitData["slug"] = "lb" AbstractUnitData["lbm"]
AbstractUnitData["slugUS"] = "lb" AbstractUnitData["lbmUS"]
AbstractUnitData["me"] = 2 "R\[Infinity]" "h"/"\[Alpha]"^2/"c"
AbstractUnitData["\[Mu]0"] = 2 "\[Alpha]" "h"/"c"/"e"^2
AbstractUnitData["\[HBar]"] = "h"/(2 Pi)
AbstractUnitData["\[Delta]\[Mu]0"] = "\[Mu]0" - 4 Pi 10^-7
AbstractUnitData["\[Mu]pe"] = "\[Mu]pu"/"\[Mu]eu"
AbstractUnitData["Ru"] = "NA" "kB"
AbstractUnitData["\[Alpha]L"] = 1/100/"c"
AbstractUnitData["\[Alpha]G"] = (AbstractUnitData["me"]/"mP")^2
AbstractUnitData["\[HBar]1990"] = 2/"RK1990"/"KJ1990"^2/Pi
AbstractUnitData["\[HBar]2014"] = 2/"RK2014"/"KJ2014"^2/Pi
AbstractUnitData["me1990"] = 4 Pi "R\[Infinity]" AbstractUnitData["\[HBar]1990"]/"\[Alpha]"^2/"c"
AbstractUnitData["me2014"] = 4 Pi "R\[Infinity]" AbstractUnitData["\[HBar]2014"]/"\[Alpha]"^2/"c"
AbstractUnitData["GG"] = "c" AbstractUnitData["\[HBar]"]/"mP"^2
AbstractUnitData["ms"] = "GMsun"/AbstractUnitData["GG"]
AbstractUnitData["Js"] = AbstractUnitData["ms"] "au"^2/"day"^2;
AbstractUnitData["mf"] = 90/AbstractUnitData["lbm"]/Mass[AbstractUnitSystem["Metric"], AbstractUnitSystem["English"]]
AbstractUnitData["Jf"] = AbstractUnitData["mf"] (201.168/(14 "day"))^2

EvalUnitData[x] := AbstractUnitData[x]/.Normal[UnitData]

(*AppendTo[UnitData, "lbm" -> EvalUnitData["lbm"]]
AppendTo[UnitData, "lbmUS" -> EvalUnitData["lbmUS"]]
AppendTo[UnitData, "slug" -> EvalUnitData["slug"]]
AppendTo[UnitData, "slugUS" -> EvalUnitData["slugUS"]]
AppendTo[UnitData, "me" -> EvalUnitData["me"]]
AppendTo[UnitData, "\[Mu]0" -> EvalUnitData["\[Mu]0"]] (*\[TildeTilde]4\[Pi]*(1e-7+5.5e-17),exact charge*)
AppendTo[UnitData, "\[HBar]" -> EvalUnitData["\[HBar]"]]
AppendTo[UnitData, "\[Delta]\[Mu]0" -> EvalUnitData["\[Delta]\[Mu]0"]]
AppendTo[UnitData, "\[Mu]pe" -> EvalUnitData["\[Mu]pe"]]
AppendTo[UnitData, "Ru" -> EvalUnitData["Ru"]]
AppendTo[UnitData, "\[Alpha]L" -> EvalUnitData["\[Alpha]L"]]
AppendTo[UnitData, "\[Alpha]G" -> EvalUnitData["\[Alpha]G"]]
AppendTo[UnitData, "\[HBar]1990" -> EvalUnitData["\[HBar]1990"]]
AppendTo[UnitData, "\[HBar]2014" -> EvalUnitData["\[HBar]2014"]]
AppendTo[UnitData, "me1990" -> EvalUnitData["me1990"]]
AppendTo[UnitData, "me2014" -> EvalUnitData["me2014"]]*)
(*AppendTo[UnitData, "GMearth" -> Around[398600441.8 10^6, 8 10^5]]
AppendTo[UnitData, "GMjupiter" -> Around[1.26686534 10^17, 9 10^9]*)
(*AppendTo[UnitData, "LD" -> 384402 10^3]*)
(*AppendTo[UnitData, "pc" -> ("au" 648000/Pi)/.Normal[UnitData]]
AppendTo[UnitData, "ly" -> (365.25 "c" "day")/.Normal[UnitData]]
AppendTo[UnitData, "GG" -> EvalUnitData["GG"]]
AppendTo[UnitData, "ms" -> EvalUnitData["ms"]]
AppendTo[UnitData, "Js" -> EvalUnitData["Js"]]*)

Coupling["AbstractCoupling"] = Coupling["\[Alpha]G", "\[Alpha]", "\[Mu]eu", "\[Mu]pu"]
Coupling["AbstractUniverse"] = Coupling[AbstractUnitData["\[Alpha]G"], "\[Alpha]", "\[Mu]eu", "\[Mu]pu"]
Coupling["StandardModel"] = Coupling["AbstractUniverse"] /. Normal[UnitData]

AbstractUnitSystem["AbstractUnits"] = UnitSystem["kB", "\[HBar]", "c", "\[Mu]0", "me", "\[Lambda]", "\[Alpha]L"]
AbstractUnitSystem["AbstractUnits1"] = UnitSystem["kB1", "\[HBar]1", "c1", "\[Mu]01", "me1", "\[Lambda]1", "\[Alpha]L1"]
AbstractUnitSystem["AbstractUnits2"] = UnitSystem["kB2", "\[HBar]2", "c2", "\[Mu]02", "me2", "\[Lambda]2", "\[Alpha]L2"]

AbstractUnitSystem["Gauss"] = UnitSystem[10^10 AbstractUnitData["Ru"] AbstractUnitData["me"]/"\[Mu]eu", 10^7 AbstractUnitData["\[HBar]"], 100 "c", 1, 1000 AbstractUnitData["me"], 4 Pi, AbstractUnitData["\[Alpha]L"]]
AbstractUnitSystem["LorentzHeaviside"] = UnitSystem[10^10 AbstractUnitData["Ru"] AbstractUnitData["me"]/"\[Mu]eu", 10^7 AbstractUnitData["\[HBar]"], 100 "c", 1, 1000 AbstractUnitData["me"], 1, AbstractUnitData["\[Alpha]L"]]
AbstractUnitSystem["Thomson"] = UnitSystem[10^10 AbstractUnitData["Ru"] AbstractUnitData["me"]/"\[Mu]eu", 10^7 AbstractUnitData["\[HBar]"], 100 "c", 1, 1000 AbstractUnitData["me"], 4 Pi, 1/2];
AbstractUnitSystem["Kennelly"] = UnitSystem[1000 AbstractUnitData["Ru"] AbstractUnitData["me"]/"\[Mu]eu", AbstractUnitData["\[HBar]"], "c", 10^-7, AbstractUnitData["me"], 4 Pi]
AbstractUnitSystem["ESU"] = UnitSystem[10^10 AbstractUnitData["Ru"] AbstractUnitData["me"]/"\[Mu]eu", 10^7 AbstractUnitData["\[HBar]"], 100 "c", (100 "c")^-2, 1000 AbstractUnitData["me"], 4 Pi];
AbstractUnitSystem["ESU2019"] = UnitSystem[10^7 "kB", 10^7 AbstractUnitData["\[HBar]"], 100 "c", 10^3 AbstractUnitData["\[Mu]0"]/"c"^2, 1000 AbstractUnitData["me"]]
AbstractUnitSystem["EMU"] = UnitSystem[10^10 AbstractUnitData["Ru"] AbstractUnitData["me"]/"\[Mu]eu", 10^7 AbstractUnitData["\[HBar]"], 100 "c", 1, 1000 AbstractUnitData["me"], 4 Pi];
AbstractUnitSystem["EMU2019"] = UnitSystem[10^7 "kB", 10^7 AbstractUnitData["\[HBar]"], 100 "c", 10^7 AbstractUnitData["\[Mu]0"], 1000 AbstractUnitData["me"]]
AbstractUnitSystem["MTS"] = UnitSystem[10^6 AbstractUnitData["Ru"] AbstractUnitData["me"]/"\[Mu]eu", 1000 AbstractUnitData["\[HBar]"], "c", 4 Pi/10^4, AbstractUnitData["me"]/1000];
AbstractUnitSystem["Mixed"] = UnitSystem[1000 AbstractUnitData["Ru"] AbstractUnitData["me"]/"\[Mu]eu", AbstractUnitData["\[HBar]"], "c", AbstractUnitData["\[Mu]0"], AbstractUnitData["me"]];
AbstractUnitSystem["Metric"] = UnitSystem[1000 AbstractUnitData["Ru"] AbstractUnitData["me"]/"\[Mu]eu", AbstractUnitData["\[HBar]"], "c", 4 Pi 10^-7, AbstractUnitData["me"]];
AbstractUnitSystem["SI1976"] = UnitSystem[8314.32 AbstractUnitData["me"]/"\[Mu]eu", AbstractUnitData["\[HBar]"], "c", 4 Pi 10^-7, AbstractUnitData["me"]];
AbstractUnitSystem["SI2019"] = UnitSystem["kB", AbstractUnitData["\[HBar]"], "c", AbstractUnitData["\[Mu]0"], AbstractUnitData["me"]];

AbstractUnitSystem["CODATA"] = UnitSystem[1000 AbstractUnitData["Ru"] AbstractUnitData["me2014"]/"\[Mu]eu", AbstractUnitData["\[HBar]2014"], "c", 2 "\[Alpha]" AbstractUnitData["RK2014"]/"c", AbstractUnitData["me2014"]];
AbstractUnitSystem["Conventional"] = UnitSystem[1000 AbstractUnitData["Ru"] AbstractUnitData["me1990"]/"\[Mu]eu", AbstractUnitData["\[HBar]1990"], "c", 2 "\[Alpha]" AbstractUnitData["RK1990"]/"c", AbstractUnitData["me1990"]];
AbstractUnitSystem["English"] = UnitSystem["kB" "rankine"/AbstractUnitData["slug"]/"ft"^2, AbstractUnitData["\[HBar]"]/AbstractUnitData["slug"]/"ft"^2, "c"/"ft", 4 Pi, AbstractUnitData["me"]/AbstractUnitData["slug"]];
AbstractUnitSystem["EnglishUS"] = UnitSystem[(1000 AbstractUnitData["Ru"] AbstractUnitData["me"]/"\[Mu]eu" ) ("rankine"/AbstractUnitData["slugUS"]/"ftUS"^2),
	AbstractUnitData["\[HBar]"]/AbstractUnitData["slugUS"]/"ftUS"^2, "c"/"ftUS", 4 Pi, AbstractUnitData["me"]/AbstractUnitData["slugUS"]]
AbstractUnitSystem["IAU"] = UnitSystem[1000 AbstractUnitData["Ru"] AbstractUnitData["me"]/"\[Mu]eu"/AbstractUnitData["Js"], AbstractUnitData["\[HBar]"]/"day"/AbstractUnitData["Js"], "day" "c"/"au",
	4 Pi 10^-7 "day"^2/AbstractUnitData["Js"], AbstractUnitData["me"]/AbstractUnitData["ms"]]
AbstractUnitSystem["FFF"] = UnitSystem[1000 AbstractUnitData["Ru"] AbstractUnitData["me"]/"\[Mu]eu" "rankine"/AbstractUnitData["Jf"], AbstractUnitData["\[HBar]"]/(14 "day")/AbstractUnitData["Jf"], 14 "day" "c"/201.168, 0, AbstractUnitData["me"]/AbstractUnitData["mf"]]

AbstractUnitSystem["Planck"] = UnitSystem[1,1,1,1,Sqrt[4 Pi AbstractUnitData["\[Alpha]G"]],1,1]
AbstractUnitSystem["PlanckGauss"] = UnitSystem[1, 1, 1, 4 Pi, Sqrt[AbstractUnitData["\[Alpha]G"]]]
AbstractUnitSystem["Stoney"] = UnitSystem[1, 1/"\[Alpha]", 1, 4 Pi, Sqrt[AbstractUnitData["\[Alpha]G"]/"\[Alpha]"]]
AbstractUnitSystem["Hartree"] = UnitSystem[1,1,1/"\[Alpha]",4 Pi "\[Alpha]"^2,1]
AbstractUnitSystem["Rydberg"] = UnitSystem[1,1,2/"\[Alpha]",Pi "\[Alpha]"^2,1/2]
AbstractUnitSystem["Schrodinger"] = UnitSystem[1, 1, 1/"\[Alpha]", 4 Pi "\[Alpha]"^2, Sqrt[AbstractUnitData["\[Alpha]G"]/"\[Alpha]"]]
AbstractUnitSystem["Electronic"] = UnitSystem[1, 1/"\[Alpha]", 1, 4 Pi, 1]
AbstractUnitSystem["Natural"] = UnitSystem[1, 1, 1, 1, 1, 1, 1, "1"]
AbstractUnitSystem["NaturalGauss"] = UnitSystem[1, 1, 1, 4 Pi, 1, 1, 1, "1"]
AbstractUnitSystem["QCD"] = UnitSystem[1, 1, 1, 1, 1/AbstractUnitData["\[Mu]pe"]]
AbstractUnitSystem["QCDGauss"] = UnitSystem[1, 1, 1, 4 Pi, 1/AbstractUnitData["\[Mu]pe"]]
AbstractUnitSystem["QCDoriginal"] = UnitSystem[1, 1, 1, 4 Pi "\[Alpha]", 1/AbstractUnitData["\[Mu]pe"]]

UnitSystem["Natural"] = UnitSystem[1, 1, 1, 1, 1]
UnitSystem["NaturalGauss"] = UnitSystem[1, 1, 1, 4 Pi, 1]
Map[(UnitSystem[#] = AbstractUnitSystem[#] /. Normal[UnitData]) &,
{"AbstractUnits","AbstractUnits1","AbstractUnits2","Gauss","LorentzHeaviside","Thomson","Kennelly","ESU","ESU2019","EMU","EMU2019","MTS","Mixed","Metric","SI1976","SI2019","CODATA","Conventional","English","EnglishUS","IAU","FFF","Planck","PlanckGauss","Stoney","Hartree","Rydberg","Schrodinger","Electronic","QCD","QCDGauss","QCDoriginal"}]

(*Map[Set[Unevaluated[#],UnitData[ToString[#]]] &, {g0,ft,ftUS,lb,rankine,\[CapitalDelta]\[Nu]Cs,Kcd,mP,NA,kB,h,c,e,\[Mu]eu,\[Mu]pu,\[Alpha]inv,R\[Infinity],RK1990,RK2014,KJ1990,KJ2014}]*)
(*Map[Set[Unevaluated[#],UnitData[ToString[#]]] &, {lbm,lbmUS,slug,slugUS,kelvin,atm,kcalth,kcal4,kcla10,kcal20,kcalm,kcalit,calth,cal4,cal10,cal20,calm,calit,kcal,cal,me,\[Mu]0,\[HBar],\[Delta]\[Mu]0,\[Mu]pe,Ru,\[Alpha]L,\[Alpha]G,\[HBar]1990,\[HBar]2014,me1990,me2014}]*)
(*Map[Set[Unevaluated[#],UnitData[ToString[#]]] &, {GMsun, GMearth, GMjupiter, au, LD, day, pc, ly, GG, ms, Js, mf, Jf}]*)

StandardModel = Coupling["StandardModel"]
AbstractCoupling = Coupling["AbstractCoupling"]
AbstractUniverse = Coupling["AbstractUniverse"]

AbstractUnits = AbstractUnitSystem["AbstractUnits"]
AbstractUnits1 = AbstractUnitSystem["AbstractUnits1"]
AbstractUnits2 = AbstractUnitSystem["AbstractUnits2"]
AbstractMetric = AbstractUnitSystem["Metric"]
AbstractSI2019 = AbstractUnitSystem["SI2019"]
AbstractCGS = AbstractUnitSystem["Gauss"]
AbstractMTS = AbstractUnitSystem["MTS"]
AbstractEnglish = AbstractUnitSystem["English"]
AbstractEnglishUS = AbstractUnitSystem["EnglishUS"]
AbstractIAU = AbstractUnitSystem["IAU"]
AbstractFFF = AbstractUnitSystem["FFF"]

Gauss = UnitSystem["Gauss"]
LorentzHeaviside = UnitSystem["LorentzHeaviside"]
Thomson = UnitSystem["Thomson"]
Kennelly = UnitSystem["Kennelly"]
ESU = UnitSystem["ESU"]
ESU2019 = UnitSystem["ESU2019"]
EMU = UnitSystem["EMU"]
EMU2019 = UnitSystem["EMU2019"]
MTS = UnitSystem["MTS"]
Mixed = UnitSystem["Mixed"]
Metric = UnitSystem["Metric"]
SI1976 = UnitSystem["SI1976"]
SI2019 = UnitSystem["SI2019"]

CODATA = UnitSystem["CODATA"]
Conventional = UnitSystem["Conventional"]
English = UnitSystem["English"]
EnglishUS = UnitSystem["EnglishUS"]
IAU = UnitSystem["IAU"]
FFF = UnitSystem["FFF"]

PlanckGauss = UnitSystem["PlanckGauss"]
Stoney = UnitSystem["Stoney"]
Schrodinger = UnitSystem["Schrodinger"]
Electronic = UnitSystem["Electronic"]
Natural = UnitSystem["Natural"]
QCD = UnitSystem["QCD"]
QCDGauss = UnitSystem["QCDGauss"]
QCDoriginal = UnitSystem["QCDoriginal"]

{SI,MKS,CGS, CGS2019, CGSm, CGSe, HLU} = {SI2019, Metric, Gauss, EMU2019, EMU, ESU, LorentzHeaviside}
AbstractUnitSystem["SI"] := AbstractUnitSystem["SI2019"]
AbstractUnitSystem["MKS"] := AbstractUnitSystem["Metric"]
AbstractUnitSystem["CGS"] := AbstractUnitSystem["Gauss"]
AbstractUnitSystem["CGS2019"] := AbstractUnitSystem["EMU2019"]
AbstractUnitSystem["CGSm"] := AbstractUnitSystem["EMU"]
AbstractUnitSystem["CGSe"] := AbstractUnitSystem["ESU"]
AbstractUnitSystem["HLU"] := AbstractUnitSystem["LorentzHeaviside"]
UnitSystem["SI"] := UnitSystem["SI2019"]
UnitSystem["MKS"] := UnitSystem["Metric"]
UnitSystem["CGS"] := UnitSystem["Gauss"]
UnitSystem["CGS2019"] := UnitSystem["EMU2019"]
UnitSystem["CGSm"] := UnitSystem["EMU"]
UnitSystem["CGSe"] := UnitSystem["ESU"]
UnitSystem["HLU"] := UnitSystem["LorentzHeaviside"]

ElectronMass[UnitSystem["Planck"], c_Coupling] := Sqrt[4 Pi GravityCoupling[c]];
ElectronMass[UnitSystem["PlanckGauss"], c_Coupling] := Sqrt[GravityCoupling[c]];
ElectronMass[UnitSystem[_, _, _, _, Sqrt[EvalUnitData["\[Alpha]G"]/UnitData["\[Alpha]"]], ___], c_Coupling] := Sqrt[GravityCoupling[c]/FineStructure[c];]
ElectronMass[UnitSystem[_, _, _, _, 1/UnitData["\[Mu]pe"], ___], c_Coupling] := 1/ProtonElectron[c];
Permeability[UnitSystem[_, _, _, 4 Pi UnitData["\[Alpha]"]^2, ___], c_Coupling] := 4 Pi FineStructure[c]^2;
Permeability[UnitSystem[_, _, _, Pi UnitData["\[Alpha]"]^2, ___], c_Coupling] := Pi FineStructure[c]^2;
LightSpeed[UnitSystem[_, _, 1/UnitData["\[Alpha]"], ___], c_Coupling] := 1/FineStructure[c];
LightSpeed[UnitSystem[_, _, 2/UnitData["\[Alpha]"], ___], c_Coupling] := 2/FineStructure[c];
PlanckReduced[UnitSystem[_, 1/UnitData["\[Alpha]"], ___], c_Coupling] := 1/FineStructure[c];

ElectronMass[u : UnitSystem[_, _, UnitData["c"], _, EvalUnitData["me"], ___], c_Coupling] := ElectronMass[Planck[u], c];
ElectronMass[UnitSystem[_, _, 100 UnitData["c"], _, 1000 EvalUnitData["me"], ___], c_Coupling] := 1000 ElectronMass[UnitSystem["SI2019"], c];
ElectronMass[UnitSystem[_, _, UnitData["c"], _, EvalUnitData["me"]/1000, ___], c_Coupling] := ElectronMass[UnitSystem["SI2019"], c]/1000;
ElectronMass[u : UnitSystem[_, _, UnitData["c"], _, EvalUnitData["me2014"], ___], c_Coupling] := ElectronMass[Planck[u], c];
ElectronMass[u : UnitSystem[_, _, UnitData["c"], EvalUnitData["\[Mu]0"], EvalUnitData["me1990"], ___], c_Coupling] := ElectronMass[Planck[u], c];
ElectronMass[UnitSystem[_, _, UnitData["c"]/UnitData["ftUS"], _, EvalUnitData["me"]/EvalUnitData["slug"], ___], c_Coupling] := ElectronMass[UnitSystem["SI2019"], c]/EvalUnitData["slug"];
Permeability[UnitSystem[_, _, _, EvalUnitData["\[Mu]0"], ___], c_Coupling] := FineStructure[c] 2 UnitData["h"]/UnitData["c"]/UnitData["e"]^2;
Permeability[UnitSystem["ESU2019"], u_Coupling] := 10^3 Permeability[UnitSystem["SI2019"], u]/UnitData["c"]^2;
Permeability[UnitSystem["EMU2019"], u_Coupling] := 10^7 Permeability[UnitSystem["SI2019"], u];
Permeability[UnitSystem["CODATA"], u_Coupling] := 2 UnitData["RK2014"] FineStructure[u]/UnitData["c"];
Permeability[UnitSystem["Conventional"], u_Coupling] := 2 UnitData["RK1990"] FineStructure[u]/UnitData["c"];

MolarMass[UnitSystem[1, ___]] = 1;
MolarMass[u:UnitSystem[UnitData["kB"], ___]] := MolarMass[u, Universe[u]];
MolarMass[u:UnitSystem[UnitData["kB"], ___], c_Coupling] := UnitData["NA"] ElectronMass[u, c]/ElectronUnit[c];
MolarMass[u:UnitSystem[10^7 UnitData["kB"], ___]] := MolarMass[u, Universe[u]];
MolarMass[UnitSystem[10^7 UnitData["kB"], ___], c_Coupling] := 1000 MolarMass[UnitSystem["SI2019"], c];
MolarMass[u:UnitSystem[10^3 UnitData["kB"], ___]] := MolarMass[u, Universe[u]];
MolarMass[UnitSystem[10^3 UnitData["kB"], ___], c_Coupling] := MolarMass[UnitSystem["SI2019"], c]/1000;
MolarMass[UnitSystem[kB_, ___]] := MolarMass[UnitSystem["CGS"]]/1000;
MolarMass[UnitSystem[Boltzmann[UnitSystem["MTS"]], ___]] := MolarMass[UnitSystem["CGS"]]/10^6;
MolarMass[UnitSystem[Boltzmann[UnitSystem["CGS"]], ___]] := MolarMass[UnitSystem["Natural"]];
MolarMass[UnitSystem[Boltzmann[UnitSystem["FFF"]], ___]] := MolarMass[UnitSystem["Natural"]];
MolarMass[u:UnitSystem[Boltzmann[UnitSystem["English"]], ___]] := MolarMass[u, Universe[u]];
MolarMass[u:UnitSystem[Boltzmann[UnitSystem["English"]], ___], c_Coupling] := 1000 MolarMass[UnitSystem["SI2019"],c]
MolarMass[UnitSystem[Boltzmann[UnitSystem["EnglishUS"]], ___]] := MolarMass[UnitSystem["Natural"]];
MolarMass[UnitSystem[Boltzmann[UnitSystem["IAU"]], ___]] = EvalUnitData["ms"]/1000;

LuminousEfficacy[UnitSystem[1, ___]] = 1
LuminousEfficacy[u : UnitSystem[_?NumericQ, ___]] := Power[UnitData["Kcd"], UnitSystem["SI2019"], u]
LuminousEfficacy[u : UnitSystem[_Around, ___]] := Power[UnitData["Kcd"], UnitSystem["SI2019"], u]

Universe[_] = Coupling["AbstractUniverse"]

MolarMass[AbstractUnitSystem["AbstractUnits"]] = "Mu"
MolarMass[AbstractUnitSystem["AbstractUnits1"]] = "Mu1"
MolarMass[AbstractUnitSytem["AbstractUnits2"]] = "Mu2"
MolarMass[u : UnitSystem["kB", ___]] := MolarMass[u,Universe[u]]
MolarMass[u : UnitSystem["kB", ___], c_Coupling] := "NA" ElectronMass[AbstractUnitSystem["SI2019"],c]/ElectronUnit[c]
MolarMass[u : UnitSystem[10^7 "kB", ___]] := MolarMass[u, Universe[u]];
MolarMass[UnitSystem[10^7 "kB", ___], c_Coupling] := 1000 MolarMass[AbstractUnitSystem["SI2019"], c];
MolarMass[u : UnitSystem[10^3 "kB", ___]] := MolarMass[u, Universe[u]];
MolarMass[UnitSystem[10^3 "kB", ___], c_Coupling] := MolarMass[AbstractUnitSystem["SI2019"], c]/1000;
MolarMass[UnitSystem[Boltzmann[AbstractUnitSystem["MTS"]], ___]] := MolarMass[UnitSystem["CGS"]]/10^6;
MolarMass[UnitSystem[Boltzmann[AbstractUnitSystem["CGS"]], ___]] := MolarMass[UnitSystem["Natural"]];
MolarMass[UnitSystem[Boltzmann[AbstractUnitSystem["FFF"]], ___]] := MolarMass[UnitSystem["Natural"]];
MolarMass[u : UnitSystem[Boltzmann[AbstractUnitSystem["English"]], ___]] := MolarMass[u, Universe[u]];
MolarMass[u : UnitSystem[Boltzmann[AbstractUnitSystem["English"]], ___], c_Coupling] := 
  1000 MolarMass[AbstractUnitSystem["SI2019"], c];
MolarMass[UnitSystem[Boltzmann[AbstractUnitSystem["EnglishUS"]], ___]] := MolarMass[UnitSystem["Natural"]];
MolarMass[UnitSystem[Boltzmann[AbstractUnitSystem["IAU"]], ___]] := 1/1000 AbstractUnitData["ms"];

LuminousEfficacy[AbstractUnitSystem["AbstractUnits1"]] = "Kcd1"
LuminousEfficacy[AbstractUnitSystem["AbstractUnits2"]] = "Kcd2"
LuminousEfficacy[u_UnitSystem] := Power["Kcd", AbstractUnitSystem["SI2019"], u]

AbstractUnitData["Mu"] = MolarMass[AbstractUnitSystem["SI2019"]]

Kilograms[m_] := Kilograms[m, UnitSystem["English"]];
Kilograms[m_, u_UnitSystem] := Mass[m, UnitSystem["Metric"], u];
Slugs[m_] := Slugs[m, UnitSystem["Metric"]];
Slugs[m_, u_UnitSystems] := Mass[m, UnitSystem["English"], u];
Feet[d_] := Feet[d, UnitSystem["Metric"]];
Feet[d_, u_UnitSystem] := Length[d, UnitSystem["English"], u];
Meters[d_] := Meters[d, UnitSystem["English"]];
Meters[d_, u_UnitSystem] := Length[d, UnitSystem["Metric"], u];

(* IAU to SI *)
UnitSystem /: 
  Length[u : UnitSystem[_, _, UnitData["c"], ___],
   s : UnitSystem[_, _, UnitData["day"] UnitData["c"]/UnitData["au"], ___]] := Length[u, s, 1/UnitData["au"]];
Time[u : UnitSystem[_, _, UnitData["c"], ___],
   s : UnitSystem[_, _, UnitData["day"] UnitData["c"]/UnitData["au"], ___]] := Time[u, s, 1/UnitData["day"]];
(* SI to IAU *)
UnitSystem /: 
  Length[u : UnitSystem[_, _, UnitData["day"] UnitData["c"]/UnitData["au"], ___],
   s : UnitSystem[_, _, UnitData["c"], ___]] := Length[u, s, UnitData["au"]];
Time[u : UnitSystem[_, _, UnitData["day"] UnitData["c"]/UnitData["au"], ___],
   s : UnitSystem[_, _, UnitData["c"], ___]] := Time[u, s, UnitData["day"]];
(* IAU to CGS *)
UnitSystem /: 
  Length[u : UnitSystem[_, _, 100 UnitData["c"], ___],
   s : UnitSystem[_, _, UnitData["day"] UnitData["c"]/UnitData["au"], ___]] := Length[u, s, 1/UnitData["au"]];
Time[u : UnitSystem[_, _, 100 UnitData["c"], ___],
   s : UnitSystem[_, _, UnitData["day"] UnitData["c"]/UnitData["au"], ___]] := Time[u, s, 1/UnitData["day"]];
(* CGS to IAU *)
UnitSystem /: 
  Length[u : UnitSystem[_, _, UnitData["day"] UnitData["c"]/UnitData["au"], ___],
   s : UnitSystem[_, _, 100 UnitData["c"], ___]] := Length[u, s, UnitData["au"]];
Time[u : UnitSystem[_, _, UnitData["day"] UnitData["c"]/UnitData["au"], ___],
   s : UnitSystem[_, _, 100 UnitData["c"], ___]] := Time[u, s, UnitData["day"]];
(* IAU to English *)
UnitSystem /: 
  Length[u : UnitSystem[_, _, UnitData["c"]/UnitData["ft"], ___],
   s : UnitSystem[_, _, UnitData["day"] UnitData["c"]/UnitData["au"], ___]] := Length[u, s, UnitData["ft"]/UnitData["au"]];
Time[u : UnitSystem[_, _, UnitData["c"]/UnitData["ft"], ___],
   s : UnitSystem[_, _, UnitData["day"] UnitData["c"]/UnitData["au"], ___]] := Time[u, s, 1/UnitData["day"]];
(* English to IAU *)
UnitSystem /: 
  Length[u : UnitSystem[_, _, UnitData["day"] UnitData["c"]/UnitData["au"], ___],
   s : UnitSystem[_, _, UnitData["c"]/UnitData["ft"], ___]] := Length[u, s, UnitData["au"]/UnitData["ft"]];
Time[u : UnitSystem[_, _, UnitData["day"] UnitData["c"]/UnitData["au"], ___],
   s : UnitSystem[_, _, UnitData["c"]/UnitData["ft"], ___]] := Time[u, s, UnitData["day"]];

UnitSystem /: Length[u_UnitSystem, s_UnitSystem] := Length[u, s, 1];
UnitSystem /: Length[u_UnitSystem, s_UnitSystem, l_] := 
  Unit[(PlanckReduced[s] ElectronMass[u] LightSpeed[
       u])/(PlanckReduced[u] ElectronMass[s] LightSpeed[s]), l];
UnitSystem /: Area[u_UnitSystem, s_UnitSystem] := Unit[Length[u, s]^2];
UnitSystem /: Volume[u_UnitSystem, s_UnitSystem] := Unit[Length[u, s]^3];

WaveNumber[u_UnitSystem, s_UnitSystem] := Unit[Length[s, u]];
FuelEfficiency[u_UnitSystem, s_UnitSystem] := Area[s, u];
Time[u_UnitSystem, s_UnitSystem] := Time[u, s, 1];
Time[u_UnitSystem, s_UnitSystem, t_] := Unit[Length[u, s]/LightSpeed[u, s], 1];
Frequency[u_UnitSystem, s_UnitSystem] := Time[s, u];
FrequencyDrift[u_UnitSystem, s_UnitSystem] := Unit[Time[s, u]^2];
Speed[u_UnitSystem, s_UnitSystem] := LightSpeed[u, s];
Acceleration[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]/Time[u, s]];
Jerk[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]/Time[u, s]^2];
Snap[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]/Time[u, s]^3];
VolumeFlow[u_UnitSystem, s_UnitSystem] := Unit[Area[u, s], Speed[u, s]];
SpecificEnergy[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]^2];

Energy[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] SpecificEnergy[u, s]];
UnitSystem /: Power[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/Time[u, s]];
Force[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] Acceleration[u, s]];
Pressure[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/Length[u, s]/Time[u, s]^2];

Momentum[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] Speed[u, s]];
AngularMomentum[u_UnitSystem, s_UnitSystem] := Unit[Momentum[u, s] Length[u, s]];
Yank[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] Jerk[u, s]];
AreaDensity[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/Area[u, s]];
Density[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/Volume[u, s]];
SpecificVolume[u_UnitSystem, s_UnitSystem] := Unit[Volume[u, s]/Mass[u, s]];
Action[u_UnitSystem, s_UnitSystem] := Unit[Momentum[u, s] Length[u, s]];
Stiffness[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/Area[u, s]];
Intensity[u_UnitSystem, s_UnitSystem] := Unit[Power[u, s]/Area[u, s]];
Diffusivity[u_UnitSystem, s_UnitSystem] := 
  Unit[(PlanckReduced[s] ElectronMass[u])/(PlanckReduced[u] ElectronMass[s])];
Viscosity[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/Length[u, s]/Time[u, s]];
LinearDensity[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/Length[u, s]];
MassFlow[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/Time[u, s]];
SpectralFlux[u_UnitSystem, s_UnitSystem] := Unit[Power[u, s]/Length[u, s]];
PowerDensity[u_UnitSystem, s_UnitSystem] := Unit[Power[u, s]/Volume[u, s]];
Compressibility[u_UnitSystem, s_UnitSystem] := Pressure[s, u];
Fluence[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/Area[u, s]];
RotationalInertia[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] Area[u, s]];

SoundExposure[u_UnitSystem, s_UnitSystem] := Unit[Time[u, s] Pressure[u, s]^2];
SpecificImpedance[u_UnitSystem, s_UnitSystem] := Unit[Pressure[u, s]/Speed[u, s]];
Impedance[u_UnitSystem, s_UnitSystem] := Unit[SpecificImpedance[u, s]/Area[u, s]];
Admittance[u_UnitSystem, s_UnitSystem] := Unit[Area[u, s]/SpecificImpedance[u, s]];
Compliance[u_UnitSystem, s_UnitSystem] := Unit[Time[u, s]^2/Mass[u, s]];
Inertance[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/Length[u, s]^4];

(* Electromagnetic *)

Voltage = ElectricPotential

Charge[u_UnitSystem, s_UnitSystem] := 
  Unit[Sqrt[(PlanckReduced[s] Permeability[u] LightSpeed[
        u] Rationalization[u] Lorentz[u]^2)/(PlanckReduced[
        u] Permeability[s] LightSpeed[s] Rationalization[s] Lorentz[s]^2)]];
Current[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s]/Time[u, s]];
ElectricPotential[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/Charge[u, s]];
Capacitance[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s]/ElectricPotential[u, s]];
Resistance[u_UnitSystem, s_UnitSystem] := Unit[ElectricPotential[u, s]/Current[u, s]];
Conductance[u_UnitSystem, s_UnitSystem] := Unit[Current[u, s]/ElectricPotential[u, s]];
MagneticFlux[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/Lorentz[u, s]/Current[u, s]];
MagneticFluxDensity[u_UnitSystem, s_UnitSystem] := 
  Unit[Mass[u, s]/Lorentz[u, s]/Current[u, s]/Time[u, s]^2];
Inductance[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] Area[u, s]/Charge[u, s]^2];

ElectricFluxDensity[u_UnitSystem, s_UnitSystem] := 
  Unit[Charge[u, s] Rationalization[u, s]/Area[u, s]];
ChargeDensity[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s]/Volume[u, s]];
CurrentDensity[u_UnitSystem, s_UnitSystem] := Unit[Current[u, s]/Area[u, s]];
Conductivity[u_UnitSystem, s_UnitSystem] := Unit[Conductance[u, s]/Length[u, s]];
Permittivity[u_UnitSystem, s_UnitSystem] := 
  Unit[Capacitance[u, s] Rationalization[u, s]/Length[u, s]];
ElectricField[u_UnitSystem, s_UnitSystem] := 
  Unit[ElectricPotential[u, s]/Length[u, s]];
MagneticField[u_UnitSystem, s_UnitSystem] := 
  Unit[Current[u, s] Rationalization[u, s] Lorentz[u, s]/Length[u, s]];
Exposure[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s]/Mass[u, s]];
Resistivity[u_UnitSystem, s_UnitSystem] := Unit[Resistance[u, s] Length[u, s]];
LinearChargeDensity[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s]/Length[u, s]];
MagneticDipoleMoment[u_UnitSystem, s_UnitSystem] := 
  Unit[Current[u, s] Lorentz[u, s] Area[u, s]];
Mobility[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s] Time[u, s]/Mass[u, s]];
Reluctance[u_UnitSystem, s_UnitSystem] := 
  Unit[Rationalization[u, s] Lorentz[u, s]^2/Inductance[u, s]];
VectorPotential[u_UnitSystem, s_UnitSystem] := Unit[MagneticFlux[u, s]/Length[u, s]];
MagneticMoment[u_UnitSystem, s_UnitSystem] := Unit[MagneticFlux[u, s] Length[u, s]];
Rigidity[u_UnitSystem, s_UnitSystem] := Unit[MagneticFluxDensity[u, s] Length[u, s]];
Susceptibility[u_UnitSystem, s_UnitSystem] := Unit[Rationalization[s, u]];

(* WARNING unchecked: rigitidy, magneticmoment, vectorpotential, \
mobility, linearchargedensity, exposure *)

ElectricFlux[u_UnitSystem, s_UnitSystem] := Unit[ElectricPotential[u, s] Length[u, s]];
ElectricDipoleMoment[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s] Length[u, s]];
MagneticPotential[u_UnitSystem, s_UnitSystem] := Unit[MagneticFlux[u, s] Reluctance[u, s]];
PoleStrength[u_UnitSystem, s_UnitSystem] := Unit[MagneticDipoleMoment[u, s]/Length[u, s]];
Permeance[u_UnitSystem, s_UnitSystem] := Reluctance[s, u];
SpecificSusceptibility[u_UnitSystem, s_UnitSystem] := 
  Unit[MagneticDipoleMoment[u, s]/MagneticField[u, s]/Mass[u, s]];
Magnetizability[u_UnitSystem, s_UnitSystem] := 
  Unit[MagneticMoment[u, s]/MagneticFluxDensity[u, s]];
ElectricPolarizability[u_UnitSystem, s_UnitSystem] := 
  Unit[ElectricDipoleMoment[u, s]/ElectricField[u, s]];
MagneticPolarizability[u_UnitSystem, s_UnitSystem] := 
  Unit[MagneticDipoleMoment[u, s]/MagneticField[u, s]];
Magnetization[u_UnitSystem, s_UnitSystem] := Unit[MagneticMoment[u, s]/Volume[u, s]];

SpecificMagnetization[u_UnitSystem, s_UnitSystem] := Unit[MagneticMoment[s, u]/Mass[s, u]];
DemagnetizingFactor[u_UnitSystem, s_UnitSystem] := Unit[Rationalization[u, s]];

(* Thermodynamic *)

Moles[n_] := Moles[n, Metric];
Moles[n_, u_UnitSystem] := n/Avogadro[u];
Molecules[n_] := Molecules[n, Metric];
Molecules[n_, u_UnitSystem] := n Avogadro[u];

Temperature[u_UnitSystem, s_UnitSystem] := 
  Unit[(Boltzmann[u] ElectronMass[s] LightSpeed[s]^2)/(Boltzmann[
       s] ElectronMass[u] LightSpeed[u]^2)];
UnitSystem /: Entropy[u_UnitSystem, s_UnitSystem] := 
  Unit[Energy[u, s]/Temperature[u, s]];
SpecificEntropy[u_UnitSystem, s_UnitSystem] := 
  Unit[SpecificEnergy[u, s]/Temperature[u, s]];
VolumeHeatCapacity[u_UnitSystem, s_UnitSystem] := Unit[Entropy[u, s]/Volume[u, s]];
ThermalConductivity[u_UnitSystem, s_UnitSystem] := 
  Unit[Force[u, s]/Time[u, s]/Temperature[u, s]];
ThermalConductance[u_UnitSystem, s_UnitSystem] := 
  Unit[ThermalConductivity[u, s] Length[u, s]];
ThermalResistance[u_UnitSystem, s_UnitSystem] := ThermalConductance[s, u];
ThermalExpansion[u_UnitSystem, s_UnitSystem] := Temperature[s, u];
LapseRate[u_UnitSystem, s_UnitSystem] := Unit[Temperature[u, s]/Length[u, s]];

Molality[u_UnitSystem, s_UnitSystem] := MolarMass[s, u];
Mole[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] Molality[u, s]];
Molarity[u_UnitSystem, s_UnitSystem] := 
  Unit[Mole[u, s]/Volume[u, s]];
MolarVolume[u_UnitSystem, s_UnitSystem] := 
  Unit[Volume[u, s]/Mole[u, s]];
MolarEntropy[u_UnitSystem, s_UnitSystem] := 
  Unit[Entropy[u, s]/Mole[u, s]];
MolarEnergy[u_UnitSystem, s_UnitSystem] := 
  Unit[Energy[u, s]/Mole[u, s]];
MolarConductivity[u_UnitSystem, s_UnitSystem] := 
  Unit[Conductivity[u, s] Area[u, s]/Mole[u, s]];
MolarSusceptibility[u_UnitSystem, s_UnitSystem] := 
  Unit[SpecificSusceptibility[u, s] MolarMass[u, s]];
Catalysis[u_UnitSystem, s_UnitSystem] := Unit[Mole[u, s]/Time[u, s]];
Specificity[u_UnitSystem, s_UnitSystem] := 
  Unit[Volume[u, s]/Mole[u, s]/Time[u, s]];

LuminousFlux[u_UnitSystem, s_UnitSystem] := 
  Unit[Frequency[u, 
     s]^2 (LuminousEfficacy[s] PlanckReduced[s])/(LuminousEfficacy[
        u] PlanckReduced[u])];
Luminance[u_UnitSystem, s_UnitSystem] := 
  Unit[LuminousFlux[u, s]/Area[u, s]];
LuminousEnergy[u_UnitSystem, s_UnitSystem] := 
  Unit[Frequency[u, 
     s] (LuminousEfficacy[s] PlanckReduced[s])/(LuminousEfficacy[
        u] PlanckReduced[u])];
LuminousExposure[u_UnitSystem, s_UnitSystem] := 
  Unit[Luminance[u, s] Time[u, s]];

(* Physics *)

Hyperfine[u_UnitSystem] := Frequency[UnitData["\[CapitalDelta]\[Nu]Cs"], u];

Map[(#[u_UnitSystem] := #[u, Universe[u]]) &, {Avogadro, AtomicMass, 
   ProtonMass, Einstein, Universal, Stefan, RadiationDensity, 
   Permittivity, Coulomb, BiotSavart, Impedance, Faraday, Josephson, 
   MagneticFlux, Klitzing, Conductance, Hartree, Rydberg, Bohr, 
   BohrReduced, ElectronRadius, Magneton}];
Avogadro[u_UnitSystem, c_Coupling] := 
  MolarMass[u, c] ElectronUnit[c]/ElectronMass[u, c];
AtomicMass[u_UnitSystem, c_Coupling] := ElectronMass[u, c]/ElectronUnit[c];
ProtonMass[u_UnitSystem, c_Coupling] := ProtonElectron[c] ElectronMass[u, c];
Einstein[u_UnitSystem, c_Coupling] := (8 Pi Newton[u, c])/LightSpeed[u, c]^4;
Universal[u_UnitSystem, c_Coupling] := Boltzmann[u, c] Avogadro[u, c];
Stefan[u_UnitSystem, c_Coupling] := (2 Pi^5 Boltzmann[u, c]^4)/(
  15 Planck[u, c]^3 LightSpeed[u, c]^2);
RadiationDensity[u_UnitSystem, c_Coupling] := (4 Stefan[u, c])/LightSpeed[u, c];
Permittivity[u_UnitSystem, c_Coupling] := 
  (1/Permeability[u, c]) (LightSpeed[u, c] Lorentz[u]^2);
Coulomb[u_UnitSystem, c_Coupling] := 
  Rationalization[u]/(4 Pi)/ Permittivity[u];
BiotSavart[u_UnitSystem, c_Coupling] := 
  Permeability[u, c] Lorentz[u] (Rationalization[u]/(4 Pi));
Ampere[u_UnitSystem] := Lorentz[u] BiotSavart[u];
Impedance[u_UnitSystem, c_Coupling] :=
  
Permeability[u, c] LightSpeed[u, c] Rationalization[u] Lorentz[u]^2;
Faraday[u_UnitSystem, c_Coupling] := Charge[u, c] Avogadro[u, c];
Josephson[u_UnitSystem, c_Coupling] := 2 Charge[u, c] Lorentz[u]/Planck[u, c];
MagneticFlux[u_UnitSystem, c_Coupling] := 1/Josephson[u, c];
Klitzing[u_UnitSystem, c_Coupling] := Planck[u, c]/Charge[u, c]^2;
Conductance[u_UnitSystem, c_Coupling] := (2 Charge[u, c]^2)/Planck[u, c];
Hartree[u_UnitSystem, c_Coupling] := 
  ElectronMass[u, c] (LightSpeed[u, c] FineStructure[c])^2;
Rydberg[u_UnitSystem, c_Coupling] := Hartree[u, c]/(2 Planck[u, c])/LightSpeed[u, c];
Bohr[u_UnitSystem, c_Coupling] := 
  PlanckReduced[u, c]/ElectronMass[u, c]/LightSpeed[u, c]/FineStructure[c];
BohrReduced[u_UnitSystem, c_Coupling] := Bohr[u, c] (1 + 1/ProtonElectron[c]);
ElectronRadius[u_UnitSystem, c_Coupling] := 
  FineStructure[c] PlanckReduced[u, c]/ElectronMass[u, c]/LightSpeed[u, c];
Magneton[u_UnitSystem, c_Coupling] := 
  Charge[u, c] PlanckReduced[u, c] Lorentz[u]/2 ElectronMass[u, c];

(* more *)
(*
\[Kappa] = Einstein[UnitSystem["SI2019"]];
\[Sigma] = Stefan[UnitSystem["SI2019"]];(**)
\[Mu]B = Magneton[UnitSystem["SI2019"]];(**)
\[CurlyEpsilon]0 = Permittivity[UnitSystem["SI2019"]];(**)

ke = Coulomb[UnitSystem["SI2019"]];(**)
mp = ProtonMass[UnitSystem["SI2019"]];
mu = AtomicMass[UnitSystem["SI2019"]];
Mu = MolarMass[UnitSytem["SI2019"]];
\|01d509 = Faraday[UnitSystem["SI2019"]];(**)
\[CapitalPhi]0 = MagneticFlux[UnitSystem["SI2019"]];(**)
Z0 = Impedance[UnitSystem["SI2019"]];(**)

G0 = Conductance[UnitSystem["SI2019"]];(**)
Eh = Hartree[UnitSystem["SI2019"]];
a0 = Bohr[UnitSystem["SI2019"]];
re = ElectronRadius[UnitSystem["SI2019"]];
RK = Klitzing[UnitSystem["SI2019"]];(**)
KJ = Josephson[UnitSystem["SI2019"]];(**)
{RH, Ry} = {R\[Infinity] mp/(me + mp), h c R\[Infinity]};

\[ScriptL]P = Length[UnitSystem["PlanckGauss"], UnitSystem["SI2019"]];
tP = Time[UnitSystem["PlanckGauss"], UnitSystem["SI2019"]];
TP = Temperature[UnitSystem["PlanckGauss"], UnitSystem["SI2019"]];

lS = Length[UnitSystem["Stoney"], UnitSystem["SI2019"]];
tS = Time[UnitSystem["Stoney"], UnitSystem["SI2019"]];
mS = Mass[UnitSystem["Stoney"], UnitSystem["SI2019"]];
qS = Charge[UnitSystem["Stoney"], UnitSystem["SI2019"]];

(*lA=Length[UnitSystem["Hartree"],UnitSystem["SI2019"]];
tA=Time[UnitSystem["Hartree"],UnitSystem["SI2019"]];
mA=Mass[UnitSystem["Hartree"],UnitSystem["SI2019"]];
qA=Charge[UnitSystem["Hartree"],UnitSystem["SI2019"]];*)

lQCD = Length[UnitSystem["QCD"], UnitSystem["SI2019"]];
tQCD = Time[UnitSystem["QCD"], UnitSystem["SI2019"]];
mQCD = Mass[UnitSystem["QCD"], UnitSystem["SI2019"]];
*)
