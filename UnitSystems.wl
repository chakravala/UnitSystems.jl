(* ::Package:: *)
(* This file is part of UnitSystems. It is licensed under the MIT license *)
(* UnitSystems Copyright (C) 2021 Michael Reed *)

Unprotect[UnitSystem];
ProtectedList = {Length, Area, Volume, Power, Entropy};
UnitSystemsList = {Metric, SI2019, CODATA, Conventional, MTS, English,
    EnglishUS, IAU, SI1976, Mixed, ESU2019, EMU2019, EMU, ESU, Gauss, 
   LorentzHeaviside, Thomson, Kennelly, Planck, PlanckGauss, Stoney, 
   Hartree, Rydberg, Schrodinger, Electronic, Natural, NaturalGauss, 
   QCD, QCDGauss, QCDoriginal};
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

measure[x_] := x;
GravityCoupling[Coupling[\[Alpha]G_, ___]] := measure[\[Alpha]G];
FineStructure[Coupling[_, \[Alpha]_, ___]] := measure[\[Alpha]];
ElectronUnit[Coupling[_, _, \[Mu]eu_, ___]] := measure[\[Mu]eu];
ProtonUnit[Coupling[_, _, _, \[Mu]pu_, ___]] := measure[\[Mu]pu];
ProtonElectron[c_Coupling] := ProtonUnit[c]/ElectronUnit[c];

UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, \[Lambda]_] := 
  UnitSystem[kB, \[HBar], c, \[Mu]0, me, \[Lambda], 1];
UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_] := 
  UnitSystem[kB, \[HBar], c, \[Mu]0, me, 1];
Boltzmann[UnitSystem[kB_, ___]] := kB;
PlanckReduced[UnitSystem[_, \[HBar]_, ___]] := \[HBar];
LightSpeed[UnitSystem[_, _, c_, ___]] := c;
Permeability[UnitSystem[_, _, _, \[Mu]0_, ___]] := \[Mu]0;
ElectronMass[UnitSystem[_, _, _, _, me_, ___]] := me;
Rationalization[UnitSystem[_, _, _, _, _, \[Lambda]_, ___]] := \[Lambda];
Lorentz[UnitSystem[_, _, _, _, _, _, \[Alpha]L_, ___]] := \[Alpha]L;

RationalizedQ[u_UnitSystem] := Rationalization[u] != 4 \[Pi]

Universe[_] := StandardModel;
Unit[x_, y_ : 1] := x;
Mass[u_UnitSystem, s_UnitSystem] := ElectronMass[u, s];
ElectronMass[h_?NumberQ] := R\[Infinity] 2 h \[Alpha]inv/c;
ElectronMass[h_?NumberQ, u_Coupling] := (R\[Infinity] 2 h)/(FineStructure[u]^2 c);
PlanckMass[u_UnitSystem, c_Coupling] := ElectronMass[u, c]/Sqrt[GravityCoupling[c]];
Planck[u_UnitSystem, c_Coupling] := 2 \[Pi] PlanckReduced[u];
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
      Module[{n = #[u, s]}, If[OneQ[v], v, v/u]];
     #[v_?NumberQ, 
       u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___], 
       s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := v;
     #[u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___], 
       s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := 1;
     If[!MemberQ[Join[ProtectedList, ConstantsList, {Permittivity, Charge, MagneticFlux, 
       Impedance, Conductance}], #], #[u_UnitSystem] := #[u, Metric], Nothing];,
     UnitSystem /: #[v_?NumberQ, u_UnitSystem] := #[v, u, Metric];
     UnitSystem /: #[v_?NumberQ, u_UnitSystem, s_UnitSystem] := 
      Module[{n = #[u, s]}, If[OneQ[v], v, v/u]];
     UnitSystem /: #[v_?NumberQ, 
       u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___], 
       s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := v;
     UnitSystem /: #[
       u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___], 
       s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := 1;
     UnitSystem /: #[u_UnitSystem] := #[u, Metric];]
    ) &, ConvertList];

g0 = 9.80665;
ft = 0.3048;
ftUS = 1200/3937;
lbm = g0/ft;
lbmUS = g0/ftUS;
slug = 0.45359237 lbm;
slugUS = 0.45359237 lbmUS;
rankine = 5/9;
kelvin = 9/5;
atm = 101325;

kcalth = 4184;
kcal4 = 4204;
kcal10 = 4185.5;
kcal20 = 4182;
kcalm = 4190;
kcalit = 4186.8;
calth = kcalth/1000;
cal4 = kcal4/1000;
cal10 = kcal10/1000;
cal20 = kcal20/1000;
calm = kcalm/1000;
calit = kcalit/1000;
kcal = kcalth;
cal = kcal/1000;
(*calth = thermal calorie*)

\[CapitalDelta]\[Nu]Cs = 9192631770.0;
Kcd = 683 555.016/555;
mP = 2.176434 10^-8;
NA = 6.02214076 10^23;
kB = 1.380649 10^-23;
h = 6.62607015 10^-34;
c = 299792458.;
e = 1.602176634 10^-19;
\[Mu]eu = 1/1822.888486209;
\[Mu]pu = 1.007276466621;
\[Alpha]inv = 137.035999084;
R\[Infinity] = 10973731.5681601;
me = ElectronMass[h];
\[Mu]0 = 2 h/c/\[Alpha]inv/e^2 ;(*\[TildeTilde]4\[Pi]*(1e-7+5.5e-17),exact charge*)
ħ = h/2 \[Pi];
\[Delta]\[Mu]0 = \[Mu]0 - 4 \[Pi] 10^-7;
\[Mu]pe = \[Mu]pu/\[Mu]eu;
Ru = NA kB;
\[Alpha]L = 0.01/c;
\[Alpha]G = (me/mP)^2;
RK1990 = 25812.807;
RK2014 = 25812.8074555;
KJ1990 = 4.835979 10^14;
KJ2014 = 4.835978525 10^14;
ħ1990 = 2/RK1990/KJ1990^2/\[Pi];
ħ2014 = 2/RK2014/KJ2014^2/\[Pi];
me1990 = ElectronMass[2 \[Pi] ħ1990];
me2014 = ElectronMass[2 \[Pi] ħ2014];

StandardModel = Coupling[\[Alpha]G, 1/\[Alpha]inv, \[Mu]eu, \[Mu]pu];
Gauss = UnitSystem[10^10 Ru me/\[Mu]eu, 10^7 ħ, 100 c, 1, 1000 me, 4 \[Pi], 0.01/c];
LorentzHeaviside = UnitSystem[10^10 Ru me/\[Mu]eu, 10^7 ħ, 100 c, 1, 1000 me, 1, 0.01/c];
Thomson = UnitSystem[10^10 Ru me/\[Mu]eu, 10^7 ħ, 100 c, 1, 1000 me, 4 \[Pi], 1/2];
Kennelly = UnitSystem[Ru me/\[Mu]eu/0.001, ħ, c, 10^-7, me, 4 \[Pi]];
ESU = UnitSystem[10^10 Ru me/\[Mu]eu, 10^7 ħ, 100 c, (100 c)^-2, 1000 me, 4 \[Pi]];
ESU2019 = UnitSystem[10^7 kB, 10^7 ħ, 100 c, 10^3 \[Mu]0/c^2, 1000 me];
EMU = UnitSystem[10^10 Ru me/\[Mu]eu, 10^7 ħ, 100 c, 1, 1000 me, 4 \[Pi]];
EMU2019 = UnitSystem[10^7 kB, 10^7 ħ, 100 c, 10^7 \[Mu]0, 1000 me];
MTS = UnitSystem[10^6 Ru me/\[Mu]eu, 1000 ħ, c, 4 \[Pi]/10^4, me/1000];
Mixed = UnitSystem[Ru me/\[Mu]eu/0.001, ħ, c, \[Mu]0, me];
Metric = UnitSystem[Ru me/\[Mu]eu/0.001, ħ, c, 4 \[Pi] 10^-7, me];
SI1976 = UnitSystem[8.31432 me/\[Mu]eu/0.001, ħ, c, 4 \[Pi] 10^-7, me];
SI2019 = UnitSystem[kB, ħ, c, \[Mu]0, me];
CODATA = UnitSystem[Ru me2014/\[Mu]eu/0.001, ħ2014, c, 2 RK2014/c/\[Alpha]inv, me2014];
Conventional = 
  UnitSystem[Ru me1990/\[Mu]eu/0.001, ħ1990, c, 2 RK1990/c/\[Alpha]inv, me1990];
English = 
  UnitSystem[kB rankine/slug/ft^2, ħ/slug/ft^2, c/ft, 4 \[Pi], me/slug];
EnglishUS = 
  UnitSystem[(1000 Ru me/\[Mu]eu ) (rankine/slug/ftUS^2), 
   ħ/slug/ftUS^2, c/ftUS, 4 \[Pi], me/slug];

GMsun = 1.32712442099 10^20;
GMearth = 398600441.8 10^6;
GMjupiter = 1.26686534 10^17;
au = 149597870.7 10^3;
LD = 384402 10^3;
day = 60^2 24;
pc = au 648000/\[Pi];
ly = 365.25 c day;
GG = c ħ/mP^2;
ms = GMsun/GG;
Js = ms au^2/day^2;
IAU = UnitSystem[Ru me/\[Mu]eu/0.001/Js, ħ/day/Js, day c/au, 
   4 \[Pi] 10^-7 day^2/Js, me/ms];

mf = Mass[90/lbm, Metric, English];
Jf = mf (201.168/14 day)^2;
FFF = UnitSystem[1000 Ru me/\[Mu]eu rankine/Jf, ħ/14 day/Jf, 14 day c/201.168, 0, me/mf];
SI = SI2019;
MKS = Metric;
{CGS, CGS2019, CGSm, CGSe, HLU} = {Gauss, EMU2019, EMU, ESU, LorentzHeaviside};

mf = Mass[90/lbm, Metric, English];
Jf = mf (201.168/14 day)^2;
FFF = UnitSystem[1000 Ru me/\[Mu]eu rankine/Jf, ħ/14 day/Jf, 14 day c/201.168, 0, me/mf];
SI = SI2019;
MKS = Metric;
{CGS, CGS2019, CGSm, CGSe, HLU} = {Gauss, EMU2019, EMU, ESU, LorentzHeaviside};

(*Planck=UnitSystem[1,1,1,1,Sqrt[4\[Pi] \[Alpha]G]];*)

PlanckGauss = UnitSystem[1, 1, 1, 4 \[Pi], Sqrt[\[Alpha]G]];
Stoney = UnitSystem[1, \[Alpha]inv, 1, 4 \[Pi], Sqrt[\[Alpha]G \[Alpha]inv]];
(*Hartree=UnitSystem[1,1,\[Alpha]inv,4\[Pi]/\[Alpha]inv^2,1];
Rydberg=UnitSystem[1,1,2\[Alpha]inv,\[Pi]/\[Alpha]inv^2,1/2];*)

Schrodinger = 
  UnitSystem[1, 1, \[Alpha]inv, 4 \[Pi]/\[Alpha]inv^2, 
   Sqrt[\[Alpha]G \[Alpha]inv]];
Electronic = UnitSystem[1, \[Alpha]inv, 1, 4 \[Pi], 1];
Natural = UnitSystem[1, 1, 1, 1, 1];
NaturalGauss = UnitSystem[1, 1, 1, 4 \[Pi], 1];
QCD = UnitSystem[1, 1, 1, 1, 1/\[Mu]pe];
QCDGauss = UnitSystem[1, 1, 1, 4 \[Pi], 1/\[Mu]pe];
QCDoriginal = UnitSystem[1, 1, 1, 4 \[Pi]/\[Alpha]inv, 1/\[Mu]pe];

ElectronMass[Planck, c_Coupling] := Sqrt[4 \[Pi] GravityCoupling[c]];
ElectronMass[PlanckGauss, c_Coupling] := Sqrt[GravityCoupling[c]];
ElectronMass[UnitSystem[_, _, _, _, Sqrt[\[Alpha]G \[Alpha]inv], ___],
   c_Coupling] := Sqrt[GravityCoupling[c]/FineStructure[c];]
ElectronMass[UnitSystem[_, _, _, _, 1/\[Mu]pe, ___], c_Coupling] := 1/ProtonElectron[c];
Permeability[UnitSystem[_, _, _, 4 \[Pi]/\[Alpha]inv^2, ___], 
   c_Coupling] := 4 \[Pi] FineStructure[c]^2;
Permeability[UnitSystem[_, _, _, \[Pi]/\[Alpha]inv^2, ___], 
   c_Coupling] := \[Pi] FineStructure[c]^2;
LightSpeed[UnitSystem[_, _, \[Alpha]inv, ___], c_Coupling] := 1/FineStructure[c];
LightSpeed[UnitSystem[_, _, 2 \[Alpha]inv, ___], c_Coupling] := 2/FineStructure[c];
PlanckReduced[UnitSystem[_, \[Alpha]inv, ___], c_Coupling] := 1/FineStructure[c];

ElectronMass[u : UnitSystem[_, _, c, _, me, ___], c_Coupling] := 
  ElectronMass[Planck[u], c];
ElectronMass[UnitSystem[_, _, 100 c, _, 1000 me, ___], c_Coupling] := 
  1000 ElectronMass[SI, c];
ElectronMass[UnitSystem[_, _, c, _, me/1000, ___], c_Coupling] := ElectronMass[SI, c]/1000;
ElectronMass[u : UnitSystem[_, _, c, _, me2014, ___], c_Coupling] := 
  ElectronMass[Planck[u], c];
ElectronMass[u : UnitSystem[_, _, c, \[Mu]0, me1990, ___], 
   c_Coupling] := ElectronMass[Planck[u], c];
ElectronMass[UnitSystem[_, _, c/ftUS, _, me/slug, ___], c_Coupling] :=
   ElectronMass[SI, c]/slug;
Permeability[UnitSystem[_, _, _, \[Mu]0, ___], c_Coupling] := FineStructure[c] 2 h/c/e^2;
Permeability[ESU2019, u_Coupling] := 10^3 Permeability[SI, u]/c^2;
Permeability[EMU2019, c_Coupling] := 10^7 Permeability[SI, c];
Permeability[CODATA, u_Coupling] := 2 RK2014 FineStructure[u]/c;
Permeability[Conventional, u_Coupling] := 2 RK1990 FineStructure[u]/c;

MolarMass[UnitSystem[1, ___]] = 1;
MolarMass[u : UnitSystem[kB, ___]] := MolarMass[u, Universe[u]];
MolarMass[u : UnitSystem[kB, ___], c_Coupling] := NA ElectronMass[u, c]/ElectronUnit[c];
MolarMass[u : UnitSystem[10^7 kB, ___]] := MolarMass[u, Universe[u]];
MolarMass[UnitSystem[10^7 kB, ___], c_Coupling] := 1000 MolarMass[SI2019, c];
MolarMass[u : UnitSystem[10^3 kB, ___]] := MolarMass[u, Universe[u]];
MolarMass[UnitSystem[10^3 kB, ___], c_Coupling] := MolarMass[SI2019, c]/1000;
MolarMass[UnitSystem[kB_, ___]] := MolarMass[CGS]/1000;
MolarMass[UnitSystem[Boltzmann[MTS], ___]] := MolarMass[CGS]/10^6;
MolarMass[UnitSystem[Boltzmann[CGS], ___]] := MolarMass[Natural];
MolarMass[UnitSystem[Boltzmann[FFF], ___]] := MolarMass[Natural];
MolarMass[u : UnitSystem[Boltzmann[English], ___]] := MolarMass[u, Universe[u]];
MolarMass[u : UnitSystem[Boltzmann[English], ___], c_Coupling] := 
  1000 MolarMass[SI2019, c];
MolarMass[UnitSystem[Boltzmann[EnglishUS], ___]] := MolarMass[Natural];
MolarMass[UnitSystem[Boltzmann[IAU], ___]] := 1/1000 ms;

LuminousEfficacy[UnitSystem[1, ___]] = 1;
LuminousEfficacy[u : UnitSystem] := Power[Kcd, SI2019, u];

Kilograms[m_] := Kilograms[m, English];
Kilograms[m_, u_UnitSystem] := Mass[m, Metric, u];
Slugs[m_] := Slugs[m, Metric];
Slugs[m_, u_UnitSystems] := Mass[m, English, u];
Feet[d_] := Feet[d, Metric];
Feet[d_, u_UnitSystem] := Length[d, English, u];
Meters[d_] := Meters[d, English];
Meters[d_, u_UnitSystem] := Length[d, Metric, u];

(* IAU to SI *)
UnitSystem /: 
  Length[u : UnitSystem[_, _, c, ___], 
   s : UnitSystem[_, _, day c/au, ___]] := Length[u, s, 1/au];
Time[u : UnitSystem[_, _, c, ___], 
   s : UnitSystem[_, _, day c/au, ___]] := Time[u, s, 1/day];
(* SI to IAU *)
UnitSystem /: 
  Length[u : UnitSystem[_, _, day c/au, ___], 
   s : UnitSystem[_, _, c, ___]] := Length[u, s, au];
Time[u : UnitSystem[_, _, day c/au, ___], 
   s : UnitSystem[_, _, c, ___]] := Time[u, s, day];
(* IAU to CGS *)
UnitSystem /: 
  Length[u : UnitSystem[_, _, 100 c, ___], 
   s : UnitSystem[_, _, day c/au, ___]] := Length[u, s, 1/au];
Time[u : UnitSystem[_, _, 100 c, ___], 
   s : UnitSystem[_, _, day c/au, ___]] := Time[u, s, 1/day];
(* CGS to IAU *)
UnitSystem /: 
  Length[u : UnitSystem[_, _, day c/au, ___], 
   s : UnitSystem[_, _, 100 c, ___]] := Length[u, s, au];
Time[u : UnitSystem[_, _, day c/au, ___], 
   s : UnitSystem[_, _, 100 c, ___]] := Time[u, s, day];
(* IAU to English *)
UnitSystem /: 
  Length[u : UnitSystem[_, _, c/ft, ___], 
   s : UnitSystem[_, _, day c/au, ___]] := Length[u, s, ft/au];
Time[u : UnitSystem[_, _, c/ft, ___], 
   s : UnitSystem[_, _, day c/au, ___]] := Time[u, s, 1/day];
(* English to IAU *)
UnitSystem /: 
  Length[u : UnitSystem[_, _, day c/au, ___], 
   s : UnitSystem[_, _, c/ft, ___]] := Length[u, s, au/ft];
Time[u : UnitSystem[_, _, day c/au, ___], 
   s : UnitSystem[_, _, c/ft, ___]] := Time[u, s, day];

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

Hyperfine[u_UnitSystem] := Frequency[\[CapitalDelta]\[Nu]Cs, u];

Map[(#[u_UnitSystem] := #[u, Universe[u]]) &, {Avogadro, AtomicMass, 
   ProtonMass, Einstein, Universal, Stefan, RadiationDensity, 
   Permittivity, Coulomb, BiotSavart, Impedance, Faraday, Josephson, 
   MagneticFlux, Klitzing, Conductance, Hartree, Rydberg, Bohr, 
   BohrReduced, ElectronRadius, Magneton}];
Avogadro[u_UnitSystem, c_Coupling] := 
  MolarMass[u, c] ElectronUnit[c]/ElectronMass[u, c];
AtomicMass[u_UnitSystem, c_Coupling] := ElectronMass[u, c]/ElectronUnit[c];
ProtonMass[u_UnitSystem, c_Coupling] := ProtonElectron[c] ElectronMass[u, c];
Einstein[u_UnitSystem, c_Coupling] := (8 \[Pi] Newton[u, c])/LightSpeed[u, c]^4;
Universal[u_UnitSystem, c_Coupling] := Boltzmann[u, c] Avogadro[u, c];
Stefan[u_UnitSystem, c_Coupling] := (2 \[Pi]^5 Boltzmann[u, c]^4)/(
  15 Planck[u, c]^3 LightSpeed[u, c]^2);
RadiationDensity[u_UnitSystem, c_Coupling] := (4 Stefan[u, c])/LightSpeed[u, c];
Permittivity[u_UnitSystem, c_Coupling] := 
  (1/Permeability[u, c]) (LightSpeed[u, c] Lorentz[u]^2);
Coulomb[u_UnitSystem, c_Coupling] := 
  Rationalization[u]/4 \[Pi]/ Permittivity[u];
BiotSavart[u_UnitSystem, c_Coupling] := 
  Permeability[u, c] Lorentz[u] (Rationalization[u]/4 \[Pi]);
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

\[Kappa] = Einstein[SI2019];
\[Sigma] = Stefan[SI2019];(**)
\[Mu]B = Magneton[SI2019];(**)
\[CurlyEpsilon]0 = Permittivity[SI2019];(**)

ke = Coulomb[SI2019];(**)
mp = ProtonMass[SI2019];
mu = AtomicMass[SI2019];
Mu = MolarMass[SI2019];
\|01d509 = Faraday[SI2019];(**)
\[CapitalPhi]0 = MagneticFlux[SI2019];(**)
Z0 = Impedance[SI2019];(**)

G0 = Conductance[SI2019];(**)
Eh = Hartree[SI2019];
a0 = Bohr[SI2019];
re = ElectronRadius[SI2019];
RK = Klitzing[SI2019];(**)
KJ = Josephson[SI2019];(**)
{RH, Ry} = {R\[Infinity] mp/(me + mp), h c R\[Infinity]};

\[ScriptL]P = Length[PlanckGauss, SI2019];
tP = Time[PlanckGauss, SI2019];
TP = Temperature[PlanckGauss, SI2019];

lS = Length[Stoney, SI2019];
tS = Time[Stoney, SI2019];
mS = Mass[Stoney, SI2019];
qS = Charge[Stoney, SI2019];

(*lA=Length[Hartree,SI2019];
tA=Time[Hartree,SI2019];
mA=Mass[Hartree,SI2019];
qA=Charge[Hartree,SI2019];*)

lQCD = Length[QCD, SI2019];
tQCD = Time[QCD, SI2019];
mQCD = Mass[QCD, SI2019];
