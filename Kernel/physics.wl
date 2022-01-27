
(* This file is part of UnitSystems. It is licensed under the MIT license *)
(* UnitSystems Copyright (C) 2022 Michael Reed *)

ElectronMass[h_?NumberQ] := UnitData["R\[Infinity]"] 2 h/UnitData["\[Alpha]"]^2/UnitData["c"];
ElectronMass[h_?NumberQ, u_Coupling] := (UnitData["R\[Infinity]"] 2 h)/(FineStructureConstant[u]^2 UnitData["c"]);
PlanckMass[u_UnitSystem, c_Coupling] := PowerExpand[ElectronMass[u, c]/Sqrt[GravitationalCouplingConstantElectronElectron[c]]]
PlanckConstant[u_UnitSystem, c_Coupling] := 2 Pi ReducedPlanckConstant[u];
GravitationalConstant[u_UnitSystem, c_Coupling] := PowerExpand[(SpeedOfLight[u, c] ReducedPlanckConstant[u, c])/PlanckMass[u, c]^2]
ElementaryCharge[u_UnitSystem, c_Coupling] := PowerExpand[Sqrt[2 PlanckConstant[u]/(MagneticConstant[u]/FineStructureConstant[u])/(SpeedOfLight[u] RationalizationConstant[u] LorentzConstant[u]^2)]]

ElectronMass[UnitSystem["Planck"], c_Coupling] := PowerExpand[Sqrt[4 Pi GravityCoupling[c]]]
ElectronMass[UnitSystem["PlanckGauss"], c_Coupling] := Sqrt[GravitationalCouplingConstantElectronElectron[c]];
ElectronMass[UnitSystem[_, _, _, _, Sqrt[EvalUnitData["\[Alpha]G"]/UnitData["\[Alpha]"]], ___], c_Coupling] := PowerExpand[Sqrt[GravitationalCouplingConstantElectronElectron[c]/FineStructureConstant[c]]]
ElectronMass[UnitSystem[_, _, _, _, 1/UnitData["\[Mu]pe"], ___], c_Coupling] := 1/ProtonElectronMassRatio[c];
MagneticConstant[UnitSystem[_, _, _, 4 Pi UnitData["\[Alpha]"]^2, ___], c_Coupling] := 4 Pi FineStructureConstant[c]^2;
MagneticConstant[UnitSystem[_, _, _, Pi UnitData["\[Alpha]"]^2, ___], c_Coupling] := Pi FineStructureConstant[c]^2;
SpeedOfLight[UnitSystem[_, _, 1/UnitData["\[Alpha]"], ___], c_Coupling] := 1/FineStructureConstant[c];
SpeedOfLight[UnitSystem[_, _, 2/UnitData["\[Alpha]"], ___], c_Coupling] := 2/FineStructureConstant[c];
PlanckReduced[UnitSystem[_, 1/UnitData["\[Alpha]"], ___], c_Coupling] := 1/FineStructure[c];

ElectronMass[u : UnitSystem[_, _, UnitData["c"], _, EvalUnitData["me"], ___], c_Coupling] := ElectronMass[PlanckConstant[u], c];
ElectronMass[UnitSystem[_, _, 100 UnitData["c"], _, 1000 EvalUnitData["me"], ___], c_Coupling] := 1000 ElectronMass[UnitSystem["SI2019"], c];
ElectronMass[UnitSystem[_, _, UnitData["c"], _, EvalUnitData["me"]/1000, ___], c_Coupling] := ElectronMass[UnitSystem["SI2019"], c]/1000;
ElectronMass[u : UnitSystem[_, _, UnitData["c"], _, EvalUnitData["me2014"], ___], c_Coupling] := ElectronMass[PlanckConstant[u], c];
ElectronMass[u : UnitSystem[_, _, UnitData["c"], EvalUnitData["\[Mu]0"], EvalUnitData["me1990"], ___], c_Coupling] := ElectronMass[PlanckConstant[u], c];
ElectronMass[UnitSystem[_, _, UnitData["c"]/UnitData["ftUS"], _, EvalUnitData["me"]/EvalUnitData["slug"], ___], c_Coupling] := ElectronMass[UnitSystem["SI2019"], c]/EvalUnitData["slug"];

MagneticConstant[UnitSystem[_, _, _, EvalUnitData["\[Mu]0"], ___], c_Coupling] := FineStructureConstant[c] 2 UnitData["h"]/UnitData["c"]/UnitData["e"]^2;
MagneticConstant[UnitSystem["ESU2019"], u_Coupling] := 10^3 MagneticConstant[UnitSystem["SI2019"], u]/UnitData["c"]^2;
MagneticConstant[UnitSystem["EMU2019"], u_Coupling] := 10^7 MagneticConstant[UnitSystem["SI2019"], u];
MagneticConstant[UnitSystem["CODATA"], u_Coupling] := 2 UnitData["RK2014"] FineStructureConstant[u]/UnitData["c"];
MagneticConstant[UnitSystem["Conventional"], u_Coupling] := 2 UnitData["RK1990"] FineStructureConstant[u]/UnitData["c"];

MolarMassConstant[UnitSystem[1, ___]] = 1;
MolarMassConstant[u:UnitSystem[UnitData["kB"], ___]] := MolarMassConstant[u, Universe[u]];
MolarMassConstant[u:UnitSystem[UnitData["kB"], ___], c_Coupling] := UnitData["NA"] ElectronMass[u, c]/ElectronRelativeAtomicMass[c];
MolarMassConstant[u:UnitSystem[10^7 UnitData["kB"], ___]] := MolarMassConstant[u, Universe[u]];
MolarMassConstant[UnitSystem[10^7 UnitData["kB"], ___], c_Coupling] := 1000 MolarMassConstant[UnitSystem["SI2019"], c];
MolarMassConstant[u:UnitSystem[10^3 UnitData["kB"], ___]] := MolarMassConstant[u, Universe[u]];
MolarMassConstant[UnitSystem[10^3 UnitData["kB"], ___], c_Coupling] := MolarMassConstant[UnitSystem["SI2019"], c]/1000;
MolarMassConstant[UnitSystem[kB_, ___]] := MolarMassConstant["CGS"]/1000;
MolarMassConstant[UnitSystem[BoltzmannConstant[UnitSystem["MTS"]], ___]] := MolarMassConstant["CGS"]/10^6;
MolarMassConstant[UnitSystem[BoltzmannConstant[UnitSystem["CGS"]], ___]] := MolarMassConstant["Natural"];
MolarMassConstant[UnitSystem[BoltzmannConstant[UnitSystem["FFF"]], ___]] := MolarMassConstant["Natural"];
MolarMassConstant[u:UnitSystem[BoltzmannConstant[UnitSystem["English"]], ___]] := MolarMassConstant[u, Universe[u]];
MolarMassConstant[u:UnitSystem[BoltzmannConstant[UnitSystem["English"]], ___], c_Coupling] := 1000 MolarMassConstant[UnitSystem["SI2019"],c]
MolarMassConstant[UnitSystem[BoltzmannConstant[UnitSystem["EnglishUS"]], ___]] := MolarMassConstant["Natural"];
MolarMassConstant[UnitSystem[BoltzmannConstant[UnitSystem["IAU"]], ___]] = EvalUnitData["ms"]/1000;

MolarMassConstant[AbstractUnitSystem["AbstractUnits"]] = "Mu"
MolarMassConstant[AbstractUnitSystem["AbstractUnits1"]] = "Mu1"
MolarMassConstant[AbstractUnitSystem["AbstractUnits2"]] = "Mu2"
MolarMassConstant[u : UnitSystem["kB", ___]] := MolarMassConstant[u,Universe[u]]
MolarMassConstant[u : UnitSystem["kB", ___], c_Coupling] := "NA" ElectronMass[AbstractUnitSystem["SI2019"],c]/ElectronRelativeAtomicMass[c]
MolarMassConstant[u : UnitSystem[10^7 "kB", ___]] := MolarMassConstant[u, Universe[u]];
MolarMassConstant[UnitSystem[10^7 "kB", ___], c_Coupling] := 1000 MolarMassConstant[AbstractUnitSystem["SI2019"], c];
MolarMassConstant[u : UnitSystem[10^3 "kB", ___]] := MolarMassConstant[u, Universe[u]];
MolarMassConstant[UnitSystem[10^3 "kB", ___], c_Coupling] := MolarMassConstant[AbstractUnitSystem["SI2019"], c]/1000;
MolarMassConstant[UnitSystem[BoltzmannConstant[AbstractUnitSystem["MTS"]], ___]] := MolarMassConstant["CGS"]/10^6;
MolarMassConstant[UnitSystem[BoltzmannConstant[AbstractUnitSystem["CGS"]], ___]] := MolarMassConstant["Natural"];
MolarMassConstant[UnitSystem[BoltzmannConstant[AbstractUnitSystem["FFF"]], ___]] := MolarMassConstant["Natural"];
MolarMassConstant[u : UnitSystem[BoltzmannConstant[AbstractUnitSystem["English"]], ___]] := MolarMassConstant[u, Universe[u]];
MolarMassConstant[u : UnitSystem[BoltzmannConstant[AbstractUnitSystem["English"]], ___], c_Coupling] := 1000 MolarMassConstant[AbstractUnitSystem["SI2019"], c];
MolarMassConstant[UnitSystem[BoltzmannConstant[AbstractUnitSystem["EnglishUS"]], ___]] := MolarMassConstant["Natural"];
MolarMassConstant[UnitSystem[BoltzmannConstant[AbstractUnitSystem["IAU"]], ___]] := 1/1000 AbstractUnitData["ms"];

MonochromaticRadiation540THzLuminousEfficacy[UnitSystem[1, ___]] = 1
MonochromaticRadiation540THzLuminousEfficacy[u : UnitSystem[_?NumericQ, ___]] := Power[UnitData["Kcd"], UnitSystem["SI2019"], u]
MonochromaticRadiation540THzLuminousEfficacy[u : UnitSystem[_Around, ___]] := Power[UnitData["Kcd"], UnitSystem["SI2019"], u]
MonochromaticRadiation540THzLuminousEfficacy[AbstractUnitSystem["AbstractUnits1"]] = "Kcd1"
MonochromaticRadiation540THzLuminousEfficacy[AbstractUnitSystem["AbstractUnits2"]] = "Kcd2"
MonochromaticRadiation540THzLuminousEfficacy[u_UnitSystem] := Power["Kcd", AbstractUnitSystem["SI2019"], u]

Kilograms[m_] := Kilograms[m, "English"];
Kilograms[m_, u_UnitSystem] := Mass[m, UnitSystem["Metric"], u];
Kilograms[m_, u_String] := Mass[m, UnitSystem[u]];
Slugs[m_] := Slugs[m, "Metric"];
Slugs[m_, u_UnitSystems] := Mass[m, UnitSystem["English"], u];
Slugs[m_, u_String] := Mass[m, UnitSystem[u]];
Feet[d_] := Feet[d, "Metric"];
Feet[d_, u_UnitSystem] := Length[d, UnitSystem["English"], u];
Feet[d_, u_String] := Length[d, UnitSystem[u]];
Meters[d_] := Meters[d, "English"];
Meters[d_, u_UnitSystem] := Length[d, UnitSystem["Metric"], u];
Meters[d_, u_String] := Length[d, UnitSystem[u]];

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
  Unit[(ReducedPlanckConstant[s] ElectronMass[u] SpeedOfLight[
       u])/(ReducedPlanckConstant[u] ElectronMass[s] SpeedOfLight[s]), l];
UnitSystem /: Area[u_UnitSystem, s_UnitSystem] := Unit[Length[u, s]^2];
UnitSystem /: Volume[u_UnitSystem, s_UnitSystem] := Unit[Length[u, s]^3];

WaveNumber[u_UnitSystem, s_UnitSystem] := Unit[Length[s, u]];
FuelEfficiency[u_UnitSystem, s_UnitSystem] := Area[s, u];
Time[u_UnitSystem, s_UnitSystem] := Time[u, s, 1];
Time[u_UnitSystem, s_UnitSystem, t_] := Unit[Length[u, s]/SpeedOfLight[u, s], 1];
Frequency[u_UnitSystem, s_UnitSystem] := Time[s, u];
FrequencyDrift[u_UnitSystem, s_UnitSystem] := Unit[Time[s, u]^2];
Speed[u_UnitSystem, s_UnitSystem] := SpeedOfLight[u, s];
Acceleration[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]/Time[u, s]];
Jerk[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]/Time[u, s]^2];
Snap[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]/Time[u, s]^3];
VolumeFlow[u_UnitSystem, s_UnitSystem] := Unit[Area[u, s], Speed[u, s]];
SpecificEnergy[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]^2];

Mass[u_UnitSystem, s_UnitSystem] := ElectronMass[u, s];
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
Diffusivity[u_UnitSystem, s_UnitSystem] := Unit[(ReducedPlanckConstant[s] ElectronMass[u])/(ReducedPlanckConstant[u] ElectronMass[s])];
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
  Unit[Sqrt[(ReducedPlanckConstant[s] MagneticConstant[u] SpeedOfLight[
        u] RationalizationConstant[u] LorentzConstant[u]^2)/(ReducedPlanckConstant[
        u] MagneticConstant[s] SpeedOfLight[s] RationalizationConstant[s] LorentzConstant[s]^2)]];
Current[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s]/Time[u, s]];
ElectricPotential[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/Charge[u, s]];
Capacitance[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s]/ElectricPotential[u, s]];
Resistance[u_UnitSystem, s_UnitSystem] := Unit[ElectricPotential[u, s]/Current[u, s]];
Conductance[u_UnitSystem, s_UnitSystem] := Unit[Current[u, s]/ElectricPotential[u, s]];
MagneticFlux[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/LorentzConstant[u, s]/Current[u, s]];
MagneticFluxDensity[u_UnitSystem, s_UnitSystem] :=
  Unit[Mass[u, s]/LorentzConstant[u, s]/Current[u, s]/Time[u, s]^2];
Inductance[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] Area[u, s]/Charge[u, s]^2];

ElectricFluxDensity[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s] RationalizationConstant[u, s]/Area[u, s]];
ChargeDensity[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s]/Volume[u, s]];
CurrentDensity[u_UnitSystem, s_UnitSystem] := Unit[Current[u, s]/Area[u, s]];
Conductivity[u_UnitSystem, s_UnitSystem] := Unit[Conductance[u, s]/Length[u, s]];
Permittivity[u_UnitSystem, s_UnitSystem] := Unit[Capacitance[u, s] RationalizationConstant[u, s]/Length[u, s]];
ElectricField[u_UnitSystem, s_UnitSystem] := Unit[ElectricPotential[u, s]/Length[u, s]];
MagneticField[u_UnitSystem, s_UnitSystem] := Unit[Current[u, s] RationalizationConstant[u, s] LorentzConstant[u, s]/Length[u, s]];
Exposure[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s]/Mass[u, s]];
Resistivity[u_UnitSystem, s_UnitSystem] := Unit[Resistance[u, s] Length[u, s]];
LinearChargeDensity[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s]/Length[u, s]];
MagneticDipoleMoment[u_UnitSystem, s_UnitSystem] := Unit[Current[u, s] LorentzConstant[u, s] Area[u, s]];
Mobility[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s] Time[u, s]/Mass[u, s]];
Reluctance[u_UnitSystem, s_UnitSystem] := Unit[RationalizationConstant[u, s] LorentzConstant[u, s]^2/Inductance[u, s]];
VectorPotential[u_UnitSystem, s_UnitSystem] := Unit[MagneticFlux[u, s]/Length[u, s]];
MagneticMoment[u_UnitSystem, s_UnitSystem] := Unit[MagneticFlux[u, s] Length[u, s]];
Rigidity[u_UnitSystem, s_UnitSystem] := Unit[MagneticFluxDensity[u, s] Length[u, s]];
Susceptibility[u_UnitSystem, s_UnitSystem] := Unit[RationalizationConstant[s, u]];

(* WARNING unchecked: rigitidy, magneticmoment, vectorpotential, \
mobility, linearchargedensity, exposure *)

ElectricFlux[u_UnitSystem, s_UnitSystem] := Unit[ElectricPotential[u, s] Length[u, s]];
ElectricDipoleMoment[u_UnitSystem, s_UnitSystem] := Unit[Charge[u, s] Length[u, s]];
MagneticPotential[u_UnitSystem, s_UnitSystem] := Unit[MagneticFlux[u, s] Reluctance[u, s]];
PoleStrength[u_UnitSystem, s_UnitSystem] := Unit[MagneticDipoleMoment[u, s]/Length[u, s]];
Permeance[u_UnitSystem, s_UnitSystem] := Reluctance[s, u];
SpecificSusceptibility[u_UnitSystem, s_UnitSystem] := Unit[MagneticDipoleMoment[u, s]/MagneticField[u, s]/Mass[u, s]];
Magnetizability[u_UnitSystem, s_UnitSystem] := Unit[MagneticMoment[u, s]/MagneticFluxDensity[u, s]];
ElectricPolarizability[u_UnitSystem, s_UnitSystem] := Unit[ElectricDipoleMoment[u, s]/ElectricField[u, s]];
MagneticPolarizability[u_UnitSystem, s_UnitSystem] := Unit[MagneticDipoleMoment[u, s]/MagneticField[u, s]];
Magnetization[u_UnitSystem, s_UnitSystem] := Unit[MagneticMoment[u, s]/Volume[u, s]];

SpecificMagnetization[u_UnitSystem, s_UnitSystem] := Unit[MagneticMoment[s, u]/Mass[s, u]];
DemagnetizingFactor[u_UnitSystem, s_UnitSystem] := Unit[Rationalization[u, s]];

(* Thermodynamic *)

Moles[n_] := Moles[n, "Metric"];
Moles[n_, u_UnitSystem] := n/AvogadroConstant[u];
Moles[n_, u_String] := Moles[n, UnitSystem[u]];
Molecules[n_] := Molecules[n, "Metric"];
Molecules[n_, u_UnitSystem] := n AvogadroConstant[u];
Molecules[n_, u_String] := Molecules[n, UnitSystem[u]];

Temperature[u_UnitSystem, s_UnitSystem] :=
  Unit[(BoltzmannConstant[u] ElectronMass[s] SpeedOfLight[s]^2)/(BoltzmannConstant[
       s] ElectronMass[u] SpeedOfLight[u]^2)];
UnitSystem /: Entropy[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/Temperature[u, s]];
SpecificEntropy[u_UnitSystem, s_UnitSystem] := Unit[SpecificEnergy[u, s]/Temperature[u, s]];
VolumeHeatCapacity[u_UnitSystem, s_UnitSystem] := Unit[Entropy[u, s]/Volume[u, s]];
ThermalConductivity[u_UnitSystem, s_UnitSystem] := Unit[Force[u, s]/Time[u, s]/Temperature[u, s]];
ThermalConductance[u_UnitSystem, s_UnitSystem] := Unit[ThermalConductivity[u, s] Length[u, s]];
ThermalResistance[u_UnitSystem, s_UnitSystem] := ThermalConductance[s, u];
ThermalExpansion[u_UnitSystem, s_UnitSystem] := Temperature[s, u];
LapseRate[u_UnitSystem, s_UnitSystem] := Unit[Temperature[u, s]/Length[u, s]];

Molality[u_UnitSystem, s_UnitSystem] := MolarMassConstant[s, u];
Mole[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] Molality[u, s]];
Molarity[u_UnitSystem, s_UnitSystem] := Unit[Mole[u, s]/Volume[u, s]];
MolarVolume[u_UnitSystem, s_UnitSystem] := Unit[Volume[u, s]/Mole[u, s]];
MolarEntropy[u_UnitSystem, s_UnitSystem] := Unit[Entropy[u, s]/Mole[u, s]];
MolarEnergy[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/Mole[u, s]];
MolarConductivity[u_UnitSystem, s_UnitSystem] := Unit[Conductivity[u, s] Area[u, s]/Mole[u, s]];
MolarSusceptibility[u_UnitSystem, s_UnitSystem] := Unit[SpecificSusceptibility[u, s] MolarMass[u, s]];
Catalysis[u_UnitSystem, s_UnitSystem] := Unit[Mole[u, s]/Time[u, s]];
Specificity[u_UnitSystem, s_UnitSystem] := Unit[Volume[u, s]/Mole[u, s]/Time[u, s]];

LuminousFlux[u_UnitSystem, s_UnitSystem] := Unit[Frequency[u, s]^2 (MonochromaticRadiation540THzLuminousEfficacy[s] ReducedPlanckConstant[s])/(MonochromaticRadiation540THzLuminousEfficacy[u] ReducedPlanckConstant[u])];
Luminance[u_UnitSystem, s_UnitSystem] := Unit[LuminousFlux[u, s]/Area[u, s]];
LuminousEnergy[u_UnitSystem, s_UnitSystem] := Unit[Frequency[u, s] (MonochromaticRadiation540THzLuminousEfficacy[s] ReducedPlanckConstant[s])/(MonochromaticRadiation540THzLuminousEfficacy[u] ReducedPlanckConstant[u])];
LuminousExposure[u_UnitSystem, s_UnitSystem] := Unit[Luminance[u, s] Time[u, s]];

(* Physics *)

Cesium133HyperfineSplittingFrequency[u_UnitSystem] := Frequency[UnitData["\[CapitalDelta]\[Nu]Cs"], u];
HubbleParameter[u_UnitSystem] := Time[u, MatchSystem[u,"Hubble"]]

Map[(#[u_UnitSystem] := #[u, Universe[u]]) &, {CosmologicalConstant, AvogadroConstant, AtomicMassConstant, ProtonMass, EinsteinConstantSpeedOfLightToTheFourth, MolarGasConstant, StefanBoltzmannConstant, RadiationConstant, ElectricConstant, CoulombConstant, BiotSavartConstant, VacuumImpedance, FaradayConstant, JosephsonConstant, MagneticFluxQuantum, VonKlitzingConstant, ConductanceQuantum, HartreeEnergy, RydbergConstant, BohrRadius, RelativisticBohrRadius, ClassicalElectronRadius, BohrMagneton}];
CosmologicalConstant[u_UnitSystem, c_Coupling] := 3 UniverseDarkEnergyMassDensity[c] (HubbleParameter[u]/SpeedOfLight[u, c])^2
AvogadroConstant[u_UnitSystem, c_Coupling] := MolarMassConstant[u, c] ElectronRelativeAtomicMass[c]/ElectronMass[u, c];
AtomicMassConstant[u_UnitSystem, c_Coupling] := ElectronMass[u, c]/ElectronRelativeAtomicMass[c];
ProtonMass[u_UnitSystem, c_Coupling] := ProtonElectronMassRatio[c] ElectronMass[u, c];
EinsteinConstantSpeedOfLightToTheFourth[u_UnitSystem, c_Coupling] := (8 Pi GravitationalConstant[u, c])/SpeedOfLight[u, c]^4;
MolarGasConstant[u_UnitSystem, c_Coupling] := BoltzmannConstant[u, c] AvogadroConstant[u, c];
StefanBoltzmannConstant[u_UnitSystem, c_Coupling] := (2 Pi^5 BoltzmannConstant[u, c]^4)/(15 PlanckConstant[u, c]^3 SpeedOfLight[u, c]^2);
RadiationConstant[u_UnitSystem, c_Coupling] := (4 StefanBoltzmannConstant[u, c])/SpeedOfLight[u, c];
ElectricConstant[u_UnitSystem, c_Coupling] := (1/MagneticConstant[u, c]) (SpeedOfLight[u, c] LorentzConstant[u]^2);
CoulombConstant[u_UnitSystem, c_Coupling] := RationalizationConstant[u]/(4 Pi)/ ElectricConstant[u];
BiotSavartConstant[u_UnitSystem, c_Coupling] := MagneticConstant[u, c] LorentzConstant[u] (RationalizationConstant[u]/(4 Pi));
AmpereConstant[u_UnitSystem] := LorentzConstant[u] BiotSavartConstant[u];
VacuumImpedance[u_UnitSystem, c_Coupling] := MagneticConstant[u, c] SpeedOfLight[u, c] RationalizationConstant[u] LorentzConstant[u]^2;
FaradayConstant[u_UnitSystem, c_Coupling] := ElementaryCharge[u, c] AvogadroConstant[u, c];
JosephsonConstant[u_UnitSystem, c_Coupling] := 2 ElementaryCharge[u, c] LorentzConstant[u]/PlanckConstant[u, c];
MagneticFluxQuantum[u_UnitSystem, c_Coupling] := 1/JosephsonConstant[u, c];
VonKlitzingConstant[u_UnitSystem, c_Coupling] := PlanckConstant[u, c]/ElementaryCharge[u, c]^2;
ConductanceQuantum[u_UnitSystem, c_Coupling] := (2 ElementaryCharge[u, c]^2)/PlanckConstant[u, c];
HartreeEnergy[u_UnitSystem, c_Coupling] := ElectronMass[u, c] (SpeedOfLight[u, c] FineStructureConstant[c])^2;
RydbergConstant[u_UnitSystem, c_Coupling] := HartreeEnergy[u, c]/(2 PlanckConstant[u, c])/SpeedOfLight[u, c];
BohrRadius[u_UnitSystem, c_Coupling] := ReducedPlanckConstant[u, c]/ElectronMass[u, c]/SpeedOfLight[u, c]/FineStructureConstant[c];
RelativisticBohrRadius[u_UnitSystem, c_Coupling] := BohrRadius[u, c] (1 + 1/ProtonElectronMassRatio[c]);
ClassicalElectronRadius[u_UnitSystem, c_Coupling] := FineStructureConstant[c] ReducedPlanckConstant[u, c]/ElectronMass[u, c]/SpeedOfLight[u, c];
BohrMagneton[u_UnitSystem, c_Coupling] := ElementaryCharge[u, c] ReducedPlanckConstant[u, c] LorentzConstant[u]/2 ElectronMass[u, c];
