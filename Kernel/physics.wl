
(*   This file is part of UnitSystems                   *)
(*   It is licensed under the MIT license               *)
(*   UnitSystems Copyright (C) 2022 Michael Reed        *)
(*       _           _                         _        *)
(*      | |         | |                       | |       *)
(*   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _  *)
(*  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` | *)
(* | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| | *)
(*  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_| *)

UnitSystem /: SolidAngle[u_UnitSystem,s_UnitSystem] := Unit[Angle[u,s]^2]
Spat[u_UnitSystem] := Unit[2 Turn[u]/Angle[u] (Turn[u]/MM[Turn[u]])]

ElectronMass[h_?NumberQ] := UnitData["R\[Infinity]"] 2 h/UnitData["\[Alpha]"]^2/UnitData["c"];
ElectronMass[h_?NumberQ, u_Coupling] := (UnitData["R\[Infinity]"] 2 h)/(FineStructureConstant[u]^2 UnitData["c"]);
PlanckMass[u_UnitSystem, c_Coupling] := PowerExpand[ElectronMass[u, c]/Sqrt[GravitationalCouplingConstantElectronElectron[c]]]
PlanckConstant[u_UnitSystem, c_Coupling] := Turn[u] ReducedPlanckConstant[u];
GravitationalConstant[u_UnitSystem, c_Coupling] := PowerExpand[(SpeedOfLight[u, c] PlanckConstant[u, c])/(2 Pi PlanckMass[u, c]^2)]
ElementaryCharge[u_UnitSystem, c_Coupling] := PowerExpand[Sqrt[2 PlanckConstant[u]/(MagneticConstant[u]/FineStructureConstant[c])/(SpeedOfLight[u] RationalizationConstant[u] LorentzConstant[u]^2)]]

ElectronMass[u:UnitSystem["Planck"], c_Coupling] := PowerExpand[Sqrt[Spat[u] GravitationalCouplingConstantElectronElectron[c]]]
ElectronMass[UnitSystem["PlanckGauss"], c_Coupling] := Sqrt[GravitationalCouplingConstantElectronElectron[c]];
ElectronMass[UnitSystem[_, _, _, _, Sqrt[EvalUnitData["\[Alpha]G"]/UnitData["\[Alpha]"]], ___], c_Coupling] := PowerExpand[Sqrt[GravitationalCouplingConstantElectronElectron[c]/FineStructureConstant[c]]]
ElectronMass[UnitSystem[_, _, _, _, 1/UnitData["\[Mu]pe"], ___], c_Coupling] := 1/ProtonElectronMassRatio[c];
MagneticConstant[u:UnitSystem[_, _, _, 4 Pi UnitData["\[Alpha]"]^2, ___], c_Coupling] := Spat[u] FineStructureConstant[c]^2;
MagneticConstant[u:UnitSystem[_, _, _, Pi UnitData["\[Alpha]"]^2, ___], c_Coupling] := Spat[u]/4 FineStructureConstant[c]^2;
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

(*(* IAU to SI *)
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
   s : UnitSystem[_, _, UnitData["c"]/UnitData["ft"], ___]] := Time[u, s, UnitData["day"]];*)

(*UnitSystem /: Length[u_UnitSystem, s_UnitSystem] := Length[u, s, 1];*)
UnitSystem /: Length[u_UnitSystem, s_UnitSystem] :=
  Unit[measuratio[Turn[s] ReducedPlanckConstant[s] GravityConstant[s]/(ElectronMass[s] SpeedOfLight[s]),Turn[u] ReducedPlanckConstant[u] GravityConstant[u]/(ElectronMass[u] SpeedOfLight[u])]];
UnitSystem /: Area[u_UnitSystem, s_UnitSystem] := Unit[Length[u, s]^2];
UnitSystem /: Volume[u_UnitSystem, s_UnitSystem] := Unit[Length[u, s]^3];

Force = dF
Mass = dM
Time = dT

AngularLength = dL*dA
AngularArea = dA^2/dA^2
Wavenumber = d1/dL
AngularWavenumber = dA/dL
Etendue = dL^2*dA^2
FuelEconomy = dL^2
NumberDensity = d1/dL^3
AngularTime = dT/dA
Frequency = d1/Time
AngularFrequency = dA/Time
FrequencyDrift = Frequency^2
Speed = dL/Time
Stagnance = Time/dL
PhotonIntensity = Frequency/dA^2
PhotonIrradiance = Wavenumber/Speed
PhotonRadiance = PhotonIrradiance/dA^2
Acceleration = Speed/Time
Jerk = Speed/Time^2
Snap = Speed/Time^3
Crackle = Speed/Time^4
Pop = Speed/Time^5
VolumeFlowRate = dL^2 Speed

GravityForce = Mass*dL/Time^2/Force
SpecificEnergy = Speed^2/GravityForce
Inertia = Mass/GravityForce
Energy = Mass*SpecificEnergy
SpecificForce = Acceleration/GravityForce
Pressure = Force/dL^2

(*Wavenumber[u_UnitSystem, s_UnitSystem] := Unit[Length[s, u]];
AngularWavenumber[u_UnitSystem, s_UnitSystem] := Unit[Angle[u,s] Length[s, u]];
FuelEconomy[u_UnitSystem, s_UnitSystem] := Area[s, u];
NumberDensity[u_UnitSystem, s_UnitSystem] := Volume[s, u];
Time[u_UnitSystem, s_UnitSystem] := Time[u, s, 1];
Time[u_UnitSystem, s_UnitSystem, t_] := Unit[Length[u, s]/SpeedOfLight[u, s], 1];
Frequency[u_UnitSystem, s_UnitSystem] := Time[s, u];
AngularFrequency[u_UnitSystem, s_UnitSystem] := Unit[Angle[u, s] Time[s, u]];
FrequencyDrift[u_UnitSystem, s_UnitSystem] := Unit[Time[s, u]^2];
Speed[u_UnitSystem, s_UnitSystem] := SpeedOfLight[u, s];
Acceleration[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]/Time[u, s]];
Jerk[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]/Time[u, s]^2];
Snap[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]/Time[u, s]^3];
Crackle[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]/Time[u, s]^4];
Pop[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]/Time[u, s]^5];
VolumeFlowRate[u_UnitSystem, s_UnitSystem] := Unit[Area[u, s] Speed[u, s]];
SpecificEnergy[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s]^2/GravityConstant[u, s]];*)

(*Inertia[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/GravityConstant[u, s]];
Mass[u_UnitSystem, s_UnitSystem] := ElectronMass[u, s];
Energy[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] SpecificEnergy[u, s]];*)
UnitSystem /: Power[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/Time[u, s]];
(*Force[u_UnitSystem, s_UnitSystem] := Unit[Inertia[u, s] Acceleration[u, s]];
SpecificForce[u_UnitSystem, s_UnitSystem] := Unit[Acceleration[u, s]/GravityConstant[u, s]];
GravityForce[u_UnitSystem, s_UnitSystem] := Unit[GravityConstant[u, s]];
Pressure[u_UnitSystem, s_UnitSystem] := Unit[Force[u, s]/Area[u, s]];*)

Impulse = Force*Time
Momentum = Mass*Speed
AngularMomentum = Momentum*dL*dA
ForceOnsetRate = Mass*Jerk
MassPerArea = Mass/dL^2
MassDensity = Mass/dL^3
SpecificWeight = Force/dL^3
SpecificVolume = dL^3/Mass
Action = Energy*Time
Irradiance = Energy/Time/dL^2
Radiance = Irradiance/dA^2
RadiantIntensity = Energy/Time/dA^2
SpectralExposure = Force/Speed
KinematicViscosity = Speed*dL
DynamicViscosity = Force/Speed/dL
LinearMassDensity = Mass/dL
MassFlowRate = Mass/Time
PowerGradient = Energy/Time/dL
PowerDensity = PowerGradient/dL^2
Compressibility = d1/Pressure
EnergyPerArea = Energy/dL^2

(*Impulse[u_UnitSystem, s_UnitSystem] := Unit[Force[u, s] Time[u, s]];
Momentum[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] Speed[u, s]];
AngularMomentum[u_UnitSystem, s_UnitSystem] := Unit[Momentum[u, s] Length[u, s] Angle[u, s]];
ForceOnsetRate[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] Jerk[u, s]];
MassPerArea[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/Area[u, s]];
MassDensity[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/Volume[u, s]];
SpecificWeight[u_UnitSystem, s_UnitSystem] := Unit[Force[u, s]/Volume[u, s]];
SpecificVolume[u_UnitSystem, s_UnitSystem] := Unit[Volume[u, s]/Mass[u, s]];
Action[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s] Time[u, s]];
(*Stiffness[u_UnitSystem, s_UnitSystem] := Unit[Force[u, s]/Length[u, s]];*)
Irradiance[u_UnitSystem, s_UnitSystem] := Unit[Power[u, s]/Area[u, s]];
KinematicViscosity[u_UnitSystem, s_UnitSystem] := Unit[Speed[u, s] Length[u, s]];
DynamicViscosity[u_UnitSystem, s_UnitSystem] := Unit[Force[u, s]/Speed[u, s]/Length[u, s]];
LinearMassDensity[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/Length[u, s]];
MassFlowRate[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/Time[u, s]];
PowerGradient[u_UnitSystem, s_UnitSystem] := Unit[Power[u, s]/Length[u, s]];
PowerDensity[u_UnitSystem, s_UnitSystem] := Unit[Power[u, s]/Volume[u, s]];
Compressibility[u_UnitSystem, s_UnitSystem] := Pressure[s, u];
EnergyPerArea[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/Area[u, s]];*)
UnitSystem /: MomentOfInertia[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] Area[u, s]];

SoundExposure = Time*Pressure^2
SpecificAcousticImpedance = Pressure/Speed
AcousticImpedance = SpecificAcousticImpedance/dL^2
Admittance = dL^2/SpecificAcousticImpedance
Compliance = Time^2/Mass
Inertance = Mass/dL^4

(*SoundExposure[u_UnitSystem, s_UnitSystem] := Unit[Time[u, s] Pressure[u, s]^2];
SpecificAcousticImpedance[u_UnitSystem, s_UnitSystem] := Unit[Pressure[u, s]/Speed[u, s]];
AcousticImpedance[u_UnitSystem, s_UnitSystem] := Unit[SpecificAcousticImpedance[u, s]/Area[u, s]];
Admittance[u_UnitSystem, s_UnitSystem] := Unit[Area[u, s]/SpecificAcousticImpedance[u, s]];
Compliance[u_UnitSystem, s_UnitSystem] := Unit[Time[u, s]^2/Mass[u, s]];
Inertance[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s]/Length[u, s]^4];*)

(* Electromagnetic *)

ElectricCharge = dQ
(*Rationalization = dR*)
MagneticPermeability = Force*Time^2/ElectricCharge^2/Rationalization*dC^2
ElectricCurrent = ElectricCharge/Time
ElectricPotential = Energy/ElectricCharge
ElectricCapacitance = ElectricCharge/ElectricPotential
ElectricResistance = ElectricPotential/ElectricCurrent
ElectricConductance = ElectricCurrent/ElectricPotential
MagneticFlux = Energy*dC
MagneticFluxDensity = MagneticFlux/dL^2
MagneticInductance = MagneticFlux/ElectricCurrent

(*MagneticPermeability[u_UnitSystem, s_UnitSystem] := Unit[MagneticConstant[u, s]];
ElectricCharge[u_UnitSystem, s_UnitSystem] :=
  Unit[Sqrt[measuratio[Turn[s] ReducedPlanckConstant[s]/(MagneticConstant[s] SpeedOfLight[s] RationalizationConstant[s] LorentzConstant[s]^2),
	Turn[u] ReducedPlanckConstant[u]/(MagneticConstant[u] SpeedOfLight[u] RationalizationConstant[u] LorentzConstant[u]^2)]]];
ElectricCurrent[u_UnitSystem, s_UnitSystem] := Unit[ElectricCharge[u, s]/Time[u, s]];
ElectricPotential[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/ElectricCharge[u, s]];
ElectricCapacitance[u_UnitSystem, s_UnitSystem] := Unit[ElectricCharge[u, s]/ElectricPotential[u, s]];
ElectricResistance[u_UnitSystem, s_UnitSystem] := Unit[ElectricPotential[u, s]/ElectricCurrent[u, s]];
ElectricConductance[u_UnitSystem, s_UnitSystem] := Unit[ElectricCurrent[u, s]/ElectricPotential[u, s]];
MagneticFlux[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/LorentzConstant[u, s]/ElectricCurrent[u, s]];
MagneticFluxDensity[u_UnitSystem, s_UnitSystem] :=
  Unit[MagneticFlux[u, s]/Area[u, s]];
MagneticInductance[u_UnitSystem, s_UnitSystem] := Unit[MagneticFlux[u, s]/ElectricCurrent[u, s]];*)

LinearElectricChargeDensity = ElectricCharge/dL
(*ElectricDisplacement?ElectricFluxDensity = ElectricCharge*dR/dL^2*)
ElectricChargeDensity = ElectricCharge/dL^3
ElectricCurrentDensity = ElectricCurrent/dL^2
Conductivity = ElectricConductance/dL
ElectricPermittivity = ElectricCapacitance*dR/dL
ElectricFieldStrength = ElectricPotential/dL
MagneticFieldStrength = ElectricCurrent*dR/dC/dL
ElectricChargePerMass = ElectricCharge/Mass
Resistivity = ElectricResistance*dL
MagneticDipoleMoment = ElectricCurrent/dC*dL^2/GravityForce/dA
ElectricalMobility = dL*Speed*ElectricPotential
MagneticReluctance = dR/dC^2/MagneticInductance
MagneticVectorPotential = MagneticFlux/dL
MagneticMoment = MagneticFlux*dL
MagneticSusceptibility = d1/dR

(*ElectricFluxDensity[u_UnitSystem, s_UnitSystem] := Unit[ElectricCharge[u, s] RationalizationConstant[u, s]/Area[u, s]];
ElectricChargeDensity[u_UnitSystem, s_UnitSystem] := Unit[ElectricCharge[u, s]/Volume[u, s]];
ElectricCurrentDensity[u_UnitSystem, s_UnitSystem] := Unit[ElectricCurrent[u, s]/Area[u, s]];
Conductivity[u_UnitSystem, s_UnitSystem] := Unit[ElectricConductance[u, s]/Length[u, s]];
ElectricPermittivity[u_UnitSystem, s_UnitSystem] := Unit[ElectricCapacitance[u, s] RationalizationConstant[u, s]/Length[u, s]];
ElectricFieldStrength[u_UnitSystem, s_UnitSystem] := Unit[ElectricPotential[u, s]/Length[u, s]];
MagneticFieldStrength[u_UnitSystem, s_UnitSystem] := Unit[ElectricCurrent[u, s] RationalizationConstant[u, s] LorentzConstant[u, s]/Length[u, s]];
ElectricChargePerMass[u_UnitSystem, s_UnitSystem] := Unit[ElectricCharge[u, s]/Mass[u, s]];
Resistivity[u_UnitSystem, s_UnitSystem] := Unit[ElectricResistance[u, s] Length[u, s]];
LinearElectricChargeDensity[u_UnitSystem, s_UnitSystem] := Unit[ElectricCharge[u, s]/Length[u, s]];
MagneticDipoleMoment[u_UnitSystem, s_UnitSystem] := Unit[ElectricCurrent[u, s] LorentzConstant[u, s] Area[u, s]/GravityCsontant[u, s]/Angle[u, s]];
ElectricalMobility[u_UnitSystem, s_UnitSystem] := Unit[Length[u, s] Speed[u, s] ElectricPotential[u, s]];
MagneticReluctance[u_UnitSystem, s_UnitSystem] := Unit[RationalizationConstant[u, s] LorentzConstant[u, s]^2/MagneticInductance[u, s]];
MagneticVectorPotential[u_UnitSystem, s_UnitSystem] := Unit[MagneticFlux[u, s]/Length[u, s]];
MagneticMoment[u_UnitSystem, s_UnitSystem] := Unit[MagneticFlux[u, s] Length[u, s]];
(*MagneticRigidity[u_UnitSystem, s_UnitSystem] := Unit[MagneticFluxDensity[u, s] Length[u, s]];*)
MagneticSusceptibility[u_UnitSystem, s_UnitSystem] := Unit[RationalizationConstant[s, u]];*)

(* WARNING unchecked: rigitidy, magneticmoment, vectorpotential, \
mobility, linearchargedensity, exposure *)

ElectricFlux = ElectricPotential*dL
ElectricDipoleMoment = ElectricCharge*dL
MagnetomotiveForce = MagneticFlux*MagneticReluctance
MagneticPoleStrength = MagneticDipoleMoment/dL
MagneticPermeance = d1/MagneticReluctance
SpecificSusceptibility = MagneticDipoleMoment/MagneticFieldStrength/Mass
ElectricPolarizability = ElectricDipoleMoment/ElectricFieldStrength
MagneticPolarizability = MagneticDipoleMoment/MagneticFieldStrength

SpecificMagnetization = MagneticMoment/Mass
DemagnetizingFactor = dR

(*ElectricFlux[u_UnitSystem, s_UnitSystem] := Unit[ElectricPotential[u, s] Length[u, s]];
ElectricDipoleMoment[u_UnitSystem, s_UnitSystem] := Unit[ElectricCharge[u, s] Length[u, s]];
MagnetomotiveForce[u_UnitSystem, s_UnitSystem] := Unit[MagneticFlux[u, s] MagneticReluctance[u, s]];
MagneticPoleStrength[u_UnitSystem, s_UnitSystem] := Unit[MagneticDipoleMoment[u, s]/Length[u, s]];
MagneticPermeance[u_UnitSystem, s_UnitSystem] := MagneticReluctance[s, u];
SpecificSusceptibility[u_UnitSystem, s_UnitSystem] := Unit[MagneticDipoleMoment[u, s]/MagneticFieldStrength[u, s]/Mass[u, s]];
(*Magnetizability[u_UnitSystem, s_UnitSystem] := Unit[MagneticMoment[u, s]/MagneticFluxDensity[u, s]];*)
ElectricPolarizability[u_UnitSystem, s_UnitSystem] := Unit[ElectricDipoleMoment[u, s]/ElectricFieldStrength[u, s]];
MagneticPolarizability[u_UnitSystem, s_UnitSystem] := Unit[MagneticDipoleMoment[u, s]/MagneticFieldStrength[u, s]];
(*Magnetization[u_UnitSystem, s_UnitSystem] := Unit[MagneticMoment[u, s]/Volume[u, s]];*)

SpecificMagnetization[u_UnitSystem, s_UnitSystem] := Unit[MagneticMoment[s, u]/Mass[s, u]];
DemagnetizationFactor[u_UnitSystem, s_UnitSystem] := Unit[RationalizationConstant[u, s]];*)

(* Thermodynamic *)

Moles[n_] := Moles[n, "Metric"];
Moles[n_, u_UnitSystem] := n/AvogadroConstant[u];
Moles[n_, u_String] := Moles[n, UnitSystem[u]];
Molecules[n_] := Molecules[n, "Metric"];
Molecules[n_, u_UnitSystem] := n AvogadroConstant[u];
Molecules[n_, u_String] := Molecules[n, UnitSystem[u]];

Temperature = d0
SpecificEntropy = SpecificEnergy/Temperature
EntropyPerVolume = Energy/Temperature/dL^3
ThermalConductivity = Force/Time/Temperature
ThermalConductance = ThermalConductivity*dL
ThermalResistivity = d1/ThermalConductivity
ThermalResistance = d1/ThermalConductance
ThermalExpansion = d1/Temperature
TemperatureGradient = Temperature/dL

(*Temperature[u_UnitSystem, s_UnitSystem] :=
  Unit[measuratio[ElectronMass[s] SpeedOfLight[s]^2/BoltzmannConstant[s]/GravityConstant[s],
	ElectronMass[u] SpeedOfLight[u]^2/BoltzmannConstant[u]/GravityConstant[u]]];*)
UnitSystem /: Entropy[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/Temperature[u, s]];
(*SpecificEntropy[u_UnitSystem, s_UnitSystem] := Unit[SpecificEnergy[u, s]/Temperature[u, s]];
EntropyPerVolume[u_UnitSystem, s_UnitSystem] := Unit[Entropy[u, s]/Volume[u, s]];
ThermalConductivity[u_UnitSystem, s_UnitSystem] := Unit[Force[u, s]/Time[u, s]/Temperature[u, s]];
ThermalConductance[u_UnitSystem, s_UnitSystem] := Unit[ThermalConductivity[u, s] Length[u, s]];
ThermalResistivity[u_UnitSystem, s_UnitSystem] := ThermalConductivity[s, u];
ThermalResistance[u_UnitSystem, s_UnitSystem] := ThermalConductance[s, u];
ThermalExpansion[u_UnitSystem, s_UnitSystem] := Temperature[s, u];
TemperatureGradient[u_UnitSystem, s_UnitSystem] := Unit[Temperature[u, s]/Length[u, s]];*)

Amount = dN
MolarMass = Amount/Mass
Molality = Mass/Amount
AmountConcentration = Amount/dL^3
MolarVolume = dL^3/Amount
MolarEntropy = Energy/Temperature/Amount
MolarEnergy = Energy/Amount
MolarConductivity = Conductivity*dL^2/Amount
MolarMagneticSusceptibility = SpecificSusceptibility*MolarMass
Catalysis = Amount/Time
Specificity = dL^3/Amount/Time
DiffusionFlux = Amount*PhotonIrradiance

(*MolarMass[u_UnitSystem, s_UnitSystem] := Unit[MolarMassConstant[u, s]];
Molality[u_UnitSystem, s_UnitSystem] := MolarMassConstant[s, u];
Amount[u_UnitSystem, s_UnitSystem] := Unit[Mass[u, s] Molality[u, s]];
AmountConcentration[u_UnitSystem, s_UnitSystem] := Unit[Amount[u, s]/Volume[u, s]];
MolarVolume[u_UnitSystem, s_UnitSystem] := Unit[Volume[u, s]/Amount[u, s]];
MolarEntropy[u_UnitSystem, s_UnitSystem] := Unit[Entropy[u, s]/Amount[u, s]];
MolarEnergy[u_UnitSystem, s_UnitSystem] := Unit[Energy[u, s]/Amount[u, s]];
MolarConductivity[u_UnitSystem, s_UnitSystem] := Unit[Conductivity[u, s] Area[u, s]/Amount[u, s]];
MolarMagneticSusceptibility[u_UnitSystem, s_UnitSystem] := Unit[SpecificSusceptibility[u, s] MolarMass[u, s]];
Catalysis[u_UnitSystem, s_UnitSystem] := Unit[Amount[u, s]/Time[u, s]];
Specificity[u_UnitSystem, s_UnitSystem] := Unit[Volume[u, s]/Amount[u, s]/Time[u, s]];*)

LuminousFlux = dJ
LuminousEfficacyOfRadiation = Time*LuminousFlux/Force/dL
LuminousIntensity = LuminousFlux/dA^2
Illuminance = LuminousFlux/dL^2
Luminance = LuminousIntensity/dL^2
LuminousEnergy = Time*LuminousFlux
LuminousExposure = Luminance*Time

(*LuminousEfficacyOfRadiation[u_UnitSystem, s_UnitSystem] := Unit[MonochromaticRadiation540THzLuminousEfficacy[u, s]];
LuminousFlux[u_UnitSystem, s_UnitSystem] := Unit[Frequency[u, s] LuminousEnergy[u, s]];
LuminousIntensity[u_UnitSystem, s_UnitSystem] := Unit[LuminousFlux[u, s]/SolidAngle[u, s]];
Illuminance[u_UnitSystem, s_UnitSystem] := Unit[LuminousFlux[u, s]/Area[u, s]];
Luminance[u_UnitSystem, s_UnitSystem] := Unit[LuminousIntensity[u, s]/Area[u, s]];
LuminousEnergy[u_UnitSystem, s_UnitSystem] := Unit[Frequency[u, s] measuratio[MonochromaticRadiation540THzLuminousEfficacy[s] ReducedPlanckConstant[s], MonochromaticRadiation540THzLuminousEfficacy[u] ReducedPlanckConstant[u]]];
LuminousExposure[u_UnitSystem, s_UnitSystem] := Unit[Luminance[u, s] Time[u, s]];*)

(* Physics *)

Cesium133HyperfineSplittingFrequency[u : UnitSystem[_?NumberQ, ___]] := Frequency[If[AbstractUnitSystemQ[u], "\[CapitalDelta]\[Nu]Cs", UnitData["\[CapitalDelta]\[Nu]Cs"]], u];
Cesium133HyperfineSplittingFrequency[u : UnitSystem[_Around, ___]] := Frequency[UnitData["\[CapitalDelta]\[Nu]Cs"], u];
Cesium133HyperfineSplittingFrequency[u_UnitSystem] := Frequency["\[CapitalDelta]\[Nu]Cs", u];
HubbleParameter[u_UnitSystem] := Time[u, MatchSystem[u,"Hubble"]]
CosmologicalConstant[u_UnitSystem, c_Coupling] := 3 UniverseDarkEnergyMassDensity[c] (HubbleParameter[u]/SpeedOfLight[u, c])^2

Map[(#[u_UnitSystem] := #[u, Universe[u]]) &, {CosmologicalConstant, AvogadroConstant, AtomicMassConstant, ProtonMass, EinsteinConstantSpeedOfLightToTheFourth, MolarGasConstant, StefanBoltzmannConstant, RadiationConstant, ElectricConstant, CoulombConstant, BiotSavartConstant, VacuumImpedance, FaradayConstant, JosephsonConstant, MagneticFluxQuantum, VonKlitzingConstant, ConductanceQuantum, HartreeEnergy, RydbergConstant, BohrRadius, RelativisticBohrRadius, ClassicalElectronRadius, BohrMagneton}];
AvogadroConstant[u_UnitSystem, c_Coupling] := MolarMassConstant[u, c] ElectronRelativeAtomicMass[c]/ElectronMass[u, c];
AtomicMassConstant[u_UnitSystem, c_Coupling] := ElectronMass[u, c]/ElectronRelativeAtomicMass[c];
ProtonMass[u_UnitSystem, c_Coupling] := ProtonElectronMassRatio[c] ElectronMass[u, c];
GaussianGravitationalConstant[u_UnitSystem, c_Couling] = Measure[MM@Sqrt[GravitationalConstant[u,c]],dA/dT,"IAU"][u]
EinsteinConstantSpeedOfLightSquared[u_UnitSystem, c_Coupling] := (2 Spat[u] GravitationalConstant[u, c])/SpeedOfLight[u, c]^2;
EinsteinConstantSpeedOfLightToTheFourth[u_UnitSystem, c_Coupling] := (2 Spat[u] GravitationalConstant[u, c])/SpeedOfLight[u, c]^4;
MolarGasConstant[u_UnitSystem, c_Coupling] := BoltzmannConstant[u, c] AvogadroConstant[u, c];
StefanBoltzmannConstant[u_UnitSystem, c_Coupling] := (MM[Turn[u]]^4/2^5 Spat[u] BoltzmannConstant[u, c]^4)/(15 PlanckConstant[u, c]^3 SpeedOfLight[u, c]^2);
RadiationConstant[u_UnitSystem, c_Coupling] := (4 StefanBoltzmannConstant[u, c])/SpeedOfLight[u, c];
ElectricConstant[u_UnitSystem, c_Coupling] := 1/(MagneticConstant[u, c] (SpeedOfLight[u, c] LorentzConstant[u])^2);
CoulombConstant[u_UnitSystem, c_Coupling] := RationalizationConstant[u]/Spat[u]/ ElectricConstant[u];
BiotSavartConstant[u_UnitSystem, c_Coupling] := MagneticConstant[u, c] LorentzConstant[u] (RationalizationConstant[u]/Spat[u]);
AmpereConstant[u_UnitSystem] := LorentzConstant[u] BiotSavartConstant[u];
VacuumImpedance[u_UnitSystem, c_Coupling] := MagneticConstant[u, c] SpeedOfLight[u, c] RationalizationConstant[u] LorentzConstant[u]^2;
FaradayConstant[u_UnitSystem, c_Coupling] := ElementaryCharge[u, c] AvogadroConstant[u, c];
JosephsonConstant[u_UnitSystem, c_Coupling] := 2 ElementaryCharge[u, c] LorentzConstant[u]/PlanckConstant[u, c];
MagneticFluxQuantum[u_UnitSystem, c_Coupling] := 1/JosephsonConstant[u, c];
VonKlitzingConstant[u_UnitSystem, c_Coupling] := PlanckConstant[u, c]/ElementaryCharge[u, c]^2;
ConductanceQuantum[u_UnitSystem, c_Coupling] := (2 ElementaryCharge[u, c]^2)/PlanckConstant[u, c];
HartreeEnergy[u_UnitSystem, c_Coupling] := ElectronMass[u, c]/GravityConstant[u] (SpeedOfLight[u, c] FineStructureConstant[c])^2;
RydbergConstant[u_UnitSystem, c_Coupling] := HartreeEnergy[u, c]/(2 PlanckConstant[u, c])/SpeedOfLight[u, c];
BohrRadius[u_UnitSystem, c_Coupling] := ReducedPlanckConstant[u, c] GravityConstant[u]/ElectronMass[u, c]/SpeedOfLight[u, c]/FineStructureConstant[c];
ClassicalElectronRadius[u_UnitSystem, c_Coupling] := FineStructureConstant[c] ReducedPlanckConstant[u, c]/ElectronMass[u, c]/SpeedOfLight[u, c];
BohrMagneton[u_UnitSystem, c_Coupling] := ReducedPlanckConstant[u, c] GravityConstant[u]/ElectronMass[u, c]/SpeedOfLight[u]/FineStructureConstant[c];
