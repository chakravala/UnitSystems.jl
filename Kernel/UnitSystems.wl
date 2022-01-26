(* ::Package:: *)
(* This file is part of UnitSystems. It is licensed under the MIT license *)
(* UnitSystems Copyright (C) 2021 Michael Reed *)

Unprotect[UnitSystem];
ProtectedList = {Length, Area, Volume, Power, Entropy};
UnitSystemsList = {"Metric", "SI2019", "CODATA", "Conventional", "MTS", "English",
    "EnglishUS", "FFF", "IAU", "SI1976", "Mixed", "ESU2019", "EMU2019", "EMU", "ESU",
	"Gauss", "LorentzHeaviside", "Thomson", "Kennelly", "Planck", "PlanckGauss", "Stoney",
	"Hartree", "Rydberg", "Schrodinger", "Electronic", "Natural", "NaturalGauss",
	"QCD", "QCDGauss", "QCDoriginal", "Hubble", "Cosmological", "CosmologicalQuantum"};
ConstantsList = {Cesium133HyperfineSplittingFrequency, SpeedOfLight,
	PlanckConstant, ReducedPlanckConstant, ElectronMass,
	MolarMassConstant, BoltzmannConstant, MagneticConstant, RationalizationConstant,
	LorentzConstant, MonochromaticRadiation540THzLuminousEfficacy};
PhysicsList = {AtomicMassConstant, ProtonMass, PlanckMass,
	GravitationalConstant, EinsteinConstantSpeedOfLightToTheFourth,
	HartreeEnergy, RydbergConstant, BohrRadius,
	RelativisticBohrRadius, ClassicalElectronRadius,
	AvogadroConstant, MolarGasConstant, StefanBoltzmannConstant,
	RadiationConstant, ElectricConstant, CoulombConstant, AmpereConstant,
	BiotSavartConstant, ElementaryCharge, FaradayConstant, VacuumImpedance,
	ConductanceQuantum, VonKlitzingConstant, JosephsonConstant,
	MagneticFluxQuantum, BohrMagneton};
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
    Permittivity, Permeability, Susceptibility, Current, Conductance, 
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
GravitationalCouplingConstantElectronElectron[Coupling[\[Alpha]G_, ___]] := measure[\[Alpha]G];
FineStructureConstant[Coupling[_, \[Alpha]_, ___]] := measure[\[Alpha]];
ElectronRelativeAtomicMass[Coupling[_, _, \[Mu]eu_, ___]] := measure[\[Mu]eu];
ProtonRelativeAtomicMass[Coupling[_, _, _, \[Mu]pu_, ___]] := measure[\[Mu]pu];
ProtonElectronMassRatio[c_Coupling] := ProtonRelativeAtomicMass[c]/ElectronRelativeAtomicMass[c];

BoltzmannConstant[UnitSystem[kB_, ___]] := kB;
ReducedPlanckConstant[UnitSystem[_, \[HBar]_, ___]] := \[HBar];
SpeedOfLight[UnitSystem[_, _, c_, ___]] := c;
MagneticConstant[UnitSystem[_, _, _, \[Mu]0_, ___]] := \[Mu]0;
ElectronMass[UnitSystem[_, _, _, _, me_, ___]] := me;
RationalizationConstant[UnitSystem[_, _, _, _, _, \[Lambda]_, ___]] := \[Lambda];
LorentzConstant[UnitSystem[_, _, _, _, _, _, \[Alpha]L_, ___]] := \[Alpha]L;

LorentzConstant[_UnitSystem] := 1
RationalizationConstant[_UnitSystem] := 1
RationalizedQ[u_UnitSystem] := RationalizationConstant[u] != 4 Pi
RationalizedQ[u_String] := RationalizedQ[UnitSystem[u]]

Universe[UnitSystem[_?NumericQ, ___]] := StandardModel
Universe[UnitSystem[_Around, ___]] := StandardModel
Unit[x_, y_ : 1] := PowerExpand[x];

UnitData = <||>

AppendTo[UnitData, "g0" -> 980665/10^5]
AppendTo[UnitData, "ft" -> 3048/10000]
AppendTo[UnitData, "ftUS" -> 1200/3937]
AppendTo[UnitData, "lb" -> 45359237/10^8]
AppendTo[UnitData, "rankine" -> 5/9]
AppendTo[UnitData, "\[CapitalDelta]\[Nu]Cs" -> 9192631770]
AppendTo[UnitData, "Kcd" -> 683 555016/555000]
AppendTo[UnitData, "mP" -> Around[2.176434/10^8, 24/10^14]]
AppendTo[UnitData, "NA" -> 602214076 10^15]
AppendTo[UnitData, "kB" -> 1380649/10^29]
AppendTo[UnitData, "h" -> 662607015/10^42]
AppendTo[UnitData, "c" -> 299792458]
AppendTo[UnitData, "e" -> 1602176634/10^28]
AppendTo[UnitData, "\[Mu]eu" -> 1/Around[1822.888486209, 53/10^9]]
AppendTo[UnitData, "\[Mu]pu" -> Around[1.007276466621, 53/10^12]]
AppendTo[UnitData, "\[Alpha]" -> 1/Around[137.035999084, 21/10^9]]
AppendTo[UnitData, "R\[Infinity]" -> Around[10973731.5681601, 21/10^6]]
AppendTo[UnitData, "RK1990" -> 25812807/1000]
AppendTo[UnitData, "RK2014" -> Around[25812.8074555, 59/10^7]]
AppendTo[UnitData, "KJ1990" -> 4835979 10^8]
AppendTo[UnitData, "KJ2014" -> Around[4.835978525 10^14, 3 10^6]]
AppendTo[UnitData, "GMsun" -> Around[1.32712442099 10^20, 9 10^9]]
AppendTo[UnitData, "day" -> 60^2 24]
AppendTo[UnitData, "au" -> 149597870700]
AppendTo[UnitData, "fur" -> 201168/1000]
AppendTo[UnitData, "H0" -> Around[67.66, 0.42]]
AppendTo[UnitData, "\[CapitalOmega]\[CapitalLambda]" -> Around[0.6889, 0.0056]]

AbstractUnitData["slug"] = "lb" "g0"/"ft"
AbstractUnitData["slugUS"] = "lb" "g0"/"ftUS"
AbstractUnitData["Ru"] = "NA" "kB"
AbstractUnitData["me"] = 2 "R\[Infinity]" "h"/"\[Alpha]"^2/"c"
AbstractUnitData["\[Mu]0"] = 2 "\[Alpha]" "h"/"c"/"e"^2 (*\[TildeTilde]4\[Pi]*(1e-7+5.5e-17),exact charge*)
AbstractUnitData["\[HBar]"] = "h"/(2 Pi)
AbstractUnitData["\[Mu]pe"] = "\[Mu]pu"/"\[Mu]eu"
AbstractUnitData["\[Alpha]G"] = (AbstractUnitData["me"]/"mP")^2
AbstractUnitData["GG"] = "c" AbstractUnitData["\[HBar]"]/"mP"^2
AbstractUnitData["pc"] = "au" 3*60^3/Pi
AbstractUnitData["th"] = 1000 AbstractUnitData["pc"]/"H0"
AbstractUnitData["\[CapitalLambda]"] = 3 (1/AbstractUnitData["th"]/"c")^2 "\[CapitalOmega]\[CapitalLambda]"
AbstractUnitData["lc"] = PowerExpand[2 Sqrt[2 Pi/AbstractUnitData["\[CapitalLambda]"]]]
AbstractUnitData["mc"] = PowerExpand["c"^2/(2 Sqrt[2 Pi AbstractUnitData["\[CapitalLambda]"]] AbstractUnitData["GG"])]
AbstractUnitData["\[Rho]\[CapitalLambda]"] = "c"^4 AbstractUnitData["\[CapitalLambda]"]/(8 Pi AbstractUnitData["GG"])
AbstractUnitData["lcq"] = PowerExpand[("c" AbstractUnitData["\[HBar]"]/AbstractUnitData["\[Rho]\[CapitalLambda]"])^(1/4)]
AbstractUnitData["mcq"] = PowerExpand[(AbstractUnitData["\[HBar]"]^3 AbstractUnitData["\[Rho]\[CapitalLambda]"]/"c"^5)^(1/4)]
AbstractUnitData["ecq"] = PowerExpand[("c"^3 AbstractUnitData["\[HBar]"]^3 AbstractUnitData["\[Rho]\[CapitalLambda]"])^(1/4)]
AbstractUnitData["tcq"] = PowerExpand[AbstractUnitData["lcq"] Sqrt[AbstractUnitData["mcq"]/AbstractUnitData["ecq"]]]

EvalUnitData[x] := AbstractUnitData[x]/.Normal[UnitData]

Coupling["AbstractCoupling"] = Coupling["\[Alpha]G", "\[Alpha]", "\[Mu]eu", "\[Mu]pu"]
Coupling["AbstractUniverse"] = Coupling[AbstractUnitData["\[Alpha]G"], "\[Alpha]", "\[Mu]eu", "\[Mu]pu"]
Coupling["StandardModel"] = Coupling["AbstractUniverse"] /. Normal[UnitData]

AbstractUnitSystem["AbstractUnits"] = UnitSystem["kB", "\[HBar]", "c", "\[Mu]0", "me", "\[Lambda]", "\[Alpha]L"]
AbstractUnitSystem["AbstractUnits1"] = UnitSystem["kB1", "\[HBar]1", "c1", "\[Mu]01", "me1", "\[Lambda]1", "\[Alpha]L1"]
AbstractUnitSystem["AbstractUnits2"] = UnitSystem["kB2", "\[HBar]2", "c2", "\[Mu]02", "me2", "\[Lambda]2", "\[Alpha]L2"]
AbstractUnitSystem["SI2019"] = UnitSystem["kB", AbstractUnitData["\[HBar]"], "c", AbstractUnitData["\[Mu]0"], AbstractUnitData["me"]];

DeriveMetric[ru_, perm_] := DeriveMetric[ru, perm, AbstractUnitData["\[HBar]"], AbstractUnitData["me"]];
DeriveMetric[ru_, perm_, hbar_, mass_] := UnitSystem[ru mass/"\[Mu]eu", hbar, "c", perm, mass];

AbstractUnitSystem["Metric"] = DeriveMetric[1000 AbstractUnitData["Ru"], 4 Pi/10^7]
AbstractUnitSystem["SI1976"] = DeriveMetric[831432/100, 4 Pi/10^7]

DeriveCODATA[klitz_, joseph_] := DeriveCODATA[klitz, Nothing, 2/klitz/joseph^2/Pi]
DeriveCODATA[klitz_, _, hbar_] := DeriveCODATA[klitz, Nothing, hbar, 4 Pi "R\[Infinity]" hbar/"\[Alpha]"^2/"c"]
DeriveCODATA[klitz_, _, hbar_, mass_] := DeriveMetric[1000 AbstractUnitData["Ru"], 2 "\[Alpha]" klitz/"c", hbar, mass]

AbstractUnitSystem["CODATA"] = DeriveCODATA["RK2014", "KJ2014"];
AbstractUnitSystem["Conventional"] = DeriveCODATA["RK1990", "KJ1990"];

DeriveGaussSystem[u_, perm_, ratio_] := DeriveGaussSystem[u, perm, ratio, 1]
DeriveGaussSystem[u_, perm_, ratio_, lorentz_] := DeriveGaussSystem[u, perm, ratio, lorentz, 1/100, 1/1000]
DeriveGaussSystem[u_, perm_, ratio_, lorentz_, length_, mass_] :=
	DeriveGaussSystem[u, perm, ratio, lorentz, length, mass, PowerExpand[mass length^2]]
DeriveGaussSystem[u_, perm_, ratio_, 1, length_, mass_, energy_] :=
	UnitSystem[BoltzmannConstant[u]/energy, ReducedPlanckConstant[u]/energy, SpeedOfLight[u]/length, perm, ElectronMass[u]/mass, ratio]
DeriveGaussSystem[u_, perm_, ratio_, lorentz_, length_, mass_, energy_] :=
	UnitSystem[BoltzmannConstant[u]/energy, ReducedPlanckConstant[u]/energy, SpeedOfLight[u]/length, perm, ElectronMass[u]/mass, ratio, lorentz]

AbstractUnitSystem["Gauss"] = DeriveGaussSystem["AbstractMetric",1,4 Pi,1/100/"c"]
AbstractUnitSystem["LorentzHeaviside"] = DeriveGaussSystem["AbstractMetric",1,1,1/100/"c"]
AbstractUnitSystem["Thomson"] = DeriveGaussSystem["AbstractMetric",1,4 Pi,1/2]
AbstractUnitSystem["Kennelly"] = DeriveGaussSystem["AbstractMetric",10^-7,4 Pi,1,1,1]
AbstractUnitSystem["ESU"] = DeriveGaussSystem["AbstractMetric",(100 "c")^-2,4 Pi]
AbstractUnitSystem["EMU"] = DeriveGaussSystem["AbstractMetric",1,4 Pi]

DeriveEnergySystem[u_, time_, length_, mass_] := DeriveTempSystem[u, time, length, mass, 1]
DeriveEnergySystem[u_, time_, length_, mass_, energy_] :=
	DeriveTempSystem[u, time, length, mass, 1, PowerExpand[MagneticConstant[u] time^2/energy], energy]
DeriveTempSystem[u_, time_, length_, mass_, temp_] := Module[{energy = PowerExpand[mass length^2/time^2]},
	DeriveTempSystem[u, time, length, mass, temp, PowerExpand[MagneticConstant[u] time^2/energy]]]
DeriveTempSystem[u_, time_, length_, mass_, temp_, perm_] :=
	DeriveTempSystem[u, time, length, mass, temp, perm, PowerExpand[mass length^2/time^2]]
DeriveTempSystem[u_, time_, length_, mass_, temp_, perm_, energy_] :=
	UnitSystem[BoltzmannConstant[u] temp/energy, ReducedPlanckConstant[u]/time/energy, time SpeedOfLight[u]/length, perm, ElectronMass[u]/mass]

AbstractUnitSystem["MTS"] = DeriveEnergySystem["AbstractMetric", 1, 1, 1000]
AbstractUnitSystem["EMU2019"] = DeriveEnergySystem["AbstractSI2019", 1, 1/100, 1/1000]
AbstractUnitSystem["ESU2019"] = DeriveTempSystem["AbstractSI2019", 1, 1/100, 1/1000, 1, 1000 AbstractUnitData["\[Mu]0"]/"c"^2]
AbstractUnitSystem["Mixed"] = DeriveTempSystem["AbstractMetric", 1, 1, 1, 1, AbstractUnitData["\[Mu]0"]]
AbstractUnitSystem["English"] = DeriveTempSystem["AbstractSI2019", 1, "ft", AbstractUnitData["slug"], "rankine", 4 Pi]
AbstractUnitSystem["EnglishUS"] = DeriveTempSystem["AbstractMetric", 1, "ftUS", AbstractUnitData["slugUS"], "rankine", 4 Pi]
AbstractUnitSystem["FFF"] = DeriveTempSystem["AbstractMetric", 14 "day", "fur", 90 "lb", "rankine", 0]
AbstractUnitSystem["IAU"] = DeriveEnergySystem["AbstractMetric", "day", "au", "GMsun"/AbstractUnitData["GG"]]
AbstractUnitSystem["Hubble"] = DeriveEnergySystem["AbstractMetric", AbstractUnitData["th"], "c" AbstractUnitData["th"], 1]
AbstractUnitSystem["Cosmological"] = DeriveEnergySystem["AbstractMetric", AbstractUnitData["lc"]/"c", AbstractUnitData["lc"], AbstractUnitData["mc"]]
AbstractUnitSystem["CosmologicalQuantum"] = DeriveEnergySystem["AbstractMetric", AbstractUnitData["tcq"], AbstractUnitData["lcq"], AbstractUnitData["mcq"], AbstractUnitData["ecq"]]

AbstractUnitSystem["Planck"] = UnitSystem[1, 1, 1, 1, PowerExpand[Sqrt[4 Pi AbstractUnitData["\[Alpha]G"]]]]
AbstractUnitSystem["PlanckGauss"] = UnitSystem[1, 1, 1, 4 Pi, PowerExpand[Sqrt[AbstractUnitData["\[Alpha]G"]]]]
AbstractUnitSystem["Stoney"] = UnitSystem[1, 1/"\[Alpha]", 1, 4 Pi, PowerExpand[Sqrt[AbstractUnitData["\[Alpha]G"]/"\[Alpha]"]]]
AbstractUnitSystem["Hartree"] = UnitSystem[1,1,1/"\[Alpha]",4 Pi "\[Alpha]"^2,1]
AbstractUnitSystem["Rydberg"] = UnitSystem[1,1,2/"\[Alpha]",Pi "\[Alpha]"^2,1/2]
AbstractUnitSystem["Schrodinger"] = UnitSystem[1, 1, 1/"\[Alpha]", 4 Pi "\[Alpha]"^2, PowerExpand[Sqrt[AbstractUnitData["\[Alpha]G"]/"\[Alpha]"]]]
AbstractUnitSystem["Electronic"] = UnitSystem[1, 1/"\[Alpha]", 1, 4 Pi, 1]
AbstractUnitSystem["Natural"] = UnitSystem[1, 1, 1, 1, 1, 1, 1, "1"]
AbstractUnitSystem["NaturalGauss"] = UnitSystem[1, 1, 1, 4 Pi, 1, 1, 1, "1"]
AbstractUnitSystem["QCD"] = UnitSystem[1, 1, 1, 1, 1/AbstractUnitData["\[Mu]pe"]]
AbstractUnitSystem["QCDGauss"] = UnitSystem[1, 1, 1, 4 Pi, 1/AbstractUnitData["\[Mu]pe"]]
AbstractUnitSystem["QCDoriginal"] = UnitSystem[1, 1, 1, 4 Pi "\[Alpha]", 1/AbstractUnitData["\[Mu]pe"]]

UnitSystem["Natural"] = UnitSystem[1, 1, 1, 1, 1]
UnitSystem["NaturalGauss"] = UnitSystem[1, 1, 1, 4 Pi, 1]
Map[(UnitSystem[#] = AbstractUnitSystem[#] /. Normal[UnitData]) &,
{"AbstractUnits","AbstractUnits1","AbstractUnits2","Gauss","LorentzHeaviside","Thomson","Kennelly","ESU","ESU2019","EMU","EMU2019","MTS","Mixed","Metric","SI1976","SI2019","CODATA","Conventional","English","EnglishUS","IAU","FFF","Planck","PlanckGauss","Stoney","Hartree","Rydberg","Schrodinger","Electronic","QCD","QCDGauss","QCDoriginal","Hubble","Cosmological","CosmologicalQuantum"}]

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
Hubble = UnitSystem["Hubble"]
Cosmological = UnitSystem["Cosmological"]
CosmologicalQuantum = UnitSystem["CosmologicalQuantum"]

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

UnitSystem[u_String] := AbstractUnitSystem[StringDelete[u, "Abstract"]] /; StringStartsQ[u, "Abstract"]
Coupling[c_String] := Universe[UnitSystem[c]]

Map[(#[u_UnitSystem] := #[u, Universe[u]]) &, {PlanckMass, PlanckConstant, GravitationalConstant, ElementaryCharge}];
Map[(#[u_UnitSystem] := #[Universe[u]]) &, {GravitationalCouplingConstantElectronElectron, FineStructureConstant, ElectronRelativeAtomicMass, ProtonRelativeAtomicMass, ProtonElectronMassRatio}];
Map[(#[u_UnitSystem, c_Coupling] := #[u]) &, {BoltzmannConstant, ReducedPlanckConstant, SpeedOfLight, MagneticConstant, ElectronMass, MolarMassConstant}];
Map[(#[u_UnitSystem, s_UnitSystem] := Unit[#[s]/#[u]]) &, ConstantsList];
Map[Unprotect, ProtectedList]
Map[(If[!MemberQ[ProtectedList, #],
	#[v_?NumberQ, u_String] := #[v, UnitSystem[u]];
	#[v_?NumberQ, u_String, s_String] := #[v, UnitSystem[u], UnitSystem[s]];
	#[u_String] := #[UnitSystem[u]];
	#[u_String, s_String] := #[UnitSystem[u], UnitSystem[s]];
	#[v_?NumberQ, u_UnitSystem] := #[v, u, UnitSystem["Metric"]];
	#[v_?NumberQ, u_UnitSystem, s_UnitSystem] :=
		Module[{n = #[u, s]}, If[OneQ[n], v, v/n]];
	#[v_?NumberQ,
		u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___],
		s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := v;
	#[
		u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___],
		s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := 1;
	#[u_UnitSystem] := #[u, UnitSystem["Metric"]];,
	#[v_?NumberQ, u_String] := #[v, UnitSystem[u]];
	#[v_?NumberQ, u_String, s_String] := #[v, UnitSystem[u], UnitSystem[s]];
	#[u_String] := #[UnitSystem[u]];
	#[u_String, s_String] := #[UnitSystem[u], UnitSystem[s]];
	UnitSystem /: #[v_?NumberQ, u_UnitSystem] := #[v, u, UnitSystem["Metric"]];
	UnitSystem /: #[v_?NumberQ, u_UnitSystem, s_UnitSystem] :=
		Module[{n = #[u, s]}, If[OneQ[n], v, v/n]];
	UnitSystem /: #[v_?NumberQ,
		u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___],
		s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := v;
	UnitSystem /: #[
		u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___],
		s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := 1;
	UnitSystem /: #[u_UnitSystem] := #[u, UnitSystem["Metric"]];]
) &, ConvertList];
Map[Protect, ProtectedList]

Mass[u_UnitSystem, s_UnitSystem] := ElectronMass[u, s];
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

MonochromaticRadiation540THzLuminousEfficacy[UnitSystem[1, ___]] = 1
MonochromaticRadiation540THzLuminousEfficacy[u : UnitSystem[_?NumericQ, ___]] := Power[UnitData["Kcd"], UnitSystem["SI2019"], u]
MonochromaticRadiation540THzLuminousEfficacy[u : UnitSystem[_Around, ___]] := Power[UnitData["Kcd"], UnitSystem["SI2019"], u]

Universe[_] = Coupling["AbstractUniverse"]
Universe[c_String] := Coupling[c]
Map[(Universe[AbstractUnitSystem[#]] := Coupling["AbstractUniverse"]) &,
{"Planck","PlanckGauss","Stoney","Hartree","Rydberg","Schrodinger","Electronic","Natural","NaturalGauss","QCD","QCDGauss","QCDoriginal"}]

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

MonochromaticRadiation540THzLuminousEfficacy[AbstractUnitSystem["AbstractUnits1"]] = "Kcd1"
MonochromaticRadiation540THzLuminousEfficacy[AbstractUnitSystem["AbstractUnits2"]] = "Kcd2"
MonochromaticRadiation540THzLuminousEfficacy[u_UnitSystem] := Power["Kcd", AbstractUnitSystem["SI2019"], u]

AbstractUnitData["Mu"] = MolarMassConstant[AbstractUnitSystem["SI2019"]]

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

Map[(#[u_UnitSystem] := #[u, Universe[u]]) &, {AvogadroConstant, AtomicMassConstant, ProtonMass, EinsteinConstantSpeedOfLightToTheFourth, MolarGasConstant, StefanBoltzmannConstant, RadiationConstant, ElectricConstant, CoulombConstant, BiotSavartConstant, VacuumImpedance, FaradayConstant, JosephsonConstant, MagneticFluxQuantum, VonKlitzingConstant, ConductanceQuantum, HartreeEnergy, RydbergConstant, BohrRadius, RelativisticBohrRadius, ClassicalElectronRadius, BohrMagneton}];
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

Map[(
	#[u_String, c_String] := #[UnitSystem[u], Coupling[c]];
	#[u_String] := #[UnitSystem[u]]
) &, Join[ConstantsList, PhysicsList]]

Map[(#[c_String] := #[Coupling[c]]) &, {GravitationalCouplingConstantElectronElectron, FineStructureConstant, ElectronRelativeAtomicMass, ProtonRelativeAtomicMass, ProtonElectronMassRatio}]

Map[(symbol[ToString[#]] = #) &, Join[ConstantsList, PhysicsList, {GravitationalCouplingConstantElectronElectron, FineStructureConstant, ElectronRelativeAtomicMass, ProtonRelativeAtomicMass, ProtonElectronMassRatio}]]

Coupling /: Entity["PhysicalConstant", x_][c_Coupling] := symbol[x][c]
UnitSystem /: Entity["PhysicalConstant", x_][u_UnitSystem] := symbol[x][u]
UnitSystem /: Entity["PhysicalConstant", x_][u_UnitSystem, c_Coupling] := symbol[x][u, c]
(u_UnitSystem)[Entity["PhysicalConstant", x_]] := symbol[x][u]
(u_UnitSystem)[x_String] := symbol[x][u]
(c_Coupling)[Entity["PhysicalConstant", x_]] := symbol[x][c]
(c_Coupling)[x_String] := symbol[x][c]

UnitConstant[x_Symbol, u_] := x[u]
UnitConstant[x_String, u_] := UnitConstant[symbol[x], u]
UnitConstant[x_Symbol, u_, c_] := x[u, c]
UnitConstant[x_String, u_, c_] := UnitConstant[symbol[x], u, c]

Map[(convert[ToString[#]] = #) &, ConvertList]

ConvertUnit[x_Symbol, u_] := x[u]
ConvertUnit[x_String, u_] := convert[x][u]
ConvertUnit[x_Symbol, a_, b_] := x[a, b]
ConvertUnit[x_String, a_, b_] := convert[x][a, b]
ConvertUnit[x_Symbol, v_, u_, s_] := x[v, u, s]
ConvertUnit[x_String, v_, u_, s_] := convert[x][v, u, s]

(* more *)

DerivedUnits = <|
"NaturalUnitOfEnergy" -> Energy["AbstractNatural","AbstractSI2019"],
"NaturalUnitOfLength" -> Length["AbstractNatural","AbstractSI2019"],
"NaturalUnitOfMomentum" -> Momentum["AbstractNatural","AbstractSI2019"],
"NaturalUnitOfTime" -> Time["AbstractNatural","AbstractSI2019"],
"PlanckArea" -> Area[UnitSystem["AbstractPlanckGauss"],UnitSystem["AbstractSI2019"]]
"PlanckFrequency" -> Frequency["AbstractPlanckGauss","AbstractSI2019"],
"PlanckLength" -> Length["AbstractPlanckGauss","AbstractSI2019"],
"PlanckMassDensity" -> Density["AbstractPlanckGauss","AbstractSI2019"],
"PlanckTemperature" -> Temperature["AbstractPlanckGauss","AbstractSI2019"],
"PlanckTime" -> Time["AbstractPlanckGauss","AbstractSI2019"],
"PlanckVolume" -> Volume[UnitSystem["AbstractPlanckGauss"],UnitSystem["AbstractSI2019"]],
"AtomicUnitOfElectricConductance" -> Conductance["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfElectricChargeDensity" -> ChargeDensity["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfElectricCurrent" -> Current["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfElectricFieldStrength" -> ElectricField["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfElectricPermittivity" -> ElectricConstant["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfElectricPolarizability" -> ElectricPolarizability["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfElectricPotential" -> ElectricPotential["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfForce" -> Force["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfFrequency" -> Frequency["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfMagneticFlux" -> MagneticFlux["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfMagneticFluxDensity" -> MagneticFluxDensity["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfMagneticMoment" -> MagneticDipoleMoment["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfMomentum" -> Momentum["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfPressure" -> Pressure["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfTemperature" -> Temperature["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfTime" -> Time["AbstractHartree","AbstractSI2019"],
"AtomicUnitOfVelocity" -> Speed["AbstractHartree","AbstractSI2019"],
"RydbergEnergy" -> Energy["AbstractRydberg","AbstractSI2019"],
"SolarMass" -> Mass["AbstractIAU","AbstractSI2019"],
"AstronomicalUnit" -> Length["AbstractIAU","AbstractSI2019"],
"CosmologicalNaturalLength" -> Length["AbstractCosmological","AbstractSI2019"],
"CosmologicalNaturalMass" -> Mass["AbstractCosmological","AbstractSI2019"],
"CosmologicalNaturalTime" -> Time["AbstractCosmological","AbstractSI2019"],
"CosmologicalQuantumPointLength" -> Length["AbstractCosmologicalQuantum","AbstractSI2019"],
"CosmologicalQuantumPointMass" -> Mass["AbstractCosmologicalQuantum","AbstractSI2019"],
"CosmologicalQuantumPointEnergy" -> Energy["AbstractCosmologicalQuantum","AbstractSI2019"],
"HubbleLength" -> Length["AbstractHubble","AbstractSI2019"],
"HubbleTime" -> Time["AbstractHubble","AbstractSI2019"],
"HubbleVolume" -> Volume[UnitSystem["AbstractHubble"],UnitSystem["AbstractSI2019"]]
|>

(*
"AtomicUnitOfLength" -> Length[UnitSystem["Hartree"],UnitSystem["SI2019"]];
"lS" -> Length[UnitSystem["Stoney"], UnitSystem["SI2019"]];
"tS" -> Time[UnitSystem["Stoney"], UnitSystem["SI2019"]];
"mS" -> Mass[UnitSystem["Stoney"], UnitSystem["SI2019"]];
"lQCD" -> Length[UnitSystem["QCD"], UnitSystem["SI2019"]];
"tQCD" -> Time[UnitSystem["QCD"], UnitSystem["SI2019"]];
"mQCD" -> Mass[UnitSystem["QCD"], UnitSystem["SI2019"]];
*)

(*ConstantsList = {Hyperfine, LightSpeed, Planck, PlanckReduced,
	ElectronMass, MolarMass, Boltzmann, Permeability, Rationalization,
	Lorentz, LuminousEfficacy};
PhysicsList = {AtomicMass, ProtonMass, PlanckMass, Newton, Einstein,
	Hartree, Rydberg, Bohr, BohrReduced, ElectronRadius, Avogadro,
	Universal, Stefan, RadiationDensity, Permittivity, Coulomb, Ampere,
	BiotSavart, Charge, Faraday, Impedance, Conductance, Klitzing,
	Josephson, MagneticFlux, Magneton};*)

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
AppendTo[UnitData, "cal" -> UnitData[["kcal"]]/1000]
AppendTo[UnitData, "calth" -> thermal calorie
AppendTo[UnitData, "GMearth" -> Around[398600441.8 10^6, 8 10^5]]
AppendTo[UnitData, "GMjupiter" -> Around[1.26686534 10^17, 9 10^9]
AppendTo[UnitData, "LD" -> 384402 10^3]
AppendTo[UnitData, "ly" -> (365.25 "c" "day")/.Normal[UnitData]]*)

UnitSystem["Natural"] = UnitSystem[1, 1, 1, 1, 1]
UnitSystem["NaturalGauss"] = UnitSystem[1, 1, 1, 4 Pi, 1]
Map[(UnitSystem[#] = AbstractUnitSystem[#] /. Normal[UnitData]) &,
{"AbstractUnits","AbstractUnits1","AbstractUnits2","Gauss","LorentzHeaviside","Thomson","Kennelly","ESU","ESU2019","EMU","EMU2019","MTS","Mixed","Metric","SI1976","SI2019","CODATA","Conventional","English","EnglishUS","IAU","FFF","Planck","PlanckGauss","Stoney","Hartree","Rydberg","Schrodinger","Electronic","QCD","QCDGauss","QCDoriginal","Hubble","Cosmological","CosmologicalQuantum"}];
