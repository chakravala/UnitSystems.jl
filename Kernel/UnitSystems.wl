(* ::Package:: *)
(* This file is part of UnitSystems. It is licensed under the MIT license *)
(* UnitSystems Copyright (C) 2021 Michael Reed *)

Unprotect[UnitSystem];
ProtectedList = {Length, Area, Volume, Power, Entropy, MomentOfInertia};
UnitSystemsList = {"Metric", "SI2019", "CODATA", "Conventional", "MTS", "English",
    "EnglishUS", "FFF", "IAU", "SI1976", "Mixed", "ESU2019", "EMU2019", "EMU", "ESU",
	"Gauss", "LorentzHeaviside", "Thomson", "Kennelly", "Planck", "PlanckGauss", "Stoney",
	"Hartree", "Rydberg", "Schrodinger", "Electronic", "Natural", "NaturalGauss",
	"QCD", "QCDGauss", "QCDoriginal", "Hubble", "Cosmological", "CosmologicalQuantum"};
DimensionlessList = {GravitationalCouplingConstantElectronElectron, FineStructureConstant, ElectronRelativeAtomicMass, ProtonRelativeAtomicMass, ProtonElectronMassRatio, UniverseDarkEnergyMassDensity}
ConstantsList = {SpeedOfLight, PlanckConstant, ReducedPlanckConstant, ElectronMass,
	MolarMassConstant, BoltzmannConstant, MagneticConstant, RationalizationConstant,
	LorentzConstant, MonochromaticRadiation540THzLuminousEfficacy};
PhysicsList = {Cesium133HyperfineSplittingFrequency, HubbleParameter,
	CosmologicalConstant, AtomicMassConstant, ProtonMass, PlanckMass,
	GravitationalConstant, EinsteinConstantSpeedOfLightToTheFourth,
	HartreeEnergy, RydbergConstant, BohrRadius,
	RelativisticBohrRadius, ClassicalElectronRadius,
	AvogadroConstant, MolarGasConstant, StefanBoltzmannConstant,
	RadiationConstant, ElectricConstant, CoulombConstant, AmpereConstant,
	BiotSavartConstant, ElementaryCharge, FaradayConstant, VacuumImpedance,
	ConductanceQuantum, VonKlitzingConstant, JosephsonConstant,
	MagneticFluxQuantum, BohrMagneton};
KinematicList = {Time, Length, Area, Volume, Wavenumber,
   FuelEconomy, Frequency, FrequencyDrift, Speed, Acceleration,
   Jerk, Snap, VolumeFlowRate};
MechanicalList = {Mass, MassFlowRate, LinearMassDensity, MassPerArea, MassDensity,
    SpecificVolume, Force, Stiffness, Pressure, Compressibility,
   DynamicViscosity, KinematicViscosity, MomentOfInertia, Momentum,
   AngularMomentum, ForceOnsetRate, Energy, SpecificEnergy, Action, EnergyPerArea,
   Power, PowerDensity, Irradiance, PowerGradient, SoundExposure,
   AcousticImpedance, SpecificAcousticImpedance, Admittance, Compliance, Inertance};
ElectromagneticList = {ElectricCharge, ElectricChargeDensity, LinearElectricChargeDensity,
   ElectricChargePerMass, ElectricalMobility, ElectricCapacitance, MagneticInductance, MagneticReluctance, MagneticPermeance,
    ElectricPermittivity, MagneticPermeability, MagneticSusceptibility, ElectricCurrent, ElectricConductance,
   SpecificSusceptibility, DemagnetizationFactor, MagneticVectorPotential,
   ElectricPotential, MagnetomotiveForce, ElectricFieldStrength, MagneticFieldStrength,
    ElectricFlux, MagneticFlux, ElectricFluxDensity,
   MagneticFluxDensity, ElectricDipoleMoment, MagneticDipoleMoment,
   ElectricPolarizability, MagneticPolarizability, MagneticMoment,
   Magnetizability, Magnetization, SpecificMagnetization, MagneticRigidity,
   MagneticPoleStrength};
ThermodynamicList = {Temperature, Entropy, SpecificEntropy,
   EntropyPerVolume, ThermalConductivity, ThermalConductance,
   ThermalResistance, ThermalExpansion, TemperatureGradient};
MolarList = {MolarMass, Molality, Amount, AmountConcentration, MolarVolume,
   MolarEntropy, MolarEnergy, MolarConductivity, MolarMagneticSusceptibility,
   Catalysis, Specificity};
PhotometricList = {LuminousFlux, Luminance, LuminousEnergy,
   LuminousExposure, LuminousEfficacyOfRadiation};
MechanicsList = Join[KinematicList, MechanicalList];
ConvertList =
  Join[MechanicsList, ElectromagneticList, ThermodynamicList,
   MolarList, PhotometricList];

measure[x_] := x;
GravitationalCouplingConstantElectronElectron[Coupling[\[Alpha]G_, ___]] := measure[\[Alpha]G];
FineStructureConstant[Coupling[_, \[Alpha]_, ___]] := measure[\[Alpha]];
ElectronRelativeAtomicMass[Coupling[_, _, \[Mu]eu_, ___]] := measure[\[Mu]eu];
ProtonRelativeAtomicMass[Coupling[_, _, _, \[Mu]pu_, ___]] := measure[\[Mu]pu];
UniverseDarkEnergyMassDensity[Coupling[_, _, _, _, OL_, ___]] := measure[OL];
ProtonElectronMassRatio[c_Coupling] := ProtonRelativeAtomicMass[c]/ElectronRelativeAtomicMass[c];

UnitSystem[u_UnitSystem] := u
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

Universe[UnitSystem[_?NumericQ, ___]] := Coupling["StandardModel"]
Universe[UnitSystem[_Around, ___]] := Coupling["StandardModel"]
Unit[x_, y_ : 1] := PowerExpand[x];

UnitSystem[u_String] := If[AbstractUnitSystemQ[u], AbstractUnitSystem[StringDelete[u, "Abstract"]], DimensionSystem[StringDelete[u, "Dimension"]]] /; AbstractUnitSystemQ[u] || StringStartsQ[u, "Dimension"]
Coupling[c_String] := Universe[UnitSystem[c]]

AbstractUnitSystemQ[u_String] := StringStartsQ[u, "Abstract"]
AbstractUnitSystemQ[u:UnitSystem[_,__]] := MemberQ[Values[AbstractUnitSystem],u]
AbstractUnitSystemQ[_] := False
MatchSystem[u_,s_] := If[AbstractUnitSystemQ[u], AbstractUnitSystem[s], UnitSystem[s]]
DefaultSystem[u_] := MatchSystem[u, "SI2019"]

Map[(
	#[u_String, c_String] := #[UnitSystem[u], Coupling[c]];
	#[u_String] := #[UnitSystem[u]]
) &, Join[ConstantsList, PhysicsList]]

Map[(#[c_String] := #[Coupling[c]]) &, DimensionlessList]

Get["UnitSystems`systems`"]

Universe[_] = Coupling["AbstractUniverse"]
Universe[c_String] := Coupling[c]
Map[(Universe[AbstractUnitSystem[#]] := Coupling["AbstractUniverse"]) &,
{"Planck","PlanckGauss","Stoney","Hartree","Rydberg","Schrodinger","Electronic","Natural","NaturalGauss","QCD","QCDGauss","QCDoriginal"}]

OneQ[1] := True; OneQ[-1] := False; OneQ[True] := True
OneQ[0] := False; OneQ[_] := False; OneQ[False] := False
OneQ[x_?NumberQ] := x == 1

CouplingSystem[u_String] := CouplingSystem[UnitSystem[u]]
CouplingSystem[UnitSystem[u_String]] := Coupling[u]
CouplingSystem[u_UnitSystem] := u

Map[(#[u_UnitSystem] := #[u, Universe[u]]) &, {PlanckMass, PlanckConstant, GravitationalConstant, ElementaryCharge}];
Map[(#[u_UnitSystem] := #[Universe[u]]) &, DimensionlessList];
Map[(#[u_UnitSystem, c_Coupling] := #[u]) &, {BoltzmannConstant, ReducedPlanckConstant, SpeedOfLight, MagneticConstant, ElectronMass, MolarMassConstant}];
Map[(
	#[u_UnitSystem, s_UnitSystem] := Unit[#[s]/#[u]];
	#[u_String, s_String] := #[UnitSystem[u], CouplingSystem[s]];
) &, Join[ConstantsList,PhysicsList]]
Map[Unprotect, ProtectedList]
Map[(If[!MemberQ[ProtectedList, #],
	#[v_, u_String] := #[v, UnitSystem[u]];
	#[v_, u_String, s_String] := #[v, UnitSystem[u], UnitSystem[s]];
	#[u_String] := #[UnitSystem[u]];
	#[u_String, s_String] := #[UnitSystem[u], UnitSystem[s]];
	#[v_, u_UnitSystem] := #[v, u, DefaultSystem[u]];
	#[v_, u_UnitSystem, s_UnitSystem] :=
		Module[{n = #[u, s]}, If[OneQ[n], v, v/n]];
	#[v_,
		u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___],
		s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := v;
	#[
		u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___],
		s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := 1;
	#[u_UnitSystem] := #[u, DefaultSystem[u]];,
	#[v_, u_String] := #[v, UnitSystem[u]];
	#[v_, u_String, s_String] := #[v, UnitSystem[u], UnitSystem[s]];
	#[u_String] := #[UnitSystem[u]];
	#[u_String, s_String] := #[UnitSystem[u], UnitSystem[s]];
	UnitSystem /: #[v_, u_UnitSystem] := #[v, u, DefaultSystem[u]];
	UnitSystem /: #[v_, u_UnitSystem, s_UnitSystem] :=
		Module[{n = #[u, s]}, If[OneQ[n], v, v/n]];
	UnitSystem /: #[v_,
		u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___],
		s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := v;
	UnitSystem /: #[
		u : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___],
		s : UnitSystem[kB_, \[HBar]_, c_, \[Mu]0_, me_, ___]] := 1;
	UnitSystem /: #[u_UnitSystem] := #[u, DefaultSystem[u]];]
) &, ConvertList];
Map[Protect, ProtectedList]

Get["UnitSystems`physics`"]

symbol = Association[Map[Rule[ToString[#],#] &, Join[ConstantsList, PhysicsList, DimensionlessList]]]

UnitConstant[x_Symbol, u_] := x[u]
UnitConstant[x_String, u_] := UnitConstant[symbol[x], u]
UnitConstant[x_Symbol, u_, c_] := x[u, c]
UnitConstant[x_String, u_, c_] := UnitConstant[symbol[x], u, c]

(u_UnitSystem)[Entity["PhysicalConstant", x_]] := UnitConstant[x, u]
(u_UnitSystem)[x_String] := UnitConstant[x, u]
(c_Coupling)[Entity["PhysicalConstant", x_]] := UnitConstant[x, c]
(c_Coupling)[x_String] := UnitConstant[x, c]
Coupling /: Entity["PhysicalConstant", x_][c_Coupling] := UnitConstant[x, c]
UnitSystem /: Entity["PhysicalConstant", x_][u_UnitSystem] := UnitConstant[x, u]
UnitSystem /: Entity["PhysicalConstant", x_][u_UnitSystem, c_Coupling] := UnitConstant[x, u, c]

Unprotect[Entity]
Entity["UnitSystem", u_][Entity["PhysicalConstant", x_]] := UnitConstant[x, u]
Entity["PhysicalConstant", x_][Entity["UnitSystem", u_]] := UnitConstant[x, u]
Entity["PhysicalConstant", x_][Entity["UnitSystem", u_], c_] := UnitConstant[x, u, c]
Entity["UnitSystem", u_][Entity["UnitSystem", s_]] := ComparePhysics[AbstractUnitSystem[u], AbstractUnitSystem[s]]
Protect[Entity]

convert = Association[Map[Rule[ToString[#],#] &, ConvertList]]

ConvertUnit[x_Symbol, u_] := x[u]
ConvertUnit[x_String, u_] := ConvertUnit[convert[x], u]
ConvertUnit[x_Symbol, a_, b_] := x[a, b]
ConvertUnit[x_String, a_, b_] := ConvertUnit[convert[x], a, b]
ConvertUnit[x_Symbol, v_, u_, s_] := x[v, u, s]
ConvertUnit[x_String, v_, u_, s_] := ConvertUnit[convert[x], v, u, s]
ConvertUnit[Area, u_] := Area[UnitSystem[u],DefaultSystem[u]]
ConvertUnit[Area, u_String, s_String] := Area[UnitSystem[u],UnitSystem[s]]
ConvertUnit[Area, u_UnitSystem, s_UnitSystem] := Area[u ,s]
ConvertUnit[Area, v_, u_String] := Area[v, UnitSystem[u],DefaultSystem[u]]
ConvertUnit[Area, v_, u_String, s_String] := Area[v, UnitSystem[u],UnitSystem[s]]
ConvertUnit[Volume, u_] := Volume[UnitSystem[u],DefaultSystem[u]]
ConvertUnit[Volume, u_String, s_String] := Volume[UnitSystem[u],UnitSystem[s]]
ConvertUnit[Volume, u_UnitSystem, s_UnitSystem] := Volume[u, s]
ConvertUnit[Volume, v_, u_] := Volume[v, UnitSystem[u],DefaultSystem[u]]
ConvertUnit[Volume, v_, u_String, s_String] := Volume[v, UnitSystem[u],UnitSystem[s]]
ConvertUnit[Entropy, u_] := Entropy[UnitSystem[u],DefaultSystem[u]]
ConvertUnit[Entropy, u_String, s_String] := Entropy[UnitSystem[u], UnitSystem[s]]
ConvertUnit[Entropy, u_UnitSystem, s_UnitSystem] := Entropy[u, s]
ConvertUnit[Entropy, v_, u_] := Entropy[v, UnitSystem[u], DefaultSystem[u]]
ConvertUnit[Entropy, v_, u_String, s_String] := Entropy[v, UnitSystem[u], UnitSystem[s]]
ConvertUnit[MomentOfInertia, u_] := MomentOfInertia[UnitSystem[u],DefaultSystem[u]]
ConvertUnit[MomentOfInertia, u_String, s_String] := MomentOfInertia[UnitSystem[u],UnitSystem[s]]
ConvertUnit[MomentOfInertia, u_UnitSystem, s_UnitSystem] := MomentOfInertia[u ,s]
ConvertUnit[MomentOfInertia, v_, u_String] := MomentOfInertia[v, UnitSystem[u],DefaultSystem[u]]
ConvertUnit[MomentOfInertia, v_, u_String, s_String] := MomentOfInertia[v, UnitSystem[u],UnitSystem[s]]

CompareBase[a_, b_] := CompareUnits[a, b, ConstantsList]
CompareDerived[a_, b_] := CompareUnits[a, b, PhysicsList]
CompareConstants[a_, b_] := CompareUnits[a, b, Join[ConstantsList, PhysicsList]]
ComparePhysics[a_, b_] := CompareUnits[a, b, ConvertList]
CompareUnits[a_, b_] := CompareUnits[a, b, Join[ConstantsList, PhysicsList, ConvertList]]
CompareUnits[a_String, b_, c_] := CompareUnits[UnitSystem[a], b, c]
CompareUnits[a_, b_String, c_] := CompareUnits[a, UnitSystem[b], c]
CompareUnits[a_String, b_String, c_] := CompareUnits[UnitSystem[a], UnitSystem[b], c]
CompareUnits[a_, b_, c_] := Association[Map[ToString[#] -> #[a, b] &, c]]
DimensionUnits[u_, l_] := CompareUnits["AbstractNatural",DimensionSystem[u], l]
DimensionUnits[u_] := DimensionUnits[u, Join[ConstantsList, PhysicsList, ConvertList]]
DimensionPhysics[u_] := DimensionUnits[u, ConvertList]
DimensionConstants[u_] := DimensionUnits[u, Join[ConstantsList, PhysicsList]]
DimensionDerived[u_] := DimensionUnits[u, PhysicsList]
DimensionBase[u_] := DimensionUnits[u, ConstantsList]

DerivedConstants[u_] := DerivedConstants[u, Join[ConstantsList, PhysicsList]]
DerivedConstants[u_, l_] := Association[Map[ToString[#] -> #[u] &, l]]
DerivedConstants[u_String, l_] := If[AbstractUnitSystemQ[u], DerivedConstants[UnitSystem[u],l], Association@Normal[DerivedConstants[StringJoin["Abstract",u],l]/.Normal[UnitData]]]

EntityUnregister["UnitSystem"]
CreateEntity[u_String] := u -> <|
	"Abstract" :> AbstractUnitSystem[u],
	"Dimensions" :> DimensionSystem[u],
	"Value" :> UnitSystem[u],
	"Formulas" :> DerivedConstants[AbstractUnitSystem[u]],
	"Constants" :> DerivedConstants[u],
	"Quantities" :> DimensionPhysics[u]|>
EntityRegister[EntityStore["UnitSystem" -> <|
	"Entities" -> Association[Map[CreateEntity, UnitSystemsList]]|>]]

(* more *)

DerivedUnits = <|
"NaturalUnitOfEnergy" -> Energy["AbstractNatural"],
"NaturalUnitOfLength" -> Length["AbstractNatural"],
"NaturalUnitOfMomentum" -> Momentum["AbstractNatural"],
"NaturalUnitOfTime" -> Time["AbstractNatural"],
"PlanckArea" -> Area[UnitSystem["AbstractPlanckGauss"]],
"PlanckFrequency" -> Frequency["AbstractPlanckGauss"],
"PlanckLength" -> Length["AbstractPlanckGauss"],
"PlanckMassDensity" -> MassDensity["AbstractPlanckGauss"],
"PlanckTemperature" -> Temperature["AbstractPlanckGauss"],
"PlanckTime" -> Time["AbstractPlanckGauss"],
"PlanckVolume" -> Volume[UnitSystem["AbstractPlanckGauss"]],
"AtomicUnitOfElectricConductance" -> ElectricConductance["AbstractHartree"],
"AtomicUnitOfElectricElectricChargeDensity" -> ElectricChargeDensity["AbstractHartree"],
"AtomicUnitOfElectricElectricCurrent" -> ElectricCurrent["AbstractHartree"],
"AtomicUnitOfElectricFieldStrengthStrength" -> ElectricFieldStrength["AbstractHartree"],
"AtomicUnitOfElectricElectricPermittivity" -> ElectricPermittivity["AbstractHartree"],
"AtomicUnitOfElectricPolarizability" -> ElectricPolarizability["AbstractHartree"],
"AtomicUnitOfElectricPotential" -> ElectricPotential["AbstractHartree"],
"AtomicUnitOfForce" -> Force["AbstractHartree"],
"AtomicUnitOfFrequency" -> Frequency["AbstractHartree"],
"AtomicUnitOfMagneticFlux" -> MagneticFlux["AbstractHartree"],
"AtomicUnitOfMagneticFluxDensity" -> MagneticFluxDensity["AbstractHartree"],
"AtomicUnitOfMagneticMoment" -> MagneticDipoleMoment["AbstractHartree"],
"AtomicUnitOfMomentum" -> Momentum["AbstractHartree"],
"AtomicUnitOfPressure" -> Pressure["AbstractHartree"],
"AtomicUnitOfTemperature" -> Temperature["AbstractHartree"],
"AtomicUnitOfTime" -> Time["AbstractHartree"],
"AtomicUnitOfVelocity" -> Speed["AbstractHartree"],
"RydbergEnergy" -> Energy["AbstractRydberg"],
"SolarMass" -> Mass["AbstractIAU"],
"AstronomicalUnit" -> Length["AbstractIAU"],
"CosmologicalNaturalLength" -> Length["AbstractCosmological"],
"CosmologicalNaturalMass" -> Mass["AbstractCosmological"],
"CosmologicalNaturalTime" -> Time["AbstractCosmological"],
"CosmologicalQuantumPointLength" -> Length["AbstractCosmologicalQuantum"],
"CosmologicalQuantumPointMass" -> Mass["AbstractCosmologicalQuantum"],
"CosmologicalQuantumPointEnergy" -> Energy["AbstractCosmologicalQuantum"],
"HubbleLength" -> Length["AbstractHubble"],
"HubbleTime" -> Time["AbstractHubble"],
"HubbleVolume" -> Volume[UnitSystem["AbstractHubble"]]|>;

(*AppendTo[DerivedUnits, {
"AtomicUnitOfLength" -> Length["AbstractHartree"],
"lS" -> Length["AbstractStoney"],
"tS" -> Time["AbstractStoney"],
"mS" -> Mass["AbstractStoney"],
"lQCD" -> Length["AbstractQCD"],
"tQCD" -> Time["AbstractQCD"],
"mQCD" -> Mass["AbstractQCD"]}]*)

(*AppendTo[UnitData, {
"GMearth" -> Around[398600441.8 10^6, 8 10^5],
"GMjupiter" -> Around[1.26686534 10^17, 9 10^9],
"LD" -> 384402 10^3,
"kelvin" -> 9/5,
"atm" -> 101325,
"kcalth" -> 4184,
"kcal4" -> 4204,
"kcal10" -> 4185+1/2,
"kcal20" -> 4182,
"kcalm" -> 4190,
"kcalit" -> 4186+8/10}]
AppendTo[AbstractUnitData, {
"cal4" -> "kcal4"/1000,
"cal10" -> "kcal10"/1000,
"cal20" -> "kcal20"/1000,
"calm" -> "kcalm"/1000,
"calit" -> "kcalit"/1000,
"calth" -> "kcalth"/1000,
"kcal" -> "kcalth",
"cal" -> "kcalth"/1000,
"calth" -> "kcalth"/1000,
"ly" -> 36525/100 "c" "day"}]*)

(*ConstantsList = {Hyperfine, LightSpeed, Planck, PlanckReduced,
	ElectronMass, MolarMass, Boltzmann, Permeability, Rationalization,
	Lorentz, LuminousEfficacy};
PhysicsList = {AtomicMass, ProtonMass, PlanckMass, Newton, Einstein,
	Hartree, Rydberg, Bohr, BohrReduced, ElectronRadius, Avogadro,
	Universal, Stefan, RadiationDensity, Permittivity, Coulomb, Ampere,
	BiotSavart, Charge, Faraday, Impedance, Conductance, Klitzing,
	Josephson, MagneticFlux, Magneton};*)

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

{SI,MKS,CGS, CGS2019, CGSm, CGSe, HLU} = {SI2019, Metric, Gauss, EMU2019, EMU, ESU, LorentzHeaviside};
