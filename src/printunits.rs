
// Copyright (C) 2020 Michael Reed

pub mod lib;
pub use lib as us;
pub use lib::UnitSystem;

pub fn print_unitsystem(u: &UnitSystem) {
    println!("boltzmann: {}",us::boltzmann(u));
    println!("planckreduced: {}",us::planckreduced(u));
    println!("lightspeed: {}",us::lightspeed(u));
    println!("vacuumpermeability: {}",us::vacuumpermeability(u));
    println!("electronmass: {}",us::electronmass(u));
    println!("molarmass: {}",us::molarmass(u));
    println!("rationalization: {}",us::rationalization(u));
    println!("lorentz: {}",us::lorentz(u));
}

pub fn print_constants(u: &UnitSystem) {
    print_unitsystem(u);
    println!("hyperfine: {}",us::hyperfine(u));
    println!("luminousefficacy: {}",us::luminousefficacy(u));
    println!("avogadro: {}",us::avogadro(u));
    println!("planck: {}",us::planck(u));
    println!("atomicmass: {}",us::atomicmass(u));
    println!("protonmass: {}",us::protonmass(u));
    println!("planckmass: {}",us::planckmass(u));
    println!("newton: {}",us::newton(u));
    println!("einstein: {}",us::einstein(u));
    println!("universal: {}",us::universal(u));
    println!("stefan: {}",us::stefan(u));
    println!("radiationdensity: {}",us::radiationdensity(u));
    println!("vacuumpermittivity: {}",us::vacuumpermittivity(u));
    println!("coulomb: {}",us::coulomb(u));
    println!("biotsavart: {}",us::biotsavart(u));
    println!("ampere: {}",us::ampere(u));
    println!("vacuumimpedance: {}",us::vacuumimpedance(u));
    println!("elementarycharge: {}",us::elementarycharge(u));
    println!("faraday: {}",us::faraday(u));
    println!("josephson: {}",us::josephson(u));
    println!("magneticfluxquantum: {}",us::magneticfluxquantum(u));
    println!("klitzing: {}",us::klitzing(u));
    println!("conductancequantum: {}",us::conductancequantum(u));
    println!("hartree: {}",us::hartree(u));
    println!("rydberg: {}",us::rydberg(u));
    println!("bohr: {}",us::bohr(u));
    println!("bohrreduced: {}",us::bohrreduced(u));
    println!("electronradius: {}",us::electronradius(u));
    println!("magneton: {}",us::magneton(u));
}

pub fn print_kinematic(u: &UnitSystem, s: &UnitSystem) {
    println!("length: {}",us::length(u,s));
    println!("area: {}",us::area(u,s));
    println!("volume: {}",us::volume(u,s));
    println!("wavenumber: {}",us::wavenumber(u,s));
    println!("fuelefficiency: {}",us::fuelefficiency(u,s));
    println!("time: {}",us::time(u,s));
    println!("frequency: {}",us::frequency(u,s));
    println!("frequencydrift: {}",us::frequencydrift(u,s));
    println!("speed: {}",us::speed(u,s));
    println!("acceleration: {}",us::acceleration(u,s));
    println!("jerk: {}",us::jerk(u,s));
    println!("snap: {}",us::snap(u,s));
    println!("volumeflow: {}",us::volumeflow(u,s));
    println!("specificenergy: {}",us::specificenergy(u,s));
    println!("mass: {}",us::mass(u,s));
    println!("energy: {}",us::energy(u,s));
    println!("power: {}",us::power(u,s));
    println!("force: {}",us::force(u,s));
    println!("pressure: {}",us::pressure(u,s));
}

pub fn print_mechanical(u: &UnitSystem, s: &UnitSystem) {
    println!("momentum: {}",us::momentum(u,s));
    println!("angularmomentum: {}",us::angularmomentum(u,s));
    println!("yank: {}",us::yank(u,s));
    println!("areadensity: {}",us::areadensity(u,s));
    println!("density: {}",us::density(u,s));
    println!("specificvolume: {}",us::specificvolume(u,s));
    println!("action: {}",us::action(u,s));
    println!("stiffness: {}",us::stiffness(u,s));
    println!("intensity: {}",us::intensity(u,s));
    println!("diffusivity: {}",us::diffusivity(u,s));
    println!("viscosity: {}",us::viscosity(u,s));
    println!("lineardensity: {}",us::lineardensity(u,s));
    println!("massflow: {}",us::massflow(u,s));
    println!("spectralflux: {}",us::spectralflux(u,s));
    println!("powerdensity: {}",us::powerdensity(u,s));
    println!("compressibility: {}",us::compressibility(u,s));
    println!("fluence: {}",us::fluence(u,s));
    println!("rotationalinertia: {}",us::rotationalinertia(u,s));
    println!("soundexposure: {}",us::soundexposure(u,s));
    println!("specificimpedance: {}",us::specificimpedance(u,s));
    println!("acousticimpedance: {}",us::acousticimpedance(u,s));
    println!("admittance: {}",us::admittance(u,s));
    println!("compliance: {}",us::compliance(u,s));
    println!("inertance: {}",us::inertance(u,s));
}

pub fn print_electromagnetic(u: &UnitSystem, s: &UnitSystem) {
    println!("charge: {}",us::charge(u,s));
    println!("current: {}",us::current(u,s));
    println!("electricpotential: {}",us::electricpotential(u,s));
    println!("capacitance: {}",us::capacitance(u,s));
    println!("resistance: {}",us::resistance(u,s));
    println!("conductance: {}",us::conductance(u,s));
    println!("magneticflux: {}",us::magneticflux(u,s));
    println!("magneticfluxdensity: {}",us::magneticfluxdensity(u,s));
    println!("inductance: {}",us::inductance(u,s));
    println!("electricfluxdensity: {}",us::electricfluxdensity(u,s));
    println!("chargedensity: {}",us::chargedensity(u,s));
    println!("currentdensity: {}",us::currentdensity(u,s));
    println!("conductivity: {}",us::conductivity(u,s));
    println!("permittivity: {}",us::permittivity(u,s));
    println!("permeability: {}",us::permeability(u,s));
    println!("electricfield: {}",us::electricfield(u,s));
    println!("magneticfield: {}",us::magneticfield(u,s));
    println!("exposure: {}",us::exposure(u,s));
    println!("resistivity: {}",us::resistivity(u,s));
    println!("linearchargedensity: {}",us::linearchargedensity(u,s));
    println!("magneticdipolemoment: {}",us::magneticdipolemoment(u,s));
    println!("mobility: {}",us::mobility(u,s));
    println!("reluctance: {}",us::reluctance(u,s));
    println!("vectorpotential: {}",us::vectorpotential(u,s));
    println!("magneticmoment: {}",us::magneticmoment(u,s));
    println!("rigidity: {}",us::rigidity(u,s));
    println!("susceptibility: {}",us::susceptibility(u,s));
    println!("electricflux: {}",us::electricflux(u,s));
    println!("electricdipolemoment: {}",us::electricdipolemoment(u,s));
    println!("magneticpotential: {}",us::magneticpotential(u,s));
    println!("polestrength: {}",us::polestrength(u,s));
    println!("permeance: {}",us::permeance(u,s));
    println!("specificsusceptibility: {}",us::specificsusceptibility(u,s));
    println!("magnetizability: {}",us::magnetizability(u,s));
    println!("electricpolarizability: {}",us::electricpolarizability(u,s));
    println!("magneticpolarizability: {}",us::magneticpolarizability(u,s));
    println!("magnetization: {}",us::magnetization(u,s));
    println!("specificmagnetization: {}",us::specificmagnetization(u,s));
    println!("demagnetizingfactor: {}",us::demagnitizingfactor(u,s));
}

pub fn print_thermodynamic(u: &UnitSystem, s: &UnitSystem) {
    println!("temperature: {}",us::temperature(u,s));
    println!("entropy: {}",us::entropy(u,s));
    println!("specificentropy: {}",us::specificentropy(u,s));
    println!("volumeheatcapacity: {}",us::volumeheatcapacity(u,s));
    println!("thermalconductivity: {}",us::thermalconductivity(u,s));
    println!("thermalconductance: {}",us::thermalconductance(u,s));
    println!("thermalresistance: {}",us::thermalresistance(u,s));
    println!("thermalexpansion: {}",us::thermalexpansion(u,s));
    println!("lapserate: {}",us::lapserate(u,s));
}

pub fn print_molar(u: &UnitSystem, s: &UnitSystem) {
    println!("molality: {}",us::molality(u,s));
    println!("mole: {}",us::mole(u,s));
    println!("molarity: {}",us::molarity(u,s));
    println!("molarvolume: {}",us::molarvolume(u,s));
    println!("molarenergy: {}",us::molarenergy(u,s));
    println!("molarconductivity: {}",us::molarconductivity(u,s));
    println!("molarsusceptibility: {}",us::molarsusceptibility(u,s));
    println!("catalysis: {}",us::catalysis(u,s));
    println!("specificity: {}",us::specificity(u,s));
}

pub fn print_photometric(u: &UnitSystem, s: &UnitSystem) {
    println!("luminousflux: {}",us::luminousflux(u,s));
    println!("luminance: {}",us::luminance(u,s));
    println!("luminousenergy: {}",us::luminousenergy(u,s));
    println!("luminousexposure: {}",us::luminousexposure(u,s));
}

pub fn print_conversions(u: &UnitSystem, s: &UnitSystem) {
    print_kinematic(u,s);
    print_mechanical(u,s);
    print_electromagnetic(u,s);
    print_thermodynamic(u,s);
    print_molar(u,s);
    print_photometric(u,s);
}

fn main() {
    print_constants(&us::SI2019);
}
