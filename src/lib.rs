
// Copyright (C) 2020 Michael Reed

pub struct UnitSystem {
    kb: f64, // boltzmann
    hbar: f64, // planckreduced
    c: f64, // lightspeed
    m0: f64, // vacuumpermeability
    me: f64, // electronmass
    mu: f64, // molarmass
    lambda: f64, // rationalization
    al: f64 // lorentz
}

pub use std::f64::consts::PI;
use std::f64::consts::FRAC_2_SQRT_PI;

// conversions
pub const G0: f64 = 9.80665;
pub const FT: f64 = 0.3048;
pub const FTUS: f64 = 1200./3937.;
pub const LBM: f64 = G0/FT;
pub const LBMUS: f64 = G0/FTUS;
pub const SLUG: f64 = 0.45359237*LBM;
pub const SLUGUS: f64 = 0.45359237*LBMUS;
pub const RANKINE: f64 = 5./9.;
pub const KELVIN: f64 = 1.8;
pub const ATM: f64 = 101325.0;

// constants
pub const DVCS: f64 = 9192631770.0;
pub const KCD: f64 = 683.002;
pub const MP: f64 = 2.176434e-8;
pub const NA: f64 = 6.02214076e23;
pub const KB: f64 = 1.380649e-23;
pub const H: f64 = 6.62607015e-34;
pub const C: f64 = 299792458.;
pub const E: f64 = 1.602176634e-19;
pub const MEU: f64 = 1./1822.888486209;
pub const MPU: f64 = 1.007276466621;
pub const AINV: f64 = 137.035999084;
pub const RINF: f64 = 10973731.5681601;
pub const AL: f64 = 0.01/C;
pub const M0: f64 = 2.*H/C/AINV/(E*E);
pub const HBAR: f64 = H/2./PI;
pub const DM0: f64 = M0-4.*PI*1e-7;
pub const MPE: f64 = MPU/MEU;
pub const RU: f64 = NA*KB;
pub const ME: f64 = AINV*AINV*RINF*2.*H/C;
pub const MU: f64 = ME*NA/MEU;
pub const RK1990: f64 = 25812.807;
pub const RK2014: f64 = 25812.8074555;
pub const KJ1990: f64 = 4.835979e14;
pub const KJ2014: f64 = 4.835978525e14;
pub const HBAR1990: f64 = 2./RK1990/(KJ1990*KJ1990)/PI;
pub const HBAR2014: f64 = 2./RK2014/(KJ2014*KJ2014)/PI;
pub const ME1990: f64 = AINV*AINV*RINF*4.*PI*HBAR1990/C;
pub const ME2014: f64 = AINV*AINV*RINF*4.*PI*HBAR2014/C;

// engineering units

pub const GAUSS: UnitSystem = UnitSystem{kb: 1e10*RU*ME/MEU, hbar: 1e7*HBAR, c: 1e2*C, m0: 1., me: 1e3*ME, mu: 1., lambda: 4.*PI, al: 0.01/C};
pub const LORENTZHEAVISIDE: UnitSystem = UnitSystem{kb: 1e10*RU*ME/MEU, hbar: 1e7*HBAR, c: 1e2*C, m0: 1., me: 1e3*ME, mu: 1., lambda: 1., al: 0.01/C};
pub const THOMSON: UnitSystem = UnitSystem{kb: 1e10*RU*ME/MEU, hbar: 1e7*HBAR, c: 1e2*C, m0: 1., me: 1e3*ME, mu: 1., lambda: 4.*PI, al: 0.5};
pub const KENNELLY: UnitSystem = UnitSystem{kb: RU*ME/MEU/0.001, hbar: HBAR, c: C, m0: 1e-7, me: ME, mu: 0.001, lambda: 4.*PI, al: 1.};
pub const ESU: UnitSystem = UnitSystem{kb: 1e10*RU*ME/MEU, hbar: 1e7*HBAR, c: 1e2*C, m0: 1./((1e2*C)*(1e2*C)), me: 1e3*ME, mu: 1., lambda: 4.*PI, al: 1.};
pub const ESU2019: UnitSystem = UnitSystem{kb: 1e7*KB, hbar: 1e7*HBAR, c: 1e2*C, m0: 1e3*M0/(C*C), me: 1e3*ME, mu: 1e3*MU, lambda: 1., al: 1.};
pub const EMU: UnitSystem = UnitSystem{kb: 1e10*RU*ME/MEU, hbar: 1e7*HBAR, c: 1e2*C, m0: 1., me: 1e3*ME, mu: 1., lambda: 4.*PI, al: 1.};
pub const EMU2019: UnitSystem = UnitSystem{kb: 1e7*KB, hbar: 1e7*HBAR, c: 1e2*C, m0: 1e7*M0, me: 1e3*ME, mu: 1e3*MU, lambda: 1., al: 1.};
pub const MTS: UnitSystem = UnitSystem{kb: 1e6*RU*ME/MEU, hbar: 1e3*HBAR, c: C, m0: 4.*PI/1e4, me: ME/1e3, mu: MU/1e3, lambda: 1e-6, al: 1.};
pub const MIXED: UnitSystem = UnitSystem{kb: RU*ME/MEU/0.001, hbar: HBAR, c: C, m0: M0, me: ME, mu: 0.001, lambda: 1., al: 1.};
pub const METRIC: UnitSystem = UnitSystem{kb: RU*ME/MEU/0.001, hbar: HBAR, c: C, m0: 4.*PI*1e-7, me: ME, mu: 0.001, lambda: 1., al: 1.};
pub const SI1976: UnitSystem = UnitSystem{kb: 8.31432*ME/MEU/0.001, hbar: HBAR, c: C, m0: 4.*PI*1e-7, me: ME, mu: 0.001, lambda: 1., al: 1.};
pub const SI2019: UnitSystem = UnitSystem{kb: KB, hbar: HBAR, c: C, m0: M0, me: ME, mu: MU, lambda: 1., al: 1.};
pub const CODATA: UnitSystem = UnitSystem{kb: RU*ME2014/MEU/0.001, hbar: HBAR2014, c: C, m0: 2.*RK2014/C/AINV, me: ME2014, mu: 0.001, lambda: 1., al: 1.};
pub const CONVENTIONAL: UnitSystem = UnitSystem{kb: RU*ME1990/MEU/0.001, hbar: HBAR1990, c: C, m0: 2.*RK1990/C/AINV, me: ME1990, mu: 0.001, lambda: 1., al: 1.};
pub const ENGLISH: UnitSystem = UnitSystem{kb: KB*RANKINE/SLUG/(FT*FT), hbar: HBAR/SLUG/(FT*FT), c: C/FT, m0: 4.*PI, me: ME/SLUG, mu: 1e3*MU, lambda: 1., al: 1.};
pub const ENGLISHUS: UnitSystem = UnitSystem{kb: KB*RANKINE/SLUG/(FTUS*FTUS), hbar: HBAR/SLUG/(FTUS*FTUS), c: C/FTUS, m0: 4.*PI, me: ME/SLUG, mu: 1., lambda: 1., al: 1.};

// astronomical units

pub const GMSUN: f64 =  1.32712442099e20;
pub const GMEARTH: f64 =  398600441.8e6;
pub const GMJUPITER: f64 =  1.26686534e17;
pub const AU: f64 = 149597870.7e3;
pub const LD: f64 = 384402e3;
pub const DAY: f64 = 36e2*24.;
pub const PC: f64 = AU*648e3/PI;
pub const LY: f64 = 365.25*C*DAY;
pub const GG: f64 = C*HBAR/(MP*MP);
pub const MS: f64 = GMSUN/GG;
pub const JS: f64 = MS*AU*AU/(DAY*DAY);
pub const IAU: UnitSystem = UnitSystem{kb: RU*ME/MEU/0.001/JS, hbar: HBAR/DAY/JS, c: DAY*C/AU, m0: 4.*PI*1e-7*DAY*DAY/JS, me: ME/MS, mu: 0.001, lambda: 1., al: 1.};

// aliased & humorous units

pub const MF: f64 = 9e1*SLUG/LBM;
pub const JF: f64 = MF*(201.168/14./DAY)*(201.168/14./DAY);
pub const FFF: UnitSystem = UnitSystem{kb: 1e3*RU*ME/MEU*RANKINE/JF, hbar: HBAR/14./DAY/JF, c: 14.*DAY*C/201.168, m0: 0., me: ME/MF, mu: 1., lambda: 1., al: 1.};

// pub const units, US, temp = UnitSystem, UnitSystem, temperature
pub const SI: UnitSystem = SI2019;
pub const MKS: UnitSystem = METRIC;
pub const CGS: UnitSystem = GAUSS;
pub const CGS2019: UnitSystem = EMU2019;
pub const CGSM: UnitSystem = EMU;
pub const CGSE: UnitSystem = ESU;
pub const HLU: UnitSystem = LORENTZHEAVISIDE;

// natural units

pub const AG: f64 = (ME/MP)*(ME/MP);
pub const SQRT_PI: f64 = 2./FRAC_2_SQRT_PI;
pub const SQRT_AINV: f64 = 11.706237614366112;
pub const PLANCK: UnitSystem = UnitSystem{kb: 1., hbar: 1., c: 1., m0: 1., me: 2.*SQRT_PI*(ME/MP), mu: 1., lambda: 1., al: 1.};
pub const PLANCKGAUSS: UnitSystem = UnitSystem{kb: 1., hbar: 1., c: 1., m0: 4.*PI, me: ME/MP, mu: 1., lambda: 1., al: 1.};
pub const STONEY: UnitSystem = UnitSystem{kb: 1., hbar: AINV, c: 1., m0: 4.*PI, me: (ME/MP)*SQRT_AINV, mu: 1., lambda: 1., al: 1.};
pub const HARTREE: UnitSystem = UnitSystem{kb: 1., hbar: 1., c: AINV, m0: 4.*PI/(AINV*AINV), me: 1., mu: 1., lambda: 1., al: 1.};
pub const RYDBERG: UnitSystem = UnitSystem{kb: 1., hbar: 1., c: 2.*AINV, m0: PI/(AINV*AINV), me: 0.5, mu: 1., lambda: 1., al: 1.};
pub const SCHRODINGER: UnitSystem = UnitSystem{kb: 1., hbar: 1., c: AINV, m0: 4.*PI/(AINV*AINV), me: (ME/MP)*SQRT_AINV, mu: 1., lambda: 1., al: 1.};
pub const ELECTRONIC: UnitSystem = UnitSystem{kb: 1., hbar: AINV, c: 1., m0: 4.*PI, me: 1., mu: 1., lambda: 1., al: 1.};
pub const NATURAL: UnitSystem = UnitSystem{kb: 1., hbar: 1., c: 1., m0: 1., me: 1., mu: 1., lambda: 1., al: 1.};
pub const NATURALGAUSS: UnitSystem = UnitSystem{kb: 1., hbar: 1., c: 1., m0: 4.*PI, me: 1., mu: 1., lambda: 1., al: 1.};
pub const QCD: UnitSystem = UnitSystem{kb: 1., hbar: 1., c: 1., m0: 1., me: 1./MPE, mu: 1., lambda: 1., al: 1.};
pub const QCDGAUSS: UnitSystem = UnitSystem{kb: 1., hbar: 1., c: 1., m0: 4.*PI, me: 1./MPE, mu: 1., lambda: 1., al: 1.};
pub const QCDORIGINAL: UnitSystem = UnitSystem{kb: 1., hbar: 1., c: 1., m0: 4.*PI/AINV, me: 1./MPE, mu: 1., lambda: 1., al: 1.};

// constants

pub const fn boltzmann(u: &UnitSystem) -> f64 {return u.kb}
pub const fn planckreduced(u: &UnitSystem) -> f64 {return u.hbar}
pub const fn lightspeed(u: &UnitSystem) -> f64 {return u.c}
pub const fn vacuumpermeability(u: &UnitSystem) -> f64 {return u.m0}
pub const fn electronmass(u: &UnitSystem) -> f64 {return u.me}
pub const fn molarmass(u: &UnitSystem) -> f64 {return u.mu}
pub const fn rationalization(u: &UnitSystem) -> f64 {return u.lambda}
pub const fn lorentz(u: &UnitSystem) -> f64 {return u.al}

pub fn hyperfine(u: &UnitSystem) -> f64 {
    return DVCS/frequency(u,&SI2019)
}

pub fn luminousefficacy(u: &UnitSystem) -> f64 {
    return KCD*power(u,&SI2019)
}

pub fn avogadro(u: &UnitSystem) -> f64 {
    return MEU*molarmass(u)/electronmass(u)
}

pub fn planck(u: &UnitSystem) -> f64 {
    return 2.*PI*planckreduced(u)
}

pub fn atomicmass(u: &UnitSystem) -> f64 {
    return electronmass(u)/MEU
}

pub fn protonmass(u: &UnitSystem) -> f64 {
    return MPE*electronmass(u)
}

pub fn planckmass(u: &UnitSystem) -> f64 {
    return MP/mass(u,&SI2019)
}

pub fn newton(u: &UnitSystem) -> f64 {
    return lightspeed(u)*planckreduced(u)/planckmass(u).powi(2)
}

pub fn einstein(u: &UnitSystem) -> f64 {
    return 8.*PI*newton(u)/lightspeed(u).powi(4)
}

pub fn universal(u: &UnitSystem) -> f64 {
    return boltzmann(u)*avogadro(u)
}

pub fn stefan(u: &UnitSystem) -> f64 {
    return 2.*PI.powi(5)*boltzmann(u).powi(4)/(15.*planck(u).powi(3)*lightspeed(u).powi(2))
}

pub fn radiationdensity(u: &UnitSystem) -> f64 {
    return 4.*stefan(u)/lightspeed(u)
}

pub fn vacuumpermittivity(u: &UnitSystem) -> f64 {
    return 1./(vacuumpermeability(u)*(lightspeed(u)*lorentz(u)).powi(2))
}

pub fn coulomb(u: &UnitSystem) -> f64 {
    return rationalization(u)/(4.*PI*vacuumpermittivity(u))
}

pub fn biotsavart(u: &UnitSystem) -> f64 {
    return vacuumpermeability(u)*lorentz(u)*(rationalization(u)/(4.*PI))
}

pub fn ampere(u: &UnitSystem) -> f64 {
    return lorentz(u)*biotsavart(u)
}

pub fn vacuumimpedance(u: &UnitSystem) -> f64 {
    return vacuumpermeability(u)*lightspeed(u)*rationalization(u)*lorentz(u).powi(2)
}

pub fn elementarycharge(u: &UnitSystem) -> f64 {
    return (2.*planck(u)/vacuumimpedance(u)/AINV).sqrt()
}

pub fn faraday(u: &UnitSystem) -> f64 {
    return elementarycharge(u)*avogadro(u)
}

pub fn josephson(u: &UnitSystem) -> f64 {
    return 2.*elementarycharge(u)*lorentz(u)/planck(u)
}

pub fn magneticfluxquantum(u: &UnitSystem) -> f64 {
    return 1./josephson(u)
}

pub fn klitzing(u: &UnitSystem) -> f64 {
    return planck(u)/elementarycharge(u).powi(2)
}

pub fn conductancequantum(u: &UnitSystem) -> f64 {
    return 2.*elementarycharge(u).powi(2)/planck(u)
}

pub fn hartree(u: &UnitSystem) -> f64 {
    return electronmass(u)*(lightspeed(u)/AINV).powi(2)
}

pub fn rydberg(u: &UnitSystem) -> f64 {
    return hartree(u)/(2.*planck(u)*lightspeed(u))
}

pub fn bohr(u: &UnitSystem) -> f64 {
    return AINV*planckreduced(u)/(electronmass(u)*lightspeed(u))
}

pub fn bohrreduced(u: &UnitSystem) -> f64 {
    return bohr(u)*(1.+1./MPE)
}

pub fn electronradius(u: &UnitSystem) -> f64 {
    return planckreduced(u)/(electronmass(u)*lightspeed(u)*AINV)
}

pub fn magneton(u: &UnitSystem) -> f64 {
    return elementarycharge(u)*planckreduced(u)*lorentz(u)/(2.*electronmass(u))
}

// conversions

pub fn length(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return (planckreduced(s)*electronmass(u)*lightspeed(u))/(planckreduced(u)*electronmass(s)*lightspeed(s))
}

pub fn area(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return length(u,s).powi(2)
}

pub fn volume(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return length(u,s).powi(3)
}

pub fn wavenumber(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return length(s,u)
}

pub fn fuelefficiency(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return area(s,u)
}

pub fn time(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return length(u,s)/speed(u,s)
}

pub fn frequency(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return time(s,u)
}

pub fn frequencydrift(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return time(s,u).powi(2)
}

pub fn speed(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return lightspeed(s)/lightspeed(u)
}

pub fn acceleration(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return speed(u,s)/time(u,s)
}

pub fn jerk(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return speed(u,s)/time(u,s).powi(2)
}

pub fn snap(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return speed(u,s)/time(u,s).powi(3)
}

pub fn volumeflow(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return area(u,s)*speed(u,s)
}

pub fn specificenergy(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return speed(u,s).powi(2)
}

// kinematic

pub fn mass(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return electronmass(s)/electronmass(u)
}

pub fn energy(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)*specificenergy(u,s)
}

pub fn power(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return energy(u,s)/time(u,s)
}

pub fn force(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)*acceleration(u,s)
}

pub fn pressure(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)/length(u,s)/time(u,s).powi(2)
}

// mechanical

pub fn momentum(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)*speed(u,s)
}

pub fn angularmomentum(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return momentum(u,s)*length(u,s)
}

pub fn yank(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)*jerk(u,s)
}

pub fn areadensity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)/area(u,s)
}

pub fn density(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)/volume(u,s)
}

pub fn specificvolume(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return volume(u,s)/mass(u,s)
}

pub fn action(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return momentum(u,s)*length(u,s)
}

pub fn stiffness(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return energy(u,s)/area(u,s)
}

pub fn intensity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return power(u,s)/area(u,s)
}

pub fn diffusivity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return (planck(s)/planck(u))/mass(u,s)
}

pub fn viscosity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)/length(u,s)/time(u,s)
}

pub fn lineardensity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)/length(u,s)
}

pub fn massflow(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)/time(u,s)
}

pub fn spectralflux(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return power(u,s)/length(u,s)
}

pub fn powerdensity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return power(u,s)/volume(u,s)
}

pub fn compressibility(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return pressure(s,u)
}

pub fn fluence(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return energy(u,s)/area(u,s)
}

pub fn rotationalinertia(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)/area(u,s)
}

// acoustic

pub fn soundexposure(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return time(u,s)*pressure(u,s).powi(2)
}

pub fn specificimpedance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return pressure(u,s)/speed(u,s)
}

pub fn acousticimpedance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return specificimpedance(u,s)/area(u,s)
}

pub fn admittance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return area(u,s)/specificimpedance(u,s)
}

pub fn compliance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return time(u,s).powi(2)/mass(u,s)
}

pub fn inertance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)/length(u,s).powi(4)
}

// electromagnetic

pub fn charge(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return ((planckreduced(s)*vacuumpermeability(u)*lightspeed(u)*rationalization(u)*lorentz(u).powi(2))/(planckreduced(u)*vacuumpermeability(s)*lightspeed(s)*rationalization(s)*lorentz(s).powi(2))).sqrt()
}

pub fn current(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return charge(u,s)/time(u,s)
}

pub fn electricpotential(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return energy(u,s)/charge(u,s)
}

pub fn capacitance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return charge(u,s)/electricpotential(u,s)
}

pub fn resistance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return electricpotential(u,s)/current(u,s)
}

pub fn conductance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return current(u,s)/electricpotential(u,s)
}

pub fn magneticflux(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return energy(u,s)/(lorentz(s)/lorentz(u))/current(u,s)
}

pub fn magneticfluxdensity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)/(lorentz(s)/lorentz(u))/current(u,s)/time(u,s).powi(2)
}

pub fn inductance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)*area(u,s)/charge(u,s)
}

// electromagnetics

pub fn electricfluxdensity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return charge(u,s)*(rationalization(s)/rationalization(u))/area(u,s)
}

pub fn chargedensity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return charge(u,s)/volume(u,s)
}

pub fn currentdensity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return current(u,s)/area(u,s)
}

pub fn conductivity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return conductance(u,s)/length(u,s)
}

pub fn permittivity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return capacitance(u,s)*(rationalization(s)/rationalization(s))/length(u,s)
}

pub fn permeability(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return vacuumpermeability(s)/vacuumpermeability(u)
}

pub fn electricfield(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return electricpotential(u,s)/length(u,s)
}

pub fn magneticfield(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return current(u,s)*(rationalization(s)/rationalization(u))*(lorentz(s)/lorentz(u))/length(u,s)
}

pub fn exposure(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return charge(u,s)/mass(u,s)
}

pub fn resistivity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return resistance(u,s)*length(u,s)
}

pub fn linearchargedensity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return charge(u,s)/length(u,s)
}

pub fn magneticdipolemoment(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return current(u,s)*(lorentz(s)/lorentz(u))/area(u,s)
}

pub fn mobility(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return charge(u,s)*time(u,s)/mass(u,s)
}

pub fn reluctance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return (rationalization(s)/rationalization(u))*(lorentz(s)*lorentz(u)).powi(2)/inductance(u,s)
}

pub fn vectorpotential(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return magneticflux(u,s)/length(u,s)
}

pub fn magneticmoment(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return magneticflux(u,s)*length(u,s)
}

pub fn rigidity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return magneticfluxdensity(u,s)*length(u,s)
}

pub fn susceptibility(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return rationalization(u)/rationalization(s)
}

pub fn electricflux(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return electricpotential(u,s)*length(u,s)
}

pub fn electricdipolemoment(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return charge(u,s)*length(u,s)
}

pub fn magneticpotential(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return magneticflux(u,s)*reluctance(u,s)
}

pub fn polestrength(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return magneticdipolemoment(u,s)/length(u,s)
}

pub fn permeance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return reluctance(s,u)
}

pub fn specificsusceptibility(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return magneticdipolemoment(u,s)/magneticfield(u,s)/mass(u,s)
}

pub fn magnetizability(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return magneticmoment(u,s)/magneticfluxdensity(u,s)
}

pub fn electricpolarizability(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return electricdipolemoment(u,s)/electricfield(u,s)
}

pub fn magneticpolarizability(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return magneticdipolemoment(u,s)/magneticfield(u,s)
}

pub fn magnetization(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return magneticmoment(u,s)/volume(u,s)
}

pub fn specificmagnetization(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return magneticmoment(u,s)/mass(u,s)
}

pub fn demagnitizingfactor(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return susceptibility(s,u)
}

// magneticfluxrate
// magneticfluxratedensity
// magneticresistance

// thermodynamics

pub fn temperature(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return (boltzmann(u)*electronmass(s)*lightspeed(s).powi(2))/(boltzmann(s)*electronmass(u)*lightspeed(u).powi(2))
}

pub fn entropy(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return energy(u,s)/temperature(u,s)
}

pub fn specificentropy(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return specificenergy(u,s)/temperature(u,s)
}

pub fn volumeheatcapacity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return entropy(u,s)/volume(u,s)
}

pub fn thermalconductivity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return force(u,s)/time(u,s)/temperature(u,s)
}

pub fn thermalconductance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return thermalconductivity(u,s)*length(u,s)
}

pub fn thermalresistance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return thermalconductance(s,u)
}

pub fn thermalexpansion(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return temperature(s,u)
}

pub fn lapserate(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return temperature(u,s)/length(u,s)
}

// molar

pub fn molality(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return molarmass(u)/molarmass(s)
}

pub fn mole(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mass(u,s)*molality(u,s)
}

pub fn molarity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mole(u,s)/volume(u,s)
}

pub fn molarvolume(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return volume(u,s)/mole(u,s)
}

pub fn molarentropy(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return entropy(u,s)/mole(u,s)
}

pub fn molarenergy(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return energy(u,s)/mole(u,s)
}

pub fn molarconductivity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return conductivity(u,s)*area(u,s)/mole(u,s)
}

pub fn molarsusceptibility(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return specificsusceptibility(u,s)*molality(s,u)
}

pub fn catalysis(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return mole(u,s)/time(u,s)
}

pub fn specificity(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return volume(u,s)/mole(u,s)/time(u,s)
}

// photometrics

pub fn luminousflux(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return ((luminousefficacy(s)*planckreduced(s))/(luminousefficacy(u)*planckreduced(u)))*frequency(u,s).powi(2)
}

pub fn luminance(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return luminousflux(u,s)/area(u,s)
}

pub fn luminousenergy(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return luminousflux(u,s)*time(u,s)
}

pub fn luminousexposure(u: &UnitSystem, s: &UnitSystem) -> f64 {
    return luminance(u,s)*time(u,s)
}
