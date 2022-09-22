//
//  Catalogs.swift
//

/// Returns the mass (in Solar Masses) given a relative "Earth Mass" value (1.0 = Earth)
/// - Parameter x: The earth mass to convert
func EM(_ x: Double) -> Double {
    x / SUN_MASS_IN_EARTH_MASSES
}

func AVG(_ x: Double, _ y: Double) -> Double {
    (x + y) / 2.0
}

// Earth Moons
let luna = Planet(planet_no: 1, a: 2.571e-3, e: 0.055, axial_tilt: 1.53, mass: EM(0.01229), gas_giant: false,
                  dust_mass: EM(0.01229), gas_mass: 0, next_planet: nil)

// Jupiter Moons
let callisto = Planet(planet_no: 4, a: 1.259e-2, e: 0, axial_tilt: 0, mass: EM(1.62e-2), gas_giant: false, dust_mass: EM(1.62e-2), gas_mass: 0, next_planet: nil)
let ganymede = Planet(planet_no: 3, a: 7.16e-3, e: 0.0796, axial_tilt: 0, mass: EM(2.6e-2), gas_giant: false, dust_mass: EM(2.6e-2), gas_mass: 0, next_planet: callisto)
let europa = Planet(planet_no: 2, a: 4.49e-3, e: 0.0075, axial_tilt: 0, mass: EM(7.9e-3), gas_giant: false, dust_mass: EM(7.9e-3), gas_mass: 0, next_planet: ganymede)
let io = Planet(planet_no: 1, a: 2.82e-3, e: 0.0006, axial_tilt: 0, mass: EM(1.21e-2), gas_giant: false, dust_mass: EM(1.21e-2), gas_mass: 0, next_planet: europa)

// Saturn Moons

let iapetus = Planet(planet_no: 6, a: 2.38e-2, e: 0.029, axial_tilt: 0, mass: EM(8.4e-4), gas_giant: false, dust_mass: EM(8.4e-4), gas_mass: 0, next_planet: nil)
let hyperion = Planet(planet_no: 5, a: 9.89e-3, e: 0.110, axial_tilt: 0, mass: EM(1.82e-5), gas_giant: false, dust_mass: EM(1.82e-5), gas_mass: 0, next_planet: iapetus)
let titan = Planet(planet_no: 4, a: 8.17e-3, e: 0.0289, axial_tilt: 0, mass: EM(2.3e-2), gas_giant: false, dust_mass: EM(2.3e-2), gas_mass: 0, next_planet: hyperion)
let rhea = Planet(planet_no: 3, a: 3.52e-3, e: 0.0009, axial_tilt: 0, mass: EM(3.85e-4), gas_giant: false, dust_mass: EM(3.85e-4), gas_mass: 0, next_planet: titan)
let dione = Planet(planet_no: 2, a: 2.52e-3, e: 0.0021, axial_tilt: 0, mass: EM(1.74e-4), gas_giant: false, dust_mass: EM(1.74e-4), gas_mass: 0, next_planet: rhea)
let tethys = Planet(planet_no: 1, a: 1.97e-3, e: 0, axial_tilt: 0, mass: EM(1.09e-4), gas_giant: false, dust_mass: EM(1.09e-4), gas_mass: 0, next_planet: dione)

let triton = Planet(planet_no: 1, a: 2.36e-3, e: 0, axial_tilt: 0, mass: EM(2.31e-2), gas_giant: false, dust_mass: EM(2.31e-2), gas_mass: 0, next_planet: nil)
let charon = Planet(planet_no: 1, a: 19571 / KM_PER_AU, e: 0, axial_tilt: 0, mass: EM(2.54e-4), gas_giant: false, dust_mass: EM(2.54e-4), gas_mass: 0, next_planet: nil)

let xena = Planet(planet_no: 11, a: 67.6681, e: 0.44177, axial_tilt: 0, mass: EM(0.0025), gas_giant: false, dust_mass: EM(0.0025), gas_mass: 0, next_planet: nil)
let pluto = Planet(planet_no: 10, a: 39.529, e: 0.248, axial_tilt: 122.5, mass: EM(0.002), gas_giant: false, dust_mass: EM(0.002), gas_mass: 0, next_planet: xena)
let neptune = Planet(planet_no: 9, a: 30.061, e: 0.010, axial_tilt: 29.6, mass: EM(17.14), gas_giant: true, dust_mass: 0, gas_mass: EM(17.14), first_moon: triton, next_planet: pluto)
let uranus = Planet(planet_no: 8, a: 19.191, e: 0.046, axial_tilt: 97.9, mass: EM(14.530), gas_giant: true, dust_mass: 0, gas_mass: EM(14.530), next_planet: neptune)
let saturn = Planet(planet_no: 7, a: 9.539, e: 0.056, axial_tilt: 26.7, mass: EM(95.18), gas_giant: true, dust_mass: 0, gas_mass: EM(95.18), first_moon: tethys, next_planet: uranus)
let jupiter = Planet(planet_no: 6, a: 5.203, e: 0.048, axial_tilt: 3.1, mass: EM(317.9), gas_giant: true, dust_mass: 0, gas_mass: EM(317.9), first_moon: io, next_planet: saturn)
let ceres = Planet(planet_no: 5, a: 2.766, e: 0.08, axial_tilt: 0, mass: 9.5e20 / SOLAR_MASS_IN_KILOGRAMS, gas_giant: false, dust_mass: 9.5e20 / SOLAR_MASS_IN_KILOGRAMS, gas_mass: 0, next_planet: jupiter)
let mars = Planet(planet_no: 4, a: 1.524, e: 0.093, axial_tilt: 25.2, mass: EM(0.1074), gas_giant: false, dust_mass: EM(0.1074), gas_mass: 0, next_planet: ceres)
let earth = Planet(planet_no: 3, a: 1.0, e: 0.017, axial_tilt: 23.5, mass: EM(1.0), gas_giant: false, dust_mass: EM(1.0), gas_mass: 0, first_moon: luna, next_planet: mars)
let venus = Planet(planet_no: 2, a: 0.723, e: 0.007, axial_tilt: 177.3, mass: EM(0.815), gas_giant: false, dust_mass: EM(0.815), gas_mass: 0, next_planet: earth)
let mercury = Planet(planet_no: 1, a: 0.387, e: 0.206, axial_tilt: 2, mass: EM(0.055), gas_giant: false, dust_mass: EM(0.055), gas_mass: 0, next_planet: venus)
// planet_pointer solar_system = &mercury;

let eriEpsI = Planet(planet_no: 1, a: 3.3, e: 0.608, axial_tilt: 0, mass: 0, gas_giant: false, dust_mass: 0, gas_mass: 0, next_planet: nil)

let UMa47II = Planet(planet_no: 2, a: 3.73, e: 0, axial_tilt: 0, mass: 0, gas_giant: false, dust_mass: 9, gas_mass: 9, next_planet: nil)
let UMa47I = Planet(planet_no: 1, a: 2.11, e: 0.096, axial_tilt: 0, mass: 0, gas_giant: false, dust_mass: 0, gas_mass: 0, next_planet: UMa47II)
let horIotI = Planet(planet_no: 1, a: 0.925, e: 0.161, axial_tilt: 0, mass: 0, gas_giant: false, dust_mass: 0, gas_mass: 0, next_planet: nil)

//
//            /*    No Orbit Eccen. Tilt   Mass    Gas Giant? Dust Mass   Gas */
// planets    smallest={0, 0.0, 0.0,    20.0,    EM(0.4),   FALSE,  EM(0.4),   0.0, ZEROES,0,NULL, NULL};
// planets    average    ={0, 0.0, 0.0,    20.0,    EM(1.0),   FALSE,  EM(1.0),    0.0, ZEROES,0,NULL, NULL};
// planets    largest    ={0, 0.0, 0.0,    20.0,    EM(1.6),   FALSE,  EM(1.6),   0.0, ZEROES,0,NULL, NULL};
let perdole: [Star] = [
    Star(0, 1.00, 0, 0, 0, mercury, "Sol", 1, "The Solar System"),
    Star(0, 1.08, 0.88, 0.52, 23.2, nil, "ALF Cen A", 1, "Alpha Centauri A"),
    Star(0, 0.88, 1.08, 0.52, 23.2, nil, "ALF Cen B", 1, "Alpha Centauri B"),
    Star(0, 0.80, 0, 0, 0, eriEpsI, "EPS Eri", 1, "Epsilon Eridani"),
    Star(0, 0.82, 0, 0, 0, nil, "TAU Cet", 1, "Tau Ceti"),
    Star(0, 0.90, 0.65, 0.50, AVG(22.8, 24.3), nil, "70 Oph", 1, "70 Ophiuchi A"),
    Star(0, 0.94, 0.58, 0.53, AVG(69.0, 71.0), nil, "ETA Cas", 1, "Eta Cassiopeiae A"),
    Star(0, 0.82, 0, 0, 0, nil, "SIG Dra", 1, "Sigma Draconis"),
    Star(0, 0.77, 0.76, 0, 22.0, nil, "36 Oph", 1, "36 Ophiuchi A"),
    Star(0, 0.76, 0.77, 0, 22.0, nil, "36 Oph B", 0, "36 Ophiuchi B"),
    /* Fake up a B just to clip the distances -- need the real data */
    Star(0, 0.76, 0, 0, 46.0, nil, "HD 191408", 1, "HR7703 A"),
    Star(0, 0.76, 0.5, 0.5, 46.0, nil, "HD 191408", 1, "HR7703 A"),
    Star(0, 0.98, 0, 0, 0, nil, "DEL Pav", 1, "Delta Pavonis"),
    Star(0, 0.91, 0, 0, 0, nil, "82 Eri", 1, "82 Eridani"),
    Star(0, 1.23, 0, 0, 0, nil, "BET Hyi", 1, "Beta Hydri"),
    Star(0, 0.74, 0, 0, 0, nil, "HD 219134", 1, "HR8832"),
    Star(0, 0.725, 0, 0, 1100.0, nil, "HD 16160", 1, "HR753 A"),
]

/*
 *    The following values were taken from: http://www.solstation.com/stars.htm
 */
let web: [Star] = [
    //// L            Mass            Mass2            Eccen.    SMAxis     Planets    Designation    Name
    Star(1.00, 1.00, 0, 0, 0, mercury, "Sol", 1, "The Solar System"), // 0
    Star(1.60, 1.09, 0.90, 0.519, 23.7, nil, "ALF Cen A", 1, "Alpha Centauri A"), // 4.4
    Star(0.45, 0.90, 1.09, 0.519, 23.7, nil, "ALF Cen B", 1, "Alpha Centauri B"), // 4.4
    Star(0.34, 0.85, 0, 0, 0, eriEpsI, "EPS Eri", 1, "Epsilon Eridani"), // 10.5
    // Star(AVE(6.3,8.9),0.59,            0.5,            .48,    85.2,     nil,        "61 Cyg A",     1, "61 Cygni A"),                // 11.4
    Star(0.085, 0.59, 0.5, 0.48, 85.2, nil, "61 Cyg A", 1, "61 Cygni A"), // 11.4
    Star(0.59, 0.82, 0, 0, 0, nil, "TAU Cet", 1, "Tau Ceti"), // 11.9
    Star(0.38, 0.75, 0.16 + 0.43, 0, 418.0, nil, "40 Eri", 1, "40 Eridani A"), // 16.5
    Star(AVG(0.44, 0.47), 0.924, 0.701, 0.495, 23.3, nil, "70 Oph", 1, "70 Ophiuchi A"), // 16.6
    Star(0.39, 0.82, 0, 0, 0, nil, "SIG Dra", 1, "Sigma Draconis"), // 18.8
    Star(0.156, 0.76, 0.55 + 0.35, 0.20, 190.0, nil, "33 g Lib", 1, "HR 5568"), // 19.3
    Star(AVG(1.0, 1.29), 0.91, 0.56, 0.497, 71.0, nil, "ETA Cas", 1, "Eta Cassiopeiae A"), // 19.4
    Star(0.23, 0.82, 0.20, 0.5, 43.0, nil, "HD 191408", 1, "HR 7703 (HJ 5173) A"), // 19.7
    Star(0.65, 0.97, 0, 0, 0, nil, "82 Eri", 1, "82 Eridani"), // 19.8
    Star(1.2, 0.98, 0, 0, 0, nil, "DEL Pav", 1, "Delta Pavonis"), // 19.9
    Star(0, 0.74, 0, 0, 0, nil, "HD 219134", 1, "HR 8832"), // 21.3
    Star(0.52, 0.90, 0.76, 0.51, 33.6, nil, "XI Boo", 1, "Xi Bootis A"), // 21.8
    Star(0.21, 0.81, 0.082, 0.75, 15.0, nil, "HD 16160", 1, "HR 753 A"), // 23.5
    Star(0.24, 0.83, 0, 0, 0, nil, "HD 4628", 1, "BD+04 123 (HR 222)"), // 24.3
    Star(3.6, 1.1, 0, 0, 0, nil, "BET Hyi", 1, "Beta Hydri"), // 24.4
    Star(0.37, 0.89, 0, 0, 0, nil, "107 Psc", 1, "107 Piscium"), // 24.4
    // 107 Psc p1 = Klotho in Celestia's imagined.ssc
    Star(3.0, 1.3, 0, 0, 0, nil, "PI3 Ori", 1, "Pi3 Orionis A"), // 26.2
    Star(0.28, 0.88, 0.86, 0.534, 63.7, nil, "RHO1 Eri", 1, "Rho Eridani A"), // 26.6
    Star(0.25, 0.86, 0.88, 0.534, 63.7, nil, "RHO2 Eri", 1, "Rho Eridani B"), // 26.6
    Star(1.2, 1.07, 0, 0, 0, nil, "BET CVn", 1, "Chara"), // 27.3
    Star(2.9, 0.90, 1.45, 0.412, 21.2, nil, "XI UMa", 1, "Xi Ursae Majoris Ba"), // 27.3
//                                                                    Xi Urs Maj aka Alula Australis
//                    55203:Alula Australis:XI UMa:53 UMa defined in Celestia starnames, but no data
    Star(0.80, 0.96, 0, 0, 0, nil, "61 Vir", 1, "61 Virginis"), // 27.8
    Star(1.3, 0.98, 0, 0, 0, nil, "ZET Tuc", 1, "Zeta Tucanae"), // 28.0
    Star(1.08, 1.0, 0.15, 0.45, 6.4, nil, "CHI1 Ori", 1, "Chi1 Orionis A"), // 28.3
//                    41 Arae masses are Wieth-Knudsen's 1957 estimates,
    Star(0.41, 0.9, 0.6, 0.779, 91.5, nil, "41 Ari", 1, "41 Arae A"), // 28.7
    Star(0.21, 0.845, 0, 0, 0, nil, "HR 1614", 0, "BD-05 1123 (HR 1614) A"), // 28.8
    Star(0.33, 0.87, 0, 0, 0, nil, "HR 7722", 0, "CD-27 14659 (HR 7722)"), // 28.8
    Star(2.6, 1.2, 0.63, 0.5, 864.0, nil, "GAM Lep", 1, "Gamma Leporis A"), // 29.3
    Star(1.4, 1.05, 0, 0, 0, nil, "BET Com", 1, "Beta Comae Berenices"), // 29.9
    Star(0.85, 1.0, 0, 0, 0, nil, "KAP1 Cet", 1, "Kappa Ceti"), // 29.9
    Star(1.5, 0.8, 0, 0, 0, nil, "GAM Pav", 1, "Gamma Pavonis"), // 30.1
    Star(0.82, 0.8, 0.07, 0.6, 235.0, nil, "HD 102365", 1, "HR 4523"), // 30.1
    Star(0.588, 0.81, 0, 0, 0, nil, "61 UMa", 1, "61 Ursae Majoris"), // 31.1
    Star(0.31, 0.87, 0, 0.5, 80.5, nil, "HR 4458", 0, "CD-32 8179 (HR 4458)"), // 31.1
    Star(AVG(0.39, 0.41), 0.90, 0, 0, 0, nil, "12 Oph", 1, "12 Ophiuchi"), // 31.9
    Star(0.46, 0.92, 0, 0, 0, nil, "HR 511", 0, "BD+63 238 (HR 511)"), // 32.5
    Star(0.83, 0.87, 0, 0, 0, nil, "ALF Men", 1, "Alpha Mensae"), // 33.1
    Star(0.93, 0.79, 1.02, 0.5, 9000.0, nil, "ZET1 Ret", 1, "Zeta 1 Reticuli"), // 39.4-39.5
    Star(0.99, 1.02, 0.79, 0.5, 9000.0, nil, "ZET2 Ret", 1, "Zeta 2 Reticuli"), // 39.4-39.5
    Star(1.14, 1.05, 2.0, 0.55, 48.5, nil, "44 Boo", 1, "44 Bootis A"), // 41.6
    Star(1.7, 1.03, 0, 0, 0, UMa47I, "47 UMa", 1, "47 Ursae Majoris"), // 45.9
    Star(1.8, 1.03, 0, 0, 0, horIotI, "IOT Hor", 1, "Iota Horologii"), // 56.2
    Star(AVG(0.13, 0.15), AVG(0.59, 0.71), 0, 0, 0, nil, "EPS Ind", 1, "Epsilon Indi"), // 11.8
    Star(AVG(0.083, 0.09), 0.701, 0.924, 0.495, 23.3, nil, "70 Oph", 1, "70 Ophiuchi B"), // 16.6
    Star(0.28, 0.85, 0.85, 0.922, 88.0, nil, "36 Oph", 1, "36 Ophiuchi A"), // 19.5
    Star(0.27, 0.85, 0.85, 0.922, 88.0, nil, "36 Oph B", 0, "36 Ophiuchi B"), // 19.5
    Star(0.12, 0.75, 0.65, 0.58, 12.6, nil, "HR 6426", 0, "MLO 4 (HR 6426) A"), // 22.7
    Star(0.146, 0.80, 0.50, 0.5, 500.0, nil, "BD-05 1844 A", 0, "BD-05 1844 A"), // 28.3
]

// BD-05 1123 A:     HR 1614, Gl 183 A, Hip 23311, HD 32147, SAO 131688, LHS 200, LTT 2412, LFT 382, and LPM 200.
// CD-27 14659:         HR 7722, Gl 785, Hip 99825, HD 192310, CP(D)-27 6972, SAO 189065, LHS 488, LTT 8009, LFT 1535, and LPM 731
// CD-32 8179 A:     HR 4458, Gl 432 A, Hip 56452, HD 100623, CP(D)-32 3122, SAO 202583, LHS 308, LTT 4280, LFT 823, LPM 389, and E 439-246.
// BD+63 238:         HR 511*, Gl 75, Hip 8362, HD 10780, SAO 11983, LHS 1297, LTT 10619, and LFT 162.
// 36 Ophiuchi B:     HR 6401, Gl 663 B, HD 155885, SAO 185199, LHS 438, and ADS 10417 B.
// MLO 4 A:             HR 6426*, Gl 667 A, Hip 84709, HD 156384, CD-34 11626 A, CP-34 6803, SAO 208670, LHS 442, LTT 6888, LFT 1336, LPM 638, and UGPMF 433.
// BD-05 1844 A:     Gl 250 A, Hip 32984, HD 50281, SAO 133805, LHS 1875, LTT 2662, LFT 494, and LPM 244.
//
// {.00006,        0.105,            0.1,            0.62,    5.5,     NULL,        "",    "Luyten 726-8 A"},        // 8.7
// {0.039,        0.5,            0.59,            .48,    85.2,     NULL,        "",    "61 Cygni B"},            // 11.4
// {0.05,        0.65,            0.75,            0.58,    12.6,     NULL,        "",    "MLO 4 (HR 6426) B"},    // 22.7
// {1.1,        1.05,            0.4,            0.53,    0.6,     NULL,        "",    "Xi Ursae Majoris Aa"},    // 27.3
// {0,            0.4,            1.05,            0.53,    0.6,     NULL,        "",    "Xi Ursae Majoris Ab"},    // 27.3
// {0.064,        0.76,            0.90,            0.51,    33.0,     NULL,        "",    "Xi Bootis B"},            // 21.8
//
// star    various[] =
// {
//// L            Mass            Mass2            Eccen.    SMAxis     Planets    Designation    Name
// {1.00,            1.00,            0,                0,        0,         &mercury,    "Sol",         1, "The Solar System"},        // 0
// {14800.,        8,                0,                0,        0,         NULL,        "ALF Car",     1, "Canopus"}
// };
//
let solstation = Catalog(count: web.count, arg: "w", stars: web)
let dole = Catalog(count: perdole.count, arg: "d", stars: perdole)
// catalog jimb        = {sizeof(various) / sizeof (star), "F",    &various};
