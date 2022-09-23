//
//  Types.swift
//

// #pragma once
//
// typedef struct dust_record    *dust_pointer;
// typedef struct planets_record  *planet_pointer;
// typedef struct gen *gen_pointer;
import Foundation

enum PlanetType {
    case unknown
    case rock
    case venusian
    case terrestrial
    case gasgiant
    case martian
    case water
    case ice
    case subgasgiant
    case subsubgasgiant
    case asteroids
    case face
}

// typedef enum planet_type {
//    tUnknown,
//    tRock,
//    tVenusian,
//    tTerrestrial,
//    tGasGiant,
//    tMartian,
//    tWater,
//    tIce,
//    tSubGasGiant,
//    tSubSubGasGiant,
//    tAsteroids,
//    t1Face
// } planet_type;

final class Gas {
    var type: ChemTableEntry
    var surf_pressure: Double /* units of millibars (mb)             */

    init(_ type: ChemTableEntry, surf_pressure: Double) {
        self.type = type
        self.surf_pressure = surf_pressure
    }
}

// typedef struct gas {
//      int         num;
//    long double    surf_pressure;        /* units of millibars (mb)             */
//    } gas;

final class Sun {
    var luminosity: Double
    var mass: Double
    var life: Double
    var age: Double
    var r_ecosphere: Double
    var name: String

    init(luminosity: Double, mass: Double, life: Double, age: Double, r_ecosphere: Double, name: String) {
        self.luminosity = luminosity
        self.mass = mass
        self.life = life
        self.age = age
        self.r_ecosphere = r_ecosphere
        self.name = name
    }
}

// typedef struct sun {
//    long double    luminosity;
//    long double    mass;
//    long double life;
//    long double age;
//    long double r_ecosphere;
//    char        *name;
//    } sun;

// For the seeding systems, the initializer that makes the most sense takes an "orbit #", "AU" distance, and "E" eccentricity, axial tilt (degrees), and mass.
// From "SolarSystem.swift" conversion space: /*       No Orbit   Eccen. Tilt Mass    Gas Giant? Dust Mass   Gas */
// No reason we couldn't splat in a name as well ;-)
final class Planet: Equatable {
    // NOTE(heckj): I slapped in a quick unique identifier for each generated planet to easily
    // compare linked list references...
    static func == (lhs: Planet, rhs: Planet) -> Bool {
        lhs.id == rhs.id
    }

    var id: UUID = .init()

    var planet_no: Int
    var a: Double /* semi-major axis of solar orbit (in AU)*/
    var e: Double /* eccentricity of solar orbit         */
    var axial_tilt: Double /* units in degrees */
    var mass: Double /* mass (in solar masses)             */
    var gas_giant: Bool /* TRUE if the planet is a gas giant */
    var dust_mass: Double /* mass, ignoring gas                 */
    var gas_mass: Double /* mass, ignoring dust                 */

    var core_radius: Double = 0 /* radius of the rocky core (in km)     */
    var radius: Double = 0 /* equatorial radius (in km)         */
    var orbit_zone: Int = 0 /* the 'zone' of the planet             */
    var density: Double = 0 /* density (in g/cc)                 */
    var orb_period: Double = 0 /* length of the local year (days)     */
    var day: Double = 0 /* length of the local day (hours)     */
    var resonant_period: Bool = false /* TRUE if in resonant rotation         */
    var esc_velocity: Double = 0 /* units of cm/sec                     */
    var surf_accel: Double = 0 /* units of cm/sec2                     */
    var surf_grav: Double = 0 /* units of Earth gravities             */
    var rms_velocity: Double = 0 /* units of cm/sec                     */
    var molec_weight: Double = 0 /* smallest molecular weight retained*/
    var volatile_gas_inventory: Double = 0
    var surf_pressure: Double = 0 /* units of millibars (mb)             */
    var greenhouse_effect: Bool = false /* runaway greenhouse effect?         */
    var boil_point: Double = 0 /* the boiling point of water (Kelvin)*/
    var albedo: Double = 0 /* albedo of the planet                 */
    var exospheric_temp: Double = 0 /* units of degrees Kelvin             */
    var estimated_temp: Double = 0 /* quick non-iterative estimate (K)  */
    var estimated_terr_temp: Double = 0 /* for terrestrial moons and the like*/
    var surf_temp: Double = 0 /* surface temperature in Kelvin     */
    var greenhs_rise: Double = 0 /* Temperature rise due to greenhouse */
    var high_temp: Double = 0 /* Day-time temperature              */
    var low_temp: Double = 0 /* Night-time temperature             */
    var max_temp: Double = 0 /* Summer/Day                         */
    var min_temp: Double = 0 /* Winter/Night                         */
    var hydrosphere: Double = 0 /* fraction of surface covered         */
    var cloud_cover: Double = 0 /* fraction of surface covered         */
    var ice_cover: Double = 0 /* fraction of surface covered         */

    var sun: Sun?
    var gases: Int = 0 /* Count of gases in the atmosphere: */
    var atmosphere: [Gas] = []

    var planet_type: PlanetType = .unknown
    var minor_moons: Int = 0
    var first_moon: Planet?
    var next_planet: Planet?

    init(planet_no: Int, a: Double, e: Double, axial_tilt: Double, mass: Double, gas_giant: Bool, dust_mass: Double, gas_mass: Double, first_moon: Planet? = nil, next_planet: Planet?) {
        self.planet_no = planet_no
        self.a = a
        self.e = e
        self.axial_tilt = axial_tilt
        self.mass = mass
        self.gas_giant = gas_giant
        self.dust_mass = dust_mass
        self.gas_mass = gas_mass
        self.first_moon = first_moon
        self.next_planet = next_planet
    }
}

// typedef struct planets_record {
//    int            planet_no;
//    long double a;                    /* semi-major axis of solar orbit (in AU)*/
//    long double e;                    /* eccentricity of solar orbit         */
//    long double    axial_tilt;            /* units of degrees                     */
//    long double mass;                /* mass (in solar masses)             */
//    int         gas_giant;            /* TRUE if the planet is a gas giant */
//    long double    dust_mass;            /* mass, ignoring gas                 */
//    long double    gas_mass;            /* mass, ignoring dust                 */
//                                    /*   ZEROES start here               */
//    long double moon_a;                /* semi-major axis of lunar orbit (in AU)*/
//    long double moon_e;                /* eccentricity of lunar orbit         */
//    long double    core_radius;        /* radius of the rocky core (in km)     */
//    long double radius;                /* equatorial radius (in km)         */
//    int         orbit_zone;            /* the 'zone' of the planet             */
//    long double density;            /* density (in g/cc)                 */
//    long double orb_period;            /* length of the local year (days)     */
//    long double day;                /* length of the local day (hours)     */
//    int         resonant_period;    /* TRUE if in resonant rotation         */
//    long double    esc_velocity;        /* units of cm/sec                     */
//    long double    surf_accel;            /* units of cm/sec2                     */
//    long double    surf_grav;            /* units of Earth gravities             */
//    long double    rms_velocity;        /* units of cm/sec                     */
//    long double molec_weight;        /* smallest molecular weight retained*/
//    long double    volatile_gas_inventory;
//    long double    surf_pressure;        /* units of millibars (mb)             */
//    int             greenhouse_effect;    /* runaway greenhouse effect?         */
//    long double    boil_point;            /* the boiling point of water (Kelvin)*/
//    long double    albedo;                /* albedo of the planet                 */
//    long double    exospheric_temp;    /* units of degrees Kelvin             */
//    long double estimated_temp;     /* quick non-iterative estimate (K)  */
//    long double estimated_terr_temp;/* for terrestrial moons and the like*/
//    long double    surf_temp;            /* surface temperature in Kelvin     */
//    long double    greenhs_rise;        /* Temperature rise due to greenhouse */
//    long double high_temp;            /* Day-time temperature              */
//    long double low_temp;            /* Night-time temperature             */
//    long double max_temp;            /* Summer/Day                         */
//    long double min_temp;            /* Winter/Night                         */
//    long double    hydrosphere;        /* fraction of surface covered         */
//    long double    cloud_cover;        /* fraction of surface covered         */
//    long double    ice_cover;            /* fraction of surface covered         */
//    sun*        sun;
//    int            gases;                /* Count of gases in the atmosphere: */
//    gas*        atmosphere;
//    planet_type type;                /* Type code                         */
//    int            minor_moons;
//    planet_pointer first_moon;
//                                    /*   ZEROES end here               */
//    planet_pointer next_planet;
//    } planets;

/*    Define the solar system for comparisons, etc. */
// #define ZEROES 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,NULL,0,NULL,tUnknown

final class Dust {
    var inner_edge: Double
    var outer_edge: Double
    var dust_present: Bool
    var gas_present: Bool
    var next_band: Dust?

    init(inner_edge: Double, outer_edge: Double, dust_present: Bool, gas_present: Bool, next_band: Dust?) {
        self.inner_edge = inner_edge
        self.outer_edge = outer_edge
        self.dust_present = dust_present
        self.gas_present = gas_present
        self.next_band = next_band
    }
}

// typedef struct dust_record {
//    long double inner_edge;
//    long double outer_edge;
//    int         dust_present;
//    int         gas_present;
//    dust_pointer next_band;
//     } dust;

final class Star {
    var luminosity: Double
    var mass: Double
    var m2: Double
    var e: Double
    var a: Double
    var known_planets: Planet?
    var desig: String
    var in_celestia: Bool
    var name: String

    /*    L  Mass    Mass2    Eccen.    SemiMajorAxis    Designation    Name    */
    init(luminosity: Double, mass: Double, m2: Double, e: Double, a: Double, known_planets: Planet?, desig: String, in_celestia: Bool, name: String) {
        self.luminosity = luminosity
        self.mass = mass
        self.m2 = m2
        self.e = e
        self.a = a
        self.known_planets = known_planets
        self.desig = desig
        self.in_celestia = in_celestia
        self.name = name
    }

    //// L            Mass            Mass2            Eccen.    SMAxis     Planets    Designation    Name
    // {{1.00,            1.00,            0,                0,        0,         &mercury,    "Sol",         1, "The Solar System"},        // 0
    init(_ luminosity: Double, _ mass: Double, _ m2: Double, _ e: Double, _ a: Double, _ known_planets: Planet?, _ desig: String, _ in_celestia: Int, _ name: String) {
        self.luminosity = luminosity
        self.mass = mass
        self.m2 = m2
        self.e = e
        self.a = a
        self.known_planets = known_planets
        self.desig = desig
        self.in_celestia = in_celestia > 0
        self.name = name
    }
}

// typedef struct star {
//    long double        luminosity;
//    long double        mass;
//    long double        m2;
//    long double        e;
//    long double        a;
//    planet_pointer    known_planets;
//    char            *desig;
//    int                in_celestia;
//    char            *name;
//    } star;

final class Catalog {
    var count: Int
    var arg: String
    var stars: [Star]

    init(count: Int, arg: String, stars: [Star]) {
        self.count = count
        self.arg = arg
        self.stars = stars
    }
}

// typedef struct catalog {
//    int                count;
//    char*            arg;
//    star            (*stars)[];
//    } catalog;

final class Generation {
    var dust: Dust?
    var planet: Planet?
    var next: Generation?

    init(dust: Dust?, planet: Planet?, next: Generation?) {
        self.dust = dust
        self.planet = planet
        self.next = next
    }
}

// typedef    struct gen
// {
//    dust_pointer    dusts;
//    planet_pointer    planets;
//    gen_pointer        next;
// } generation;
