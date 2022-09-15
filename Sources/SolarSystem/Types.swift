//
//  Types.swift
//  

//#pragma once
//
//typedef struct dust_record    *dust_pointer;
//typedef struct planets_record  *planet_pointer;
//typedef struct gen *gen_pointer;

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

//typedef enum planet_type {
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
//} planet_type;

final class Gas {
    var num: Int
    var surf_pressure: Double /* units of millibars (mb)             */
    
    init(num: Int, surf_pressure: Double) {
        self.num = num
        self.surf_pressure = surf_pressure
    }
}
//typedef struct gas {
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
//typedef struct sun {
//    long double    luminosity;
//    long double    mass;
//    long double life;
//    long double age;
//    long double r_ecosphere;
//    char        *name;
//    } sun;

final class Planet {
    var planet_no: Int
    var a: Double /* semi-major axis of solar orbit (in AU)*/
    var e: Double /* eccentricity of solar orbit         */
    var axial_tilt: Double /* units in degrees */
    var mass: Double /* mass (in solar masses)             */
    var gas_giant: Bool /* TRUE if the planet is a gas giant */
    var dust_mass: Double /* mass, ignoring gas                 */
    var gas_mass: Double /* mass, ignoring dust                 */
    
    var moon_a: Double /* semi-major axis of lunar orbit (in AU)*/
    var moon_e: Double /* eccentricity of lunar orbit         */
    var core_radius: Double /* radius of the rocky core (in km)     */
    var radius: Double /* equatorial radius (in km)         */
    var orbit_zone: Int /* the 'zone' of the planet             */
    var density: Double /* density (in g/cc)                 */
    var orb_period: Double /* length of the local year (days)     */
    var day: Double /* length of the local day (hours)     */
    var resonant_period: Bool /* TRUE if in resonant rotation         */
    var esc_velocity: Double /* units of cm/sec                     */
    var surf_accel: Double /* units of cm/sec2                     */
    var surf_grav: Double /* units of Earth gravities             */
    var    rms_velocity: Double        /* units of cm/sec                     */
    var molec_weight: Double        /* smallest molecular weight retained*/
    var    volatile_gas_inventory: Double
    var    surf_pressure: Double        /* units of millibars (mb)             */
    var             greenhouse_effect: Bool    /* runaway greenhouse effect?         */
    var    boil_point: Double            /* the boiling point of water (Kelvin)*/
    var    albedo: Double                /* albedo of the planet                 */
    var    exospheric_temp: Double    /* units of degrees Kelvin             */
    var estimated_temp: Double     /* quick non-iterative estimate (K)  */
    var estimated_terr_temp: Double /* for terrestrial moons and the like*/
    var    surf_temp: Double            /* surface temperature in Kelvin     */
    var    greenhs_rise: Double        /* Temperature rise due to greenhouse */
    var high_temp: Double            /* Day-time temperature              */
    var low_temp: Double            /* Night-time temperature             */
    var max_temp: Double            /* Summer/Day                         */
    var min_temp: Double            /* Winter/Night                         */
    var    hydrosphere: Double        /* fraction of surface covered         */
    var    cloud_cover: Double        /* fraction of surface covered         */
    var    ice_cover: Double            /* fraction of surface covered         */

    var sun: Sun
    var gases: Int /* Count of gases in the atmosphere: */
    var atmosphere: Gas?
    
    var planet_type: PlanetType
    var minor_moons: Int
    var first_moon: Planet?
    var next_planet: Planet?
    
    init(planet_no: Int, a: Double, e: Double, axial_tilt: Double, mass: Double, gas_giant: Bool, dust_mass: Double, gas_mass: Double, moon_a: Double, moon_e: Double, core_radius: Double, radius: Double, orbit_zone: Int, density: Double, orb_period: Double, day: Double, resonant_period: Bool, esc_velocity: Double, surf_accel: Double, surf_grav: Double, rms_velocity: Double, molec_weight: Double, volatile_gas_inventory: Double, surf_pressure: Double, greenhouse_effect: Bool, boil_point: Double, albedo: Double, exospheric_temp: Double, estimated_temp: Double, estimated_terr_temp: Double, surf_temp: Double, greenhs_rise: Double, high_temp: Double, low_temp: Double, max_temp: Double, min_temp: Double, hydrosphere: Double, cloud_cover: Double, ice_cover: Double, sun: Sun, gases: Int, atmosphere: Gas?, planet_type: PlanetType, minor_moons: Int, first_moon: Planet?, next_planet: Planet?) {
        self.planet_no = planet_no
        self.a = a
        self.e = e
        self.axial_tilt = axial_tilt
        self.mass = mass
        self.gas_giant = gas_giant
        self.dust_mass = dust_mass
        self.gas_mass = gas_mass
        self.moon_a = moon_a
        self.moon_e = moon_e
        self.core_radius = core_radius
        self.radius = radius
        self.orbit_zone = orbit_zone
        self.density = density
        self.orb_period = orb_period
        self.day = day
        self.resonant_period = resonant_period
        self.esc_velocity = esc_velocity
        self.surf_accel = surf_accel
        self.surf_grav = surf_grav
        self.rms_velocity = rms_velocity
        self.molec_weight = molec_weight
        self.volatile_gas_inventory = volatile_gas_inventory
        self.surf_pressure = surf_pressure
        self.greenhouse_effect = greenhouse_effect
        self.boil_point = boil_point
        self.albedo = albedo
        self.exospheric_temp = exospheric_temp
        self.estimated_temp = estimated_temp
        self.estimated_terr_temp = estimated_terr_temp
        self.surf_temp = surf_temp
        self.greenhs_rise = greenhs_rise
        self.high_temp = high_temp
        self.low_temp = low_temp
        self.max_temp = max_temp
        self.min_temp = min_temp
        self.hydrosphere = hydrosphere
        self.cloud_cover = cloud_cover
        self.ice_cover = ice_cover
        self.sun = sun
        self.gases = gases
        self.atmosphere = atmosphere
        self.planet_type = planet_type
        self.minor_moons = minor_moons
        self.first_moon = first_moon
        self.next_planet = next_planet
    }
}

//typedef struct planets_record {
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
//#define ZEROES 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,NULL,0,NULL,tUnknown

final class Dust {
    var inner_edge: Double
    var outer_edge: Double
    var dust_present: Bool
    var gas_present: Bool
    var next_band: Dust? = nil
    
    init(inner_edge: Double, outer_edge: Double, dust_present: Bool, gas_present: Bool, next_band: Dust?) {
        self.inner_edge = inner_edge
        self.outer_edge = outer_edge
        self.dust_present = dust_present
        self.gas_present = gas_present
        self.next_band = next_band
    }
}

//typedef struct dust_record {
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
    var known_planets: Planet
    var desig: String
    var in_celestia: Bool
    var name: String
    
    init(luminosity: Double, mass: Double, m2: Double, e: Double, a: Double, known_planets: Planet, desig: String, in_celestia: Bool, name: String) {
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
}
//typedef struct star {
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
//typedef struct catalog {
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

//typedef    struct gen
//{
//    dust_pointer    dusts;
//    planet_pointer    planets;
//    gen_pointer        next;
//} generation;

// From Keris

final class ChemTable {
    var num: Int
    var symbol: String
    var html_symbol: String
    var name: String
    var weight: Double
    var melt: Double
    var boil: Double
    var density: Double
    var abunde: Double
    var abunds: Double
    var reactivity: Double
    var max_ipp: Double // Max inspired partial pressure im millibars
    
    init(num: Int, symbol: String, html_symbol: String, name: String, weight: Double, melt: Double, boil: Double, density: Double, abunde: Double, abunds: Double, reactivity: Double, max_ipp: Double) {
        self.num = num
        self.symbol = symbol
        self.html_symbol = html_symbol
        self.name = name
        self.weight = weight
        self.melt = melt
        self.boil = boil
        self.density = density
        self.abunde = abunde
        self.abunds = abunds
        self.reactivity = reactivity
        self.max_ipp = max_ipp
    }
}
//typedef struct ChemTableS
//{
//  int             num;
//  char           *symbol;
//  char            *html_symbol;
//  char           *name;
//  long double    weight;
//  long double    melt;
//  long double    boil;
//  long double    density;
//  long double    abunde;
//  long double    abunds;
//  long double    reactivity;
//  long double    max_ipp;    // Max inspired partial pressure im millibars
//} ChemTable;

