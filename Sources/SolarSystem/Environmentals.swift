//
//  Environmentals.swift
//

import Foundation

func pow3(_ x: Double) -> Double {
    x * x * x
}

func pow2(_ x: Double) -> Double {
    x * x
}

func pow1_4(_ x: Double) -> Double {
    sqrt(sqrt(x))
}

enum Breathability: Int {
    case none = 0
    case breathable = 1
    case unbreathable = 2
    case poisonous = 3
}

extension Breathability: CustomStringConvertible {
    var description: String {
        switch self {
        case .none:
            return "none"
        case .breathable:
            return "breathable"
        case .unbreathable:
            return "unbreathable"
        case .poisonous:
            return "poisonous"
        }
    }
}

public func orbital_range_string(a: Double, e: Double) -> String {
    "\((a * (1 - e)).formatted(AUFormat))AU-\((a * (1 + e)).formatted(AUFormat))AU"
}

/// Determines an approximate luminosity for a star based on stellar mass
/// - Parameter mass_ratio: The mass in solar masses.
/// - Returns: The luminosity relative to solar luminosity.
public func luminosity(mass_ratio: Double) -> Double {
    let n: Double
    if mass_ratio < 1.0 {
        n = 1.75 * (mass_ratio - 0.1) + 3.325 // 3.325 (0.1) -> ~4.9 (->1.0)
    } else if mass_ratio < 2.0 {
        n = 0.5 * (2.0 - mass_ratio) + 4.4 // 4.0 (1) -> 4.4 (2)
    } else if mass_ratio < 16.0 {
        n = 0.7 * (16.0 - mass_ratio) + 3.7 // 4.4 (2) -> 3.7 (16)
    } else {
        n = 3.5
    }
    return (pow(mass_ratio, n))
}

/*--------------------------------------------------------------------------*/
/*     This function, given the orbital radius of a planet in AU, returns        */
/*     the orbital 'zone' of the particle.                                    */
/*--------------------------------------------------------------------------*/
func orbital_zone(luminosity: Double, orbital_radius: Double) -> Int {
    if orbital_radius < (4.0 * sqrt(luminosity)) {
        return (1)
    } else if orbital_radius < (15.0 * sqrt(luminosity)) {
        return (2)
    } else {
        return (3)
    }
}

/// *--------------------------------------------------------------------------*/
/// *     The mass is in units of solar masses, and the density is in units        */
/// *     of grams/cc.  The radius returned is in units of km.                    */
/// *--------------------------------------------------------------------------*/
func volume_radius(mass: Double, density: Double) -> Double {
    let mass_g = mass * SOLAR_MASS_IN_GRAMS
    let volume = mass_g / density
    return (pow((3.0 * volume) / (4.0 * Double.pi), 1.0 / 3.0) / CM_PER_KM)
}

/// *--------------------------------------------------------------------------*/
/// *     Returns the radius of the planet in kilometers.                        */
/// *     The mass passed in is in units of solar masses.                        */
/// *     This formula is listed as eq.9 in Fogg's article, although some typos    */
/// *     crop up in that eq.  See "The Internal Constitution of Planets", by    */
/// *     Dr. D. S. Kothari, Mon. Not. of the Royal Astronomical Society, vol 96 */
/// *     pp.833-843, 1936 for the derivation.  Specifically, this is Kothari's    */
/// *     eq.23, which appears on page 840.                                        */
/// *--------------------------------------------------------------------------*/
func kothari_radius(mass: Double, giant: Bool, zone: Int) -> Double {
    let atomic_weight: Double
    let atomic_num: Double

    if zone == 1 {
        if giant {
            atomic_weight = 9.5
            atomic_num = 4.5
        } else {
            atomic_weight = 15.0
            atomic_num = 8.0
        }
    } else if zone == 2 {
        if giant {
            atomic_weight = 2.47
            atomic_num = 2.0
        } else {
            atomic_weight = 10.0
            atomic_num = 5.0
        }
    } else { // zone == 3
        if giant {
            atomic_weight = 7.0
            atomic_num = 4.0
        } else {
            atomic_weight = 10.0
            atomic_num = 5.0
        }
    }

    let temp1 = atomic_weight * atomic_num
    var temp = (2.0 * BETA_20 * pow(SOLAR_MASS_IN_GRAMS, 1.0 / 3.0))
        / (A1_20 * pow(temp1, 1.0 / 3.0))

    var temp2 = A2_20 * pow(atomic_weight, 4.0 / 3.0) * pow(SOLAR_MASS_IN_GRAMS, 2.0 / 3.0)
    temp2 = temp2 * pow(mass, 2.0 / 3.0)
    temp2 = temp2 / (A1_20 * (atomic_num * atomic_num))
    temp2 = 1.0 + temp2
    temp = temp / temp2
    temp = (temp * pow(mass, 1.0 / 3.0)) / CM_PER_KM

    temp /= JIMS_FUDGE /* Make Earth = actual earth */
    return (temp)
}

/// *--------------------------------------------------------------------------*/
/// *    The mass passed in is in units of solar masses, and the orbital radius    */
/// *    is in units of AU.    The density is returned in units of grams/cc.        */
/// *--------------------------------------------------------------------------*/
func empirical_density(mass: Double, orb_radius: Double,
                       r_ecosphere: Double, gas_giant: Bool) -> Double
{
    var temp = pow(mass * SUN_MASS_IN_EARTH_MASSES, 1.0 / 8.0)
    temp = temp * pow1_4(r_ecosphere / orb_radius)
    if gas_giant {
        return (temp * 1.2)
    }
    return (temp * 5.5)
}

/// *--------------------------------------------------------------------------*/
/// *    The mass passed in is in units of solar masses, and the equatorial        */
/// *    radius is in km.  The density is returned in units of grams/cc.            */
/// *--------------------------------------------------------------------------*/
func volume_density(mass: Double, equat_radius: Double) -> Double {
    let mass_grams = mass * SOLAR_MASS_IN_GRAMS
    let equat_radius_in_cm = equat_radius * CM_PER_KM
    let volume = (4.0 * Double.pi * pow3(equat_radius_in_cm)) / 3.0
    return (mass_grams / volume)
}

/// *--------------------------------------------------------------------------*/
/// *    The separation is in units of AU, and both masses are in units of solar */
/// *    masses.     The period returned is in terms of Earth days.                    */
/// *--------------------------------------------------------------------------*/
func period(separation: Double, small_mass: Double, large_mass: Double) -> Double {
    let period_in_years = sqrt(pow3(separation) / (small_mass + large_mass))
    return (period_in_years * DAYS_IN_A_YEAR)
}

/// *--------------------------------------------------------------------------*/
/// *     Fogg's information for this routine came from Dole "Habitable Planets    */
/// * for Man", Blaisdell Publishing Company, NY, 1964.  From this, he came    */
/// * up with his eq.12, which is the equation for the 'base_angular_velocity' */
/// * below.  He then used an equation for the change in angular velocity per    */
/// * time (dw/dt) from P. Goldreich and S. Soter's paper "Q in the Solar        */
/// * System" in Icarus, vol 5, pp.375-389 (1966).     Using as a comparison the    */
/// * change in angular velocity for the Earth, Fogg has come up with an        */
/// * approximation for our new planet (his eq.13) and take that into account. */
/// * This is used to find 'change_in_angular_velocity' below.                    */
/// *                                                                            */
/// *     Input parameters are mass (in solar masses), radius (in Km), orbital    */
/// * period (in days), orbital radius (in AU), density (in g/cc),                */
/// * eccentricity, and whether it is a gas giant or not.                        */
/// *     The length of the day is returned in units of hours.                    */
/// *--------------------------------------------------------------------------*/
//
func day_length(planet: Planet) -> Double {
    let planetary_mass_in_grams = planet.mass * SOLAR_MASS_IN_GRAMS
    let equatorial_radius_in_cm = planet.radius * CM_PER_KM
    let year_in_hours = planet.orb_period * 24.0
    let giant: Bool = planet.planet_type == .subgasgiant || planet.planet_type == .gasgiant || planet.planet_type == .subsubgasgiant

    let k2: Double
    let day_in_hours: Double

    var stopped = false
    planet.resonant_period = false /* Warning: Modify the planet */

    if giant {
        k2 = 0.24
    } else {
        k2 = 0.33
    }

    let base_angular_velocity = sqrt(2.0 * J * planetary_mass_in_grams /
        (k2 * equatorial_radius_in_cm * equatorial_radius_in_cm))

    /*    This next calculation determines how much the planet's rotation is     */
    /*    slowed by the presence of the star.                                 */

    let change_in_angular_velocity = CHANGE_IN_EARTH_ANG_VEL *
        (planet.density / EARTH_DENSITY) *
        (equatorial_radius_in_cm / EARTH_RADIUS) *
        (EARTH_MASS_IN_GRAMS / planetary_mass_in_grams) *
        pow(planet.sun!.mass, 2.0) *
        (1.0 / pow(planet.a, 6.0))

    let ang_velocity = base_angular_velocity + (change_in_angular_velocity * planet.sun!.age)

    /* Now we change from rad/sec to hours/rotation. */

    if ang_velocity <= 0.0 {
        stopped = true
        day_in_hours = INCREDIBLY_LARGE_NUMBER
    } else {
        day_in_hours = RADIANS_PER_ROTATION / (SECONDS_PER_HOUR * ang_velocity)
    }

    if (day_in_hours >= year_in_hours) || stopped {
        if planet.e > 0.1 {
            let spin_resonance_factor = (1.0 - planet.e) / (1.0 + planet.e)
            planet.resonant_period = true
            return (spin_resonance_factor * year_in_hours)
        } else {
            return (year_in_hours)
        }
    }
    return (day_in_hours)
}

/// *--------------------------------------------------------------------------*/
/// *     The orbital radius is expected in units of Astronomical Units (AU).    */
/// *     Inclination is returned in units of degrees.                            */
/// *--------------------------------------------------------------------------*/
func inclination(orb_radius: Double, prng: RNGWrapper<Xoshiro>) -> Double {
    pow(orb_radius, 0.2) * prng.about(EARTH_AXIAL_TILT, variation: 0.4) / 360.0
}

/// *--------------------------------------------------------------------------*/
/// *     This function implements the escape velocity calculation.    Note that    */
/// *    it appears that Fogg's eq.15 is incorrect.                                */
/// *    The mass is in units of solar mass, the radius in kilometers, and the    */
/// *    velocity returned is in cm/sec.                                            */
/// *--------------------------------------------------------------------------*/
func escape_vel(mass: Double, radius: Double) -> Double {
    let mass_in_grams = mass * SOLAR_MASS_IN_GRAMS
    let radius_in_cm = radius * CM_PER_KM
    return (sqrt(2.0 * GRAV_CONSTANT * mass_in_grams / radius_in_cm))
}

/// *--------------------------------------------------------------------------*/
/// *    This is Fogg's eq.16.  The molecular weight (usually assumed to be N2)    */
/// *    is used as the basis of the Root Mean Square (RMS) velocity of the        */
/// *    molecule or atom.  The velocity returned is in cm/sec.                    */
/// *    Orbital radius is in A.U.(ie: in units of the earth's orbital radius).    */
/// *--------------------------------------------------------------------------*/
//
func rms_vel(molecular_weight: Double, exospheric_temp: Double) -> Double {
    sqrt((3.0 * MOLAR_GAS_CONST * exospheric_temp) / molecular_weight)
        * CM_PER_METER
}

/// *--------------------------------------------------------------------------*/
/// *     This function returns the smallest molecular weight retained by the    */
/// *    body, which is useful for determining the atmosphere composition.        */
/// *    Mass is in units of solar masses, and equatorial radius is in units of    */
/// *    kilometers.                                                                */
/// *--------------------------------------------------------------------------*/
//
func molecule_limit(mass: Double, equat_radius: Double, exospheric_temp: Double) -> Double {
    let esc_velocity = escape_vel(mass: mass, radius: equat_radius)

    return (3.0 * MOLAR_GAS_CONST * exospheric_temp) / pow2((esc_velocity / GAS_RETENTION_THRESHOLD) / CM_PER_METER)
}

/// *--------------------------------------------------------------------------*/
/// *     This function calculates the surface acceleration of a planet.     The    */
/// *    mass is in units of solar masses, the radius in terms of km, and the    */
/// *    acceleration is returned in units of cm/sec2.                            */
/// *--------------------------------------------------------------------------*/
//
func acceleration(mass: Double, radius: Double) -> Double {
    GRAV_CONSTANT * (mass * SOLAR_MASS_IN_GRAMS) /
        pow2(radius * CM_PER_KM)
}

/// *--------------------------------------------------------------------------*/
/// *     This function calculates the surface gravity of a planet.    The            */
/// *    acceleration is in units of cm/sec2, and the gravity is returned in        */
/// *    units of Earth gravities.                                                */
/// *--------------------------------------------------------------------------*/
//
func gravity(acceleration: Double) -> Double {
    acceleration / EARTH_ACCELERATION
}

/// *--------------------------------------------------------------------------*/
/// *    This implements Fogg's eq.17.  The 'inventory' returned is unitless.    */
/// *--------------------------------------------------------------------------*/

func vol_inventory(mass: Double, escape_vel: Double, rms_vel: Double, stellar_mass: Double, zone: Int,
                   greenhouse_effect: Bool, accreted_gas: Bool, prng: RNGWrapper<Xoshiro>) -> Double
{
    // proportion_const, temp1, temp2, earth_units;
    //
    let velocity_ratio = escape_vel / rms_vel
    let proportion_const: Double
    if velocity_ratio >= GAS_RETENTION_THRESHOLD {
        switch zone {
        case 1:
            proportion_const = 140_000.0 /* 100 -> 140 JLB */
        case 2:
            proportion_const = 75000.0
        case 3:
            proportion_const = 250.0
        default:
            proportion_const = 0.0
            print("Error: orbital zone not initialized correctly!\n")
        }

        let earth_units = mass * SUN_MASS_IN_EARTH_MASSES
        let temp1 = (proportion_const * earth_units) / stellar_mass
        let temp2 = prng.about(temp1, variation: 0.2)

        if greenhouse_effect || accreted_gas {
            return (temp2)
        } else {
            return (temp2 / 140.0) /* 100 -> 140 JLB */
        }
    }
    return (0.0)
}

/// *--------------------------------------------------------------------------*/
/// *    This implements Fogg's eq.18.  The pressure returned is in units of        */
/// *    millibars (mb).     The gravity is in units of Earth gravities, the radius */
/// *    in units of kilometers.                                                    */
/// *                                                                            */
/// *  JLB: Apparently this assumed that earth pressure = 1000mb. I've added a    */
/// *    fudge factor (EARTH_SURF_PRES_IN_MILLIBARS / 1000.) to correct for that    */
/// *--------------------------------------------------------------------------*/
//
func pressure(volatile_gas_inventory: Double, equat_radius: Double, gravity: Double) -> Double {
    let equat_radius_in_earths = KM_EARTH_RADIUS / equat_radius
    return (volatile_gas_inventory * gravity * (EARTH_SURF_PRES_IN_MILLIBARS / 1000.0) / pow2(equat_radius_in_earths))
}

/// *--------------------------------------------------------------------------*/
/// *     This function returns the boiling point of water in an atmosphere of    */
/// *     pressure 'surf_pressure', given in millibars.    The boiling point is    */
/// *     returned in units of Kelvin.  This is Fogg's eq.21.                    */
/// *--------------------------------------------------------------------------*/
func boiling_point(surf_pressure: Double) -> Double {
    let surface_pressure_in_bars = surf_pressure / MILLIBARS_PER_BAR
    return (1.0 / ((log(surface_pressure_in_bars) / -5050.5) +
            (1.0 / 373.0)))
}

/// *--------------------------------------------------------------------------*/
/// *     This function is Fogg's eq.22.     Given the volatile gas inventory and    */
/// *     planetary radius of a planet (in Km), this function returns the        */
/// *     fraction of the planet covered with water.                                */
/// *     I have changed the function very slightly:     the fraction of Earth's    */
/// *     surface covered by water is 71%, not 75% as Fogg used.                    */
/// *--------------------------------------------------------------------------*/
func hydro_fraction(volatile_gas_inventory: Double, planet_radius: Double) -> Double {
    let temp = (0.71 * volatile_gas_inventory / 1000.0)
        * pow2(KM_EARTH_RADIUS / planet_radius)
    if temp >= 1.0 {
        return (1.0)
    }
    return (temp)
}

/// *--------------------------------------------------------------------------*/
/// *     Given the surface temperature of a planet (in Kelvin), this function    */
/// *     returns the fraction of cloud cover available.     This is Fogg's eq.23.    */
/// *     See Hart in "Icarus" (vol 33, pp23 - 39, 1978) for an explanation.        */
/// *     This equation is Hart's eq.3.                                            */
/// *     I have modified it slightly using constants and relationships from        */
/// *     Glass's book "Introduction to Planetary Geology", p.46.                */
/// *     The 'CLOUD_COVERAGE_FACTOR' is the amount of surface area on Earth        */
/// *     covered by one Kg. of cloud.                                            */
/// *--------------------------------------------------------------------------*/
//
func cloud_fraction(surf_temp: Double, smallest_MW_retained: Double, equat_radius: Double, hydro_fraction: Double) -> Double {
    if smallest_MW_retained > WATER_VAPOR {
        return 0.0
    }
    let surf_area = 4.0 * Double.pi * pow2(equat_radius)
    let hydro_mass = hydro_fraction * surf_area * EARTH_WATER_MASS_PER_AREA
    let water_vapor_in_kg = (0.00000001 * hydro_mass) * exp(Q2_36 * (surf_temp - EARTH_AVERAGE_KELVIN))
    let fraction = CLOUD_COVERAGE_FACTOR * water_vapor_in_kg / surf_area
    if fraction >= 1.0 {
        return 1.0
    }
    return fraction
}

/// *--------------------------------------------------------------------------*/
/// *     Given the surface temperature of a planet (in Kelvin), this function    */
/// *     returns the fraction of the planet's surface covered by ice.  This is    */
/// *     Fogg's eq.24.    See Hart[24] in Icarus vol.33, p.28 for an explanation. */
/// *     I have changed a constant from 70 to 90 in order to bring it more in    */
/// *     line with the fraction of the Earth's surface covered with ice, which    */
/// *     is approximatly .016 (=1.6%).                                            */
/// *--------------------------------------------------------------------------*/
//
func ice_fraction(hydro_fraction: Double, surf_temp: Double) -> Double {
    let clamped_surface_temp = (surf_temp > 328.0) ? 328.0 : surf_temp

    var temp = pow((328.0 - clamped_surface_temp) / 90.0, 5.0)
    if temp > (1.5 * hydro_fraction) {
        temp = (1.5 * hydro_fraction)
    }
    if temp >= 1.0 {
        return (1.0)
    }

    return (temp)
}

/// *--------------------------------------------------------------------------*/
/// *    This is Fogg's eq.19.  The ecosphere radius is given in AU, the orbital */
/// *    radius in AU, and the temperature returned is in Kelvin.                */
/// *--------------------------------------------------------------------------*/
public func est_temp(ecosphere_radius: Double, orb_radius: Double, albedo: Double) -> Double {
    sqrt(ecosphere_radius / orb_radius)
        * pow1_4((1.0 - albedo) / (1.0 - EARTH_ALBEDO))
        * EARTH_AVERAGE_KELVIN
}

/// *--------------------------------------------------------------------------*/
/// * Old grnhouse:                                                            */
/// *    Note that if the orbital radius of the planet is greater than or equal    */
/// *    to R_inner, 99% of it's volatiles are assumed to have been deposited in */
/// *    surface reservoirs (otherwise, it suffers from the greenhouse effect).    */
/// *--------------------------------------------------------------------------*/
/// *    if ((orb_radius < r_greenhouse) && (zone == 1)) */
//
/// *--------------------------------------------------------------------------*/
/// *    The new definition is based on the inital surface temperature and what    */
/// *    state water is in. If it's too hot, the water will never condense out    */
/// *    of the atmosphere, rain down and form an ocean. The albedo used here    */
/// *    was chosen so that the boundary is about the same as the old method        */
/// *    Neither zone, nor r_greenhouse are used in this version                JLB    */
/// *--------------------------------------------------------------------------*/
func grnhouse(r_ecosphere: Double, orb_radius: Double) -> Bool {
    let temp = est_temp(ecosphere_radius: r_ecosphere, orb_radius: orb_radius, albedo: GREENHOUSE_TRIGGER_ALBEDO)
    return temp > FREEZING_POINT_OF_WATER
}

/// *--------------------------------------------------------------------------*/
/// *    This is Fogg's eq.20, and is also Hart's eq.20 in his "Evolution of        */
/// *    Earth's Atmosphere" article.  The effective temperature given is in        */
/// *    units of Kelvin, as is the rise in temperature produced by the            */
/// *    greenhouse effect, which is returned.                                    */
/// *    I tuned this by changing a pow(x,.25) to pow(x,.4) to match Venus - JLB    */
/// *--------------------------------------------------------------------------*/
func green_rise(optical_depth: Double, effective_temp: Double, surf_pressure: Double) -> Double {
    let convection_factor = EARTH_CONVECTION_FACTOR * pow(surf_pressure / EARTH_SURF_PRES_IN_MILLIBARS, 0.4)
    var rise = (pow1_4(1.0 + 0.75 * optical_depth) - 1.0) * effective_temp * convection_factor

    if rise < 0.0 {
        rise = 0.0
    }
    return rise
}

/// *--------------------------------------------------------------------------*/
/// *     The surface temperature passed in is in units of Kelvin.                */
/// *     The cloud adjustment is the fraction of cloud cover obscuring each        */
/// *     of the three major components of albedo that lie below the clouds.        */
/// *--------------------------------------------------------------------------*/
//
func planet_albedo(water_fraction: Double, cloud_fraction: Double, ice_fraction: Double, surf_pressure: Double) -> Double {
    let rock_fraction = 1.0 - water_fraction - ice_fraction
    var components = 0.0
    if water_fraction > 0.0 {
        components = components + 1.0
    }
    if ice_fraction > 0.0 {
        components = components + 1.0
    }
    if rock_fraction > 0.0 {
        components = components + 1.0
    }

    let cloud_adjustment = cloud_fraction / components
    let adj_rock_fraction: Double
    let adj_water_fraction: Double
    let adj_ice_fraction: Double

    if rock_fraction >= cloud_adjustment {
        adj_rock_fraction = rock_fraction - cloud_adjustment
    } else {
        adj_rock_fraction = 0.0
    }

    if water_fraction > cloud_adjustment {
        adj_water_fraction = water_fraction - cloud_adjustment
    } else {
        adj_water_fraction = 0.0
    }

    if ice_fraction > cloud_adjustment {
        adj_ice_fraction = ice_fraction - cloud_adjustment
    } else {
        adj_ice_fraction = 0.0
    }

    let cloud_part = cloud_fraction * CLOUD_ALBEDO /* about(...,0.2); */
    let rock_part: Double
    let water_part: Double
    let ice_part: Double

    if surf_pressure == 0.0 {
        rock_part = adj_rock_fraction * ROCKY_AIRLESS_ALBEDO /* about(...,0.3); */
        ice_part = adj_ice_fraction * AIRLESS_ICE_ALBEDO /* about(...,0.4); */
        water_part = 0
    } else {
        rock_part = adj_rock_fraction * ROCKY_ALBEDO /* about(...,0.1); */
        water_part = adj_water_fraction * WATER_ALBEDO /* about(...,0.2); */
        ice_part = adj_ice_fraction * ICE_ALBEDO /* about(...,0.1); */
    }

    return (cloud_part + rock_part + water_part + ice_part)
}

/// *--------------------------------------------------------------------------*/
/// *     This function returns the dimensionless quantity of optical depth,        */
/// *     which is useful in determining the amount of greenhouse effect on a    */
/// *     planet.                                                                */
/// *--------------------------------------------------------------------------*/
func opacity(molecular_weight: Double, surf_pressure: Double) -> Double {
    var optical_depth = 0.0
    if molecular_weight >= 0.0, molecular_weight < 10.0 {
        optical_depth = optical_depth + 3.0
    }
    if molecular_weight >= 10.0, molecular_weight < 20.0 {
        optical_depth = optical_depth + 2.34
    }
    if molecular_weight >= 20.0, molecular_weight < 30.0 {
        optical_depth = optical_depth + 1.0
    }
    if molecular_weight >= 30.0, molecular_weight < 45.0 {
        optical_depth = optical_depth + 0.15
    }
    if molecular_weight >= 45.0, molecular_weight < 100.0 {
        optical_depth = optical_depth + 0.05
    }

    if surf_pressure >= (70.0 * EARTH_SURF_PRES_IN_MILLIBARS) {
        optical_depth = optical_depth * 8.333
    } else if surf_pressure >= (50.0 * EARTH_SURF_PRES_IN_MILLIBARS) {
        optical_depth = optical_depth * 6.666
    } else if surf_pressure >= (30.0 * EARTH_SURF_PRES_IN_MILLIBARS) {
        optical_depth = optical_depth * 3.333
    } else if surf_pressure >= (10.0 * EARTH_SURF_PRES_IN_MILLIBARS) {
        optical_depth = optical_depth * 2.0
    } else if surf_pressure >= (5.0 * EARTH_SURF_PRES_IN_MILLIBARS) {
        optical_depth = optical_depth * 1.5
    }

    return (optical_depth)
}

/*
 *    calculates the number of years it takes for 1/e of a gas to escape
 *    from a planet's atmosphere.
 *    Taken from Dole p. 34. He cites Jeans (1916) & Jones (1923)
 */
func gas_life(molecular_weight: Double, planet: Planet) -> Double {
    let v = rms_vel(molecular_weight: molecular_weight, exospheric_temp: planet.exospheric_temp)
    let g = planet.surf_grav * EARTH_ACCELERATION
    let r = (planet.radius * CM_PER_KM)
    let t = (pow3(v) / (2.0 * pow2(g) * r)) * exp((3.0 * g * r) / pow2(v))
    var years = t / (SECONDS_PER_HOUR * 24.0 * DAYS_IN_A_YEAR)

    // earlier algorithm:
//    long double ve = planet->esc_velocity;
//    long double k = 2;
//    long double t2 = ((k * pow3(v) * r) / pow4(ve)) * exp((3.0 * pow2(ve)) / (2.0 * pow2(v)));
//    long double years2 = t2 / (SECONDS_PER_HOUR * 24.0 * DAYS_IN_A_YEAR);
//    if (flag_verbose & 0x0040)
//        fprintf (stderr, "gas_life: %LGs, V ratio: %Lf\n",
//                years, ve / v);

    if years > 2.0e10 {
        years = INCREDIBLY_LARGE_NUMBER
    }

    return years
}

func min_molec_weight(planet: Planet) -> Double {
    let mass = planet.mass
    let radius = planet.radius
    let temp = planet.exospheric_temp
    var target = 5.0e9

    var guess_1 = molecule_limit(mass: mass, equat_radius: radius, exospheric_temp: temp)
    var guess_2 = guess_1

    var life = gas_life(molecular_weight: guess_1, planet: planet)
    if let age = planet.sun?.age {
        target = age
    }
    var loops = 0
    if life > target {
        while life > target, loops < 25 {
            guess_1 = guess_1 / 2.0
            life = gas_life(molecular_weight: guess_1, planet: planet)
            loops += 1
        }
    } else {
        while life < target, loops < 25 {
            guess_2 = guess_2 * 2.0
            life = gas_life(molecular_weight: guess_2, planet: planet)
        }
    }

    loops = 0

    while (guess_2 - guess_1) > 0.1, loops < 25 {
        let guess_3 = (guess_1 + guess_2) / 2.0
        life = gas_life(molecular_weight: guess_3, planet: planet)

        if life < target {
            guess_1 = guess_3
        } else {
            guess_2 = guess_3
        }
        loops += 1
    }

    life = gas_life(molecular_weight: guess_2, planet: planet)

    return (guess_2)
}

/// *--------------------------------------------------------------------------*/
/// *     The temperature calculated is in degrees Kelvin.                        */
/// *     Quantities already known which are used in these calculations:            */
/// *         planet->molec_weight                                                */
/// *         planet->surf_pressure                                                */
/// *         R_ecosphere                                                        */
/// *         planet->a                                                            */
/// *         planet->volatile_gas_inventory                                        */
/// *         planet->radius                                                        */
/// *         planet->boil_point                                                    */
/// *--------------------------------------------------------------------------*/
func calculate_surface_temp(planet: Planet,
                            first: Bool,
                            last_water: Double,
                            last_clouds: Double,
                            last_ice: Double,
                            last_temp: Double,
                            last_albedo: Double,
                            prng: RNGWrapper<Xoshiro>)
{
    var effective_temp: Double
//    var water_raw: Double
//    var clouds_raw: Double
    var greenhouse_temp: Double
    var boil_off = false

    if first {
        planet.albedo = EARTH_ALBEDO

        effective_temp = est_temp(ecosphere_radius: planet.sun!.r_ecosphere, orb_radius: planet.a, albedo: planet.albedo)
        let planet_opacity = opacity(molecular_weight: planet.molec_weight, surf_pressure: planet.surf_pressure)
        greenhouse_temp = green_rise(optical_depth: planet_opacity,
                                     effective_temp: effective_temp,
                                     surf_pressure: planet.surf_pressure)
        planet.surf_temp = effective_temp + greenhouse_temp
        set_temp_range(planet: planet)
    }

    if planet.greenhouse_effect, planet.max_temp < planet.boil_point {
//        if (flag_verbose & 0x0010)
//            fprintf (stderr, "Deluge: %s %d max (%Lf) < boil (%Lf)\n",
//                    planet->sun->name,
//                    planet->planet_no,
//                    planet->max_temp,
//                    planet->boil_point);

        planet.greenhouse_effect = false
        planet.volatile_gas_inventory = vol_inventory(mass: planet.mass,
                                                      escape_vel: planet.esc_velocity,
                                                      rms_vel: planet.rms_velocity,
                                                      stellar_mass: planet.sun!.mass,
                                                      zone: planet.orbit_zone,
                                                      greenhouse_effect: planet.greenhouse_effect,
                                                      accreted_gas: (planet.gas_mass / planet.mass) > 0.000001,
                                                      prng: prng)
        planet.surf_pressure = pressure(volatile_gas_inventory: planet.volatile_gas_inventory,
                                        equat_radius: planet.radius,
                                        gravity: planet.surf_grav)
        planet.boil_point = boiling_point(surf_pressure: planet.surf_pressure)
    }

//    water_raw                 =
    planet.hydrosphere = hydro_fraction(volatile_gas_inventory: planet.volatile_gas_inventory,
                                        planet_radius: planet.radius)

//    clouds_raw                 =
    planet.cloud_cover = cloud_fraction(surf_temp: planet.surf_temp,
                                        smallest_MW_retained: planet.molec_weight,
                                        equat_radius: planet.radius,
                                        hydro_fraction: planet.hydrosphere)
    planet.ice_cover = ice_fraction(hydro_fraction: planet.hydrosphere,
                                    surf_temp: planet.surf_temp)

    if planet.greenhouse_effect, planet.surf_pressure > 0.0 {
        planet.cloud_cover = 1.0
    }

    if planet.high_temp >= planet.boil_point, !first, !(Int(planet.day) == Int(planet.orb_period * 24.0) ||
        (planet.resonant_period))
    {
        planet.hydrosphere = 0.0
        boil_off = true

        if planet.molec_weight > WATER_VAPOR {
            planet.cloud_cover = 0.0
        } else {
            planet.cloud_cover = 1.0
        }
    }

    if planet.surf_temp < (FREEZING_POINT_OF_WATER - 3.0) {
        planet.hydrosphere = 0.0
    }

    planet.albedo = planet_albedo(water_fraction: planet.hydrosphere,
                                  cloud_fraction: planet.cloud_cover,
                                  ice_fraction: planet.ice_cover,
                                  surf_pressure: planet.surf_pressure)

    effective_temp = est_temp(ecosphere_radius: planet.sun!.r_ecosphere, orb_radius: planet.a, albedo: planet.albedo)
    let planet_opacity = opacity(molecular_weight: planet.molec_weight, surf_pressure: planet.surf_pressure)
    greenhouse_temp = green_rise(optical_depth: planet_opacity,
                                 effective_temp: effective_temp,
                                 surf_pressure: planet.surf_pressure)
    planet.surf_temp = effective_temp + greenhouse_temp

    if !first {
        if !boil_off {
            planet.hydrosphere = (planet.hydrosphere + (last_water * 2)) / 3
        }
        planet.cloud_cover = (planet.cloud_cover + (last_clouds * 2)) / 3
        planet.ice_cover = (planet.ice_cover + (last_ice * 2)) / 3
        planet.albedo = (planet.albedo + (last_albedo * 2)) / 3
        planet.surf_temp = (planet.surf_temp + (last_temp * 2)) / 3
    }

    set_temp_range(planet: planet)

//    if (flag_verbose & 0x0020)
//        fprintf (stderr, "%5.1Lf AU: %5.1Lf = %5.1Lf ef + %5.1Lf gh%c "
//                "(W: %4.2Lf (%4.2Lf) C: %4.2Lf (%4.2Lf) I: %4.2Lf A: (%4.2Lf))\n",
//                planet->a,
//                planet->surf_temp - FREEZING_POINT_OF_WATER,
//                effective_temp - FREEZING_POINT_OF_WATER,
//                greenhouse_temp,
//                (planet->greenhouse_effect) ? '*' :' ',
//                planet->hydrosphere, water_raw,
//                planet->cloud_cover, clouds_raw,
//                planet->ice_cover,
//                planet->albedo);
}

func iterate_surface_temp(planet: Planet, prng: RNGWrapper<Xoshiro>) {
    let initial_temp = est_temp(ecosphere_radius: planet.sun!.r_ecosphere, orb_radius: planet.a, albedo: planet.albedo)
//    let h2_life  = gas_life (molecular_weight: MOL_HYDROGEN,    planet: planet)
//    let  h2o_life = gas_life (molecular_weight: WATER_VAPOR,     planet: planet)
//    let  n2_life  = gas_life (molecular_weight: MOL_NITROGEN,    planet: planet)
//    let  n_life   = gas_life (molecular_weight: ATOMIC_NITROGEN, planet: planet)

//    if (flag_verbose & 0x20000)
//        fprintf (stderr, "%d:                     %5.1Lf it [%5.1Lf re %5.1Lf a %5.1Lf alb]\n",
//                planet->planet_no,
//                initial_temp,
//                planet->sun->r_ecosphere, planet->a, planet->albedo
//                );

//    if (flag_verbose & 0x0040)
//        fprintf (stderr, "\nGas lifetimes: H2 - %Lf, H2O - %Lf, N - %Lf, N2 - %Lf\n",
//                h2_life, h2o_life, n_life, n2_life);

    calculate_surface_temp(planet: planet, first: true, last_water: 0, last_clouds: 0, last_ice: 0, last_temp: 0, last_albedo: 0, prng: prng)

    for _ in 0 ... 25 {
        let last_water = planet.hydrosphere
        let last_clouds = planet.cloud_cover
        let last_ice = planet.ice_cover
        let last_temp = planet.surf_temp
        let last_albedo = planet.albedo

        calculate_surface_temp(planet: planet, first: false,
                               last_water: last_water, last_clouds: last_clouds, last_ice: last_ice,
                               last_temp: last_temp, last_albedo: last_albedo, prng: prng)

        if fabs(planet.surf_temp - last_temp) < 0.25 {
            break
        }
    }

    planet.greenhs_rise = planet.surf_temp - initial_temp

//    if (flag_verbose & 0x20000)
//        fprintf (stderr, "%d: %5.1Lf gh = %5.1Lf (%5.1Lf C) st - %5.1Lf it [%5.1Lf re %5.1Lf a %5.1Lf alb]\n",
//                planet->planet_no,
//                planet->greenhs_rise,
//                planet->surf_temp,
//                planet->surf_temp - FREEZING_POINT_OF_WATER,
//                initial_temp,
//                planet->sun->r_ecosphere, planet->a, planet->albedo
//                );
}

/// *--------------------------------------------------------------------------*/
/// *     Inspired partial pressure, taking into account humidification of the    */
/// *     air in the nasal passage and throat This formula is on Dole's p. 14    */
/// *--------------------------------------------------------------------------*/
//
func inspired_partial_pressure(surf_pressure: Double, gas_pressure: Double) -> Double {
    let pH2O = H20_ASSUMED_PRESSURE
    let fraction = gas_pressure / surf_pressure
    return (surf_pressure - pH2O) * fraction
}

/// *--------------------------------------------------------------------------*/
/// *     This function uses figures on the maximum inspired partial pressures   */
/// *   of Oxygen, other atmospheric and traces gases as laid out on pages 15, */
/// *   16 and 18 of Dole's Habitable Planets for Man to derive breathability  */
/// *   of the planet's atmosphere.                                       JLB  */
/// *--------------------------------------------------------------------------*/

func breathability(planet: Planet) -> Breathability {
    var oxygen_ok = false

    if planet.gases == 0 {
        return .none
    }

    for a_gas in planet.atmosphere {
        let ipp = inspired_partial_pressure(surf_pressure: planet.surf_pressure, gas_pressure: a_gas.surf_pressure)
        if ipp > a_gas.type.max_ipp {
            return .poisonous
        }
        if a_gas.type.num == AN_O {
            oxygen_ok = ipp > MIN_O2_IPP && ipp <= MAX_O2_IPP
        }
    }

    if oxygen_ok {
        return .breathable
    }

    return .unbreathable
}

/// * function for 'soft limiting' temperatures */
func lim(_ x: Double) -> Double {
    x / sqrt(sqrt(1 + x * x * x * x))
}

func soft(_ v: Double, _ max: Double, _ min: Double) -> Double {
    let dv = v - min
    let dm = max - min
    return (lim(2 * dv / dm - 1) + 1) / 2 * dm + min
}

func set_temp_range(planet: Planet) {
    let pressmod = 1 / sqrt(1 + 20 * planet.surf_pressure / 1000.0)
    let ppmod = 1 / sqrt(10 + 5 * planet.surf_pressure / 1000.0)
    let tiltmod = fabs(cos(planet.axial_tilt * Double.pi / 180) * pow(1 + planet.e, 2))
    let daymod = 1 / (200 / planet.day + 1)
    let mh = pow(1 + daymod, pressmod)
    let ml = pow(1 - daymod, pressmod)
    let hi = mh * planet.surf_temp
    var lo = ml * planet.surf_temp
    let sh = hi + pow((100 + hi) * tiltmod, sqrt(ppmod))
    var wl = lo - pow((150 + lo) * tiltmod, sqrt(ppmod))
    let max = planet.surf_temp + sqrt(planet.surf_temp) * 10
    let min = planet.surf_temp / sqrt(planet.day + 24)

    if lo < min {
        lo = min
    }
    if wl < 0 {
        wl = 0
    }

    planet.high_temp = soft(hi, max, min)
    planet.low_temp = soft(lo, max, min)
    planet.max_temp = soft(sh, max, min)
    planet.min_temp = soft(wl, max, min)
}
