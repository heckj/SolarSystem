//
//  Accrete.swift
//

import Combine
import Foundation

/*------------------------------------------------------------------------*/
/*                             BIBLIOGRAPHY                               */
/*    Dole, Stephen H.  "Formation of Planetary Systems by Aggregation:   */
/*        a Computer Simulation"    October 1969,  Rand Corporation Paper */
/*        P-4226.                                                         */
/*------------------------------------------------------------------------*/

// NOTE(heckj): throughout this ported code, 'a' represents the AU distance from the core of the disk,
// 'mass' is the accumulated mass at that AU, and 'e' represents the eccentricity of the mass being calculated.

// NOTE(heckj): Porting Notes:
// Accrete was a variety of free functions that manipulated a set of global variables.
// I'm encapsulating those into a struct, which seems a bit awkward, but better than slapping
// around global variables and pointers to me.
// To support the linked-list data structures, I'm using reference types (classes) and passing them
// around as need be in place of passing "pointers" in the original C code. This is NOT the most efficient
// mechanism for Swift, and quite a lot could be refactored - I suspect most of what's happening in linked
// lists could be replaced with Array structures and updates, especially capturing the the dust lanes and
// collecting dust from them. The planet structure that accumulates ends up being a graph rather than a list
// when you include `moons`, so that may well best sit as a linked-list structure.
public struct AccretionState {
    public let dustlanes: [Dust]
    public let planets: [Planet]
    public let dust_left: Bool

    public static let example = AccretionState(dustlanes: [
        Dust(inner_edge: 0, outer_edge: 0.02963, dust_present: true, gas_present: true, next_band: nil),
        Dust(inner_edge: 0.02963, outer_edge: 0.3544, dust_present: false, gas_present: true, next_band: nil),
        Dust(inner_edge: 0.3544, outer_edge: 0.9784, dust_present: false, gas_present: false, next_band: nil),
        Dust(inner_edge: 0.9784, outer_edge: 29.59, dust_present: false, gas_present: false, next_band: nil),
        Dust(inner_edge: 29.59, outer_edge: 78.19, dust_present: false, gas_present: false, next_band: nil),
        Dust(inner_edge: 78.19, outer_edge: 212.5, dust_present: true, gas_present: true, next_band: nil),
    ], planets: [
        Planet(planet_no: 1, a: 0.3188, e: 0, axial_tilt: 0, mass: 0.9438, gas_giant: false, dust_mass: 0.9438, gas_mass: 0, next_planet: nil),
        Planet(planet_no: 2, a: 0.735, e: 0, axial_tilt: 0, mass: 1956, gas_giant: true, dust_mass: 1956, gas_mass: 0, next_planet: nil),
        Planet(planet_no: 3, a: 1.663, e: 0, axial_tilt: 0, mass: 1288, gas_giant: true, dust_mass: 1288, gas_mass: 0, next_planet: nil),
        Planet(planet_no: 4, a: 3.11, e: 0, axial_tilt: 0, mass: 882.2, gas_giant: true, dust_mass: 882.2, gas_mass: 0, next_planet: nil),
        Planet(planet_no: 5, a: 7.425, e: 0, axial_tilt: 0, mass: 504.2, gas_giant: true, dust_mass: 504.2, gas_mass: 0, next_planet: nil),
        Planet(planet_no: 6, a: 16.08, e: 0, axial_tilt: 0, mass: 1917, gas_giant: true, dust_mass: 1017, gas_mass: 0, next_planet: nil),
        Planet(planet_no: 7, a: 48.98, e: 0, axial_tilt: 0, mass: 1870, gas_giant: true, dust_mass: 1870, gas_mass: 0, next_planet: nil),
    ], dust_left: false)
}

public struct AccretionDisk {
    let FPStyle: FloatingPointFormatStyle<Double> = .number.precision(.significantDigits(1 ... 4))
    let tight: FloatingPointFormatStyle<Double> = .number.precision(.integerAndFractionLength(integerLimits: 1..., fractionLimits: 0 ... 3))

    var dust_left: Bool
    var cloud_eccentricity: Double
    var planet_inner_bound: Double
    var planet_outer_bound: Double

    var r_inner: Double = 0 // the inner radius of the computed band/width that a mass accumulates
    var r_outer: Double = 0 // the outer radius of the computed band/width that a mass accumulates
    var reduced_mass: Double = 0 // a reduced mass computed during accretion to collect dust and gas
    var dust_density: Double = 0

    public var stellar_mass_ratio: Double // the mass of the star, in solar masses
    public var stellar_luminosity_ratio: Double // the luminosity of the star, comparable to sol

    var dust_density_coeff: Double
    // NOTE(heckj): create a proper initializer to set this value, aligned somewhere between the initializer and calling
    // dist_planetary_masses, which is the primarily call site into this...
    var prng: RNGWrapper<Xoshiro>

    var dust_lanes: [Dust]
    var planets: [Planet]

    var dust_head: Dust?
    var planet_head: Planet?
    var do_moons: Bool

    var current_seed: Planet?
    public let updater: PassthroughSubject<AccretionState, Never>
    public let msgs: PassthroughSubject<String, Never>

    public init(prng: RNGWrapper<Xoshiro>,
                inner_limit_of_dust: Double,
                outer_limit_of_dust: Double,
                stellar_mass_ratio: Double,
                stellar_luminosity_ratio: Double,
                outer_planet_limit: Double? = nil,
                do_moons: Bool = true,
                dust_density_multipler: Double = 1,
                cloud_eccentricity: Double = 0.2,
                seed_system: Planet? = nil)
    {
        self.prng = prng
        self.do_moons = do_moons
        self.stellar_mass_ratio = stellar_mass_ratio
        self.stellar_luminosity_ratio = stellar_luminosity_ratio
        current_seed = seed_system

        // initial dust
        dust_density_coeff = DUST_DENSITY_COEFF * dust_density_multipler
        self.cloud_eccentricity = cloud_eccentricity
        dust_left = true
        if outer_limit_of_dust != 0 {
            let initialDustLane = Dust(
                inner_edge: inner_limit_of_dust,
                outer_edge: outer_limit_of_dust,
                dust_present: true, gas_present: true, next_band: nil
            )
            dust_head = initialDustLane
            dust_lanes = [initialDustLane]
        } else {
            let initialDustLane = Dust(
                inner_edge: inner_limit_of_dust,
                outer_edge: AccretionDisk.stellar_dust_limit(stellar_mass: stellar_mass_ratio),
                dust_present: true, gas_present: true, next_band: nil
            )
            dust_head = initialDustLane
            dust_lanes = [initialDustLane]
        }
        // initial planets
        planets = []
        planet_head = nil

        // Combine publishers to send updates about accretion status
        updater = PassthroughSubject<AccretionState, Never>()
        msgs = PassthroughSubject<String, Never>()

        // set or compute the range within which planets can be generated
        planet_inner_bound = AccretionDisk.nearest_planet(stellar_mass: stellar_mass_ratio)
        if let outer_planet_limit = outer_planet_limit, outer_planet_limit != 0.0 {
            planet_outer_bound = outer_planet_limit
        } else {
            planet_outer_bound = AccretionDisk.farthest_planet(stellar_mass: stellar_mass_ratio)
        }
    }

    /// The outer limit, in AU, of a planetary dust cloud for a star.
    /// - Parameter stellar_mass: The mass of the star, in solar masses
    /// - Returns: The distance, in AU, of the outermost range of dust for planetary accretion.
    public static func stellar_dust_limit(stellar_mass: Double) -> Double {
        200.0 * pow(stellar_mass, 1.0 / 3.0)
    }

    /// Computes the distance from the start to the closest interior planet that would likely form, in AU.
    /// - Parameter stellar_mass: The mass of the star, in solar masses
    /// - Returns: The distance to the nearest possilble planet, in AU.
    public static func nearest_planet(stellar_mass: Double) -> Double {
        0.3 * pow(stellar_mass, 1.0 / 3.0)
    }

    /// The outer limit, in AU, of a where a planet will likely accrete for a star.
    /// - Parameter stellar_mass: The mass of the star, in solar masses
    /// - Returns: The distance, in AU, of the outermost range of a planet for the given stellar mass.
    public static func farthest_planet(stellar_mass: Double) -> Double {
        50.0 * pow(stellar_mass, 1.0 / 3.0)
    }

    /// Returns the mass, in solar masses, at which a planet begins to accrete gas as well as dust.
    /// - Parameters:
    ///   - orb_radius: The orbital radius, in AU
    ///   - eccentricity: The eccentricity of the orbit.
    ///   - stell_luminosity_ratio: The luminosity of the star.
    public static func critical_limit(orbital_radius: Double, eccentricity: Double, stell_luminosity_ratio: Double) -> Double {
        let perihelion_dist = (orbital_radius - orbital_radius * eccentricity)
        let temp = perihelion_dist * sqrt(stell_luminosity_ratio)
        return B * pow(temp, -0.75)
    }

    /// Returns the inner most effect range for a given distance (a) in AU, and eccentricity that would accumulate dust and gas.
    /// - Parameters:
    ///   - a: The distance of the planetesimal mass from the star, in AU.
    ///   - e: The eccentricity of the orbit of the planetesimal mass.
    ///   - mass: The mass of the planetesimal.
    static func inner_effect_limit(a: Double, e: Double, mass: Double, cloud_eccentricity: Double) -> Double {
        // Given an orbital distance, mass, and eccentricity, calculates the
        // inner-most distance (in AU) that the orbit "sweeps up" dust.
        let minAU = a * (1.0 - e) * (1.0 - mass) / (1.0 + cloud_eccentricity)
        // clamp minimum value to 0 AU
        return minAU < 0.0 ? 0 : minAU
    }

    /// Returns the outer most effect range for a given distance (a) in AU, and eccentricity that would accumulate dust and gas.
    /// - Parameters:
    ///   - a: The distance of the planetesimal mass from the star, in AU.
    ///   - e: The eccentricity of the orbit of the planetesimal mass.
    ///   - mass: The mass of the planetesimal.
    static func outer_effect_limit(a: Double, e: Double, mass: Double, cloud_eccentricity: Double) -> Double {
        // Given an orbital distance, mass, and eccentricity, calculates the
        // outer-most distance (in AU) that the orbit "sweeps up" dust.
        a * (1.0 + e) * (1.0 + mass) / (1.0 - cloud_eccentricity)
    }

    /// Returns a Boolean value indicating if dust is available to accumulate within the given range.
    /// - Parameter sweeprange: The range, in AU, to check for available dust.
    func dust_available(_ sweeprange: ClosedRange<Double>) -> Bool {
        // called from `accrete_dust`
        dust_lanes.reduce(into: false) { partialResult, dustlane in
            if dustlane.range.overlaps(sweeprange) {
                partialResult = partialResult || dustlane.dust_present
            }
        }
    }

    func updated_dust_lanes(accretion_effect_range: ClosedRange<Double>, mass: Double, crit_mass: Double) -> [Dust] {
        // called from `accrete_dust`

        // With the collected mass accumulated earlier in sweeping through the dust lanes, we need to deplete the dust
        // stored within our data structures (the list of dust lanes) between the inner_bound and outer_bound
        // distances in AU.

        let min_effect = accretion_effect_range.lowerBound
        let max_effect = accretion_effect_range.upperBound
        let gas: Bool = (mass <= crit_mass) // Boolean to indicate we to deplete gas as well as dust

        let depleted_lanes = dust_lanes.flatMap { dust in
            if dust.inner_edge < min_effect && dust.outer_edge > max_effect {
                //  dust.inner ------------------------- dust.outer
                //                     min ---- max
                //  [  lane1          ][   lane2   ][      lane3   ]

                let lane1 = Dust(inner_edge: dust.inner_edge, outer_edge: min_effect, dust_present: dust.dust_present, gas_present: dust.gas_present, next_band: nil)
                let lane2 = Dust(inner_edge: min_effect, outer_edge: max_effect, dust_present: false, gas_present: dust.gas_present ? gas : false, next_band: nil)
                let lane3 = Dust(inner_edge: max_effect, outer_edge: dust.outer_edge, dust_present: dust.dust_present, gas_present: dust.gas_present, next_band: nil)
                return [lane1, lane2, lane3]
            } else if dust.inner_edge < max_effect && dust.outer_edge > max_effect {
                //         dust.inner ---------- dust.outer
                //   min -------------- max
                //         [    lane1     ][   lane2      ]

                let lane1 = Dust(inner_edge: dust.inner_edge, outer_edge: max_effect, dust_present: false, gas_present: dust.gas_present ? gas : false, next_band: nil)
                let lane2 = Dust(inner_edge: max_effect, outer_edge: dust.outer_edge, dust_present: dust.dust_present, gas_present: dust.gas_present, next_band: nil)
                return [lane1, lane2]
            } else if dust.inner_edge < min_effect && dust.outer_edge > min_effect {
                //        dust.inner ---------- dust.outer
                //                        min ----------------- max
                //        [    lane1     ][     lane2    ]

                let lane1 = Dust(inner_edge: dust.inner_edge, outer_edge: min_effect, dust_present: dust.dust_present, gas_present: dust.gas_present, next_band: nil)
                let lane2 = Dust(inner_edge: min_effect, outer_edge: dust.outer_edge, dust_present: false, gas_present: dust.gas_present ? gas : false, next_band: nil)
                return [lane1, lane2]
            } else if dust.inner_edge >= min_effect && dust.outer_edge <= max_effect {
                //       dust.inner --- dust.outer
                //  min ------------------------------ max
                //       [         lane1         ]
                let lane1 = Dust(inner_edge: dust.inner_edge, outer_edge: dust.outer_edge, dust_present: false, gas_present: dust.gas_present ? gas : false, next_band: nil)
                return [lane1]
            } else if dust.outer_edge < min_effect || dust.inner_edge > max_effect {
                //       node1.inner --- node1.outer
                //                                     min ----- max
                // OR
                //                       node1.inner --- node1.outer
                //       min ----- max
                // no overlap or effect on this lane:
                let lane1 = Dust(inner_edge: dust.inner_edge, outer_edge: dust.outer_edge, dust_present: dust.dust_present, gas_present: dust.gas_present, next_band: nil)
                return [lane1]
            }
            return []
        }

        let collapsed: [Dust] = depleted_lanes.reduce(into: []) { partialResult, dust in
            if partialResult.last != nil {
                let prev = partialResult.removeLast()
                if prev.dust_present == dust.dust_present, prev.gas_present == dust.gas_present {
                    // equivalent values, collapse them together
                    let newDust = Dust(inner_edge: min(prev.inner_edge, dust.inner_edge),
                                       outer_edge: max(prev.outer_edge, dust.outer_edge),
                                       dust_present: prev.dust_present, gas_present: prev.gas_present, next_band: nil)
                    partialResult.append(newDust)
                } else {
                    partialResult.append(prev)
                    partialResult.append(dust)
                }

            } else {
                // empty partialResult, so just pop in the current value - first in the list
                partialResult.append(dust)
            }
        }
        return collapsed
    }

    struct DustAndGas {
        var dust: Double
        var gas: Double
        func accretion_effect_range(a: Double, e: Double, cloud_eccentricity: Double) -> ClosedRange<Double> {
            // calculate the effect limits of dust collection based on a reduced mass
            let reduced_mass = pow(mass / (1.0 + mass), 1.0 / 4.0)
            return AccretionDisk.inner_effect_limit(a: a, e: e, mass: reduced_mass, cloud_eccentricity: cloud_eccentricity) ... AccretionDisk.outer_effect_limit(a: a, e: e, mass: reduced_mass, cloud_eccentricity: cloud_eccentricity)
        }

        var mass: Double {
            dust + gas
        }

        var reduced_mass: Double {
            pow(mass / (1.0 + mass), 1.0 / 4.0)
        }

        public static func + (lhs: DustAndGas, rhs: DustAndGas) -> DustAndGas {
            DustAndGas(dust: lhs.dust + rhs.dust, gas: lhs.gas + rhs.gas)
        }
    }

    /// Computes the mass collected during an accretion pass at the distance, and eccentricity you provide.
    /// - Parameters:
    ///   - mass: The mass, in solar masses, accreting dust and/or gas.
    ///   - critical_mass: The mass at which gas accretes as well as dust.
    ///   - a: The distance, in AU, of the orbit.
    ///   - e: The eccentricity of the orbit.
    /// - Returns: The mass, in solar masses, available to be gathered.
    func collect_dust(mass: Double, critical_mass: Double, a: Double, e: Double, dust_density: Double) -> DustAndGas {
        // sweeps through the existing dust lanes, collecting all the dust defined by the initial
        // orbit position (a: in AU) with an eccentricity (e), the boundary
        // markers of which are calculated by inner_effect_limit() and outer_effect_limit().

        let seed_dg = DustAndGas(dust: mass, gas: 0)
        let accretion_effect_range = seed_dg.accretion_effect_range(a: a, e: e, cloud_eccentricity: cloud_eccentricity)

        let massByLane: [DustAndGas] = dust_lanes.map { dustlane in
            let temp_density: Double
            let mass_density: Double
            let gas_density: Double

            if dustlane.dust_present == false {
                temp_density = 0.0
            } else {
                temp_density = dust_density
            }

            if mass < critical_mass || dustlane.gas_present == false {
                gas_density = 0.0
                mass_density = temp_density
            } else {
                mass_density = K * temp_density / (1.0 + sqrt(critical_mass / mass) * (K - 1.0))
                gas_density = mass_density - temp_density
            }

            if accretion_effect_range.overlaps(dustlane.inner_edge ... dustlane.outer_edge) {
                let bandwidth = accretion_effect_range.upperBound - accretion_effect_range.lowerBound
                var width = bandwidth

                var outer_width_reduction = accretion_effect_range.upperBound - dustlane.outer_edge
                if outer_width_reduction < 0 {
                    outer_width_reduction = 0
                }
                var inner_width_reduction = dustlane.inner_edge - accretion_effect_range.lowerBound
                if inner_width_reduction < 0 {
                    inner_width_reduction = 0
                }

                width -= outer_width_reduction
                width -= inner_width_reduction

                let volume = 4.0 * Double.pi * pow(a, 2.0) * mass * (1.0 - e * (outer_width_reduction - inner_width_reduction) / bandwidth) * width

                let accumulated_mass = volume * mass_density
                let gas_mass = volume * gas_density
                let dust_mass = accumulated_mass - gas_mass
                return DustAndGas(dust: dust_mass, gas: gas_mass)
            }
            return DustAndGas(dust: 0, gas: 0)
        }

        let accumulatedMass: DustAndGas = massByLane.reduce(into: DustAndGas(dust: 0, gas: 0)) { partialResult, dg in
            partialResult = partialResult + dg
        }

        return accumulatedMass
    }

    func updated_orbit(planet: Planet, planetesimal: DustAndGas, a: Double, e: Double) -> (a: Double, e: Double) {
        let new_a = (planet.mass + planetesimal.mass) / ((planet.mass / planet.a) + (planetesimal.mass / a))

        var temp = planet.mass * sqrt(planet.a) * sqrt(1.0 - pow(planet.e, 2.0))
        temp += (planetesimal.mass * sqrt(a) * sqrt(sqrt(1.0 - pow(e, 2.0))))
        temp /= ((planet.mass + planetesimal.mass) * sqrt(new_a))
        temp = 1.0 - pow(temp, 2.0)
        if temp < 0.0 || temp >= 1.0 {
            temp = 0.0
        }
        let new_e = sqrt(temp)
        return (new_a, new_e)
    }

    mutating func accrete_dust(seed_mass: DustAndGas, a: Double, e: Double, crit_mass: Double, dust_density: Double) -> DustAndGas {
        // called from `coalesce_planetesimals` and `dist_planetary_masses`/`advance`

        // This function is creating planetesimals (represented by a DustAndGas instance),
        // sweeping through dust lanes using 'collect_dust' until it stops growing at a notable rate.
        // Growth from the last sweep being less than 0.0001 times larger than the previous sweep.
        var temp_mass = seed_mass
        var new_mass: DustAndGas = temp_mass
        var growth: Double = 0
        repeat {
            temp_mass = new_mass
            new_mass = collect_dust(mass: temp_mass.mass, critical_mass: crit_mass, a: a, e: e, dust_density: dust_density)
            growth = new_mass.mass - temp_mass.mass
        } while !(growth < (0.0001 * temp_mass.mass))

        let grownMass = seed_mass + new_mass
        let accretion_effect_range = grownMass.accretion_effect_range(a: a, e: e, cloud_eccentricity: cloud_eccentricity)
        // MUTATION PARTS BELOW
        dust_lanes = updated_dust_lanes(accretion_effect_range: accretion_effect_range, mass: grownMass.mass, crit_mass: crit_mass)
        dust_left = dust_available(planet_inner_bound ... planet_outer_bound)
        return grownMass
    }

    func attempted_moon_capture(planet: Planet, planetesimal: DustAndGas, a: Double, e: Double, crit_mass: Double) -> (captured: Bool, Planet) {
        let existing_mass: Double = planet.moons.reduce(into: 0.0) { partialResult, moon in
            partialResult += moon.mass
        }
        // only planetesimals below critical mass can become moons
        if planetesimal.mass < crit_mass {
            let moon_capture = planetesimal.mass * SUN_MASS_IN_EARTH_MASSES < 2.5 &&
                planetesimal.mass * SUN_MASS_IN_EARTH_MASSES > 0.001 &&
                existing_mass < planet.mass * 0.05
            if moon_capture {
                if planetesimal.mass > planet.mass {
                    // captured planet out-masses the existing planet, so the existing
                    // planet *becomes* the moon.
                    let new_moon = Planet(
                        planet_no: 0,
                        a: a,
                        e: e,
                        axial_tilt: planet.axial_tilt,
                        mass: planet.mass,
                        gas_giant: false,
                        dust_mass: planet.dust_mass,
                        gas_mass: planet.gas_mass,
                        next_planet: nil
                    )
                    planet.dust_mass = planetesimal.dust
                    planet.gas_mass = planetesimal.gas
                    planet.mass = planetesimal.mass
                    planet.moons.append(new_moon)
                } else {
                    let new_moon = Planet(
                        planet_no: 0,
                        a: a,
                        e: e,
                        axial_tilt: 0,
                        mass: planetesimal.mass,
                        gas_giant: false,
                        dust_mass: planetesimal.dust,
                        gas_mass: planetesimal.gas,
                        next_planet: nil
                    )
                    planet.moons.append(new_moon)
                }
                print("Moon captured: \(planet.a) AU \(planet.mass * SUN_MASS_IN_EARTH_MASSES)EM <- \(planetesimal.mass * SUN_MASS_IN_EARTH_MASSES)EM")
                return (true, planet)
            } else {
                // moon escape
                var escape_reason = ""
                if (planetesimal.mass * SUN_MASS_IN_EARTH_MASSES) >= 2.5 {
                    escape_reason = ", too big"
                } else if (planetesimal.mass * SUN_MASS_IN_EARTH_MASSES) <= 0.0001 {
                    escape_reason = ", too small"
                }
                print("Moon escapes: \(planet.a) AU (\(planet.mass * SUN_MASS_IN_EARTH_MASSES))EM \(planetesimal.mass * SUN_MASS_IN_EARTH_MASSES) EM\(escape_reason)")
                return (false, planet)
            }
        }
        // planetesimal has achieved critical mass - gathers gas as well as dust
        return (false, planet)
    }

    mutating func collide(planet: Planet, planetesimal: DustAndGas, a: Double, e: Double) -> (Planet) {
        let crit_mass = AccretionDisk.critical_limit(orbital_radius: a, eccentricity: e, stell_luminosity_ratio: stellar_luminosity_ratio)

        print("Collision between planet and planetesimal! \(planet.a) AU (\(planet.mass * SUN_MASS_IN_EARTH_MASSES)EM) + \(a) AU (\(planetesimal.mass * SUN_MASS_IN_EARTH_MASSES)EM = \(planetesimal.dust * SUN_MASS_IN_EARTH_MASSES)EMd + \(planetesimal.gas * SUN_MASS_IN_EARTH_MASSES)EMg [\(crit_mass * SUN_MASS_IN_EARTH_MASSES)EM]) -> \(a) AU (\(e))")

        let combined = DustAndGas(dust: planet.dust_mass + planetesimal.dust, gas: planet.gas_mass + planetesimal.gas)

        let local_dust_density = dust_density_coeff * sqrt(stellar_mass_ratio) * exp(-ALPHA * pow(a, 1.0 / N))
        // NOTE(heckj): accrete_dust has the side effect of updating the dust lanes...
        let updated_mass = accrete_dust(seed_mass: combined, a: a, e: e, crit_mass: crit_mass, dust_density: local_dust_density)

        planet.a = a
        planet.e = e
        planet.dust_mass = updated_mass.dust
        planet.gas_mass = updated_mass.gas
        planet.mass = updated_mass.mass
        if updated_mass.mass >= crit_mass {
            planet.gas_giant = true
        }
        return planet
    }

    mutating func coalesce_planetesimals(planetesimal: DustAndGas, a: Double, e: Double, crit_mass: Double, do_moons: Bool) -> [Planet] {
        // Takes the provided planetesimal and iterates through, colliding or merging it into the list of existing planets
        // as appropriate.
        var planets = planets
        var finished = false

        // First we try to find an existing planet with an over-lapping orbit.
        if let planet_with_overlaping_orbit = planets.first(where: { planet in
            planet.a > a
        }) {
            if let index_first_overlapping = planets.firstIndex(of: planet_with_overlaping_orbit),
               index_first_overlapping > 0
            {
                let previous_planet = planets[index_first_overlapping - 1]

                // prev planet calculations

                // planetesimal orbit with extreme eccentricity
                let max_dist_to_prev = (a * (1.0 + e) * (1.0 - planetesimal.reduced_mass))

                /* previous planet's perihelion based on eccentricity */
                let min_dist_to_prev = (previous_planet.a * (1.0 + previous_planet.e) * (1.0 + previous_planet.reduced_mass)) - previous_planet.a

                let mean_dist_to_prev = a - previous_planet.a

                // fabs(_) returns the absolute value of the provided floating point number
                if fabs(mean_dist_to_prev) <= fabs(max_dist_to_prev) || fabs(mean_dist_to_prev) <= fabs(min_dist_to_prev) {
                    // collision with orbital flow of the previous planet

                    let (new_a, new_e) = updated_orbit(planet: previous_planet, planetesimal: planetesimal, a: a, e: e)
                    // NOW DO SOMETHING WITH THE NEW a and e values
                    if do_moons {
                        let (captured, _) = attempted_moon_capture(planet: planet_with_overlaping_orbit, planetesimal: planetesimal, a: new_a, e: new_e, crit_mass: crit_mass)
                        if captured {
                            finished = true
                        }
                    }

                    if !finished {
                        // collide previous planet and planetesimal
                        // collide has the side effect of potential additional accretion, and updating the dust lanes
                        _ = collide(planet: previous_planet, planetesimal: planetesimal, a: new_a, e: new_e)
                        finished = true
                    }
                }
                // no intersection with previous planets, continue looking for overlap and possible collision
            }

            if !finished {
                // if finished, then we collided and interacted with the previous planet and shouldn't check or attempt
                // capture or collision with the next planet in the sequence.

                // next planet calculations
                let mean_dist_to_next = planet_with_overlaping_orbit.a - a

                // planetesimals' maximum orbital swing with eccentricity
                let max_dist_to_next = (a * (1.0 + e) * (1.0 - planetesimal.reduced_mass))

                /* next planets aphelion based on eccentricity */
                let min_dist_to_next = planet_with_overlaping_orbit.a - (planet_with_overlaping_orbit.a * (1.0 - planet_with_overlaping_orbit.e) * (1.0 - planet_with_overlaping_orbit.reduced_mass))

                // fabs(_) returns the absolute value of the provided floating point number
                if fabs(mean_dist_to_next) <= fabs(max_dist_to_next) || fabs(mean_dist_to_next) <= fabs(min_dist_to_next) {
                    // collision with orbital flow of the next planet
                    let (new_a, new_e) = updated_orbit(planet: planet_with_overlaping_orbit, planetesimal: planetesimal, a: a, e: e)
                    // NOW DO SOMETHING WITH THE NEW a and e values

                    if do_moons {
                        let (captured, _) = attempted_moon_capture(planet: planet_with_overlaping_orbit, planetesimal: planetesimal, a: new_a, e: new_e, crit_mass: crit_mass)
                        if captured {
                            finished = true
                        }
                    }

                    if !finished {
                        // collide next planet and planetesimal
                        // collide has the side effect of potential additional accretion, and updating the dust lanes
                        _ = collide(planet: planet_with_overlaping_orbit, planetesimal: planetesimal, a: new_a, e: new_e)
                        finished = true
                    }
                }
                // no intersection/collision with the next planet in the sequence
            } else {
                // No planet exists with an overlapping orbit
                // Planetesimals didn't collide, so we make it into a planet.
                let the_planet = Planet(planet_no: 0, a: a, e: e, axial_tilt: 0, mass: planetesimal.mass,
                                        gas_giant: planetesimal.mass >= crit_mass,
                                        dust_mass: planetesimal.dust, gas_mass: planetesimal.gas,
                                        next_planet: nil)
                if let insert_index = planets.firstIndex(where: { planet in
                    planet.a > a
                }) {
                    planets.insert(the_planet, at: insert_index)
                } else {
                    planets.append(the_planet)
                }
            }
        } else {
            // No planet exists with an overlapping orbit
            // Planetesimals didn't collide, so we make it into a planet.
            let the_planet = Planet(planet_no: 0, a: a, e: e, axial_tilt: 0, mass: planetesimal.mass,
                                    gas_giant: planetesimal.mass >= crit_mass,
                                    dust_mass: planetesimal.dust, gas_mass: planetesimal.gas,
                                    next_planet: nil)
            planets.append(the_planet)
        }

        return planets
    }

    public func currentState() -> AccretionState {
        var current_dust_lanes: [Dust] = []
        var current_planets: [Planet] = []

        if let firstDustLane = dust_head {
            for dust_lane in sequence(first: firstDustLane, next: \.next_band) {
                current_dust_lanes.append(dust_lane)
            }
        }

        if let firstPlanet = planet_head {
            for some_planet in sequence(first: firstPlanet, next: \.next_planet) {
                current_planets.append(some_planet)
            }
        }
        return AccretionState(dustlanes: current_dust_lanes, planets: current_planets, dust_left: dust_left)
    }

    public mutating func advance() {
        let a: Double // distance, in AU
        let e: Double // eccentricity of orbit
        let seed_dg = DustAndGas(dust: PROTOPLANET_MASS, gas: 0)

        // The general flow is to seed planetesimals that start gravitational accretion
        // (called 'gravitational instability' in current (~2020) astrophysics research)
        // and continue to deploy "seeds" until the dust of the disk is "consumed".
        if dust_left {
            if let definitely_seed = current_seed {
                a = definitely_seed.a
                e = definitely_seed.e
                current_seed = definitely_seed.next_planet
            } else {
                a = prng.random_number(in: planet_inner_bound ... planet_outer_bound)
                e = prng.random_eccentricity()
            }
            // print("Checking \(a.formatted(FPStyle)) AU")
            msgs.send("Checking \(a.formatted(FPStyle)) AU")

            let accretion_effect_range = seed_dg.accretion_effect_range(a: a, e: e, cloud_eccentricity: cloud_eccentricity)
            if dust_available(accretion_effect_range) {
                // print("Injecting protoplanet at \(a.formatted(FPStyle)) AU")
                msgs.send("Injecting protoplanet at \(a.formatted(FPStyle)) AU")

                // N = 3, ALPHA = 5
                dust_density = dust_density_coeff * sqrt(stellar_mass_ratio) * exp(-ALPHA * pow(a, 1.0 / N))

                // Determine the mass (in solar masses) at which a body will start accumulating gasses
                let crit_mass = AccretionDisk.critical_limit(orbital_radius: a, eccentricity: e, stell_luminosity_ratio: stellar_luminosity_ratio)

                let planetesimal = accrete_dust(seed_mass: seed_dg, a: a, e: e, crit_mass: crit_mass, dust_density: dust_density)

                if planetesimal.mass > PROTOPLANET_MASS {
                    let updated_planets = coalesce_planetesimals(planetesimal: planetesimal, a: a, e: e, crit_mass: crit_mass, do_moons: do_moons)
                    planets = updated_planets.sorted()
                } else {
                    msgs.send("Insufficient mass additional from accretion.")
                }
            } // dust available for relevant distance (a) and eccentricity (e)
            else {
                msgs.send("No dust available at \(a.formatted(FPStyle)) AU")
            }
            updater.send(currentState())
//            for dust in final_state.dustlanes {
//                let dust_symbols = "\(dust.dust_present ? "+" : " ")\(dust.gas_present ? "." : " ")"
//                print("  \(dust_symbols) \(dust.inner_edge.formatted(FPStyle)) - \(dust.outer_edge.formatted(FPStyle))")
//            }
//            for planet in final_state.planets {
//                print(" \(planet.a.formatted(FPStyle)) AU : \(planet.id) \( (planet.mass * SUN_MASS_IN_EARTH_MASSES).formatted(FPStyle)) EM")
//            }
        }
    }

    // primary entry point? - called from SolarSystem
    mutating func dist_planetary_masses() -> Planet? {
        let seed_dg = DustAndGas(dust: PROTOPLANET_MASS, gas: 0)

        // The general flow is to seed planetesimals that start gravitational accretion
        // (called 'gravitational instability' in current (~2020) astrophysics research)
        // and continue to deploy "seeds" until the dust of the disk is "consumed".
        while dust_left {
            let a: Double // distance, in AU
            let e: Double // eccentricity of orbit

            print("DUST STATUS")
            updater.send(currentState())

            if let definitely_seed = current_seed {
                a = definitely_seed.a
                e = definitely_seed.e
                current_seed = definitely_seed.next_planet
            } else {
                a = prng.random_number(in: planet_inner_bound ... planet_outer_bound)
                e = prng.random_eccentricity()
            }

            print("Checking \(a.formatted(FPStyle)) AU")
            let accretion_effect_range = seed_dg.accretion_effect_range(a: a, e: e, cloud_eccentricity: cloud_eccentricity)
            if dust_available(accretion_effect_range) {
                print("Injecting protoplanet at \(a.formatted(FPStyle)) AU")

                dust_density = dust_density_coeff * sqrt(stellar_mass_ratio) * exp(-ALPHA * pow(a, 1.0 / N))

                // Determine the mass (in solar masses) at which a body will start accumulating gasses
                let crit_mass = AccretionDisk.critical_limit(orbital_radius: a, eccentricity: e, stell_luminosity_ratio: stellar_luminosity_ratio)

                let planetesimal = accrete_dust(seed_mass: seed_dg, a: a, e: e, crit_mass: crit_mass, dust_density: dust_density)

                if planetesimal.mass > PROTOPLANET_MASS {
                    let updated_planets = coalesce_planetesimals(planetesimal: planetesimal, a: a, e: e, crit_mass: crit_mass, do_moons: do_moons)
                    planets = updated_planets.sorted()
                } else {
                    msgs.send("Insufficient mass additional from accretion.")
                }
            }
        }

        let final_state = currentState()
        print("DUST CONSUMED")
        for dust in final_state.dustlanes {
            let dust_symbols = "\(dust.dust_present ? "+" : " ")\(dust.gas_present ? "." : " ")"
            print("  \(dust_symbols) \(dust.inner_edge.formatted(FPStyle)) - \(dust.outer_edge.formatted(FPStyle))")
        }
        for planet in final_state.planets {
            print(" \(planet.a.formatted(FPStyle)) AU : \(planet.id) \((planet.mass * SUN_MASS_IN_EARTH_MASSES).formatted(FPStyle)) EM")
        }
        return (planet_head)
    }
}
