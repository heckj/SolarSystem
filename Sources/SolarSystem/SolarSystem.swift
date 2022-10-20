import Foundation

public struct SolarSystem {
    var sun: Sun
    var planets: [Planet]
    var seed: UInt64
    
    func text_describe_system() {
        var planet: Planet?
        var counter = 1

        print("SolarSystem Generation seed: \(seed)")
        print("                          SYSTEM  CHARACTERISTICS")
        print("Stellar mass: \(sun.mass.formatted(FPStyle)) solar masses")
        print("Stellar luminosity: \(sun.luminosity.formatted(FPStyle))")
        print("Age: \((sun.age / 1.0e9).formatted(FPStyle)) billion years (\(((sun.life - sun.age) / 1.0e9).formatted(FPStyle)) billion left on main sequence)")
        print("Habitable ecosphere radius: \(sun.r_ecosphere.formatted(FPStyle)) AU")
        print()
        print("Planets present at:")
        for Aplanet in planets {
            let textSymbol: String
            if Aplanet.gas_giant {
                textSymbol = "O"
            } else if Aplanet.greenhouse_effect, Aplanet.surf_pressure > 0.0 {
                textSymbol = "+"
            } else if Aplanet.hydrosphere > 0.05, Aplanet.hydrosphere < 0.95 {
                textSymbol = "*"
            } else if (Aplanet.mass * SUN_MASS_IN_EARTH_MASSES) > 0.1 {
                textSymbol = "o"
            } else {
                textSymbol = "."
            }
            print("\(counter)\t\(textSymbol)\t\(Aplanet.a.formatted(AUFormat)) AU (\(Aplanet.e.formatted(eFormat))) \(Aplanet.orbital_range())\t\((Aplanet.mass * SUN_MASS_IN_EARTH_MASSES).formatted(massFormat)) EM\n")
            counter += 1
            planet = planet?.next_planet
        }

        print()

        for Aplanet in planets {
            print("Planet \(Aplanet.id) \(Aplanet.gas_giant ? "*gas giant*" : "")")

            if Int(Aplanet.day) == Int(Aplanet.orb_period * 24.0) {
                print("Planet is tidally locked with one face to star.")
            }
            if Aplanet.resonant_period {
                print("Planet's rotation is in a resonant spin lock with the star")
            }
            print("   Distance from primary star:\t\(Aplanet.a.formatted(FPStyle))\tAU")
            print("   Mass:\t\t\t\((Aplanet.mass * SUN_MASS_IN_EARTH_MASSES).formatted(FPStyle))\tEarth masses")
            if !(Aplanet.gas_giant) {
                print("   Surface gravity:\t\t\(Aplanet.surf_grav.formatted(FPStyle))\tEarth Gs")
                print("   Surface pressure:\t\t\((Aplanet.surf_pressure / 1000.0).formatted(FPStyle))\tEarth atmospheres")
                if Aplanet.greenhouse_effect, Aplanet.surf_pressure > 0.0 {
                    print("\tGREENHOUSE EFFECT")
                }
                print("   Surface temperature:\t\t\((Aplanet.surf_temp - FREEZING_POINT_OF_WATER).formatted(FPStyle))\tdegrees Celcius\n")
            }
            print("   Equatorial radius:\t\t\(Aplanet.radius.formatted(FPStyle))\tKm")
            print("   Density:\t\t\t\(Aplanet.density.formatted(FPStyle))\tgrams/cc")
            print("   Eccentricity of orbit:\t\(Aplanet.e.formatted(FPStyle))")
            print("   Escape Velocity:\t\t\((Aplanet.esc_velocity / CM_PER_KM).formatted(FPStyle))\tKm/sec")
            print("   Molecular weight retained:\t\(Aplanet.molec_weight.formatted(FPStyle)) and above")
            print("   Surface acceleration:\t\(Aplanet.surf_accel.formatted(FPStyle))\tcm/sec2")
            print("   Axial tilt:\t\t\t\(Aplanet.axial_tilt.formatted(FPStyle))\tdegrees")
            print("   Planetary albedo:\t\t\(Aplanet.albedo.formatted(FPStyle))")
            print("   Length of year:\t\t\(Aplanet.orb_period.formatted(FPStyle))\tdays")
            print("   Length of day:\t\t\(Aplanet.day.formatted(FPStyle))\thours")
            if !(Aplanet.gas_giant) {
                print("   Boiling point of water:\t\((Aplanet.boil_point - FREEZING_POINT_OF_WATER).formatted(FPStyle))\tdegrees Celcius")
                print("   Hydrosphere percentage:\t\((Aplanet.hydrosphere * 100.0).formatted(FPStyle))")
                print("   Cloud cover percentage:\t\((Aplanet.cloud_cover * 100).formatted(FPStyle))")
                print("   Ice cover percentage:\t\((Aplanet.ice_cover * 100).formatted(FPStyle))")
            }
            print()

//            if do_gases,
//               Aplanet.planet_type != .gasgiant,
//               Aplanet.planet_type != .subgasgiant,
//               Aplanet.planet_type != .subsubgasgiant
//            {
//                // gases?
//            }
        }
    }
}

let FPStyle: FloatingPointFormatStyle<Double> = .number.precision(.significantDigits(1 ... 4))

// called from stargen()
func generate_stellar_system(sun: Sun,
                             use_seed_system: Bool,
                             seed_system: Planet?,
                             sys_no: Int,
                             system_name: String,
                             outer_planet_limit: Double,
                             dust_density_coeff _: Double,
                             do_gases: Bool,
                             do_moons: Bool,
                             prng: RNGWrapper<Xoshiro>,
                             counts: inout InterestingCounts) -> SolarSystem
{
    let outer_dust_limit = AccretionDisk.stellar_dust_limit(stellar_mass: sun.mass)

    // NOTE(heckj): This only invokes the generation sequence IF `use_seed_system` is false,
    // which I suspect means that it was an indicator that planetary refinement should be processed
    // against a known set (such as our solar system) without generating a sequence of new planets
    // in order to explore and verify those algorithms.
    if use_seed_system {
        sun.age = 5.0e9
    } else {
        let min_age = 1.0e9
        let max_age = 6.0e9

        if sun.life < max_age {
            sun.life = max_age
        }
        sun.age = prng.random_number(in: min_age ... max_age)
    }

    var accretionDisk = AccretionDisk(prng: prng, inner_limit_of_dust: 0.0, outer_limit_of_dust: outer_dust_limit, stellar_mass_ratio: sun.mass, stellar_luminosity_ratio: sun.luminosity, outer_planet_limit: outer_planet_limit, do_moons: do_moons, dust_density_multipler: 1, seed_system: seed_system)

    let planets = accretionDisk.dist_planetary_masses()
    generate_planet_details(sun: sun,
                            planets: planets,
                            random_tilt: !use_seed_system,
                            sys_no: sys_no,
                            system_name: system_name,
                            do_gases: do_gases,
                            do_moons: do_moons,
                            prng: prng,
                            counts: &counts)

    print(" .. GENERATION FINISHED")
    // At this point, we have the star for the system defined in sun,
    // the linked-list of planets in innermost_planet,
    // and the counts of planet types, number of breathable atmospheres, habitable, etc in counts.
    return SolarSystem(sun: sun, planets: planets, seed: prng.seed)
}

func calculate_gases(sun: Sun, planet: Planet, planet_id: String) {
    if planet.surf_pressure > 0 {
        var total_amount: Double = 0
        let pressure: Double = planet.surf_pressure / MILLIBARS_PER_BAR
        var n: Int = 0
        for gas in gases {
            let yp = gas.boil / (373.0 * ((log(pressure + 0.001) / -5050.5) + (1.0 / 373.0)))
            if yp >= 0, yp < planet.low_temp, gas.weight >= planet.molec_weight {
                let vrms = rms_vel(molecular_weight: gas.weight, exospheric_temp: planet.exospheric_temp)
                let pvrms = pow(1 / (1 + vrms / planet.esc_velocity), sun.age / 1e9)
                var abund = gas.abunds
                var react = 1.0
                var fract = 1.0
                var pres2 = 1.0

                if gas.symbol == "Ar" {
                    react = 0.15 * sun.age / 4e9
                } else if gas.symbol == "He" {
                    abund = abund * (0.001 + (planet.gas_mass / planet.mass))
                    pres2 = (0.75 + pressure)
                    react = pow(1 / (1 + gas.reactivity), sun.age / 2e9 * pres2)
                } else if gas.symbol == "O" || gas.symbol == "O2", sun.age > 2e9, planet.surf_temp > 270, planet.surf_temp < 400 {
                    /*    pres2 = (0.65 + pressure/2);            Breathable - M: .55-1.4     */
                    pres2 = (0.89 + pressure / 4) /*    Breathable - M: .6 -1.8     */
                    react = pow(1 / (1 + gas.reactivity), pow(sun.age / 2e9, 0.25) * pres2)

                } else if gas.symbol == "CO2", sun.age > 2e9, planet.surf_temp > 270, planet.surf_temp < 400 {
                    pres2 = (0.75 + pressure)
                    react = pow(1 / (1 + gas.reactivity), pow(sun.age / 2e9, 0.5) * pres2)
                    react *= 1.5
                } else {
                    pres2 = (0.75 + pressure)
                    react = pow(1 / (1 + gas.reactivity), sun.age / 2e9 * pres2)
                }

                fract = (1 - (planet.molec_weight / gas.weight))

                let amount = abund * pvrms * react * fract
                if amount > 0 {
                    n += 1
                }
                total_amount += amount
                planet.atmosphere.append(Gas(gas, surf_pressure: amount)) // stash the amount into surf_pressure for the
                // gas, to be updated after we've computed the total amount of gas to arrange the appropriate surface
                // pressures.
                if ["O", "N", "Ar", "He", "CO2"].contains(gas.symbol) {
                    print(" .. .. .. \((planet.mass * SUN_MASS_IN_EARTH_MASSES).formatted(FPStyle)) \(gas.symbol), \(amount.formatted(FPStyle)) = \(abund.formatted(FPStyle)) * \(pvrms.formatted(FPStyle)) * \(react.formatted(FPStyle)) * \(pres2.formatted(FPStyle)) * \(fract.formatted(FPStyle)) (ratio \((planet.gas_mass / planet.mass).formatted(.percent)))")
                }
            }
        } // foreach gas

        // update
        for gas in planet.atmosphere {
            let placeholder_amount = gas.surf_pressure
            gas.surf_pressure = planet.surf_pressure * placeholder_amount / total_amount

            if gas.type.num == AN_O, inspired_partial_pressure(surf_pressure: planet.surf_pressure, gas_pressure: gas.surf_pressure) > gas.type.max_ipp {
                print("\(planet_id) Poisoned by O2")
            }
        }
        planet.atmosphere.sort { gasA, gasB in
            gasA.surf_pressure > gasB.surf_pressure
        }
        print("\(planet_id), (\(planet.a.formatted(FPStyle)) AU) gases:")
        for gas in planet.atmosphere {
            print("    \(gas.type.symbol): \(gas.surf_pressure.formatted(FPStyle)), \((gas.surf_pressure / planet.surf_pressure).formatted(.percent))")
        }
    }
}

func generate_planet(planet: Planet,
                     planet_no: Int,
                     sun: Sun,
                     random_tilt: Bool,
                     planet_id: String,
                     do_gases: Bool,
                     do_moons: Bool,
                     is_moon: Bool,
                     prng: RNGWrapper<Xoshiro>,
                     counts: inout InterestingCounts)
{
    planet.planet_no = planet_no
    planet.sun = sun

    planet.orbit_zone = orbital_zone(luminosity: sun.luminosity, orbital_radius: planet.a)
    planet.orb_period = period(separation: planet.a, small_mass: planet.mass, large_mass: sun.mass)
    if random_tilt {
        planet.axial_tilt = inclination(orb_radius: planet.a, prng: prng)
    }

    planet.exospheric_temp = EARTH_EXOSPHERE_TEMP / pow2(planet.a / sun.r_ecosphere)
    planet.rms_velocity = rms_vel(molecular_weight: MOL_NITROGEN, exospheric_temp: planet.exospheric_temp)
    planet.core_radius = kothari_radius(mass: planet.dust_mass, giant: false, zone: planet.orbit_zone)

    // Calculate the radius as a gas giant, to verify it will retain gas.
    // Then if mass > Earth, it's at least 5% gas and retains He, it's
    // some flavor of gas giant.

    planet.density = empirical_density(mass: planet.mass, orb_radius: planet.a, r_ecosphere: sun.r_ecosphere, gas_giant: true)
    planet.radius = volume_radius(mass: planet.mass, density: planet.density)

    planet.surf_accel = acceleration(mass: planet.mass, radius: planet.radius)
    planet.surf_grav = gravity(acceleration: planet.surf_accel)

    planet.molec_weight = min_molec_weight(planet: planet)

    if (planet.mass * SUN_MASS_IN_EARTH_MASSES) > 1.0,
       (planet.gas_mass / planet.mass) > 0.05,
       min_molec_weight(planet: planet) <= 4.0
    {
        print(" .. .. Gas Giant with mass ratio (gas/mass) of \((planet.gas_mass / planet.mass).formatted(FPStyle))")
        if (planet.gas_mass / planet.mass) < 0.20 {
            planet.planet_type = .subsubgasgiant
        } else if (planet.mass * SUN_MASS_IN_EARTH_MASSES) < 20.0 {
            planet.planet_type = .subgasgiant
        } else {
            planet.planet_type = .gasgiant
        }
    } else {
        // If not, it's rocky.
        print(" .. .. Rocky Planet")
        planet.radius = kothari_radius(mass: planet.mass, giant: false, zone: planet.orbit_zone)
        planet.density = volume_density(mass: planet.mass, equat_radius: planet.radius)

        planet.surf_accel = acceleration(mass: planet.mass, radius: planet.radius)
        planet.surf_grav = gravity(acceleration: planet.surf_accel)
        if (planet.gas_mass / planet.mass) > 0.000001 {
            let h2_mass = planet.gas_mass * 0.85
            let he_mass = (planet.gas_mass - h2_mass) * 0.999

            var h2_loss = 0.0
            var he_loss = 0.0
            print(" .. .. atmosphere (gas/mass ratio: \((planet.gas_mass / planet.mass).formatted(FPStyle)))")
            let h2_life = gas_life(molecular_weight: MOL_HYDROGEN, planet: planet)
            let he_life = gas_life(molecular_weight: HELIUM, planet: planet)

            if h2_life < sun.age {
                h2_loss = ((1.0 - (1.0 / exp(sun.age / h2_life))) * h2_mass)

                planet.gas_mass -= h2_loss
                planet.mass -= h2_loss

                planet.surf_accel = acceleration(mass: planet.mass, radius: planet.radius)
                planet.surf_grav = gravity(acceleration: planet.surf_accel)
            }

            if he_life < sun.age {
                he_loss = ((1.0 - (1.0 / exp(sun.age / he_life))) * he_mass)

                planet.gas_mass -= he_loss
                planet.mass -= he_loss

                planet.surf_accel = acceleration(mass: planet.mass, radius: planet.radius)
                planet.surf_grav = gravity(acceleration: planet.surf_accel)
            }
            if (h2_loss + he_loss) > 0.000001 {
                print(" .. .. \(planet_id) losing gas: H2 \((h2_loss * SUN_MASS_IN_EARTH_MASSES).formatted(FPStyle)) EM, He \((he_loss * SUN_MASS_IN_EARTH_MASSES).formatted(FPStyle)) EM")
            }
        }

        planet.day = day_length(planet: planet) /* Modifies planet->resonant_period */
        planet.esc_velocity = escape_vel(mass: planet.mass, radius: planet.radius)

        if (planet.planet_type == .gasgiant)
            || (planet.planet_type == .subgasgiant)
            || (planet.planet_type == .subsubgasgiant)
        {
            planet.greenhouse_effect = false
            planet.volatile_gas_inventory = INCREDIBLY_LARGE_NUMBER
            planet.surf_pressure = INCREDIBLY_LARGE_NUMBER

            planet.boil_point = INCREDIBLY_LARGE_NUMBER

            planet.surf_temp = INCREDIBLY_LARGE_NUMBER
            planet.greenhs_rise = 0
            planet.albedo = prng.about(GAS_GIANT_ALBEDO, variation: 0.1)
            planet.hydrosphere = 1.0
            planet.cloud_cover = 1.0
            planet.ice_cover = 0.0
            planet.surf_grav = gravity(acceleration: planet.surf_accel)
            planet.molec_weight = min_molec_weight(planet: planet)
            planet.surf_grav = INCREDIBLY_LARGE_NUMBER
            planet.estimated_temp = est_temp(ecosphere_radius: sun.r_ecosphere, orb_radius: planet.a, albedo: planet.albedo)
            planet.estimated_terr_temp = est_temp(ecosphere_radius: sun.r_ecosphere, orb_radius: planet.a, albedo: EARTH_ALBEDO)

            if planet.estimated_terr_temp >= FREEZING_POINT_OF_WATER,
               planet.estimated_terr_temp <= EARTH_AVERAGE_KELVIN + 10.0,
               sun.age > 2.0e9
            {
                counts.habitable_jovians += 1

                print(" .. .. \(planet_id) (\((planet.mass * SUN_MASS_IN_EARTH_MASSES).formatted(FPStyle))EM, age \((sun.age / 1.0e9).formatted(FPStyle))) at \(planet.estimated_terr_temp - FREEZING_POINT_OF_WATER)°C")
            }
        } else {
            // rocky world, not gas giant
            planet.estimated_temp = est_temp(ecosphere_radius: sun.r_ecosphere, orb_radius: planet.a, albedo: EARTH_ALBEDO)
            planet.estimated_terr_temp = est_temp(ecosphere_radius: sun.r_ecosphere, orb_radius: planet.a, albedo: EARTH_ALBEDO)

            planet.surf_grav = gravity(acceleration: planet.surf_accel)
            planet.molec_weight = min_molec_weight(planet: planet)

            planet.greenhouse_effect = grnhouse(r_ecosphere: sun.r_ecosphere, orb_radius: planet.a)
            planet.volatile_gas_inventory = vol_inventory(mass: planet.mass,
                                                          escape_vel: planet.esc_velocity,
                                                          rms_vel: planet.rms_velocity,
                                                          stellar_mass: sun.mass,
                                                          zone: planet.orbit_zone,
                                                          greenhouse_effect: planet.greenhouse_effect,
                                                          accreted_gas: (planet.gas_mass
                                                              / planet.mass) > 0.000001,
                                                          prng: prng)
            planet.surf_pressure = pressure(volatile_gas_inventory: planet.volatile_gas_inventory,
                                            equat_radius: planet.radius,
                                            gravity: planet.surf_grav)

            if planet.surf_pressure == 0.0 {
                planet.boil_point = 0.0
            } else {
                planet.boil_point = boiling_point(surf_pressure: planet.surf_pressure)
            }

            iterate_surface_temp(planet: planet, prng: prng) /*    Sets:
             *        planet->surf_temp
             *        planet->greenhs_rise
             *        planet->albedo
             *        planet->hydrosphere
             *        planet->cloud_cover
             *        planet->ice_cover
             */

            if do_gases,
               planet.max_temp >= FREEZING_POINT_OF_WATER,
               planet.min_temp <= planet.boil_point
            {
                print(" .. CALCULATE GASES")
                calculate_gases(sun: sun, planet: planet, planet_id: planet_id)
                print(" .. CALCULATE GASES (done)")
            }

            /*
             *    Next we assign a type to the planet.
             */

            if planet.surf_pressure < 1.0 {
                if !is_moon,
                   (planet.mass * SUN_MASS_IN_EARTH_MASSES) < ASTEROID_MASS_LIMIT
                {
                    planet.planet_type = .asteroids
                } else {
                    planet.planet_type = .rock
                }
            } else if planet.surf_pressure > 6000.0, planet.molec_weight <= 2.0 {
                // Retains Hydrogen
                planet.planet_type = .subsubgasgiant
                planet.gases = 0
                planet.atmosphere = []
            } else {
                // Atmospheres:
                if Int(planet.day) == Int(planet.orb_period * 24.0) || planet.resonant_period {
                    planet.planet_type = .face
                } else if planet.hydrosphere >= 0.95 {
                    planet.planet_type = .water // >95% water
                } else if planet.ice_cover >= 0.95 {
                    planet.planet_type = .ice // >95% ice
                } else if planet.hydrosphere > 0.05 {
                    planet.planet_type = .terrestrial // Terrestrial
                }
                // else <5% water
                else if planet.max_temp > planet.boil_point {
                    planet.planet_type = .venusian // Hot = Venusian
                } else if (planet.gas_mass / planet.mass) > 0.0001 {
                    // Accreted gas
                    planet.planet_type = .ice // But no Greenhouse
                    planet.ice_cover = 1.0 // or liquid water
                } // Make it an Ice World
                else if planet.surf_pressure <= 250.0 {
                    // Thin air = Martian
                    planet.planet_type = .martian
                } else if planet.surf_temp < FREEZING_POINT_OF_WATER {
                    planet.planet_type = .ice
                } else {
                    planet.planet_type = .unknown
                    //                        if (flag_verbose & 0x0001)
                    //                            fprintf (stderr, "%12s\tp=%4.2Lf\tm=%4.2Lf\tg=%4.2Lf\tt=%+.1Lf\t%s\t Unknown %s\n",
                    //                                            type_string (planet->type),
                    //                                            planet->surf_pressure,
                    //                                            planet->mass * SUN_MASS_IN_EARTH_MASSES,
                    //                                            planet->surf_grav,
                    //                                            planet->surf_temp  - EARTH_AVERAGE_KELVIN,
                    //                                            planet_id,
                    //                                            ((int)planet->day == (int)(planet->orb_period * 24.0) ||
                    //                                             (planet->resonant_period)) ? "(1-Face)" : ""
                    //                                     );
                }
                print(" .. .. .. -> planet type determined: \(planet.planet_type)")
            }
        }

        if do_moons, !is_moon {
            print(" .. .. .. PROCESSING MOONS")
            if planet.first_moon != nil {
                var n: Int = 0
                var ptr: Planet? = planet.first_moon

                while ptr != nil {
                    if ptr!.mass * SUN_MASS_IN_EARTH_MASSES > 0.000001 {
                        var moon_id: String
                        var roche_limit: Double = 0.0
                        var hill_sphere: Double = 0.0
                        ptr!.a = planet.a
                        ptr!.e = planet.e
                        n += 1
                        moon_id = "\(planet_id).\(n)"
                        generate_planet(planet: ptr!, planet_no: n, sun: sun, random_tilt: random_tilt, planet_id: moon_id, do_gases: do_gases, do_moons: do_moons, is_moon: true, prng: prng, counts: &counts)
                        // Adjusts ptr->density

                        roche_limit = 2.44 * planet.radius * pow(planet.density / ptr!.density, 1.0 / 3.0)
                        hill_sphere = planet.a * KM_PER_AU * pow(planet.mass / (3.0 * sun.mass), 1.0 / 3.0)

                        if (roche_limit * 3.0) < hill_sphere {
                            ptr!.a = prng.random_number(in: (roche_limit * 1.5) ... (hill_sphere / 2.0)) / KM_PER_AU
                            ptr!.e = prng.random_eccentricity()
                        } else {
                            ptr!.a = 0
                            ptr!.e = 0
                        }
                        //                            if (flag_verbose & 0x40000)
                        //                            {
                        //                                fprintf (stderr,
                        //                                            "   Roche limit: R = %4.2Lg, rM = %4.2Lg, rm = %4.2Lg -> %.0Lf km\n"
                        //                                            "   Hill Sphere: a = %4.2Lg, m = %4.2Lg, M = %4.2Lg -> %.0Lf km\n"
                        //                                            "%s Moon orbit: a = %.0Lf km, e = %.0Lg\n",
                        //                                            planet->radius, planet->density, ptr->density,
                        //                                            roche_limit,
                        //                                            planet->a * KM_PER_AU, planet->mass * SOLAR_MASS_IN_KILOGRAMS, sun->mass * SOLAR_MASS_IN_KILOGRAMS,
                        //                                            hill_sphere,
                        //                                            moon_id,
                        //                                            ptr->moon_a * KM_PER_AU, ptr->moon_e
                        //                                        );
                        //                            }
                        //
                        //                            if (flag_verbose & 0x1000)
                        //                            {
                        //                                fprintf (stderr, "  %s: (%7.2LfEM) %d %4.2LgEM\n",
                        //                                    planet_id,
                        //                                    planet->mass * SUN_MASS_IN_EARTH_MASSES,
                        //                                    n,
                        //                                    ptr->mass * SUN_MASS_IN_EARTH_MASSES);
                        //                            }
                    }
                    ptr = ptr?.next_planet
                }
            }
        }
    }
}

// NOTE(heckj): implementation/porting notes: mostly just collects counts and updates an overall count of what's been generated
// for the cases of searching for generations that result in habitable planets.
func check_planet(planet: Planet, planet_id: String, is_moon: Bool, counts: inout InterestingCounts) {
    if let a_count = counts.type_counts[planet.planet_type] {
        counts.type_counts[planet.planet_type] = a_count + 1
    } else {
        counts.type_counts[planet.planet_type] = 1
    }

    var min_breathable_terrestrial_g = 1000.0
    var min_breathable_g = 1000.0
    var max_breathable_terrestrial_g = 0.0
    var max_breathable_g = 0.0
    var min_breathable_temp = 1000.0
    var max_breathable_temp = 0.0
    var min_breathable_p = 100_000.0
    var max_breathable_p = 0.0
    var min_breathable_terrestrial_l = 1000.0
    var min_breathable_l = 1000.0
    var max_breathable_terrestrial_l = 0.0
    var max_breathable_l = 0.0
    var max_moon_mass = 0.0

    /* Check for and list planets with breathable atmospheres */

    let breathe: Breathability = breathability(planet: planet)
    var habitable = false

    if (breathe == .breathable) &&
        (!planet.resonant_period) && // Option needed?
        (Int(planet.day) != Int(planet.orb_period * 24.0))
    {
        var list_it = false
        let illumination = pow2(1.0 / planet.a) * planet.sun!.luminosity

        counts.habitable += 1
        habitable = true

        if min_breathable_temp > planet.surf_temp {
            min_breathable_temp = planet.surf_temp
            list_it = true
        }

        if max_breathable_temp < planet.surf_temp {
            max_breathable_temp = planet.surf_temp
            list_it = true
        }

        if min_breathable_g > planet.surf_grav {
            min_breathable_g = planet.surf_grav
            list_it = true
        }

        if max_breathable_g < planet.surf_grav {
            max_breathable_g = planet.surf_grav
            list_it = true
        }

        if min_breathable_l > illumination {
            min_breathable_l = illumination
            list_it = true
        }

        if max_breathable_l < illumination {
            max_breathable_l = illumination
            list_it = true
        }

        if planet.planet_type == .terrestrial {
            if min_breathable_terrestrial_g > planet.surf_grav {
                min_breathable_terrestrial_g = planet.surf_grav
                list_it = true
            }

            if max_breathable_terrestrial_g < planet.surf_grav {
                max_breathable_terrestrial_g = planet.surf_grav
                list_it = true
            }

            if min_breathable_terrestrial_l > illumination {
                min_breathable_terrestrial_l = illumination
                list_it = true
            }

            if max_breathable_terrestrial_l < illumination {
                max_breathable_terrestrial_l = illumination
                list_it = true
            }
        }

        if min_breathable_p > planet.surf_pressure {
            min_breathable_p = planet.surf_pressure
            list_it = true
        }

        if max_breathable_p < planet.surf_pressure {
            max_breathable_p = planet.surf_pressure
            list_it = true
        }

        if list_it {
            print("\(planet.planet_type)\tp=\(planet.surf_pressure.formatted(FPStyle))\tm=\((planet.mass * SUN_MASS_IN_EARTH_MASSES).formatted(FPStyle))\tg=\(planet.surf_grav.formatted(FPStyle))\tt=\((planet.surf_temp - EARTH_AVERAGE_KELVIN).formatted(FPStyle))\tl=\(illumination.formatted(FPStyle))\t\(planet_id)")
        }
    }

    if is_moon && max_moon_mass < planet.mass {
        max_moon_mass = planet.mass
        print("moon \(planet.planet_type)\tp=\(planet.surf_pressure.formatted(FPStyle))\tm=\((planet.mass * SUN_MASS_IN_EARTH_MASSES).formatted(FPStyle))\tg=\(planet.surf_grav.formatted(FPStyle))\tt=\((planet.surf_temp - EARTH_AVERAGE_KELVIN).formatted(FPStyle))\t\(planet_id)")
    }

    if planet.dust_mass * SUN_MASS_IN_EARTH_MASSES >= 0.0006,
       planet.gas_mass * SUN_MASS_IN_EARTH_MASSES >= 0.0006,
       planet.planet_type != .gasgiant,
       planet.planet_type != .subgasgiant
    {
        let core_size = Int((50.0 * planet.core_radius) / planet.radius)

        if core_size <= 49 {
            print("\(planet.planet_type)\tp=\(planet.core_radius.formatted(FPStyle))\tr=\(planet.radius.formatted(FPStyle))\tm=\((planet.mass * SUN_MASS_IN_EARTH_MASSES).formatted(FPStyle))\t\(planet_id)\t\(50 - core_size)")
        }
    }

    let rel_temp = (planet.surf_temp - FREEZING_POINT_OF_WATER) - EARTH_AVERAGE_CELSIUS
    let seas = (planet.hydrosphere * 100.0)
    let clouds = (planet.cloud_cover * 100.0)
    let pressure = (planet.surf_pressure / EARTH_SURF_PRES_IN_MILLIBARS)
    let ice = (planet.ice_cover * 100.0)
    let gravity = planet.surf_grav

    if gravity >= 0.8,
       gravity <= 1.2,
       rel_temp >= -2.0,
       rel_temp <= 3.0,
       ice <= 10.0,
       pressure >= 0.5,
       pressure <= 2.0,
       clouds >= 40.0,
       clouds <= 80.0,
       seas >= 50.0,
       seas <= 80.0,
       planet.planet_type != .water,
       breathe == .breathable
    {
        counts.earthlike += 1
        print("EARTHLIKE: \(planet.planet_type)\tp=\(planet.surf_pressure.formatted(FPStyle))\tm=\((planet.mass * SUN_MASS_IN_EARTH_MASSES).formatted(FPStyle))\tg=\(planet.surf_grav.formatted(FPStyle))\tt=\((planet.surf_temp - EARTH_AVERAGE_KELVIN).formatted(FPStyle))\t\(planet_id)")

    } else if breathe == .breathable,
              gravity > 1.3,
              habitable,
              (rel_temp < -2.0) || (ice > 10.0)
    {
        print("Sphinx-like: \(planet.planet_type)\tp=\(planet.surf_pressure.formatted(FPStyle))\tm=\((planet.mass * SUN_MASS_IN_EARTH_MASSES).formatted(FPStyle))\tg=\(planet.surf_grav.formatted(FPStyle))\tt=\((planet.surf_temp - EARTH_AVERAGE_KELVIN).formatted(FPStyle))\t\(planet_id)")
    }
}

func generate_planet_details(sun: Sun,
                             planets: [Planet],
                             random_tilt: Bool,
                             sys_no _: Int,
                             system_name: String,
                             do_gases: Bool,
                             do_moons: Bool,
                             prng: RNGWrapper<Xoshiro>,
                             counts: inout InterestingCounts)
{
    var planet_no = 1
    print(" .. GENERATE PLANET DETAILS")
    for planet in planets {
        let planet_id = "\(system_name) \(planet_no)"
        //                "%s (-s%ld -%c%d) %d",
        //                system_name, flag_seed, flag_char, sys_no, planet_no);
        print(" .. SYNTHESIZING DETAILS FOR \(planet_id) (\(planet.a.formatted(FPStyle)) AU)")
        generate_planet(planet: planet,
                        planet_no: planet_no,
                        sun: sun,
                        random_tilt: random_tilt,
                        planet_id: planet_id,
                        do_gases: do_gases,
                        do_moons: do_moons,
                        is_moon: false,
                        prng: prng,
                        counts: &counts)

        /*
         *    Now we're ready to test for habitable planets,
         *    so we can count and log them and such
         */

        print(" .. CHECK PLANET FOR \(planet_id) (\(planet.a.formatted(FPStyle)) AU)")
        check_planet(planet: planet, planet_id: planet_id, is_moon: false, counts: &counts)

        var moons = 1
        for moon in planet.moons {
            let moon_id = "\(planet_id).\(moons)"
            print(" .. CHECK PLANET FOR \(moon_id) (\(moon.a.formatted(FPStyle)) AU)")
            check_planet(planet: moon, planet_id: moon_id, is_moon: true, counts: &counts)
            moons += 1
        }
        planet_no += 1
    }
}

public struct FunctionFlags {
    public enum OutputFormat {
        case html
        case text
        case csv
        case svg
    }

    public enum GraphicFormat {
        case gif
        case png
        case jpeg
        case svg
    }

    public var do_catalog: Bool
    public var do_moons: Bool
    public var do_gases: Bool
    public var use_solar_system: Bool
    public var reuse_solar_system: Bool
    public var use_known_planets: Bool
    public var dont_generate: Bool
    // reporting results flags
    public var only_habitable: Bool
    public var only_multi_habitable: Bool
    public var only_jovian_habitable: Bool
    public var only_earthlike: Bool

    public var output_path: String
    public var filename_argument: String
    public var output_format: OutputFormat
    public var graphic_format: GraphicFormat

    public var system_name: String
    // public var sys_no_argument: Int
    public var mass_argument: Double // number of solar masses to use as a basis for RNG generation
    public var seed_argument: UInt64 // RNG seed
    public var count_argument: Int // number of star systems to generate
    public var catalog_argument: Catalog? // the star catalog to use for stellar luminosity and mass
    public var ratio_argument: Double // multiplier to the dust density coefficient

    public init(do_catalog: Bool, do_moons: Bool, do_gases: Bool, use_solar_system: Bool, reuse_solar_system: Bool, use_known_planets: Bool, dont_generate: Bool, only_habitable: Bool, only_multi_habitable: Bool, only_jovian_habitable: Bool, only_earthlike: Bool, output_path: String, filename_argument: String, output_format: OutputFormat, graphic_format: GraphicFormat, system_name: String, mass_argument: Double, seed_argument: UInt64, count_argument: Int = 1, catalog_argument: Catalog? = nil, ratio_argument: Double = 1.0) {
        self.do_catalog = do_catalog
        self.do_moons = do_moons
        self.do_gases = do_gases
        self.use_solar_system = use_solar_system
        self.reuse_solar_system = reuse_solar_system
        self.use_known_planets = use_known_planets
        self.dont_generate = dont_generate
        self.only_habitable = only_habitable
        self.only_multi_habitable = only_multi_habitable
        self.only_jovian_habitable = only_jovian_habitable
        self.only_earthlike = only_earthlike
        self.output_path = output_path
        self.filename_argument = filename_argument
        self.output_format = output_format
        self.graphic_format = graphic_format
        self.system_name = system_name
        self.mass_argument = mass_argument
        self.seed_argument = seed_argument
        self.count_argument = count_argument
        self.catalog_argument = catalog_argument
        self.ratio_argument = ratio_argument
    }
}

struct InterestingCounts {
    var earthlike = 0
    var habitable = 0
    var habitable_jovians = 0
    var type_counts: [PlanetType: Int] = [:]
}

//  main entrance point to invoking generation or exploration
public func stargen(flags: FunctionFlags) {
    let sun = Sun(luminosity: 0, mass: 0, life: 0, age: 0, r_ecosphere: 0, name: "")
    let min_mass = 0.4
    let inc_mass = 0.05
    let max_mass = 2.35
    var system_count = 1

    let prng = RNGWrapper(Xoshiro(seed: flags.seed_argument))

        sun.mass = flags.mass_argument
        system_count = flags.count_argument
        var seed_planets: Planet?
        let dust_density_coeff: Double
        var use_seed_system = false

        if flags.ratio_argument > 0.0 {
            dust_density_coeff = DUST_DENSITY_COEFF * flags.ratio_argument
        } else {
            dust_density_coeff = DUST_DENSITY_COEFF
        }

        if flags.reuse_solar_system {
            system_count = 1 + Int((max_mass - min_mass) / inc_mass)
            sun.luminosity = 1
            sun.mass = 1
            sun.life = 1e10
            sun.age = 5e9
            sun.r_ecosphere = 1
        } else if flags.do_catalog {
            system_count = flags.catalog_argument?.count ?? 0
        }

        for iteration in 1 ... system_count {
            let system_name: String
//            let designation: String
            let outer_limit: Double

            if flags.do_catalog {
                guard let catalog = flags.catalog_argument else {
                    break
                }
                sun.mass = catalog.stars[iteration].mass
                sun.luminosity = catalog.stars[iteration].luminosity
                sun.name = catalog.stars[iteration].name

                system_name = catalog.stars[iteration].desig
//                designation = catalog.stars[iteration].desig

                if catalog.stars[iteration].m2 > 0.001 {
                    /*
                     *    The following is Holman & Wiegert's equation 1 from
                     *    Long-Term Stability of Planets in Binary Systems
                     *    The Astronomical Journal, 117:621-628, Jan 1999
                     */
                    let m1 = sun.mass
                    let m2 = catalog.stars[iteration].m2
                    let mu = m2 / (m1 + m2)
                    let e = catalog.stars[iteration].e
                    let a = catalog.stars[iteration].a
                    outer_limit = (0.464 + (-0.380 * mu) + (-0.631 * e) +
                        (0.586 * mu * e) + (0.150 * pow2(e)) +
                        (-0.198 * mu * pow2(e))) * a
                } else {
                    outer_limit = 0.0
                }

            } else if flags.reuse_solar_system {
                system_name = "Earth-M\(earth.mass * SUN_MASS_IN_EARTH_MASSES)"
//                designation = system_name
                outer_limit = 0.0
            } else if !flags.system_name.isEmpty {
                system_name = flags.system_name
//                designation = system_name
                outer_limit = 0.0
            } else {
                system_name = "SolarSystem \(flags.seed_argument)-\(sun.mass)"
//                designation = system_name
                outer_limit = 0.0
            }
            sun.name = system_name

            // Sun constraints

            if sun.mass < 0.2 || sun.mass > 1.5 {
                sun.mass = prng.random_number(in: 0.7 ... 1.4)
            }
            if sun.luminosity == 0 {
                sun.luminosity = luminosity(mass_ratio: sun.mass)
            }
            sun.r_ecosphere = sqrt(sun.luminosity)
            sun.life = 1.0e10 * (sun.mass / sun.luminosity)

            var counts = InterestingCounts()
            if flags.reuse_solar_system || flags.use_solar_system {
                seed_planets = mercury
                use_seed_system = true
            }

            let solarsystem = generate_stellar_system(sun: sun, use_seed_system: use_seed_system, seed_system: seed_planets, sys_no: iteration, system_name: system_name, outer_planet_limit: outer_limit, dust_density_coeff: dust_density_coeff, do_gases: flags.do_gases, do_moons: flags.do_moons, prng: prng, counts: &counts)
            solarsystem.text_describe_system()
        }
}
