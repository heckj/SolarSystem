//
//  Accrete.swift
//

import Foundation

/*------------------------------------------------------------------------*/
/*                             BIBLIOGRAPHY                               */
/*    Dole, Stephen H.  "Formation of Planetary Systems by Aggregation:   */
/*        a Computer Simulation"    October 1969,  Rand Corporation Paper */
/*        P-4226.                                                         */
/*------------------------------------------------------------------------*/

// NOTE(heckj): throughout this ported code, 'a' represents the AU distance from the core of the disk,
// 'mass' is the accumulated mass at that AU, and 'e' represents the eccentricity of the mass being calculated.

// Original C functions in header:
//
//void set_initial_conditions(long double, long double );
//long double stellar_dust_limit(long double);
//long double nearest_planet(long double);
//long double farthest_planet(long double);
//long double inner_effect_limit(long double, long double, long double );
//long double outer_effect_limit(long double, long double, long double );
//int dust_available(long double, long double );
//void update_dust_lanes(long double, long double, long double, long double, long double, long double );
//long double collect_dust(long double, long double *, long double *, long double, long double, long double, dust_pointer);
//long double critical_limit(long double, long double, long double );
//void accrete_dust(long double *, long double *, long double *, long double, long double, long double, long double, long double );
//void coalesce_planetesimals(long double, long double, long double, long double, long double, long double, long double, long double, long double, int );
//planet_pointer dist_planetary_masses(long double, long double, long double, long double, long double, long double, planet_pointer, int);
//void free_dust (dust_pointer);
//void free_planet (planet_pointer);
//void free_atmosphere(planet_pointer);
//void free_generations();

//#include    <stdio.h>
//#include    <stdlib.h>
//#include    <math.h>
//#include    "const.h"
//#include    "structs.h"
//#include     "accrete.h"
//#include     "stargen.h"

// Accrete was a variety of free functions that manipulated a set of global variables.
// I'm encapsulating those into a struct, which seems a bit awkward, but better than slapping
// around global variables and pointers to me.
struct AccretionDisk {
    var dust_left: Bool = true
    var r_inner: Double = 0 //? unsure about initial value here - set within code, but not at initial call sites
    var r_outer: Double = 0 //? unsure about initial value here - set within code, but not at initial call sites
    var reduced_mass: Double = 0 //? unsure about initial value here - set within code?, but not at initial call sites
    var dust_density: Double = 0 // initial value set within `dist_planetary_masses` and used while accumulating dust
    // ^^ could probably be a constant, as I don't think it's changed after it's set, if we collapse initial settings
    // and 'dist_planetary_masses' into some sort of initializer
    var cloud_eccentricity: Double = 0.2
    
    // NOTE(heckj): create a proper initializer to set this value, aligned somewhere between the initializer and calling
    // dist_planetary_masses, which is the primarily call site into this...
    var prng = RNGWrapper(Xoshiro(seed: 541))
    
    // In the earlier versions of accrete, these were implemented as a linked-list of pointers
    // in C, so I've transferred that structure to reference types (final class), originally
    // defined in structs.h to the types included within Types.swift.
    var dust_head: Dust? = nil
    var planet_head: Planet? = nil
    var generation_head: Generation? = nil
    
    ///* Now for some variables global to the accretion process:        */
    //int             dust_left;
    //long double        r_inner;
    //long double        r_outer;
    //long double        reduced_mass;
    //long double        dust_density;
    //long double        cloud_eccentricity;
    //dust_pointer    dust_head    = NULL;
    //planet_pointer    planet_head    = NULL;
    //gen_pointer        hist_head    = NULL;
    
    // called from dist_planetary_masses, which seems to be the primary entry point to these functions.
    // NOTE(heckj): This is completely setup related, so this might be WAY more sane in an
    // intializer during refactoring. For the initial porting, I'm trying to maintain
    // the same basic flow as the C code.
    mutating func set_initial_conditions(inner_limit_of_dust: Double, outer_limit_of_dust: Double) {
        
        // original code sets up history _before_ adding explicit dust cloud and details,
        // so I think the first iteration of history is intended to be empty...
        var hist = Generation(dust: nil, planet: nil, next: nil)
        hist.dust = dust_head
        hist.planet = planet_head
        generation_head = hist
        
        let dust = Dust(inner_edge: inner_limit_of_dust, outer_edge: outer_limit_of_dust, dust_present: true, gas_present: true, next_band: nil)
        dust_head = dust
        planet_head = nil
        cloud_eccentricity = 0.2
        dust_left = true
    }
    // called from dist_planetary_masses, which seems to be the primary entry point to these functions.
    //void set_initial_conditions(long double inner_limit_of_dust,
    //                            long double outer_limit_of_dust)
    //{
    //    gen_pointer hist;
    //    hist = (gen_pointer)malloc(sizeof(generation));
    //    hist->dusts = dust_head;
    //    hist->planets = planet_head;
    //    hist->next = hist_head;
    //    hist_head = hist;
    //
    //    dust_head = (dust *)malloc(sizeof(dust));
    //    planet_head = NULL;
    //    dust_head->next_band = NULL;
    //    dust_head->outer_edge = outer_limit_of_dust;
    //    dust_head->inner_edge = inner_limit_of_dust;
    //    dust_head->dust_present = TRUE;
    //    dust_head->gas_present = TRUE;
    //    dust_left = TRUE;
    //    cloud_eccentricity = 0.2;
    //}
    
    func stellar_dust_limit(stell_mass_ratio: Double) -> Double {
        200.0 * pow(stell_mass_ratio, (1.0/3.0))
    }
    //    long double stellar_dust_limit(long double stell_mass_ratio)
    //    {
    //        return(200.0 * pow(stell_mass_ratio,(1.0 / 3.0)));
    //    }
    
    // Uses the stellar mass ratio to determine the closest "interior" planet
    // that would likely form. The result is in "AU".
    func nearest_planet(stell_mass_ratio: Double) -> Double {
        0.3 * pow(stell_mass_ratio, (1.0/3.0))
    }
    //    long double nearest_planet(long double stell_mass_ratio)
    //    {
    //        return(0.3 * pow(stell_mass_ratio,(1.0 / 3.0)));
    //    }
    
    // Uses the stellar mass ratio to determine the farthest planet
    // that would likely form. The result is in "AU".
    func farthest_planet(stell_mass_ratio: Double) -> Double {
        50.0 * pow(stell_mass_ratio, (1.0/3.0))
    }
    //    long double farthest_planet(long double stell_mass_ratio)
    //    {
    //        return(50.0 * pow(stell_mass_ratio,(1.0 / 3.0)));
    //    }
    
    func inner_effect_limit(a: Double, e: Double, mass: Double) -> Double {
        // Given an orbital distance, mass, and eccentricity, calculates the
        // inner-most distance (in AU) that the orbit "sweeps up" dust.
        a * (1.0 - e) * (1.0 - mass) / (1.0 + cloud_eccentricity)
    }
    //    long double inner_effect_limit(long double a, long double e, long double mass)
    //    {
    //        return (a * (1.0 - e) * (1.0 - mass) / (1.0 + cloud_eccentricity));
    //    }
    
    func outer_effect_limit(a: Double, e: Double, mass: Double) -> Double {
        // Given an orbital distance, mass, and eccentricity, calculates the
        // outer-most distance (in AU) that the orbit "sweeps up" dust.
        a * (1.0 + e) * (1.0 + mass) / (1.0 - cloud_eccentricity)
    }
    //    long double outer_effect_limit(long double a, long double e, long double mass)
    //    {
    //        return (a * (1.0 + e) * (1.0 + mass) / (1.0 - cloud_eccentricity));
    //    }
    
    func dust_available(inside_range: Double, outside_range: Double) -> Bool {
        // This iterates through the linked-list of all the "bands of dust"
        // (represented by Dust.next -> Dust.next -> nil) to search for the areas
        // that match the region defined between inside_range and outside_range.
        // If any dust band exists, and has dust, within those ranges, this will return
        // true.
        
        // NOTE(heckj): refactoring note - I'll have to see how these dust bands are
        // created and updated, but it might make a lot more sense to represent them in
        // an array, maybe using range values to determine what elements of the array
        // are relevant, and using a functional map or filter to return a boolean.
        // At a minimum, doing something that allows us to use Sequence and iterator would
        // make this a much more understandable bit of code.
        var current_dust_band = dust_head
        var dust_here: Bool = false
        
        while current_dust_band != nil && current_dust_band!.outer_edge < inside_range {
            current_dust_band = current_dust_band?.next_band
        }
        
        if let definitely_current_dust_band = current_dust_band {
            dust_here = definitely_current_dust_band.dust_present
        } else {
            dust_here = false
        }
        
        while current_dust_band != nil && current_dust_band!.inner_edge < outside_range {
            dust_here = dust_here || current_dust_band!.dust_present
            current_dust_band = current_dust_band?.next_band
        }
        
        return dust_here
    }
    //    int dust_available(long double inside_range, long double outside_range)
    //    {
    //        dust_pointer current_dust_band;
    //        int dust_here;
    //
    //        current_dust_band = dust_head;
    
    //        while ((current_dust_band != NULL)
    //            && (current_dust_band->outer_edge < inside_range))
    //            current_dust_band = current_dust_band->next_band;
    
    //        if (current_dust_band == NULL)
    //            dust_here = FALSE;
    //        else dust_here = current_dust_band->dust_present;
    
    //        while ((current_dust_band != NULL)
    //            && (current_dust_band->inner_edge < outside_range)) {
    //                dust_here = dust_here || current_dust_band->dust_present;
    //                current_dust_band = current_dust_band->next_band;
    //            }
    //        return(dust_here);
    //    }
    
    
    // called from `accrete_dust`
    mutating func update_dust_lanes(min: Double, max: Double, mass: Double, crit_mass: Double, body_inner_bound: Double, body_outer_bound: Double) {
        // With the collected mass accumulated earlier in sweeping through the dust lanes, we need to deplete the dust
        // stored within our data structures (the linked list of dust lanes) between the inner_bound and outer_bound
        // distances in AU.
        
        // This does so by updating the linked list, breaking up the dust into separate lanes and they get "swept through"
        // and accumulated into masses.
        
        var gas: Bool = (mass <= crit_mass)
        var node1: Dust? = nil
        var node2: Dust? = nil
        var node3: Dust? = nil
        dust_left = false
        node1 = dust_head
        while (node1 != nil) {
            if let definitely_node1 = node1 {
                if definitely_node1.inner_edge < min && definitely_node1.outer_edge > max {
                    //  node1.inner ------------------------- node1.outer
                    //                     min ---- max
                    // make a new dust lane and insert it into the linked list, marking
                    // that dust lane as consumed, squeezing down the first lane to 'min'
                    // and creating a new dust lane to cover 'max' to the outer edge
                    node2 = Dust(inner_edge: min, outer_edge: max, dust_present: false, gas_present: definitely_node1.gas_present ? gas : false, next_band: nil)
                    node3 = Dust(inner_edge: max, outer_edge: definitely_node1.outer_edge, dust_present: definitely_node1.dust_present, gas_present: definitely_node1.gas_present, next_band: definitely_node1.next_band)
                    definitely_node1.next_band = node2
                    node2?.next_band = node3
                    node1?.outer_edge = min
                    node1 = node3?.next_band
                } else if definitely_node1.inner_edge < max && definitely_node1.outer_edge > max {
                    //        node1.inner ------ node1.outer
                    //   min -------------- max
                    
                    node2 = Dust(inner_edge: max, outer_edge: definitely_node1.outer_edge, dust_present: definitely_node1.dust_present, gas_present: definitely_node1.gas_present, next_band: definitely_node1.next_band)
                    definitely_node1.next_band = node2
                    definitely_node1.outer_edge = max
                    definitely_node1.gas_present = definitely_node1.gas_present ? gas : false
                    definitely_node1.dust_present = false
                    node1 = node2?.next_band
                } else if definitely_node1.inner_edge < min && definitely_node1.outer_edge > min {
                    //        node1.inner -------- node1.outer
                    //                       min -------------- max
                    node2 = Dust(inner_edge: min, outer_edge: definitely_node1.outer_edge, dust_present: false, gas_present: definitely_node1.gas_present ? gas : false, next_band: definitely_node1.next_band)
                    definitely_node1.next_band = node2
                    definitely_node1.outer_edge = min
                    node1 = node2?.next_band
                } else if definitely_node1.inner_edge >= min && definitely_node1.outer_edge <= max {
                    //       node1.inner --- node1.outer
                    //  min ------------------------------ max
                    definitely_node1.dust_present = false
                    definitely_node1.gas_present = definitely_node1.gas_present ? gas : false
                    node1 = node1?.next_band
                } else if definitely_node1.outer_edge < min || definitely_node1.inner_edge > max {
                    //       node1.inner --- node1.outer
                    //                                     min ----- max
                    // OR
                    //                       node1.inner --- node1.outer
                    //       min ----- max
                    node1 = node1?.next_band
                }
            }
        }
        // walk through the tree after the updates, and collapse bands that are equivalent
        // together
        node1 = dust_head
        while (node1 != nil) {
            if let definitely_node1 = node1 {
                if definitely_node1.dust_present && definitely_node1.outer_edge >= body_inner_bound && definitely_node1.inner_edge <= body_outer_bound {
                    dust_left = true
                }
                node2 = definitely_node1.next_band
                if let definitely_node2 = node2 {
                    if definitely_node1.dust_present == definitely_node2.dust_present && definitely_node1.gas_present == definitely_node2.gas_present {
                        definitely_node1.outer_edge = definitely_node2.outer_edge
                        definitely_node1.next_band = definitely_node2.next_band
                        // free(node2)
                    }
                }
            }
            node1 = node1?.next_band
        }
    }
    
    
    // called from `accrete_dust`
    //void update_dust_lanes(long double min, long double max, long double mass,
    //                       long double crit_mass, long double body_inner_bound,
    //                       long double body_outer_bound)
    //{
    //    int             gas;
    //    dust_pointer    node1;
    //    dust_pointer    node2;
    //    dust_pointer    node3;
    //
    //    dust_left = FALSE;
    //    if ((mass > crit_mass))
    //        gas = FALSE;
    //    else
    //        gas = TRUE;
    //    node1 = dust_head;
    //    while ((node1 != NULL))
    //    {
    //        if (((node1->inner_edge < min) && (node1->outer_edge > max)))
    //        {
    //            node2 = (dust *)malloc(sizeof(dust));
    //            node2->inner_edge = min;
    //            node2->outer_edge = max;
    //            if ((node1->gas_present == TRUE))
    //                node2->gas_present = gas;
    //            else
    //                node2->gas_present = FALSE;
    //            node2->dust_present = FALSE;
    //            node3 = (dust *)malloc(sizeof(dust));
    //            node3->inner_edge = max;
    //            node3->outer_edge = node1->outer_edge;
    //            node3->gas_present = node1->gas_present;
    //            node3->dust_present = node1->dust_present;
    //            node3->next_band = node1->next_band;
    
    //            node1->next_band = node2;
    //            node2->next_band = node3;
    //            node1->outer_edge = min;
    //            node1 = node3->next_band;
    //        }
    //        else
    //            if (((node1->inner_edge < max) && (node1->outer_edge > max)))
    //            {
    //                node2 = (dust *)malloc(sizeof(dust));
    //                node2->next_band = node1->next_band;
    //                node2->dust_present = node1->dust_present;
    //                node2->gas_present = node1->gas_present;
    //                node2->outer_edge = node1->outer_edge;
    //                node2->inner_edge = max;
    //                node1->next_band = node2;
    //                node1->outer_edge = max;
    //                if ((node1->gas_present == TRUE))
    //                    node1->gas_present = gas;
    //                else
    //                    node1->gas_present = FALSE;
    //                node1->dust_present = FALSE;
    //                node1 = node2->next_band;
    //            }
    //            else
    //                if (((node1->inner_edge < min) && (node1->outer_edge > min)))
    //                {
    //                    node2 = (dust *)malloc(sizeof(dust));
    //                    node2->next_band = node1->next_band;
    //                    node2->dust_present = FALSE;
    //                    if ((node1->gas_present == TRUE))
    //                        node2->gas_present = gas;
    //                    else
    //                        node2->gas_present = FALSE;
    //                    node2->outer_edge = node1->outer_edge;
    //                    node2->inner_edge = min;
    //                    node1->next_band = node2;
    //                    node1->outer_edge = min;
    //                    node1 = node2->next_band;
    //                }
    //                else
    //                    if (((node1->inner_edge >= min) && (node1->outer_edge <= max)))
    //                    {
    //                        if ((node1->gas_present == TRUE))
    //                            node1->gas_present = gas;
    //                        node1->dust_present = FALSE;
    //                        node1 = node1->next_band;
    //                    }
    //                    else
    //                        if (((node1->outer_edge < min) || (node1->inner_edge > max)))
    //                            node1 = node1->next_band;
    //    }
    //    node1 = dust_head;
    //    while ((node1 != NULL))
    //    {
    //        if (((node1->dust_present)
    //            && (((node1->outer_edge >= body_inner_bound)
    //                && (node1->inner_edge <= body_outer_bound)))))
    //            dust_left = TRUE;
    //        node2 = node1->next_band;
    //        if ((node2 != NULL))
    //        {
    //            if (((node1->dust_present == node2->dust_present)
    //                && (node1->gas_present == node2->gas_present)))
    //            {
    //                node1->outer_edge = node2->outer_edge;
    //                node1->next_band = node2->next_band;
    //                free(node2);
    //            }
    //        }
    //        node1 = node1->next_band;
    //    }
    //}
    
    
    // called from `accrete_dust`
    func collect_dust(last_mass: Double, new_dust: inout Double, new_gas: inout Double, a: Double, e: Double, crit_mass: Double, dust_band: Dust?) -> Double {
        // recursively sweeps through the existing dust lanes, recursively walking down the linked list of dust lanes
        // and collecting all the dust defined by the initial orbit position (a: in AU) with an eccentricity (e), the boundary
        // markers of which are calculated by inner_effect_limit() and outer_effect_limit().
        
        var mass_density: Double
        var    temp1: Double
        var    temp2: Double
        
        var    temp_density: Double
        var    bandwidth: Double
        var    width: Double
        var    volume: Double
        var    gas_density: Double = 0.0
        var    new_mass: Double
        var    next_mass: Double
        var    next_dust: Double = 0
        var    next_gas: Double = 0
        
        var temp = last_mass / (1.0 + last_mass)
        let reduced_mass = pow(temp,(1.0 / 4.0))
        
        var r_inner = inner_effect_limit(a: a, e: e, mass: reduced_mass);
        let r_outer = outer_effect_limit(a: a, e: e, mass: reduced_mass);
        if r_inner < 0.0 {
            r_inner = 0.0
        }
        guard let dust_band = dust_band else {
            return 0.0
        }
        
        if (dust_band.dust_present == false) {
            temp_density = 0.0
        } else {
            temp_density = dust_density
        }
        
        if last_mass < crit_mass || dust_band.gas_present == false {
            mass_density = temp_density
        } else {
            mass_density = K * temp_density / (1.0 + sqrt(crit_mass / last_mass) * (K - 1.0))
            gas_density = mass_density - temp_density
        }
        
        if dust_band.outer_edge <= r_inner || dust_band.inner_edge >= r_outer {
            // nothing in this band, walk down and process the next...
            return collect_dust(last_mass: last_mass, new_dust: &new_dust, new_gas: &new_gas, a: a, e: e, crit_mass: crit_mass, dust_band: dust_band.next_band)
        } else {
            bandwidth = (r_outer - r_inner)
            
            temp1 = r_outer - dust_band.outer_edge
            if (temp1 < 0.0) {
                temp1 = 0.0
            }
            width = bandwidth - temp1
            
            temp2 = dust_band.inner_edge - r_inner
            if (temp2 < 0.0) {
                temp2 = 0.0
            }
            width = width - temp2
            
            temp = 4.0 * Double.pi * pow(a,2.0) * reduced_mass * (1.0 - e * (temp1 - temp2) / bandwidth)
            volume = temp * width
            
            new_mass  = volume * mass_density
            new_gas  = volume * gas_density
            new_dust = new_mass - new_gas
            
            next_mass = collect_dust(last_mass: last_mass, new_dust: &next_dust, new_gas: &next_gas, a: a, e: e, crit_mass: crit_mass, dust_band: dust_band.next_band)
            
            new_gas  +=  next_gas
            new_dust +=  next_dust
            
            return(new_mass + next_mass)
        }
    }
    
    //// recursive function called originally from `accrete_dust`
    //long double collect_dust(long double last_mass, long double *new_dust,
    //                         long double *new_gas,
    //                         long double a, long double e,
    //                         long double crit_mass, dust_pointer dust_band)
    //{
    //    long double    mass_density;
    //    long double    temp1;
    //    long double    temp2;
    //    long double    temp;
    //    long double    temp_density;
    //    long double    bandwidth;
    //    long double    width;
    //    long double    volume;
    //    long double    gas_density = 0.0;
    //    long double    new_mass;
    //    long double    next_mass;
    //    long double    next_dust = 0;
    //    long double    next_gas = 0;
    //
    //
    //    temp = last_mass / (1.0 + last_mass);
    //    reduced_mass = pow(temp,(1.0 / 4.0));
    //    r_inner = inner_effect_limit(a, e, reduced_mass);
    //    r_outer = outer_effect_limit(a, e, reduced_mass);
    //
    //    if ((r_inner < 0.0))
    //        r_inner = 0.0;
    //
    //    if ((dust_band == NULL))
    //        return(0.0);
    //    else
    //    {
    //        if ((dust_band->dust_present == FALSE))
    //            temp_density = 0.0;
    //        else
    //            temp_density = dust_density;
    //
    //        if (((last_mass < crit_mass) || (dust_band->gas_present == FALSE)))
    //            mass_density = temp_density;
    //        else
    //        {
    //            mass_density = K * temp_density / (1.0 + sqrt(crit_mass / last_mass)
    //                                        * (K - 1.0));
    //            gas_density = mass_density - temp_density;
    //        }
    //
    //        if (((dust_band->outer_edge <= r_inner)
    //          || (dust_band->inner_edge >= r_outer)))
    //        {
    //            return(collect_dust(last_mass, new_dust, new_gas,
    //                                a,e,crit_mass, dust_band->next_band));
    //        }
    //        else
    //        {
    //            bandwidth = (r_outer - r_inner);
    //
    //            temp1 = r_outer - dust_band->outer_edge;
    //            if (temp1 < 0.0)
    //                temp1 = 0.0;
    //            width = bandwidth - temp1;
    //
    //            temp2 = dust_band->inner_edge - r_inner;
    //            if (temp2 < 0.0)
    //                temp2 = 0.0;
    //            width = width - temp2;
    //
    //            temp = 4.0 * PI * pow(a,2.0) * reduced_mass
    //                * (1.0 - e * (temp1 - temp2) / bandwidth);
    //            volume = temp * width;
    //
    //            new_mass  = volume * mass_density;
    //            *new_gas  = volume * gas_density;
    //            *new_dust = new_mass - *new_gas;
    //
    //            next_mass = collect_dust(last_mass, &next_dust, &next_gas,
    //                                     a,e,crit_mass, dust_band->next_band);
    //
    //            *new_gas  = *new_gas + next_gas;
    //            *new_dust = *new_dust + next_dust;
    //
    //            return(new_mass + next_mass);
    //        }
    //    }
    //}
    
    
    /*--------------------------------------------------------------------------*/
    /*     Orbital radius is in AU, eccentricity is unitless, and the stellar        */
    /*    luminosity ratio is with respect to the sun.  The value returned is the */
    /*    mass at which the planet begins to accrete gas as well as dust, and is    */
    /*    in units of solar masses.                                                */
    /*--------------------------------------------------------------------------*/
    
    func critical_limit(orb_radius: Double, eccentricity: Double, stell_luminosity_ratio: Double) -> Double {
        let perihelion_dist = (orb_radius - orb_radius * eccentricity)
        let temp = perihelion_dist * sqrt(stell_luminosity_ratio)
        return B * pow(temp, -0.75)
    }
    //long double critical_limit(long double orb_radius, long double eccentricity,
    //                           long double stell_luminosity_ratio)
    //{
    //    long double    temp;
    //    long double    perihelion_dist;
    //
    //    perihelion_dist = (orb_radius - orb_radius * eccentricity);
    //    temp = perihelion_dist * sqrt(stell_luminosity_ratio);
    //    return(B * pow(temp,-0.75));
    //}
    
    // called from `coalesce_planetesimals` and `dist_planetary_masses`
    mutating func accrete_dust(seed_mass: inout Double, new_dust: inout Double, new_gas: inout Double, a: Double, e: Double, crit_mass: Double, body_inner_bound: Double, body_outer_bound: Double) {
        // This function is effectively creating planetesimals, sweeping through dust lanes using 'collect dust' until
        // it stops growing at a notable rate. That's set by the comparison of "mass collected" and the last sweep through
        // being less thatn 0.0001 times larger than the previous sweep.
        var new_mass = seed_mass
        var temp_mass: Double
        repeat {
            temp_mass = new_mass
            new_mass = collect_dust(last_mass: new_mass, new_dust: &new_dust, new_gas: &new_gas, a: a, e: e, crit_mass: crit_mass, dust_band: dust_head)
            
        } while !((new_mass - temp_mass) < (0.0001 * temp_mass))
        
        seed_mass += new_mass
        update_dust_lanes(min: r_inner, max: r_outer, mass: seed_mass, crit_mass: crit_mass, body_inner_bound: body_inner_bound, body_outer_bound: body_outer_bound)
    }
    // called from `coalesce_planetesimals` and `dist_planetary_masses`
    //    void accrete_dust(long double *seed_mass, long double *new_dust, long double *new_gas,
    //                      long double a, long double e, long double crit_mass,
    //                      long double body_inner_bound, long double body_outer_bound)
    //    {
    //        long double    new_mass = (*seed_mass);
    //        long double    temp_mass;
    //
    //        do
    //        {
    //            temp_mass = new_mass;
    //            new_mass = collect_dust(new_mass, new_dust, new_gas,
    //                                    a,e,crit_mass, dust_head);
    //        }
    //        while (!(((new_mass - temp_mass) < (0.0001 * temp_mass))));
    //
    //        (*seed_mass) = (*seed_mass) + new_mass;
    //        update_dust_lanes(r_inner,r_outer,(*seed_mass),crit_mass,body_inner_bound,body_outer_bound);
    //    }
    
    // called from `dist_planetary_masses`
    // NOTE(heckj): Refactoring note - it looks like we might be able to make this a lot easier
    // by implementing a set of planets in an array and iterating through them, updating the ordering as needed. The
    // current code is straight from the original accrete, and follows the linked-list mechanisms that were used there.
    mutating func coalesce_planetesimals(a: Double, e: Double, mass: Double, crit_mass: Double, dust_mass: Double, gas_mass: Double, stell_luminosity_ratio: Double, body_inner_bound: Double, body_outer_bound: Double, do_moons: Bool) {
        
        var the_planet: Planet? = nil
        var next_planet: Planet? = nil
        var prev_planet: Planet? = nil
        var finished: Bool = false
        let dist1: Double
        let dist2: Double
        var temp: Double = 0
        // First we try to find an existing planet with an over-lapping orbit.
        
        the_planet = planet_head
        while let this_planet = the_planet {
            
            let diff = this_planet.a - a
            the_planet = the_planet?.next_planet
            if diff > 0.0 {
                dist1 = (a * (1.0 + e) * (1.0 + reduced_mass)) - a
                /* x aphelion     */
                reduced_mass = pow((this_planet.mass / (1.0 + this_planet.mass)),(1.0 / 4.0))
                dist2 = this_planet.a - (this_planet.a * (1.0 - this_planet.e) * (1.0 - reduced_mass))
            } else {
                dist1 = a - (a * (1.0 - e) * (1.0 - reduced_mass))
                /* x perihelion */
                reduced_mass = pow((this_planet.mass / (1.0 + this_planet.mass)),(1.0 / 4.0));
                dist2 = (this_planet.a * (1.0 + this_planet.e) * (1.0 + reduced_mass)) - this_planet.a;
            }
            // fabs(_) returns the absolute value of the provided floating point number
            if fabs(diff) <= fabs(dist1) || fabs(diff) <= fabs(dist2) {
                var new_dust: Double = 0
                var new_gas: Double = 0
                var new_a = (this_planet.mass + mass) / ((this_planet.mass / this_planet.a) + (mass / a))
                temp = this_planet.mass * sqrt(this_planet.a) * sqrt(1.0 - pow(this_planet.e, 2.0))
                temp += (mass * sqrt(a) * sqrt(sqrt(1.0 - pow(e, 2.0))))
                temp /= ((this_planet.mass + mass) * sqrt(new_a))
                temp = 1.0 - pow(temp,2.0)
                if temp < 0.0 || temp >= 1.0 {
                    temp = 0.0
                }
                let e = sqrt(temp)
                
                if do_moons {
                    var existing_mass: Double = 0.0
                    
                    // walk through the linked list of moons and get the summed mass of them
                    var current_moon: Planet? = this_planet.first_moon
                    while (current_moon != nil) {
                        existing_mass += current_moon?.mass ?? 0.0
                        current_moon = current_moon?.next_planet
                    }
                    if mass < crit_mass {
                        if mass * SUN_MASS_IN_EARTH_MASSES < 2.5 &&
                            mass * SUN_MASS_IN_EARTH_MASSES > 0.001 &&
                            existing_mass < this_planet.mass * 0.05 {
                            //NOTE(heckj): Refactoring note - all we're really caring about with Planets (planetisimals)
                            // at this point is mass (total, dust, gas) and distance (au), and to a lesser extent we
                            // set gas_giant = true when it exceeds a critical mass point. So the rest of the initializer
                            // bits are just excessive noise at this point, most of which can default to 0 to be calculated
                            // at a later stage.
                            let the_moon = Planet(
                                planet_no: 0,
                                a: a,
                                e: e,
                                axial_tilt: 0,
                                mass: mass,
                                gas_giant: false,
                                dust_mass: dust_mass, gas_mass: gas_mass,
                                moon_a: 0, moon_e: 0,
                                core_radius: 0, radius: 0,
                                orbit_zone: 0,
                                density: 0,
                                orb_period: 0,
                                day: 0,
                                resonant_period: false,
                                esc_velocity: 0,
                                surf_accel: 0,
                                surf_grav: 0,
                                rms_velocity: 0,
                                molec_weight: 0,
                                volatile_gas_inventory: 0,
                                surf_pressure: 0,
                                greenhouse_effect: false,
                                boil_point: 0,
                                albedo: 0, exospheric_temp: 0, estimated_temp: 0, estimated_terr_temp: 0,
                                surf_temp: 0,
                                greenhs_rise: 0,
                                high_temp: 0, low_temp: 0, max_temp: 0, min_temp: 0,
                                hydrosphere: 0, cloud_cover: 0, ice_cover: 0, // default to 0
                                sun: nil, gases: 0, atmosphere: nil, planet_type: .unknown,
                                minor_moons: 0, first_moon: nil, next_planet: nil)
                            if (the_moon.dust_mass + the_moon.gas_mass) > this_planet.dust_mass + this_planet.gas_mass {
                                // if the moon has more mass than the planet, switch them around
                                
                                let temp_dust = this_planet.dust_mass
                                let temp_gas = this_planet.gas_mass;
                                let temp_mass = this_planet.mass;
                                
                                this_planet.dust_mass = the_moon.dust_mass
                                this_planet.gas_mass = the_moon.gas_mass
                                this_planet.mass = the_moon.mass
                                
                                the_moon.dust_mass = temp_dust
                                the_moon.gas_mass = temp_gas
                                the_moon.mass = temp_mass
                            }
                            
                            if (this_planet.first_moon == nil) {
                                // if this planet doesn't already have a first moon, make this one it's first moon
                                this_planet.first_moon = the_moon
                            } else {
                                // otherwise insert it into the list of moons
                                the_moon.next_planet = this_planet.first_moon
                                this_planet.first_moon = the_moon
                            }
                            finished = true
                            print("Moon captured: \(this_planet.a) AU \(this_planet.mass * SUN_MASS_IN_EARTH_MASSES) <- \(mass * SUN_MASS_IN_EARTH_MASSES)")
                            
                        } else {
                            var escape_reason = ""
                            if (mass * SUN_MASS_IN_EARTH_MASSES) >= 2.5 {
                                escape_reason = ", too big"
                            } else if ((mass * SUN_MASS_IN_EARTH_MASSES) <= 0.0001) {
                                escape_reason = ", too small"
                            }
                            print("Moon escapes: \(this_planet.a) AU (\(this_planet.mass * SUN_MASS_IN_EARTH_MASSES))\(mass * SUN_MASS_IN_EARTH_MASSES) \(escape_reason)")
                        }
                    }
                }
                if !finished {
                    print("Collision between two planetesimals! \(this_planet.a) AU (\(this_planet.mass * SUN_MASS_IN_EARTH_MASSES)EM) + \(a) AU (\(mass * SUN_MASS_IN_EARTH_MASSES)EM = \(dust_mass * SUN_MASS_IN_EARTH_MASSES)EMd + \(gas_mass * SUN_MASS_IN_EARTH_MASSES)EMg [\(crit_mass * SUN_MASS_IN_EARTH_MASSES)EM]) -> \(new_a) AU (\(e))")
                    temp = this_planet.mass + mass
                    accrete_dust(seed_mass: &temp, new_dust: &new_dust, new_gas: &new_gas, a: new_a, e: e, crit_mass: stell_luminosity_ratio, body_inner_bound: body_inner_bound, body_outer_bound: body_outer_bound)
                    this_planet.a = new_a
                    this_planet.e = e
                    this_planet.mass = temp
                    this_planet.dust_mass += dust_mass + new_dust
                    this_planet.gas_mass += gas_mass + new_gas
                    if (temp >= crit_mass) {
                        this_planet.gas_giant = true
                    }
                    
                    // Walk the list of planets and move the linked list around to re-order the planetisimals
                    while (the_planet!.next_planet != nil && the_planet!.next_planet!.a < new_a) {
                        next_planet = the_planet!.next_planet
                        if (the_planet == planet_head) {
                            planet_head = next_planet
                        } else {
                            prev_planet?.next_planet = next_planet
                        }
                        the_planet?.next_planet = next_planet?.next_planet
                        next_planet?.next_planet = the_planet
                        prev_planet = next_planet
                    }
                }
                finished = true
                break
            } else {
                prev_planet = the_planet
            }
        }
    }
    
    
    //void coalesce_planetesimals(long double a, long double e, long double mass, long double crit_mass,
    //                            long double dust_mass, long double gas_mass,
    //                            long double stell_luminosity_ratio,
    //                            long double body_inner_bound, long double body_outer_bound,
    //                            int            do_moons)
    //{
    //    planet_pointer    the_planet;
    //    planet_pointer    next_planet;
    //    planet_pointer    prev_planet;
    //    int             finished;
    //    long double     temp;
    //    long double     diff;
    //    long double     dist1;
    //    long double     dist2;
    //
    //    finished = FALSE;
    //    prev_planet = NULL;
    //
    //// First we try to find an existing planet with an over-lapping orbit.
    //
    //    for (the_planet = planet_head;
    //         the_planet != NULL;
    //         the_planet = the_planet->next_planet)
    //    {
    //        diff = the_planet->a - a;
    //
    //        if ((diff > 0.0))
    //        {
    //            dist1 = (a * (1.0 + e) * (1.0 + reduced_mass)) - a;
    //            /* x aphelion     */
    //            reduced_mass = pow((the_planet->mass / (1.0 + the_planet->mass)),(1.0 / 4.0));
    //            dist2 = the_planet->a
    //                - (the_planet->a * (1.0 - the_planet->e) * (1.0 - reduced_mass));
    //        }
    //        else
    //        {
    //            dist1 = a - (a * (1.0 - e) * (1.0 - reduced_mass));
    //            /* x perihelion */
    //            reduced_mass = pow((the_planet->mass / (1.0 + the_planet->mass)),(1.0 / 4.0));
    //            dist2 = (the_planet->a * (1.0 + the_planet->e) * (1.0 + reduced_mass))
    //                - the_planet->a;
    //        }
    //
    //        if (((fabs(diff) <= fabs(dist1)) || (fabs(diff) <= fabs(dist2))))
    //        {
    //            long double new_dust = 0;
    //            long double    new_gas = 0;
    //            long double new_a = (the_planet->mass + mass) /
    //                                ((the_planet->mass / the_planet->a) + (mass / a));
    //
    //            temp = the_planet->mass * sqrt(the_planet->a) * sqrt(1.0 - pow(the_planet->e,2.0));
    //            temp = temp + (mass * sqrt(a) * sqrt(sqrt(1.0 - pow(e,2.0))));
    //            temp = temp / ((the_planet->mass + mass) * sqrt(new_a));
    //            temp = 1.0 - pow(temp,2.0);
    //            if (((temp < 0.0) || (temp >= 1.0)))
    //                temp = 0.0;
    //            e = sqrt(temp);
    //
    //            if (do_moons)
    //            {
    //                long double existing_mass = 0.0;
    //
    //                if (the_planet->first_moon != NULL)
    //                {
    //                    planet_pointer    m;
    //
    //                    for (m = the_planet->first_moon;
    //                         m != NULL;
    //                         m = m->next_planet)
    //                    {
    //                        existing_mass += m->mass;
    //                    }
    //                }
    //
    //                if (mass < crit_mass)
    //                {
    //                    if ((mass * SUN_MASS_IN_EARTH_MASSES) < 2.5
    //                     && (mass * SUN_MASS_IN_EARTH_MASSES) > .0001
    //                     && existing_mass < (the_planet->mass * .05)
    //                       )
    //                    {
    //                        planet_pointer    the_moon = (planets *)malloc(sizeof(planets));
    //
    //                        the_moon->type             = tUnknown;
    //    /*                     the_moon->a             = a; */
    //    /*                     the_moon->e             = e; */
    //                        the_moon->mass             = mass;
    //                        the_moon->dust_mass     = dust_mass;
    //                        the_moon->gas_mass         = gas_mass;
    //                        the_moon->atmosphere     = NULL;
    //                        the_moon->next_planet     = NULL;
    //                        the_moon->first_moon     = NULL;
    //                        the_moon->gas_giant     = FALSE;
    //                        the_moon->atmosphere    = NULL;
    //                        the_moon->albedo        = 0;
    //                        the_moon->gases            = 0;
    //                        the_moon->surf_temp        = 0;
    //                        the_moon->high_temp        = 0;
    //                        the_moon->low_temp        = 0;
    //                        the_moon->max_temp        = 0;
    //                        the_moon->min_temp        = 0;
    //                        the_moon->greenhs_rise    = 0;
    //                        the_moon->minor_moons     = 0;
    //
    //                        if ((the_moon->dust_mass + the_moon->gas_mass)
    //                          > (the_planet->dust_mass + the_planet->gas_mass))
    //                        {
    //                            long double    temp_dust = the_planet->dust_mass;
    //                            long double temp_gas  = the_planet->gas_mass;
    //                            long double temp_mass = the_planet->mass;
    //
    //                            the_planet->dust_mass = the_moon->dust_mass;
    //                            the_planet->gas_mass  = the_moon->gas_mass;
    //                            the_planet->mass      = the_moon->mass;
    //
    //                            the_moon->dust_mass   = temp_dust;
    //                            the_moon->gas_mass    = temp_gas;
    //                            the_moon->mass        = temp_mass;
    //                        }
    //
    //                        if (the_planet->first_moon == NULL)
    //                            the_planet->first_moon = the_moon;
    //                        else
    //                        {
    //                            the_moon->next_planet = the_planet->first_moon;
    //                            the_planet->first_moon = the_moon;
    //                        }
    //
    //                        finished = TRUE;
    //
    //                        if (flag_verbose & 0x0100)
    //                            fprintf (stderr, "Moon Captured... "
    //                                     "%5.3Lf AU (%.2LfEM) <- %.2LfEM\n",
    //                                    the_planet->a, the_planet->mass * SUN_MASS_IN_EARTH_MASSES,
    //                                    mass * SUN_MASS_IN_EARTH_MASSES
    //                                    );
    //                    }
    //                    else
    //                    {
    //                        if (flag_verbose & 0x0100)
    //                            fprintf (stderr, "Moon Escapes... "
    //                                     "%5.3Lf AU (%.2LfEM)%s %.2LfEM%s\n",
    //                                    the_planet->a, the_planet->mass * SUN_MASS_IN_EARTH_MASSES,
    //                                    existing_mass < (the_planet->mass * .05) ? "" : " (big moons)",
    //                                    mass * SUN_MASS_IN_EARTH_MASSES,
    //                                    (mass * SUN_MASS_IN_EARTH_MASSES) >= 2.5 ? ", too big" :
    //                                      (mass * SUN_MASS_IN_EARTH_MASSES) <= .0001 ? ", too small" : ""
    //                                    );
    //                    }
    //                }
    //            }
    //
    //            if (!finished)
    //            {
    //                if (flag_verbose & 0x0100)
    //                        fprintf (stderr, "Collision between two planetesimals! "
    //                                "%4.2Lf AU (%.2LfEM) + %4.2Lf AU (%.2LfEM = %.2LfEMd + %.2LfEMg [%.3LfEM])-> %5.3Lf AU (%5.3Lf)\n",
    //                                the_planet->a, the_planet->mass * SUN_MASS_IN_EARTH_MASSES,
    //                                a, mass * SUN_MASS_IN_EARTH_MASSES,
    //                                dust_mass * SUN_MASS_IN_EARTH_MASSES, gas_mass * SUN_MASS_IN_EARTH_MASSES,
    //                                crit_mass * SUN_MASS_IN_EARTH_MASSES,
    //                                new_a, e);
    //
    //                temp = the_planet->mass + mass;
    //                accrete_dust(&temp, &new_dust, &new_gas,
    //                             new_a,e,stell_luminosity_ratio,
    //                             body_inner_bound,body_outer_bound);
    //
    //                the_planet->a = new_a;
    //                the_planet->e = e;
    //                the_planet->mass = temp;
    //                the_planet->dust_mass += dust_mass + new_dust;
    //                the_planet->gas_mass += gas_mass + new_gas;
    //                if (temp >= crit_mass)
    //                    the_planet->gas_giant = TRUE;
    //
    //                while (the_planet->next_planet != NULL && the_planet->next_planet->a < new_a)
    //                {
    //                    next_planet = the_planet->next_planet;
    //
    //                    if (the_planet == planet_head)
    //                        planet_head = next_planet;
    //                    else
    //                        prev_planet->next_planet = next_planet;
    //
    //                    the_planet->next_planet = next_planet->next_planet;
    //                    next_planet->next_planet = the_planet;
    //                    prev_planet = next_planet;
    //                }
    //            }
    //
    //            finished = TRUE;
    //            break;
    //        }
    //        else
    //        {
    //            prev_planet = the_planet;
    //        }
    //    }
    //
    //    if (!(finished))            // Planetesimals didn't collide. Make it a planet.
    //    {
    //        the_planet = (planets *)malloc(sizeof(planets));
    //
    //        the_planet->type             = tUnknown;
    //        the_planet->a                 = a;
    //        the_planet->e                 = e;
    //        the_planet->mass             = mass;
    //        the_planet->dust_mass         = dust_mass;
    //        the_planet->gas_mass         = gas_mass;
    //        the_planet->atmosphere         = NULL;
    //        the_planet->first_moon         = NULL;
    //        the_planet->atmosphere        = NULL;
    //        the_planet->albedo            = 0;
    //        the_planet->gases            = 0;
    //        the_planet->surf_temp        = 0;
    //        the_planet->high_temp        = 0;
    //        the_planet->low_temp        = 0;
    //        the_planet->max_temp        = 0;
    //        the_planet->min_temp        = 0;
    //        the_planet->greenhs_rise    = 0;
    //        the_planet->minor_moons     = 0;
    //
    //        if ((mass >= crit_mass))
    //            the_planet->gas_giant = TRUE;
    //        else
    //            the_planet->gas_giant = FALSE;
    //
    //        if ((planet_head == NULL))
    //        {
    //            planet_head = the_planet;
    //            the_planet->next_planet = NULL;
    //        }
    //        else if ((a < planet_head->a))
    //        {
    //            the_planet->next_planet = planet_head;
    //            planet_head = the_planet;
    //        }
    //        else if ((planet_head->next_planet == NULL))
    //        {
    //            planet_head->next_planet = the_planet;
    //            the_planet->next_planet = NULL;
    //        }
    //        else
    //        {
    //            next_planet = planet_head;
    //            while (((next_planet != NULL) && (next_planet->a < a)))
    //            {
    //                prev_planet = next_planet;
    //                next_planet = next_planet->next_planet;
    //            }
    //            the_planet->next_planet = next_planet;
    //            prev_planet->next_planet = the_planet;
    //        }
    //    }
    //}
    
    // primary entry point? - called from SolarSystem
    mutating func dist_planetary_masses(stell_mass_ratio: Double, stell_luminosity_ratio: Double, inner_dust: Double, outer_dust:Double, outer_planet_limit: Double, dust_density_coeff: Double, seed_system: Planet?, do_moons: Bool) -> Planet? {
        var a: Double // distance, in AU
        var e: Double // eccentricity of orbit
        var mass: Double = PROTOPLANET_MASS // units of Solar Mass
        var dust_mass: Double = 0
        var gas_mass: Double = 0
        var crit_mass: Double
        var planet_inner_bound: Double
        var planet_outer_bound: Double
        var seeds: Planet? = seed_system
        
        set_initial_conditions(inner_limit_of_dust: inner_dust, outer_limit_of_dust: outer_dust)
        planet_inner_bound = nearest_planet(stell_mass_ratio: stell_mass_ratio)
        if planet_outer_bound == 0 {
            planet_outer_bound = farthest_planet(stell_mass_ratio: stell_mass_ratio)
        } else {
            planet_outer_bound = outer_planet_limit
        }
        
        // The general flow is to seed planetesimals that start gravitational accretion
        // (called 'gravitational instability' in current (~2020) astrophysics research)
        // and continue to deploy "seeds" until the dust of the disk is "consumed".
        while(dust_left) {
            if let definitely_seed = seeds {
                a = definitely_seed.a
                e = definitely_seed.e
                seeds = definitely_seed.next_planet
            } else {
                a = prng.random_number(in: planet_inner_bound...planet_outer_bound)
                e = prng.random_eccentricity()
            }

            print("Checking \(a) AU")

            if dust_available(inside_range: inner_effect_limit(a: a, e: e, mass: mass),
                              outside_range: outer_effect_limit(a: a, e: e, mass: mass)) {
                print("Injecting protoplanet at \(a) AU")
                
                // NOTE(heckj): this is used in collect_dist, which is called from inside accrete_dust - and doesn't
                // appear to be changing, so this calculation likely doesn't need to be within the loop
                // and could be part of setting up initial conditions.
                dust_density = dust_density_coeff * sqrt(stell_mass_ratio) * exp(-ALPHA * pow(a,(1.0 / N)))
                
                // Determine the mass (in solar masses) at which a body will start accumulating gasses
                crit_mass = critical_limit(orb_radius: a, eccentricity: e, stell_luminosity_ratio: stell_luminosity_ratio)
                
                accrete_dust(seed_mass: &mass, new_dust: &dust_mass, new_gas: &gas_mass, a: a, e: e, crit_mass: crit_mass, body_inner_bound: planet_inner_bound, body_outer_bound: planet_outer_bound)
                dust_mass += PROTOPLANET_MASS
                
                if mass > PROTOPLANET_MASS {
                    coalesce_planetesimals(a: a, e: e, mass: mass, crit_mass: crit_mass, dust_mass: dust_mass, gas_mass: gas_mass, stell_luminosity_ratio: stell_luminosity_ratio, body_inner_bound: planet_inner_bound, body_outer_bound: planet_outer_bound, do_moons: do_moons)
                } else {
                    print("failed..")
                }

            }
        }
        return(planet_head)
    }
    
    //planet_pointer dist_planetary_masses(long double stell_mass_ratio,
    //                                     long double stell_luminosity_ratio,
    //                                     long double inner_dust,
    //                                     long double outer_dust,
    //                                     long double outer_planet_limit,
    //                                     long double dust_density_coeff,
    //                                     planet_pointer seed_system,
    //                                     int         do_moons)
    //{
    //    long double     a;
    //    long double     e;
    //    long double     mass;
    //    long double        dust_mass;
    //    long double        gas_mass;
    //    long double     crit_mass;
    //    long double     planet_inner_bound;
    //    long double     planet_outer_bound;
    //    planet_pointer     seeds = seed_system;
    //
    //    set_initial_conditions(inner_dust,outer_dust);
    //    planet_inner_bound = nearest_planet(stell_mass_ratio);
    //
    //    if (outer_planet_limit == 0)
    //        planet_outer_bound = farthest_planet(stell_mass_ratio);
    //    else
    //        planet_outer_bound = outer_planet_limit;
    //
    //    while (dust_left)
    //    {
    //        if (seeds != NULL)
    //        {
    //            a = seeds->a;
    //            e = seeds->e;
    //            seeds = seeds->next_planet;
    //        }
    //        else
    //        {
    //            a = random_number(planet_inner_bound,planet_outer_bound);
    //            e = random_eccentricity( );
    //        }
    //
    //        mass      = PROTOPLANET_MASS;
    //        dust_mass = 0;
    //        gas_mass  = 0;
    //
    //        if (flag_verbose & 0x0200)
    //            fprintf (stderr, "Checking %Lg AU.\n",a);
    //
    //        if (dust_available(inner_effect_limit(a, e, mass),
    //                           outer_effect_limit(a, e, mass)))
    //        {
    //            if (flag_verbose & 0x0100)
    //                fprintf (stderr, "Injecting protoplanet at %Lg AU.\n", a);
    //
    //            dust_density = dust_density_coeff * sqrt(stell_mass_ratio)
    //                           * exp(-ALPHA * pow(a,(1.0 / N)));
    //            crit_mass = critical_limit(a,e,stell_luminosity_ratio);
    //            accrete_dust(&mass, &dust_mass, &gas_mass,
    //                         a,e,crit_mass,
    //                         planet_inner_bound,
    //                         planet_outer_bound);
    //
    //            dust_mass += PROTOPLANET_MASS;
    //
    //            if (mass > PROTOPLANET_MASS)
    //                coalesce_planetesimals(a,e,mass,crit_mass,
    //                                       dust_mass, gas_mass,
    //                                       stell_luminosity_ratio,
    //                                       planet_inner_bound,planet_outer_bound,
    //                                       do_moons);
    //            else if (flag_verbose & 0x0100)
    //                fprintf (stderr, ".. failed due to large neighbor.\n");
    //        }
    //        else if (flag_verbose & 0x0200)
    //            fprintf (stderr, ".. failed.\n");
    //    }
    //    return(planet_head);
    //}

}

//
//void free_dust (dust_pointer head)
//{
//    dust_pointer    node;
//    dust_pointer    next;
//
//    for(node = head;
//        node != NULL;
//        node = next)
//    {
//        next = node->next_band;
//        free (node);
//    }
//
//}
//
//void free_planet (planet_pointer head)
//{
//    planet_pointer    node;
//    planet_pointer    next;
//
//    for(node = head;
//        node != NULL;
//        node = next)
//    {
//        next = node->next_planet;
//
//        free (node);
//    }
//}
//
//void free_generations()
//{
//    gen_pointer    node;
//    gen_pointer    next;
//
//    for(node = hist_head;
//        node != NULL;
//        node = next)
//    {
//        next = node->next;
//
//        if (node->dusts)
//            free_dust (node->dusts);
//
//        if (node->planets)
//            free_planet (node->planets);
//
//        free (node);
//    }
//
//    if (dust_head != NULL)
//        free_dust (dust_head);
//
//    if (planet_head != NULL)
//        free_planet (planet_head);
//
//    dust_head = NULL;
//    planet_head = NULL;
//    hist_head = NULL;
//}
//
//void free_atmosphere(planet_pointer head)
//{
//    planet_pointer    node;
//
//    for (node = head;
//         node != NULL;
//         node = node->next_planet)
//    {
//        if (node->atmosphere != NULL)
//        {
//            free(node->atmosphere);
//
//            node->atmosphere = NULL;
//        }
//
//        if (node->first_moon != NULL)
//        {
//            free_atmosphere(node->first_moon);
//        }
//    }
//}

