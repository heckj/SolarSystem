//  Constants.swift

import Foundation

//#define PI                        (3.1415926536)
// Using Double.pi

//#define RADIANS_PER_ROTATION    (2.0 * PI)
let RADIANS_PER_ROTATION = Double.pi * 2.0

let ECCENTRICITY_COEFF  = 0.077            /* Dole's was 0.077            */
let PROTOPLANET_MASS  = 1.0E-15         /* Units of solar masses    */
let CHANGE_IN_EARTH_ANG_VEL = -1.3E-15        /* Units of radians/sec/year*/
let SOLAR_MASS_IN_GRAMS = 1.989E33        /* Units of grams            */
let SOLAR_MASS_IN_KILOGRAMS    = 1.989E30        /* Units of kg                */
let EARTH_MASS_IN_GRAMS        = 5.977E27        /* Units of grams            */
let EARTH_RADIUS             = 6.378E8        /* Units of cm                */
let EARTH_DENSITY            = 5.52            /* Units of g/cc            */
let KM_EARTH_RADIUS            = 6378.0        /* Units of km                */
let EARTH_ACCELERATION       = 980.7            /* Units of cm/sec2            */
let EARTH_AXIAL_TILT        = 23.4            /* Units of degrees            */
let EARTH_EXOSPHERE_TEMP    = 1273.0        /* Units of degrees Kelvin    */
let SUN_MASS_IN_EARTH_MASSES = 332775.64
let ASTEROID_MASS_LIMIT        = 0.001            /* Units of Earth Masses    */
let EARTH_EFFECTIVE_TEMP    = 250.0            /* Units of degrees Kelvin (was 255)    */
let CLOUD_COVERAGE_FACTOR    = 1.839E-8        /* Km2/kg                    */
let EARTH_WATER_MASS_PER_AREA     = 3.83E15    /* grams per square km        */
let EARTH_SURF_PRES_IN_MILLIBARS = 1013.25
let EARTH_SURF_PRES_IN_MMHG    = 760.0            /* Dole p. 15                */
let EARTH_SURF_PRES_IN_PSI    = 14.696        /* Pounds per square inch    */
let MMHG_TO_MILLIBARS = (EARTH_SURF_PRES_IN_MILLIBARS / EARTH_SURF_PRES_IN_MMHG)
let PSI_TO_MILLIBARS = (EARTH_SURF_PRES_IN_MILLIBARS / EARTH_SURF_PRES_IN_PSI)
let H20_ASSUMED_PRESSURE   = (47.0 * MMHG_TO_MILLIBARS) /* Dole p. 15      */
let MIN_O2_IPP  =  (72.0 * MMHG_TO_MILLIBARS)    /* Dole, p. 15                */
let MAX_O2_IPP  =  (400.0 * MMHG_TO_MILLIBARS)    /* Dole, p. 15                */
let MAX_HE_IPP  =  (61000.0 * MMHG_TO_MILLIBARS)    /* Dole, p. 16            */
let MAX_NE_IPP  =  (3900.0 * MMHG_TO_MILLIBARS)    /* Dole, p. 16                */
let MAX_N2_IPP  =  (2330.0 * MMHG_TO_MILLIBARS)    /* Dole, p. 16                */
let MAX_AR_IPP  =  (1220.0 * MMHG_TO_MILLIBARS)    /* Dole, p. 16                */
let MAX_KR_IPP  =  (350.0 * MMHG_TO_MILLIBARS)    /* Dole, p. 16                */
let MAX_XE_IPP  =  (160.0 * MMHG_TO_MILLIBARS)    /* Dole, p. 16                */
let MAX_CO2_IPP = (7.0 * MMHG_TO_MILLIBARS)    /* Dole, p. 16                */
let MAX_HABITABLE_PRESSURE = (118 * PSI_TO_MILLIBARS)    /* Dole, p. 16        */
// The next gases are listed as poisonous in parts per million by volume at 1 atm:
let PPM_PRSSURE = (EARTH_SURF_PRES_IN_MILLIBARS / 1000000.0)
let MAX_F_IPP    = (0.1 * PPM_PRSSURE)            /* Dole, p. 18                */
let MAX_CL_IPP     = (1.0 * PPM_PRSSURE)            /* Dole, p. 18                */
let MAX_NH3_IPP    = (100.0 * PPM_PRSSURE)        /* Dole, p. 18                */
let MAX_O3_IPP    = (0.1 * PPM_PRSSURE)            /* Dole, p. 18                */
let MAX_CH4_IPP    = (50000.0 * PPM_PRSSURE)        /* Dole, p. 18                */



let EARTH_CONVECTION_FACTOR = 0.43            /* from Hart, eq.20            */
//      FREEZING_POINT_OF_WATER (273.0)            /* Units of degrees Kelvin    */
let FREEZING_POINT_OF_WATER = 273.15        /* Units of degrees Kelvin    */
//      EARTH_AVERAGE_CELSIUS   (15.5)            /* Average Earth Temperature */
let EARTH_AVERAGE_CELSIUS   = 14.0            /* Average Earth Temperature */
let EARTH_AVERAGE_KELVIN    = (EARTH_AVERAGE_CELSIUS + FREEZING_POINT_OF_WATER)
let DAYS_IN_A_YEAR           = 365.256        /* Earth days per Earth year*/
//        gas_retention_threshold = 5.0;          /* ratio of esc vel to RMS vel */
let GAS_RETENTION_THRESHOLD = 6.0            /* ratio of esc vel to RMS vel */

let ICE_ALBEDO               = 0.7
let CLOUD_ALBEDO            = 0.52
let GAS_GIANT_ALBEDO        = 0.5            /* albedo of a gas giant    */
let AIRLESS_ICE_ALBEDO        = 0.5
let EARTH_ALBEDO            = 0.3            /* was .33 for a while */
let GREENHOUSE_TRIGGER_ALBEDO = 0.20
let ROCKY_ALBEDO            = 0.15
let ROCKY_AIRLESS_ALBEDO    = 0.07
let WATER_ALBEDO            = 0.04

let SECONDS_PER_HOUR        = 3600.0
let CM_PER_AU                = 1.495978707E13/* number of cm in an AU    */
let CM_PER_KM                = 1.0E5             /* number of cm in a km        */
let KM_PER_AU               = (CM_PER_AU / CM_PER_KM)
let CM_PER_METER            = 100.0
let MILLIBARS_PER_BAR       = 1000.00

let GRAV_CONSTANT            = 6.672E-8        /* units of dyne cm2/gram2    */
let MOLAR_GAS_CONST           = 8314.41        /* units: g*m2/(sec2*K*mol) */
let K                        = 50.0             /* K = gas/dust ratio        */
let B                        = 1.2E-5        /* Used in Crit_mass calc    */
let DUST_DENSITY_COEFF       = 2.0E-3        /* A in Dole's paper        */
let ALPHA                    = 5.0           /* Used in density calcs    */
let N                        = 3.0            /* Used in density calcs    */
let J                        = 1.46E-19        /* Used in day-length calcs (cm2/sec2 g) */
let INCREDIBLY_LARGE_NUMBER = Double.max

/*    Now for a few molecular weights (used for RMS velocity calcs):       */
/*    This table is from Dole's book "Habitable Planets for Man", p. 38  */

let ATOMIC_HYDROGEN            = 1.0    /* H   */
let MOL_HYDROGEN            = 2.0    /* H2  */
let HELIUM                    = 4.0    /* He  */
let ATOMIC_NITROGEN            = 14.0    /* N   */
let ATOMIC_OXYGEN            = 16.0    /* O   */
let METHANE                    = 16.0    /* CH4 */
let AMMONIA                    = 17.0    /* NH3 */
let WATER_VAPOR                = 18.0    /* H2O */
let NEON                    = 20.2    /* Ne  */
let MOL_NITROGEN            = 28.0    /* N2  */
let CARBON_MONOXIDE            = 28.0    /* CO  */
let NITRIC_OXIDE            = 30.0    /* NO  */
let MOL_OXYGEN                = 32.0    /* O2  */
let HYDROGEN_SULPHIDE        = 34.1    /* H2S */
let ARGON                    = 39.9    /* Ar  */
let CARBON_DIOXIDE            = 44.0    /* CO2 */
let NITROUS_OXIDE            = 44.0    /* N2O */
let NITROGEN_DIOXIDE        = 46.0    /* NO2 */
let OZONE                    = 48.0    /* O3  */
let SULPH_DIOXIDE            = 64.1    /* SO2 */
let SULPH_TRIOXIDE            = 80.1    /* SO3 */
let KRYPTON                    = 83.8    /* Kr  */
let XENON                    = 131.3 /* Xe  */

//    And atomic numbers, for use in ChemTable indexes
let AN_H   = 1
let AN_HE  =    2
let AN_N   =   7
let AN_O   =   8
let AN_F   =   9
let AN_NE  =    10
let AN_P   =   15
let AN_CL  =    17
let AN_AR  =    18
let AN_BR  =    35
let AN_KR  =    36
let AN_I   =   53
let AN_XE  =    54
let AN_HG  =    80
let AN_AT  =    85
let AN_RN  =    86
let AN_FR  =    87

let AN_NH3 =     900
let AN_H2O =     901
let AN_CO2 =     902
let AN_O3  =    903
let AN_CH4 =     904
let AN_CH3CH2OH =   905

/*    The following defines are used in the kothari_radius function in    */
/*    file enviro.c.                                                        */
let A1_20                    = 6.485E12        /* All units are in cgs system.     */
let A2_20                    = 4.0032E-8        /*     ie: cm, g, dynes, etc.         */
let BETA_20                    = 5.71E12

let JIMS_FUDGE                = 1.004

/*     The following defines are used in determining the fraction of a planet     */
/*    covered with clouds in function cloud_fraction in file enviro.c.         */
let Q1_36                    = 1.258E19        /* grams    */
let Q2_36                    = 0.0698        /* 1/Kelvin */
