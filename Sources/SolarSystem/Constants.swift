//  Constants.swift

// #define PI                        (3.1415926536)
// Using Double.pi

// #define RADIANS_PER_ROTATION    (2.0 * PI)
public let RADIANS_PER_ROTATION = Double.pi * 2.0

public let ECCENTRICITY_COEFF = 0.077 /* Dole's was 0.077            */
public let PROTOPLANET_MASS = 1.0e-15 /* Units of solar masses    */
public let CHANGE_IN_EARTH_ANG_VEL = -1.3e-15 /* Units of radians/sec/year*/
public let SOLAR_MASS_IN_GRAMS = 1.989e33 /* Units of grams            */
public let SOLAR_MASS_IN_KILOGRAMS = 1.989e30 /* Units of kg                */
public let EARTH_MASS_IN_GRAMS = 5.977e27 /* Units of grams            */
public let EARTH_RADIUS = 6.378e8 /* Units of cm                */
public let EARTH_DENSITY = 5.52 /* Units of g/cc            */
public let KM_EARTH_RADIUS = 6378.0 /* Units of km                */
public let EARTH_ACCELERATION = 980.7 /* Units of cm/sec2            */
public let EARTH_AXIAL_TILT = 23.4 /* Units of degrees            */
public let EARTH_EXOSPHERE_TEMP = 1273.0 /* Units of degrees Kelvin    */
public let SUN_MASS_IN_EARTH_MASSES = 332_775.64
public let ASTEROID_MASS_LIMIT = 0.001 /* Units of Earth Masses    */
public let EARTH_EFFECTIVE_TEMP = 250.0 /* Units of degrees Kelvin (was 255)    */
public let CLOUD_COVERAGE_FACTOR = 1.839e-8 /* Km2/kg                    */
public let EARTH_WATER_MASS_PER_AREA = 3.83e15 /* grams per square km        */
public let EARTH_SURF_PRES_IN_MILLIBARS = 1013.25
public let EARTH_SURF_PRES_IN_MMHG = 760.0 /* Dole p. 15                */
public let EARTH_SURF_PRES_IN_PSI = 14.696 /* Pounds per square inch    */
public let MMHG_TO_MILLIBARS = (EARTH_SURF_PRES_IN_MILLIBARS / EARTH_SURF_PRES_IN_MMHG)
public let PSI_TO_MILLIBARS = (EARTH_SURF_PRES_IN_MILLIBARS / EARTH_SURF_PRES_IN_PSI)
public let H20_ASSUMED_PRESSURE = (47.0 * MMHG_TO_MILLIBARS) /* Dole p. 15      */
public let MIN_O2_IPP = (72.0 * MMHG_TO_MILLIBARS) /* Dole, p. 15                */
public let MAX_O2_IPP = (400.0 * MMHG_TO_MILLIBARS) /* Dole, p. 15                */
public let MAX_HE_IPP = (61000.0 * MMHG_TO_MILLIBARS) /* Dole, p. 16            */
public let MAX_NE_IPP = (3900.0 * MMHG_TO_MILLIBARS) /* Dole, p. 16                */
public let MAX_N2_IPP = (2330.0 * MMHG_TO_MILLIBARS) /* Dole, p. 16                */
public let MAX_AR_IPP = (1220.0 * MMHG_TO_MILLIBARS) /* Dole, p. 16                */
public let MAX_KR_IPP = (350.0 * MMHG_TO_MILLIBARS) /* Dole, p. 16                */
public let MAX_XE_IPP = (160.0 * MMHG_TO_MILLIBARS) /* Dole, p. 16                */
public let MAX_CO2_IPP = (7.0 * MMHG_TO_MILLIBARS) /* Dole, p. 16                */
public let MAX_HABITABLE_PRESSURE = (118 * PSI_TO_MILLIBARS) /* Dole, p. 16        */
// The next gases are listed as poisonous in parts per million by volume at 1 atm:
public let PPM_PRSSURE = (EARTH_SURF_PRES_IN_MILLIBARS / 1_000_000.0)
public let MAX_F_IPP = (0.1 * PPM_PRSSURE) /* Dole, p. 18                */
public let MAX_CL_IPP = (1.0 * PPM_PRSSURE) /* Dole, p. 18                */
public let MAX_NH3_IPP = (100.0 * PPM_PRSSURE) /* Dole, p. 18                */
public let MAX_O3_IPP = (0.1 * PPM_PRSSURE) /* Dole, p. 18                */
public let MAX_CH4_IPP = (50000.0 * PPM_PRSSURE) /* Dole, p. 18                */

public let EARTH_CONVECTION_FACTOR = 0.43 /* from Hart, eq.20            */
//      FREEZING_POINT_OF_WATER (273.0)            /* Units of degrees Kelvin    */
public let FREEZING_POINT_OF_WATER = 273.15 /* Units of degrees Kelvin    */
//      EARTH_AVERAGE_CELSIUS   (15.5)            /* Average Earth Temperature */
public let EARTH_AVERAGE_CELSIUS = 14.0 /* Average Earth Temperature */
public let EARTH_AVERAGE_KELVIN = (EARTH_AVERAGE_CELSIUS + FREEZING_POINT_OF_WATER)
public let DAYS_IN_A_YEAR = 365.256 /* Earth days per Earth year*/
//        gas_retention_threshold = 5.0;          /* ratio of esc vel to RMS vel */
public let GAS_RETENTION_THRESHOLD = 6.0 /* ratio of esc vel to RMS vel */

public let ICE_ALBEDO = 0.7
public let CLOUD_ALBEDO = 0.52
public let GAS_GIANT_ALBEDO = 0.5 /* albedo of a gas giant    */
public let AIRLESS_ICE_ALBEDO = 0.5
public let EARTH_ALBEDO = 0.3 /* was .33 for a while */
public let GREENHOUSE_TRIGGER_ALBEDO = 0.20
public let ROCKY_ALBEDO = 0.15
public let ROCKY_AIRLESS_ALBEDO = 0.07
public let WATER_ALBEDO = 0.04

public let SECONDS_PER_HOUR = 3600.0
public let CM_PER_AU = 1.495978707e13 /* number of cm in an AU    */
public let CM_PER_KM = 1.0e5 /* number of cm in a km        */
public let KM_PER_AU = (CM_PER_AU / CM_PER_KM)
public let CM_PER_METER = 100.0
public let MILLIBARS_PER_BAR = 1000.00

public let GRAV_CONSTANT = 6.672e-8 /* units of dyne cm2/gram2    */
public let MOLAR_GAS_CONST = 8314.41 /* units: g*m2/(sec2*K*mol) */
public let K = 50.0 /* K = gas/dust ratio        */
public let B = 1.2e-5 /* Used in Crit_mass calc    */
public let DUST_DENSITY_COEFF = 2.0e-3 /* A in Dole's paper        */
public let ALPHA = 5.0 /* Used in density calcs    */
public let N = 3.0 /* Used in density calcs    */
public let J = 1.46e-19 /* Used in day-length calcs (cm2/sec2 g) */
public let INCREDIBLY_LARGE_NUMBER = Double.greatestFiniteMagnitude

/*    Now for a few molecular weights (used for RMS velocity calcs):       */
/*    This table is from Dole's book "Habitable Planets for Man", p. 38  */

public let ATOMIC_HYDROGEN = 1.0 /* H   */
public let MOL_HYDROGEN = 2.0 /* H2  */
public let HELIUM = 4.0 /* He  */
public let ATOMIC_NITROGEN = 14.0 /* N   */
public let ATOMIC_OXYGEN = 16.0 /* O   */
public let METHANE = 16.0 /* CH4 */
public let AMMONIA = 17.0 /* NH3 */
public let WATER_VAPOR = 18.0 /* H2O */
public let NEON = 20.2 /* Ne  */
public let MOL_NITROGEN = 28.0 /* N2  */
public let CARBON_MONOXIDE = 28.0 /* CO  */
public let NITRIC_OXIDE = 30.0 /* NO  */
public let MOL_OXYGEN = 32.0 /* O2  */
public let HYDROGEN_SULPHIDE = 34.1 /* H2S */
public let ARGON = 39.9 /* Ar  */
public let CARBON_DIOXIDE = 44.0 /* CO2 */
public let NITROUS_OXIDE = 44.0 /* N2O */
public let NITROGEN_DIOXIDE = 46.0 /* NO2 */
public let OZONE = 48.0 /* O3  */
public let SULPH_DIOXIDE = 64.1 /* SO2 */
public let SULPH_TRIOXIDE = 80.1 /* SO3 */
public let KRYPTON = 83.8 /* Kr  */
public let XENON = 131.3 /* Xe  */

//    And atomic numbers, for use in ChemTable indexes
public let AN_H = 1
public let AN_HE = 2
public let AN_N = 7
public let AN_O = 8
public let AN_F = 9
public let AN_NE = 10
public let AN_P = 15
public let AN_CL = 17
public let AN_AR = 18
public let AN_BR = 35
public let AN_KR = 36
public let AN_I = 53
public let AN_XE = 54
public let AN_HG = 80
public let AN_AT = 85
public let AN_RN = 86
public let AN_FR = 87

public let AN_NH3 = 900
public let AN_H2O = 901
public let AN_CO2 = 902
public let AN_O3 = 903
public let AN_CH4 = 904
public let AN_CH3CH2OH = 905

/*    The following defines are used in the kothari_radius function in    */
/*    file enviro.c.                                                        */
public let A1_20 = 6.485e12 /* All units are in cgs system.     */
public let A2_20 = 4.0032e-8 /*     ie: cm, g, dynes, etc.         */
public let BETA_20 = 5.71e12

public let JIMS_FUDGE = 1.004

/*     The following defines are used in determining the fraction of a planet     */
/*    covered with clouds in function cloud_fraction in file enviro.c.         */
public let Q1_36 = 1.258e19 /* grams    */
public let Q2_36 = 0.0698 /* 1/Kelvin */
