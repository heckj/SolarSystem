//
//  ChamTable.swift
//

import Foundation

// From Keris

final class ChemTableEntry {
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

    //   An   sym   HTML symbol                      name                 Aw      melt    boil    dens       ABUNDe       ABUNDs         Rea    Max inspired pp
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

let gases: [ChemTableEntry] = [
    ChemTableEntry(num: AN_H, symbol: "H", html_symbol: "H<SUB><SMALL>2</SMALL></SUB>", name: "Hydrogen",
                   weight: 1.0079, melt: 14.06, boil: 20.40, density: 8.99e-05,
                   abunde: 0.00125893, abunds: 27925.4, reactivity: 1, max_ipp: 0.0),
    ChemTableEntry(num: AN_HE, symbol: "He", html_symbol: "He", name: "Helium",
                   weight: 4.0026, melt: 3.46, boil: 4.20, density: 0.0001787,
                   abunde: 7.94328e-09, abunds: 2722.7, reactivity: 0, max_ipp: MAX_HE_IPP),
    ChemTableEntry(num: AN_N, symbol: "N", html_symbol: "N<SUB><SMALL>2</SMALL></SUB>", name: "Nitrogen",
                   weight: 14.0067, melt: 63.34, boil: 77.40, density: 0.0012506,
                   abunde: 1.99526e-05, abunds: 3.13329, reactivity: 0, max_ipp: MAX_N2_IPP),
    ChemTableEntry(num: AN_O, symbol: "O", html_symbol: "O<SUB><SMALL>2</SMALL></SUB>", name: "Oxygen",
                   weight: 15.9994, melt: 54.80, boil: 90.20, density: 0.001429,
                   abunde: 0.501187, abunds: 23.8232, reactivity: 10, max_ipp: MAX_O2_IPP),
    ChemTableEntry(num: AN_NE, symbol: "Ne", html_symbol: "Ne", name: "Neon",
                   weight: 20.1700, melt: 24.53, boil: 27.10, density: 0.0009,
                   abunde: 5.01187e-09, abunds: 3.4435e-5, reactivity: 0, max_ipp: MAX_NE_IPP),
    ChemTableEntry(num: AN_AR, symbol: "Ar", html_symbol: "Ar", name: "Argon",
                   weight: 39.9480, melt: 84.00, boil: 87.30, density: 0.0017824,
                   abunde: 3.16228e-06, abunds: 0.100925, reactivity: 0, max_ipp: MAX_AR_IPP),
    ChemTableEntry(num: AN_KR, symbol: "Kr", html_symbol: "Kr", name: "Krypton",
                   weight: 83.8000, melt: 116.60, boil: 119.70, density: 0.003708,
                   abunde: 1e-10, abunds: 4.4978e-05, reactivity: 0, max_ipp: MAX_KR_IPP),
    ChemTableEntry(num: AN_XE, symbol: "Xe", html_symbol: "Xe", name: "Xenon",
                   weight: 131.3000, melt: 161.30, boil: 165.00, density: 0.00588,
                   abunde: 3.16228e-11, abunds: 4.69894e-06, reactivity: 0, max_ipp: MAX_XE_IPP),
    ChemTableEntry(num: AN_NH3, symbol: "NH3", html_symbol: "NH<SUB><SMALL>3</SMALL></SUB>", name: "Ammonia",
                   weight: 17.0000, melt: 195.46, boil: 239.66, density: 0.001,
                   abunde: 0.002, abunds: 0.0001, reactivity: 1, max_ipp: MAX_NH3_IPP),
    ChemTableEntry(num: AN_H2O, symbol: "H20", html_symbol: "H<SUB><SMALL>2</SMALL></SUB>O", name: "Water",
                   weight: 18.0000, melt: 273.16, boil: 373.16, density: 1.000,
                   abunde: 0.03, abunds: 0.001, reactivity: 0, max_ipp: 0.0),
    ChemTableEntry(num: AN_CO2, symbol: "CO2", html_symbol: "CO<SUB><SMALL>2</SMALL></SUB>", name: "CarbonDioxide",
                   weight: 44.0000, melt: 194.66, boil: 194.66, density: 0.001,
                   abunde: 0.01, abunds: 0.0005, reactivity: 0, max_ipp: MAX_CO2_IPP),
    ChemTableEntry(num: AN_O3, symbol: "O3", html_symbol: "O<SUB><SMALL>3</SMALL></SUB>", name: "Ozone",
                   weight: 48.0000, melt: 80.16, boil: 161.16, density: 0.001,
                   abunde: 0.001, abunds: 0.000001, reactivity: 2, max_ipp: MAX_O3_IPP),
    ChemTableEntry(num: AN_CH4, symbol: "CH4", html_symbol: "CH<SUB><SMALL>4</SMALL></SUB>", name: "Methane",
                   weight: 16.0000, melt: 90.16, boil: 109.16, density: 0.010,
                   abunde: 0.005, abunds: 0.0001, reactivity: 1, max_ipp: MAX_CH4_IPP),
    ChemTableEntry(num: AN_F, symbol: "F", html_symbol: "F", name: "Fluorine",
                   weight: 18.9984, melt: 53.58, boil: 85.10, density: 0.001696,
                   abunde: 0.000630957, abunds: 0.000843335, reactivity: 50, max_ipp: MAX_F_IPP),
    ChemTableEntry(num: AN_CL, symbol: "Cl", html_symbol: "Cl", name: "Chlorine",
                   weight: 35.4530, melt: 172.22, boil: 239.20, density: 0.003214,
                   abunde: 0.000125893, abunds: 0.005236, reactivity: 40, max_ipp: MAX_CL_IPP),
    ChemTableEntry(num: AN_CH3CH2OH, symbol: "CH3CH2OH", html_symbol: "CH3CH2OH", name: "Ethanol",
                   weight: 46.0000, melt: 159.06, boil: 351.66, density: 0.895,
                   abunde: 0.001, abunds: 0.001, reactivity: 0, max_ipp: 0),
]
