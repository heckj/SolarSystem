//
//  AccretionModel.swift
//  SolarSystemGenerator
//
//  Created by Joseph Heck on 10/5/22.
//

import Foundation
import SolarSystem
public class AccretionModel: ObservableObject {
    
    var disk: AccretionDisk
    let solar_masses: Double
    @Published var state: AccretionState
    
    public init(mass: Double) {
        solar_masses = mass
        disk = AccretionDisk(prng: RNGWrapper(Xoshiro(seed: 23)), inner_limit_of_dust: 0.0, outer_limit_of_dust: 0.0, stellar_mass_ratio: mass, stellar_luminosity_ratio: luminosity(mass_ratio: mass))
        state = disk.currentState()
    }
    
    public func adv() {
        disk.advance()
        state = disk.currentState()
    }
}
