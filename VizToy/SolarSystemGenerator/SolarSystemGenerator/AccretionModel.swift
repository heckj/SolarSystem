//
//  AccretionModel.swift
//  SolarSystemGenerator
//
//  Created by Joseph Heck on 10/5/22.
//

import Foundation
import Combine
import SolarSystem

public class AccretionModel: ObservableObject {
    
    var disk: AccretionDisk
    let solar_masses: Double
    @Published var state: AccretionState
    
    public var msgs: AnyPublisher<String, Never>
    
    public init(mass: Double, seed: UInt64 = 23) {
        solar_masses = mass
        disk = AccretionDisk(prng: RNGWrapper(Xoshiro(seed: seed)), inner_limit_of_dust: 0.0, outer_limit_of_dust: 0.0, stellar_mass_ratio: mass, stellar_luminosity_ratio: luminosity(mass_ratio: mass))
        state = disk.currentState()
        msgs = disk.msgs.eraseToAnyPublisher()
    }
    
    public func adv() {
        disk.advance()
        state = disk.currentState()
    }
}
