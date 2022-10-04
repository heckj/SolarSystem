//
//  SwiftUIView.swift
//  SolarSystemGenerator
//
//  Created by Joseph Heck on 10/3/22.
//

import SwiftUI
import SolarSystem

struct StarMassSelectorView: View {
    let min_mass = 0.4
    let inc_mass = 0.05
    let max_mass = 2.35

    @State private var mass: Double = 0.4
    
    var body: some View {
        VStack {
            Text("\(mass) Solar Masses")
            Text("luminosity: \(luminosity(mass_ratio: mass))")
            Slider(value: $mass, in: min_mass...max_mass) {
                Text("Mass")
            } minimumValueLabel: {
                Text(min_mass.formatted())
            } maximumValueLabel: {
                Text(max_mass.formatted())
            }
        }
    }
}

struct StarMassSelectorView_Previews: PreviewProvider {
    static var previews: some View {
        StarMassSelectorView()
    }
}
