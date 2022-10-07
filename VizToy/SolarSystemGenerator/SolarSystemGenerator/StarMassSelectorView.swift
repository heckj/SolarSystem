//
//  SwiftUIView.swift
//  SolarSystemGenerator
//
//  Created by Joseph Heck on 10/3/22.
//

import SolarSystem
import SwiftUI

struct StarMassSelectorView: View {
    let min_mass = 0.4
    let inc_mass = 0.05
    let max_mass = 2.35

    @State private var mass: Double = 0.4

    let FPStyle: FloatingPointFormatStyle<Double> = .number.precision(.integerAndFractionLength(integerLimits: 1..., fractionLimits: 0 ... 3))

    var body: some View {
        VStack {
            Text("\(mass.formatted(FPStyle)) \u{2609} Solar Masses")
            Text("luminosity: \(luminosity(mass_ratio: mass).formatted(FPStyle))")
            Text("r_ecosphere: \(sqrt(luminosity(mass_ratio: mass)).formatted(FPStyle))")
            Text("lifetime: \((1.0e10 * (mass / luminosity(mass_ratio: mass))).formatted(.number.precision(.integerAndFractionLength(integerLimits: 1..., fractionLimits: 0 ... 3))))")

            Slider(value: $mass, in: min_mass ... max_mass) {
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
