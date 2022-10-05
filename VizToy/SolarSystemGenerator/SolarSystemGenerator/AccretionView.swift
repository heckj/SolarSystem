//
//  AccretionView.swift
//  SolarSystemGenerator
//
//  Created by Joseph Heck on 10/4/22.
//

import SwiftUI
import SolarSystem


struct AccretionView: View {
    func dustSymbol(_ d:Dust) -> String {
        if d.dust_present && d.gas_present {
            return "*"
        } else if d.dust_present {
            return "+"
        } else if d.gas_present {
            return "."
        } else {
            return "_"
        }
        
    }
    
    let tight: FloatingPointFormatStyle<Double> = .number.precision(.integerAndFractionLength(integerLimits: 1..., fractionLimits: 0...3))
    
    let accretionDisk: AccretionDisk
    @State var accreteState: AccretionState? = nil
    var body: some View {
        VStack {
            Text("\(accretionDisk.stellar_mass_ratio) \u{2609} Solar Masses")
            if let accreteState = accreteState {
                ForEach(accreteState.dustlanes, id: \.inner_edge) { dustlane in
                    Text("\(dustlane.inner_edge.formatted(tight)) - \(dustlane.outer_edge.formatted(tight)) \(dustSymbol(dustlane))")
                }
            }
        }
        .onAppear() {
            accreteState = accretionDisk.currentState()
        }
    }
    
    init(mass: Double) {
        self.accretionDisk = AccretionDisk(prng: RNGWrapper(Xoshiro(seed: 23)), inner_limit_of_dust: 0.0, outer_limit_of_dust: 0.0, stellar_mass_ratio: mass, stellar_luminosity_ratio: luminosity(mass_ratio: mass))
    }
}

struct AccretionView_Previews: PreviewProvider {
    static var previews: some View {
        AccretionView(mass: 1.0).frame(width: 300, height: 300)
    }
}
