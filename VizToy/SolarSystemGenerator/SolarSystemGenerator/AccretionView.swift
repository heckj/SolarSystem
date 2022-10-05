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
    
    @ObservedObject var accretionDisk: AccretionModel
    
    var body: some View {
        VStack {
            Text("\(accretionDisk.solar_masses) \u{2609} Solar Masses")
            Button {
                print("hi")
                accretionDisk.adv()
            } label: {
                Image(systemName: "play")
            }
            AccretionStateView(accretionState: accretionDisk.state)
        }
    }
    
    init(model: AccretionModel) {
        self.accretionDisk = model
    }
}

struct AccretionView_Previews: PreviewProvider {
    static var previews: some View {
        AccretionView(model: AccretionModel(mass: 1.1)).frame(width: 300, height: 300)
    }
}
