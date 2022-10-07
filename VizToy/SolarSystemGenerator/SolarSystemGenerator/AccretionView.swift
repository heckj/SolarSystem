//
//  AccretionView.swift
//  SolarSystemGenerator
//
//  Created by Joseph Heck on 10/4/22.
//

import SolarSystem
import SwiftUI

struct AccretionView: View {
    func dustSymbol(_ d: Dust) -> String {
        if d.dust_present, d.gas_present {
            return "*"
        } else if d.dust_present {
            return "+"
        } else if d.gas_present {
            return "."
        } else {
            return "_"
        }
    }

    let tight: FloatingPointFormatStyle<Double> = .number.precision(.integerAndFractionLength(integerLimits: 1..., fractionLimits: 0 ... 3))

    @ObservedObject var accretionDisk: AccretionModel
    @State private var msgStrings: [String] = []

    var body: some View {
        VStack {
            Text("\(accretionDisk.solar_masses) \u{2609} Solar Masses")
            if accretionDisk.state.dust_left {
                Button {
                    msgStrings = []
                    accretionDisk.adv()
                } label: {
                    Image(systemName: "play")
                }
                ForEach(msgStrings.suffix(3), id: \.self) { str in
                    Text(str)
                }
            } else {
                Text("accretion complete")
            }
            AccretionStateView(accretionState: accretionDisk.state)
        }
        .onReceive(accretionDisk.msgs) { msgStr in
            msgStrings.append(msgStr)
        }
    }

    init(model: AccretionModel) {
        accretionDisk = model
    }
}

struct AccretionView_Previews: PreviewProvider {
    static var previews: some View {
        AccretionView(model: AccretionModel(mass: 1.1)).frame(width: 300, height: 300)
    }
}
