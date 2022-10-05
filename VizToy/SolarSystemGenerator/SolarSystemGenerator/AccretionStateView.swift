//
//  AccretionStateView.swift
//  SolarSystemGenerator
//
//  Created by Joseph Heck on 10/5/22.
//

import SwiftUI
import SolarSystem
import SwiftVizScale

struct AccretionStateView: View {
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
    func dustColor(_ d:Dust) -> Color {
        if d.dust_present && d.gas_present {
            return .brown
        } else if d.dust_present {
            return .red
        } else if d.gas_present {
            return .yellow
        } else {
            return .clear
        }
    }

    
    let tight: FloatingPointFormatStyle<Double> = .number.precision(.integerAndFractionLength(integerLimits: 1..., fractionLimits: 0...3))

    let accretionState: AccretionState
    let scale: SwiftVizScale.ContinuousScale<CGFloat>
    let planetScale: SwiftVizScale.ContinuousScale<CGFloat>
    
    var body: some View {
        VStack {
            ForEach(accretionState.dustlanes, id: \.inner_edge) { dustlane in
                Text("\(dustlane.inner_edge.formatted(tight)) - \(dustlane.outer_edge.formatted(tight)) \(dustSymbol(dustlane))")
//                Text("log: \(scale.scale(dustlane.inner_edge, from: 5.0, to: 90.0) ?? 1.1) to \(scale.scale(dustlane.outer_edge, from: 5.0, to: 90.0) ?? 2.2)")
            }
            Canvas { ctx, size in
//                ctx.draw(Text("\(size.height) \(size.width)"), in: CGRect(origin: .zero, size: size))
//
                let y = size.height/2
                let insetWidth = size.width * 0.9
                let leftInsetMargin = size.width * 0.05
                let right = leftInsetMargin + insetWidth
                //ctx.draw(Text("\(leftInsetMargin) \(right)"), in: CGRect(origin: .zero, size: size))
                
                for d in accretionState.dustlanes {
//                    let basePath = Path { p in
//                        p.move(to: CGPoint(x: d.inner_edge, y: y))
//                        p.addLine(to: CGPoint(x: d.outer_edge, y: y))
//                    }
//                    ctx.stroke(basePath, with: .color(dustColor(d)))
                    
                    let dPath = Path { p in
                        p.move(to: CGPoint(x: scale.scale(d.inner_edge, from: leftInsetMargin, to: right) ?? 0, y: y))
                        p.addLine(to: CGPoint(x: scale.scale(d.outer_edge, from: leftInsetMargin, to: right) ?? 10, y: y))
                    }
                    ctx.stroke(dPath, with: GraphicsContext.Shading.color(dustColor(d)), lineWidth: 4.0)
                }
                for plnt in accretionState.planets {
                    let planetCircle = Path { p in
                        let xy: CGFloat = planetScale.scale(plnt.mass) ?? 1.0
                        if let xValue = scale.scale(plnt.a, from: leftInsetMargin, to: right) {
                            p.addEllipse(in: CGRect(x: xValue - xy/2.0, y: y-xy/2.0, width: xy, height: xy))
                        }
                    }
                    ctx.stroke(planetCircle, with: .color(.blue))
                }
            }.frame(height: 50)
        }.frame(width: 300, height: 500)
    }
    
    public init(accretionState: AccretionState) {
        self.accretionState = accretionState
//        let min: Double = accretionState.dustlanes.reduce(into: 0.01) { partialResult, dustlane in
//            if dustlane.inner_edge < partialResult {
//                partialResult = dustlane.inner_edge
//            }
//        }
        let max: Double = accretionState.dustlanes.reduce(into: 0.0) { partialResult, dustlane in
            if dustlane.outer_edge > partialResult {
                partialResult = dustlane.outer_edge
            }
        }
        let maxPlanetMass: Double = accretionState.planets.reduce(into: 0.0) { partialResult, planet in
            if planet.mass > partialResult {
                partialResult = planet.mass
            }
        }
        scale = ContinuousScale(lower: 0.01, higher: max, type: .log, transform: .clamp)
        planetScale = ContinuousScale(lower: 0.1, higher: maxPlanetMass, type: .linear, rangeLower: 1, rangeHigher: 15)
    }
}

struct AccretionStateView_Previews: PreviewProvider {
    static var previews: some View {
        AccretionStateView(accretionState: AccretionState.example)
    }
}
