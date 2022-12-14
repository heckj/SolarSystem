//
//  ContentView.swift
//  SolarSystemGenerator
//
//  Created by Joseph Heck on 9/27/22.
//

import SwiftUI

struct ContentView: View {
    var body: some View {
        VStack {
            Image(systemName: "globe")
                .imageScale(.large)
                .foregroundColor(.accentColor)
            Text("Hello, world!")
            AccretionView(model: AccretionModel(mass: 1.1))
        }
        .padding()
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
