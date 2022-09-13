// swift-tools-version: 5.6

import PackageDescription

let package = Package(
    name: "SolarSystem",
    products: [
        .library(
            name: "SolarSystem",
            targets: ["SolarSystem"]
        ),
    ],
    dependencies: [
        .package(url: "https://github.com/apple/swift-docc-plugin", from: "1.0.0"),
    ],
    targets: [
        .target(
            name: "SolarSystem",
            dependencies: []
        ),
        .testTarget(
            name: "SolarSystemTests",
            dependencies: ["SolarSystem"]
        ),
    ]
)
