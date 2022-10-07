// swift-tools-version: 5.6

import PackageDescription

let package = Package(
    name: "SolarSystem",
    platforms: [
        .iOS(.v15),
        .macOS(.v12),
        .tvOS(.v15),
        .watchOS(.v8),
    ],
    products: [
        .library(
            name: "SolarSystem",
            targets: ["SolarSystem"]
        ),
        .executable(name: "Stargen", targets: ["Stargen"]),
    ],
    dependencies: [
        .package(url: "https://github.com/apple/swift-docc-plugin", from: "1.0.0"),
        .package(url: "https://github.com/apple/swift-argument-parser", from: "1.0.0"),
        .package(url: "https://github.com/apple/swift-system", from: "1.0.0"),
        .package(url: "https://github.com/johnsundell/plot.git", from: "0.9.0"),
    ],
    targets: [
        .executableTarget(
            name: "Stargen",
            dependencies: [
                "SolarSystem",
                .product(name: "ArgumentParser", package: "swift-argument-parser"),
                .product(name: "SystemPackage", package: "swift-system"),
            ]
        ),
        .target(
            name: "SolarSystem",
            dependencies: [
                .product(name: "Plot", package: "plot"),
            ],
            swiftSettings: [
                //                .unsafeFlags(["-Xfrontend", "-warn-concurrency", "-Xfrontend", "-enable-actor-data-race-checks"])
                .unsafeFlags(["-Xfrontend", "-strict-concurrency=complete"]),
                /*
                    Summation from https://www.donnywals.com/enabling-concurrency-warnings-in-xcode-14/
                    Set `strict-concurrency` to `targeted` to enforce Sendable and actor-isolation checks in your code.
                      This explicitly verifies that `Sendable` constraints are met when you mark one of your types as `Sendable`.
                      This mode is essentially a bit of a hybrid between the behavior that's intended in Swift 6, and the default in Swift 5.7.
                      Use this mode to have a bit of checking on your code that uses Swift concurrency without too many warnings and / or errors in your current codebase.

                    Set `strict-concurrency` to `complete` to get the full suite of concurrency constraints, essentially as they will work in Swift 6.
                    */
            ]
        ),
        .testTarget(
            name: "SolarSystemTests",
            dependencies: ["SolarSystem"]
        ),
    ]
)
