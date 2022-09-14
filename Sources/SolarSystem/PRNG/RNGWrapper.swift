//
//  RNGWrapper.swift
//
//
//  Created by Joseph Heck on 12/15/21.
//

import Foundation

/// A class that provides probabilistic functions based on the seedable psuedo-random number generator used to create it.
///
/// ## Topics
///
/// ### Creating a Random Number Generator Wrapper
///
/// - ``RNGWrapper/init(_:)``
///
/// ### Inspecting a Random Number Generator Wrapper
///
/// - ``RNGWrapper/seed``
/// - ``RNGWrapper/position``
///
/// ### Resetting the Seed of a Random Number Generator Wrapper
///
/// - ``RNGWrapper/resetRNG(seed:)``
///
public final class RNGWrapper<PRNG> where PRNG: SeededRandomNumberGenerator {
    private var _prng: PRNG
    #if DEBUG
        var _invokeCount: UInt64 = 0
    #endif
    // access to the underlying PRNG state

    public var seed: UInt64 {
        _prng.seed
    }

    public var position: UInt64 {
        _prng.position
    }

    public func resetRNG(seed: UInt64) {
        _prng = PRNG(seed: seed)
        #if DEBUG
            _invokeCount = 0
        #endif
    }

    /// Creates a new random number generator wrapper class with the random number generator you provide.
    /// - Parameter prng: A random number generator.
    public init(_ prng: PRNG) {
        _prng = prng
    }

    /// Returns a random double value within the range you provide.
    /// - Parameter range: The range of possible values for the double.
    func random_number(in range: ClosedRange<Double>) -> Double {
        #if DEBUG
            _invokeCount += 1
        #endif
        return Double.random(in: range, using: &_prng)
    }

    /// Returns a random double value within the range you provide.
    /// - Parameter range: The range of possible values for the double.
    func about(_ value: Double, variation: Double) -> Double {
        #if DEBUG
            _invokeCount += 1
        #endif
        let adjustment = Double.random(in: -variation...variation, using: &_prng)
        return value + adjustment
    }

    func random_eccentricity() -> Double {
        let e = 1.0 - pow(random_number(in: 0.0...1.0), ECCENTRICITY_COEFF);
        if e > 0.99 {
            return 0.99
        }
        return e
    }
    
}
