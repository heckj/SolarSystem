@testable import SolarSystem
import XCTest

final class SolarSystemTests: XCTestCase {
    let mass = 1.2
    let seed: UInt64 = 2346
    var disk: AccretionDisk?
    
    let seed_sequence: [(Double,Double)] = [
        (1.2885, 0.0401084),    // 0
        (34.6525, 0.000524515), // 1
        (46.3695, 0.0294519),   // 2
        (19.6841, 0.0365445),   // 3
        (5.99764, 0.125626),    // 4
        (49.6731, 0.000530321), // 5
        (45.1105, 0.223403),    // 6
        (9.30144, 0.0464744),   // 7
        (51.7443, 0.173246),    // 8
        (0.43533, 0.173706),    // 9
        (2.10938, 0.0155614),   // 10
        (3.45524, 0.157992),    // 11
        (0.776317, 0.039044)    // 12
    ]
    
    override func setUp() async throws {
        self.disk = AccretionDisk(prng: RNGWrapper(Xoshiro(seed: seed)),
                                  inner_limit_of_dust: 0.0,
                                  outer_limit_of_dust: 0.0,
                                  stellar_mass_ratio: mass,
                                  stellar_luminosity_ratio: luminosity(mass_ratio: mass))
    }

    func testInitialParameters() throws {
        guard let disk = disk else {
            XCTFail("Accretion disk failed to set up")
            return
        }
        let state = disk.currentState()
        XCTAssertEqual(state.dust_left, true)
        XCTAssertEqual(state.planets.count, 0)
        XCTAssertEqual(state.dustlanes.count, 1)
        XCTAssertEqual(state.dustlanes[0].inner_edge, 0.0)
        XCTAssertEqual(state.dustlanes[0].outer_edge, 212.53171, accuracy: 0.0001)
        XCTAssertEqual(disk.dust_density, 9.49948e-06, accuracy: 0.0001)
    }
    
    func testAccretionStep() throws {
        guard let disk = disk else {
            XCTFail("Accretion disk failed to set up")
            return
        }
        var diskCopy = disk
        diskCopy.advance(distance: 1.2885, eccentricity: 0.0401084)
        let state = diskCopy.currentState()
        XCTAssertEqual(state.dustlanes.count, 3)
        XCTAssertEqual(state.planets.count, 1)
        guard let planet = state.planets.first else {
            XCTFail("Accretion disk failed to generate intial planet")
            return
        }
        XCTAssertEqual(planet.a, 1.29, accuracy: 0.01)
        XCTAssertEqual(planet.e, 0.04, accuracy: 0.01)
        XCTAssertEqual(planet.mass, 1.1528e-5, accuracy: 1.0e-6)
        XCTAssertEqual(planet.dust_mass, 9.263e-06, accuracy: 1.0e-8)
        XCTAssertEqual(planet.gas_mass, 2.264e-06, accuracy: 1.0e-8)
    }
//    Injecting protoplanet (mass=1e-15) at 1.2885 AU (e=0.0401084).
//    Dust density at 1.2885 AU is 9.49948e-06.
//    accrete_dust
//      . mass of 1e-15 collecting dust between 1.0305 and 1.67552
//       . new_mass 2.27328e-08 = vol(0.00239306) * density(9.49948e-06)
//      . mass of 1e-15 collecting dust between 1.0305 and 1.67552
//       . accumulated_mass 2.27328e-08 new_mass(2.27328e-08) from dust(3.04713e-314) gas(3.04713e-314)
//    Accreted 3.08291EM dust and 0.753819EM gas.
//    Added planet 1.29 AU (0.040) (3.84EM = 3.08EMd + 0.75EMg [2.452EM])

    func testMultipleSteps() throws {
        guard let disk = disk else {
            XCTFail("Accretion disk failed to set up")
            return
        }
        var diskCopy = disk
        
        let state0 = diskCopy.currentState()
        XCTAssertEqual(state0.planets.count, 0)

        diskCopy.advance(distance: seed_sequence[0].0, eccentricity: seed_sequence[0].1)
        XCTAssertEqual(diskCopy.currentState().planets.count, 1)

        diskCopy.advance(distance: seed_sequence[1].0, eccentricity: seed_sequence[1].1)
        XCTAssertEqual(diskCopy.currentState().planets.count, 2)
        diskCopy.advance(distance: seed_sequence[2].0, eccentricity: seed_sequence[2].1)
        XCTAssertEqual(diskCopy.currentState().planets.count, 3)
        diskCopy.advance(distance: seed_sequence[3].0, eccentricity: seed_sequence[3].1)
        XCTAssertEqual(diskCopy.currentState().planets.count, 4)
//        diskCopy.advance(distance: seed_sequence[4].0, eccentricity: seed_sequence[4].1)
//        XCTAssertEqual(diskCopy.currentState().planets.count, 5)
//        diskCopy.advance(distance: seed_sequence[5].0, eccentricity: seed_sequence[5].1)
//        XCTAssertEqual(diskCopy.currentState().planets.count, 6)
//        diskCopy.advance(distance: seed_sequence[6].0, eccentricity: seed_sequence[6].1)
//        XCTAssertEqual(diskCopy.currentState().planets.count, 7)
//        diskCopy.advance(distance: seed_sequence[7].0, eccentricity: seed_sequence[7].1)
//        XCTAssertEqual(diskCopy.currentState().planets.count, 8)
//        diskCopy.advance(distance: seed_sequence[8].0, eccentricity: seed_sequence[8].1)
//        XCTAssertEqual(diskCopy.currentState().planets.count, 9)
//        diskCopy.advance(distance: seed_sequence[9].0, eccentricity: seed_sequence[9].1)
//        diskCopy.advance(distance: seed_sequence[10].0, eccentricity: seed_sequence[10].1)
//        diskCopy.advance(distance: seed_sequence[11].0, eccentricity: seed_sequence[11].1)
//        diskCopy.advance(distance: seed_sequence[12].0, eccentricity: seed_sequence[12].1)
//
//        let state = diskCopy.currentState()
//        XCTAssertEqual(state.planets.count, 11)
    }
}
