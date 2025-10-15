#pragma once

namespace Physics {
constexpr double m_tau = 1.77686; // Tau mass in GeV
constexpr double m_mu = 0.105658; // Muon mass in GeV
constexpr double m_pi = 0.139570; // Charged pion mass in GeV
constexpr double m_rho = 0.77526; // Rho meson mass in GeV
}

namespace EWParameters {
constexpr double Mz = 91.188;    // Z-boson mass in GeV
constexpr double Gz = 2.4955;     // Z-boson width in GeV
constexpr double Mw = 80.3692;     // W-boson mass in GeV
// constexpr double sw2 = 0.22336;   // sin^2(theta_W)
// constexpr double cw2 = 1.0 - sw2; // cos^2(theta_W)
constexpr double GF = 1.1663788e-5; // Fermi constant in GeV^{-2}
}

namespace PDG {
    enum class ParticleID {
        electron = 11,
        positron = -11,
        muon     = 13,
        tau      = 15
    };
}