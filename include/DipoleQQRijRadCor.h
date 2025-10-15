#ifndef DIPOLE_QQRIJ_RADCOR_H
#define DIPOLE_QQRIJ_RADCOR_H

#include <array>
#include <complex>

class DipoleQQRijRadCor {
public:
    DipoleQQRijRadCor(int iqed_enabled);

    // Calculate the spin-correlation matrix R given energy, angle, and dipole form factors
    std::array<std::array<double,4>,4> calculate(
        double energy,
        double theta,
        double ReA, double ImA,
        double ReB, double ImB,
        double ReX, double ImX,
        double ReY, double ImY,
        double channel
    );

private:
    int iqed_;

    // Physical constants
    const double pi_ = 3.141592653589793238;
    const double alpha_ = 7.2973525693e-3;
    const double Mz_ = 91.1876;     // Z-boson mass (GeV)
    const double Gz_ = 2.4952;      // Z-boson width (GeV)
    double v0_ = -0.03783;    // Vector coupling for tau
    double a0_ = -0.50123;    // Axial coupling for tau
    const double mtau_ = 1.77686;   // Tau mass (GeV)
    const double GF_ = 1.1663788e-5; // Fermi constant (GeV^{-2})
    const double GMu_ = 1.66389e-5; // Muon decay constant (GeV^{-2})
    const double Mw_ = 80.353;     // W-boson mass (GeV)
};

#endif
