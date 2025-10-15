#pragma once
#include <array>

struct EventInitilizers {
    bool pvTaum_found=false, pvTaul_found=false;
    int heTaum=0, heTaul=0;
    float wt=0.0f, wt_approx=0.0f;
    std::array<double,4> m_HvCloneTaum, m_HvCloneTaul;
    std::array<double,4> P1, P2;
    std::array<double,4> PIPM, PIPL, PIZM, PIZL;
    double beamEnergy=0.0;
};