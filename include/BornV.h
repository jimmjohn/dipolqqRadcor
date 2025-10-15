#pragma once
#include <vector>
#include <complex>
#include <string>
#include <stdexcept>
#include <cstddef>
#include "PhysicsConstants.h"

class BornV {
public:
    // ---- sizes (match your 'BornV.h' params) ----
    // energy grid sizes (each loop is 0..poinX)
    int m_poin1 = 100;                        // LEP1 LOG(SQRT(s)) spacing
    int m_poin2 = 120;                        // Z pole sqrt(s)spacing
    int m_poin3 = 145;                        // LEP2 interval sqrt(s) spacing
    int m_poin4 = 80;                         // NLC range sqrt(s) spacing

    // cos(theta) grid sizes (each loop is 0..poThX)
    int m_poTh2 = 14;                       // Z pole Cost(theta) spacing
    int m_poTh3 = 30;                       // LEP2 Cost(theta) spacing
    int m_poTh4 = 14;                       // NLC Cost(theta) spacing

    int m_poinG = 7;                        // number of EW form factors (complex)
    int m_poinQ = 4;                        // number of QCD correction entries (real)

    // ---- physics parameters read/checked from file header ----
    double m_MZ = EWParameters::Mz;                        // Z mass (expected, set by you before reading; also read from file for check)

    //Set these before proceeding - set from the file header of table you use
    double m_amh   = 0.0;                   // Higgs mass (expected)
    double m_amtop = 0.0;                   // Top mass (expected)
    double m_swsq  = 0.0;                   // sin^2(thetaW)  (read from file)
    double m_gammz = 0.0;                   // Gamma_Z        (read from file)
    double m_MW    = 0.0;                   // W mass         (read from file)
    double m_GammW = 0.0;                   // Gamma_W        (read from file)

    // ---- interpolation windows (set by you) ----
    double m_WdelZ      = 60.0;             // Z range (amz + m_WdelZ)
    double m_WminLEP1   = 0.01;             // LEP1 basic range (m_WminLEP1,m_WmaxLEP1)
    double m_WmaxLEP1   = 95.0;             // LEP1 basic range (m_WminLEP1,m_WmaxLEP1)
    double m_WmaxLEP2   = 240.001;          // LEP2 interval (m_WmaxLEP1,m_WmaxLEP2)
    double m_WmaxNLC    = 1040.001;         // NLC range (m_WmaxLEP2,m_WmaxNLC)

    // ---- outputs of interpolation ----
    std::vector<std::complex<double>> m_GSW;  // size m_poinG (complex)
    std::vector<double>               m_QCDcorR; // size m_poinQ (real)
    double                            m_QCDcor = 0.0; // legacy: m_QCDcorR[0] - 1

    // ---- EW/QCD lookup tables (flattened) ----
    // basic (LEP1 outside Z)
    std::vector<std::complex<double>> m_cyy; // (m_poin1+1, m_poinG, 16)
    std::vector<double>               m_syy; // (m_poin1+1, m_poinQ, 16)
    // Z pole
    std::vector<std::complex<double>> m_czz; // (m_poin2+1, m_poTh2+1, m_poinG, 16)
    std::vector<double>               m_szz; // (m_poin2+1, m_poinQ, 16)
    // LEP2
    std::vector<std::complex<double>> m_ctt; // (m_poin3+1, m_poTh3+1, m_poinG, 16)
    std::vector<double>               m_stt; // (m_poin3+1, m_poinQ, 16)
    // NLC
    std::vector<std::complex<double>> m_clc; // (m_poin4+1, m_poTh4+1, m_poinG, 16)
    std::vector<double>               m_slc; // (m_poin4+1, m_poinQ, 16)

    // ---- lifecycle ----
    BornV();

    // set expected masses (for file consistency check) and region edges
    void set_expected_masses(double MZ, double amh, double amtop) {
        m_MZ = MZ; m_amh = amh; m_amtop = amtop;
    }

    // ---- I/O + interpolation ----
    void BornV_ReadAll(const std::string& dizetDir);
    void BornV_ReadFile(const std::string& DiskFile, int KFfin);
    void BornV_InterpoGSW(int KFf, double svar, double CosThe);


private:
    // index helpers (zero-based everywhere)
    inline std::size_t idx_cyy(int i, int k, int kff) const {
        return (static_cast<std::size_t>(i) * m_poinG + k) * 16u + kff;
    }
    inline std::size_t idx_syy(int i, int k, int kff) const {
        return (static_cast<std::size_t>(i) * m_poinQ + k) * 16u + kff;
    }
    inline std::size_t idx_czz(int i, int j, int k, int kff) const {
        const std::size_t A = static_cast<std::size_t>(m_poTh2 + 1);
        const std::size_t B = static_cast<std::size_t>(m_poinG);
        return (((static_cast<std::size_t>(i) * A) + j) * B + k) * 16u + kff;
    }
    inline std::size_t idx_szz(int i, int k, int kff) const {
        return (static_cast<std::size_t>(i) * m_poinQ + k) * 16u + kff;
    }
    inline std::size_t idx_ctt(int i, int j, int k, int kff) const {
        const std::size_t A = static_cast<std::size_t>(m_poTh3 + 1);
        const std::size_t B = static_cast<std::size_t>(m_poinG);
        return (((static_cast<std::size_t>(i) * A) + j) * B + k) * 16u + kff;
    }
    inline std::size_t idx_stt(int i, int k, int kff) const {
        return (static_cast<std::size_t>(i) * m_poinQ + k) * 16u + kff;
    }
    inline std::size_t idx_clc(int i, int j, int k, int kff) const {
        const std::size_t A = static_cast<std::size_t>(m_poTh4 + 1);
        const std::size_t B = static_cast<std::size_t>(m_poinG);
        return (((static_cast<std::size_t>(i) * A) + j) * B + k) * 16u + kff;
    }
    inline std::size_t idx_slc(int i, int k, int kff) const {
        return (static_cast<std::size_t>(i) * m_poinQ + k) * 16u + kff;
    }

    void allocate_all_();
};

