#include "BornV.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cmath>

static bool read_fortran_double(std::istream& in, double& out) {
    std::string tok;
    if (!(in >> tok)) return false;
    for (char& c : tok) if (c == 'D' || c == 'd') c = 'E'; // D exponent -> E
    out = std::stod(tok);
    return true;
}

// ---- all-files loader (BornV_ReadAll) ----
void BornV::BornV_ReadAll(const std::string& dizetDir) {
    // PDG codes used in your Fortran
    const int KFdown=1, KFup=2, KFstran=3, KFcharm=4, KFbotom=5, KFtop=6;
    const int KFel=11, KFelnu=12, KFmu=13, KFmunu=14, KFtau=15, KFtaunu=16;

    std::cout << "BornV: Loading Dizet EW tables from directory: " << dizetDir << std::endl;
    auto path = [&](const char* name){ return dizetDir + "/" + std::string(name); };
    std::cout << " BornV: Reading file: " << path("table.mu") << std::endl;
    BornV_ReadFile(path("table.mu"),  KFmu);
    std::cout << " BornV: Reading file: " << path("table.up") << std::endl;
    BornV_ReadFile(path("table.up"),    KFup);
    std::cout << " BornV: Reading file: " << path("table.down") << std::endl;
    BornV_ReadFile(path("table.down"),  KFstran);
    std::cout << " BornV: Reading file: " << path("table.up") << std::endl;
    BornV_ReadFile(path("table.up"),    KFcharm);
    //BornV_ReadFile(path("table.botom"), KFbotom);
    std::cout << " BornV: Reading file: " << path("table.down") << std::endl;
    BornV_ReadFile(path("table.down"),    KFdown);
    //BornV_ReadFile(path("table.tau"),   KFtau);
    //BornV_ReadFile(path("table.nue"),   KFelnu);
    //BornV_ReadFile(path("table.numu"),  KFmunu);
    //BornV_ReadFile(path("table.nutau"), KFtaunu);

    // (No electron file in the original Fortran sequence.)
}


void BornV::BornV_ReadFile(const std::string& DiskFile, int KFfin) {
    int KFf = std::abs(KFfin);
    if (KFf < 1 || KFf > 16) {
        throw std::runtime_error("BornV_ReadFile: wrong KFf = " + std::to_string(KFfin));
    }
    const int kf = KFf - 1; // zero-based flavour
    std::ifstream in(DiskFile);

    if (!in) throw std::runtime_error("Cannot open file: " + DiskFile);
    // header: amz, amh, amtop, m_swsq, m_gammz, m_MW, m_GammW
    if (!read_fortran_double(in, m_MZ) ||
        !read_fortran_double(in, m_amh) ||
        !read_fortran_double(in, m_amtop) ||
        !read_fortran_double(in, m_swsq) ||
        !read_fortran_double(in, m_gammz) ||
        !read_fortran_double(in, m_MW) ||
        !read_fortran_double(in, m_GammW)) {
        throw std::runtime_error("BornV_ReadFile: header read failed");
    }
    std::cout << std::fixed << std::setprecision(6);

    std::cout << " BornV_ReadFile: read header mz,mh,mtop = "
              << m_MZ << ", " << m_amh << ", " << m_amtop << std::endl;
    std::cout << "                 swsq, gammz, MW, GammW = "
              << m_swsq << ", " << m_gammz << ", "
              << m_MW << ", " << m_GammW << std::endl;

    // // Optional consistency check (matches Fortran STOPs)
    // auto bad = [](double expected, double got){
    //     if (expected == 0.0) return false; // skip check if not preset
    //     return std::abs(1.0 - expected / got) > 1e-5;
    // };
    // if (amz == 0.0 || amh == 0.0 || amtop == 0.0
    //     || bad(m_, amz) || bad(m_amh, amh) || bad(m_amtop, amtop)) {
    //     std::ostringstream oss;
    //     oss << "BornV_ReadFile: CHECK mz,mh,mtop (expected/read) = "
    //         << m_MZ << "/" << amz << "  "
    //         << m_amh << "/" << amh << "  "
    //         << m_amtop << "/" << amtop;
    //     throw std::runtime_error(oss.str());
    // }

    // ---- helpers to read header cards ----
    // 2) Read header line "(a, i4, f11.5 ...)" → chr, n, ww
    auto read_line_header3 = [&](char& chr, int& n, double& ww){
        // (a, i4, f11.5, i4, f11.5) but we only need first 3 here
        in >> std::ws;
        if (!(in >> chr)) throw std::runtime_error("header3: chr");
        if (!(in >> n))   throw std::runtime_error("header3: n");
        if (!read_fortran_double(in, ww)) throw std::runtime_error("header3: ww");
    };
    auto read_line_header4 = [&](char& chr, int& n1, double& ww, int& n2){
        // (a, i4, f11.5, i4, f11.5) but we only need first 3 here
        in >> std::ws;
        if (!(in >> chr)) throw std::runtime_error("header4: chr");
        if (!(in >> n1))  throw std::runtime_error("header4: n1");
        if (!read_fortran_double(in, ww)) throw std::runtime_error("header4: ww");
        if (!(in >> n2))  throw std::runtime_error("header4: n2");
    };
    auto read_line_header5 = [&](char& chr, int& n1, double& ww, int& n2, double& cosi){
        in >> std::ws;
        if (!(in >> chr)) throw std::runtime_error("header5: chr");
        if (!(in >> n1))  throw std::runtime_error("header5: n1");
        if (!read_fortran_double(in, ww)) throw std::runtime_error("header5: ww");
        if (!(in >> n2))  throw std::runtime_error("header5: n2");
        if (!read_fortran_double(in, cosi)) throw std::runtime_error("header5: cosi");
    };

    // 1) Read exactly one g13.7 field (13 chars), tolerating missing '1x' and CR/LF
    auto read_g13_number = [&]() -> double {
        while (in && (in.peek()=='\n' || in.peek()=='\r')) in.get(); // before
        char buf[13];
        in.read(buf, 13);
        if (!in) throw std::runtime_error("EOF while reading g13.7 field");
        std::string s(buf, 13);
        // trim
        auto l = s.find_first_not_of(" \t");
        auto r = s.find_last_not_of(" \t");
        s = (l == std::string::npos) ? std::string() : s.substr(l, r - l + 1);
        for (char& c : s) if (c=='D' || c=='d') c='E';
        double val = s.empty() ? 0.0 : std::stod(s);
        if (in && in.peek()==' ') in.get();                       // optional 1x
        while (in && (in.peek()=='\n' || in.peek()=='\r')) in.get(); // after
        return val;
    };


    char chr; int n, n1, n2; double ww, cosi;
    std::vector<std::complex<double>> rowC;
    std::vector<double> rowR;

    // ---- basic range (cyy, syy) ----
    for (int i = 0; i <= m_poin1; ++i) {
        read_line_header3(chr, n, ww);
        //std::cout << " read_line_header3: chr = " << chr << ", n = " << n << ", ww = " << ww << std::endl;
        for (int k = 0; k < m_poinG; ++k) {
            double re = read_g13_number();
            double im = read_g13_number();
            cyy[i][k][kf] = {re, im};
            //std::cout << " read_g13_number: re = " << re << ", im = " << im << std::endl;
        }
        for (int k = 0; k < m_poinQ; ++k) {
            syy[i][k][kf] = read_g13_number();
            //std::cout << " read_g13_number: syy[" << i << "][" << k << "] = " << syy[i][k][kf] << std::endl;
        }
    }

    // ---- Z pole (czz, szz) ----
    for (int i = 0; i <= m_poin2; ++i) {
        for (int j = 0; j <= m_poTh2; ++j) {
            read_line_header4(chr, n1, ww, n2); // trailing j on the line is ignored by the parser
            //std::cout << " read_line_header4: chr = " << chr << ", n1 = " << n1 << ", ww = " << ww << ", n2 = " << n2 << std::endl;
            for (int k = 0; k < m_poinG; ++k) {
                double re = read_g13_number();
                double im = read_g13_number();
                czz[i][j][k][kf] = {re, im};
                //std::cout << " read_g13_number: re = " << re << ", im = " << im << std::endl;
            }
        }
        for (int k = 0; k < m_poinQ; ++k) {
            szz[i][k][kf] = read_g13_number();
            //std::cout << " read_g13_number: szz[" << i << "][" << k << "] = " << szz[i][k][kf] << std::endl;
        }

    }

    // ---- LEP2 (ctt, stt) ----
    for (int i = 0; i <= m_poin3; ++i) {
        for (int j = 0; j <= m_poTh3; ++j) {
            read_line_header5(chr, n1, ww, n2, cosi);
            //std::cout << " read_line_header5: chr = " << chr << ", n1 = " << n1 << ", ww = " << ww << ", n2 = " << n2 << ", cosi = " << cosi << std::endl;
            for (int k = 0; k < m_poinG; ++k) {
                double re = read_g13_number();
                double im = read_g13_number();
                ctt[i][j][k][kf] = {re, im};
                //std::cout << " read_g13_number: re = " << re << ", im = " << im << std::endl;
            }
        }
        for (int k = 0; k < m_poinQ; ++k) {
            stt[i][k][kf] = read_g13_number();
            //std::cout << " read_g13_number: stt[" << i << "][" << k << "] = " << stt[i][k][kf] << std::endl;
        }
    }

    // ---- NLC (clc, slc) ----
    for (int i = 0; i <= m_poin4; ++i) {
        for (int j = 0; j <= m_poTh4; ++j) {
            read_line_header5(chr, n1, ww, n2, cosi);
            //std::cout << " read_line_header5: chr = " << chr << ", n1 = " << n1 << ", ww = " << ww << ", n2 = " << n2 << ", cosi = " << cosi << std::endl;
            for (int k = 0; k < m_poinG; ++k) {
                double re = read_g13_number();
                double im = read_g13_number();
                clc[i][j][k][kf] = {re, im};
                //std::cout << " read_g13_number: re = " << re << ", im = " << im << std::endl;
            }
        }
        for (int k = 0; k < m_poinQ; ++k) {
            slc[i][k][kf] = read_g13_number();
            //std::cout << " read_g13_number: slc[" << i << "][" << k << "] = " << slc[i][k][kf] << std::endl;
        }
    }
}

void BornV::BornV_InterpoGSW(int KFf, double svar, double CosThe) {
    if (KFf < 1 || KFf > 16) throw std::runtime_error("BornV_InterpoGSW: KFf out of [1,16]");
    const int kf = KFf - 1;

    const double ww = std::sqrt(std::max(0.0, svar));
    const double WminZ = m_MZ - m_WdelZ;
    const double WmaxZ = m_MZ + m_WdelZ;

    auto lerp = [](auto a, auto b, double t){ return a*(1.0 - t) + b*t; };

     // clear outputs
    for (int kk = 0; kk < m_poinG; ++kk) m_GSW[kk] = {0.0, 0.0};
    for (int kk = 0; kk < m_poinQ; ++kk) m_QCDcorR[kk] = 0.0;
    m_QCDcor = 0.0;

    if (ww >= WminZ && ww <= WmaxZ) {
        // LEP1 near Z
        const double x = (ww - WminZ) / (WmaxZ - WminZ);
        int ib = std::min(static_cast<int>(m_poin2 * x), m_poin2 - 1);
        double h = x * m_poin2 - ib;

        if (m_poTh2 == 0) {
            for (int kk = 0; kk < m_poinG; ++kk) {
                auto a = czz[ib  ][0][kk][kf];
                auto b = czz[ib+1][0][kk][kf];
                m_GSW[kk] = lerp(a, b, h);
            }
        } else {
            double y  = 0.5 * (1.0 + CosThe);
            int jb    = std::min(static_cast<int>(m_poTh2 * y), m_poTh2 - 1);
            double hy = y * m_poTh2 - jb;
            for (int kk=0; kk<m_poinG; ++kk) {
                auto c00 = czz[ib  ][jb  ][kk][kf];
                auto c10 = czz[ib+1][jb  ][kk][kf];
                auto c01 = czz[ib  ][jb+1][kk][kf];
                auto c11 = czz[ib+1][jb+1][kk][kf];
                m_GSW[kk] = lerp( lerp(c00, c10, h), lerp(c01, c11, h), hy );
            }
        }
        for (int kk=0; kk<m_poinQ; ++kk) {
            double a = szz[ib  ][kk][kf];
            double b = szz[ib+1][kk][kf];
            m_QCDcorR[kk] = lerp(a,b,h);
        }

    } else if (ww >= m_WminLEP1 && ww <= m_WmaxLEP1) {
        // LEP1 outside Z and low energies (log spacing)
        const double x = std::log(ww / m_WminLEP1) / std::log(m_WmaxLEP1 / m_WminLEP1);
        int ib = std::min(static_cast<int>(m_poin1 * x), m_poin1 - 1);
        double h = x * m_poin1 - ib;

        for (int kk=0; kk<m_poinG; ++kk) {
            auto a = cyy[ib ][kk][kf];
            auto b = cyy[ib+1][kk][kf];
            m_GSW[kk] = lerp(a,b,h);
        }
        for (int kk=0; kk<m_poinQ; ++kk) {
            double a = syy[ib  ][kk][kf];
            double b = syy[ib+1][kk][kf];
            m_QCDcorR[kk] = lerp(a,b,h);
        }

    } else if (ww >= m_WmaxLEP1 && ww <= m_WmaxLEP2) {
        // LEP2 region
        const double x = (ww - m_WmaxLEP1) / (m_WmaxLEP2 - m_WmaxLEP1);
        int ib = std::min(static_cast<int>(m_poin3 * x), m_poin3 - 1);
        double h = x * m_poin3 - ib;

        double y  = 0.5 * (1.0 + CosThe);
        int jb    = std::min(static_cast<int>(m_poTh3 * y), m_poTh3 - 1);
        double hy = y * m_poTh3 - jb;

        for (int kk=0; kk<m_poinG; ++kk) {
            auto c00 = ctt[ib  ][jb  ][kk][kf];
            auto c10 = ctt[ib+1][jb  ][kk][kf];
            auto c01 = ctt[ib  ][jb+1][kk][kf];
            auto c11 = ctt[ib+1][jb+1][kk][kf];
            m_GSW[kk] = lerp( lerp(c00, c10, h), lerp(c01, c11, h), hy );
        }
        for (int kk=0; kk<m_poinQ; ++kk) {
            double a = stt[ib  ][kk][kf];
            double b = stt[ib+1][kk][kf];
            m_QCDcorR[kk] = lerp(a,b,h);
        }

    } else if (ww >= m_WmaxLEP2 && ww <= m_WmaxNLC) {
        // NLC region
        const double x = (ww - m_WmaxLEP2) / (m_WmaxNLC - m_WmaxLEP2);
        int ib = std::min(static_cast<int>(m_poin4 * x), m_poin4 - 1);
        double h = x * m_poin4 - ib;

        double y  = 0.5 * (1.0 + CosThe);
        int jb    = std::min(static_cast<int>(m_poTh4 * y), m_poTh4 - 1);
        double hy = y * m_poTh4 - jb;

        for (int kk=0; kk<m_poinG; ++kk) {
            auto c00 = clc[ib  ][jb  ][kk][kf];
            auto c10 = clc[ib+1][jb  ][kk][kf];
            auto c01 = clc[ib  ][jb+1][kk][kf];
            auto c11 = clc[ib+1][jb+1][kk][kf];
            m_GSW[kk] = lerp( lerp(c00, c10, h), lerp(c01, c11, h), hy );
        }
        for (int kk=0; kk<m_poinQ; ++kk) {
            double a = slc[ib  ][kk][kf];
            double b = slc[ib+1][kk][kf];
            m_QCDcorR[kk] = lerp(a,b,h);
        }

    } else {
        // Out of predefined range
        std::cerr << "BornV_InterpoGSW WARNING: s out of predefined range, ww=" << ww << "\n";
    }

    // legacy field
    m_QCDcor = m_QCDcorR[0] - 1.0;
}
