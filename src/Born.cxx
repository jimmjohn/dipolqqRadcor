#include "BornV.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

static bool read_fortran_double(std::istream& in, double& out) {
    std::string tok;
    if (!(in >> tok)) return false;
    for (char& c : tok) if (c == 'D' || c == 'd') c = 'E'; // D exponent -> E
    out = std::stod(tok);
    return true;
}

BornV::BornV() {
    if (m_poin1<0 || m_poin2<0 || m_poin3<0 || m_poin4<0
        || m_poTh2<0 || m_poTh3<0 || m_poTh4<0
        || m_poinG<=0 || m_poinQ<=0) {
        throw std::invalid_argument("BornV: invalid grid sizes");
    }
    allocate_all_();
}

void BornV::allocate_all_() {
    // outputs
    m_GSW.assign(static_cast<std::size_t>(m_poinG), {0.0,0.0});
    m_QCDcorR.assign(static_cast<std::size_t>(m_poinQ), 0.0);
    m_QCDcor = 0.0;

    // tables
    m_cyy.assign(static_cast<std::size_t>((m_poin1+1) * m_poinG * 16), {0.0,0.0});
    m_syy.assign(static_cast<std::size_t>((m_poin1+1) * m_poinQ * 16), 0.0);

    m_czz.assign(static_cast<std::size_t>((m_poin2+1) * (m_poTh2+1) * m_poinG * 16), {0.0,0.0});
    m_szz.assign(static_cast<std::size_t>((m_poin2+1) * m_poinQ * 16), 0.0);

    m_ctt.assign(static_cast<std::size_t>((m_poin3+1) * (m_poTh3+1) * m_poinG * 16), {0.0,0.0});
    m_stt.assign(static_cast<std::size_t>((m_poin3+1) * m_poinQ * 16), 0.0);

    m_clc.assign(static_cast<std::size_t>((m_poin4+1) * (m_poTh4+1) * m_poinG * 16), {0.0,0.0});
    m_slc.assign(static_cast<std::size_t>((m_poin4+1) * m_poinQ * 16), 0.0);
}

// ---- all-files loader (BornV_ReadAll) ----
void BornV::BornV_ReadAll(const std::string& dizetDir) {
    // PDG codes used in your Fortran
    const int KFdown=1, KFup=2, KFstran=3, KFcharm=4, KFbotom=5, KFtop=6;
    const int KFel=11, KFelnu=12, KFmu=13, KFmunu=14, KFtau=15, KFtaunu=16;

    auto path = [&](const char* name){ return dizetDir + "/" + std::string(name); };

    BornV_ReadFile(path("table.down"),  KFdown);
    BornV_ReadFile(path("table.up"),    KFup);
    BornV_ReadFile(path("table.down"),  KFstran);
    BornV_ReadFile(path("table.up"),    KFcharm);
    BornV_ReadFile(path("table.botom"), KFbotom);
    BornV_ReadFile(path("table.mu"),    KFmu);
    BornV_ReadFile(path("table.tau"),   KFtau);
    BornV_ReadFile(path("table.nue"),   KFelnu);
    BornV_ReadFile(path("table.numu"),  KFmunu);
    BornV_ReadFile(path("table.nutau"), KFtaunu);

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
    double amz=0, amh=0, amtop=0;
    if (!read_fortran_double(in, amz) ||
        !read_fortran_double(in, amh) ||
        !read_fortran_double(in, amtop) ||
        !read_fortran_double(in, m_swsq) ||
        !read_fortran_double(in, m_gammz) ||
        !read_fortran_double(in, m_MW) ||
        !read_fortran_double(in, m_GammW)) {
        throw std::runtime_error("BornV_ReadFile: header read failed");
    }

    // Optional consistency check (matches Fortran STOPs)
    auto bad = [](double expected, double got){
        if (expected == 0.0) return false; // skip check if not preset
        return std::abs(1.0 - expected / got) > 1e-5;
    };
    if (amz == 0.0 || amh == 0.0 || amtop == 0.0
        || bad(m_MZ, amz) || bad(m_amh, amh) || bad(m_amtop, amtop)) {
        std::ostringstream oss;
        oss << "BornV_ReadFile: CHECK mz,mh,mtop (expected/read) = "
            << m_MZ << "/" << amz << "  "
            << m_amh << "/" << amh << "  "
            << m_amtop << "/" << amtop;
        throw std::runtime_error(oss.str());
    }

    // ---- helpers to read header cards ----
    auto read_header_3 = [&](char& chr, int& n, double& ww){
        // (a, i4, f11.5, i4, f11.5) but we only need first 3 here
        in >> std::ws;
        if (!(in >> chr)) throw std::runtime_error("header3: chr");
        if (!(in >> n))   throw std::runtime_error("header3: n");
        if (!read_fortran_double(in, ww)) throw std::runtime_error("header3: ww");
    };
    auto read_header_5 = [&](char& chr, int& n1, double& ww, int& n2, double& cosi){
        in >> std::ws;
        if (!(in >> chr)) throw std::runtime_error("header5: chr");
        if (!(in >> n1))  throw std::runtime_error("header5: n1");
        if (!read_fortran_double(in, ww)) throw std::runtime_error("header5: ww");
        if (!(in >> n2))  throw std::runtime_error("header5: n2");
        if (!read_fortran_double(in, cosi)) throw std::runtime_error("header5: cosi");
    };

    // The complex numbers are sometimes written with no space after the exponent. So we read them one by one.
    // tolerating missing 1x space
    auto read_g13_fields = [&](int need, std::vector<double>& out) {
    out.clear(); out.reserve(need);
      for (int t = 0; t < need; ++t) {
        // skip line breaks
        while (in && (in.peek() == '\n' || in.peek() == '\r')) in.get();
        char buf[13];
        in.read(buf, 13);
        if (!in) throw std::runtime_error("EOF while reading g13 field");
        std::string s(buf, 13);
        // trim spaces
        auto l = s.find_first_not_of(" \t");
        auto r = s.find_last_not_of(" \t");
        s = (l == std::string::npos) ? std::string() : s.substr(l, r - l + 1);
        for (char& c : s) if (c=='D' || c=='d') c='E';
        out.push_back(s.empty() ? 0.0 : std::stod(s));
        // optional 1x separator
        if (in && in.peek() == ' ') in.get();
        // swallow any CR/LF before next field
        while (in && (in.peek() == '\n' || in.peek() == '\r')) in.get();
      }
    };

    auto read_complex_row = [&](std::vector<std::complex<double>>& dst, int count){
        dst.resize(count);
        std::vector<double> vals;
        read_g13_fields(2*count, vals);
        dst.resize(count);
        for (int i=0;i<count;++i) {
            dst[i] = {vals[2*i], vals[2*i+1]};
        }
    };
    auto read_real_row = [&](std::vector<double>& dst, int count){
        dst.resize(count);
        for (int i=0;i<count;++i) {
            if (!read_fortran_double(in, dst[i])) throw std::runtime_error("real row");
        }
    };

    char chr; int n, n1, n2; double ww, cosi;
    std::vector<std::complex<double>> rowC;
    std::vector<double> rowR;

    // ---- basic range (cyy, syy) ----
    for (int i=0; i<=m_poin1; ++i) {
        read_header_3(chr, n, ww);
        read_complex_row(rowC, m_poinG);
        for (int k=0; k<m_poinG; ++k) m_cyy[idx_cyy(i,k,kf)] = rowC[k];

        read_real_row(rowR, m_poinQ);
        for (int k=0; k<m_poinQ; ++k) m_syy[idx_syy(i,k,kf)] = rowR[k];
    }

    // ---- Z pole (czz, szz) ----
    for (int i=0; i<=m_poin2; ++i) {
        for (int j=0; j<=m_poTh2; ++j) {
            read_header_3(chr, n, ww);
            read_complex_row(rowC, m_poinG);
            for (int k=0; k<m_poinG; ++k) m_czz[idx_czz(i,j,k,kf)] = rowC[k];
        }
        read_real_row(rowR, m_poinQ);
        for (int k=0; k<m_poinQ; ++k) m_szz[idx_szz(i,k,kf)] = rowR[k];
    }

    // ---- LEP2 (ctt, stt) ----
    for (int i=0; i<=m_poin3; ++i) {
        for (int j=0; j<=m_poTh3; ++j) {
            read_header_5(chr, n1, ww, n2, cosi);
            read_complex_row(rowC, m_poinG);
            for (int k=0; k<m_poinG; ++k) m_ctt[idx_ctt(i,j,k,kf)] = rowC[k];
        }
        read_real_row(rowR, m_poinQ);
        for (int k=0; k<m_poinQ; ++k) m_stt[idx_stt(i,k,kf)] = rowR[k];
    }

    // ---- NLC (clc, slc) ----
    for (int i=0; i<=m_poin4; ++i) {
        for (int j=0; j<=m_poTh4; ++j) {
            read_header_5(chr, n1, ww, n2, cosi);
            read_complex_row(rowC, m_poinG);
            for (int k=0; k<m_poinG; ++k) m_clc[idx_clc(i,j,k,kf)] = rowC[k];
        }
        read_real_row(rowR, m_poinQ);
        for (int k=0; k<m_poinQ; ++k) m_slc[idx_slc(i,k,kf)] = rowR[k];
    }
}

void BornV::BornV_InterpoGSW(int KFf, double svar, double CosThe) {
    if (KFf < 1 || KFf > 16) throw std::runtime_error("BornV_InterpoGSW: KFf out of [1,16]");
    const int kf = KFf - 1;

    const double ww = std::sqrt(std::max(0.0, svar));
    const double WminZ = m_MZ - m_WdelZ;
    const double WmaxZ = m_MZ + m_WdelZ;

    auto lerp = [](auto a, auto b, double t){ return a*(1.0 - t) + b*t; };

    // Clear outputs by default
    std::fill(m_GSW.begin(), m_GSW.end(), std::complex<double>{0.0,0.0});
    std::fill(m_QCDcorR.begin(), m_QCDcorR.end(), 0.0);
    m_QCDcor = 0.0;

    if (ww >= WminZ && ww <= WmaxZ) {
        // LEP1 near Z
        const double x = (ww - WminZ) / (WmaxZ - WminZ);
        int ib = std::min(static_cast<int>(m_poin2 * x), m_poin2 - 1);
        double h = x * m_poin2 - ib;

        if (m_poTh2 == 0) {
            for (int kk=0; kk<m_poinG; ++kk) {
                auto a = m_czz[idx_czz(ib,   0, kk, kf)];
                auto b = m_czz[idx_czz(ib+1, 0, kk, kf)];
                m_GSW[kk] = lerp(a,b,h);
            }
        } else {
            double y  = 0.5 * (1.0 + CosThe);
            int jb    = std::min(static_cast<int>(m_poTh2 * y), m_poTh2 - 1);
            double hy = y * m_poTh2 - jb;
            for (int kk=0; kk<m_poinG; ++kk) {
                auto c00 = m_czz[idx_czz(ib,   jb,   kk, kf)];
                auto c10 = m_czz[idx_czz(ib+1, jb,   kk, kf)];
                auto c01 = m_czz[idx_czz(ib,   jb+1, kk, kf)];
                auto c11 = m_czz[idx_czz(ib+1, jb+1, kk, kf)];
                m_GSW[kk] = lerp(lerp(c00, c10, h), lerp(c01, c11, h), hy);
            }
        }
        for (int kk=0; kk<m_poinQ; ++kk) {
            double a = m_szz[idx_szz(ib,   kk, kf)];
            double b = m_szz[idx_szz(ib+1, kk, kf)];
            m_QCDcorR[kk] = lerp(a,b,h);
        }

    } else if (ww >= m_WminLEP1 && ww <= m_WmaxLEP1) {
        // LEP1 outside Z and low energies (log spacing)
        const double x = std::log(ww / m_WminLEP1) / std::log(m_WmaxLEP1 / m_WminLEP1);
        int ib = std::min(static_cast<int>(m_poin1 * x), m_poin1 - 1);
        double h = x * m_poin1 - ib;

        for (int kk=0; kk<m_poinG; ++kk) {
            auto a = m_cyy[idx_cyy(ib,   kk, kf)];
            auto b = m_cyy[idx_cyy(ib+1, kk, kf)];
            m_GSW[kk] = lerp(a,b,h);
        }
        for (int kk=0; kk<m_poinQ; ++kk) {
            double a = m_syy[idx_syy(ib,   kk, kf)];
            double b = m_syy[idx_syy(ib+1, kk, kf)];
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
            auto c00 = m_ctt[idx_ctt(ib,   jb,   kk, kf)];
            auto c10 = m_ctt[idx_ctt(ib+1, jb,   kk, kf)];
            auto c01 = m_ctt[idx_ctt(ib,   jb+1, kk, kf)];
            auto c11 = m_ctt[idx_ctt(ib+1, jb+1, kk, kf)];
            m_GSW[kk] = lerp(lerp(c00, c10, h), lerp(c01, c11, h), hy);
        }
        for (int kk=0; kk<m_poinQ; ++kk) {
            double a = m_stt[idx_stt(ib,   kk, kf)];
            double b = m_stt[idx_stt(ib+1, kk, kf)];
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
            auto c00 = m_clc[idx_clc(ib,   jb,   kk, kf)];
            auto c10 = m_clc[idx_clc(ib+1, jb,   kk, kf)];
            auto c01 = m_clc[idx_clc(ib,   jb+1, kk, kf)];
            auto c11 = m_clc[idx_clc(ib+1, jb+1, kk, kf)];
            m_GSW[kk] = lerp(lerp(c00, c10, h), lerp(c01, c11, h), hy);
        }
        for (int kk=0; kk<m_poinQ; ++kk) {
            double a = m_slc[idx_slc(ib,   kk, kf)];
            double b = m_slc[idx_slc(ib+1, kk, kf)];
            m_QCDcorR[kk] = lerp(a,b,h);
        }

    } else {
        // Out of predefined range
        std::cerr << "BornV_InterpoGSW WARNING: s out of predefined range, ww=" << ww << "\n";
    }

    // legacy field
    if (!m_QCDcorR.empty()) m_QCDcor = m_QCDcorR[0] - 1.0;
}
