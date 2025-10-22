#include "DipoleQQRijRadCor.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "KKdizet.h"

DipoleQQRijRadCor::DipoleQQRijRadCor(int iqed_enabled)
    : iqed_(iqed_enabled)
{}

// Main calculation function
std::array<std::array<double,4>,4> DipoleQQRijRadCor::calculate(
    double energy,
    double theta,
    double ReA, double ImA,
    double ReB, double ImB,
    double ReX, double ImX,
    double ReY, double ImY,
    double channel
) {
    double e = sqrt(4.0 * pi_ * alpha_);
    double f = Mz_ * sqrt(GF_ * sqrt(2.0));
    double V = sqrt(1.0 - pow(mtau_/energy, 2));
    double gam = energy / mtau_;
    double m = mtau_;

    std::complex<double>  A, B, X, Y, Pg, Pz, Den, vi, vf;
    std::complex<double>  K1, K2, K3, K11, K12, K13, Gammavp, Rho11, Rho12, Rho13;
    std::complex<double>  GSW[100] {};
    double re[100]{}, im[100]{};
    double Svar, CosThetaD;
    std::complex<double>  RhoEW, VPgamma, CorEle, CorFin, CorEleFin, VVCef;
    double Qi, Qf, ai, af, sw2, sw, cw2, cw, t, Fvivf;
    int KFi, KFf, IFGSW, KeyElw, IfPrint;


    if(iqed_>=5) {
        v0_ = 0.0;
        a0_ = 0.0;
    }

    double s = 4.0 * energy * energy;
    double cth = cos(theta);
    double sth = sin(theta);
    double s2th = sin(2.0 * theta);
    double c2th = cos(2.0 * theta);


    //Setting the electroweak corrections using KKdizet class.
    // The dizet table from dizet6.45 is copied from the workMini to the
    // dizet folder before running this code.
    // The dizet table contain four different energy ranges:
    // LEP1 low energies, LEP1 around Z, LEP2, NLC/LHC energies. (a, b, c, d)
    // The seven different form-factors are interpolated from the table.
    // This is done based on centre of mass energy and angle (cos(theta)).

    KFi = 1 ; // initial state electron
    KFf = 15 ; // tau final state
    Svar = s;
    CosThetaD = cth;
    IFGSW = 1;
    IfPrint = 0; // Switch for printing debug information
    if(IFGSW==1) {
        KKdizet::instance().InterpoGSW(KFi, KFf, Svar, CosThetaD);
        KKdizet::instance().GetGSWxy(re, im);
        std::transform(re, re + 7, im,  GSW, [](double r, double i) { return std::complex<double>(r, i); });
        RhoEW         = GSW[0];
        //VPgamma       = GSW[5];
        VPgamma       = 1.0/(2.0- GSW[5]);
        CorEle        = GSW[1];
        CorFin        = GSW[2];
        CorEleFin     = GSW[3];
        sw2           = KKdizet::instance().D_swsq;//0.23113; //sin^2(theta_W)
        Gz_           = Gz_*Svar/(Mz_*Mz_); // Running width
        // std::cout << "GSW[0] = " << GSW[0] << ", GSW[1] = " << GSW[1] << std::endl;
        // std::cout << "GSW[2] = " << GSW[2] << ", GSW[3] = " << GSW[3] << std::endl;
        // std::cout << "GSW[4] = " << GSW[4] << ", GSW[5] = " << GSW[5] << std::endl;
        // std::cout << "GSW[6] = " << GSW[6] << std::endl;
    } else {
        RhoEW     = 1.0;
        VPgamma   = 1.0;
        CorEle    = 1.0;
        CorFin    = 1.0;
        CorEleFin = 1.0;
        sw2       = 0.22351946; //sin(theta_w)^2 from  Eur.Phys.J. C (2019) 79, 480
    }

    K3      = CorEle;
    K1      = CorEle;
    K13     = CorEleFin;
    K2      = CorEle;
    K12     = CorEleFin;
    K11     = CorEleFin;
    Gammavp = VPgamma;
    Rho11   = RhoEW;
    Rho12   = RhoEW;
    Rho13   = RhoEW;



    // QED loop corrections (optional)
    double ArQED = 0.0, AiQED = 0.0;
    double Ar1 = 0.0, Ai1 = 0.0;
    if (iqed_==10) { //  ! Warning: temporarily this contribution is blocked
        ArQED = -alpha_ * m * m / (pi_ * V * s) * log((1 + V) / (1 - V));
        AiQED = alpha_ * m * m / (V * s);
    }

    double ReA1 = ArQED + ReA;
    double ImA1 = AiQED + ImA;

    sw  = sqrt(sw2);
    cw2 = 1.0 - sw2;
    cw  = sqrt(cw2);
    // Mandelstam variables s and t (mass of tau is included in t). s is previously defined.
    t   = pow(m,2) - 0.5 * s * (1.0 - V * cth);

    // std::cout<< "ReA1 = " << ReA1 << ", ImA1 = " << ImA1 << std::endl;
    // std::cout<< "ArQED = " << ArQED << ", AiQED = " << AiQED << std::endl;

    std::array<std::array<double,4>,4> RSM{};
    std::array<std::array<double,4>,4> RDM{};
    std::array<std::array<double,4>,4> R{};

    //Reduce the compiler overhead
    double a0_4 = pow(a0_, 4);
    double a0_2 = pow(a0_, 2);
    double v0_4 = pow(v0_, 4);
    double v0_2 = pow(v0_, 2);
    double Mz_4 = pow(Mz_, 4);
    double Mz_2 = pow(Mz_, 2);
    double Gz_4 = pow(Gz_, 4);
    double Gz_2 = pow(Gz_, 2);
    double gam6 = pow(gam, 6);
    double gam5 = pow(gam, 5);
    double gam4 = pow(gam, 4);
    double gam3 = pow(gam, 3);
    double gam2 = pow(gam, 2);
    double m4   = pow(m, 4);
    double m2   = pow(m, 2);
    double e4   = pow(e, 4);
    double e2   = pow(e, 2);
    double f4   = pow(f, 4);
    double f2   = pow(f, 2);
    double V2   = pow(V, 2);
    double cth2 = pow(cth, 2);
    double sth2 = pow(sth, 2);
    double sw4  = pow(sw, 4);

    double denom = (16.0*gam4*m4+(Gz_2-8.0*gam2*m2)*Mz_2+Mz_4);

    // Denominators of Z propagator with running width
    Den = s - Mz_2 + std::complex<double>(0.0, 1.0) * Gz_ * Mz_;
    // Constant for effective photon like couplings correction coming from Z boson
    Fvivf = 4.0 * (f2/e2) * sw4; //! WARNING: the *sw**4 factor come from vector couplings numerators Fvivf
    // Note: Fvivf can also be written as Fvivf= sw**2/cw**2, or through the constant G_mu:
    // Fvivf= sqrt(2) *G_mu *Mz**2 *sw**4 /(pi*alpha)  This agrees with definition below.

    // FINAL FERMION IS TAU LEPTON (f=lepton) with couplings:
    Qf = -1.0;
    vf = -0.03783;  //! from PDG, used in previous code (v0 coupling)
    af = -0.50123;  //! from PDG, used in previous code (a0 coupling)
    if(IFGSW==1) {
        vf = -0.5 -2.0 * Qf * sw2 * K1;  // is function for lepton
        af = -0.5;
    }

    // Couplings for Initial Fermions e, u or d
    if(channel == 1) { // Leptons
        Qi = -1.0;
        vi = -0.03783;  //! from PDG, used in the previous code
        ai = -0.50123;  //! from PDG, used in the previous code
        // Effective Z propagator  (reduces to 1/Den without rad. corr.):
        Pz = Rho11/Den;
        if(IFGSW==1) {
            vi = -0.5 -2.0 * Qi * sw2 * K1;  // is function for lepton
            ai = -0.5;
            Fvivf = GMu_ * Mz_2 / (alpha_ * sqrt(2.0) * 8.0*pi_) *
                    (sw2*(1.0-sw2))*16.;
            Fvivf = Fvivf * (sw2/(1.0-sw2)); //! why this factor? Because  vi ai couplings in KKMC
                 // are divided by deno= 4 sqrt(m_swsq*(1d0-m_swsq)) in addition multiplied by 2 and we are not using it here
                 // Ve = (2*T3e -4*Qe*m_swsq)/deno and Ae = 2*T3e/deno.
        }
        //Effective photon propagator with v_{if}-v_i*v_f correction to Z exchange
        //  (NOTE: Pz reduces to 1/s if  rad. corr. are off):
        Pg = 1.0/s*(Gammavp + (Fvivf * s/Den) * (Rho11*(K11-K1*K1)));
    } else if(channel == 2) {;}//up quark
    else if(channel == 3) {;}//down quark
    else{
        Qi = 0.0;
        vi = {0.0,0.0};
        ai = 0.0;
        Pz = {0.0,0.0};
        Pg = {0.0,0.0};
    }

    // Complex magnetic and electric form-factors
    A = std::complex<double>(ReA1, ImA1);
    B = std::complex<double>(ReB, ImB);
    X = std::complex<double>(ReX, ImX);
    Y = std::complex<double>(ReY, ImY);

    // The complete translation for each RSM(i,j) and RDM(i,j)

    // CONTRIBUTIONS IN STANDARD MODEL FOR ALL DIPOLE MOMENTS ZERO

    // Mostly appearing terms are pre-calculated to reduce the compiler overhead
    double Qi2 = pow(Qi, 2);
    double Qf2 = pow(Qf, 2);
    double ai2 = pow(ai, 2);
    double af2 = pow(af, 2);
    std::complex<double> vi2 = vi * vi;
    std::complex<double> vf2 = vf * vf;

    std::complex<double> PgC = std::conj(Pg);
    std::complex<double> PzC = std::conj(Pz);
    std::complex<double> viC = std::conj(vi);
    std::complex<double> aiC = std::conj(ai);
    std::complex<double> vfC = std::conj(vf);

    // //Print things out for debugging
    // std::cout << "----------------------------------------" << std::endl;
    // std::cout << "Denominator: " << Den << std::endl;
    // std::cout << "Pg: " << Pg << ", Pz: " << Pz << std::endl;
    // std::cout << "V:" << V << std::endl;
    // std::cout << "Vi: " << vi << ", Ai: " << ai << std::endl;
    // std::cout << "Vf: " << vf << ", Af: " << af << std::endl;
    // std::cout << "Qi: " << Qi << ", Qf: " << Qf << std::endl;
    // std::cout << "A: " << A << ", B: " << B << ", X: " << X << ", Y: " << Y << std::endl;
    // std::cout << "Fvivf: " << Fvivf << std::endl;
    // std::cout << "Svar: " << Svar << std::endl;
    // std::cout << "t: " << t << std::endl;
    // std::cout << "s: " << s << std::endl;
    // std::cout << "theta: " << theta << std::endl;
    // std::cout << "cth: " << cth << ", sth: " << sth << std::endl;
    // std::cout << "gam: " << gam << std::endl;
    // std::cout << "m: " << m << std::endl;
    // std::cout << "e: " << e << ", f: " << f << std::endl;
    // std::cout << "sw2: " << sw2 << ", cw2: " << cw2 << std::endl;
    // std::cout << "alpha: " << alpha_ << ", GMu: " << GMu_ << std::endl;
    // std::cout << "Mz: " << Mz_ << ", Gz: " << Gz_ << std::endl;
    // std::cout << "iqed: " << iqed_ << std::endl;
    // std::cout << "----------------------------------------" << std::endl;

    RSM[0][0] = 4.0*gam2*m4 * std::real(e2*(1.0 + gam2)*Qf*Qi*
      ((e2*Pg*Qf*Qi + f2*Pz*vf*vi) * std::conj(Pg))
    + f2*std::conj(Pz) *
      ( -(af2*f2*(-1.0 + gam2)*Pz*(ai2 + vi*std::conj(vi))) +
        (1.0 + gam2)*std::conj(vf) *
          ( ai2*f2*Pz*vf +
            (e2*Pg*Qf*Qi + f2*Pz*vf*vi) * std::conj(vi) )
      )
    ) * sth2;

    RSM[0][1] = 4.0 * af * f2 * gam4 * m4 * V *
    std::imag(e2 * Pz * Qf * Qi * vi * std::conj(Pg)
        + std::conj(Pz) * (-(ai2 * f2 * Pz * vf)
            - (e2*Pg*Qf*Qi + f2*Pz*vf*vi) * std::conj(vi)
            + f2 * Pz * std::conj(vf) * (ai2 + vi*std::conj(vi))
        )
    ) * sth2;

    RSM[1][0] = RSM[0][1];

    RSM[1][1] = 4.0 * gam2 * (-1.0 + gam2) * m4 *
    std::real(-( e2 * Qf * Qi * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(Pg) )
        + std::conj(Pz) * ( af2 * f4 * Pz * ( ai2 + vi * std::conj(vi) )
        - f2 * std::conj(vf) * (
                ai2 * f2 * Pz * vf
            + ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(vi)
            )
        )
    ) * sth2;


    RSM[0][2] = 4.0 * (gam2*gam) * m4 *
    std::real( e2 * Qf * Qi * std::conj(Pg) *
        ( af * ai * f2 * Pz * V
            - 2.0 * gam2 * (-1.0 + V*V)
            * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * cth )
        +  f2 * std::conj(Pz) *
           ( af * ( ai * V * ( e2 * Pg * Qf * Qi
                            + f2 * Pz * vf * ( vi + std::conj(vi) ) )
                + 2.0 * af * f2 * Pz * ( 1.0 + gam2 * (-1.0 + V2) )
                * ( ai2 + vi * std::conj(vi) ) * cth )
        +
        std::conj(vf) * (
            ai * f2 * Pz * ( af * V * vi
                            - 2.0 * ai * gam2 * (-1.0 + V2) * vf * cth )
            +
            std::conj(vi) * ( af * ai * f2 * Pz * V
                            - 2.0 * gam2 * (-1.0 + V2)
                                * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi )
                                * cth )
        )
        )
    ) * sth;

    RSM[2][0] = RSM[0][2];

    RSM[1][2] = 2.0 * af * f2 * (gam3) * m4 * V *
    std::imag( e2 * Pz * Qf * Qi * vi * std::conj(Pg)
        + std::conj(Pz) * ( -(ai2 * f2 * Pz * vf)
            - (e2*Pg*Qf*Qi + f2*Pz*vf*vi) * std::conj(vi)
            + f2 * Pz * std::conj(vf) * (ai2 + vi*std::conj(vi))
        )
    ) * s2th;

    RSM[2][1] = RSM[1][2];

    RSM[2][2] = 2.0 * gam2 * m4 *
    std::real( -( e4 * Pg * (Qf2) * (Qi2) * std::conj(Pg) *
        ( 1.0 - 3.0*gam2
            + ( -1.0 + 3.0*gam2 + 4.0*gam4 * (-1.0 + V2) ) * c2th ) )
        - e2 * f2 * Qf * Qi * ( Pz * std::conj(Pg) *
            ( -4.0*af*ai*gam2*V*cth
                + vf*vi * ( 1.0 - 3.0*gam2
                            + ( -1.0 + 3.0*gam2 + 4.0*gam4 * (-1.0 + V2) )
                            * c2th ) )
        + Pg * std::conj(Pz) *
            ( -4.0*af*ai*gam2*V*cth
                + std::conj(vf)*std::conj(vi) *
                ( 1.0 - 3.0*gam2
                    + ( -1.0 + 3.0*gam2 + 4.0*gam4 * (-1.0 + V2) )
                    * c2th ) ) )
        + f4 * Pz * std::conj(Pz) * (
            af * ( 4.0*ai*gam2*V*vf * ( vi + std::conj(vi) ) * cth
                + af * ( ai2 + vi*std::conj(vi) ) *
                    ( 1.0 + gam2 * ( -1.0 + 4.0*V2 )
                    + ( -1.0 + 5.0*gam2 + 4.0*gam4 * ( -1.0 + V2 ) )
                        * c2th ) )
        + std::conj(vf) * (
                std::conj(vi) * ( 4.0*af*ai*gam2*V*cth
                                + vf*vi * ( -1.0 + 3.0*gam2
                                    + ( 1.0 - 3.0*gam2 - 4.0*gam4 * ( -1.0 + V2 ) )
                                    * c2th ) )
            + ai * ( 4.0*af*gam2*V*vi*cth
                    - ai*vf * ( 1.0 - 3.0*gam2
                        + ( -1.0 + 3.0*gam2 + 4.0*gam4 * ( -1.0 + V2 ) )
                        * c2th ) )
            )
        )
    );


    RSM[0][3] = -4.0 * f2 * gam3 * m4 *
    std::real( e2 * Pz * Qf * Qi * std::conj(Pg) *
        ( 2.0*ai*vf + af*V*vi*cth )
        + std::conj(Pz) * ( af*V *
            ( ai2 * f2 * Pz * vf
                + ( e2*Pg*Qf*Qi + f2*Pz*vf*vi ) * std::conj(vi) )
            * cth
        + std::conj(vf) *
            ( 2.0*ai * ( e2*Pg*Qf*Qi + f2*Pz*vf*( vi + std::conj(vi) ) )
                + af * f2 * Pz * V * ( ai2 + vi*std::conj(vi) ) * cth )
        )
    ) * sth;

    RSM[3][0] = RSM[0][3];

    RSM[1][3] = -4.0 * af * ai * f2 * gam3 * m4 * V *
    std::imag( e2 * Pz * Qf * Qi * std::conj(Pg)
        + std::conj(Pz) * ( -(e2 * Pg * Qf * Qi)
            - f2 * Pz * vf * (vi + std::conj(vi))
            + f2 * Pz * std::conj(vf) * (vi + std::conj(vi))
        )
    ) * sth;

    RSM[3][1] = RSM[1][3];

    RSM[2][3] = -4.0 * f2 * gam2 * m4 *
    std::real( e2 * gam2 * Pz * Qf * Qi * std::conj(Pg) *
        ( 2.0*ai*vf*cth + af*V*vi*(1.0 + cth2) )
        + (std::conj(Pz) * ( std::conj(vi) *
          ( 4.0*af2*ai*f2 * (-1.0 + gam2) * Pz * cth
            + 4.0*ai*f2*gam2 * Pz * vf * std::conj(vf) * cth
            + af*gam2*V *
                ( e2*Pg*Qf*Qi + f2*Pz*vi*( vf + std::conj(vf) ) ) *
                ( 3.0 + c2th ) )
        + ai * ( gam2 * std::conj(vf) *
                ( 4.0*( e2*Pg*Qf*Qi + f2*Pz*vf*vi ) * cth
                + af*ai*f2*Pz*V * ( 3.0 + c2th ) )
            + af*f2*Pz *
                ( 4.0*af*(-1.0 + gam2) * vi * cth
                + ai*gam2*V * vf * ( 3.0 + c2th ) )
        )
        )
    ) / 2.0 );

    RSM[3][2] = RSM[2][3];

    RSM[3][3] = 2.0 * gam2 * m4 *
    std::real( e4 * Pg * Qf2 * Qi2 * std::conj(Pg) *
        ( 1.0 + 3.0*gam2 + (-1.0 + gam2) * c2th )
        + e2 * f2 * Qf * Qi * ( Pz * std::conj(Pg) *
            ( 4.0*af*ai*gam2*V*cth
                + vf*vi * ( 1.0 + 3.0*gam2 + (-1.0 + gam2) * c2th ) )
        + Pg * std::conj(Pz) *
            ( 4.0*af*ai*gam2*V*cth
                + std::conj(vf)*std::conj(vi) *
                ( 1.0 + 3.0*gam2 + (-1.0 + gam2) * c2th ) ) )
        - f4 * Pz * std::conj(Pz) *
        ( -( af * ( 4.0*ai*gam2*V*vf * ( vi + std::conj(vi) ) * cth
                    + af * (-1.0 + gam2) * ( ai2 + vi*std::conj(vi) )
                    * ( 3.0 + c2th ) ) )
        - std::conj(vf) *
            ( ai * ( 4.0*af*gam2*V*vi*cth
                    + ai*vf * ( 1.0 + 3.0*gam2 + (-1.0 + gam2) * c2th ) )
            + std::conj(vi) *
                ( 4.0*af*ai*gam2*V*cth
                + vf*vi * ( 1.0 + 3.0*gam2 + (-1.0 + gam2) * c2th ) )
            )
        )
    );


    // CONTRIBUTION LINEAR IN DIPOLE MOMENTS A, B, X, Y

        // RDM entries
    RDM[0][0] = 8.0 * gam4 * m4 *
    std::real( e2 * Qf * Qi *
        ( e2 * Pg * Qf * Qi * ( A + std::conj(A) )
            + f2 * Pz * vi * ( X + vf * std::conj(A) ) ) * std::conj(Pg)
        + f2 * std::conj(Pz) *
        ( std::conj(vf) * ( ai2 * f2 * Pz * X
                + ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * X ) * std::conj(vi) )
            + ( ai2 * f2 * Pz * vf
                + ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(vi) )
            * std::conj(X) ) ) * sth2;

    RDM[0][1] = 4.0 * gam4 * m4 * V *
    std::real( e2 * Qf * Qi *
        ( e2 * Pg * Qf * Qi * ( B + std::conj(B) )
            + f2 * Pz * vi * ( Y
                            - std::complex<double>(0.0, 1.0) * af * std::conj(A)
                            + vf * std::conj(B) ) ) * std::conj(Pg)
        +  f2 * std::conj(Pz) *
        ( std::conj(vf) *  ( ai2 * f2 * Pz * Y
                + ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * Y ) * std::conj(vi) )
            + std::complex<double>(0.0, 1.0) * af * ( std::conj(vi) *
                    ( A * e2 * Pg * Qf * Qi
                    + f2 * Pz * vi * ( X - std::conj(X) ) )
                + ai2 * f2 * Pz * ( X - std::conj(X) ) )
            + ( ai2 * f2 * Pz * vf
                + ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(vi) )
                * std::conj(Y) ) ) * sth2;

    RDM[1][0] =  -4.0 * gam4 * m4 * V *
    std::real( e2 * Qf * Qi *
        ( e2 * Pg * Qf * Qi * ( B + std::conj(B) )
            + f2 * Pz * vi * ( Y
                                + std::complex<double>(0.0, 1.0) * af * std::conj(A)
                                + vf * std::conj(B) ) ) * std::conj(Pg)
        + f2 * std::conj(Pz) *
        ( std::conj(vf) *
            ( ai2 * f2 * Pz * Y
                + ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * Y ) * std::conj(vi) )
            - std::complex<double>(0.0, 1.0) * af *  ( std::conj(vi) *
                    ( A * e2 * Pg * Qf * Qi
                    + f2 * Pz * vi * ( X - std::conj(X) ) )
                + ai2 * f2 * Pz * ( X - std::conj(X) ) )
            + ( ai2 * f2 * Pz * vf
                + ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(vi) )
                * std::conj(Y) ) ) * sth2;

    RDM[1][1] = 0.0;

    RDM[0][2] =  4.0 * (gam4 * gam) * m4 *
    std::real( -( e2 * Qf * Qi * std::conj(Pg) *
        ( ( -std::complex<double>(0.0, 1.0) ) * ai * f2 * Pz * V *
            ( Y - std::complex<double>(0.0, 1.0) * af * std::conj(A)
                - vf * std::conj(B) )
            + ( e2 * Pg * Qf * Qi * ( -2.0 + V2 ) * ( A + std::conj(A) )
            + f2 * Pz * vi *
                ( ( -2.0 + V2 ) * ( X + vf * std::conj(A) )
                    + std::complex<double>(0.0, 1.0) * af * (V2) * std::conj(B) ))
                    * cth ) )
            + f2 * std::conj(Pz) *
            ( std::conj(vf) * ( ai *
                ( std::complex<double>(0.0, 1.0) * V *
                    ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * Y )
                - ai * f2 * Pz * ( -2.0 + V2 ) * X * cth )
                + std::conj(vi) *
                ( std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V * Y
                    - ( -2.0 + V2 ) * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * X ) *
                        cth ))
            + std::conj(vi) *
            ( ( -std::complex<double>(0.0, 1.0) ) * ai * f2 * Pz * V * vf * std::conj(Y)
                - ( -2.0 + V2 ) * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) *
                    std::conj(X) * cth
                + af * V *  ( ai * f2 * Pz * ( X + std::conj(X) )
                    + std::complex<double>(0.0, 1.0) * V *
                        ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y - std::conj(Y) ) )
                        * cth ) )
            + ai * ( ( -std::complex<double>(0.0, 1.0) ) * V *
                ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(Y)
                + ai * f2 * Pz * ( 2.0 - V2 ) * vf * std::conj(X) * cth
                + af * V * ( A * e2 * Pg * Qf * Qi
                    + f2 * Pz * vi * ( X + std::conj(X) )
                    + std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V
                        * ( Y - std::conj(Y) ) * cth ) ) ) ) * sth;

    RDM[2][0] = 4.0 * (gam4 * gam) * m4 *
    std::real( e2 * Qf * Qi * std::conj(Pg) *
        ( ai * f2 * Pz * V * ( af * std::conj(A)
                                - std::complex<double>(0.0, 1.0) * ( Y - vf * std::conj(B) ) )
            - ( e2 * Pg * Qf * Qi * ( -2.0 + V2 ) * ( A + std::conj(A) )
                + f2 * Pz * vi * ( ( -2.0 + V2 ) * ( X + vf * std::conj(A) )
                - std::complex<double>(0.0, 1.0) * af * (V2) * std::conj(B) ) ) * cth )
        + f2 * std::conj(Pz) *
        ( std::conj(vf) * (
                -( ai * ( std::complex<double>(0.0, 1.0) * V
                            * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * Y )
                        + ai * f2 * Pz * ( -2.0 + V2 ) * X * cth ) )
                + std::conj(vi) *
                    ( -std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V * Y
                    - ( -2.0 + V2 ) * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * X )
                        * cth ) )
            + std::conj(vi) *
            ( std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V * vf * std::conj(Y)
                - ( -2.0 + V2 ) * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi )
                    * std::conj(X) * cth
                + af * V * ( ai * f2 * Pz * ( X + std::conj(X) )
                    - std::complex<double>(0.0, 1.0) * V
                        * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y - std::conj(Y) ) )
                        * cth ) )
            +ai * ( std::complex<double>(0.0, 1.0) * V
                * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(Y)
                + ai * f2 * Pz * ( 2.0 - V2 ) * vf * std::conj(X) * cth
                + af * V * ( A * e2 * Pg * Qf * Qi
                    + f2 * Pz * vi * ( X + std::conj(X) )
                    - std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V
                        * ( Y - std::conj(Y) ) * cth ) ) ) ) * sth;


    RDM[1][2] =  -4.0 * (gam4 * gam) * m4 * V *
    std::imag( e2 * Qf * Qi * std::conj(Pg) *
        ( ai * f2 * Pz * V * ( X - vf * std::conj(A)
        + std::complex<double>(0.0, 1.0) * af * std::conj(B) )
        + std::complex<double>(0.0, 1.0) * ( e2 * Pg * Qf * Qi * ( B + std::conj(B) )
        + f2 * Pz * vi * ( Y + std::complex<double>(0.0, 1.0) * af * std::conj(A)
        + vf * std::conj(B) ) ) * cth )
        + f2 * std::conj(Pz) * ( std::conj(vf) *
            ( ai * ( A * e2 * Pg * Qf * Qi * V + f2 * Pz * V * vi * X
        + std::complex<double>(0.0, 1.0) * ai * f2 * Pz * Y * cth )
        + std::conj(vi) * ( ai * f2 * Pz * V * X
        + std::complex<double>(0.0, 1.0) * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * Y )
                                    * cth ) )
        + std::complex<double>(0.0, 1.0) * ( std::conj(vi) *
            ( std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V * vf * std::conj(X)
        + ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(Y) * cth
        + af * ( ai * f2 * Pz * V * ( Y + std::conj(Y) )
        - std::complex<double>(0.0, 1.0) * ( A * e2 * Pg * Qf * Qi
        + f2 * Pz * vi * ( X - std::conj(X) ) ) * cth ) )
        + ai * ( std::complex<double>(0.0, 1.0) * V * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi )
                        * std::conj(X)
        + ai * f2 * Pz * vf * std::conj(Y) * cth
        + af * ( V * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y + std::conj(Y) ) )
        - std::complex<double>(0.0, 1.0) * ai * f2 * Pz * ( X - std::conj(X) )
                                    * cth ) ) ) ) ) * sth;


    RDM[2][1] =  4.0 * gam5 * m4 * V *
    std::real( e2 * Qf * Qi * std::conj(Pg) *
        (
            ai * f2 * Pz * V *
            ( std::complex<double>(0.0, 1.0) * ( X - vf * std::conj(A) ) + af * std::conj(B) )
            +
            ( e2 * Pg * Qf * Qi * ( B + std::conj(B) )
            + f2 * Pz * vi * ( Y - std::complex<double>(0.0, 1.0) * af * std::conj(A) + vf * std::conj(B) ) )
            * cth
        )
        +
        f2 * std::conj(Pz) *
        (
            std::conj(vf) *
            (
                ai * ( std::complex<double>(0.0, 1.0) * V * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * X )
                    + ai * f2 * Pz * Y * cth )
                +
                std::conj(vi) *
                ( std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V * X
                    + ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * Y ) * cth )
            )
            +
            std::conj(vi) *
            (
                ( -std::complex<double>(0.0, 1.0) ) * ai * f2 * Pz * V * vf * std::conj(X)
                + ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(Y) * cth
                + af * ( ai * f2 * Pz * V * ( Y + std::conj(Y) )
                        + std::complex<double>(0.0, 1.0)
                            * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X - std::conj(X) ) )
                            * cth )
            )
            +
            ai *
            (
                ( -std::complex<double>(0.0, 1.0) ) * V * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(X)
                + ai * f2 * Pz * vf * std::conj(Y) * cth
                + af * ( V * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y + std::conj(Y) ) )
                        + std::complex<double>(0.0, 1.0) * ai * f2 * Pz * ( X - std::conj(X) ) * cth )
            )
        )
    ) * sth;

    RDM[2][2] =
    -8.0 * gam6 * m4 * (-1.0 + V2) * cth *
    std::real(
        e2 * Qf * Qi * std::conj(Pg) *
        ( ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * X ) * cth
            + std::conj(A) * ( af * ai * f2 * Pz * V
                            + ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * cth ) )
        +
        f2 * std::conj(Pz) *
        ( ai * ( af * V * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X + std::conj(X) ) )
                + ai * f2 * Pz * X * std::conj(vf) * cth
                + ai * f2 * Pz * vf * std::conj(X) * cth )
            + std::conj(vi) *
            ( af * ai * f2 * Pz * V * ( X + std::conj(X) )
                + ( ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * X ) * std::conj(vf)
                    + ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(X) )
                * cth ) )
    );


    RDM[0][3] =
    -4.0 * gam3 * m4 *
    std::real(
        e2 * Qf * Qi * std::conj(Pg) *
        (
            ai * f2 * Pz * ( (1.0 + gam2) * ( X + vf * std::conj(A) )
                            - std::complex<double>(0.0, 1.0) * af * (-1.0 + gam2) * std::conj(B) )
            + std::complex<double>(0.0, 1.0) * gam2 * V *
                ( e2 * Pg * Qf * Qi * ( B - std::conj(B) )
                + f2 * Pz * vi * ( Y - std::complex<double>(0.0, 1.0) * af * std::conj(A) - vf * std::conj(B) ) )
                * cth
        )
        +
        f2 * std::conj(Pz) *
        (
            std::conj(vf) *
            (
                ai * ( (1.0 + gam2) * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * X )
                    + std::complex<double>(0.0, 1.0) * ai * f2 * gam2 * Pz * V * Y * cth )
                + std::conj(vi) *
                    ( ai * f2 * (1.0 + gam2) * Pz * X
                    + std::complex<double>(0.0, 1.0) * gam2 * V
                        * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * Y )
                        * cth )
            )
            + ai *
            (
                (1.0 + gam2) * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(X)
                - std::complex<double>(0.0, 1.0) * ai * f2 * gam2 * Pz * V * vf * std::conj(Y) * cth
                + af * ( std::complex<double>(0.0, 1.0) * (-1.0 + gam2)
                        * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y - std::conj(Y) ) )
                        + ai * f2 * gam2 * Pz * V * ( X + std::conj(X) ) * cth )
            )
            + std::conj(vi) *
            (
                ai * f2 * (1.0 + gam2) * Pz * vf * std::conj(X)
                - std::complex<double>(0.0, 1.0) * gam2 * V
                    * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(Y) * cth
                + af * ( std::complex<double>(0.0, 1.0) * ai * f2 * (-1.0 + gam2) * Pz * ( Y - std::conj(Y) )
                        + gam2 * V * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X + std::conj(X) ) )
                        * cth )
            )
        )
    ) * sth;


    RDM[3][0] =
    -4.0 * gam3 * m4 *
    std::real(
        e2 * Qf * Qi * std::conj(Pg) *
        (
            ai * f2 * Pz * ( (1.0 + gam2) * ( X + vf * std::conj(A) )
                            + std::complex<double>(0.0, 1.0) * af * (-1.0 + gam2) * std::conj(B) )
            - std::complex<double>(0.0, 1.0) * gam2 * V *
                ( e2 * Pg * Qf * Qi * ( B - std::conj(B) )
                + f2 * Pz * vi * ( Y + std::complex<double>(0.0, 1.0) * af * std::conj(A) - vf * std::conj(B) ) )
                * cth
        )
        +
        f2 * std::conj(Pz) *
        (
            std::conj(vf) *
            (
                ai * ( (1.0 + gam2) * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * X )
                    - std::complex<double>(0.0, 1.0) * ai * f2 * gam2 * Pz * V * Y * cth )
                + std::conj(vi) *
                    ( ai * f2 * (1.0 + gam2) * Pz * X
                    - std::complex<double>(0.0, 1.0) * gam2 * V
                        * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * Y )
                        * cth )
            )
            + ai *
            (
                (1.0 + gam2) * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(X)
                + std::complex<double>(0.0, 1.0) * ai * f2 * gam2 * Pz * V * vf * std::conj(Y) * cth
                + af * ( -std::complex<double>(0.0, 1.0) * (-1.0 + gam2)
                        * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y - std::conj(Y) ) )
                        + ai * f2 * gam2 * Pz * V * ( X + std::conj(X) ) * cth )
            )
            + std::conj(vi) *
            (
                ai * f2 * (1.0 + gam2) * Pz * vf * std::conj(X)
                + std::complex<double>(0.0, 1.0) * gam2 * V
                    * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(Y) * cth
                + af * ( -std::complex<double>(0.0, 1.0) * ai * f2 * (-1.0 + gam2) * Pz * ( Y - std::conj(Y) )
                        + gam2 * V * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X + std::conj(X) ) )
                        * cth )
            )
        )
    ) * sth;


    RDM[1][3] =
    4.0 * gam3 * m4 *
    std::real(
        e2 * Qf * Qi * std::conj(Pg) *
        ( ai * f2 * gam2 * Pz * V *
            ( Y + std::complex<double>(0.0, 1.0) * af * std::conj(A) + vf * std::conj(B) )
            - std::complex<double>(0.0, 1.0) * (-1.0 + gam2) *
            ( e2 * Pg * Qf * Qi * ( A - std::conj(A) )
                + f2 * Pz * vi * ( X - vf * std::conj(A)
                                    + std::complex<double>(0.0, 1.0) * af * std::conj(B) ) )
            * cth
        )
        +
        f2 * std::conj(Pz) *
        ( std::conj(vf) *
            ( ai * ( gam2 * V * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * Y )
                    - std::complex<double>(0.0, 1.0) * ai * f2 * (-1.0 + gam2) * Pz * X * cth )
                + std::conj(vi) *
                    ( ai * f2 * gam2 * Pz * V * Y
                    - std::complex<double>(0.0, 1.0) * (-1.0 + gam2)
                        * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * X ) * cth )
            )
            + ai *
            ( gam2 * V * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(Y)
                + std::complex<double>(0.0, 1.0) * ai * f2 * (-1.0 + gam2) * Pz * vf * std::conj(X) * cth
                + af * ( -std::complex<double>(0.0, 1.0) * gam2 * V
                        * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X - std::conj(X) ) )
                        + ai * f2 * (-1.0 + gam2) * Pz * ( Y + std::conj(Y) ) * cth )
            )
            + std::conj(vi) *
            ( std::complex<double>(0.0, 1.0) *
                ( -std::complex<double>(0.0, 1.0) * ai * f2 * gam2 * Pz * V * vf * std::conj(Y)
                    + (-1.0 + gam2) * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi )
                        * std::conj(X) * cth )
                + af * ( -std::complex<double>(0.0, 1.0) * ai * f2 * gam2 * Pz * V * ( X - std::conj(X) )
                        + (-1.0 + gam2)
                            * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y + std::conj(Y) ) )
                            * cth )
            )
        )
    ) * sth;


    RDM[3][1] =
    4.0 * gam3 * m4 *
    std::real(
        e2 * Qf * Qi * std::conj(Pg) *
        (
            -( ai * f2 * gam2 * Pz * V * ( Y - std::complex<double>(0.0,1.0) * af * std::conj(A) + vf * std::conj(B) ) )
            - std::complex<double>(0.0,1.0) * (-1.0 + gam2) *
                ( e2 * Pg * Qf * Qi * ( A - std::conj(A) )
                + f2 * Pz * vi * ( X - vf * std::conj(A) - std::complex<double>(0.0,1.0) * af * std::conj(B) ) )
                * std::cos(theta)
        )
        - f2 * std::conj(Pz) *
        (
            std::conj(vf) *
            (
                ai * ( gam2 * V * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * Y )
                    + std::complex<double>(0.0,1.0) * ai * f2 * (-1.0 + gam2) * Pz * X * cth )
                + std::conj(vi) *
                    ( ai * f2 * gam2 * Pz * V * Y
                    + std::complex<double>(0.0,1.0) * (-1.0 + gam2)
                        * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * X )
                        * cth )
            )
            + ai *
            ( gam2 * V * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(Y)
                - std::complex<double>(0.0,1.0) * ai * f2 * (-1.0 + gam2) * Pz * vf * std::conj(X) * cth
                + af * ( std::complex<double>(0.0,1.0) * gam2 * V
                        * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X - std::conj(X) ) )
                        + ai * f2 * (-1.0 + gam2) * Pz * ( Y + std::conj(Y) ) * cth )
            )
            + std::conj(vi) *
            (
                ( -std::complex<double>(0.0,1.0) ) *
                ( std::complex<double>(0.0,1.0) * ai * f2 * gam2 * Pz * V * vf * std::conj(Y)
                    + (-1.0 + gam2) * ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi )
                        * std::conj(X) * cth )
                + af * ( std::complex<double>(0.0,1.0) * ai * f2 * gam2 * Pz * V * ( X - std::conj(X) )
                        + (-1.0 + gam2)
                            * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y + std::conj(Y) ) )
                            * cth )
            )
        )
    ) * sth;



    RDM[2][3] =
    -2.0 * gam4 * m4 *
    std::real(
        f2 * std::conj(Pz) *
        (
            -( af * V *
                ( ai2 * f2 * Pz * ( X + std::conj(X) )
                + std::conj(vi) * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X + std::conj(X) ) ) )
            * ( -2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) )
            + std::complex<double>(0.0, 1.0) *
            (
                std::complex<double>(0.0, -4.0) * ai *
                ( e2 * Pg * Qf * Qi + f2 * Pz * vf * ( vi + std::conj(vi) ) )
                * std::conj(X) * std::cos(theta)
                + V * ( ai2 * f2 * Pz * vf
                        + ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(vi) )
                    * std::conj(Y)
                    * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th )
            )
            + std::conj(vf) *
            (
                std::complex<double>(0.0, 1.0) * ai *
                ( std::complex<double>(0.0, -4.0) * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * X ) * cth
                    + ai * f2 * Pz * V * Y
                        * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) )
                + std::conj(vi) *
                ( 4.0 * ai * f2 * Pz * X * cth
                    - std::complex<double>(0.0, 1.0) * V
                        * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * Y )
                        * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) )
            )
        )
        + e2 * Qf * Qi * std::conj(Pg) *
        (
            4.0 * ai * f2 * Pz * ( X + vf * std::conj(A) ) * cth
            - std::complex<double>(0.0, 1.0) * V *
            ( e2 * Pg * Qf * Qi * ( B - std::conj(B) )
                * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th )
                + f2 * Pz * vi *
                    ( std::complex<double>(0.0, -1.0) * af * std::conj(A)
                        * ( -2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th )
                    + ( Y - vf * std::conj(B) )
                        * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) ) )
        )
    );


    RDM[3][2] =
    -2.0 * gam4 * m4 *
    std::real(
        f2 * std::conj(Pz) *
        (
            -( af * V *
                ( ai2 * f2 * Pz * ( X + std::conj(X) )
                + std::conj(vi) * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X + std::conj(X) ) ) )
            * ( -2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) )
            - std::complex<double>(0.0, 1.0) *
            (
                std::complex<double>(0.0, 4.0) * ai *
                ( e2 * Pg * Qf * Qi + f2 * Pz * vf * ( vi + std::conj(vi) ) )
                * std::conj(X) * cth
                + V * ( ai2 * f2 * Pz * vf
                        + ( e2 * Pg * Qf * Qi + f2 * Pz * vf * vi ) * std::conj(vi) )
                    * std::conj(Y)
                    * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th )
            )
            + std::conj(vf) *
            (
                std::complex<double>(0.0, 1.0) * ai *
                ( std::complex<double>(0.0, -4.0) * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * X )
                    * cth
                    + ai * f2 * Pz * V * Y
                        * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) )
                + std::conj(vi) *
                ( 4.0 * ai * f2 * Pz * X * cth
                    + std::complex<double>(0.0, 1.0) * V
                        * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * Y )
                        * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) )
            )
        )
        + e2 * Qf * Qi * std::conj(Pg) *
        (
            4.0 * ai * f2 * Pz * ( X + vf * std::conj(A) ) * cth
            + std::complex<double>(0.0, 1.0) * V *
            ( e2 * Pg * Qf * Qi * ( B - std::conj(B) )
                * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th )
                + f2 * Pz * vi *
                    ( std::complex<double>(0.0, 1.0) * af * std::conj(A)
                        * ( -2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th )
                    + ( Y - vf * std::conj(B) )
                        * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) ) )
        )
    );


    RDM[3][3] =
    8.0 * gam4 * m4 *
    std::real(
        e4 * Pg * Qf2 * Qi2 * ( A + std::conj(A) ) * std::conj(Pg)
        +
        e2 * f2 * Qf * Qi * (
            Pz * vi * X * std::conj(Pg)
        + A * Pg * std::conj(Pz) * std::conj(vf) * std::conj(vi)
        + Pg * std::conj(Pz) * std::conj(vi) * std::conj(X)
        + A * af * ai * Pg * V * std::conj(Pz) * cth
        + Pz * std::conj(A) * std::conj(Pg) * ( vf * vi + af * ai * V * cth )
        )
        +
        f4 * Pz * std::conj(Pz) * (
            ai * ( ai * X * std::conj(vf)
                + ai * vf * std::conj(X)
                + af * V * vi * X * cth
                + af * V * vi * std::conj(X) * cth )
        + std::conj(vi) * (
                vi * ( X * std::conj(vf) + vf * std::conj(X) )
            + af * ai * V * ( X + std::conj(X) ) * cth
            )
        )
    );

    // std::cout<< "RSM03: " << RSM[0][3] << std::endl;
    // std::cout<< "RDM03: " << RDM[0][3] << std::endl;
    // std::cout<< "X+vf: " << X+vf << std::endl;
    // Combine RSM and RDM to form R
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            R[i][j] = RSM[i][j] + RDM[i][j];
        }
    }

    return R;
}
