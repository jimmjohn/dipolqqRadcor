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
    double a0_2 = a0_*a0_;
    double a0_4 = a0_2*a0_2;
    double v0_2 = v0_*v0_;
    double v0_4 = v0_2*v0_2;
    double Mz_2 = Mz_*Mz_;
    double Mz_4 = Mz_2*Mz_2;
    double Gz_2 = Gz_*Gz_;
    double Gz_4 = Gz_2*Gz_2;
    double gam2 = gam*gam;
    double gam3 = gam*gam2;
    double gam4 = gam2*gam2;
    double gam5 = gam2*gam3;
    double gam6 = gam2*gam4;
    double m2   = m*m;
    double m4   = m2*m2;
    double e2   = e*e;
    double e4   = e2*e2;
    double f2   = f*f;
    double f4   = f2*f2;

    double V2   = V*V;
    double cth2 = cth*cth;
    double sth2 = sth*sth;
    double sw4  = sw*sw*sw*sw;

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

    // Mostly appearing terms are pre-calculated to reduce the overhead
    double Qi2 = Qi * Qi;
    double Qf2 = Qf * Qf;
    double ai2 = ai * ai;
    double af2 = af * af;
    std::complex<double> vi2 = vi * vi;
    std::complex<double> vf2 = vf * vf;

    std::complex<double> PgC = std::conj(Pg);
    std::complex<double> PzC = std::conj(Pz);
    std::complex<double> viC = std::conj(vi);
    std::complex<double> vfC = std::conj(vf);
    std::complex<double> AC  = std::conj(A);
    std::complex<double> BC  = std::conj(B);
    std::complex<double> XC  = std::conj(X);
    std::complex<double> YC  = std::conj(Y);

    // Did profiling and found that the complex multiplications are taking time.
    // Consider optimizing the following lines if they are a bottleneck:
    std::complex<double>  e2PgQfQiPlusf2Pzvfvi = e2*Pg*Qf*Qi + f2*Pz*vf*vi;
    std::complex<double>  Ae2PgQfQiPlusf2PzvfviX =  A*e2*Pg*Qf*Qi + f2*Pz*vi*X;
    std::complex<double>  Be2PgQfQiPlusf2PzvfviY =  B*e2*Pg*Qf*Qi + f2*Pz*vi*Y;



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
      (e2PgQfQiPlusf2Pzvfvi * PgC)
    + f2*PzC *
      ( -(af2*f2*(-1.0 + gam2)*Pz*(ai2 + vi*viC)) +
        (1.0 + gam2)*vfC *
          ( ai2*f2*Pz*vf +
            e2PgQfQiPlusf2Pzvfvi * viC )
      )
    ) * sth2;

    RSM[0][1] = 4.0 * af * f2 * gam4 * m4 * V *
    std::imag(e2 * Pz * Qf * Qi * vi * PgC
        + PzC * (-(ai2 * f2 * Pz * vf)
            - e2PgQfQiPlusf2Pzvfvi * viC
            + f2 * Pz * vfC * (ai2 + vi*viC)
        )
    ) * sth2;

    RSM[1][0] = RSM[0][1];

    RSM[1][1] = 4.0 * gam2 * (-1.0 + gam2) * m4 *
    std::real(-( e2 * Qf * Qi * e2PgQfQiPlusf2Pzvfvi * PgC )
        + PzC * ( af2 * f4 * Pz * ( ai2 + vi * viC )
        - f2 * vfC * (
                ai2 * f2 * Pz * vf
            + e2PgQfQiPlusf2Pzvfvi * viC
            )
        )
    ) * sth2;


    RSM[0][2] = 4.0 * gam3 * m4 *
    std::real( e2 * Qf * Qi * PgC *
        ( af * ai * f2 * Pz * V
            - 2.0 * gam2 * (-1.0 + V*V)
            * e2PgQfQiPlusf2Pzvfvi * cth )
        +  f2 * PzC *
           ( af * ( ai * V * (e2*Pg*Qf*Qi + f2*Pz*vf*(vi + viC) )
                + 2.0 * af * f2 * Pz * ( 1.0 + gam2 * (-1.0 + V2) )
                * ( ai2 + vi * viC ) * cth )
        +
        vfC * (
            ai * f2 * Pz * ( af * V * vi
                            - 2.0 * ai * gam2 * (-1.0 + V2) * vf * cth )
            +
            viC * ( af * ai * f2 * Pz * V
                            - 2.0 * gam2 * (-1.0 + V2)
                                * e2PgQfQiPlusf2Pzvfvi
                                * cth )
        )
        )
    ) * sth;

    RSM[2][0] = RSM[0][2];

    RSM[1][2] = 2.0 * af * f2 * gam3 * m4 * V *
    std::imag( e2 * Pz * Qf * Qi * vi * PgC
        + PzC * ( -(ai2 * f2 * Pz * vf)
            - e2PgQfQiPlusf2Pzvfvi * viC
            + f2 * Pz * vfC * (ai2 + vi*viC)
        )
    ) * s2th;

    RSM[2][1] = RSM[1][2];

    RSM[2][2] = 2.0 * gam2 * m4 *
    std::real( -( e4 * Pg * (Qf2) * (Qi2) * PgC *
        ( 1.0 - 3.0*gam2
            + ( -1.0 + 3.0*gam2 + 4.0*gam4 * (-1.0 + V2) ) * c2th ) )
        - e2 * f2 * Qf * Qi * ( Pz * PgC *
            ( -4.0*af*ai*gam2*V*cth
                + vf*vi * ( 1.0 - 3.0*gam2
                            + ( -1.0 + 3.0*gam2 + 4.0*gam4 * (-1.0 + V2) )
                            * c2th ) )
        + Pg * PzC *
            ( -4.0*af*ai*gam2*V*cth
                + vfC*viC *
                ( 1.0 - 3.0*gam2
                    + ( -1.0 + 3.0*gam2 + 4.0*gam4 * (-1.0 + V2) )
                    * c2th ) ) )
        + f4 * Pz * PzC * (
            af * ( 4.0*ai*gam2*V*vf * ( vi + viC ) * cth
                + af * ( ai2 + vi*viC ) *
                    ( 1.0 + gam2 * ( -1.0 + 4.0*V2 )
                    + ( -1.0 + 5.0*gam2 + 4.0*gam4 * ( -1.0 + V2 ) )
                        * c2th ) )
        + vfC * (
                viC * ( 4.0*af*ai*gam2*V*cth
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
    std::real( e2 * Pz * Qf * Qi * PgC *
        ( 2.0*ai*vf + af*V*vi*cth )
        + PzC * ( af*V *
            ( ai2 * f2 * Pz * vf
                + e2PgQfQiPlusf2Pzvfvi * viC )
            * cth
        + vfC *
            ( 2.0*ai * ( e2*Pg*Qf*Qi + f2*Pz*vf*( vi + viC ) )
                + af * f2 * Pz * V * ( ai2 + vi*viC ) * cth )
        )
    ) * sth;

    RSM[3][0] = RSM[0][3];

    RSM[1][3] = -4.0 * af * ai * f2 * gam3 * m4 * V *
    std::imag( e2 * Pz * Qf * Qi * PgC
        + PzC * ( -(e2 * Pg * Qf * Qi)
            - f2 * Pz * vf * (vi + viC)
            + f2 * Pz * vfC * (vi + viC)
        )
    ) * sth;

    RSM[3][1] = RSM[1][3];

    RSM[2][3] = -4.0 * f2 * gam2 * m4 *
    std::real( e2 * gam2 * Pz * Qf * Qi * PgC *
        ( 2.0*ai*vf*cth + af*V*vi*(1.0 + cth2) )
        + (PzC * ( viC *
          ( 4.0*af2*ai*f2 * (-1.0 + gam2) * Pz * cth
            + 4.0*ai*f2*gam2 * Pz * vf * vfC * cth
            + af*gam2*V *
                ( e2*Pg*Qf*Qi + f2*Pz*vi*( vf + vfC ) ) *
                ( 3.0 + c2th ) )
        + ai * ( gam2 * vfC *
                ( 4.0*e2PgQfQiPlusf2Pzvfvi * cth
                + af*ai*f2*Pz*V * ( 3.0 + c2th ) )
            + af*f2*Pz *
                ( 4.0*af*(-1.0 + gam2) * vi * cth
                + ai*gam2*V * vf * ( 3.0 + c2th ) )
        )
        )
    ) / 2.0 );

    RSM[3][2] = RSM[2][3];

    RSM[3][3] = 2.0 * gam2 * m4 *
    std::real( e4 * Pg * Qf2 * Qi2 * PgC *
        ( 1.0 + 3.0*gam2 + (-1.0 + gam2) * c2th )
        + e2 * f2 * Qf * Qi * ( Pz * PgC *
            ( 4.0*af*ai*gam2*V*cth
                + vf*vi * ( 1.0 + 3.0*gam2 + (-1.0 + gam2) * c2th ) )
        + Pg * PzC *
            ( 4.0*af*ai*gam2*V*cth
                + vfC*viC *
                ( 1.0 + 3.0*gam2 + (-1.0 + gam2) * c2th ) ) )
        - f4 * Pz * PzC *
        ( -( af * ( 4.0*ai*gam2*V*vf * ( vi + viC ) * cth
                    + af * (-1.0 + gam2) * ( ai2 + vi*viC )
                    * ( 3.0 + c2th ) ) )
        - vfC *
            ( ai * ( 4.0*af*gam2*V*vi*cth
                    + ai*vf * ( 1.0 + 3.0*gam2 + (-1.0 + gam2) * c2th ) )
            + viC *
                ( 4.0*af*ai*gam2*V*cth
                + vf*vi * ( 1.0 + 3.0*gam2 + (-1.0 + gam2) * c2th ) )
            )
        )
    );


    // CONTRIBUTION LINEAR IN DIPOLE MOMENTS A, B, X, Y

        // RDM entries
    RDM[0][0] = 8.0 * gam4 * m4 *
    std::real( e2 * Qf * Qi *
        ( e2 * Pg * Qf * Qi * ( A + AC )
            + f2 * Pz * vi * ( X + vf * AC ) ) * PgC
        + f2 * PzC *
        ( vfC * ( ai2 * f2 * Pz * X
                + Ae2PgQfQiPlusf2PzvfviX * viC )
            + ( ai2 * f2 * Pz * vf
                + e2PgQfQiPlusf2Pzvfvi * viC )
            * XC ) ) * sth2;

    RDM[0][1] = 4.0 * gam4 * m4 * V *
    std::real( e2 * Qf * Qi *
        ( e2 * Pg * Qf * Qi * ( B + BC )
            + f2 * Pz * vi * ( Y
                            - std::complex<double>(0.0, 1.0) * af * AC
                            + vf * BC ) ) * PgC
        +  f2 * PzC *
        ( vfC *  ( ai2 * f2 * Pz * Y
                + Be2PgQfQiPlusf2PzvfviY * viC )
            + std::complex<double>(0.0, 1.0) * af * ( viC *
                    ( A * e2 * Pg * Qf * Qi
                    + f2 * Pz * vi * ( X - XC ) )
                + ai2 * f2 * Pz * ( X - XC ) )
            + ( ai2 * f2 * Pz * vf
                + e2PgQfQiPlusf2Pzvfvi * viC )
                * YC ) ) * sth2;

    RDM[1][0] =  -4.0 * gam4 * m4 * V *
    std::real( e2 * Qf * Qi *
        ( e2 * Pg * Qf * Qi * ( B + BC )
            + f2 * Pz * vi * ( Y
                                + std::complex<double>(0.0, 1.0) * af * AC
                                + vf * BC ) ) * PgC
        + f2 * PzC *
        ( vfC *
            ( ai2 * f2 * Pz * Y
                + Be2PgQfQiPlusf2PzvfviY * viC )
            - std::complex<double>(0.0, 1.0) * af *  ( viC *
                    ( A * e2 * Pg * Qf * Qi
                    + f2 * Pz * vi * ( X - XC ) )
                + ai2 * f2 * Pz * ( X - XC ) )
            + ( ai2 * f2 * Pz * vf
                + e2PgQfQiPlusf2Pzvfvi * viC )
                * YC ) ) * sth2;

    RDM[1][1] = 0.0;

    RDM[0][2] =  4.0 * (gam4 * gam) * m4 *
    std::real( -( e2 * Qf * Qi * PgC *
        ( ( -std::complex<double>(0.0, 1.0) ) * ai * f2 * Pz * V *
            ( Y - std::complex<double>(0.0, 1.0) * af * AC
                - vf * BC )
            + ( e2 * Pg * Qf * Qi * ( -2.0 + V2 ) * ( A + AC )
            + f2 * Pz * vi *
                ( ( -2.0 + V2 ) * ( X + vf * AC )
                    + std::complex<double>(0.0, 1.0) * af * (V2) * BC ))
                    * cth ) )
            + f2 * PzC *
            ( vfC * ( ai *
                ( std::complex<double>(0.0, 1.0) * V *
                    Be2PgQfQiPlusf2PzvfviY
                - ai * f2 * Pz * ( -2.0 + V2 ) * X * cth )
                + viC *
                ( std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V * Y
                    - ( -2.0 + V2 ) * Ae2PgQfQiPlusf2PzvfviX *
                        cth ))
            + viC *
            ( ( -std::complex<double>(0.0, 1.0) ) * ai * f2 * Pz * V * vf * YC
                - ( -2.0 + V2 ) * e2PgQfQiPlusf2Pzvfvi *
                    XC * cth
                + af * V *  ( ai * f2 * Pz * ( X + XC )
                    + std::complex<double>(0.0, 1.0) * V *
                        ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y - YC ) )
                        * cth ) )
            + ai * ( ( -std::complex<double>(0.0, 1.0) ) * V *
                e2PgQfQiPlusf2Pzvfvi * YC
                + ai * f2 * Pz * ( 2.0 - V2 ) * vf * XC * cth
                + af * V * ( A * e2 * Pg * Qf * Qi
                    + f2 * Pz * vi * ( X + XC )
                    + std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V
                        * ( Y - YC ) * cth ) ) ) ) * sth;

    RDM[2][0] = 4.0 * (gam4 * gam) * m4 *
    std::real( e2 * Qf * Qi * PgC *
        ( ai * f2 * Pz * V * ( af * AC
                                - std::complex<double>(0.0, 1.0) * ( Y - vf * BC ) )
            - ( e2 * Pg * Qf * Qi * ( -2.0 + V2 ) * ( A + AC )
                + f2 * Pz * vi * ( ( -2.0 + V2 ) * ( X + vf * AC )
                - std::complex<double>(0.0, 1.0) * af * (V2) * BC ) ) * cth )
        + f2 * PzC *
        ( vfC * (
                -( ai * ( std::complex<double>(0.0, 1.0) * V
                            * Be2PgQfQiPlusf2PzvfviY
                        + ai * f2 * Pz * ( -2.0 + V2 ) * X * cth ) )
                + viC *
                    ( -std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V * Y
                    - ( -2.0 + V2 ) * Ae2PgQfQiPlusf2PzvfviX
                        * cth ) )
            + viC *
            ( std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V * vf * YC
                - ( -2.0 + V2 ) * e2PgQfQiPlusf2Pzvfvi
                    * XC * cth
                + af * V * ( ai * f2 * Pz * ( X + XC )
                    - std::complex<double>(0.0, 1.0) * V
                        * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y - YC ) )
                        * cth ) )
            +ai * ( std::complex<double>(0.0, 1.0) * V
                * e2PgQfQiPlusf2Pzvfvi * YC
                + ai * f2 * Pz * ( 2.0 - V2 ) * vf * XC * cth
                + af * V * ( A * e2 * Pg * Qf * Qi
                    + f2 * Pz * vi * ( X + XC )
                    - std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V
                        * ( Y - YC ) * cth ) ) ) ) * sth;


    RDM[1][2] =  -4.0 * (gam4 * gam) * m4 * V *
    std::imag( e2 * Qf * Qi * PgC *
        ( ai * f2 * Pz * V * ( X - vf * AC
        + std::complex<double>(0.0, 1.0) * af * BC )
        + std::complex<double>(0.0, 1.0) * ( e2 * Pg * Qf * Qi * ( B + BC )
        + f2 * Pz * vi * ( Y + std::complex<double>(0.0, 1.0) * af * AC
        + vf * BC ) ) * cth )
        + f2 * PzC * ( vfC *
            ( ai * ( A * e2 * Pg * Qf * Qi * V + f2 * Pz * V * vi * X
        + std::complex<double>(0.0, 1.0) * ai * f2 * Pz * Y * cth )
        + viC * ( ai * f2 * Pz * V * X
        + std::complex<double>(0.0, 1.0) * Be2PgQfQiPlusf2PzvfviY
                                    * cth ) )
        + std::complex<double>(0.0, 1.0) * ( viC *
            ( std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V * vf * XC
        + e2PgQfQiPlusf2Pzvfvi * YC * cth
        + af * ( ai * f2 * Pz * V * ( Y + YC )
        - std::complex<double>(0.0, 1.0) * ( A * e2 * Pg * Qf * Qi
        + f2 * Pz * vi * ( X - XC ) ) * cth ) )
        + ai * ( std::complex<double>(0.0, 1.0) * V * e2PgQfQiPlusf2Pzvfvi
                        * XC
        + ai * f2 * Pz * vf * YC * cth
        + af * ( V * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y + YC ) )
        - std::complex<double>(0.0, 1.0) * ai * f2 * Pz * ( X - XC )
                                    * cth ) ) ) ) ) * sth;


    RDM[2][1] =  4.0 * gam5 * m4 * V *
    std::real( e2 * Qf * Qi * PgC *
        (
            ai * f2 * Pz * V *
            ( std::complex<double>(0.0, 1.0) * ( X - vf * AC ) + af * BC )
            +
            ( e2 * Pg * Qf * Qi * ( B + BC )
            + f2 * Pz * vi * ( Y - std::complex<double>(0.0, 1.0) * af * AC + vf * BC ) )
            * cth
        )
        +
        f2 * PzC *
        (
            vfC *
            (
                ai * ( std::complex<double>(0.0, 1.0) * V * Ae2PgQfQiPlusf2PzvfviX
                    + ai * f2 * Pz * Y * cth )
                +
                viC *
                ( std::complex<double>(0.0, 1.0) * ai * f2 * Pz * V * X
                    + Be2PgQfQiPlusf2PzvfviY * cth )
            )
            +
            viC *
            (
                ( -std::complex<double>(0.0, 1.0) ) * ai * f2 * Pz * V * vf * XC
                + e2PgQfQiPlusf2Pzvfvi * YC * cth
                + af * ( ai * f2 * Pz * V * ( Y + YC )
                        + std::complex<double>(0.0, 1.0)
                            * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X - XC ) )
                            * cth )
            )
            +
            ai *
            (
                ( -std::complex<double>(0.0, 1.0) ) * V * e2PgQfQiPlusf2Pzvfvi * XC
                + ai * f2 * Pz * vf * YC * cth
                + af * ( V * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y + YC ) )
                        + std::complex<double>(0.0, 1.0) * ai * f2 * Pz * ( X - XC ) * cth )
            )
        )
    ) * sth;

    RDM[2][2] =
    -8.0 * gam6 * m4 * (-1.0 + V2) * cth *
    std::real(
        e2 * Qf * Qi * PgC *
        ( Ae2PgQfQiPlusf2PzvfviX * cth
            + AC * ( af * ai * f2 * Pz * V
                            + e2PgQfQiPlusf2Pzvfvi * cth ) )
        +
        f2 * PzC *
        ( ai * ( af * V * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X + XC ) )
                + ai * f2 * Pz * X * vfC * cth
                + ai * f2 * Pz * vf * XC * cth )
            + viC *
            ( af * ai * f2 * Pz * V * ( X + XC )
                + ( Ae2PgQfQiPlusf2PzvfviX * vfC
                    + e2PgQfQiPlusf2Pzvfvi * XC )
                * cth ) )
    );


    RDM[0][3] =
    -4.0 * gam3 * m4 *
    std::real(
        e2 * Qf * Qi * PgC *
        (
            ai * f2 * Pz * ( (1.0 + gam2) * ( X + vf * AC )
                            - std::complex<double>(0.0, 1.0) * af * (-1.0 + gam2) * BC )
            + std::complex<double>(0.0, 1.0) * gam2 * V *
                ( e2 * Pg * Qf * Qi * ( B - BC )
                + f2 * Pz * vi * ( Y - std::complex<double>(0.0, 1.0) * af * AC - vf * BC ) )
                * cth
        )
        +
        f2 * PzC *
        (
            vfC *
            (
                ai * ( (1.0 + gam2) * Ae2PgQfQiPlusf2PzvfviX
                    + std::complex<double>(0.0, 1.0) * ai * f2 * gam2 * Pz * V * Y * cth )
                + viC *
                    ( ai * f2 * (1.0 + gam2) * Pz * X
                    + std::complex<double>(0.0, 1.0) * gam2 * V
                        * Be2PgQfQiPlusf2PzvfviY
                        * cth )
            )
            + ai *
            (
                (1.0 + gam2) * e2PgQfQiPlusf2Pzvfvi * XC
                - std::complex<double>(0.0, 1.0) * ai * f2 * gam2 * Pz * V * vf * YC * cth
                + af * ( std::complex<double>(0.0, 1.0) * (-1.0 + gam2)
                        * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y - YC ) )
                        + ai * f2 * gam2 * Pz * V * ( X + XC ) * cth )
            )
            + viC *
            (
                ai * f2 * (1.0 + gam2) * Pz * vf * XC
                - std::complex<double>(0.0, 1.0) * gam2 * V
                    * e2PgQfQiPlusf2Pzvfvi * YC * cth
                + af * ( std::complex<double>(0.0, 1.0) * ai * f2 * (-1.0 + gam2) * Pz * ( Y - YC )
                        + gam2 * V * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X + XC ) )
                        * cth )
            )
        )
    ) * sth;


    RDM[3][0] =
    -4.0 * gam3 * m4 *
    std::real(
        e2 * Qf * Qi * PgC *
        (
            ai * f2 * Pz * ( (1.0 + gam2) * ( X + vf * AC )
                            + std::complex<double>(0.0, 1.0) * af * (-1.0 + gam2) * BC )
            - std::complex<double>(0.0, 1.0) * gam2 * V *
                ( e2 * Pg * Qf * Qi * ( B - BC )
                + f2 * Pz * vi * ( Y + std::complex<double>(0.0, 1.0) * af * AC - vf * BC ) )
                * cth
        )
        +
        f2 * PzC *
        (
            vfC *
            (
                ai * ( (1.0 + gam2) * Ae2PgQfQiPlusf2PzvfviX
                    - std::complex<double>(0.0, 1.0) * ai * f2 * gam2 * Pz * V * Y * cth )
                + viC *
                    ( ai * f2 * (1.0 + gam2) * Pz * X
                    - std::complex<double>(0.0, 1.0) * gam2 * V
                        * Be2PgQfQiPlusf2PzvfviY
                        * cth )
            )
            + ai *
            (
                (1.0 + gam2) * e2PgQfQiPlusf2Pzvfvi * XC
                + std::complex<double>(0.0, 1.0) * ai * f2 * gam2 * Pz * V * vf * YC * cth
                + af * ( -std::complex<double>(0.0, 1.0) * (-1.0 + gam2)
                        * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y - YC ) )
                        + ai * f2 * gam2 * Pz * V * ( X + XC ) * cth )
            )
            + viC *
            (
                ai * f2 * (1.0 + gam2) * Pz * vf * XC
                + std::complex<double>(0.0, 1.0) * gam2 * V
                    * e2PgQfQiPlusf2Pzvfvi * YC * cth
                + af * ( -std::complex<double>(0.0, 1.0) * ai * f2 * (-1.0 + gam2) * Pz * ( Y - YC )
                        + gam2 * V * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X + XC ) )
                        * cth )
            )
        )
    ) * sth;


    RDM[1][3] =
    4.0 * gam3 * m4 *
    std::real(
        e2 * Qf * Qi * PgC *
        ( ai * f2 * gam2 * Pz * V *
            ( Y + std::complex<double>(0.0, 1.0) * af * AC + vf * BC )
            - std::complex<double>(0.0, 1.0) * (-1.0 + gam2) *
            ( e2 * Pg * Qf * Qi * ( A - AC )
                + f2 * Pz * vi * ( X - vf * AC
                                    + std::complex<double>(0.0, 1.0) * af * BC ) )
            * cth
        )
        +
        f2 * PzC *
        ( vfC *
            ( ai * ( gam2 * V * Be2PgQfQiPlusf2PzvfviY
                    - std::complex<double>(0.0, 1.0) * ai * f2 * (-1.0 + gam2) * Pz * X * cth )
                + viC *
                    ( ai * f2 * gam2 * Pz * V * Y
                    - std::complex<double>(0.0, 1.0) * (-1.0 + gam2)
                        * Ae2PgQfQiPlusf2PzvfviX * cth )
            )
            + ai *
            ( gam2 * V * e2PgQfQiPlusf2Pzvfvi * YC
                + std::complex<double>(0.0, 1.0) * ai * f2 * (-1.0 + gam2) * Pz * vf * XC * cth
                + af * ( -std::complex<double>(0.0, 1.0) * gam2 * V
                        * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X - XC ) )
                        + ai * f2 * (-1.0 + gam2) * Pz * ( Y + YC ) * cth )
            )
            + viC *
            ( std::complex<double>(0.0, 1.0) *
                ( -std::complex<double>(0.0, 1.0) * ai * f2 * gam2 * Pz * V * vf * YC
                    + (-1.0 + gam2) * e2PgQfQiPlusf2Pzvfvi
                        * XC * cth )
                + af * ( -std::complex<double>(0.0, 1.0) * ai * f2 * gam2 * Pz * V * ( X - XC )
                        + (-1.0 + gam2)
                            * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y + YC ) )
                            * cth )
            )
        )
    ) * sth;


    RDM[3][1] =
    4.0 * gam3 * m4 *
    std::real(
        e2 * Qf * Qi * PgC *
        (
            -( ai * f2 * gam2 * Pz * V * ( Y - std::complex<double>(0.0,1.0) * af * AC + vf * BC ) )
            - std::complex<double>(0.0,1.0) * (-1.0 + gam2) *
                ( e2 * Pg * Qf * Qi * ( A - AC )
                + f2 * Pz * vi * ( X - vf * AC - std::complex<double>(0.0,1.0) * af * BC ) )
                * std::cos(theta)
        )
        - f2 * PzC *
        (
            vfC *
            (
                ai * ( gam2 * V * Be2PgQfQiPlusf2PzvfviY
                    + std::complex<double>(0.0,1.0) * ai * f2 * (-1.0 + gam2) * Pz * X * cth )
                + viC *
                    ( ai * f2 * gam2 * Pz * V * Y
                    + std::complex<double>(0.0,1.0) * (-1.0 + gam2)
                        * Ae2PgQfQiPlusf2PzvfviX
                        * cth )
            )
            + ai *
            ( gam2 * V * e2PgQfQiPlusf2Pzvfvi * YC
                - std::complex<double>(0.0,1.0) * ai * f2 * (-1.0 + gam2) * Pz * vf * XC * cth
                + af * ( std::complex<double>(0.0,1.0) * gam2 * V
                        * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X - XC ) )
                        + ai * f2 * (-1.0 + gam2) * Pz * ( Y + YC ) * cth )
            )
            + viC *
            (
                ( -std::complex<double>(0.0,1.0) ) *
                ( std::complex<double>(0.0,1.0) * ai * f2 * gam2 * Pz * V * vf * YC
                    + (-1.0 + gam2) * e2PgQfQiPlusf2Pzvfvi
                        * XC * cth )
                + af * ( std::complex<double>(0.0,1.0) * ai * f2 * gam2 * Pz * V * ( X - XC )
                        + (-1.0 + gam2)
                            * ( B * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( Y + YC ) )
                            * cth )
            )
        )
    ) * sth;



    RDM[2][3] =
    -2.0 * gam4 * m4 *
    std::real(
        f2 * PzC *
        (
            -( af * V *
                ( ai2 * f2 * Pz * ( X + XC )
                + viC * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X + XC ) ) )
            * ( -2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) )
            + std::complex<double>(0.0, 1.0) *
            (
                std::complex<double>(0.0, -4.0) * ai *
                ( e2 * Pg * Qf * Qi + f2 * Pz * vf * ( vi + viC ) )
                * XC * std::cos(theta)
                + V * ( ai2 * f2 * Pz * vf
                        + e2PgQfQiPlusf2Pzvfvi * viC )
                    * YC
                    * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th )
            )
            + vfC *
            (
                std::complex<double>(0.0, 1.0) * ai *
                ( std::complex<double>(0.0, -4.0) * Ae2PgQfQiPlusf2PzvfviX * cth
                    + ai * f2 * Pz * V * Y
                        * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) )
                + viC *
                ( 4.0 * ai * f2 * Pz * X * cth
                    - std::complex<double>(0.0, 1.0) * V
                        * Be2PgQfQiPlusf2PzvfviY
                        * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) )
            )
        )
        + e2 * Qf * Qi * PgC *
        (
            4.0 * ai * f2 * Pz * ( X + vf * AC ) * cth
            - std::complex<double>(0.0, 1.0) * V *
            ( e2 * Pg * Qf * Qi * ( B - BC )
                * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th )
                + f2 * Pz * vi *
                    ( std::complex<double>(0.0, -1.0) * af * AC
                        * ( -2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th )
                    + ( Y - vf * BC )
                        * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) ) )
        )
    );


    RDM[3][2] =
    -2.0 * gam4 * m4 *
    std::real(
        f2 * PzC *
        (
            -( af * V *
                ( ai2 * f2 * Pz * ( X + XC )
                + viC * ( A * e2 * Pg * Qf * Qi + f2 * Pz * vi * ( X + XC ) ) )
            * ( -2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) )
            - std::complex<double>(0.0, 1.0) *
            (
                std::complex<double>(0.0, 4.0) * ai *
                ( e2 * Pg * Qf * Qi + f2 * Pz * vf * ( vi + viC ) )
                * XC * cth
                + V * ( ai2 * f2 * Pz * vf
                        + e2PgQfQiPlusf2Pzvfvi * viC )
                    * YC
                    * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th )
            )
            + vfC *
            (
                std::complex<double>(0.0, 1.0) * ai *
                ( std::complex<double>(0.0, -4.0) * Ae2PgQfQiPlusf2PzvfviX
                    * cth
                    + ai * f2 * Pz * V * Y
                        * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) )
                + viC *
                ( 4.0 * ai * f2 * Pz * X * cth
                    + std::complex<double>(0.0, 1.0) * V
                        * Be2PgQfQiPlusf2PzvfviY
                        * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) )
            )
        )
        + e2 * Qf * Qi * PgC *
        (
            4.0 * ai * f2 * Pz * ( X + vf * AC ) * cth
            + std::complex<double>(0.0, 1.0) * V *
            ( e2 * Pg * Qf * Qi * ( B - BC )
                * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th )
                + f2 * Pz * vi *
                    ( std::complex<double>(0.0, 1.0) * af * AC
                        * ( -2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th )
                    + ( Y - vf * BC )
                        * ( 2.0 + gam2 * ( -1.0 + V2 ) + gam2 * ( -1.0 + V2 ) * c2th ) ) )
        )
    );


    RDM[3][3] =
    8.0 * gam4 * m4 *
    std::real(
        e4 * Pg * Qf2 * Qi2 * ( A + AC ) * PgC
        +
        e2 * f2 * Qf * Qi * (
            Pz * vi * X * PgC
        + A * Pg * PzC * vfC * viC
        + Pg * PzC * viC * XC
        + A * af * ai * Pg * V * PzC * cth
        + Pz * AC * PgC * ( vf * vi + af * ai * V * cth )
        )
        +
        f4 * Pz * PzC * (
            ai * ( ai * X * vfC
                + ai * vf * XC
                + af * V * vi * X * cth
                + af * V * vi * XC * cth )
        + viC * (
                vi * ( X * vfC + vf * XC )
            + af * ai * V * ( X + XC ) * cth
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
