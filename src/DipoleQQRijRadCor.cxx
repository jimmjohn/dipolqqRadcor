#include "DipoleQQRijRadCor.h"
#include <cmath>
#include <iostream>
#include <iomanip>
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
    int KFi, KFf, IfGSW, KeyElw, IfPrint;


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
    cosThetaD = cth;
    IFGSW = 1;
    IfPrint = 0; // Switch for printing debug information
    KKdizet kkDizet;
    kkDizet.ReadEWtabs();  // Read the EW tables from file
    if(IFGSW==1) {
        kkDizet.InterpoGSW(KFi, KFf, Svar, cosThetaD);
        kkDizet.GetGSWxy(re, im);
        std::transform(re, re + 7, im,  GSW, [](double r, double i) { return std::complex<double>(r, i); });
        RhoEW         = GSW[1];
        //VPgamma       = GSW[6];
        VPgamma       = 1.0/(2.0- GSW[6]);
        CorEle        = GSW[2];
        CorFin        = GSW[3];
        CorEleFin     = GSW[4];
        sw2           = kkDizet.D_swsq;//0.23113; //sin^2(theta_W)
        Gz_           = Gz_*Svar/(Mz_*Mz_); // Running width
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
    t   = pow(m,2) - 0.5 * s * (1.0 - V * cth)

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
    double gam4 = pow(gam, 4);
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
                    (m_sw2*(1.0-m_sw2))*16.;
            Fvivf = Fvivf * (m_sw2/(1.0-m_sw2)) //! why this factor? Because  vi ai couplings in KKMC
                 // are divided by deno= 4 sqrt(m_swsq*(1d0-m_swsq)) in addition multiplied by 2 and we are not using it here
                 // Ve = (2*T3e -4*Qe*m_swsq)/deno and Ae = 2*T3e/deno.
        }
        //Effective photon propagator with v_{if}-v_i*v_f correction to Z exchange
        //  (NOTE: Pz reduces to 1/s if  rad. corr. are off):
        Pg = Gammavp/s + (Fvivf * s/Den) * (Rho11*(K11-K1*K1));
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

    RSM[0][0] = ((-16.0*a0_4*f4*gam4*(-1.0+gam2)*m4
                + 32.0*a0_2*f4*gam4*m4*v0_2
                + (1.0+gam2)*(e4*Mz_4
                + 16.0*gam4*m4*pow(e2+f2*v0_2,2)
                + Mz_2*(e4*(Gz_2-8.0*gam2*m2)-8.0*e2*f2*gam2*m2*v0_2)))
                *sth2)
                /(4.0*gam2*denom);

    RSM[0][1] = (-2.0*a0_*e2*f2*gam2*Gz_*m2*Mz_*V*v0_*sth2)/denom;

    RSM[1][0] = RSM[0][1];

    RSM[1][1] = -((-1.0 + gam2)*(-16.0*a0_4*f4*gam4*m4
                + e4*Mz_4 + 16.0*gam4*m4*pow(e2 + f2*v0_2,2)
                + Mz_2*(e4*(Gz_2 - 8.0*gam2*m2)
                - 8.0*e2*f2*gam2*m2*v0_2))*sth2)
                / (4.0*gam2*denom);

    RSM[0][2] = ((e4*Mz_4 + 16.0*gam4*m4*pow(e2+f2*v0_2,2)
                + Mz_2*(e4*(Gz_2-8.0*gam2*m2)
                - 8.0*e2*f2*gam2*m2*v0_2))*cth
                + 4.0*a0_2*f2*gam2*m2*(-(e2*Mz_2*V)
                + 4.0*gam2*m2*(e2*V + f2*v0_2*(2.0*V+cth))))*sth
                / (2.0*gam*denom);

    RSM[2][0] = RSM[0][2];

    RSM[1][2] = (-a0_*e2*f2*gam*Gz_*m2*Mz_*V*v0_*s2th)/denom;

    RSM[2][1] = RSM[1][2];

    RSM[2][2] = (8.0*a0_4*f4*gam4*(-1.0+gam2)*m4*(3.0+c2th)
                + ((e4*Mz_4 + 16.0*gam4*m4*pow(e2+f2*v0_2,2)
                + Mz_2*(e4*(Gz_2-8.0*gam2*m2)
                - 8.0*e2*f2*gam2*m2*v0_2))
                *(-1.0+3.0*gam2+(1.0+gam2)*c2th))/2.0
                + 16.0*a0_2*f2*gam4*m2*(-(e2*Mz_2*V*cth)
                + m2*(4.0*e2*gam2*V*cth
                + f2*v0_2*(-2.0+3.0*gam2+8.0*gam2*V*cth+gam2*c2th))))
                /(4.0*gam2*denom);

    RSM[0][3] = (-2.0*a0_*f2*gam*m2*v0_
                *(4.0*a0_2*f2*gam2*m2*V*cth
                + (-(e2*Mz_2) + 4.0*gam2*m2*(e2 + f2*v0_2))
                *(2.0 + V*cth))*sth)/denom;

    RSM[3][0] = RSM[0][3];

    RSM[1][3] = (2.0*a0_2*e2*f2*gam*Gz_*m2*Mz_*V*sth)/denom;

    RSM[3][1] = RSM[1][3];

    RSM[2][3] = -((a0_*f2*gam2*m2*v0_
                *((-(e2*Mz_2) + 4.0*gam2*m2*(e2 + f2*v0_2))
                *(4.0*cth + V*(3.0 + c2th)) + 4.0*a0_2*f2*m2
                *(4.0*(-1.0 + gam2)*cth + gam2*V*(3.0 + c2th))))/denom);

    RSM[3][2] = RSM[2][3];

    RSM[3][3] = (e4*(1.0 + gam2 + (-1.0 + gam2)*cth2)
                + (4.0*e2*f2*gam2*m2*(4.0*gam2*m2 - Mz_2)
                *(4.0*a0_2*gam2*V*cth + v0_2*(1.0 + 3.0*gam2 + (-1.0 + gam2)*c2th)))
                /denom + (8.0*f4*gam4*m4
                *(a0_4*(-1.0 + gam2)*(3.0 + c2th) + v0_4*(1.0 + 3.0*gam2 + (-1.0 + gam2)*c2th)
                + 2.0*a0_2*v0_2*(-1.0 + 3.0*gam2 + 8.0*gam2*V*cth + (-1.0 + gam2)*c2th)))
                /denom) / (4.0*gam2);



    // CONTRIBUTION LINEAR IN DIPOLE MOMENTS A, B, X, Y

        // RDM entries
    RDM[0][0] = ((e2*(e2*Mz_4 + 16.0*gam4*m4*(e2 + f2*v0_2)
                + Mz_2*(e2*(Gz_2 - 8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2)) * ReA1
                + 4.0*f2*gam2*m2*v0_*(e2*Gz_*Mz_*(-v0_*ImA1 + ImX)
                + (4.0*a0_2*f2*gam2*m2 - e2*Mz_2
                + 4.0*gam2*m2*(e2 + f2*v0_2)) * ReX))*sth2) / denom;

    RDM[0][1] = (V*(e2*(e2*Mz_4
                + 16.0*gam4*m4*(e2 + f2*v0_2)
                + Mz_2*(e2*(Gz_2 - 8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2)) * ReB
                - 4.0*f2*gam2*m2
                *(4.0*a0_*f2*gam2*m2*(a0_2 + v0_2) * ImX
                + v0_*(e2*(a0_*(4.0*gam2*m2 - Mz_2) * ImA1
                + Gz_*Mz_*(v0_*ImB - ImY + a0_*ReA1))
                + (-4.0*a0_2*f2*gam2*m2 + e2*Mz_2
                -4.0*gam2*m2*(e2 + f2*v0_2)) * ReY)))*sth2)/(2.0*denom);

    RDM[1][0] = -(V * (e2*(e2*Mz_4+16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2 - 8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2)) * ReB
                + 4.0*f2*gam2*m2*(4.0*a0_*f2*gam2*m2*(a0_2 + v0_2) * ImX
                + v0_*(e2*(a0_*(4.0*gam2*m2 - Mz_2) * ImA1
                + Gz_*Mz_*(-(v0_*ImB) + ImY + a0_*ReA1))+(4.0*a0_2*f2*gam2*m2-e2*Mz_2
                + 4.0*gam2*m2*(e2 + f2*v0_2)) * ReY)))*sth2)/(2.0*denom);

    RDM[1][1] = 0.0;

    RDM[0][2] = -(gam * (e2*(4.0*a0_2*f2*gam2*m2
                *(-4.0*gam2*m2+Mz_2)*V + (-2.0 + V2) * (e2*Mz_4
                + 16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2 - 8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2)) * cth)*ReA1
                + 4.0*f2*gam2*m2*(a0_*V*(-(e2*Mz_2)
                + 4.0*a0_2*f2*gam2*m2*V*cth
                + 4.0*gam2*m2*(e2+f2*v0_2*(2.0+V*cth)))*ImY
                + v0_*((-2.0 + V*V)*(-(e2*Mz_2)+4.0*gam2*m2*(e2+f2*v0_2))*cth
                + 4.0*a0_2*f2*gam2*m2*(-2.0*V+(-2.0+V2)*cth))*ReX
                + e2*(a0_*(4.0*gam2*m2-Mz_2)*V*v0_*(1.0+V*cth)*ImB
                + Gz_*Mz_*((a0_2*V-(-2.0+V2)*v0_2*cth)*ImA1
                + v0_*((-2.0+V2)*cth*ImX + a0_*V*(1.0+V*cth)*ReB)
                - a0_*V*ReY))))*sth)/(2.0 * denom);

    RDM[2][0] = -(gam*(e2*(4.0*a0_2*f2*gam2*m2
                *(-4.0*gam2*m2+Mz_2)*V + (-2.0+V2)*(e2*Mz_4
                + 16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2-8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2))*cth)*ReA1
                - 4.0*f2*gam2*m2*(a0_*V*(-(e2*Mz_2)
                + 4.0*a0_2*f2*gam2*m2*V*cth
                + 4.0*gam2*m2*(e2+f2*v0_2*(2.0+V*cth)))*ImY
                - v0_*((-2.0 + V2)*(-(e2*Mz_2)
                + 4.0*gam2*m2*(e2+f2*v0_2))*cth
                + 4.0*a0_2*f2*gam2*m2*(-2.0*V+(-2.0+V2)*cth))*ReX
                + e2*(a0_*(4.0*gam2*m2-Mz_2)*V*v0_*(1.0+V*cth)*ImB
                + Gz_*Mz_*((-(a0_*a0_*V)+(-2.0+V2)*v0_2*cth)*ImA1
                + v0_*(-((-2.0+V2)*cth*ImX)+a0_*V*(1.0+V*cth)*ReB)
                - a0_*V*ReY))))*sth)/(2.0 * denom);


//Check this one more time
    RDM[1][2] = -(gam*V*(e2*(4.0*a0_2*f2*gam2*m2*(4.0*gam2*m2-Mz_2)*V
                + (e2*Mz_4+16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2-8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2))*cth)*ReB
                + 4.0*f2*gam2*m2*(a0_*(-(e2*Mz_2*V)
                + 4.0*a0_2*f2*gam2*m2*cth
                + 4.0*gam2*m2*(e2*V+f2*v0_2*(2.0*V+cth)))*ImX
                + e2*(a0_*(4.0*gam2*m2-Mz_2)*v0_*(V+cth)*ImA1
                + Gz_*Mz_*(-((a0_2*V+v0_2*cth)*ImB)+v0_*(cth*ImY+a0_*(V+cth)*ReA1)
                - a0_*V*ReX))+v0_*((-(e2*Mz_2)+4.0*gam2*m2*(e2+f2*v0_2))*cth
                + 4.0*a0_2*f2*gam2*m2*(2.0*V+cth))*ReY))*sth)/(2.0*denom);



    RDM[2][1] = (gam*V*(e2*(4.0*a0_2*f2*gam2*m2
                *(4.0*gam2*m2-Mz_2)*V
                + (e2*Mz_4+16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2-8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2))*cth)*ReB
                - 4.0*f2*gam2*m2*(a0_*(-(e2*Mz_2*V)
                + 4.0*a0_2*f2*gam2*m2*cth
                + 4.0*gam2*m2*(e2*V+f2*v0_2*(2.0*V+cth)))*ImX
                + e2*(a0_*(4.0*gam2*m2-Mz_2)*v0_*(V+cth)*ImA1
                + Gz_*Mz_*((a0_2*V+v0_2*cth)*ImB+v0_*(-(cth*ImY)
                + a0_*(V+cth)*ReA1)-a0_*V*ReX))-v0_*((-(e2*Mz_2)
                + 4.0*gam2*m2*(e2+f2*v0_2))*cth
                + 4.0*a0_2*f2*gam2*m2*(2.0*V+cth))*ReY))*sth)/(2.0*denom);

    RDM[2][2] = -(gam2*(-1.0+V2)*cth*(e2*(4.0*a0_2*f2*gam2*m2
                *(4.0*gam2*m2-Mz_2)*V
                + (e2*Mz_4+16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2-8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2))*cth)*ReA1
                + 4.0*f2*gam2*m2*(-(e2*Gz_*Mz_*((a0_2*V
                + v0_2*cth)*ImA1- v0_*cth*ImX))
                + v0_*((-(e2*Mz_2)+4.0*gam2*m2*(e2+f2*v0_2))*cth
                + 4.0*a0_2*f2*gam2*m2*(2.0*V+cth))*ReX)))/denom;


    RDM[0][3] = (gam*(e2*(4.0*a0_2*f2*(-1.0+gam2)*m2
                *(4.0*gam2*m2-Mz_2) + V*(e2*Mz_4
                + 16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2-8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2))*cth)*ImB
                - 4.0*f2*m2*(-(gam2*v0_*(V*(-(e2*Mz_2)
                + 4.0*gam2*m2*(e2+f2*v0_2))*cth
                + 4.0*a0_2*f2*m2*(-2.0+2.0*gam2+gam2*V*cth))*ImY)
                + a0_*(-(e2*(1.0+gam2)*Mz_2)+4.0*a0_2*f2*gam4*m2*V*cth
                + 4.0*gam2*m2*(e2*(1.0+gam2)
                + f2*v0_2*(2.0+2.0*gam2+gam2*V*cth)))*ReX
                + e2*(a0_*(4.0*gam2*m2-Mz_2)*v0_*(1.0+gam2+gam2*V*cth)*ReA1
                + Gz_*Mz_*(a0_*(1.0+gam2)*ImX +(a0_2*(1.0-gam2)
                - gam2*V*v0_2*cth)*ReB
                + v0_*(-(a0_*(1.0+gam2+gam2*V*cth)*ImA1)
                + gam2*V*cth*ReY)))))*sth)/(2.0*denom);

    RDM[3][0] = -(gam*(e2*(4.0*a0_2*f2*(-1.0+gam2)*m2*(4.0*gam2*m2-Mz_2)
                + V*(e2*Mz_4+16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2-8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2))*cth)*ImB
                + 4.0*f2*m2*(gam2*v0_*(V*(-(e2*Mz_2)
                + 4.0*gam2*m2*(e2+f2*v0_2))*cth
                + 4.0*a0_2*f2*m2*(-2.0+2.0*gam2+gam2*V*cth))*ImY
                + a0_*(-(e2*(1.0+gam2)*Mz_2)+4.0*a0_2*f2*gam4*m2*V*cth
                + 4.0*gam2*m2*(e2*(1.0+gam2)
                + f2*v0_2*(2.0+2.0*gam2+gam2*V*cth)))*ReX
                + e2*(a0_*(4.0*gam2*m2-Mz_2)*v0_*(1.0+gam2+gam2*V*cth)*ReA1
                + Gz_*Mz_*(a0_*(1.0+gam2)*ImX+(a0_2*(-1.0+gam2)+gam2*V*v0_2*cth)*ReB
                - v0_*(a0_*(1.0+gam2+gam2*V*cth)*ImA1+gam2*V*cth*ReY)))))*sth)/(2.0*denom);

    RDM[1][3] = (gam*V*(e2*(4.0*a0_2*f2*gam2*m2*(4.0*gam2*m2-Mz_2)
                + V*(e2*Mz_4+16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2-8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2))*cth)*ImA1
                + 4.0*f2*gam2*m2*(v0_*(V*(-(e2*Mz_2)
                + 4.0*gam2*m2*(e2+f2*v0_2))*cth
                + 4.0*a0_2*f2*gam2*m2*(2.0+V*cth))*ImX
                + e2*(a0_*(4.0*gam2*m2-Mz_2)*v0_*(1.0+V*cth)*ReB
                + Gz_*Mz_*(a0_*ImY+(a0_2+V*v0_2*cth)*ReA1-v0_*(a0_*(1.0+V*cth)*ImB+V*cth*ReX)))
                + a0_*(-(e2*Mz_2)+4.0*a0_2*f2*gam2*m2*V*cth
                + 4.0*gam2*m2*(e2+f2*v0_2*(2.0+V*cth)))*ReY))*sth)/(2.0*denom);

    RDM[3][1] = (gam*V*(e2*(4.0*a0_2*f2*gam2*m2*(4.0*gam2*m2-Mz_2)
                + V*(e2*Mz_4+16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2-8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2))*cth)*ImA1
                - 4.0*f2*gam2*m2*(-(v0_*(V*(-(e2*Mz_2)
                + 4.0*gam2*m2*(e2+f2*v0_2))*cth
                + 4.0*a0_2*f2*gam2*m2*(2.0+V*cth))*ImX)
                + e2*(a0_*(4.0*gam2*m2-Mz_2)*v0_*(1.0+V*cth)*ReB
                + Gz_*Mz_*(a0_*ImY-(a0_2+V*v0_2*cth)*ReA1
                + v0_*(-(a0_*(1.0+V*cth)*ImB)+V*cth*ReX)))
                + a0_*(-(e2*Mz_2) + 4.0*a0_2*f2*gam2*m2*V*cth
                + 4.0*gam2*m2*(e2+f2*v0_2*(2.0+V*cth)))*ReY))*sth)/(2.0*denom);

    RDM[2][3] = -(2.0*e2*V*(e2*Mz_4+16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2-8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2))*ImB*sth2
                + 4.0*f2*gam2*m2*(-4.0*a0_*(e2*Mz_2*cth
                - a0_2*f2*gam2*m2*V*(3.0+c2th)
                - gam2*m2*(4.0*e2*cth+f2*v0_2*(8.0*cth+V*(3.0+c2th))))*ReX
                + 2.0*V*v0_*(4.0*a0_2*f2*gam2*m2-e2*Mz_2
                + 4.0*gam2*m2*(e2+f2*v0_2))*ImY*sth2
                + e2*(a0_*(4.0*gam2*m2-Mz_2)*v0_*(4.0*cth+V*(3.0+c2th))*ReA1
                + 2.0*Gz_*Mz_*(2.0*a0_*cth*ImX+v0_*(-(a0_*(V+2.0*cth+V*cth2)*ImA1)
                + V*v0_*ReB*sth2-V*ReY*sth2)))))/(4.0*denom);

    RDM[3][2] = (2.0*e2*V*(e2*Mz_4+16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2-8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2))*ImB*sth2
                + 4.0*f2*gam2*m2*(4.0*a0_*(e2*Mz_2*cth
                - a0_2*f2*gam2*m2*V*(3.0+c2th)
                - gam2*m2*(4.0*e2*cth+f2*v0_2*(8.0*cth+V*(3.0+c2th))))*ReX
                + 2.0*V*v0_*(4.0*a0_2*f2*gam2*m2-e2*Mz_2
                + 4.0*gam2*m2*(e2+f2*v0_2))*ImY*sth2
                + e2*(-(a0_*(4.0*gam2*m2-Mz_2)*v0_*(4.0*cth+V*(3.0+c2th))*ReA1)
                + 2.0*Gz_*Mz_*(-2.0*a0_*cth*ImX+v0_*(a0_*(V+2.0*cth+V*cth2)*ImA1
                + V*v0_*ReB*sth2-V*ReY*sth2)))))/(4.0*denom);

    RDM[3][3] = (e2*(e2*Mz_4+16.0*gam4*m4*(e2+f2*v0_2)
                + Mz_2*(e2*(Gz_2-8.0*gam2*m2)
                - 4.0*f2*gam2*m2*v0_2)
                + 4.0*a0_2*f2*gam2*m2*(4.0*gam2*m2-Mz_2)*V*cth)*ReA1
                + 4.0*f2*gam2*m2*(-(e2*Gz_*Mz_*((v0_2+a0_2*V*cth)*ImA1-v0_*ImX))
                + v0_*(-(e2*Mz_2)+4.0*gam2*m2*(e2+f2*v0_2)
                + 4.0*a0_2*f2*gam2*m2*(1.0+2.0*V*cth))*ReX))/denom;



    // Combine RSM and RDM to form R
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            R[i][j] = RSM[i][j] + RDM[i][j];
        }
    }

    return R;
}
