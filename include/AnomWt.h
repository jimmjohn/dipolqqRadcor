#ifndef ANOMWT_CLASS_H
#define ANOMWT_CLASS_H

#include <array>
#include <cmath>
#include <iostream>
#include "KinLib.h"
#include "DipoleEERij.h"
#include "DipoleQQRijRadCor.h"
#include "EventInitilizers.h"
#include <TVector3.h>
#include <TRandom3.h>
#include <TRotation.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>


class AnomWt {
public:
    // ---- Configuration / inputs ----
    void setCouplings(double Ar, double Ai, double Br, double Bi) {
        Ar0 = Ar; Ai0 = Ai; Br0 = Br; Bi0 = Bi;
    }

    void setKinematics(EventInitilizers& evtIn) {
        H1 = evtIn.m_HvCloneTaum;
        H2 = evtIn.m_HvCloneTaul;
        P1 = evtIn.P1;
        P2 = evtIn.P2;
        beamEnergy = evtIn.beamEnergy;
        kkmc_wt = evtIn.wt;
    }

    // frameOption: 1 = Collins–Soper, 2 = use PB1, 3 = use PB2.
    void setOptions(int frameOptionI, int nevtsI = 0, int printEvtsI = 0) {
        frameOption = frameOptionI;
        nevts = nevtsI;
        printEvts = printEvtsI;
    }

    // ---- Do the computation ----
    void compute();

    // ---- Outputs (getters) ----
    double wtME()        const { return m_wtME; }
    double wtSPIN()      const { return m_wtSPIN; }
    double wtSPIN0()     const { return m_wtSPIN0; }
    double theta()       const { return m_theta; }
    double phi()         const { return m_phi; }
    double Invariant()         const { return E; }
    double hardSoft()    const { return m_hardSoft; }
    double prob()        const { return m_prob; }
    double spinApproxWt()  const { return m_spinWeight_approx; }
    int    approxHel1()        const { return m_hel1; }
    int    approxHel2()        const { return m_hel2; }
    bool   selected()    const { return m_selected; }
    double Rtt()         const { return Rtt_val; }
    double AFB()         const { return A;}
    double tauMomentum() const { return m_tauMomentum; }
    double tauEnergy()   const { return m_tauEnergy; }
    double directThetaVal() const { return directTheta; }
    double TransverseMomentum() const { return m_transverseMomentum; }
    double taumomentaz() const { return m_tau_momentum_z; }  // added
    double (&getwtw())[2][2] { return wtwProb;  }
    TVector3 getH1vec() const { return H1_vec; }
    TVector3 getH2vec() const { return H2_vec; }

    // If you prefer public fields, you can expose these directly instead of getters.

private:
    // Inputs
    double Ar0{0.0}, Ai0{0.0}, Br0{0.0}, Bi0{0.0};
    std::array<double,4> H1{0,0,0,0}, H2{0,0,0,0};
    std::array<double,4> P1{0,0,0,0}, P2{0,0,0,0};
    double beamEnergy{0.0};
    double kkmc_wt{0.0};
    int frameOption{1};
    int nevts{0}, printEvts{0};

    // Outputs
    double m_wtME{1.0}, m_wtSPIN{1.0}, m_wtSPIN0{1.0};
    double m_theta{0.0}, m_phi{0.0}, m_hardSoft{0.0}, m_prob{1.0};
    double m_spinWeight_approx{0.0};
    int m_hel1{0}, m_hel2{0};
    bool m_selected{true};
    double E{0.0};
    double Rtt_val{0.0};
    double A{0.0};
    double m_tauMomentum{0.0};
    double m_tauEnergy{0.0};
    double directTheta{0.0}; // for debugging
    double m_transverseMomentum{0.0};
    double m_tau_momentum_z{0.0};  // added
    double wtw[2][2] = {0.0};
    double wtwProb[2][2] = {0.0};
    TVector3 H1_vec;
    TVector3 H2_vec;

    // Build a right-handed beam basis: z_b = PB1, x_b = (PBB ⟂ z_b), y_b = z_b × x_b
    static inline TRotation BuildBeamRotation(const TVector3& pb1, const TVector3& pbb)
    {
        TVector3 z = pb1.Unit();

        // x from plane, orthogonalize to z
        TVector3 x = pbb - (pbb.Dot(z))*z;
        if (x.Mag2() < 1e-30) { // fallback if degenerate plane
            TVector3 a(1,0,0); if (std::fabs(z.X())>0.9) a.SetXYZ(0,1,0);
            x = a - (a.Dot(z))*z;
        }
        x = x.Unit();

        TVector3 y = z.Cross(x).Unit();   // right-handed

        TRotation R;                 // rows or cols? We'll verify below.
        R.SetXAxis(x);               // define the target axes
        R.SetYAxis(y);
        R.SetZAxis(z);

        // Ensure it maps pb1 to z (choose R or R^-1 accordingly)
        //TVector3 test = pb1; test.Transform(R);
        //if (std::fabs(test.X()) + std::fabs(test.Y()) > 1e-10) R.Invert();
        return R;
    }

    // Convert a 3x3 TRotation into a 4x4 block matrix S = diag(1, R3)
    static inline TMatrixD MakeS4(const TRotation& R)
    {
        TMatrixD S4(4,4); S4.Zero();
        S4(3,3) = 1.0; // time untouched
        for (int i=0;i<3;++i)
        for (int j=0;j<3;++j)
            S4(i, j) = R(i,j);  // (x,y,z) -> rows 1..3
        return S4;
    }


    // Rotate H with TRotation (no TMatrixD needed)
    static inline void RotateH_TR(const TRotation& R, std::array<double,4>& H) {
        TLorentzVector v(H[0], H[1], H[2], H[3]);  // (x,y,z,t)
        v.Transform(R);                            // rotates spatial part only
        H[0]=v.X(); H[1]=v.Y(); H[2]=v.Z(); H[3]=v.T();
    }

    // R' = S R S^T   (4x4 matrices)
    static inline TMatrixD RotateR_S4(const TMatrixD& S4, const TMatrixD& R4)
    {
        TMatrixD T = S4 * R4;
        TMatrixD ST = TMatrixD(TMatrixD::kTransposed, S4);
        return T * ST;
    }

    // Convenience: copy your R0 (any [4][4]-like) into TMatrixD
    template <class RArr>
    static inline TMatrixD ToMat4(const RArr& R)
    {
        TMatrixD M(4,4);
        for (int i=0;i<4;++i) for (int j=0;j<4;++j) M(i,j) = R[i][j];
        return M;
    }

    // And back, if you need the rotated numbers in your original container type
    template <class RArr>
    static inline void FromMat4(const TMatrixD& M, RArr& R)
    {
        for (int i=0;i<4;++i) for (int j=0;j<4;++j) R[i][j] = M(i,j);
    }


};

#endif // ANOMWT_CLASS_H
