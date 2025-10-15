#ifndef KinLib_H
#define KinLib_H

#include <TVector3.h>
#include <TLorentzVector.h>
#include <array>
#include <cmath>
#include <stdexcept>

namespace KinLib {


/**
 * @brief Replace Fortran bostdq: boost along arbitrary timelike 4-vector
 * @param idir  1 for lab→rest, -1 for rest→lab
 * @param vv    Timelike 4-vector {px,py,pz,E}
 * @param pp    Input 4-vector
 * @param q     Output boosted 4-vector
 */
inline void BostDQ(int idir,
                   const std::array<double,4>& vv,
                   const std::array<double,4>& pp,
                   std::array<double,4>& q)
{
    std::array<double,4> v = vv;
    std::array<double,4> p = pp;
    double wsp = 0.0;

    double mass2 = v[0]*v[0] - v[1]*v[1] - v[2]*v[2] - v[3]*v[3];
    if (mass2 <= 0.0) {
        throw std::invalid_argument("BostDQ: vv is not timelike");
    }

    double mass = std::sqrt(std::abs(mass2));
    if(idir == -1) {
        // boost in the direction of vv
        q[0] = (p[1]*v[1] + p[2]*v[2] + p[3]*v[3] + p[0]*v[0]) / mass;
        wsp = (q[0] + p[0]) / (v[0] + mass);
    } else if(idir == 1) {
        // boost in the opposite direction to vv
        q[0] = (-p[1]*v[1] - p[2]*v[2] - p[3]*v[3] + p[0]*v[0]) / mass;
        wsp = -(q[0] + p[0]) / (v[0] + mass);
    } else {
        throw std::invalid_argument("BostDQ: idir must be -1 or +1");
    }

    // Spatial components
    q[1] = p[1] + wsp * v[1];
    q[2] = p[2] + wsp * v[2];
    q[3] = p[3] + wsp * v[3];

    // TLorentzVector QQ(vv[1], vv[2], vv[3], vv[0]);
    // TVector3 beta = QQ.BoostVector();
    // TLorentzVector PP(pp[1], pp[2], pp[3], pp[0]);
    // idir==-1 ? PP.Boost(beta) : PP.Boost(-beta);
    // q = {{ PP.E(), PP.X(), PP.Y(), PP.Z() }};

}

/**
 * @brief Rotate 4-vector about z-axis (ROTOD3)
 */
inline void ROTOD3(double phi,
                   const std::array<double,4>& pvec,
                   std::array<double,4>& qvec)
{
    double cs =  std::cos(phi);
    double sn =  std::sin(phi);
    std::array<double,4> pvec_copy = pvec; // Copy to avoid modifying original
    qvec[1] = cs * pvec_copy[1] - sn * pvec_copy[2];
    qvec[2] = sn * pvec_copy[1] + cs * pvec_copy[2];
    qvec[3] = pvec_copy[3];
    qvec[0] = pvec_copy[0];
    //  TVector3 sp(pvec[1],pvec[2],pvec[3]);
    //  sp.RotateZ(phi);
    //  qvec = {{ pvec[0], sp.X(), sp.Y(), sp.Z() }};
}

/**
 * @brief Rotate 4-vector about y-axis (ROTOD2)
 */
inline void ROTOD2(double phi,
                   const std::array<double,4>& pvec,
                   std::array<double,4>& qvec)
{
    double cs =  std::cos(phi);
    double sn =  std::sin(phi);
    std::array<double,4> pvec_copy = pvec; // Copy to avoid modifying original
    qvec[1] = cs * pvec_copy[1] + sn * pvec_copy[3];
    qvec[2] = pvec_copy[2];
    qvec[3] = -sn * pvec_copy[1] + cs * pvec_copy[3];
    qvec[0] = pvec_copy[0];

    // TVector3 sp(pvec[1],pvec[2],pvec[3]);
    // sp.RotateY(phi);
    // qvec = {{ pvec[0], sp.X(), sp.Y(), sp.Z() }};
}


/**
 * @brief Compute the angle in [0, 2π) from Cartesian (x,y).
 */
inline double ANGFI(double x, double y) {
    constexpr double PI = 3.141592653589793238462643;
    double theta;

    if (std::abs(y) < std::abs(x)) {
        // small |y|: use arctan
        theta = std::atan(std::abs(y/x));
        if (x <= 0.0)
            theta = PI - theta;
    } else {
        // large |y|: use arccos for better numerical stability
        theta = std::acos(x / std::sqrt(x*x + y*y));
    }

    // reflect into [0,2π)
    if (y < 0.0)
        theta = 2.0 * PI - theta;

    return theta;
}

/**
* It computes the angle θ (theta) between the positive x-axis and the point (x, y) in the XY-plane,
* but with a nonstandard quadrant handling (unlike atan2).
*/
inline double ANGXY(double x, double y) {
    const double PI = TMath::Pi();
    double the = 0.0;

    if (std::abs(y) < std::abs(x)) {
        the = std::atan(std::abs(y / x));
        if (x <= 0.0) the = PI - the;
    } else {
        the = std::acos(x / std::sqrt(x * x + y * y));
    }
    return the;
}


} // namespace KinLib

#endif // KinLib_H


