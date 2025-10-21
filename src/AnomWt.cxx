#include "AnomWt.h"
#include "PhysicsConstants.h"


void AnomWt::compute() {
  // 1) Build the tau‐pair rest frame
  std::array<double,4> QQ, PB1{0,0,0,1}, PB2{0,0,0,1}, PBB ;
  for(int k=0; k<4; ++k) {
      QQ[k] = P1[k] + P2[k];
      // beam directions along ±z in tau‐pair frame
      PB1[k] = PB2[k] = 0.0;
  }

  PB1[3] = beamEnergy/2.0; PB2[3] = -beamEnergy/2.0;
  PB1[0] = PB2[0] = beamEnergy/2.0;
  double initial_theta = std:: acos(P1[3]/(std::sqrt(P1[1]*P1[1] + P2[2]*P2[2] + P2[3]*P2[3])));
  m_tau_momentum_z = P1[3];
  m_transverseMomentum = sqrt(pow((P1[1] + P2[1]), 2.) + pow((P1[2] + P2[2]), 2.));

//   double p1_z = P1[3];
//   double p2_z = P2[3];

  for(auto& arr : { &H1, &H2, &P1, &P2, &PB1, &PB2}) {
      KinLib::BostDQ(1, QQ, *arr, *arr);
  }
  //Invariant mass of tau pair system (s^\prime = (pbeam1 + pbeam2 - k)**2)
  E = sqrt(pow((P1[0] + P2[0]), 2.) - pow((P1[1] + P2[1]), 2.) - pow((P1[2] + P2[2]), 2.) - pow((P1[3] + P2[3]), 2.)) / 2.0;   //Since we are adding the photon

  // 2) Rotate about z to eliminate y‐components of P2
  double fi;
  double randSave = gRandom->Uniform();
  fi = KinLib::ANGFI(P1[1], P1[2]);
  for(auto& arr : { &H1, &H2, &P1, &P2, &PB1, &PB2}) {
      KinLib::ROTOD3(-fi, *arr, *arr);
  }

  // 3) Rotate about y to align taus along z
  //double thet = std::acos(-PB1[3]/std::sqrt(PB1[1]*PB1[1] + PB1[2]*PB1[2] + PB1[3]*PB1[3]));
  double thet;
  thet = KinLib::ANGXY(P1[3], P1[1]);
  for(auto& arr : { &H1, &H2, &P1, &P2, &PB1, &PB2}) {
      KinLib::ROTOD2(-thet, *arr, *arr);
  }

  double thetataupair = std::acos(PB1[3]/(std::sqrt(PB1[1]*PB1[1] + PB1[2]*PB1[2] +PB1[3]*PB1[3])));

  // 4) Align beam‐difference to x–z plane
  switch(frameOption){
      case 1:
          for(int k=0; k<4; ++k) PBB[k] = PB1[k] - PB2[k];//Collins-Soper frame
          m_phi = std::atan2(PB1[2], PB1[1]);
          if(m_phi<0) m_phi += 2.0*TMath::Pi();   // to have phi in [0, 2pi)
          break;
      case 2:
          for(int k=0; k<4; ++k) PBB[k] = PB1[k];
          m_phi = std::atan2(PB1[2], PB1[1]);
          if(m_phi<0) m_phi += 2.0*TMath::Pi();   // to have phi in [0, 2pi)
          break;
      case 3:
          for(int k=0; k<4; ++k) PBB[k] = -PB2[k];
          m_phi = std::atan2(-PB2[2], -PB2[1]);
          if(m_phi<0) m_phi += 2.0*TMath::Pi();   // to have phi in [0, 2pi)
          break;
  }



  //double fi1 = KinLib::ANGFI(PBB[1], PBB[2]);
  double fi1 = KinLib::ANGXY(PBB[1], PBB[2]);     //Modified by AT, 21.09.2025
  for(auto& arr : { &H1, &H2, &P1, &P2, &PB1, &PB2}) {
              KinLib::ROTOD3(-fi1, *arr, *arr);
  }



  TVector3 p1(P1[1],P1[2],P1[3]);
  TVector3 p2(P2[1],P2[2],P2[3]);
  TVector3 pb1(PB1[1],PB1[2],PB1[3]);
  TVector3 pb2(PB2[1],PB2[2],PB2[3]);
  double theta1 =  acos(PB1[3] /std::sqrt(PB1[1] * PB1[1] + PB1[2] * PB1[2] + PB1[3] * PB1[3]));
  double theta2 =  acos(-PB2[3] /std::sqrt(PB2[1] * PB2[1] + PB2[2] * PB2[2] + PB2[3] * PB2[3]));
  // Calculation of Asymmetry factor fromm theory
  DipoleQQRijRadCor probablity_calc(0);
  auto Rtt1 = probablity_calc.calculate(E, acos(0.5), 0., 0., 0., 0., 0., 0., 0., 0., 1.0);
  auto Rtt2 = probablity_calc.calculate(E, acos(-0.5), 0., 0., 0., 0., 0., 0., 0., 0., 1.0);
  A =  (5.0/2.0)*(Rtt1[3][3]-Rtt2[3][3])/(Rtt1[3][3]+Rtt2[3][3]);
  // New probabilities to Mustraal is added by Z. Was ---> inroduction of mass term effect and FBasymmetry.
  double T1 = PB1[0] * PB1[0] * (1.0 + cos(theta1) * cos(theta1)+ pow(Physics::m_tau,2)/(P1[0]*P1[0])*sin(theta1)*sin(theta1)+ A*cos(theta1));
  double T2 = PB2[0] * PB2[0] * (1.0 + cos(theta2) * cos(theta2)+ pow(Physics::m_tau,2)/(P1[0]*P1[0])*sin(theta2)*sin(theta2)+ A*cos(theta2));
  switch(frameOption){
      case 1:
          m_prob = 1.0;
          break;
      case 2:
          m_prob = T1 / (T1 + T2);
          break;
      case 3:
          m_prob = T2 / (T1 + T2);
          break;
  }


  //Finding Hard Photons
  std::array<double,4> TR;
  for(int k=0; k<4; ++k) TR[k] = P1[k]+P2[k];
  m_hardSoft =(beamEnergy * beamEnergy - TR[0]*TR[0])/(beamEnergy*beamEnergy);

  // 8) Final boosts of polarimeters and pions to tau rest frame
  for(auto& arr : { &H1}) {
      // Boost to tau- rest frame
      KinLib::BostDQ(1, P1, *arr, *arr);
  }
  for(auto& arr : { &H2}) {
      // Boost to tau+ rest frame
      KinLib::BostDQ(1, P2, *arr, *arr);
  }


  TVector3 Fp1(P1[1],P1[2],P1[3]);
  TVector3 Fp2(P2[1],P2[2],P2[3]);
  TVector3 Fpb1(PB1[1],PB1[2],PB1[3]);
  TVector3 Fpb2(PB2[1],PB2[2],PB2[3]);

  // 5) Compute angle with respect to tau-/z axis & prepare R0,R matrices
  switch(frameOption){
      case 1:
          m_theta = std::acos(PB1[3] /std::sqrt(PB1[1]*PB1[1] + PB1[2]*PB1[2] + PB1[3]*PB1[3]));
          break;
      case 2:
          m_theta = std::acos(PB1[3]/std::sqrt(PB1[1]*PB1[1] + PB1[2]*PB1[2] + PB1[3]*PB1[3]));
          break;
      case 3:
          m_theta = std::acos(-PB2[3]/std::sqrt(PB2[1]*PB2[1] + PB2[2]*PB2[2] + PB2[3]*PB2[3]));
          break;
  }

//   if(cos(m_theta)>0.98) {
//   std::cout<<"P1Z="<< p1_z<<std::endl;
//   std::cout<<"P2Z="<< p2_z<<std::endl;
//   std::cout<<"Initial theta="<< cos(initial_theta) << " Theta tau pair="<< cos(thetataupair) << " Final theta="<< cos(m_theta) << std::endl;
//   }

  m_tauMomentum = abs(sqrt(P1[1]*P1[1] + P1[2]*P1[2] + P1[3]*P1[3]));
  m_tauEnergy = P1[0];

  double ReA, ImA, ReB, ImB, ReX, ImX, ReY, ImY;

  ReA = 0.0; ImA = 0.0;
  ReB = 0.0; ImB = 0.0;
  ReX = 0.0; ImX = 0.0;
  ReY = 0.0; ImY = 0.0;

  // baseline (no loops)
  DipoleQQRijRadCor dipole_calculator1(0);
  auto R0 = dipole_calculator1.calculate(E, m_theta, ReA, ImA, ReB, ImB, ReX, ImX, ReY, ImY, 1.0);

  // with anomalous couplings & loops
  DipoleQQRijRadCor dipole_calculator2(1);
  auto R = dipole_calculator2.calculate(E, m_theta, Ar0, Ai0, Br0, Bi0, Ar0, Ai0, Br0, Bi0, 1.0);

  // 10) Extract weights
  m_wtME    = R[3][3] / R0[3][3];

  Rtt_val = R0[3][3];


  // Change the format of polarimetric vectors to match the expected format
  std::array<double,4> H1_copy = H1;
  std::array<double,4> H2_copy = H2;
  for(int k=0; k<3; ++k) {
      H1[k] = H1_copy[k+1];
      H2[k] = H2_copy[k+1];
  }
  H1[3] = H1_copy[0]; // Time component
  H2[3] = H2_copy[0]; // Time component


  // spin‐weight numerator/denominator
  double sign[4] = {1., -1., 1., -1.};
  m_wtSPIN0 = m_wtSPIN = 0.0;
  for(int i=0; i<4; ++i) {
      for(int j=0; j<4; ++j) {
          m_wtSPIN  += R[i][j]/R[3][3]  * H1[i] * H2[j] * sign[i] * sign[j];
          m_wtSPIN0 += R0[i][j]/R0[3][3]* H1[i] * H2[j] * sign[i] * sign[j];
      }
  }

  //TVector3 pb1(PB1[1], PB1[2], PB1[3]);               // e- direction
  TVector3 pbb;                                       // plane direction
  if (frameOption==1)      pbb = TVector3(PB2[1]-PB1[1], PB2[2]-PB1[2], PB2[3]-PB1[3]); // Collins–Soper
  else if (frameOption==2) pbb = TVector3(PB1[1], PB1[2], PB1[3]);                      // beam1-based
  else                     pbb = TVector3(PB2[1], PB2[2], PB2[3]);                      // beam2-based

//   TRotation Rbeam = BuildBeamRotation(pb1, pbb);
//   TMatrixD  S4    = MakeS4(Rbeam);
//   // Rotate R0 and R:  R' = S R S^T
//   TMatrixD R0m = ToMat4(R0);
//   TMatrixD Rm  = ToMat4(R);
//   TMatrixD R0m_beam = RotateR_S4(S4, R0m);
//   TMatrixD Rm_beam  = RotateR_S4(S4, Rm);
//   // Rotate H1, H2 (remember H = (Hx,Hy,Hz,Ht) at this point)
//   RotateH_TR(Rbeam, H1);
//   RotateH_TR(Rbeam, H2);


  // Helicity approximated weight calulation
  m_spinWeight_approx= 0.0;
  for (int m=0; m<2; ++m) {
      for (int n=0; n<2; ++n) {
          double sign_m = (m == 0) ? -1 : 1;
          double sign_n = (n == 0) ? -1 : 1;
          //double factor = m_wtSPIN0 / kkmc_wt; // to normalize to KKMC weight
          double term1 = H1[3] - sign_n * H1[2];
          double term2 = (R0[3][3] - sign_n * R0[3][2] + sign_m * R0[2][3] - sign_m * sign_n * R0[2][2])/R0[3][3];
          double term3 = H2[3] + sign_m * H2[2];
          wtw[m][n] = term1 * term2 * term3;
          //wtw[m][n] = term2;
          m_spinWeight_approx += wtw[m][n];
          //std::cout<<"term2="<<term2<<std::endl;
      }
  }

  H1_vec.SetXYZ(H1[0], H1[1], H1[2]);
  H2_vec.SetXYZ(H2[0], H2[1], H2[2]);


  m_spinWeight_approx /= 4.0;  // average over 4 combinations

  // normalize to make probabilities
  double wts = wtw[0][0] + wtw[0][1] + wtw[1][0] + wtw[1][1];
  wtwProb[0][0]=wtw[0][0]/wts; wtwProb[0][1]=wtw[0][1]/wts; wtwProb[1][0]=wtw[1][0]/wts; wtwProb[1][1]=wtw[1][1]/wts;
  wtw[0][0]=wtw[0][0]/wts;
  wtw[0][1]=wtw[0][0]+wtw[0][1]/wts;
  wtw[1][0]=wtw[0][1]+wtw[1][0]/wts;
  wtw[1][1]=wtw[1][0]+wtw[1][1]/wts;

  double rndm = gRandom->Uniform();

  m_hel1=0, m_hel2=0;
  if(rndm<wtw[0][0]){m_hel1=-1; m_hel2=-1;}   //attribute approx helicities. Signs need to be fixed.
  else if(rndm<wtw[0][1]){m_hel1=-1; m_hel2= 1;}
  else if(rndm<wtw[1][0]){m_hel1= 1; m_hel2=-1;}
  else if(rndm<wtw[1][1]){m_hel1= 1; m_hel2= 1;}



}
