#include "Observables.h"
#include "KinLib.h"
#include <TVector3.h>
#include <TMath.h>
#include <array>

void Observables::Acoplanarity(EventInitilizers& evtIn) {
  // Move to di-rho rest frame
  std::array<double,4> DIRHO;
  for(int k=0; k<4; ++k) {
      DIRHO[k] = evtIn.PIPM[k] + evtIn.PIZM[k] + evtIn.PIPL[k] + evtIn.PIZL[k];
  }
  for(auto& arr : { &evtIn.PIPM, &evtIn.PIZM, &evtIn.PIPL, &evtIn.PIZL}) {
      KinLib::BostDQ(1, DIRHO, *arr, *arr);
  }

  y1=0.0, y2=0.0;
  y1=(evtIn.PIPM[0] - evtIn.PIZM[0])/(evtIn.PIPM[0] + evtIn.PIZM[0]);
  y2=(evtIn.PIPL[0] - evtIn.PIZL[0])/(evtIn.PIPL[0] + evtIn.PIZL[0]);

  TVector3 pipm(evtIn.PIPM[1], evtIn.PIPM[2], evtIn.PIPM[3]);
  TVector3 pipl(evtIn.PIPL[1], evtIn.PIPL[2], evtIn.PIPL[3]);
  TVector3 pizm(evtIn.PIZM[1], evtIn.PIZM[2], evtIn.PIZM[3]);
  TVector3 pizl(evtIn.PIZL[1], evtIn.PIZL[2], evtIn.PIZL[3]);

  TVector3 n1 = pipm.Cross(pizm).Unit();
  TVector3 n2 = pipl.Cross(pizl).Unit();

  acoplanarity_angle = n1.Angle(n2);//Here the main issue is it will go from [0, pi].
  // Define a reference axis (e.g., the direction of one of the rhos, or their sum)
  TVector3 ref = (pipm + pizm).Unit();
  double sign_ang = ref.Dot(n2.Cross(n1)); // Scalar triple product
  if (sign_ang < 0) acoplanarity_angle = 2 * TMath::Pi() - acoplanarity_angle; // Now acc_angle in [0, 2pi]

  //Checking direction of polarimetric vectors
  // std::cout<<"Direction of H1: "<<evtIn.m_HvCloneTaum[1]<<" "<<evtIn.m_HvCloneTaum[2]<<" "<<evtIn.m_HvCloneTaum[3]<<std::endl;
  // std::cout<<"Direction of pipm+pizm: "<<(pipm+pizm).Unit()[0]<<" "<<(pipm+pizm).Unit()[1]<<" "<<(pipm+pizm).Unit()[2]<<std::endl;
  // std::cout<<"Direction of H2: "<<evtIn.m_HvCloneTaul[1]<<" "<<evtIn.m_HvCloneTaul[2]<<" "<<evtIn.m_HvCloneTaul[3]<<std::endl;
  // std::cout<<"Direction of pipl+pizl: "<<(pipl+pizl).Unit()[0]<<" "<<(pipl+pizl).Unit()[1]<<" "<<(pipl+pizl).Unit()[2]<<std::endl;



}
