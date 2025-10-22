#include <iostream>
#include <fstream>
#include <sstream>

#include <TLatex.h>
#include <TVector3.h>
#include <TApplication.h>
#include <TMath.h>
#include <TCutG.h>

#include "AnomWt.h"
#include "BornV.h"
#include "HistManager.h"
#include "EventInitilizers.h"
#include "HepMCReader.h"
#include "Observables.h"
#include "PhysicsConstants.h"


#include "HepMC3/ReaderAscii.h"
#include "HepMC3/GenEvent.h"

#include "DipoleEERij.h"
#include "DipoleQQRijRadCor.h"
#include "KKdizet.h"



int main(int argc, char** argv) {
    TApplication theApp("App", &argc, argv);  // removes ROOT args like -b, etc.
    // std::string filename;
    // int maxEvents = 1000000; // default Modify to parse from command line if needed

    const char* filename = theApp.Argv()[1];//argv[1] will not work as TApplication eats ROOT options (like -b, -q, etc.), and rewrites argc/argv
    std::ifstream infile(filename);
    std::cout << "============================================================\n";
    std::cout << "   Demo program associated with preprint:\n";
    std::cout << "   'arXiv:2509.04400 [hep-ph]'\n";
    std::cout << "   Authors: Jim John (jim.john@ifj.edu.pl)\n";
    std::cout << "            Ananya Tapadar (ananya.tapadar@ifj.edu.pl)\n";
    std::cout << "------------------------------------------------------------\n";
    std::cout << "   Purpose:\n";
    std::cout << "   Using polarimetric vectors and helicity-like quantities\n";
    std::cout << "   from HepMC3 files filled in by KKMCee.\n";
    std::cout << "------------------------------------------------------------\n";
    std::cout << "   Version: Sept 28, 2025\n";
    std::cout << "------------------------------------------------------------\n";
    std::cout << "   Notes:\n";
    std::cout << "   * Algorithm works for any tau decay modes.\n";
    std::cout << "   * Demo plots at present are prepared only for:\n";
    std::cout << "       tau -> pi+ pi0 nu\n";
    std::cout << "   * Other decay modes require distinct plots.\n";
    std::cout << "============================================================\n";
    if (!infile) {
        std::cout << "Error: cannot open '" << filename << "'.\n";
        return 1;
    }
    std::string line;

    EventInitilizers evtIn;
    HistManager* histSave = new HistManager();

    int nevts = 0; // Event counter
    int printEvts = 0;
    int countSelected = 0;
    double beamEnergy;
    double invariantCut = 3.0; // GeV//10.58 GeV
    //double invariantCut = 25.0; // GeV//91.18 GeV
    double sum_m2OverE2CS = 0.0;
    double sum_AFB= 0.0;

    //Try reading using HepMC3 Reader
    HepMC3::ReaderAscii reader(filename);
    if (reader.failed()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 1;
    }
    HepMC3::GenEvent evt;
    HepMCReader readEvent;

    // Load all the dizet tables
    std::cout << "Loading Dizet EW tables from directory: " << std::endl;
    BornV bornV;
    std::cout << "sizeof(BornV) = " << sizeof(BornV) << " bytes\n";
    std::string dizetDir = "../dizet"; // Adjust path as needed

    bornV.BornV_ReadAll(dizetDir);



   TCutG *cutgD = new TCutG("Diagonal_CUTG",6);
   cutgD->SetVarX("");
   cutgD->SetVarY("");
   cutgD->SetTitle("Diagonal_CUTG");
   cutgD->SetFillStyle(1000);
   cutgD->SetPoint(0,1.13405,1.25);
   cutgD->SetPoint(1,1.18301,1.08333);
   cutgD->SetPoint(2,1.55162,1.51042);
   cutgD->SetPoint(3,1.50554,1.63021);
   cutgD->SetPoint(4,1.13693,1.23437);
   cutgD->SetPoint(5,1.13405,1.25);

   TCutG *cutgH = new TCutG("Horizontal_CUTG",6);
   cutgH->SetVarX("");
   cutgH->SetVarY("");
   cutgH->SetTitle("Horizontal_CUTG");
   cutgH->SetFillStyle(1000);
   cutgH->SetPoint(0,1.20029,1.06771);
   cutgH->SetPoint(1,1.20605,0.916667);
   cutgH->SetPoint(2,1.74744,0.932292);
   cutgH->SetPoint(3,1.74168,1.08333);
   cutgH->SetPoint(4,1.20029,1.06771);
   cutgH->SetPoint(5,1.20029,1.06771);

    KKdizet::instance().ReadEWtabs();  // Read the EW tables from file


    while (reader.read_event(evt)) {
        if (reader.failed()) break;
        // Skip empty lines

        bool selected = true;;

        if(nevts > 10000000) break; // Limit to 1000000 events
        evtIn = {};
        if(nevts<printEvts){std::cout<< "\n"<< "Evt number: "<<nevts<<std::endl;}
        if(nevts%10000==0){std::cout<< "\n"<< "Evt number: "<<nevts<<std::endl;}
        nevts++;

        readEvent.ReadAndFillEvent(evtIn, evt);

        int frameOption = 1;
        double acc_angle, y1y2;

        if(evtIn.P1[0]==0 || evtIn.P2[0]==0) continue; //skip if tau+ or tau- not found
        if(evtIn.PIPL[0]==0 || evtIn.PIPM[0]==0) continue; //skip if pi+ or pi- not found
        if(evtIn.heTaum==0 || evtIn.heTaul==0) continue; //skip if helicities not found
        if(evtIn.pvTaum_found==false || evtIn.pvTaul_found==false) continue; //skip if polarimetric vectors not found

        histSave->tauZMomenta->Fill(evtIn.P1[3]);
        histSave->tauZMomenta->Fill(evtIn.P2[3]);

        //Observables calculation
        Observables obs;
        obs.Acoplanarity(evtIn);
        acc_angle = obs.acoplanarity_angle;
        y1y2 = obs.y1 * obs.y2;

        AnomWt anomwtCS1, anomwtCS2, anomwtCS3;
        anomwtCS1.setCouplings(0.04, 0.0, 0.0, 0.0);   //Left Plot
        anomwtCS2.setCouplings(0.0, 0.0, 0.04, 0.0);   //Middle Plot
        anomwtCS3.setCouplings(0.02828, 0.0, 0.02828, 0.0); // Right Plot
        anomwtCS1.setKinematics(evtIn);
        anomwtCS2.setKinematics(evtIn);
        anomwtCS3.setKinematics(evtIn);
        anomwtCS1.setOptions(frameOption, nevts, printEvts);
        anomwtCS2.setOptions(frameOption, nevts, printEvts);
        anomwtCS3.setOptions(frameOption, nevts, printEvts);
        anomwtCS1.compute();
        anomwtCS2.compute();
        anomwtCS3.compute();

        //if(anomwtCS1.hardSoft()<0.02) {
        if(1) {
        histSave->wt_corr_KT_KA->Fill(evtIn.wt, evtIn.wt_approx);
        histSave->wt_corr_KT_ST->Fill(evtIn.wt, anomwtCS1.wtSPIN0());
        histSave->wt_corr_KT_SA->Fill(evtIn.wt, anomwtCS1.spinApproxWt());
        histSave->wt_corr_KA_SA->Fill(evtIn.wt_approx, anomwtCS1.spinApproxWt());
        histSave->wt_corr_ST_SA->Fill(anomwtCS1.wtSPIN0(), anomwtCS1.spinApproxWt());
        }

        histSave->hardSoft_Histo->Fill(anomwtCS1.hardSoft());

        if(y1y2 > 0.0 && anomwtCS1.hardSoft()<0.02 ) {
            histSave->acoplanarity_angle1->Fill(acc_angle, anomwtCS1.wtSPIN0()/evtIn.wt);
            histSave->h_angle_01->Fill(acc_angle, anomwtCS1.wtME());
            histSave->h_angle1->Fill(acc_angle, anomwtCS1.wtSPIN()/anomwtCS1.wtSPIN0()*anomwtCS1.wtME());
        }
        if(y1y2 > 0.0 && anomwtCS2.hardSoft()<0.02) {
            histSave->acoplanarity_angle2->Fill(acc_angle, anomwtCS2.wtSPIN0()/evtIn.wt);
            histSave->h_angle_02->Fill(acc_angle, anomwtCS2.wtME());
            histSave->h_angle2->Fill(acc_angle, anomwtCS2.wtSPIN()/anomwtCS2.wtSPIN0()*anomwtCS2.wtME());
        }
        if(y1y2 > 0.0 && anomwtCS3.hardSoft()<0.02) {
            histSave->acoplanarity_angle3->Fill(acc_angle, anomwtCS3.wtSPIN0()/evtIn.wt);
            histSave->h_angle_03->Fill(acc_angle, anomwtCS3.wtME());
            histSave->h_angle3->Fill(acc_angle, anomwtCS3.wtSPIN()/anomwtCS3.wtSPIN0()*anomwtCS3.wtME());
        }

        double E2CS = pow(anomwtCS1.Invariant(),2);
        if (E2CS > 0 && anomwtCS1.Invariant() < invariantCut) {
            sum_m2OverE2CS += 1 / E2CS;
            sum_AFB += anomwtCS1.AFB();
            countSelected++;
        }

        //Mustral Frame
        frameOption=2;
        AnomWt anomwtPB1;
        anomwtPB1.setCouplings(0.0, 0.0, 0.0, 0.0);
        anomwtPB1.setKinematics(evtIn);
        anomwtPB1.setOptions(frameOption, nevts, printEvts);
        anomwtPB1.compute();

        frameOption=3;
        AnomWt anomwtPB2;
        anomwtPB2.setCouplings(0.0, 0.0, 0.0, 0.0);
        anomwtPB2.setKinematics(evtIn);
        anomwtPB2.setOptions(frameOption, nevts, printEvts);
        anomwtPB2.compute();


        double wtMustraal0 = anomwtPB1.prob()*anomwtPB1.wtSPIN0() + anomwtPB2.prob()*anomwtPB2.wtSPIN0();
        double wtMustraal = anomwtPB1.prob()*anomwtPB1.wtSPIN() + anomwtPB2.prob()*anomwtPB2.wtSPIN();

        double thetaMustraal, phiMustraal;
        double approxHel1MS, approxHel2MS;
        if(anomwtPB1.prob()>anomwtPB2.prob()) {
            thetaMustraal = anomwtPB1.theta();
            phiMustraal = anomwtPB1.phi();
            approxHel1MS = anomwtPB1.approxHel1();
            approxHel2MS = anomwtPB1.approxHel2();
        }
        else {
            thetaMustraal = anomwtPB2.theta();
            phiMustraal = anomwtPB2.phi();
            approxHel1MS = anomwtPB2.approxHel1();
            approxHel2MS = anomwtPB2.approxHel2();
        }

        if(cos(anomwtCS1.theta())<0.05 && cos(anomwtCS1.theta())>-0.05) {
        //if(cos(anomwtCS1.theta())>0.98 || cos(anomwtCS1.theta())<-0.98) {
        //if(cos(anomwtCS1.theta())>0.98) {
        //if(anomwtCS1.hardSoft()<=0.02) {
        //if(1) {
        //if(cutgH->IsInside(evtIn.wt, evtIn.wt_approx)) {
            //std::cout<<cos(anomwtCS1.theta())<<" "<<cos(thetaMustraal)<<std::endl;
            histSave->helicity_corr->Fill(evtIn.heTaum, evtIn.heTaul);
            histSave->helicity_corr_standalone->Fill(anomwtCS1.approxHel1(), anomwtCS1.approxHel2());
            histSave->wtwPP->Fill(anomwtCS1.getwtw()[0][0]);
            histSave->wtwMM->Fill(anomwtCS1.getwtw()[1][1]);
            histSave->wtwPM->Fill(anomwtCS1.getwtw()[0][1]);
            histSave->wtwMP->Fill(anomwtCS1.getwtw()[1][0]);
        }
        //if(cos(thetaMustraal)<0.05 && cos(thetaMustraal)>-0.05) {
        if(anomwtCS1.hardSoft()<=0.02) {
        //if(1) {
            histSave->helicity_corr_st_mustraal->Fill(approxHel1MS, approxHel2MS);
        }

        // histSave->H1->SetPoint(nevts, anomwtCS1.getH1vec().X(), anomwtCS1.getH1vec().Y(), anomwtCS1.getH1vec().Z());
        // histSave->H2->SetPoint(nevts, anomwtCS1.getH2vec().X(), anomwtCS1.getH2vec().Y(), anomwtCS1.getH2vec().Z());
        // histSave->H1_H2->SetPoint(nevts, anomwtCS1.getH1vec().X() - anomwtCS1.getH2vec().X(), anomwtCS1.getH1vec().Y() - anomwtCS1.getH2vec().Y(), anomwtCS1.getH1vec().Z() - anomwtCS1.getH2vec().Z());

        if(anomwtCS1.hardSoft()<=0.02) {
            histSave->collins_soper_corr_soft->Fill(evtIn.wt, anomwtCS1.wtSPIN0());
            histSave->p1_frame_corr_soft->Fill(evtIn.wt, anomwtPB1.wtSPIN0());
            histSave->p2_frame_corr_soft->Fill(evtIn.wt, anomwtPB2.wtSPIN0());
            histSave->mustraal_corr_soft->Fill(evtIn.wt, wtMustraal0);
        } else if(anomwtCS1.hardSoft()<0.8 && anomwtCS1.hardSoft()>0.02) {
            histSave->collins_soper_corr_hard->Fill(evtIn.wt, anomwtCS1.wtSPIN0());
            histSave->p1_frame_corr_hard->Fill(evtIn.wt, anomwtPB1.wtSPIN0());
            histSave->p2_frame_corr_hard->Fill(evtIn.wt, anomwtPB2.wtSPIN0());
            histSave->mustraal_corr_hard->Fill(evtIn.wt, wtMustraal0);
        }

        beamEnergy = evtIn.beamEnergy;

        if(anomwtCS1.Invariant()<invariantCut) {
            histSave->theta_dist_cs->Fill(cos(anomwtCS1.theta()));
            histSave->phi_dist_cs->Fill(anomwtCS1.phi());
        }
        if(anomwtPB1.Invariant()<invariantCut) {
            histSave->theta_dist_ms->Fill(cos(thetaMustraal));
            histSave->phi_dist_ms->Fill(phiMustraal);
        }
        // if(anomwtCS1.hardSoft()<=0.02) histSave->theta_dist_cs->Fill(cos(anomwtCS1.theta()));
        // if(anomwtCS1.hardSoft()<=0.02) histSave->theta_dist_ms->Fill(cos(thetaMustraal));

        histSave->Invariant_dist->Fill(anomwtCS1.Invariant());

        //std::cout<<"cosditrctTheta: "<<std::cos(anomwtPB1.directThetaVal())<<" mustralTheta: "<<std::cos(thetaMustraal)<<std::endl;
        histSave->directTheta_dist->Fill(std::cos(anomwtPB1.directThetaVal()));
        //if(anomwtCS1.Invariant()>beamEnergy/2.0-0



    } //Event loop Ends here


    double avg_m2OverE2 = 0.0;
    double avg_AFB = 0.0;
    if (countSelected > 0) {
        avg_m2OverE2 = pow(Physics::m_tau,2)* sum_m2OverE2CS / countSelected;
        avg_AFB = sum_AFB / countSelected;
    }
    std::cout << "<m^2/E^2> = " << avg_m2OverE2 << std::endl;
    std::cout << "<AFB> = " << avg_AFB << std::endl;

    Int_t ba=0;
    for (double costheta = -1; costheta <= 1; costheta += 0.001) {
        ba++;
        double energy = beamEnergy/2.0;
        double theta = std::acos(costheta);
        DipoleQQRijRadCor dipole_calculator3(0);
        auto RMat = dipole_calculator3.calculate(energy, theta, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        double Rtt_theory = RMat[3][3];   // no arbitrary factor
        histSave->Rtt->SetPoint(ba, costheta, Rtt_theory);
    }

    // std::cout << std::scientific << std::uppercase << std::setprecision(15);
    // // Two scientific columns similar to Fortran
    // const int idxw = 2;   // width for indices
    // const int colw = 23;  // width for numbers (tweak if you want tighter/wider)


    // std::cout<<"Jim Compare the values for E = 5.0 and theta = 0.5" <<std::endl;
    // double thetaX = 3.14159265358979323846/3.0;
    // DipoleQQRijRadCor dipole_calculator4(1);
    // auto RMat = dipole_calculator4.calculate(5.0, thetaX, 0.04, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 1.0);

    // DipoleQQRijRadCor dipole_calculator5(1);
    // auto R0Mat = dipole_calculator5.calculate(5.0, thetaX, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

    // for (int i = 1; i <= 4; ++i) {
    //     for (int j = 1; j <= 4; ++j) {
    //         std::cout << "R[" << setw(idxw) << i << "][" << setw(idxw) << j << "] = "
    //                   << setw(colw) << RMat[i-1][j-1] << " " << setw(colw) << R0Mat[i-1][j-1] << '\n';
    //     }
    // }


    histSave->ratio_angle1 = (TH1D*)histSave->h_angle1->Clone("ratio_angle1");
    histSave->ratio_angle1->Divide(histSave->h_angle_01);
    histSave->ratio_angle2 = (TH1D*)histSave->h_angle2->Clone("ratio_angle2");
    histSave->ratio_angle2->Divide(histSave->h_angle_02);
    histSave->ratio_angle3 = (TH1D*)histSave->h_angle3->Clone("ratio_angle3");
    histSave->ratio_angle3->Divide(histSave->h_angle_03);

    histSave->collins_by_mustraal = (TH2F*)histSave->collins_soper_corr_hard->Clone("collins_by_mustraal");
    histSave->collins_by_mustraal->Divide(histSave->mustraal_corr_hard);


    // Draw all canvases
    histSave->DrawHistograms(beamEnergy, avg_m2OverE2, avg_AFB);


    //theApp.Run();  // Keeps the GUI open

    return 0;
}
