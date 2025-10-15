#include "HistManager.h"

#include <TMath.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TPostScript.h>
#include <TLatex.h>
#include <TLine.h>
#include <TString.h>
#include <TLegend.h>

HistManager::~HistManager() {
}

void HistManager::DrawHistograms(double beamEnergy, double avg_m2OverE2, double avg_AFB) {
    TDirectory* currentDir = gDirectory;

    // Create a new directory for the histograms
    TDirectory* histDir = currentDir->mkdir("Histograms");
    histDir->cd();
    //TPostScript ps("output_histos.ps",111);
    //ps.Range(20,30); //ps.Range(10,20);

        // Draw the histogram
    TCanvas* c1 = new TCanvas("c1", "Helicity Correlation", 800, 600);
    helicity_corr->SetStats(0);
    helicity_corr->GetXaxis()->SetTitle("#tau^{-}");
    helicity_corr->GetYaxis()->SetTitle("#tau^{+}");
    TString title1 = Form("Helicity Correlation at %.2f GeV from KKMCee", beamEnergy);
    //helicity_corr->SetTitle(title1);
    helicity_corr->SetTitle("");
    helicity_corr->GetXaxis()->CenterTitle();
    helicity_corr->GetYaxis()->CenterTitle();
    helicity_corr->SetMarkerStyle(20);
    helicity_corr->SetMarkerColor(kBlack);
    // hide ticks/labels but KEEP titles
    helicity_corr->GetXaxis()->SetTickLength(0);
    helicity_corr->GetYaxis()->SetTickLength(0);
    helicity_corr->GetXaxis()->SetLabelSize(0);   // instead of SetNdivisions(0,...)
    helicity_corr->GetYaxis()->SetLabelSize(0);
    helicity_corr->GetXaxis()->SetTitleSize(0.06);
    helicity_corr->GetYaxis()->SetTitleSize(0.06);
    helicity_corr->GetXaxis()->SetTitleOffset(0.7);

    // Get axis limits
    double xLeft  = helicity_corr->GetXaxis()->GetBinCenter(1);
    double xRight = helicity_corr->GetXaxis()->GetBinCenter(2);
    double yBottom = helicity_corr->GetYaxis()->GetBinCenter(1);
    double yTop = helicity_corr->GetYaxis()->GetBinCenter(2);
    // Get vertical position for the text (middle of y range)

    //c1->SetGrid();
    gPad->SetRightMargin(0.15);
    helicity_corr->SetMarkerSize(4.0);
    helicity_corr->Draw("COLZ TEXT");
    c1->Update();
    TLatex tex1;
    tex1.SetTextSize(0.06);
    tex1.SetTextFont(42);
    tex1.SetTextAlign(22);  // center alignment
    tex1.DrawLatex(xLeft, -1.2, "-1");
    tex1.DrawLatex(xRight, -1.2, " 1");
    tex1.DrawLatex(-1.2, yBottom, "-1");
    tex1.DrawLatex(-1.2, yTop, " 1");
    gPad->Modified();
    gPad->Update();


    TCanvas* c2 = new TCanvas("c2", "Helicity Correlation_Standalone", 800, 600);
    helicity_corr_standalone->SetStats(0);
    helicity_corr_standalone->SetTitle("");
    helicity_corr_standalone->GetXaxis()->SetTitle("#tau^{-}");
    helicity_corr_standalone->GetYaxis()->SetTitle("#tau^{+}");
    helicity_corr_standalone->GetXaxis()->CenterTitle();
    helicity_corr_standalone->GetYaxis()->CenterTitle();
    helicity_corr_standalone->SetMarkerStyle(20);
    helicity_corr_standalone->SetMarkerColor(kBlack);
    // hide ticks/labels but KEEP titles
    helicity_corr_standalone->GetXaxis()->SetTickLength(0);
    helicity_corr_standalone->GetYaxis()->SetTickLength(0);
    helicity_corr_standalone->GetXaxis()->SetLabelSize(0);   // instead of SetNdivisions(0,...)
    helicity_corr_standalone->GetYaxis()->SetLabelSize(0);
    helicity_corr_standalone->GetXaxis()->SetTitleSize(0.06);
    helicity_corr_standalone->GetYaxis()->SetTitleSize(0.06);
    helicity_corr_standalone->GetXaxis()->SetTitleOffset(0.7);
    // Get axis limits
    xLeft  = helicity_corr_standalone->GetXaxis()->GetBinCenter(1);
    xRight = helicity_corr_standalone->GetXaxis()->GetBinCenter(2);
    yBottom = helicity_corr_standalone->GetYaxis()->GetBinCenter(1);
    yTop = helicity_corr_standalone->GetYaxis()->GetBinCenter(2);
    gPad->SetRightMargin(0.15);

    helicity_corr_standalone->SetMarkerSize(4.0);
    helicity_corr_standalone->Draw("COLZ TEXT");
    c2->Update();
    TLatex tex2;
    tex2.SetTextSize(0.06);
    tex2.SetTextFont(42);
    tex2.SetTextAlign(22);  // center alignment
    tex2.DrawLatex(xLeft, -1.2, "-1");
    tex2.DrawLatex(xRight, -1.2, " 1");
    tex2.DrawLatex(-1.2, yBottom, "-1");
    tex2.DrawLatex(-1.2, yTop, " 1");
    gPad->Modified();
    gPad->Update();

    //JMJ added to save the standalone helicity correlation
    TCanvas* c2b = new TCanvas("c2b", "Helicity Correlation_Standalone_Mustraal", 800, 600);
    helicity_corr_st_mustraal->SetStats(0);
    helicity_corr_st_mustraal->GetXaxis()->CenterTitle();
    helicity_corr_st_mustraal->GetYaxis()->CenterTitle();
    helicity_corr_st_mustraal->SetMarkerStyle(20);
    helicity_corr_st_mustraal->SetMarkerColor(kBlack);
    c2b->SetGrid();
    helicity_corr_st_mustraal->Draw("COLZ TEXT");
    c2b->Update();

    //ps.NewPage();
    TCanvas* c3 = new TCanvas("c3", "WT Correlation", 800, 600);
    wt_corr_KT_KA->SetStats(0);
    wt_corr_KT_KA->GetXaxis()->SetTitle("SpinWT");
    wt_corr_KT_KA->GetYaxis()->SetTitle("SpinWThelApprox");
    TString title3 = Form("SpinWT Correlation at %.2f GeV from KKMCee", beamEnergy);
    wt_corr_KT_KA->SetTitle("");
    wt_corr_KT_KA->GetXaxis()->CenterTitle();
    wt_corr_KT_KA->GetYaxis()->CenterTitle();
    wt_corr_KT_KA->SetMarkerStyle(20);
    wt_corr_KT_KA->SetMarkerSize(0.5);
    wt_corr_KT_KA->SetMarkerColor(kBlue);
    c3->SetGrid();
    gPad->SetRightMargin(0.15);
    wt_corr_KT_KA->Draw("COLZ");
    TLatex tex;
    tex.SetTextSize(0.04);
    //tex.DrawLatexNDC(0.15, 0.8, Form("Correlation Factor: %.2f", wt_corr_KT_KA->GetCorrelationFactor()));
    c3->Update();

    //ps.NewPage();
    // Draw acoplanarity angle
    TCanvas* c4 = new TCanvas("c4", "Acoplanarity Angle", 800, 600);
    acoplanarity_angle1->SetStats(0);
    acoplanarity_angle1->GetXaxis()->SetTitle("Acoplanarity Angle (degrees)");
    acoplanarity_angle1->GetYaxis()->SetTitle("Spin Weight");
    TString title4 = Form("Acoplanarity Angle vs Spin Weight at %.2f GeV with spin weight", beamEnergy);
    acoplanarity_angle1->SetTitle(title4);
    acoplanarity_angle1->GetXaxis()->CenterTitle();
    acoplanarity_angle1->GetYaxis()->CenterTitle();
    acoplanarity_angle1->SetMarkerStyle(20);
    acoplanarity_angle1->SetMarkerSize(0.5);
    acoplanarity_angle1->SetMarkerColor(kRed);
    c4->SetGrid();
    acoplanarity_angle1->Draw("COLZ");
    c4->Update();

    //ps.NewPage();
    TCanvas* c5 = new TCanvas("c5", "Acoplanarity Angle 0 vs Spin Weight", 800, 600);
    // h_angle_0->SetStats(0);
    // h_angle_0->GetXaxis()->SetTitle("Acoplanarity Angle (degrees)");
    // h_angle_0->GetYaxis()->SetTitle("Events");
    // h_angle_0->SetTitle("Acoplanarity Angle vs Events at 10 GeV without spin weight");
    // h_angle_0->GetXaxis()->CenterTitle();
    // h_angle_0->GetYaxis()->CenterTitle();
    h_angle_01->SetMarkerStyle(20);
    h_angle_01->SetMarkerSize(0.5);
    h_angle_01->SetMarkerColor(kBlue);
    c4->SetGrid();
    h_angle_01->Draw("E1");
    c5->Update();

    //ps.NewPage();
    TCanvas* c6 = new TCanvas("c6", "Acoplanarity Angle vs Spin Weight", 800, 600);
    //h_angle->SetStats(0);
    //h_angle->GetXaxis()->SetTitle("Acoplanarity Angle (degrees)");
    //h_angle->GetYaxis()->SetTitle("Events");
    //h_angle->SetTitle("Acoplanarity Angle vs Events at 10 GeV with spin weight");
    h_angle1->GetXaxis()->CenterTitle();
    h_angle1->GetYaxis()->CenterTitle();
    h_angle1->SetMarkerStyle(20);
    h_angle1->SetMarkerSize(0.5);
    h_angle1->SetMarkerColor(kRed);
    c6->SetGrid();
    h_angle1->Draw("E1");
    c6->Update();

    //ps.NewPage();
    TCanvas* c7 = new TCanvas("c7", "Ratio of Acoplanarity Angles", 1200, 600);
    c7->Divide(3,1);
    c7->cd(1);
    ratio_angle1->SetStats(0);
    //ratio_angle->GetXaxis()->SetTitle("Acoplanarity Angle (degrees)");
    ratio_angle1->GetXaxis()->SetTitle("#phi (radians)");
    ratio_angle1->GetYaxis()->SetTitle("wt_{spin}");
    ratio_angle1->SetTitle("");
    ratio_angle1->GetXaxis()->CenterTitle();
    ratio_angle1->GetYaxis()->CenterTitle();
    ratio_angle1->GetXaxis()->SetLabelSize(0.05);
    ratio_angle1->GetYaxis()->SetLabelSize(0.05);
    ratio_angle1->GetXaxis()->SetTitleSize(0.05);
    ratio_angle1->GetYaxis()->SetTitleSize(0.05);
    ratio_angle1->GetYaxis()->SetRangeUser(0.99,1.01);
    ratio_angle1->SetMarkerStyle(20);
    ratio_angle1->SetMarkerSize(0.5);
    ratio_angle1->SetMarkerColor(kGreen);
    //ratio_angle1->GetXaxis()->SetLabelSize(0.0);
    //ratio_angle->GetXaxis()->SetNdivisions(6);
    gPad->SetLeftMargin(0.2);
    gPad->SetGridx(true);
    ratio_angle1->GetXaxis()->SetNdivisions(4, kFALSE);
    ratio_angle1->GetXaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, "0");
    ratio_angle1->GetXaxis()->ChangeLabel(2, -1, -1, -1, -1, -1, "#pi/2");
    ratio_angle1->GetXaxis()->ChangeLabel(3, -1, -1, -1, -1, -1, "#pi");
    ratio_angle1->GetXaxis()->ChangeLabel(4, -1, -1, -1, -1, -1, "3#pi/2");
    ratio_angle1->GetXaxis()->ChangeLabel(5, -1, -1, -1, -1, -1, "2#pi");
    ratio_angle1->Draw("HIST");
    // Horizontal guide at y=1 across the full x-range
    TLine* y1 = new TLine(0, 1.0, 2*TMath::Pi(), 1.0);
    y1->SetLineStyle(2);       // dashed
    y1->SetLineWidth(1);
    // y1->SetLineColor(kGray+2); // optional
    y1->Draw("same");
    gPad->Update(); // Update the pad to reflect the new axis

    c7->cd(2);
    ratio_angle2->SetStats(0);
    //ratio_angle->GetXaxis()->SetTitle("Acoplanarity Angle (degrees)");
    ratio_angle2->GetXaxis()->SetTitle("#phi (radians)");
    ratio_angle2->GetYaxis()->SetTitle("wt_{spin}");
    ratio_angle2->SetTitle("");
    ratio_angle2->GetXaxis()->CenterTitle();
    ratio_angle2->GetYaxis()->CenterTitle();
    ratio_angle2->GetXaxis()->SetLabelSize(0.05);
    ratio_angle2->GetYaxis()->SetLabelSize(0.05);
    ratio_angle2->GetXaxis()->SetTitleSize(0.05);
    ratio_angle2->GetYaxis()->SetTitleSize(0.05);
    ratio_angle2->GetYaxis()->SetRangeUser(0.99,1.01);
    ratio_angle2->SetMarkerStyle(20);
    ratio_angle2->SetMarkerSize(0.5);
    ratio_angle2->SetMarkerColor(kGreen);
    gPad->SetGridx(true);
    gPad->SetLeftMargin(0.2);
    ratio_angle2->GetXaxis()->SetNdivisions(4, kFALSE);
    ratio_angle2->GetXaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, "0");
    ratio_angle2->GetXaxis()->ChangeLabel(2, -1, -1, -1, -1, -1, "#pi/2");
    ratio_angle2->GetXaxis()->ChangeLabel(3, -1, -1, -1, -1, -1, "#pi");
    ratio_angle2->GetXaxis()->ChangeLabel(4, -1, -1, -1, -1, -1, "3#pi/2");
    ratio_angle2->GetXaxis()->ChangeLabel(5, -1, -1, -1, -1, -1, "2#pi");
    ratio_angle2->Draw("HIST");
    y1->Draw("same");
    gPad->Update(); // Update the pad to reflect the new axis


    c7->cd(3);
    ratio_angle3->SetStats(0);
    //ratio_angle->GetXaxis()->SetTitle("Acoplanarity Angle (degrees)");
    ratio_angle3->GetXaxis()->SetTitle("#phi (radians)");
    ratio_angle3->GetYaxis()->SetTitle("wt_{spin}");
    ratio_angle3->SetTitle("");
    ratio_angle3->GetXaxis()->CenterTitle();
    ratio_angle3->GetYaxis()->CenterTitle();
    ratio_angle3->GetXaxis()->SetLabelSize(0.05);
    ratio_angle3->GetYaxis()->SetLabelSize(0.05);
    ratio_angle3->GetXaxis()->SetTitleSize(0.05);
    ratio_angle3->GetYaxis()->SetTitleSize(0.05);
    ratio_angle3->GetYaxis()->SetRangeUser(0.99,1.01);
    ratio_angle3->SetMarkerStyle(20);
    ratio_angle3->SetMarkerSize(0.5);
    ratio_angle3->SetMarkerColor(kGreen);
    gPad->SetGridx(true);
    gPad->SetLeftMargin(0.2);
    ratio_angle3->GetXaxis()->SetNdivisions(4, kFALSE);
    ratio_angle3->GetXaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, "0");
    ratio_angle3->GetXaxis()->ChangeLabel(2, -1, -1, -1, -1, -1, "#pi/2");
    ratio_angle3->GetXaxis()->ChangeLabel(3, -1, -1, -1, -1, -1, "#pi");
    ratio_angle3->GetXaxis()->ChangeLabel(4, -1, -1, -1, -1, -1, "3#pi/2");
    ratio_angle3->GetXaxis()->ChangeLabel(5, -1, -1, -1, -1, -1, "2#pi");
    ratio_angle3->Draw("HIST");
    y1->Draw("same");
    gPad->Update(); // Update the pad to reflect the new axis
    c7->Update();

    //ps.NewPage();
    TCanvas* c8 = new TCanvas("c8","weight_correlations soft photon case",1200,1200);
    c8->Divide(2,2);
    c8->cd(1);
    collins_soper_corr_soft->Draw("colz");
    c8->cd(2);
    p1_frame_corr_soft->Draw("colz");
    c8->cd(3);
    p2_frame_corr_soft->Draw("colz");
    c8->cd(4);
    mustraal_corr_soft->Draw("colz");
    c8->Update();

    TCanvas* c9 = new TCanvas("c9","weight_correlations hard photon case",1200,1200);
    c9->Divide(2,2);
    c9->cd(1);
    collins_soper_corr_hard->Draw("colz");
    c9->cd(2);
    p1_frame_corr_hard->Draw("colz");
    c9->cd(3);
    p2_frame_corr_hard->Draw("colz");
    c9->cd(4);
    mustraal_corr_hard->Draw("colz");
    c9->Update();


    //ps.NewPage();
    TCanvas* c10 = new TCanvas("c10", "hardSoft" ,800,600);
    c10->cd();
    hardSoft_Histo->Draw();
    c10->Update();

    //ps.NewPage();
    TCanvas* c11 = new TCanvas("c11", "collins_by_mustraal" ,800,600);
    c11->cd();
    collins_by_mustraal->SetTitle("collins_by_mustraal");
    collins_by_mustraal->GetZaxis()->SetRangeUser(0,2.);
    collins_by_mustraal->Draw("colz");
    c11->Update();

    //ps.NewPage();
    TCanvas* c12 = new TCanvas("c12", "WT Correlation KT KA", 800, 600);
    wt_corr_KT_KA->SetStats(0);
    wt_corr_KT_KA->SetTitle("");
    wt_corr_KT_KA->GetXaxis()->CenterTitle();
    wt_corr_KT_KA->GetYaxis()->CenterTitle();
    wt_corr_KT_KA->SetMarkerStyle(20);
    wt_corr_KT_KA->SetMarkerSize(0.5);
    wt_corr_KT_KA->SetMarkerColor(kBlue);
    c12->SetGrid();
    gPad->SetRightMargin(0.15);
    wt_corr_KT_KA->Draw("COLZ");
    TLatex texx;
    texx.SetTextSize(0.04);
    //texx.DrawLatexNDC(0.15, 0.8, Form("Correlation Factor: %.2f", wt_corr_KT_KA->GetCorrelationFactor()));
    c12->Update();

    TCanvas* c13 = new TCanvas("c13", "WT Correlation KT ST", 800, 600);
    wt_corr_KT_ST->SetStats(0);
    wt_corr_KT_ST->GetXaxis()->CenterTitle();
    wt_corr_KT_ST->GetYaxis()->CenterTitle();
    wt_corr_KT_ST->SetMarkerStyle(20);
    wt_corr_KT_ST->SetMarkerSize(0.5);
    wt_corr_KT_ST->SetMarkerColor(kBlue);
    c13->SetGrid();
    wt_corr_KT_ST->Draw("COLZ");
    TLatex texx2;
    texx2.SetTextSize(0.04);
    texx2.DrawLatexNDC(0.15, 0.8, Form("Correlation Factor: %.2f", wt_corr_KT_ST->GetCorrelationFactor()));
    c13->Update();

    TCanvas* c14 = new TCanvas("c14", "WT Correlation KT SA", 800, 600);
    wt_corr_KT_SA->SetStats(0);
    wt_corr_KT_SA->GetXaxis()->CenterTitle();
    wt_corr_KT_SA->GetYaxis()->CenterTitle();
    wt_corr_KT_SA->SetMarkerStyle(20);
    wt_corr_KT_SA->SetMarkerSize(0.5);
    wt_corr_KT_SA->SetMarkerColor(kBlue);
    c14->SetGrid();
    wt_corr_KT_SA->Draw("COLZ");
    TLatex texx3;
    texx3.SetTextSize(0.04);
    texx3.DrawLatexNDC(0.15, 0.8, Form("Correlation Factor: %.2f", wt_corr_KT_SA->GetCorrelationFactor()));
    c14->Update();

    TCanvas* c15 = new TCanvas("c15", "WT Correlation KA SA", 800, 600);
    wt_corr_KA_SA->SetStats(0);
    wt_corr_KA_SA->GetXaxis()->CenterTitle();
    wt_corr_KA_SA->GetYaxis()->CenterTitle();
    wt_corr_KA_SA->SetMarkerStyle(20);
    wt_corr_KA_SA->SetMarkerSize(0.5);
    wt_corr_KA_SA->SetMarkerColor(kBlue);
    c15->SetGrid();
    wt_corr_KA_SA->Draw("COLZ");
    TLatex texx4;
    texx4.SetTextSize(0.04);
    texx4.DrawLatexNDC(0.15, 0.8, Form("Correlation Factor: %.2f", wt_corr_KA_SA->GetCorrelationFactor()));
    c15->Update();

    TCanvas* c16 = new TCanvas("c16", "WT Correlation KA SA", 800, 600);
    wt_corr_ST_SA->SetStats(0);
    wt_corr_ST_SA->GetXaxis()->CenterTitle();
    wt_corr_ST_SA->GetYaxis()->CenterTitle();
    wt_corr_ST_SA->SetMarkerStyle(20);
    wt_corr_ST_SA->SetMarkerSize(0.5);
    wt_corr_ST_SA->SetMarkerColor(kBlue);
    c16->SetGrid();
    wt_corr_ST_SA->Draw("COLZ");
    TLatex texx5;
    texx5.SetTextSize(0.04);
    texx5.DrawLatexNDC(0.15, 0.8, Form("Correlation Factor: %.2f", wt_corr_ST_SA->GetCorrelationFactor()));
    c16->Update();


    TCanvas* c17 = new TCanvas("c17", "Theta Distribution", 800, 600);
    TPad *pad1 = new TPad("pad1", "Top pad", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper pad has no bottom margin
    pad1->Draw();
    pad1->cd();
    theta_dist_cs->SetStats(0);
    theta_dist_cs->GetXaxis()->SetTitle("#theta (radians)");
    theta_dist_cs->GetYaxis()->SetTitle("Normalized Count");
    //theta_dist_cs->SetTitle("Theta Distribution");
    theta_dist_cs->SetTitle("");
    theta_dist_cs->GetXaxis()->CenterTitle();
    theta_dist_cs->GetYaxis()->CenterTitle();
    theta_dist_cs->SetMarkerStyle(4);
    theta_dist_cs->SetMarkerSize(0.5);
    theta_dist_cs->SetMarkerColor(kBlack);
    theta_dist_cs->GetYaxis()->SetTitleSize(0.05);
    theta_dist_cs->GetYaxis()->SetLabelSize(0.05);
    theta_dist_cs->GetYaxis()->SetTitleOffset(1.0);
    theta_dist_cs->SetLineColor(kBlack);
    theta_dist_cs->GetYaxis()->ChangeLabel(1, /*angle*/-1, /*size*/0);
    f_model_cs->SetParameters(1, 0.97,0);  // initial guesses
    //f_model_cs->FixParameter(3, avg_m2OverE2); // Fix parameter [3] to 0.97
    // Normalize to unit area
    if (theta_dist_cs->Integral() > 0)
    theta_dist_cs->Scale(1.0 / theta_dist_cs->Integral("width")); // Normalize to unit area
    // 🔹 Fit histogram with your f_model
    theta_dist_cs->Fit(f_model_cs, "R0");   // "R" = restrict fit to function range, "Q" = quiet, "0" = no drawing
    theta_dist_cs->GetXaxis()->SetNdivisions(4, kFALSE);
    //theta_dist->GetYaxis()->SetRangeUser(0.2, theta_dist->GetMaximum()*1.2);

    theta_dist_ms->SetStats(0);
    theta_dist_ms->SetLineColor(kRed);
    theta_dist_ms->Scale(1.0 / theta_dist_ms->Integral("width")); // Normalize to unit area
    f_model_ms->SetParameters(1, 0.97,0);
    //f_model_ms->FixParameter(3, avg_m2OverE2); // Fix parameter [3] to 0.97
    theta_dist_ms->Fit(f_model_ms, "R0");   // "R" = restrict fit to function range, , "Q" = quiet, "0" = no drawing
    theta_dist_ms->SetMarkerStyle(20);
    theta_dist_ms->SetMarkerSize(0.5);
    theta_dist_ms->SetMarkerColor(kRed);

    TH1D* cs_markers = (TH1D*)theta_dist_cs->Clone("cs_markers");
    for (int i = 1; i <= cs_markers->GetNbinsX(); i++) cs_markers->SetBinError(i, 0);
    TH1D* ms_markers = (TH1D*)theta_dist_ms->Clone("ms_markers");
    for (int i = 1; i <= ms_markers->GetNbinsX(); i++) ms_markers->SetBinError(i, 0);

    theta_dist_cs->Draw("hist");  // with error bars
    cs_markers->Draw("P same");  // with error bars
    f_model_norm_cs->SetParameters(f_model_cs->GetParameter(0), avg_m2OverE2, avg_AFB);
    f_model_norm_ms->SetParameters(f_model_ms->GetParameter(0), avg_m2OverE2, avg_AFB);
    f_model_cs->SetLineColor(kBlack);
    f_model_cs->SetLineWidth(1);
    f_model_cs->Draw("same");  // overlay fit
    theta_dist_ms->Draw("hist same");
    ms_markers->Draw("P same");
    f_model_ms->SetLineColor(kRed);
    f_model_ms->SetLineWidth(1);
    f_model_ms->Draw("same");  // overlay fit
    normalizeGraph(Rtt);
    //Rtt->Draw("P");  JJ enable this to see korching function
    //theta_dist_cs->GetYaxis()->SetRangeUser(0.35, 0.8); //For 42.32 GeV
    theta_dist_cs->GetYaxis()->SetRangeUser(0.35, 0.65); //For 10.58 GeV
    //theta_dist_cs->GetYaxis()->SetRangeUser(0.3, 0.9); //For 91.18 GeV

    TLatex tex17;
    tex17.SetTextSize(0.04);
    tex17.DrawLatexNDC(0.25, 0.8, Form("Expected #LT #frac{m^{2}}{E^{2}} #GT = %.3f", avg_m2OverE2));
    tex17.DrawLatexNDC(0.25, 0.70, Form("CS Fitted #LT #frac{m^{2}}{E^{2}} #GT = %.3f #pm %.3f", f_model_cs->GetParameter(1), f_model_cs->GetParError(1)));
    tex17.DrawLatexNDC(0.25, 0.60, Form("MS Fitted #LT #frac{m^{2}}{E^{2}} #GT = %.3f #pm %.3f", f_model_ms->GetParameter(1), f_model_ms->GetParError(1)));

    c17->Update();

    // Add legend
    TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.05);
    legend->AddEntry(theta_dist_cs, "E_{#tau#tau} < 6 GeV", "");
    legend->AddEntry(theta_dist_cs, "Collins-Soper frame", "lp");  // l = line, p = marker
    legend->AddEntry(theta_dist_ms, "Mustraal frame",     "lp");
    legend->SetTextFont(42);   // nice readable font
    legend->Draw();




    pad1->SetGrid();
    gPad->SetLeftMargin(0.12);
    pad1->SetFrameBorderMode(0); // or
    gPad->Modified(); gPad->Update();

    c17->cd();
    TPad *pad2 = new TPad("pad2", "Bottom pad", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3); // ratio plot needs bigger bottom margin
    pad2->Draw();
    pad2->cd();
    //TH1D* ratio1 = (TH1D*)theta_dist_cs->Clone("ratio1");
    TF1* ratio1 = MakeRatioTF1("h1", *f_model_cs, *f_model_norm_cs, -1., 1., 1e-12);
    ratio1->SetTitle(""); // Remove the title for the ratio plot
    ratio1->GetYaxis()->SetTitle("Fit/Model");
    ratio1->GetXaxis()->SetTitle("cos #theta");
    ratio1->GetXaxis()->CenterTitle();
    ratio1->GetYaxis()->CenterTitle();
    ratio1->GetYaxis()->SetTitleSize(0.15);
    ratio1->GetYaxis()->SetLabelSize(0.15);
    ratio1->GetXaxis()->SetTitleSize(0.15);
    ratio1->GetXaxis()->SetLabelSize(0.15);
    ratio1->GetYaxis()->SetTitleOffset(0.35);
    ratio1->GetXaxis()->SetNdivisions(4, kFALSE);
    ratio1->GetYaxis()->SetNdivisions(505);
    ratio1->GetYaxis()->ChangeLabel(-1, /*angle*/-1, /*size*/0); // -1 = last (top) label
    //TH1D* RttNew = graphToHist(Rtt, ratio1);
    //ratio1->Divide(RttNew);
    //ratio1->Divide(f_model_norm_cs);
    //ratio1->GetYaxis()->SetRangeUser(0.96, 1.1); //For 42.32 GeV and 91.18 GeV
    ratio1->GetYaxis()->SetRangeUser(0.96, 1.04); //For 10.58 GeV
    ratio1->SetLineColor(kBlack);
    ratio1->Draw("hist");

    //TH1D* ratio2 = (TH1D*)theta_dist_ms->Clone("ratio2");
    TF1* ratio2 = MakeRatioTF1("h2", *f_model_ms, *f_model_norm_ms, -1., 1., 1e-12);
    ratio2->SetTitle(""); // Remove the title for the ratio plot
    ratio2->GetYaxis()->SetTitle("Histo/model");
    ratio2->GetXaxis()->SetTitle("#theta (radians)");
    ratio2->GetXaxis()->CenterTitle();
    ratio2->GetYaxis()->CenterTitle();
    ratio2->GetYaxis()->SetTitleSize(0.1);
    ratio2->GetYaxis()->SetLabelSize(0.1);
    ratio2->GetXaxis()->SetTitleSize(0.1);
    ratio2->GetXaxis()->SetLabelSize(0.1);
    //ratio2->GetYaxis()->SetNdivisions(505);
    ratio2->SetLineColor(kRed);
    //ratio2->Divide(f_model_norm_ms);
    ratio2->Draw("hist same");
    c17->Update();
    // Horizontal guide at y=1 across the full x-range
    TLine* y2 = new TLine(-1.0, 1.0, 1.0, 1.0);
    y2->SetLineStyle(2);       // dashed
    y2->SetLineWidth(1);
    // y2->SetLineColor(kGray+2); // optional
    y2->Draw("same");
    pad2->SetGridx();
    gPad->SetLeftMargin(0.12);
    pad2->SetFrameBorderMode(0); // removes the frame line; pick one pads
    gPad->Modified(); gPad->Update(); // Update the pad to reflect the new axis

    TCanvas* c17b = new TCanvas("c17b", "Phi Distribution", 800, 600);
    phi_dist_cs->SetStats(0);
    phi_dist_cs->GetXaxis()->SetTitle("#phi (radians)");
    phi_dist_cs->GetYaxis()->SetTitle("Events");
    phi_dist_cs->SetTitle("Phi Distribution");
    phi_dist_cs->GetXaxis()->CenterTitle();
    phi_dist_cs->GetYaxis()->CenterTitle();
    phi_dist_cs->SetMarkerStyle(20);
    phi_dist_cs->SetMarkerSize(0.5);
    phi_dist_cs->SetLineColor(kRed);
    c17b->SetGrid();
    phi_dist_cs->Draw("hist");
    phi_dist_ms->SetLineColor(kBlue);
    phi_dist_ms->Draw("hist same");
    c17b->Update();

    TCanvas* c18 = new TCanvas("c18", "Invariant Mass Distribution", 800, 600);
    Invariant_dist->SetStats(0);
    Invariant_dist->GetXaxis()->SetTitle("Invariant Mass (GeV)");
    Invariant_dist->GetYaxis()->SetTitle("Events");
    Invariant_dist->SetTitle("Invariant Mass Distribution");
    Invariant_dist->GetXaxis()->CenterTitle();
    Invariant_dist->GetYaxis()->CenterTitle();
    Invariant_dist->SetMarkerStyle(20);
    Invariant_dist->SetMarkerSize(0.5);
    Invariant_dist->SetMarkerColor(kMagenta);
    c18->SetGrid();
    Invariant_dist->Draw("");
    c18->Update();



    TCanvas* c19 = new TCanvas("c19", "Tau Z Momenta", 800, 600);
    tauZMomenta->SetStats(0);
    tauZMomenta->GetXaxis()->SetTitle("Momentum (GeV)");
    tauZMomenta->GetYaxis()->SetTitle("Events");
    tauZMomenta->SetTitle("Tau Z Momenta Distribution");
    tauZMomenta->GetXaxis()->CenterTitle();
    tauZMomenta->GetYaxis()->CenterTitle();
    tauZMomenta->SetMarkerStyle(20);
    tauZMomenta->SetMarkerSize(0.5);
    tauZMomenta->SetMarkerColor(kMagenta);
    c19->SetGrid();
    tauZMomenta->Draw("");
    c19->Update();

    TCanvas* c20 = new TCanvas("c20", "Rtt Plot", 800, 600);
    Rtt->SetTitle("Rtt Plot");
    Rtt->GetXaxis()->SetTitle("cos(theta)");
    Rtt->GetYaxis()->SetTitle("Rtt");
    Rtt->GetXaxis()->CenterTitle();
    Rtt->GetYaxis()->CenterTitle();
    Rtt->SetMarkerStyle(20);
    Rtt->SetMarkerSize(0.5);
    Rtt->SetMarkerColor(kMagenta);
    c20->SetGrid();
    Rtt->Draw("AP");
    c20->Update();

    TCanvas* c21 = new TCanvas("c21", "wtw Plot", 800, 600);
    c21->Divide(2,2);
    c21->cd(1);
    wtwPP->SetTitle("wtwPP Plot");
    wtwPP->GetXaxis()->SetTitle("wtwPP");
    wtwPP->GetYaxis()->SetTitle("Events");
    wtwPP->GetXaxis()->CenterTitle();
    wtwPP->GetYaxis()->CenterTitle();
    wtwPP->SetMarkerStyle(20);
    wtwPP->SetMarkerSize(0.5);
    wtwPP->Draw("hist");
    c21->cd(2);
    wtwMM->SetTitle("wtwMM Plot");
    wtwMM->GetXaxis()->SetTitle("wtwMM");
    wtwMM->GetYaxis()->SetTitle("Events");
    wtwMM->GetXaxis()->CenterTitle();
    wtwMM->GetYaxis()->CenterTitle();
    wtwMM->SetMarkerStyle(20);
    wtwMM->SetMarkerSize(0.5);
    wtwMM->Draw("hist");
    c21->cd(3);
    wtwPM->SetTitle("wtwPM Plot");
    wtwPM->GetXaxis()->SetTitle("wtwPM");
    wtwPM->GetYaxis()->SetTitle("Events");
    wtwPM->GetXaxis()->CenterTitle();
    wtwPM->GetYaxis()->CenterTitle();
    wtwPM->SetMarkerStyle(20);
    wtwPM->SetMarkerSize(0.5);
    wtwPM->Draw("hist");
    c21->cd(4);
    wtwMP->SetTitle("wtwMP Plot");
    wtwMP->GetXaxis()->SetTitle("wtwMP");
    wtwMP->GetYaxis()->SetTitle("Events");
    wtwMP->GetXaxis()->CenterTitle();
    wtwMP->GetYaxis()->CenterTitle();
    wtwMP->SetMarkerStyle(20);
    wtwMP->SetMarkerSize(0.5);
    wtwMP->Draw("hist");
    c21->Update();

    TCanvas* c22 = new TCanvas("c22", "Polarimetric vector", 800, 600);
    // H1_H2->SetMarkerColor(kGreen);
    // H1_H2->Draw("P");

    //ps.Close();

    TString outFile = Form("output_histos_%.2fGeV.pdf", beamEnergy);

    // after you draw c1..c9
    c1->Print(outFile + "[");  // open
    c1->Print(outFile);   // page 1
    c2->Print(outFile);   // page 2
    c2b->Print(outFile);   // page 2b //JMJ added
    c3->Print(outFile);
    c4->Print(outFile);
    c5->Print(outFile);
    c6->Print(outFile);
    c7->Print(outFile);
    c8->Print(outFile);
    c9->Print(outFile);
    c10->Print(outFile);
    c11->Print(outFile);
    c12->Print(outFile);
    c13->Print(outFile);
    c14->Print(outFile);
    c15->Print(outFile);
    c16->Print(outFile);
    c17->Print(outFile);
    c17b->Print(outFile);
    c18->Print(outFile);
    c19->Print(outFile);
    c20->Print(outFile);
    c21->Print(outFile);
    c22->Print(outFile);
    c22->Print(outFile + "]");  // close

    c17->SaveAs("ThetaDistribution_10_58.png");
    c1->SaveAs("KKMCee_chiral_10_58_ISR_FSR.png");
    c2->SaveAs("External_Chiral_10_58_ISR_FSR.png");
    c3->SaveAs("KKMCee_WT_Correlation_10_58.png");
    c12->SaveAs("External_WT_Correlation_10_58.png");

    currentDir->cd();
}
