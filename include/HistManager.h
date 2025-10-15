#pragma once
#include <array>
#include <string>
#include <TH1D.h>
#include <TH2F.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraph2D.h>


class HistManager {
public:
  HistManager() = default;
  ~HistManager();
  // Histograms
    TH2F* helicity_corr = new TH2F("helicity_corr", "Helicity1 vs Helicity2;Helicity1;Helicity2", 2, -1.1, 1.1, 2, -1.1, 1.1);
    TH2F* helicity_corr_standalone = new TH2F("helicity_corr_standalone", "Helicity1 vs Helicity2 (External);Helicity1;Helicity2", 2, -1.1, 1.1, 2, -1.1, 1.1);
    TH2F* helicity_corr_st_mustraal = new TH2F("helicity_corr_st_mustraal", "Helicity1 vs Helicity2 (External Mustraal);Helicity1;Helicity2", 3, -1.5, 1.5, 3, -1.5, 1.5);
    TH2F* wt_corr_KT_KA = new TH2F("wt_corr_KT_KA", "wt correlation (KT-KA);wt_KT;wt_KA", 100, 0, 4.0, 100, 0, 4.0);
    TH2F* wt_corr_KT_ST = new TH2F("wt_corr_KT_ST", "wt correlation (KT-ST);wt_KT;wt_ST", 100, 0, 4.0, 100, 0, 4.0);
    TH2F* wt_corr_KT_SA = new TH2F("wt_corr_KT_SA", "wt correlation (KT-SA);wt_KT;wt_SA", 100, 0, 4.0, 100, 0, 4.0);
    TH2F* wt_corr_KA_SA = new TH2F("wt_corr_KA_SA", "wt correlation (KA-SA);wt_KA;wt_SA", 100, 0, 4.0, 100, 0, 4.0);
    TH2F* wt_corr_ST_SA = new TH2F("wt_corr_ST_SA", "wt correlation (ST-SA);wt_ST;wt_SA", 100, 0, 4.0, 100, 0, 4.0);
    TH2F* acoplanarity_angle1 = new TH2F("acoplanarity_angle1", "Acoplanarity angle;#theta_{1};#theta_{2}", 100, 0.5, TMath::Pi()/2., 100, 0, 4.0);
    TH1D* h_angle1 = new TH1D("angle1", "Acoplanarity angle1", 100, 0, 2.0*TMath::Pi());
    TH1D* h_angle_01 = new TH1D("angle_01", "Acoplanarity angle01", 100, 0, 2.0*TMath::Pi());
    TH2F* acoplanarity_angle2 = new TH2F("acoplanarity_angle2", "Acoplanarity angle;#theta_{1};#theta_{2}", 100, 0.5, TMath::Pi()/2., 100, 0, 4.0);
    TH1D* h_angle2 = new TH1D("angle2", "Acoplanarity angle2", 100, 0, 2.0*TMath::Pi());
    TH1D* h_angle_02 = new TH1D("angle_02", "Acoplanarity angle02", 100, 0, 2.0*TMath::Pi());
    TH2F* acoplanarity_angle3 = new TH2F("acoplanarity_angle3", "Acoplanarity angle;#theta_{1};#theta_{2}", 100, 0.5, TMath::Pi()/2., 100, 0, 4.0);
    TH1D* h_angle3 = new TH1D("angle3", "Acoplanarity angle3", 100, 0, 2.0*TMath::Pi());
    TH1D* h_angle_03 = new TH1D("angle_03", "Acoplanarity angle03", 100, 0, 2.0*TMath::Pi());
    TH1D* hardSoft_Histo = new TH1D("hardSoft", "hardSoft",100,0.,2.);
    TH2F* collins_soper_corr_soft = new TH2F("collins_soper_corr_soft", "Collins Soper Correlation;wt_KKMCee;wt_StandAlone", 100, 0, 4.0, 100, 0, 4.0);
    TH2F* p1_frame_corr_soft = new TH2F("p1_frame_corr_soft", "P1 Frame Correlation;wt_KKMCee;wt_StandAlone", 100, 0, 4.0, 100, 0, 4.0);
    TH2F* p2_frame_corr_soft = new TH2F("p2_frame_corr_soft", "P2 Frame Correlation;wt_KKMCee;wt_StandAlone", 100, 0, 4.0, 100, 0, 4.0);
    TH2F* mustraal_corr_soft = new TH2F("mustraal_corr_soft", "Mustraal Correlation;wt_KKMCee;wt_StandAlone", 100, 0, 4.0, 100, 0, 4.0);
    TH2F* collins_soper_corr_hard = new TH2F("collins_soper_corr_hard", "Collins Soper Correlation;wt_KKMCee;wt_StandAlone", 100, 0, 4.0, 100, 0, 4.0);
    TH2F* p1_frame_corr_hard = new TH2F("p1_frame_corr_hard", "P1 Frame Correlation;wt_KKMCee;wt_StandAlone", 100, 0, 4.0, 100, 0, 4.0);
    TH2F* p2_frame_corr_hard = new TH2F("p2_frame_corr_hard", "P2 Frame Correlation;wt_KKMCee;wt_StandAlone", 100, 0, 4.0, 100, 0, 4.0);
    TH2F* mustraal_corr_hard = new TH2F("mustraal_corr_hard", "Mustraal Correlation;wt_KKMCee;wt_StandAlone", 100, 0, 4.0, 100, 0, 4.0);

    TH1D* theta_dist_cs = new TH1D("theta_dist_cs", "Theta Distribution;#theta;Events", 100, -1.0, 1.0);
    TH1D* theta_dist_ms = new TH1D("theta_dist_ms", "Theta Distribution;#theta;Events", 100, -1.0, 1.0);
    TH1D* Invariant_dist = new TH1D("Invariant_dist", "Invariant Mass Distribution;Invariant Mass (GeV);Events", 100, 0.0, 6.0);
    TH1D* directTheta_dist = new TH1D("directTheta_dist", "Direct Theta Distribution;#theta;Events", 100, -1.0,1.0);

    TH1D* phi_dist_cs = new TH1D("phi_dist_cs", "Phi Distribution;#phi;Events", 100, 0.0, 2.0*TMath::Pi());
    TH1D* phi_dist_ms = new TH1D("phi_dist_ms", "Phi Distribution;#phi;Events", 100, 0.0, 2.0*TMath::Pi());

    TH1D* wtwPP = new TH1D("wtwPP", "wtwPP", 10000, 0.0, 1.0);
    TH1D* wtwMM = new TH1D("wtwMM", "wtwMM", 10000, 0.0, 1.0);
    TH1D* wtwPM = new TH1D("wtwPM", "wtwPM", 10000, 0.0, 1.0);
    TH1D* wtwMP = new TH1D("wtwMP", "wtwMP", 10000, 0.0, 1.0);

    double alpha = 7.2973525693e-3;
    double prefactor = 4 * TMath::Pi() * TMath::Pi() * alpha * alpha;

    TF1* f_model_cs = new TF1("f_model_cs",Form("%f * ([0]*(1  + x*x + [1]*(1-x*x) + [2]*x))", prefactor),-1, 1);
    TF1* f_model_ms = new TF1("f_model_ms",Form("%f * ([0]*(1  + x*x + [1]*(1-x*x) + [2]*x))", prefactor),-1, 1);

    TF1* f_model_norm_cs = new TF1("f_model_norm_cs",Form("%f* ([0]*(1  + x*x + [1]*(1-x*x) + [2]*x))", prefactor),-1, 1);
    TF1* f_model_norm_ms = new TF1("f_model_norm_ms",Form("%f* ([0]*(1  + x*x + [1]*(1-x*x) + [2]*x))", prefactor),-1, 1);

    TGraph* Rtt = new TGraph();

    TH1D* tauZMomenta = new TH1D("tauZMomenta", "Tau Z Momenta;Momentum (GeV);Events", 100, -6.0, 6.0);

    TGraph2D *H1= new TGraph2D();
    TGraph2D *H2= new TGraph2D();
    TGraph2D *H1_H2= new TGraph2D();





    TH1D* ratio_angle1 = nullptr;
    TH1D* ratio_angle2 = nullptr;
    TH1D* ratio_angle3 = nullptr;
    TH2F* collins_by_mustraal = nullptr;

    void DrawHistograms(double beamEnergy, double avg_m2OverE2, double avg_AFB);

    TH1D *graphToHist(const TGraph *gr, const TH1D *href) {
     TH1D *hout = (TH1D*)href->Clone("hout");
     hout->Reset();
      for (int i = 1; i <= href->GetNbinsX(); i++) {
        double x = href->GetBinCenter(i);
        double y = gr->Eval(x);   // interpolate y value at bin center
        hout->SetBinContent(i, y);
      }
      return hout;
    }

    void normalizeGraph(TGraph* g) {
      int n = g->GetN();
      double* x = g->GetX();
      double* y = g->GetY();
      double dx = x[2] - x[1];
      // Numerical integration (trapezoidal rule)
      double integral = 0.0;
      for (int i = 1; i < n ; i++) {
        double x1, y1;
        g->GetPoint(i, x1, y1);
        integral += y1 * dx;
      }

      // Scale y values
      if (integral != 0) {
        for (int i = 1; i <= n; i++) {
            double x1, y1;
            g->GetPoint(i, x1, y1);
            g->SetPoint(i, x1, y1 / integral) ; // Normalize over [-1,1]
        }
      }
    }

    TF1* MakeRatioTF1(const char* name,
                 TF1& f, TF1& g,
                 double xmin, double xmax,
                 double zero_guard = 0.0)   // e.g. 1e-12 to guard small g(x)
    {
      // restrict to overlap of domains
      double xlo = std::max({xmin, f.GetXmin(), g.GetXmin()});
      double xhi = std::min({xmax, f.GetXmax(), g.GetXmax()});
      if (xlo >= xhi) throw std::runtime_error("No overlapping x-range between f and g");

      auto ratio = [&f, &g, zero_guard](double* x, double*) -> double {
        double gx = g.EvalPar(x);
        if ((zero_guard > 0.0 && std::abs(gx) < zero_guard) || gx == 0.0)
        return std::numeric_limits<double>::quiet_NaN();
        return f.EvalPar(x) / gx;
      };

      return new TF1(name, ratio, xlo, xhi, 0);  // npar=0
    }

};