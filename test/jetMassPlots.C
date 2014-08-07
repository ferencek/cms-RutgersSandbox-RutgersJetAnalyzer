#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "exoStyle.C"


using namespace std;


void jet_mass(const string& fFile, const string& fFileDir,
              const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
              const double fXmin, const double fXmax, const string& fOutputFile, const int fRebin=2)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  //gStyle->SetOptStat("nemruoi");
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetStatX(0.94);
  gStyle->SetStatY(0.93);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);
  //gStyle->UseCurrentStyle();
  gROOT->ForceStyle();

  TFile *file = new TFile(fFile.c_str());

  TH2D *h2_nPV_JetMass = (TH2D*)file->Get((fFileDir + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH1D *h1_JetMass = h2_nPV_JetMass->ProjectionY("_py",fPVLow,fPVHigh);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  h1_JetMass->Rebin(fRebin);
  h1_JetMass->GetXaxis()->SetTitle("m_{jet} (pruned) [GeV/c^{2}]");
  h1_JetMass->GetYaxis()->SetTitle("Entries");
  h1_JetMass->GetXaxis()->SetRangeUser(fXmin,fXmax);
  h1_JetMass->SetTitleOffset(0.95,"X");
  h1_JetMass->SetTitleOffset(1.2,"Y");
  h1_JetMass->SetLineColor(kBlue+2);
  h1_JetMass->SetLineWidth(2);
  h1_JetMass->SetMaximum(1.17*h1_JetMass->GetMaximum());

  Double_t fitXmin=95.;
  Double_t fitXmax=125;
  if(fFile.find("RadionToHH")!=string::npos && fFile.find("R1p2")!=string::npos)
  {
    fitXmin=110.;
    fitXmax=155;
  }
  else if(fFile.find("RadionToHH")!=string::npos)
  {
    fitXmin=105.;
    fitXmax=135;
  }

  TF1 *f1 = new TF1("f1","gaus",fitXmin,fitXmax);
  f1->SetLineColor(kRed);
  f1->SetLineWidth(2);

  if(fTitle.find("QCD")==string::npos) h1_JetMass->Fit("f1","R");

  h1_JetMass->Draw("hists");
  if(fTitle.find("QCD")==string::npos) f1->Draw("same");

  TLatex l1;
  l1.SetTextAlign(13);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.045);
  if(fTitle.find("QCD")!=string::npos) l1.DrawLatex(0.17,0.90, fTitle.c_str());
  else                                 l1.DrawLatex(0.53,0.80, fTitle.c_str());

  l1.SetTextAlign(12);
  l1.SetTextSize(0.05);
  l1.SetTextFont(62);
  l1.DrawLatex(0.14,0.97, "CMS Simulation Preliminary, #sqrt{s} = 8 TeV");
  //l1.DrawLatex(0.14,0.97, "CMS Simulation");
  //l1.SetTextFont(42);
  //l1.DrawLatex(0.14+0.35,0.97, "#sqrt{s} = 8 TeV");

  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void makePlots()
{
  // RadionToHH_4b_M-600
  // nPV inclusive
  jet_mass("ROOT_files/jobs_6_forApproval/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root", "jetAnalyzerFatJets_PrunedSubjets",
           "Pt200toInf", 0, 52, "#splitline{H(125)#rightarrowb#bar{b}, AK R=1.2}{p_{T}>200 GeV, #DeltaR(H,jet)<0.5}",
           0, 299.5, "Jet_mass_AK12_pruned_Pt200toInf_RadionToHH_4b_M-600.eps", 1);

  // RadionToHH_4b_M-800
  // nPV inclusive
  jet_mass("ROOT_files_AK8/RadionToHH_4b_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "jetAnalyzerFatJets_PrunedSubjets",
           "Pt300toInf", 0, 52, "#splitline{H(125)#rightarrowb#bar{b}, AK R=0.8}{p_{T}>300 GeV, #DeltaR(H,jet)<0.5}",
           0, 299.5, "Jet_mass_AK8_pruned_Pt300toInf_RadionToHH_4b_M-800.eps", 1);
  
  // BprimeBprimeToBHBHinc_M-1000
  // nPV inclusive
  jet_mass("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "jetAnalyzerFatJets_PrunedSubjets",
           "Pt300toInf", 0, 52, "#splitline{H(120)#rightarrowb#bar{b}, AK R=0.8}{p_{T}>300 GeV, #DeltaR(H,jet)<0.5}",
           0, 299.5, "Jet_mass_AK8_pruned_Pt300toInf_BprimeBprimeToBHBHinc_M-1000.eps", 1);
 
  // QCDPythia6
  // nPV inclusive
  jet_mass("ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "jetAnalyzerFatJets_PrunedSubjets",
           "Pt300toInf", 0, 52, "QCD, AK R=0.8, p_{T}>300 GeV/c",
           0, 299.5, "Jet_mass_AK8_pruned_Pt300toInf_QCDPythia6.eps", 1);

}
