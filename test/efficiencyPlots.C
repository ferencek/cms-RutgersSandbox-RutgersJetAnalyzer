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
#include "TGraphAsymmErrors.h"
#include "exoStyle.C"
#include "CMS_lumi.C"


using namespace std;


void efficiency_curves_grooming(const string& fFileS, const string& fFileB,
                                const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
                                const string& fLeg1, const string& fLeg2, const string& fLeg3, const string& fLeg4,
                                const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  // signal file
  TFile *file_S = new TFile(fFileS.c_str());

  // background file
  TFile *file_B = new TFile(fFileB.c_str());

  // signal histograms
  TH2D *h2_nPV_JetMass_S          = (TH2D*)file_S->Get(("jetAnalyzerDefaultJetMass/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_Filtered = (TH2D*)file_S->Get(("jetAnalyzerFilteredJetMass/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_Pruned   = (TH2D*)file_S->Get(("jetAnalyzerPrunedJetMass/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_Trimmed  = (TH2D*)file_S->Get(("jetAnalyzerTrimmedJetMass/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S          = (TH2D*)file_S->Get(("jetAnalyzerDefaultJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_Filtered = (TH2D*)file_S->Get(("jetAnalyzerFilteredJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_Pruned   = (TH2D*)file_S->Get(("jetAnalyzerPrunedJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_Trimmed  = (TH2D*)file_S->Get(("jetAnalyzerTrimmedJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B          = (TH2D*)file_B->Get(("jetAnalyzerDefaultJetMass/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_Filtered = (TH2D*)file_B->Get(("jetAnalyzerFilteredJetMass/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_Pruned   = (TH2D*)file_B->Get(("jetAnalyzerPrunedJetMass/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_Trimmed  = (TH2D*)file_B->Get(("jetAnalyzerTrimmedJetMass/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B          = (TH2D*)file_B->Get(("jetAnalyzerDefaultJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_Filtered = (TH2D*)file_B->Get(("jetAnalyzerFilteredJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_Pruned   = (TH2D*)file_B->Get(("jetAnalyzerPrunedJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_Trimmed  = (TH2D*)file_B->Get(("jetAnalyzerTrimmedJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // signal denominator counts
  double denom_S          = h2_nPV_JetMass_S->Integral(fPVLow,fPVHigh,0,401);
  double denom_S_Filtered = h2_nPV_JetMass_S_Filtered->Integral(fPVLow,fPVHigh,0,401);
  double denom_S_Pruned   = h2_nPV_JetMass_S_Pruned->Integral(fPVLow,fPVHigh,0,401);
  double denom_S_Trimmed  = h2_nPV_JetMass_S_Trimmed->Integral(fPVLow,fPVHigh,0,401);

  // background denominator counts
  double denom_B = h2_nPV_JetMass_B->Integral(fPVLow,fPVHigh,0,401);
  double denom_B_Filtered = h2_nPV_JetMass_B_Filtered->Integral(fPVLow,fPVHigh,0,401);
  double denom_B_Pruned   = h2_nPV_JetMass_B_Pruned->Integral(fPVLow,fPVHigh,0,401);
  double denom_B_Trimmed  = h2_nPV_JetMass_B_Trimmed->Integral(fPVLow,fPVHigh,0,401);

  // Default jets
  TGraph *g_eff = new TGraph(21);
  g_eff->SetName("g_eff");
  g_eff->SetLineColor(kBlack);
  g_eff->SetLineWidth(2);
  g_eff->SetLineStyle(1);
  g_eff->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff->SetPoint(i,(num_B/denom_B),(num_S/denom_S));
  }

  // Filtered jets
  TGraph *g_eff_Filtered = new TGraph(21);
  g_eff_Filtered->SetName("g_eff_Filtered");
  g_eff_Filtered->SetLineColor(kBlue);
  g_eff_Filtered->SetLineWidth(2);
  g_eff_Filtered->SetLineStyle(7);
  g_eff_Filtered->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_Filtered->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_Filtered->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_Filtered->SetPoint(i,(num_B/denom_B_Filtered),(num_S/denom_S_Filtered));
  }

  // Pruned jets
  TGraph *g_eff_Pruned = new TGraph(21);
  g_eff_Pruned->SetName("g_eff_Pruned");
  g_eff_Pruned->SetLineColor(kRed);
  g_eff_Pruned->SetLineWidth(2);
  g_eff_Pruned->SetLineStyle(3);
  g_eff_Pruned->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_Pruned->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_Pruned->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_Pruned->SetPoint(i,(num_B/denom_B_Pruned),(num_S/denom_S_Pruned));
  }

  // Trimmed jets
  TGraph *g_eff_Trimmed = new TGraph(21);
  g_eff_Trimmed->SetName("g_eff_Trimmed");
  g_eff_Trimmed->SetLineColor(kGreen+2);
  g_eff_Trimmed->SetLineWidth(2);
  g_eff_Trimmed->SetLineStyle(8);
  g_eff_Trimmed->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_Trimmed->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_Trimmed->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_Trimmed->SetPoint(i,(num_B/denom_B_Trimmed),(num_S/denom_S_Trimmed));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg",";Mistag rate; Tagging efficiency",100,fXmin,fXmax,100,0,0.9);
  bkg->SetTitleOffset(1.05,"Y");
  bkg->Draw();

  g_eff->Draw("L");
  g_eff_Filtered->Draw("L");
  g_eff_Pruned->Draw("L");
  g_eff_Trimmed->Draw("L");

  TLegend *legend = new TLegend(.55,.25,.85,.45);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_eff_Trimmed, fLeg1.c_str(),"l");
  legend->AddEntry(g_eff_Filtered, fLeg2.c_str(),"l");
  legend->AddEntry(g_eff_Pruned, fLeg3.c_str(),"l");
  legend->AddEntry(g_eff, fLeg4.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.14,0.96, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete g_eff;
  delete g_eff_Filtered;
  delete g_eff_Pruned;
  delete g_eff_Trimmed;
  delete file_S;
  delete file_B;
}


void efficiency_curves_comp(const string& fFileS, const string& fFileB, const string& fFileDirS1, const string& fFileDirS2, const string& fFileDirB1, const string& fFileDirB2,
                            const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
                            const string& fLeg1, const string& fLeg2, const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  // signal file
  TFile *file_S  = new TFile(fFileS.c_str());

  // background file
  TFile *file_B = new TFile(fFileB.c_str());

  // signal histograms
  TH2D *h2_nPV_JetMass_S_1 = (TH2D*)file_S->Get((fFileDirS1 + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_2 = (TH2D*)file_S->Get((fFileDirS2 + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S_1 = (TH2D*)file_S->Get((fFileDirS1 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_2 = (TH2D*)file_S->Get((fFileDirS2 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B_1 = (TH2D*)file_B->Get((fFileDirB1 + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_2 = (TH2D*)file_B->Get((fFileDirB2 + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B_1 = (TH2D*)file_B->Get((fFileDirB1 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_2 = (TH2D*)file_B->Get((fFileDirB2 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // signal denominator counts
  double denom_S_1 = h2_nPV_JetMass_S_1->Integral(fPVLow,fPVHigh,0,401);
  double denom_S_2 = h2_nPV_JetMass_S_2->Integral(fPVLow,fPVHigh,0,401);

  // background denominator counts
  double denom_B_1 = h2_nPV_JetMass_B_1->Integral(fPVLow,fPVHigh,0,401);
  double denom_B_2 = h2_nPV_JetMass_B_2->Integral(fPVLow,fPVHigh,0,401);

  // Default N-subjettiness
  TGraph *g_eff_1 = new TGraph(21);
  g_eff_1->SetName("g_eff_1");
  g_eff_1->SetLineColor(kGreen+2);
  g_eff_1->SetLineWidth(2);
  g_eff_1->SetLineStyle(1);
  g_eff_1->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_1->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_1->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_1->SetPoint(i,(num_B/denom_B_1),(num_S/denom_S_1));
  }

  // Trimmed N-subjettiness
  TGraph *g_eff_2 = new TGraph(21);
  g_eff_2->SetName("g_eff_2");
  g_eff_2->SetLineColor(kRed);
  g_eff_2->SetLineWidth(2);
  g_eff_2->SetLineStyle(2);
  g_eff_2->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_2->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_2->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_2->SetPoint(i,(num_B/denom_B_2),(num_S/denom_S_2));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg",";Mistag rate; Tagging efficiency",100,fXmin,fXmax,100,0,0.9);
  bkg->SetTitleOffset(1.05,"Y");
  bkg->Draw();

  g_eff_1->Draw("L");
  g_eff_2->Draw("L");

  TLegend *legend = new TLegend(.45,.25,.75,.45);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_eff_1, fLeg1.c_str(),"l");
  legend->AddEntry(g_eff_2, fLeg2.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.035);
  l1.DrawLatex(0.06,0.96, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete g_eff_2;
  delete g_eff_1;
  delete file_S;
  delete file_B;
}


void efficiency_curve(const string& fFileS, const string& fFileB, const string& fFileDirS, const string& fFileDirB,
                            const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
                            const string& fLeg, const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  // signal file
  TFile *file_S  = new TFile(fFileS.c_str());

  // background file
  TFile *file_B = new TFile(fFileB.c_str());

  // signal histograms
  TH2D *h2_nPV_JetMass_S = (TH2D*)file_S->Get((fFileDirS + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S = (TH2D*)file_S->Get((fFileDirS + "/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B = (TH2D*)file_B->Get((fFileDirB + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B = (TH2D*)file_B->Get((fFileDirB + "/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // signal denominator counts
  double denom_S = h2_nPV_JetMass_S->Integral(fPVLow,fPVHigh,0,401);

  // background denominator counts
  double denom_B = h2_nPV_JetMass_B->Integral(fPVLow,fPVHigh,0,401);

  TGraph *g_eff_1 = new TGraph(21);
  g_eff_1->SetName("g_eff_1");
  g_eff_1->SetLineColor(kGreen+2);
  g_eff_1->SetLineWidth(2);
  g_eff_1->SetLineStyle(1);
  g_eff_1->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_1->SetPoint(i,(num_B/denom_B),(num_S/denom_S));
  }


  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg",";Mistag rate; Tagging efficiency",100,fXmin,fXmax,100,0,0.9);
  bkg->SetTitleOffset(1.05,"Y");
  bkg->Draw();

  g_eff_1->Draw("L");

  TLegend *legend = new TLegend(.45,.25,.75,.45);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_eff_1, fLeg.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.035);
  l1.DrawLatex(0.06,0.96, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete g_eff_1;
  delete file_S;
  delete file_B;
}


void efficiency_curves_nsj_massdrop(const string& fFileS, const string& fFileB, const string& fFileDir1, const string& fFileDir2,
                                    const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
                                    const string& fLeg1, const string& fLeg2, const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  // signal file
  TFile *file_S  = new TFile(fFileS.c_str());

  // background file
  TFile *file_B = new TFile(fFileB.c_str());

  // signal histograms
  TH2D *h2_nPV_JetMass_S_1 = (TH2D*)file_S->Get((fFileDir1 + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_2 = (TH2D*)file_S->Get((fFileDir2 + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S = (TH2D*)file_S->Get((fFileDir1 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_MassDrop_S = (TH2D*)file_S->Get((fFileDir2 + "/h2_nPV_MassDrop_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B_1 = (TH2D*)file_B->Get((fFileDir1 + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_2 = (TH2D*)file_B->Get((fFileDir2 + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B = (TH2D*)file_B->Get((fFileDir1 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_MassDrop_B = (TH2D*)file_B->Get((fFileDir2 + "/h2_nPV_MassDrop_" + fPtRange).c_str());

  // signal denominator counts
  double denom_S_1 = h2_nPV_JetMass_S_1->Integral(fPVLow,fPVHigh,0,401);
  double denom_S_2 = h2_nPV_JetMass_S_2->Integral(fPVLow,fPVHigh,0,401);

  // background denominator counts
  double denom_B_1 = h2_nPV_JetMass_B_1->Integral(fPVLow,fPVHigh,0,401);
  double denom_B_2 = h2_nPV_JetMass_B_2->Integral(fPVLow,fPVHigh,0,401);

  // Default N-subjettiness
  TGraph *g_eff_1 = new TGraph(21);
  g_eff_1->SetName("g_eff_1");
  g_eff_1->SetLineColor(kGreen+2);
  g_eff_1->SetLineWidth(2);
  g_eff_1->SetLineStyle(1);
  g_eff_1->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_1->SetPoint(i,(num_B/denom_B_1),(num_S/denom_S_1));
  }

  // Trimmed N-subjettiness
  TGraph *g_eff_2 = new TGraph(21);
  g_eff_2->SetName("g_eff_2");
  g_eff_2->SetLineColor(kRed);
  g_eff_2->SetLineWidth(2);
  g_eff_2->SetLineStyle(2);
  g_eff_2->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_MassDrop_S->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_MassDrop_B->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_2->SetPoint(i,(num_B/denom_B_2),(num_S/denom_S_2));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg",";Mistag rate; Tagging efficiency",100,fXmin,fXmax,100,0,0.9);
  bkg->SetTitleOffset(1.05,"Y");
  bkg->Draw();

  g_eff_1->Draw("L");
  g_eff_2->Draw("L");

  TLegend *legend = new TLegend(.45,.25,.75,.45);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_eff_1, fLeg1.c_str(),"l");
  legend->AddEntry(g_eff_2, fLeg2.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.14,0.96, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete g_eff_2;
  delete g_eff_1;
  delete file_S;
  delete file_B;
}


void efficiency_curves_comp_xrange(const string& fFileS1, const string& fFileS2, const string& fFileB1, const string& fFileB2,
                                   const string& fPlotS1, const string& fPlotS2, const string& fPlotB1, const string& fPlotB2,
                                   const double fXMin, const double fXMax, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                                   const string& fLeg1, const string& fLeg2, const double fXmin, const double fXmax, const double fYmin, const double fYmax,
                                   const string& fOutputFile, const Int_t fLogy=0, const double fOPL1=0.244, const double fOPM1=0.679, const double fOPT1=0.898,
                                   const double fOPL2=0.244, const double fOPM2=0.679, const double fOPT2=0.898)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetGridColor(kGray);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  // signal file
  TFile *file_S1  = new TFile(fFileS1.c_str());
  TFile *file_S2  = new TFile(fFileS2.c_str());

  // background file
  TFile *file_B1 = new TFile(fFileB1.c_str());
  TFile *file_B2 = new TFile(fFileB2.c_str());

  // signal histograms
  TH2 *h2_S_1 = (TH2*)file_S1->Get(fPlotS1.c_str());
  TH2 *h2_S_2 = (TH2*)file_S2->Get(fPlotS2.c_str());

  // background histograms
  TH2 *h2_B_1 = (TH2*)file_B1->Get(fPlotB1.c_str());
  TH2 *h2_B_2 = (TH2*)file_B2->Get(fPlotB2.c_str());

  // signal denominator counts
  double denom_S_1 = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),0,101);
  double denom_S_2 = h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),0,101);

  // background denominator counts
  double denom_B_1 = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),0,101);
  double denom_B_2 = h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),0,101);


  TGraph *g_eff_1 = new TGraph(29);
  g_eff_1->SetName("g_eff_1");
  g_eff_1->SetLineColor(kGreen+2);
  g_eff_1->SetLineWidth(2);
  g_eff_1->SetLineStyle(1);
  g_eff_1->SetMarkerStyle(20);

  for(int i = 0; i<5; ++i)
  {
    double num_S = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),101-i,101);
    double num_B = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),101-i,101);

    g_eff_1->SetPoint(i,(num_S/denom_S_1),(num_B/denom_B_1));
  }
  for(int i = 1; i<20; ++i)
  {
    double num_S = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),101-(i*5),101);
    double num_B = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),101-(i*5),101);

    g_eff_1->SetPoint(i+4,(num_S/denom_S_1),(num_B/denom_B_1));
  }
  for(int i = 1; i<6; ++i)
  {
    double num_S = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),6-i,101);
    double num_B = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),6-i,101);

    g_eff_1->SetPoint(i+23,(num_S/denom_S_1),(num_B/denom_B_1));
  }

  TGraph *g_eff_2 = new TGraph(29);
  g_eff_2->SetName("g_eff_2");
  g_eff_2->SetLineColor(kRed);
  g_eff_2->SetLineWidth(2);
  g_eff_2->SetLineStyle(2);
  g_eff_2->SetMarkerStyle(20);

  for(int i = 0; i<5; ++i)
  {
    double num_S = h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),101-i,101);
    double num_B = h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),101-i,101);

    g_eff_2->SetPoint(i,(num_S/denom_S_2),(num_B/denom_B_2));
  }
  for(int i = 1; i<20; ++i)
  {
    double num_S = h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),101-(i*5),101);
    double num_B = h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),101-(i*5),101);

    g_eff_2->SetPoint(i+4,(num_S/denom_S_2),(num_B/denom_B_2));
  }
  for(int i = 1; i<6; ++i)
  {
    double num_S = h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),6-i,101);
    double num_B = h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),6-i,101);

    g_eff_2->SetPoint(i+23,(num_S/denom_S_2),(num_B/denom_B_2));
  }

  // loose operating point
  TGraph *g_eff_L = new TGraph(2);
  g_eff_L->SetName("g_eff_L");
  g_eff_L->SetMarkerStyle(31);
  g_eff_L->SetMarkerSize(1.5);

  g_eff_L->SetPoint(0,(h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),h2_S_1->GetYaxis()->FindBin(fOPL1),101)/denom_S_1),(h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),h2_B_1->GetYaxis()->FindBin(fOPL1),101)/denom_B_1));
  g_eff_L->SetPoint(1,(h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),h2_S_2->GetYaxis()->FindBin(fOPL2),101)/denom_S_2),(h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),h2_B_2->GetYaxis()->FindBin(fOPL2),101)/denom_B_2));

  // medium operating point
  TGraph *g_eff_M = new TGraph(2);
  g_eff_M->SetName("g_eff_L");
  g_eff_M->SetMarkerStyle(27);
  g_eff_M->SetMarkerSize(1.5);

  g_eff_M->SetPoint(0,(h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),h2_S_1->GetYaxis()->FindBin(fOPM1),101)/denom_S_1),(h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),h2_B_1->GetYaxis()->FindBin(fOPM1),101)/denom_B_1));
  g_eff_M->SetPoint(1,(h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),h2_S_2->GetYaxis()->FindBin(fOPM2),101)/denom_S_2),(h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),h2_B_2->GetYaxis()->FindBin(fOPM2),101)/denom_B_2));

  // tight operating point
  TGraph *g_eff_T = new TGraph(2);
  g_eff_T->SetName("g_eff_L");
  g_eff_T->SetMarkerStyle(30);
  g_eff_T->SetMarkerSize(1.5);

  g_eff_T->SetPoint(0,(h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),h2_S_1->GetYaxis()->FindBin(fOPT1),101)/denom_S_1),(h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),h2_B_1->GetYaxis()->FindBin(fOPT1),101)/denom_B_1));
  g_eff_T->SetPoint(1,(h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),h2_S_2->GetYaxis()->FindBin(fOPT2),101)/denom_S_2),(h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),h2_B_2->GetYaxis()->FindBin(fOPT2),101)/denom_B_2));


  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(1.1,"Y");
  bkg->Draw();
  c->SetGridx();
  c->SetGridy();

  g_eff_1->Draw("L");
  g_eff_2->Draw("L");

  g_eff_L->Draw("P");
  g_eff_M->Draw("P");
  g_eff_T->Draw("P");

  TLegend *legend = new TLegend(.16,.64,.36,.77);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_eff_1, fLeg1.c_str(),"l");
  legend->AddEntry(g_eff_2, fLeg2.c_str(),"l");
  legend->Draw();

  TLegend *legend2 = new TLegend(.16,.45,.36,.60);
  legend2->SetBorderSize(0);
  legend2->SetFillColor(0);
  legend2->SetFillStyle(0);
  legend2->SetTextFont(42);
  legend2->SetTextSize(0.03);
  legend2->AddEntry(g_eff_L, "Loose","p");
  legend2->AddEntry(g_eff_M, "Medium","p");
  legend2->AddEntry(g_eff_T, "Tight","p");
  legend2->Draw();

  TLatex l1;
  l1.SetTextAlign(13);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.14+0.03,0.90, fTitle.c_str());

  l1.SetTextAlign(12);
  l1.SetTextSize(0.045);
  l1.SetTextFont(62);
  l1.DrawLatex(0.14,0.96, "CMS Simulation Preliminary, #sqrt{s} = 8 TeV");
  //l1.DrawLatex(0.14,0.97, "CMS Simulation");
  //l1.SetTextFont(42);
  //l1.DrawLatex(0.14+0.40,0.97, "#sqrt{s} = 8 TeV");
  l1.SetTextFont(42);
  l1.SetTextSize(0.04);
  if(fOutputFile.find("JTA")!=string::npos) l1.DrawLatex(0.48,0.18, "JTA = jet-track association");


  if(fLogy) c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete g_eff_2;
  delete g_eff_1;
  delete file_S1;
  delete file_S2;
  delete file_B1;
  delete file_B2;
}


void efficiency_vs_cut(const string& fInputFile, const string& fFileDir, const string& fVariable, const string& fPtRange,
                       const int fPVLow, const int fPVHigh, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                       const double fXmin, const double fXmax, const Double_t fYmin, const Double_t fYmax, const string& fOutputFile,
                       const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
                       const Double_t fLeftMargin=0.14, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gROOT->ForceStyle();

  // input files
  TFile *file  = new TFile(fInputFile.c_str());

  // histograms
  TH2D *h2_nPV_JetMass = (TH2D*)file->Get((fFileDir + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_variable = (TH2D*)file->Get((fFileDir + "/h2_nPV_" + fVariable + "_" + fPtRange).c_str());

  // denominator counts
  double denom = h2_nPV_JetMass->Integral(fPVLow,fPVHigh,0,401);

  // efficiency curve
  TGraph *g_eff = new TGraph(21);
  g_eff->SetName("g_eff");
  g_eff->SetLineColor(kGreen+2);
  g_eff->SetLineWidth(2);
  g_eff->SetLineStyle(1);
  g_eff->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num = h2_nPV_variable->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff->SetPoint(i,(h2_nPV_variable->GetYaxis()->GetBinUpEdge(101-(1+i*5))),(num/denom));
  }

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();

  g_eff->Draw("L");

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(fLeftMargin,0.96, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete bkg;
  delete c;
  delete g_eff;
  delete file;
}


void efficiency1D(const string& fInputFile, const string& fPlotPass, const string& fPlotTotal, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                  const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin, const Double_t fYmax,
                  const string& fOutputFile, const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
                  const Double_t fLeftMargin=0.12, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetGridColor(kGray);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gROOT->ForceStyle();

  // input files
  TFile *file  = new TFile(fInputFile.c_str());

  // histograms
  TH1D *h1_Pass = (TH1D*)file->Get(fPlotPass.c_str());

  TH1D *h1_Total = (TH1D*)file->Get(fPlotTotal.c_str());

  h1_Pass->Rebin(fRebinX);
  h1_Total->Rebin(fRebinX);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();
  c->SetGridx();
  c->SetGridy();

  TGraphAsymmErrors *g_efficiency = new TGraphAsymmErrors(h1_Pass, h1_Total,"cp");
  g_efficiency->SetLineWidth(2);
  g_efficiency->SetLineColor(kBlue+2);
  g_efficiency->SetMarkerSize(1.);
  g_efficiency->SetMarkerStyle(24);
  g_efficiency->SetMarkerColor(kBlue+2);

  g_efficiency->Draw("LP");

  TLatex l1;
  l1.SetTextAlign(13);
  l1.SetTextFont(42);
  l1.SetTextSize(0.045);
  l1.SetNDC();
  if( fOutputFile.find("Matching_eff_dR")!=string::npos ) l1.DrawLatex(fLeftMargin+0.03,0.21, fTitle.c_str());
  else                                                    l1.DrawLatex(fLeftMargin+0.03,0.90, fTitle.c_str());

  l1.SetTextAlign(12);
  l1.SetTextSize(0.05);
  l1.SetTextFont(62);
  l1.DrawLatex(fLeftMargin,0.97, "CMS Simulation Preliminary, #sqrt{s} = 8 TeV");
  //l1.DrawLatex(fLeftMargin,0.97, "CMS Simulation");
  //l1.SetTextFont(42);
  //l1.DrawLatex(fLeftMargin+0.35,0.97, "#sqrt{s} = 8 TeV");

  c->SetLogz();
  c->SaveAs(fOutputFile.c_str());

  delete g_efficiency;
  delete bkg;
  delete c;
  delete file;
}


void efficiency1D_overlay(const string& fInputFile, const string& fPlotPass, const string& fPlotTotal, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                          const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin, const Double_t fYmax,
                          const string& fOutputFile, const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
                          const Double_t fLeftMargin=0.12, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8, const string& fSubJetMode="k_{T}")
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetGridColor(kGray);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gROOT->ForceStyle();

  // input files
  TFile *file  = new TFile(fInputFile.c_str());

  // histograms
  TH1D *h1_Pass_CSVL = (TH1D*)file->Get((fPlotPass + "_CSVL").c_str());
  TH1D *h1_Pass_CSVM = (TH1D*)file->Get((fPlotPass + "_CSVM").c_str());
  TH1D *h1_Pass_SubJetCSVL = (TH1D*)file->Get((fPlotPass + "_SubJetCSVL").c_str());
  TH1D *h1_Pass_SubJetCSVM = (TH1D*)file->Get((fPlotPass + "_SubJetCSVM").c_str());
  TH1D *h1_Pass_DoubleB = (TH1D*)file->Get((fPlotPass + "_DoubleB").c_str());

  TH1D *h1_Total = (TH1D*)file->Get(fPlotTotal.c_str());

  h1_Pass_CSVL->Rebin(fRebinX);
  h1_Pass_CSVM->Rebin(fRebinX);
  h1_Pass_SubJetCSVL->Rebin(fRebinX);
  h1_Pass_SubJetCSVM->Rebin(fRebinX);
  h1_Pass_DoubleB->Rebin(fRebinX);
  h1_Total->Rebin(fRebinX);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();
  c->SetGridx();
  c->SetGridy();

  TGraphAsymmErrors *g_eff_CSVL = new TGraphAsymmErrors(h1_Pass_CSVL, h1_Total,"cp");
  g_eff_CSVL->SetLineWidth(2);
  g_eff_CSVL->SetLineColor(kBlue+2);
  g_eff_CSVL->SetMarkerSize(1.);
  g_eff_CSVL->SetMarkerStyle(20);
  g_eff_CSVL->SetMarkerColor(kBlue+2);

  TGraphAsymmErrors *g_eff_CSVM = new TGraphAsymmErrors(h1_Pass_CSVM, h1_Total,"cp");
  g_eff_CSVM->SetLineWidth(2);
  g_eff_CSVM->SetLineColor(kGreen+2);
  g_eff_CSVM->SetMarkerSize(1.);
  g_eff_CSVM->SetMarkerStyle(21);
  g_eff_CSVM->SetMarkerColor(kGreen+2);

  TGraphAsymmErrors *g_eff_SubJetCSVL = new TGraphAsymmErrors(h1_Pass_SubJetCSVL, h1_Total,"cp");
  g_eff_SubJetCSVL->SetLineWidth(2);
  g_eff_SubJetCSVL->SetLineStyle(2);
  g_eff_SubJetCSVL->SetLineColor(kBlue+2);
  g_eff_SubJetCSVL->SetMarkerSize(1.);
  g_eff_SubJetCSVL->SetMarkerStyle(24);
  g_eff_SubJetCSVL->SetMarkerColor(kBlue+2);

  TGraphAsymmErrors *g_eff_SubJetCSVM = new TGraphAsymmErrors(h1_Pass_SubJetCSVM, h1_Total,"cp");
  g_eff_SubJetCSVM->SetLineWidth(2);
  g_eff_SubJetCSVM->SetLineStyle(2);
  g_eff_SubJetCSVM->SetLineColor(kGreen+2);
  g_eff_SubJetCSVM->SetMarkerSize(1.);
  g_eff_SubJetCSVM->SetMarkerStyle(25);
  g_eff_SubJetCSVM->SetMarkerColor(kGreen+2);

  TGraphAsymmErrors *g_eff_DoubleB = new TGraphAsymmErrors(h1_Pass_DoubleB, h1_Total,"cp");
  g_eff_DoubleB->SetLineWidth(2);
  g_eff_DoubleB->SetLineColor(kRed+1);
  g_eff_DoubleB->SetMarkerSize(1.);
  g_eff_DoubleB->SetMarkerStyle(22);
  g_eff_DoubleB->SetMarkerColor(kRed+1);

  g_eff_CSVL->Draw("LP");
  g_eff_CSVM->Draw("LP");
  g_eff_SubJetCSVL->Draw("LP");
  g_eff_SubJetCSVM->Draw("LP");
  g_eff_DoubleB->Draw("LP");

  TLegend *legend = new TLegend(.15+(fLeftMargin-0.12),.67,.35+(fLeftMargin-0.12),.92);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->AddEntry(g_eff_CSVL, "CSVL","lp");
  legend->AddEntry(g_eff_CSVM, "CSVM","lp");
  legend->AddEntry(g_eff_SubJetCSVL, ("Subjet CSVL ("+ fSubJetMode +")").c_str(),"lp");
  legend->AddEntry(g_eff_SubJetCSVM, ("Subjet CSVM ("+ fSubJetMode +")").c_str(),"lp");
  legend->AddEntry(g_eff_DoubleB, "DoubleB","lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetTextSize(0.05);
  l1.SetNDC();
  l1.DrawLatex(fLeftMargin,0.97, fTitle.c_str());

  c->SetLogz();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete file;
}

void efficiency1D_overlayMulti_3(const string& fInputFile1, const string& fInputFile2, const string& fInputFile3,
                               const string& fPlotPass, const string& fPlotTotal, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                               const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin, const Double_t fYmax,
                               const string& fOutputFile, const Int_t fLogy=0, const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
                               const Double_t fLeftMargin=0.12, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetGridColor(kGray);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gROOT->ForceStyle();

  // input files
  TFile *file1  = new TFile(fInputFile1.c_str());
  TFile *file2  = new TFile(fInputFile2.c_str());
  TFile *file3  = new TFile(fInputFile3.c_str());

  // histograms
  TH1D *h1_Pass1 = (TH1D*)file1->Get(fPlotPass.c_str());
  TH1D *h1_Pass2 = (TH1D*)file2->Get(fPlotPass.c_str());
  TH1D *h1_Pass3 = (TH1D*)file3->Get(fPlotPass.c_str());

  TH1D *h1_Total1 = (TH1D*)file1->Get(fPlotTotal.c_str());
  TH1D *h1_Total2 = (TH1D*)file2->Get(fPlotTotal.c_str());
  TH1D *h1_Total3 = (TH1D*)file3->Get(fPlotTotal.c_str());

  h1_Pass1->Rebin(fRebinX);
  h1_Pass2->Rebin(fRebinX);
  h1_Pass3->Rebin(fRebinX);

  h1_Total1->Rebin(fRebinX);
  h1_Total2->Rebin(fRebinX);
  h1_Total3->Rebin(fRebinX);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();
  c->SetGridx();
  c->SetGridy();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();

  TGraphAsymmErrors *g_efficiency1 = new TGraphAsymmErrors(h1_Pass1, h1_Total1,"cp");
  g_efficiency1->SetLineWidth(2);
  g_efficiency1->SetLineStyle(1);
  g_efficiency1->SetLineColor(kGreen+2);
  g_efficiency1->SetMarkerSize(1.);
  g_efficiency1->SetMarkerStyle(24);
  g_efficiency1->SetMarkerColor(kGreen+2);

  TGraphAsymmErrors *g_efficiency2 = new TGraphAsymmErrors(h1_Pass2, h1_Total2,"cp");
  g_efficiency2->SetLineWidth(2);
  g_efficiency2->SetLineStyle(2);
  g_efficiency2->SetLineColor(kRed);
  g_efficiency2->SetMarkerSize(1.);
  g_efficiency2->SetMarkerStyle(26);
  g_efficiency2->SetMarkerColor(kRed);

  TGraphAsymmErrors *g_efficiency3 = new TGraphAsymmErrors(h1_Pass3, h1_Total3,"cp");
  g_efficiency3->SetLineWidth(2);
  g_efficiency3->SetLineStyle(3);
  g_efficiency3->SetLineColor(kBlue+2);
  g_efficiency3->SetMarkerSize(1.);
  g_efficiency3->SetMarkerStyle(20);
  g_efficiency3->SetMarkerColor(kBlue+2);

  g_efficiency1->Draw("LP");
  g_efficiency2->Draw("LPsame");
  g_efficiency3->Draw("LPsame");

  TLegend *legend = new TLegend(.15+(fLeftMargin-0.12),.67,.35+(fLeftMargin-0.12),.90);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.045);
  legend->AddEntry(g_efficiency3, "106<m_{jet}<135 GeV/c^{2} (pruned)","lp");
  legend->AddEntry(g_efficiency1, "75<m_{pruned}<135 GeV/c^{2}","lp");
  legend->AddEntry(g_efficiency2, "75<m_{jet}<106 GeV/c^{2} (pruned)","lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(13);
  l1.SetTextFont(42);
  l1.SetTextSize(0.05);
  l1.SetNDC();
  l1.DrawLatex(fLeftMargin+0.03,0.27, fTitle.c_str());

  l1.SetTextAlign(12);
  l1.SetTextSize(0.05);
  l1.SetTextFont(62);
  l1.DrawLatex(fLeftMargin,0.97, "CMS Simulation");
  l1.SetTextFont(42);
  l1.DrawLatex(fLeftMargin+0.35,0.97, "#sqrt{s} = 8 TeV");

  //c->RedrawAxis();
  c->SetLogz();
  if(fLogy) c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete g_efficiency1;
  delete g_efficiency2;
  delete bkg;
  delete c;
  delete file1;
  delete file2;
}


void efficiency1D_overlayMulti_6(const string& fInputFile1, const string& fInputFile2, const string& fInputFile3, const string& fInputFile4, const string& fInputFile5,
				 const string& fDir, const string& fPlotPass, const string& fPlotTotal, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
				 const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin, const Double_t fYmax,
				 const string& fOutputFile, const Int_t fLogy=0, const Int_t fRemoveLast=0, const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
				 const Double_t fLeftMargin=0.12, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();

  gStyle->SetGridColor(kGray);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gROOT->ForceStyle();

  // input files
  TFile *file1  = new TFile(fInputFile1.c_str());
  TFile *file2  = new TFile(fInputFile2.c_str());
  TFile *file3  = new TFile(fInputFile3.c_str());
  TFile *file4  = new TFile(fInputFile4.c_str());
  TFile *file5  = new TFile(fInputFile5.c_str());

  // histograms
  TH1D *h1_Pass1 = (TH1D*)file1->Get((fDir + "/" + fPlotPass).c_str());
  TH1D *h1_Pass2 = (TH1D*)file2->Get((fDir + "/" + fPlotPass).c_str());
  TH1D *h1_Pass3 = (TH1D*)file3->Get((fDir + "/" + fPlotPass).c_str());
  TH1D *h1_Pass4 = (TH1D*)file4->Get((fDir + "/" + fPlotPass).c_str());
  TH1D *h1_Pass5 = (TH1D*)file5->Get((fDir + "_bJetsGSP/" + fPlotPass).c_str());
  TH1D *h1_Pass6 = (TH1D*)file5->Get((fDir + "_udsgJets/" + fPlotPass).c_str());

  TH1D *h1_Total1 = (TH1D*)file1->Get((fDir + "/" + fPlotTotal).c_str());
  TH1D *h1_Total2 = (TH1D*)file2->Get((fDir + "/" + fPlotTotal).c_str());
  TH1D *h1_Total3 = (TH1D*)file3->Get((fDir + "/" + fPlotTotal).c_str());
  TH1D *h1_Total4 = (TH1D*)file4->Get((fDir + "/" + fPlotTotal).c_str());
  TH1D *h1_Total5 = (TH1D*)file5->Get((fDir + "_bJetsGSP/" + fPlotTotal).c_str());
  TH1D *h1_Total6 = (TH1D*)file5->Get((fDir + "_udsgJets/" + fPlotTotal).c_str());

  h1_Pass1->Rebin(fRebinX);
  h1_Pass2->Rebin(fRebinX);
  h1_Pass3->Rebin(fRebinX);
  h1_Pass4->Rebin(fRebinX);
  h1_Pass5->Rebin(fRebinX);
  h1_Pass6->Rebin(fRebinX);

  h1_Total1->Rebin(fRebinX);
  h1_Total2->Rebin(fRebinX);
  h1_Total3->Rebin(fRebinX);
  h1_Total4->Rebin(fRebinX);
  h1_Total5->Rebin(fRebinX);
  h1_Total6->Rebin(fRebinX);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();
  c->SetGridx();
  c->SetGridy();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();

  TGraphAsymmErrors *g_efficiency1 = new TGraphAsymmErrors(h1_Pass1, h1_Total1,"cp");
  g_efficiency1->SetLineWidth(2);
  g_efficiency1->SetLineStyle(1);
  g_efficiency1->SetLineColor(kGreen+2);
  g_efficiency1->SetMarkerSize(1.);
  g_efficiency1->SetMarkerStyle(24);
  g_efficiency1->SetMarkerColor(kGreen+2);

  TGraphAsymmErrors *g_efficiency2 = new TGraphAsymmErrors(h1_Pass2, h1_Total2,"cp");
  g_efficiency2->SetLineWidth(2);
  g_efficiency2->SetLineStyle(2);
  g_efficiency2->SetLineColor(kBlue+2);
  g_efficiency2->SetMarkerSize(1.);
  g_efficiency2->SetMarkerStyle(20);
  g_efficiency2->SetMarkerColor(kBlue+2);

  TGraphAsymmErrors *g_efficiency3 = new TGraphAsymmErrors(h1_Pass3, h1_Total3,"cp");
  g_efficiency3->SetLineWidth(2);
  g_efficiency3->SetLineStyle(3);
  g_efficiency3->SetLineColor(kOrange+1);
  g_efficiency3->SetMarkerSize(1.);
  g_efficiency3->SetMarkerStyle(22);
  g_efficiency3->SetMarkerColor(kOrange+1);

  TGraphAsymmErrors *g_efficiency4 = new TGraphAsymmErrors(h1_Pass4, h1_Total4,"cp");
  g_efficiency4->SetLineWidth(2);
  g_efficiency4->SetLineStyle(4);
  g_efficiency4->SetLineColor(kMagenta+2);
  g_efficiency4->SetMarkerSize(1.);
  g_efficiency4->SetMarkerStyle(23);
  g_efficiency4->SetMarkerColor(kMagenta+2);

  TGraphAsymmErrors *g_efficiency5 = new TGraphAsymmErrors(h1_Pass5, h1_Total5,"cp");
  g_efficiency5->SetLineWidth(2);
  g_efficiency5->SetLineStyle(5);
  g_efficiency5->SetLineColor(kRed);
  g_efficiency5->SetMarkerSize(1.);
  g_efficiency5->SetMarkerStyle(26);
  g_efficiency5->SetMarkerColor(kRed);

  TGraphAsymmErrors *g_efficiency6 = new TGraphAsymmErrors(h1_Pass6, h1_Total6,"cp");
  g_efficiency6->SetLineWidth(2);
  g_efficiency6->SetLineStyle(6);
  g_efficiency6->SetLineColor(kRed+3);
  g_efficiency6->SetMarkerSize(1.);
  g_efficiency6->SetMarkerStyle(27);
  g_efficiency6->SetMarkerColor(kRed+3);

  if(fRemoveLast)
  {
    // remove last point
    g_efficiency1->RemovePoint(g_efficiency1->GetN()-1);
    g_efficiency2->RemovePoint(g_efficiency2->GetN()-1);
    g_efficiency3->RemovePoint(g_efficiency3->GetN()-1);
    g_efficiency4->RemovePoint(g_efficiency4->GetN()-1);
    g_efficiency5->RemovePoint(g_efficiency5->GetN()-1);
    g_efficiency6->RemovePoint(g_efficiency6->GetN()-1);
  }

  g_efficiency1->Draw("LP");
  g_efficiency2->Draw("LPsame");
  g_efficiency3->Draw("LPsame");
  g_efficiency4->Draw("LPsame");
  g_efficiency5->Draw("LPsame");
  g_efficiency6->Draw("LPsame");

  TLegend *legend = new TLegend(.15+(fLeftMargin-0.12),.55,.35+(fLeftMargin-0.12),.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_efficiency1, "H(120)#rightarrowb#bar{b}","lp");
  legend->AddEntry(g_efficiency5, "g#rightarrowb#bar{b} splitting","lp");
  legend->AddEntry(g_efficiency3, "Hadronic Z","lp");
  legend->AddEntry(g_efficiency4, "Hadronic top","lp");
  legend->AddEntry(g_efficiency2, "Hadronic W","lp");
  legend->AddEntry(g_efficiency6, "udsg jets","lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(13);
  l1.SetTextFont(42);
  l1.SetTextSize(0.045);
  l1.SetNDC();
  l1.DrawLatex(fLeftMargin+0.03,0.26, fTitle.c_str());

  //l1.SetTextAlign(12);
  //l1.SetTextSize(0.05);
  //l1.SetTextFont(62);
  //l1.DrawLatex(fLeftMargin,0.97, "CMS Simulation Preliminary, #sqrt{s} = 8 TeV");
  //l1.DrawLatex(fLeftMargin,0.97, "CMS Simulation");
  //l1.SetTextFont(42);
  //l1.DrawLatex(fLeftMargin+0.35,0.97, "#sqrt{s} = 8 TeV");
  
  int iPeriod = 2;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV 

  // second parameter is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // more generally : 
  // iPos = 10*(alignment 1/2/3) + position (1/2/3 = left/center/right)

  int iPos = 0; // out of frame (in exceptional cases)

  // New Pub Comm recommendation
  CMS_lumi( c, iPeriod, iPos );

  //c->RedrawAxis();
  c->SetLogz();
  if(fLogy) c->SetLogy();

  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete g_efficiency1;
  delete g_efficiency2;
  delete bkg;
  delete c;
  delete file1;
  delete file2;
}


void makePlots()
{
  writeExtraText = true;       // if extra text
  extraText  = "Simulation Preliminary";  // default extra text is "Preliminary"
  relPosX    = 0.12;
  lumi_8TeV  = ""; // default is "19.7 fb^{-1}"

  //--------------------------------------------------------------------------------------------------------------------
  // overlay multiple backgrounds
  // Subjet IVFCSVL
  efficiency1D_overlayMulti_6("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeTprimeToBWBWinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/BprimeBprimeToBZBZinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeToTHinc_M-1700_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "jetAnalyzerFatJets_PrunedSubjets", "h1_JetPt_BosonMatched_JetMass_SubJetMinIVFCSVL", "h1_JetPt_BosonMatched_JetMass",
			      "#splitline{AK R=0.8}{75<m_{pruned}<135 GeV/c^{2}, Subjet IVFCSVL}",
			      "Fat jet p_{T} [GeV/c]", "b-tagging efficiency", 25, 0, 1000, 0.001, 1,
			      "b-tag_eff_vs_FatJetPt_SubjetIVFCSVL_PrunedJetMass.eps", 1);

  efficiency1D_overlayMulti_6("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeTprimeToBWBWinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/BprimeBprimeToBZBZinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeToTHinc_M-1700_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "jetAnalyzerFatJets_PrunedSubjets", "h1_nPV_BosonMatched_JetMass_SubJetMinIVFCSVL", "h1_nPV_BosonMatched_JetMass",
			      "#splitline{AK R=0.8}{75<m_{pruned}<135 GeV/c^{2}, Subjet IVFCSVL}",
			      "PV multiplicity", "b-tagging efficiency", 10, -0.5, 39.5, 0.001, 1,
			      "b-tag_eff_vs_nPV_SubjetIVFCSVL_PrunedJetMass.eps", 1, 1);

  // Subjet IVFCSVM
  efficiency1D_overlayMulti_6("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeTprimeToBWBWinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/BprimeBprimeToBZBZinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeToTHinc_M-1700_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "jetAnalyzerFatJets_PrunedSubjets", "h1_JetPt_BosonMatched_JetMass_SubJetMinIVFCSVM", "h1_JetPt_BosonMatched_JetMass",
			      "#splitline{AK R=0.8}{75<m_{pruned}<135 GeV/c^{2}, Subjet IVFCSVM}",
			      "Fat jet p_{T} [GeV/c]", "b-tagging efficiency", 25, 0, 1000, 0.0001, 1,
			      "b-tag_eff_vs_FatJetPt_SubjetIVFCSVM_PrunedJetMass.eps", 1);

  efficiency1D_overlayMulti_6("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeTprimeToBWBWinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/BprimeBprimeToBZBZinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeToTHinc_M-1700_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "jetAnalyzerFatJets_PrunedSubjets", "h1_nPV_BosonMatched_JetMass_SubJetMinIVFCSVM", "h1_nPV_BosonMatched_JetMass",
			      "#splitline{AK R=0.8}{75<m_{pruned}<135 GeV/c^{2}, Subjet IVFCSVM}",
			      "PV multiplicity", "b-tagging efficiency", 10, -0.5, 39.5, 0.0001, 1,
			      "b-tag_eff_vs_nPV_SubjetIVFCSVM_PrunedJetMass.eps", 1, 1);

  // Fat jet IVFCSVL
  efficiency1D_overlayMulti_6("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeTprimeToBWBWinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/BprimeBprimeToBZBZinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeToTHinc_M-1700_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "jetAnalyzerFatJets_PrunedSubjets", "h1_JetPt_BosonMatched_JetMass_IVFCSVL", "h1_JetPt_BosonMatched_JetMass",
			      "#splitline{AK R=0.8}{75<m_{pruned}<135 GeV/c^{2}, Fat jet IVFCSVL}",
			      "Fat jet p_{T} [GeV/c]", "b-tagging efficiency", 25, 0, 1000, 0.1, 1,
			      "b-tag_eff_vs_FatJetPt_FatJetIVFCSVL_PrunedJetMass.eps", 1);

  efficiency1D_overlayMulti_6("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeTprimeToBWBWinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/BprimeBprimeToBZBZinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeToTHinc_M-1700_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "jetAnalyzerFatJets_PrunedSubjets", "h1_nPV_BosonMatched_JetMass_IVFCSVL", "h1_nPV_BosonMatched_JetMass",
			      "#splitline{AK R=0.8}{75<m_{pruned}<135 GeV/c^{2}, Fat jet IVFCSVL}",
			      "PV multiplicity", "b-tagging efficiency", 10, -0.5, 39.5, 0.1, 1,
			      "b-tag_eff_vs_nPV_FatJetIVFCSVL_PrunedJetMass.eps", 1, 1);

  // Fat jet IVFCSVM
  efficiency1D_overlayMulti_6("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeTprimeToBWBWinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/BprimeBprimeToBZBZinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeToTHinc_M-1700_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "jetAnalyzerFatJets_PrunedSubjets", "h1_JetPt_BosonMatched_JetMass_IVFCSVM", "h1_JetPt_BosonMatched_JetMass",
			      "#splitline{AK R=0.8}{75<m_{pruned}<135 GeV/c^{2}, Fat jet IVFCSVM}",
			      "Fat jet p_{T} [GeV/c]", "b-tagging efficiency", 25, 0, 1000, 0.01, 1,
			      "b-tag_eff_vs_FatJetPt_FatJetIVFCSVM_PrunedJetMass.eps", 1);

  efficiency1D_overlayMulti_6("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeTprimeToBWBWinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/BprimeBprimeToBZBZinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/TprimeToTHinc_M-1700_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",
			      "jetAnalyzerFatJets_PrunedSubjets", "h1_nPV_BosonMatched_JetMass_IVFCSVM", "h1_nPV_BosonMatched_JetMass",
			      "#splitline{AK R=0.8}{75<m_{pruned}<135 GeV/c^{2}, Fat jet IVFCSVM}",
			      "PV multiplicity", "b-tagging efficiency", 10, -0.5, 39.5, 0.01, 1,
			      "b-tag_eff_vs_nPV_FatJetIVFCSVM_PrunedJetMass.eps", 1, 1);

  //--------------------------------------------------------------------------------------------------------------------
}
