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
#include "TPolyLine.h"
#include "exoStyle.C"
#include "CMS_lumi.C"


using namespace std;


void efficiency1D_overlayMulti_5(const string& fInputFile1, const string& fInputFile2, const string& fInputFile3, const string& fInputFile4, const string& fInputFile5,
				 const string& fDir, const string& fPlotPass, const string& fPlotTotal, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
				 const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin, const Double_t fYmax,
				 const string& fOutputFile, const Int_t fLogy=0, const Int_t fRemoveLast=0, const Int_t fLegendRight =0, 
				 const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0, const Double_t fLeftMargin=0.12, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
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
  TH1D *h1_Pass5 = (TH1D*)file5->Get((fDir + "/" + fPlotPass).c_str());

  TH1D *h1_Total1 = (TH1D*)file1->Get((fDir + "/" + fPlotTotal).c_str());
  TH1D *h1_Total2 = (TH1D*)file2->Get((fDir + "/" + fPlotTotal).c_str());
  TH1D *h1_Total3 = (TH1D*)file3->Get((fDir + "/" + fPlotTotal).c_str());
  TH1D *h1_Total4 = (TH1D*)file4->Get((fDir + "/" + fPlotTotal).c_str());
  TH1D *h1_Total5 = (TH1D*)file5->Get((fDir + "/" + fPlotTotal).c_str());

  h1_Pass1->Rebin(fRebinX);
  h1_Pass2->Rebin(fRebinX);
  h1_Pass3->Rebin(fRebinX);
  h1_Pass4->Rebin(fRebinX);
  h1_Pass5->Rebin(fRebinX);

  h1_Total1->Rebin(fRebinX);
  h1_Total2->Rebin(fRebinX);
  h1_Total3->Rebin(fRebinX);
  h1_Total4->Rebin(fRebinX);
  h1_Total5->Rebin(fRebinX);

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

  if(fRemoveLast)
  {
    Double_t x[5] = {fXmin,-0.5,-0.5,fXmin,fXmin};
    Double_t y[5] = {fYmin,fYmin,fYmax,fYmax,fYmin};
    TPolyLine *pline = new TPolyLine(5,x,y);
    pline->SetFillColor(kGray);
    //pline->SetFillStyle(3004);
    pline->Draw("f");
    c->RedrawAxis();
  }

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

  if(fRemoveLast)
  {
    // remove last point
    g_efficiency1->RemovePoint(g_efficiency1->GetN()-1);
    g_efficiency2->RemovePoint(g_efficiency2->GetN()-1);
    g_efficiency3->RemovePoint(g_efficiency3->GetN()-1);
    g_efficiency4->RemovePoint(g_efficiency4->GetN()-1);
    g_efficiency5->RemovePoint(g_efficiency5->GetN()-1);
  }

  g_efficiency1->Draw("LP");
  g_efficiency2->Draw("LPsame");
  g_efficiency3->Draw("LPsame");
  g_efficiency4->Draw("LPsame");
  g_efficiency5->Draw("LPsame");

  Double_t lx1 = .15+(fLeftMargin-0.12);
  Double_t lx2 = .35+(fLeftMargin-0.12);
  Double_t ly1 = .55;
  Double_t ly2 = .85;

  if(fLegendRight)
    {  
      lx1 = 1-.35+(fLeftMargin-0.1);
      lx2 = 1-.15+(fLeftMargin-0.1);
      ly1 = .55;
      ly2 = .85;
    }

  TLegend *legend = new TLegend(lx1,ly1,lx2,ly2);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_efficiency1, "H#rightarrowb#bar{b}","lp");
  legend->AddEntry(g_efficiency3, "Hadronic Z","lp");
  legend->AddEntry(g_efficiency4, "Hadronic top","lp");
  legend->AddEntry(g_efficiency2, "Hadronic W","lp");
  legend->AddEntry(g_efficiency5, "Inclusive QCD","lp");
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
  
  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

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
  lumi_sqrtS = "8 TeV"; // default is "19.7 fb^{-1}"

  //--------------------------------------------------------------------------------------------------------------------
  // overlay multiple backgrounds
  // Subjet CSVL -- double-b tagging
  efficiency1D_overlayMulti_5("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root",
			      "TprimeTprimeToBWBWinc_M-500_HiggsTagging_WBkg_R1p5_BTV-13-001_PATTuple_v3.root",
			      "TprimeTprimeToTZTZinc_M-500_HiggsTagging_ZBkg_R1p5_BTV-13-001_PATTuple_v3.root",
			      "TprimeTprimeToTHTHinc_M-600_HiggsTagging_TopBkg_R1p5_BTV-13-001_PATTuple_v3.root",
			      "QCDPythia6_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root",
			      "jetAnalyzerFatJets_FilteredSubjets", "h1_JetPt_BosonMatched_JetMass_SubJetMinCSVM", "h1_JetPt_BosonMatched_JetMass",
			      "#splitline{CA R=1.5, |#eta|<2.4}{m>60 GeV/c^{2}, Subjet CSVM}",
			      "Fat jet p_{T} [GeV/c]", "Double-b tagging efficiency", 25, 0, 500, 1E-3, 1,
			      "Double-b-tag_eff_vs_FatJetPt_SubjetCSVM.eps", 1);

  // Subjet CSVL -- Higgs tagging
  efficiency1D_overlayMulti_5("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root",
			      "TprimeTprimeToBWBWinc_M-500_HiggsTagging_WBkg_R1p5_BTV-13-001_PATTuple_v3.root",
			      "TprimeTprimeToTZTZinc_M-500_HiggsTagging_ZBkg_R1p5_BTV-13-001_PATTuple_v3.root",
			      "TprimeTprimeToTHTHinc_M-600_HiggsTagging_TopBkg_R1p5_BTV-13-001_PATTuple_v3.root",
			      "QCDPythia6_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root",
			      "jetAnalyzerFatJets_FilteredSubjets", "h1_JetPt_BosonMatched_JetMass_SubJetMinCSVM", "h1_JetPt_BosonMatched",
			      "#splitline{CA R=1.5, |#eta|<2.4}{m>60 GeV/c^{2}, Subjet CSVM}",
			      "Fat jet p_{T} [GeV/c]", "Higgs tagging efficiency", 25, 0, 500, 1E-3, 1,
			      "Higgs-tag_eff_vs_FatJetPt_SubjetCSVM.eps", 1);

  //--------------------------------------------------------------------------------------------------------------------
}
