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
#include <map>
#include "Rtypes.h"

using namespace std;


TGraph* getEfficiencyCurve(const string& fFileS1, const string& fFileB1,const string& fPlot,const double fXMin, const double fXMax)
{
  //get files and histograms
  TFile *file_S1  = new TFile(fFileS1.c_str());
  TFile *file_B1 = new TFile(fFileB1.c_str());

  TH2 *h2_S_1 = (TH2*)file_S1->Get(fPlot.c_str());
  TH2 *h2_B_1 = (TH2*)file_B1->Get(fPlot.c_str());

  //total jet count for denominator of efficiency calculation
  double denom_S_1 = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),0,101);
  double denom_B_1 = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),0,101);

  TGraph *g_eff_1 = new TGraph(29);

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
  return g_eff_1;
}

void formatGraph(TGraph* graph, int graphNum)
{
  short colors[7]={ kGreen+2, kRed, kBlue, kBlack, kMagenta, kOrange+2, kCyan };
  short graphColor = colors[graphNum % 7];
  int lineStyle = (graphNum % 11) + 1;
  graph->SetLineColor(graphColor);
  graph->SetLineStyle(lineStyle);
  graph->SetLineWidth(2);
}

std::string getHistName(std::string grooming = "Pruned",std::string algo = "CSV", bool subjetPlot = false)
{
  if (!subjetPlot)
    return "jetAnalyzerCA8FatJets_"+grooming+"Subjets/h2_JetPt_Jet"+algo+"_BosonMatched_JetMass";
  else
    return "jetAnalyzerCA8FatJets_"+grooming+"Subjets/h2_JetPt_SubJetMin"+algo+"_BosonMatched_JetMass"; 
}

void plotEfficiencyCurves(std::map< std::string,TGraph* > &graphs, const std::vector< std::string> &ordering, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                          const string& fExtraInfo, const string& fOutputFile, const double fXmin, const double fXmax, const double fYmin, const double fYmax,const Int_t fLogy=0)
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

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();
  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(1.1,"Y");
  bkg->Draw();
  c->SetGridx();
  c->SetGridy();

  TLegend *legend = new TLegend(.16,.64,.36,.77);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.021);

  int graphCounter = 0;
  for (std::vector<std::string>::const_iterator it = ordering.begin(); it != ordering.end(); ++it)
  {
    std::string label = (*it);
    TGraph* graph = graphs[label];
    legend->AddEntry(graph, label.c_str(),"l");
    formatGraph(graph,graphCounter);
    graph->Draw("L");
    graphCounter++;
  }

  if (fLogy) c->SetLogy();
  legend->Draw();
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

  l1.SetTextFont(42);
  l1.SetTextSize(0.025);
  l1.DrawLatex(0.48,0.18, fExtraInfo.c_str());
  c->SaveAs(fOutputFile.c_str());
  graphs.clear();
  delete c;
  delete legend;
  delete bkg;
}

void makePlots()
{
  // for multiple plots on the same canvas

  // vectors storing the order of legend entries
  std::vector< std::string> orderingPt300To500;
  std::vector< std::string> orderingPt700ToInf;
  // maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt300To500;
  std::map< std::string,TGraph* > graphsPt700ToInf;

  //==========================================
  // BTV-13-001 setup
  //==========================================
  
  graphsPt300To500["Fat Jet CSV (BTV-13-001)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_BTV-13-001.root","ROOT_files/QCDPythia6_HiggsTagging_BTV-13-001.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (BTV-13-001)"]  = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_BTV-13-001.root","ROOT_files/QCDPythia6_HiggsTagging_BTV-13-001.root",getHistName("Pruned","CSVL",true),300,500);
  
  orderingPt300To500.push_back("Fat Jet CSV (BTV-13-001)");
  orderingPt300To500.push_back("Subjet CSV (BTV-13-001)");
  //------------------------------------------
  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_BTV-13-001.root","ROOT_files/QCDPythia6_HiggsTagging_BTV-13-001.root",getHistName("Pruned","CSVL",false),700,1100);
  graphsPt700ToInf["Subjet CSV (BTV-13-001)"]  = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_BTV-13-001.root","ROOT_files/QCDPythia6_HiggsTagging_BTV-13-001.root",getHistName("Pruned","CSVL",true),700,1100);
  
  orderingPt700ToInf.push_back("Fat Jet CSV (BTV-13-001)");
  orderingPt700ToInf.push_back("Subjet CSV (BTV-13-001)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt300to500.eps",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt700toInf.eps",0, 1, 1E-3, 1,1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //==========================================
  // Post-BTV-13-001 setup
  //==========================================
  // Fat jets
  //------------------------------------------
  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]      = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_BTV-13-001.root","ROOT_files/QCDPythia6_HiggsTagging_BTV-13-001.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Fat Jet CSV"]                   = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","CSV",false),300,500);
  graphsPt300To500["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","IVFCSV",false),300,500);
  graphsPt300To500["Fat Jet JP"]                    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","JP",false),300,500);
  graphsPt300To500["Fat Jet JBP"]                   = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","JBP",false),300,500);
  graphsPt300To500["Fat Jet CSV (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","CSV",false),300,500);
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","IVFCSV",false),300,500);
  graphsPt300To500["Fat Jet JP (Explicit JTA)"]     = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","JP",false),300,500);
  graphsPt300To500["Fat Jet JBP (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","JBP",false),300,500);

  orderingPt300To500.push_back("Fat Jet CSV (BTV-13-001)");
  orderingPt300To500.push_back("Fat Jet CSV");
  //orderingPt300To500.push_back("Fat Jet CSV (Explicit JTA)");
  orderingPt300To500.push_back("Fat Jet IVFCSV");
  //orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt300To500.push_back("Fat Jet JP");
  //orderingPt300To500.push_back("Fat Jet JP (Explicit JTA)");
  orderingPt300To500.push_back("Fat Jet JBP");
  //orderingPt300To500.push_back("Fat Jet JBP (Explicit JTA)");
  //------------------------------------------
  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"]      = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_BTV-13-001.root","ROOT_files/QCDPythia6_HiggsTagging_BTV-13-001.root",getHistName("Pruned","CSVL",false),700,1100);
  graphsPt700ToInf["Fat Jet CSV"]                   = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","CSV",false),700,1100);
  graphsPt700ToInf["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","IVFCSV",false),700,1100);
  graphsPt700ToInf["Fat Jet JP"]                    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","JP",false),700,1100);
  graphsPt700ToInf["Fat Jet JBP"]                   = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","JBP",false),700,1100);
  graphsPt700ToInf["Fat Jet CSV (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","CSV",false),700,1100);
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","IVFCSV",false),700,1100);
  graphsPt700ToInf["Fat Jet JP (Explicit JTA)"]     = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","JP",false),700,1100);
  graphsPt700ToInf["Fat Jet JBP (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","JBP",false),700,1100);

  orderingPt700ToInf.push_back("Fat Jet CSV (BTV-13-001)");
  orderingPt700ToInf.push_back("Fat Jet CSV");
  //orderingPt700ToInf.push_back("Fat Jet CSV (Explicit JTA)");
  orderingPt700ToInf.push_back("Fat Jet IVFCSV");
  //orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt700ToInf.push_back("Fat Jet JP");
  //orderingPt700ToInf.push_back("Fat Jet JP (Explicit JTA)");
  orderingPt700ToInf.push_back("Fat Jet JBP");
  //orderingPt700ToInf.push_back("Fat Jet JBP (Explicit JTA)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt300to500_FatJets_comparison.eps",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt700toInf_FatJets_comparison.eps",0, 1, 1E-3, 1,1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //------------------------------------------
  // Fat jets: IVFCSV (Cone vs Explicit JTA)
  //------------------------------------------
  graphsPt300To500["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","IVFCSV",false),300,500);
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","IVFCSV",false),300,500);
  
  orderingPt300To500.push_back("Fat Jet IVFCSV");
  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA)");
  //------------------------------------------
  graphsPt700ToInf["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","IVFCSV",false),700,1100);
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","IVFCSV",false),700,1100);
  
  orderingPt700ToInf.push_back("Fat Jet IVFCSV");
  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt300to500_FatJets_IVFCSV_Cone_vs_Expl_JTA.eps",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt700toInf_FatJets_IVFCSV_Cone_vs_Expl_JTA.eps",0, 1, 1E-3, 1,1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //------------------------------------------
  // Subjets: CSV
  //------------------------------------------
  graphsPt300To500["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]           = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MD+Filtered)"]        = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (k_{T})"]              = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (k_{T}+Filtered)"]     = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("KtBDRSFiltered","CSV",true),300,500);

  orderingPt300To500.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt300To500.push_back("Subjet CSV (Filtered)");
  orderingPt300To500.push_back("Subjet CSV (MD+Filtered)");
  orderingPt300To500.push_back("Subjet CSV (k_{T})");
  orderingPt300To500.push_back("Subjet CSV (k_{T}+Filtered)");
  //------------------------------------------
  graphsPt700ToInf["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (Filtered)"]           = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Filtered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (MD+Filtered)"]        = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("MDBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (k_{T})"]              = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Kt","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (k_{T}+Filtered)"]     = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("KtBDRSFiltered","CSV",true),700,1100);

  orderingPt700ToInf.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt700ToInf.push_back("Subjet CSV (Filtered)");
  orderingPt700ToInf.push_back("Subjet CSV (MD+Filtered)");
  orderingPt700ToInf.push_back("Subjet CSV (k_{T})");
  orderingPt700ToInf.push_back("Subjet CSV (k_{T}+Filtered)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt300to500_Subjets_CSV.eps",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt700toInf_Subjets_CSV.eps",0, 1, 1E-3, 1,1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //------------------------------------------
  // Subjets: IVFCSV
  //------------------------------------------
  graphsPt300To500["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),300,500);
  graphsPt300To500["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Filtered","IVFCSV",true),300,500);
  graphsPt300To500["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("MDBDRSFiltered","IVFCSV",true),300,500);
  graphsPt300To500["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Kt","IVFCSV",true),300,500);
  graphsPt300To500["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("KtBDRSFiltered","IVFCSV",true),300,500);

  orderingPt300To500.push_back("Subjet IVFCSV (Pruned)");
  orderingPt300To500.push_back("Subjet IVFCSV (Filtered)");
  orderingPt300To500.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt300To500.push_back("Subjet IVFCSV (k_{T})");
  orderingPt300To500.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
  graphsPt700ToInf["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Filtered","IVFCSV",true),700,1100);
  graphsPt700ToInf["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("MDBDRSFiltered","IVFCSV",true),700,1100);
  graphsPt700ToInf["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Kt","IVFCSV",true),700,1100);
  graphsPt700ToInf["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("KtBDRSFiltered","IVFCSV",true),700,1100);

  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Filtered)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T})");
  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt300to500_Subjets_IVFCSV.eps",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt700toInf_Subjets_IVFCSV.eps",0, 1, 1E-3, 1,1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //------------------------------------------
  // Subjets: JP
  //------------------------------------------
  graphsPt300To500["Subjet JP (Pruned)"]         = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","JP",true),300,500);
  graphsPt300To500["Subjet JP (Filtered)"]       = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Filtered","JP",true),300,500);
  graphsPt300To500["Subjet JP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("MDBDRSFiltered","JP",true),300,500);
  graphsPt300To500["Subjet JP (k_{T})"]          = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Kt","JP",true),300,500);
  graphsPt300To500["Subjet JP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("KtBDRSFiltered","JP",true),300,500);

  orderingPt300To500.push_back("Subjet JP (Pruned)");
  orderingPt300To500.push_back("Subjet JP (Filtered)");
  orderingPt300To500.push_back("Subjet JP (MD+Filtered)");
  orderingPt300To500.push_back("Subjet JP (k_{T})");
  orderingPt300To500.push_back("Subjet JP (k_{T}+Filtered)");
  //------------------------------------------
  graphsPt700ToInf["Subjet JP (Pruned)"]         = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","JP",true),700,1100);
  graphsPt700ToInf["Subjet JP (Filtered)"]       = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Filtered","JP",true),700,1100);
  graphsPt700ToInf["Subjet JP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("MDBDRSFiltered","JP",true),700,1100);
  graphsPt700ToInf["Subjet JP (k_{T})"]          = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Kt","JP",true),700,1100);
  graphsPt700ToInf["Subjet JP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("KtBDRSFiltered","JP",true),700,1100);

  orderingPt700ToInf.push_back("Subjet JP (Pruned)");
  orderingPt700ToInf.push_back("Subjet JP (Filtered)");
  orderingPt700ToInf.push_back("Subjet JP (MD+Filtered)");
  orderingPt700ToInf.push_back("Subjet JP (k_{T})");
  orderingPt700ToInf.push_back("Subjet JP (k_{T}+Filtered)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt300to500_Subjets_JP.eps",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt700toInf_Subjets_JP.eps",0, 1, 1E-3, 1,1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //------------------------------------------
  // Subjets: JBP
  //------------------------------------------
  graphsPt300To500["Subjet JBP (Pruned)"]         = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","JBP",true),300,500);
  graphsPt300To500["Subjet JBP (Filtered)"]       = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Filtered","JBP",true),300,500);
  graphsPt300To500["Subjet JBP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("MDBDRSFiltered","JBP",true),300,500);
  graphsPt300To500["Subjet JBP (k_{T})"]          = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Kt","JBP",true),300,500);
  graphsPt300To500["Subjet JBP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("KtBDRSFiltered","JBP",true),300,500);

  orderingPt300To500.push_back("Subjet JBP (Pruned)");
  orderingPt300To500.push_back("Subjet JBP (Filtered)");
  orderingPt300To500.push_back("Subjet JBP (MD+Filtered)");
  orderingPt300To500.push_back("Subjet JBP (k_{T})");
  orderingPt300To500.push_back("Subjet JBP (k_{T}+Filtered)");
  //------------------------------------------
  graphsPt700ToInf["Subjet JBP (Pruned)"]         = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","JBP",true),700,1100);
  graphsPt700ToInf["Subjet JBP (Filtered)"]       = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Filtered","JBP",true),700,1100);
  graphsPt700ToInf["Subjet JBP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("MDBDRSFiltered","JBP",true),700,1100);
  graphsPt700ToInf["Subjet JBP (k_{T})"]          = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Kt","JBP",true),700,1100);
  graphsPt700ToInf["Subjet JBP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("KtBDRSFiltered","JBP",true),700,1100);

  orderingPt700ToInf.push_back("Subjet JBP (Pruned)");
  orderingPt700ToInf.push_back("Subjet JBP (Filtered)");
  orderingPt700ToInf.push_back("Subjet JBP (MD+Filtered)");
  orderingPt700ToInf.push_back("Subjet JBP (k_{T})");
  orderingPt700ToInf.push_back("Subjet JBP (k_{T}+Filtered)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt300to500_Subjets_JBP.eps",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt700toInf_Subjets_JBP.eps",0, 1, 1E-3, 1,1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //------------------------------------------
  // Subjets: Pruned comparison
  //------------------------------------------
  graphsPt300To500["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned)"]          = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]              = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","JP",true),300,500);
  graphsPt300To500["Subjet JBP (Pruned)"]             = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","JBP",true),300,500);
  
  orderingPt300To500.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned)");
  orderingPt300To500.push_back("Subjet JP (Pruned)");
  orderingPt300To500.push_back("Subjet JBP (Pruned)");
  //------------------------------------------
  graphsPt700ToInf["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","CSV",true),700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned)"]          = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),700,1100);
  graphsPt700ToInf["Subjet JP (Pruned)"]              = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","JP",true),700,1100);
  graphsPt700ToInf["Subjet JBP (Pruned)"]             = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","JBP",true),700,1100);
  
  orderingPt700ToInf.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned)");
  orderingPt700ToInf.push_back("Subjet JP (Pruned)");
  orderingPt700ToInf.push_back("Subjet JBP (Pruned)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt300to500_Subjets_Pruned_comparison.eps",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt700toInf_Subjets_Pruned_comparison.eps",0, 1, 1E-3, 1,1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //------------------------------------------
  // Subjets: Pruned IVFCSV comparison
  //------------------------------------------
  graphsPt300To500["Subjet IVFCSV (Pruned)"]                                           = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA)"]                             = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, SV Clustering)"]                            = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_SVClustering_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]              = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, Ext. PF CHS)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_ExtPFchs.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),300,500);

  orderingPt300To500.push_back("Subjet IVFCSV (Pruned)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, Ext. PF CHS)");

  //------------------------------------------
  graphsPt700ToInf["Subjet IVFCSV (Pruned)"]                                           = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA)"]                             = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, SV Clustering)"]                            = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_SVClustering_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]              = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, Ext. PF CHS)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_ExtPFchs.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),700,1100);

  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, Ext. PF CHS)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Extended PF CHS only for signal","btagperfcomp_Pt300to500_Subjets_Pruned_IVFCSV_comparison.eps",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Extended PF CHS only for signal","btagperfcomp_Pt700toInf_Subjets_Pruned_IVFCSV_comparison.eps",0, 1, 1E-3, 1,1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //------------------------------------------
  // Fat jets and subjets
  //------------------------------------------
  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]                                    = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_BTV-13-001.root","ROOT_files/QCDPythia6_HiggsTagging_BTV-13-001.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA)"]                               = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","IVFCSV",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned, BTV-13-001)"]                             = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]         = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),300,500);
  graphsPt300To500["Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",getHistName("Pruned","HybridIVFCSV",false),300,500);
  
  orderingPt300To500.push_back("Fat Jet CSV (BTV-13-001)");
  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt300To500.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt300To500.push_back("Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)");
  //------------------------------------------
  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"]                            = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_BTV-13-001.root","ROOT_files/QCDPythia6_HiggsTagging_BTV-13-001.root",getHistName("Pruned","CSVL",false),700,1100);
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root",getHistName("Pruned","IVFCSV",false),700,1100);
  graphsPt700ToInf["Subjet CSV (Pruned, BTV-13-001)"]                     = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_PATTuple_v3.root",getHistName("Pruned","CSV",true),700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",getHistName("Pruned","IVFCSV",true),700,1100);
  graphsPt700ToInf["Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root","ROOT_files/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root",getHistName("Pruned","HybridIVFCSV",false),700,1100);

  orderingPt700ToInf.push_back("Fat Jet CSV (BTV-13-001)");
  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt700ToInf.push_back("Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt300to500_FatJets_Subjets.eps",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt700toInf_FatJets_Subjets.eps",0, 1, 1E-3, 1,1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //------------------------------------------
}
