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
#include "CMS_lumi.C"

using namespace std;


TGraph* getEfficiencyCurve(const string& fFileS1, const string& fFileB1,const string& fPlot1,const string& fPlot2,const double fXMin, const double fXMax)
{
  //get files and histograms
  TFile *file_S1  = new TFile(fFileS1.c_str());
  TFile *file_B1 = new TFile(fFileB1.c_str());

  TH2 *h2_S_1 = (TH2*)file_S1->Get(fPlot1.c_str());
  TH2 *h2_B_1 = (TH2*)file_B1->Get(fPlot2.c_str());

  int nBins = h2_S_1->GetYaxis()->GetNbins();
  int overflow = nBins+1;
  int wideBin = nBins/20;
  int normalBin = nBins/100;
  int extraBins = nBins/normalBin-100;
  int fineBin = nBins/400;
  bool useFineBinning = false;
  if(nBins>400) useFineBinning = true;
  
  //total jet count for denominator of efficiency calculation 
  double denom_S_1 = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),0,overflow);
  double denom_B_1 = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),0,overflow);

  int nPoints = 29+extraBins;
  if(useFineBinning) nPoints = 32+extraBins;

  TGraph *g_eff_1 = new TGraph(nPoints);

  for(int i = 1; i<=5; ++i)
  {
    double num_S = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),1+(i-1)*normalBin,overflow);
    double num_B = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),1+(i-1)*normalBin,overflow);
  
    g_eff_1->SetPoint(i-1,(num_S/denom_S_1),(num_B/denom_B_1));
  }
  for(int i = 6; i<=23; ++i)
  {
    double num_S = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),1+5*normalBin+(i-6)*wideBin,overflow);
    double num_B = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),1+5*normalBin+(i-6)*wideBin,overflow);
  
    g_eff_1->SetPoint(i-1,(num_S/denom_S_1),(num_B/denom_B_1));
  }
  if(useFineBinning)  
  {
    for(int i = 24; i<=27; ++i)
    {
      double num_S = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),1+5*normalBin+18*wideBin+(i-24)*normalBin,overflow);
      double num_B = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),1+5*normalBin+18*wideBin+(i-24)*normalBin,overflow);

      g_eff_1->SetPoint(i-1,(num_S/denom_S_1),(num_B/denom_B_1));
    }
    for(int i = 28; i<=31; ++i)
    {
      double num_S = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),1+5*normalBin+18*wideBin+4*normalBin+(i-28)*fineBin,overflow);
      double num_B = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),1+5*normalBin+18*wideBin+4*normalBin+(i-28)*fineBin,overflow);

      g_eff_1->SetPoint(i-1,(num_S/denom_S_1),(num_B/denom_B_1));
    }
    for(int i = 32; i<=(32+extraBins); ++i)
    {
      double num_S = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),1+5*normalBin+18*wideBin+4*normalBin+4*fineBin+(i-32)*normalBin,overflow);
      double num_B = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),1+5*normalBin+18*wideBin+4*normalBin+4*fineBin+(i-32)*normalBin,overflow);

      g_eff_1->SetPoint(i-1,(num_S/denom_S_1),(num_B/denom_B_1));
    }
  }
  else
  {
    for(int i = 24; i<=(29+extraBins); ++i)
    {
      double num_S = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),1+5*normalBin+18*wideBin+(i-24)*normalBin,overflow);
      double num_B = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),1+5*normalBin+18*wideBin+(i-24)*normalBin,overflow);

      g_eff_1->SetPoint(i-1,(num_S/denom_S_1),(num_B/denom_B_1));
    }
  }

  return g_eff_1;
}

TGraph* getEfficiencyCurve(const string& fFileS1, const string& fFileB1, const string& fPlot, const double fXMin, const double fXMax)
{
  return getEfficiencyCurve(fFileS1, fFileB1, fPlot, fPlot, fXMin, fXMax);
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

std::string getHistName(const std::string & grooming = "Pruned", const std::string & algo = "JetIVFCSV", const std::string & postfix = "", const std::string & size = "")
{
    return "jetAnalyzer"+size+"FatJets_"+grooming+"Subjets" + (postfix != "" ? "_" + postfix : "") + "/h2_JetPt_"+algo+"_BosonMatched_JetMass";
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
  c->SetGridx();
  c->SetGridy();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(1.1,"Y");
  bkg->Draw();

  TLegend *legend = new TLegend(.16,.75-float(graphs.size())*(0.05),.36,.75);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);

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

  //  l1.SetTextAlign(12);
  //l1.SetTextSize(0.045);
  //l1.SetTextFont(62);
  //l1.DrawLatex(0.14,0.96, "CMS Simulation Preliminary, #sqrt{s} = 8 TeV");

  int iPeriod = 2;
  int iPos = 0;

  CMS_lumi( c, iPeriod, iPos );

  l1.SetTextFont(42);
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.48,0.2, fExtraInfo.c_str());
  c->SaveAs(fOutputFile.c_str());
  graphs.clear();
  delete c;
  delete legend;
  delete bkg;
}

void makePlots_AK_AllsubjetAllsetting()
{
  // for multiple plots on the same canvas

  // vectors storing the order of legend entries
  std::vector< std::string> orderingPt200To400;
  // maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt200To400;

  //==========================================
  // Post-BTV-13-001 setup
  //==========================================
  // Fat jets
  //------------------------------------------
  //graphsPt200To400["Fat Jet CSV"]                   = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),200,400);
  //graphsPt200To400["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  //graphsPt200To400["Fat Jet JP"]                    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJP"),200,400);
  //graphsPt200To400["Fat Jet JBP"]                   = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJBP"),200,400);
  graphsPt200To400["Fat Jet CSV (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Fat Jet JP (Explicit JTA)"]     = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJP"),200,400);
  graphsPt200To400["Fat Jet JBP (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJBP"),200,400);

  //orderingPt200To400.push_back("Fat Jet CSV");
  //orderingPt200To400.push_back("Fat Jet IVFCSV");
  //orderingPt200To400.push_back("Fat Jet JP");
  //orderingPt200To400.push_back("Fat Jet JBP");
  orderingPt200To400.push_back("Fat Jet CSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet JP (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet JBP (Explicit JTA)");


  //------------------------------------------
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt200to400_FatJets_comparison.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets - bJets QCD
  //------------------------------------------
  //graphsPt200To400["Fat Jet CSV"]                   = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","bJets"),200,400);
  //graphsPt200To400["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  //graphsPt200To400["Fat Jet JP"]                    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJP"),getHistName("Pruned","JetJP","bJets"),200,400);
  //graphsPt200To400["Fat Jet JBP"]                   = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJBP"),getHistName("Pruned","JetJBP","bJets"),200,400);
  graphsPt200To400["Fat Jet CSV (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","bJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  graphsPt200To400["Fat Jet JP (Explicit JTA)"]     = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJP"),getHistName("Pruned","JetJP","bJets"),200,400);
  graphsPt200To400["Fat Jet JBP (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJBP"),getHistName("Pruned","JetJBP","bJets"),200,400);

  //orderingPt200To400.push_back("Fat Jet CSV");
  //orderingPt200To400.push_back("Fat Jet IVFCSV");
  //orderingPt200To400.push_back("Fat Jet JP");
  //orderingPt200To400.push_back("Fat Jet JBP");
  orderingPt200To400.push_back("Fat Jet CSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet JP (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet JBP (Explicit JTA)");


  //------------------------------------------
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_FatJets_comparison_bJetsQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets - bJetsGSP QCD
  //------------------------------------------
  //graphsPt200To400["Fat Jet CSV"]                   = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","bJetsGSP"),200,400);
  //graphsPt200To400["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  //graphsPt200To400["Fat Jet JP"]                    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJP"),getHistName("Pruned","JetJP","bJetsGSP"),200,400);
  //graphsPt200To400["Fat Jet JBP"]                   = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJBP"),getHistName("Pruned","JetJBP","bJetsGSP"),200,400);
  graphsPt200To400["Fat Jet CSV (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","bJetsGSP"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Fat Jet JP (Explicit JTA)"]     = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJP"),getHistName("Pruned","JetJP","bJetsGSP"),200,400);
  graphsPt200To400["Fat Jet JBP (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJBP"),getHistName("Pruned","JetJBP","bJetsGSP"),200,400);

  //orderingPt200To400.push_back("Fat Jet CSV");
  //orderingPt200To400.push_back("Fat Jet IVFCSV");
  //orderingPt200To400.push_back("Fat Jet JP");
  //orderingPt200To400.push_back("Fat Jet JBP");
  orderingPt200To400.push_back("Fat Jet CSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet JP (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet JBP (Explicit JTA)");

  //------------------------------------------
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_FatJets_comparison_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets - cJets QCD
  //------------------------------------------
  //graphsPt200To400["Fat Jet CSV"]                   = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","cJets"),200,400);
  //graphsPt200To400["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  //graphsPt200To400["Fat Jet JP"]                    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJP"),getHistName("Pruned","JetJP","cJets"),200,400);
  //graphsPt200To400["Fat Jet JBP"]                   = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJBP"),getHistName("Pruned","JetJBP","cJets"),200,400);
  graphsPt200To400["Fat Jet CSV (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","cJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  graphsPt200To400["Fat Jet JP (Explicit JTA)"]     = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJP"),getHistName("Pruned","JetJP","cJets"),200,400);
  graphsPt200To400["Fat Jet JBP (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJBP"),getHistName("Pruned","JetJBP","cJets"),200,400);

  //orderingPt200To400.push_back("Fat Jet CSV");
  //orderingPt200To400.push_back("Fat Jet IVFCSV");
  //orderingPt200To400.push_back("Fat Jet JP");
  //orderingPt200To400.push_back("Fat Jet JBP");
  orderingPt200To400.push_back("Fat Jet CSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet JP (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet JBP (Explicit JTA)");

  //------------------------------------------
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_FatJets_comparison_cJetsQCD.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets - udsJets QCD
  //------------------------------------------
  //graphsPt200To400["Fat Jet CSV"]                   = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","udsJets"),200,400);
  //graphsPt200To400["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  //graphsPt200To400["Fat Jet JP"]                    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJP"),getHistName("Pruned","JetJP","udsJets"),200,400);
  //graphsPt200To400["Fat Jet JBP"]                   = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJBP"),getHistName("Pruned","JetJBP","udsJets"),200,400);
  graphsPt200To400["Fat Jet CSV (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","udsJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  graphsPt200To400["Fat Jet JP (Explicit JTA)"]     = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJP"),getHistName("Pruned","JetJP","udsJets"),200,400);
  graphsPt200To400["Fat Jet JBP (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJBP"),getHistName("Pruned","JetJBP","udsJets"),200,400);

  //orderingPt200To400.push_back("Fat Jet CSV");
  //orderingPt200To400.push_back("Fat Jet IVFCSV");
  //orderingPt200To400.push_back("Fat Jet JP");
  //orderingPt200To400.push_back("Fat Jet JBP");
  orderingPt200To400.push_back("Fat Jet CSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet JP (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet JBP (Explicit JTA)");


  //------------------------------------------
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_FatJets_comparison_udsJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets - gluonJets QCD
  //------------------------------------------
  //graphsPt200To400["Fat Jet CSV"]                   = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","gluonJets"),200,400);
  //graphsPt200To400["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  //graphsPt200To400["Fat Jet JP"]                    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJP"),getHistName("Pruned","JetJP","gluonJets"),200,400);
  //graphsPt200To400["Fat Jet JBP"]                   = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJBP"),getHistName("Pruned","JetJBP","gluonJets"),200,400);
  graphsPt200To400["Fat Jet CSV (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","gluonJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Fat Jet JP (Explicit JTA)"]     = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJP"),getHistName("Pruned","JetJP","gluonJets"),200,400);
  graphsPt200To400["Fat Jet JBP (Explicit JTA)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetJBP"),getHistName("Pruned","JetJBP","gluonJets"),200,400);

  //orderingPt200To400.push_back("Fat Jet CSV");
  //orderingPt200To400.push_back("Fat Jet IVFCSV");
  //orderingPt200To400.push_back("Fat Jet JP");
  //orderingPt200To400.push_back("Fat Jet JBP");
  orderingPt200To400.push_back("Fat Jet CSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet JP (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet JBP (Explicit JTA)");


  //------------------------------------------
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_FatJets_comparison_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets: IVFCSV (Cone vs Explicit JTA)
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  
  orderingPt200To400.push_back("Fat Jet IVFCSV");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD)","","btagperfcomp_Pt200to400_FatJets_IVFCSV_Cone_vs_Expl_JTA.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets: IVFCSV (Cone vs Explicit JTA) - bJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  
  orderingPt200To400.push_back("Fat Jet IVFCSV");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_FatJets_IVFCSV_Cone_vs_Expl_JTA_bJetsQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets: IVFCSV (Cone vs Explicit JTA) - bJetsGSP QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  
  orderingPt200To400.push_back("Fat Jet IVFCSV");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_FatJets_IVFCSV_Cone_vs_Expl_JTA_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets: IVFCSV (Cone vs Explicit JTA) - cJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  
  orderingPt200To400.push_back("Fat Jet IVFCSV");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_FatJets_IVFCSV_Cone_vs_Expl_JTA_cJetsQCD.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets: IVFCSV (Cone vs Explicit JTA) - udsJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  
  orderingPt200To400.push_back("Fat Jet IVFCSV");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_FatJets_IVFCSV_Cone_vs_Expl_JTA_udsJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets: IVFCSV (Cone vs Explicit JTA) - gluonJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV"]                = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  
  orderingPt200To400.push_back("Fat Jet IVFCSV");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  //------------------------------------------
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_FatJets_IVFCSV_Cone_vs_Expl_JTA_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets: IVFCSV Explicit JTA comparison
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  
  orderingPt200To400.push_back("Fat Jet IVFCSV");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA signal only)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD)","","btagperfcomp_Pt200to400_FatJets_IVFCSV_Expl_JTA_comparison.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets: IVFCSV Explicit JTA comparison - bJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  
  orderingPt200To400.push_back("Fat Jet IVFCSV");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA signal only)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_FatJets_IVFCSV_Expl_JTA_comparison_bJetsQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets: IVFCSV Explicit JTA comparison - bJetsGSP QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  
  orderingPt200To400.push_back("Fat Jet IVFCSV");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA signal only)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_FatJets_IVFCSV_Expl_JTA_comparison_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();
  
  //------------------------------------------
  // Fat jets: IVFCSV Explicit JTA comparison - cJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  
  orderingPt200To400.push_back("Fat Jet IVFCSV");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA signal only)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_FatJets_IVFCSV_Expl_JTA_comparison_cJetsQCD.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets: IVFCSV Explicit JTA comparison - udsJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  
  orderingPt200To400.push_back("Fat Jet IVFCSV");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA signal only)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_FatJets_IVFCSV_Expl_JTA_comparison_udsJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets: IVFCSV Explicit JTA comparison - gluonJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  
  orderingPt200To400.push_back("Fat Jet IVFCSV");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA signal only)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_FatJets_IVFCSV_Expl_JTA_comparison_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: CSV
  //------------------------------------------
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),200,400);
  graphsPt200To400["Subjet CSV (Filtered)"]           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinCSV"),200,400);
  //graphsPt200To400["Subjet CSV (MD+Filtered)"]        = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinCSV"),200,400);
  graphsPt200To400["Subjet CSV (k_{T})"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinCSV"),200,400);
  graphsPt200To400["Subjet CSV (k_{T}+Filtered)"]     = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinCSV"),200,400);

  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet CSV (Filtered)");
  //orderingPt200To400.push_back("Subjet CSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet CSV (k_{T})");
  orderingPt200To400.push_back("Subjet CSV (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD)","","btagperfcomp_Pt200to400_Subjets_CSV.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: CSV - bJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","bJets"),200,400);
  graphsPt200To400["Subjet CSV (Filtered)"]           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinCSV"),getHistName("Filtered","SubJetMinCSV","bJets"),200,400);
  //graphsPt200To400["Subjet CSV (MD+Filtered)"]        = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinCSV"),200,400);
  graphsPt200To400["Subjet CSV (k_{T})"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinCSV"),getHistName("Kt","SubJetMinCSV","bJets"),200,400);
  graphsPt200To400["Subjet CSV (k_{T}+Filtered)"]     = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinCSV"),getHistName("KtBDRSFiltered","SubJetMinCSV","bJets"),200,400);

  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet CSV (Filtered)");
  //orderingPt200To400.push_back("Subjet CSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet CSV (k_{T})");
  orderingPt200To400.push_back("Subjet CSV (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_Subjets_CSV_bJetsQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: CSV - bJetsGSP QCD
  //------------------------------------------
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet CSV (Filtered)"]           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinCSV"),getHistName("Filtered","SubJetMinCSV","bJetsGSP"),200,400);
  //graphsPt200To400["Subjet CSV (MD+Filtered)"]        = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinCSV"),200,400);
  graphsPt200To400["Subjet CSV (k_{T})"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinCSV"),getHistName("Kt","SubJetMinCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet CSV (k_{T}+Filtered)"]     = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinCSV"),getHistName("KtBDRSFiltered","SubJetMinCSV","bJetsGSP"),200,400);

  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet CSV (Filtered)");
  //orderingPt200To400.push_back("Subjet CSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet CSV (k_{T})");
  orderingPt200To400.push_back("Subjet CSV (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_Subjets_CSV_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();
  orderingPt200To400.clear();
  
  //------------------------------------------
  // Subjets: CSV - cJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","cJets"),200,400);
  graphsPt200To400["Subjet CSV (Filtered)"]           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinCSV"),getHistName("Filtered","SubJetMinCSV","cJets"),200,400);
  //graphsPt200To400["Subjet CSV (MD+Filtered)"]        = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinCSV"),200,400);
  graphsPt200To400["Subjet CSV (k_{T})"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinCSV"),getHistName("Kt","SubJetMinCSV","cJets"),200,400);
  graphsPt200To400["Subjet CSV (k_{T}+Filtered)"]     = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinCSV"),getHistName("KtBDRSFiltered","SubJetMinCSV","cJets"),200,400);

  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet CSV (Filtered)");
  //orderingPt200To400.push_back("Subjet CSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet CSV (k_{T})");
  orderingPt200To400.push_back("Subjet CSV (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_Subjets_CSV_cJetsQCD.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: CSV - udsJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","udsJets"),200,400);
  graphsPt200To400["Subjet CSV (Filtered)"]           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinCSV"),getHistName("Filtered","SubJetMinCSV","udsJets"),200,400);
  //graphsPt200To400["Subjet CSV (MD+Filtered)"]        = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinCSV"),200,400);
  graphsPt200To400["Subjet CSV (k_{T})"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinCSV"),getHistName("Kt","SubJetMinCSV","udsJets"),200,400);
  graphsPt200To400["Subjet CSV (k_{T}+Filtered)"]     = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinCSV"),getHistName("KtBDRSFiltered","SubJetMinCSV","udsJets"),200,400);

  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet CSV (Filtered)");
  //orderingPt200To400.push_back("Subjet CSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet CSV (k_{T})");
  orderingPt200To400.push_back("Subjet CSV (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_Subjets_CSV_udsJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: CSV - gluonJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet CSV (Filtered)"]           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinCSV"),getHistName("Filtered","SubJetMinCSV","gluonJets"),200,400);
  //graphsPt200To400["Subjet CSV (MD+Filtered)"]        = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinCSV"),200,400);
  graphsPt200To400["Subjet CSV (k_{T})"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinCSV"),getHistName("Kt","SubJetMinCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet CSV (k_{T}+Filtered)"]     = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinCSV"),getHistName("KtBDRSFiltered","SubJetMinCSV","gluonJets"),200,400);

  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet CSV (Filtered)");
  //orderingPt200To400.push_back("Subjet CSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet CSV (k_{T})");
  orderingPt200To400.push_back("Subjet CSV (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_Subjets_CSV_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: IVFCSV
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinIVFCSV"),200,400);
  //graphsPt200To400["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinIVFCSV"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Filtered)");
  //orderingPt200To400.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T})");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
 
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD)","","btagperfcomp_Pt200to400_Subjets_IVFCSV.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: IVFCSV - bJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinIVFCSV"),getHistName("Filtered","SubJetMinIVFCSV","bJets"),200,400);
  //graphsPt200To400["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinIVFCSV"),getHistName("KtBDRSFiltered","SubJetMinIVFCSV","bJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Filtered)");
  //orderingPt200To400.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T})");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
 
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_Subjets_IVFCSV_bJetsQCD.eps",0, 1, 1E-3, 1,0);
 
  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: IVFCSV - bJetsGSP QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinIVFCSV"),getHistName("Filtered","SubJetMinIVFCSV","bJetsGSP"),200,400);
  //graphsPt200To400["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinIVFCSV"),getHistName("KtBDRSFiltered","SubJetMinIVFCSV","bJetsGSP"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Filtered)");
  //orderingPt200To400.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T})");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
 
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_Subjets_IVFCSV_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);
 
  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: IVFCSV - cJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinIVFCSV"),getHistName("Filtered","SubJetMinIVFCSV","cJets"),200,400);
  //graphsPt200To400["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinIVFCSV"),getHistName("KtBDRSFiltered","SubJetMinIVFCSV","cJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Filtered)");
  //orderingPt200To400.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T})");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
 
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_Subjets_IVFCSV_cJetsQCD.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: IVFCSV - udsJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinIVFCSV"),getHistName("Filtered","SubJetMinIVFCSV","udsJets"),200,400);
  //graphsPt200To400["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinIVFCSV"),getHistName("KtBDRSFiltered","SubJetMinIVFCSV","udsJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Filtered)");
  //orderingPt200To400.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T})");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
 
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_Subjets_IVFCSV_udsJetsQCD.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: IVFCSV - gluonJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinIVFCSV"),getHistName("Filtered","SubJetMinIVFCSV","gluonJets"),200,400);
  //graphsPt200To400["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinIVFCSV"),getHistName("KtBDRSFiltered","SubJetMinIVFCSV","gluonJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Filtered)");
  //orderingPt200To400.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T})");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
 
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_Subjets_IVFCSV_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: IVFCSV - Explicit JTA, SVClustering
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinIVFCSV"),200,400);
  //graphsPt200To400["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinIVFCSV"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Filtered)");
  //orderingPt200To400.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T})");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
 
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD)","Explicit JTA, SVClustering","btagperfcomp_Pt200to400_Subjets_IVFCSV_ExplicitJTA_SVClustering.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: IVFCSV - bJets QCD - ExplicitJTA, SVClustering
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinIVFCSV"),getHistName("Filtered","SubJetMinIVFCSV","bJets"),200,400);
  //graphsPt200To400["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinIVFCSV"),getHistName("KtBDRSFiltered","SubJetMinIVFCSV","bJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Filtered)");
  //orderingPt200To400.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T})");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
 
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","Explicit JTA, SVClustering","btagperfcomp_Pt200to400_Subjets_IVFCSV_ExplicitJTA_SVClustering_bJetsQCD.eps",0, 1, 1E-3, 1,0);
 
  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: IVFCSV - bJetsGSP QCD - ExplicitJTA, SVClustering
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinIVFCSV"),getHistName("Filtered","SubJetMinIVFCSV","bJetsGSP"),200,400);
  //graphsPt200To400["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinIVFCSV"),getHistName("KtBDRSFiltered","SubJetMinIVFCSV","bJetsGSP"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Filtered)");
  //orderingPt200To400.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T})");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
 
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","Explicit JTA, SVClustering","btagperfcomp_Pt200to400_Subjets_IVFCSV_ExplicitJTA_SVClustering_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);
 
  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: IVFCSV - cJets QCD - ExplicitJTA, SVClustering
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinIVFCSV"),getHistName("Filtered","SubJetMinIVFCSV","cJets"),200,400);
  //graphsPt200To400["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinIVFCSV"),getHistName("KtBDRSFiltered","SubJetMinIVFCSV","cJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Filtered)");
  //orderingPt200To400.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T})");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
 
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","Explicit JTA, SVClustering","btagperfcomp_Pt200to400_Subjets_IVFCSV_ExplicitJTA_SVClustering_cJetsQCD.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: IVFCSV - udsJets QCD - ExplicitJTA, SVClustering
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinIVFCSV"),getHistName("Filtered","SubJetMinIVFCSV","udsJets"),200,400);
  //graphsPt200To400["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinIVFCSV"),getHistName("KtBDRSFiltered","SubJetMinIVFCSV","udsJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Filtered)");
  //orderingPt200To400.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T})");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
 
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","Explicit JTA, SVClustering","btagperfcomp_Pt200to400_Subjets_IVFCSV_ExplicitJTA_SVClustering_udsJetsQCD.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: IVFCSV - gluonJets QCD - ExplicitJTA, SVClustering
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinIVFCSV"),getHistName("Filtered","SubJetMinIVFCSV","gluonJets"),200,400);
  //graphsPt200To400["Subjet IVFCSV (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinIVFCSV"),getHistName("KtBDRSFiltered","SubJetMinIVFCSV","gluonJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Filtered)");
  //orderingPt200To400.push_back("Subjet IVFCSV (MD+Filtered)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T})");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //------------------------------------------
 
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","Explicit JTA, SVClustering","btagperfcomp_Pt200to400_Subjets_IVFCSV_ExplicitJTA_SVClustering_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: JP
  //------------------------------------------
  graphsPt200To400["Subjet JP (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJP"),200,400);
  graphsPt200To400["Subjet JP (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinJP"),200,400);
  //graphsPt200To400["Subjet JP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinJP"),200,400);
  graphsPt200To400["Subjet JP (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinJP"),200,400);
  graphsPt200To400["Subjet JP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinJP"),200,400);

  orderingPt200To400.push_back("Subjet JP (Pruned)");
  orderingPt200To400.push_back("Subjet JP (Filtered)");
  //orderingPt200To400.push_back("Subjet JP (MD+Filtered)");
  orderingPt200To400.push_back("Subjet JP (k_{T})");
  orderingPt200To400.push_back("Subjet JP (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD)","","btagperfcomp_Pt200to400_Subjets_JP.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: JP - bJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet JP (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJP"),getHistName("Pruned","SubJetMinJP","bJets"),200,400);
  graphsPt200To400["Subjet JP (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinJP"),getHistName("Filtered","SubJetMinJP","bJets"),200,400);
  //graphsPt200To400["Subjet JP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinJP"),200,400);
  graphsPt200To400["Subjet JP (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinJP"),getHistName("Kt","SubJetMinJP","bJets"),200,400);
  graphsPt200To400["Subjet JP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinJP"),getHistName("KtBDRSFiltered","SubJetMinJP","bJets"),200,400);

  orderingPt200To400.push_back("Subjet JP (Pruned)");
  orderingPt200To400.push_back("Subjet JP (Filtered)");
  //orderingPt200To400.push_back("Subjet JP (MD+Filtered)");
  orderingPt200To400.push_back("Subjet JP (k_{T})");
  orderingPt200To400.push_back("Subjet JP (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_Subjets_JP_bJetsQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: JP - bJetsGSP QCD
  //------------------------------------------
  graphsPt200To400["Subjet JP (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJP"),getHistName("Pruned","SubJetMinJP","bJetsGSP"),200,400);
  graphsPt200To400["Subjet JP (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinJP"),getHistName("Filtered","SubJetMinJP","bJetsGSP"),200,400);
  //graphsPt200To400["Subjet JP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinJP"),200,400);
  graphsPt200To400["Subjet JP (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinJP"),getHistName("Kt","SubJetMinJP","bJetsGSP"),200,400);
  graphsPt200To400["Subjet JP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinJP"),getHistName("KtBDRSFiltered","SubJetMinJP","bJetsGSP"),200,400);

  orderingPt200To400.push_back("Subjet JP (Pruned)");
  orderingPt200To400.push_back("Subjet JP (Filtered)");
  //orderingPt200To400.push_back("Subjet JP (MD+Filtered)");
  orderingPt200To400.push_back("Subjet JP (k_{T})");
  orderingPt200To400.push_back("Subjet JP (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_Subjets_JP_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: JP - cJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet JP (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJP"),getHistName("Pruned","SubJetMinJP","cJets"),200,400);
  graphsPt200To400["Subjet JP (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinJP"),getHistName("Filtered","SubJetMinJP","cJets"),200,400);
  //graphsPt200To400["Subjet JP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinJP"),200,400);
  graphsPt200To400["Subjet JP (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinJP"),getHistName("Kt","SubJetMinJP","cJets"),200,400);
  graphsPt200To400["Subjet JP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinJP"),getHistName("KtBDRSFiltered","SubJetMinJP","cJets"),200,400);

  orderingPt200To400.push_back("Subjet JP (Pruned)");
  orderingPt200To400.push_back("Subjet JP (Filtered)");
  //orderingPt200To400.push_back("Subjet JP (MD+Filtered)");
  orderingPt200To400.push_back("Subjet JP (k_{T})");
  orderingPt200To400.push_back("Subjet JP (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_Subjets_JP_cJetsQCD.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();
 
  //------------------------------------------
  // Subjets: JP - udsJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet JP (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJP"),getHistName("Pruned","SubJetMinJP","udsJets"),200,400);
  graphsPt200To400["Subjet JP (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinJP"),getHistName("Filtered","SubJetMinJP","udsJets"),200,400);
  //graphsPt200To400["Subjet JP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinJP"),200,400);
  graphsPt200To400["Subjet JP (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinJP"),getHistName("Kt","SubJetMinJP","udsJets"),200,400);
  graphsPt200To400["Subjet JP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinJP"),getHistName("KtBDRSFiltered","SubJetMinJP","udsJets"),200,400);

  orderingPt200To400.push_back("Subjet JP (Pruned)");
  orderingPt200To400.push_back("Subjet JP (Filtered)");
  //orderingPt200To400.push_back("Subjet JP (MD+Filtered)");
  orderingPt200To400.push_back("Subjet JP (k_{T})");
  orderingPt200To400.push_back("Subjet JP (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_Subjets_JP_udsJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: JP - gluonJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet JP (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJP"),getHistName("Pruned","SubJetMinJP","gluonJets"),200,400);
  graphsPt200To400["Subjet JP (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinJP"),getHistName("Filtered","SubJetMinJP","gluonJets"),200,400);
  //graphsPt200To400["Subjet JP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinJP"),200,400);
  graphsPt200To400["Subjet JP (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinJP"),getHistName("Kt","SubJetMinJP","gluonJets"),200,400);
  graphsPt200To400["Subjet JP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinJP"),getHistName("KtBDRSFiltered","SubJetMinJP","gluonJets"),200,400);

  orderingPt200To400.push_back("Subjet JP (Pruned)");
  orderingPt200To400.push_back("Subjet JP (Filtered)");
  //orderingPt200To400.push_back("Subjet JP (MD+Filtered)");
  orderingPt200To400.push_back("Subjet JP (k_{T})");
  orderingPt200To400.push_back("Subjet JP (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_Subjets_JP_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: JBP
  //------------------------------------------
  graphsPt200To400["Subjet JBP (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJBP"),200,400);
  graphsPt200To400["Subjet JBP (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinJBP"),200,400);
  //graphsPt200To400["Subjet JBP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinJBP"),200,400);
  graphsPt200To400["Subjet JBP (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinJBP"),200,400);
  graphsPt200To400["Subjet JBP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinJBP"),200,400);

  orderingPt200To400.push_back("Subjet JBP (Pruned)");
  orderingPt200To400.push_back("Subjet JBP (Filtered)");
  //orderingPt200To400.push_back("Subjet JBP (MD+Filtered)");
  orderingPt200To400.push_back("Subjet JBP (k_{T})");
  orderingPt200To400.push_back("Subjet JBP (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD)","","btagperfcomp_Pt200to400_Subjets_JBP.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //CONTINUE HERE---------------

  //------------------------------------------
  // Subjets: JBP - bJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet JBP (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJBP"),getHistName("Pruned","SubJetMinJBP","bJets"),200,400);
  graphsPt200To400["Subjet JBP (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinJBP"),getHistName("Filtered","SubJetMinJBP","bJets"),200,400);
  //graphsPt200To400["Subjet JBP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinJBP"),200,400);
  graphsPt200To400["Subjet JBP (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinJBP"),getHistName("Kt","SubJetMinJBP","bJets"),200,400);
  graphsPt200To400["Subjet JBP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinJBP"),getHistName("KtBDRSFiltered","SubJetMinJBP","bJets"),200,400);

  orderingPt200To400.push_back("Subjet JBP (Pruned)");
  orderingPt200To400.push_back("Subjet JBP (Filtered)");
  //orderingPt200To400.push_back("Subjet JBP (MD+Filtered)");
  orderingPt200To400.push_back("Subjet JBP (k_{T})");
  orderingPt200To400.push_back("Subjet JBP (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_Subjets_JBP_bJetsQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: JBP - bJetsGSP QCD
  //------------------------------------------
  graphsPt200To400["Subjet JBP (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJBP"),getHistName("Pruned","SubJetMinJBP","bJetsGSP"),200,400);
  graphsPt200To400["Subjet JBP (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinJBP"),getHistName("Filtered","SubJetMinJBP","bJetsGSP"),200,400);
  //graphsPt200To400["Subjet JBP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinJBP"),200,400);
  graphsPt200To400["Subjet JBP (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinJBP"),getHistName("Kt","SubJetMinJBP","bJetsGSP"),200,400);
  graphsPt200To400["Subjet JBP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinJBP"),getHistName("KtBDRSFiltered","SubJetMinJBP","bJetsGSP"),200,400);

  orderingPt200To400.push_back("Subjet JBP (Pruned)");
  orderingPt200To400.push_back("Subjet JBP (Filtered)");
  //orderingPt200To400.push_back("Subjet JBP (MD+Filtered)");
  orderingPt200To400.push_back("Subjet JBP (k_{T})");
  orderingPt200To400.push_back("Subjet JBP (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_Subjets_JBP_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: JBP - cJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet JBP (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJBP"),getHistName("Pruned","SubJetMinJBP","cJets"),200,400);
  graphsPt200To400["Subjet JBP (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinJBP"),getHistName("Filtered","SubJetMinJBP","cJets"),200,400);
  //graphsPt200To400["Subjet JBP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinJBP"),200,400);
  graphsPt200To400["Subjet JBP (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinJBP"),getHistName("Kt","SubJetMinJBP","cJets"),200,400);
  graphsPt200To400["Subjet JBP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinJBP"),getHistName("KtBDRSFiltered","SubJetMinJBP","cJets"),200,400);

  orderingPt200To400.push_back("Subjet JBP (Pruned)");
  orderingPt200To400.push_back("Subjet JBP (Filtered)");
  //orderingPt200To400.push_back("Subjet JBP (MD+Filtered)");
  orderingPt200To400.push_back("Subjet JBP (k_{T})");
  orderingPt200To400.push_back("Subjet JBP (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_Subjets_JBP_cJetsQCD.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: JBP - udsJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet JBP (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJBP"),getHistName("Pruned","SubJetMinJBP","udsJets"),200,400);
  graphsPt200To400["Subjet JBP (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinJBP"),getHistName("Filtered","SubJetMinJBP","udsJets"),200,400);
  //graphsPt200To400["Subjet JBP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinJBP"),200,400);
  graphsPt200To400["Subjet JBP (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinJBP"),getHistName("Kt","SubJetMinJBP","udsJets"),200,400);
  graphsPt200To400["Subjet JBP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinJBP"),getHistName("KtBDRSFiltered","SubJetMinJBP","udsJets"),200,400);

  orderingPt200To400.push_back("Subjet JBP (Pruned)");
  orderingPt200To400.push_back("Subjet JBP (Filtered)");
  //orderingPt200To400.push_back("Subjet JBP (MD+Filtered)");
  orderingPt200To400.push_back("Subjet JBP (k_{T})");
  orderingPt200To400.push_back("Subjet JBP (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_Subjets_JBP_udsJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: JBP - gluonJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet JBP (Pruned)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJBP"),getHistName("Pruned","SubJetMinJBP","gluonJets"),200,400);
  graphsPt200To400["Subjet JBP (Filtered)"]       = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Filtered","SubJetMinJBP"),getHistName("Filtered","SubJetMinJBP","gluonJets"),200,400);
  //graphsPt200To400["Subjet JBP (MD+Filtered)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("MDBDRSFiltered","SubJetMinJBP"),200,400);
  graphsPt200To400["Subjet JBP (k_{T})"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinJBP"),getHistName("Kt","SubJetMinJBP","gluonJets"),200,400);
  graphsPt200To400["Subjet JBP (k_{T}+Filtered)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("KtBDRSFiltered","SubJetMinJBP"),getHistName("KtBDRSFiltered","SubJetMinJBP","gluonJets"),200,400);

  orderingPt200To400.push_back("Subjet JBP (Pruned)");
  orderingPt200To400.push_back("Subjet JBP (Filtered)");
  //orderingPt200To400.push_back("Subjet JBP (MD+Filtered)");
  orderingPt200To400.push_back("Subjet JBP (k_{T})");
  orderingPt200To400.push_back("Subjet JBP (k_{T}+Filtered)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_Subjets_JBP_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned comparison
  //------------------------------------------
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned)"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet JP (Pruned)"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJP"),200,400);
  graphsPt200To400["Subjet JBP (Pruned)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJBP"),200,400);
  
  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet JP (Pruned)");
  orderingPt200To400.push_back("Subjet JBP (Pruned)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD)","","btagperfcomp_Pt200to400_Subjets_Pruned_comparison.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned comparison - bJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned)"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet JP (Pruned)"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJP"),getHistName("Pruned","SubJetMinJP","bJets"),200,400);
  graphsPt200To400["Subjet JBP (Pruned)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJBP"),getHistName("Pruned","SubJetMinJBP","bJets"),200,400);
  
  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet JP (Pruned)");
  orderingPt200To400.push_back("Subjet JBP (Pruned)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_Subjets_Pruned_comparison_bJetsQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned comparison - bJetsGSP QCD
  //------------------------------------------
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned)"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet JP (Pruned)"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJP"),getHistName("Pruned","SubJetMinJP","bJetsGSP"),200,400);
  graphsPt200To400["Subjet JBP (Pruned)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJBP"),getHistName("Pruned","SubJetMinJBP","bJetsGSP"),200,400);
  
  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet JP (Pruned)");
  orderingPt200To400.push_back("Subjet JBP (Pruned)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_Subjets_Pruned_comparison_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned comparison - cJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned)"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet JP (Pruned)"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJP"),getHistName("Pruned","SubJetMinJP","cJets"),200,400);
  graphsPt200To400["Subjet JBP (Pruned)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJBP"),getHistName("Pruned","SubJetMinJBP","cJets"),200,400);
  
  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet JP (Pruned)");
  orderingPt200To400.push_back("Subjet JBP (Pruned)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_Subjets_Pruned_comparison_cJetsQCD.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();
  
  //------------------------------------------
  // Subjets: Pruned comparison - udsJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned)"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet JP (Pruned)"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJP"),getHistName("Pruned","SubJetMinJP","udsJets"),200,400);
  graphsPt200To400["Subjet JBP (Pruned)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJBP"),getHistName("Pruned","SubJetMinJBP","udsJets"),200,400);
  
  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet JP (Pruned)");
  orderingPt200To400.push_back("Subjet JBP (Pruned)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_Subjets_Pruned_comparison_udsJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned comparison - gluonJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned)"]          = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet JP (Pruned)"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJP"),getHistName("Pruned","SubJetMinJP","gluonJets"),200,400);
  graphsPt200To400["Subjet JBP (Pruned)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinJBP"),getHistName("Pruned","SubJetMinJBP","gluonJets"),200,400);
  
  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet JP (Pruned)");
  orderingPt200To400.push_back("Subjet JBP (Pruned)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_Subjets_Pruned_comparison_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV comparison
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]                                           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA)"]                             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD)","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_comparison.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV comparison - bJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]                                           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA)"]                             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_comparison_bJetsQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV comparison - bJetsGSP
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]                                           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA)"]                             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_comparison_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV comparison - cJets
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]                                           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA)"]                             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_comparison_cJetsQCD.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV comparison - udsJets
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]                                           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA)"]                             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_comparison_udsJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV comparison - gluonJets
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]                                           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA)"]                             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]              = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_comparison_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV Explicit JTA comparison
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]                           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD)","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_Expl_JTA_comparison.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV Explicit JTA comparison - bJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]                           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_Expl_JTA_comparison_bJetsQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV Explicit JTA comparison - bJetsGSP QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]                           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_Expl_JTA_comparison_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV Explicit JTA comparison - cJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]                           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_Expl_JTA_comparison_cJetsQCD.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV Explicit JTA comparison - udsJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]                           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_Expl_JTA_comparison_udsJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV Explicit JTA comparison - gluonJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]                           = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_Expl_JTA_comparison_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV SV comparison
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_SV_comparison.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV SV comparison - bJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_SV_comparison_bJetsQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV SV comparison - bJetsGSP QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_SV_comparison_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();
  
  //------------------------------------------
  // Subjets: Pruned IVFCSV SV comparison - cJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_SV_comparison_cJetsQCD.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV SV comparison - udsJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_SV_comparison_udsJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV SV comparison - gluonJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, SV Clustering bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, SV Clustering bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_SV_comparison_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV, Explicit JTA, SV comparison
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD)","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_Expl_JTA_SV_comparison.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV, Explicit JTA, SV comparison - bJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_Expl_JTA_SV_comparison_bJetsQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV, Explicit JTA, SV comparison - bJetsGSP QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_Expl_JTA_SV_comparison_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV, Explicit JTA, SV comparison - cJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_Expl_JTA_SV_comparison_cJetsQCD.eps",0, 1, 1E-3, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV, Explicit JTA, SV comparison - udsJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_Expl_JTA_SV_comparison_udsJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Subjets: Pruned IVFCSV, Explicit JTA, SV comparison - gluonJets QCD
  //------------------------------------------
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA)"]                            = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)"]    = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);

  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)");
  //------------------------------------------

  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_Subjets_Pruned_IVFCSV_Expl_JTA_SV_comparison_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);

  graphsPt200To400.clear();

  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets and subjets
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"]                             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  //  graphsPt200To400["Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","HybridIVFCSV"),200,400);
   
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //orderingPt200To400.push_back("Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)");
  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD)","","btagperfcomp_Pt200to400_FatJets_Subjets.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets and subjets - bJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"]                             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  //  graphsPt200To400["Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","HybridIVFCSV"),getHistName("Pruned","HybridIVFCSV","bJets"),200,400);
   
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //orderingPt200To400.push_back("Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)");
  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJets","","btagperfcomp_Pt200to400_FatJets_Subjets_bJetsQCD.eps",0, 1, 1E-3, 1,0);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets and subjets - bJetsGSP QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"]                             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  //graphsPt200To400["Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","HybridIVFCSV"),getHistName("Pruned","HybridIVFCSV","bJetsGSP"),200,400);
   
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //orderingPt200To400.push_back("Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)");
  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - bJetsGSP","","btagperfcomp_Pt200to400_FatJets_Subjets_bJetsGSPQCD.eps",0, 1, 1E-3, 1,0);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();
   
  //------------------------------------------
  // Fat jets and subjets - cJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"]                             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  //graphsPt200To400["Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","HybridIVFCSV"),getHistName("Pruned","HybridIVFCSV","cJets"),200,400);
   
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //orderingPt200To400.push_back("Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)");
  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - cJets","","btagperfcomp_Pt200to400_FatJets_Subjets_cJetsQCD.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets and subjets - udsJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"]                             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  //graphsPt200To400["Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","HybridIVFCSV"),getHistName("Pruned","HybridIVFCSV","udsJets"),200,400);
   
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //orderingPt200To400.push_back("Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)");
  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - udsJets","","btagperfcomp_Pt200to400_FatJets_Subjets_udsJetsQCD.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets and subjets - gluonJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet CSV (Pruned, BTV-13-001)"]                             = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  //graphsPt200To400["Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","HybridIVFCSV"),getHistName("Pruned","HybridIVFCSV","gluonJets"),200,400);
   
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //orderingPt200To400.push_back("Hybrid IVFCSV (Pruned subjets, Explicit JTA, SV Clustering)");
  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "MisID prob. (QCD) - gluonJets","","btagperfcomp_Pt200to400_FatJets_Subjets_gluonJetsQCD.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

}

void makePlots_CA_AK_comparison()
{
  // for multiple plots on the same canvas

  // vectors storing the order of legend entries
  std::vector< std::string> orderingPt200To400;
  // maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt200To400;
  
  //==========================================
  // Post-BTV-13-001 setup
  //==========================================

  //------------------------------------------
  // Fat jets and subjets - Explicit JTA, SVClustering - CA&AK 
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA"]                               = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_5_CA12_ready4comparison/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA"]         = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_5_CA12_ready4comparison/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"); 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA");

  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Inclusive QCD)","","btagperfcomp_Pt200to400_FatJets_Subjets_CA_vs_AK.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets and subjets - Explicit JTA, SVClustering - CA&AK - bJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA"]                               = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_5_CA12_ready4comparison/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA"]         = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_5_CA12_ready4comparison/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK");   
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA");
 
  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (b jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_CA_vs_AK_bJets.eps",0, 1, 1E-3, 1,0);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets and subjets - Explicit JTA, SVClustering - CA&AK - bJetsGSP QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA"]                               = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_5_CA12_ready4comparison/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA"]         = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_5_CA12_ready4comparison/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK");   
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA");
 
  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (g#rightarrowb#bar{b} splitting)","","btagperfcomp_Pt200to400_FatJets_Subjets_CA_vs_AK_bJetsGSP.eps",0, 1, 1E-3, 1,0);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets and subjets - Explicit JTA, SVClustering - CA&AK - cJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA"]                               = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_5_CA12_ready4comparison/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA"]         = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_5_CA12_ready4comparison/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK");   
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA");
 
  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (c jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_CA_vs_AK_cJets.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets and subjets - Explicit JTA, SVClustering - CA&AK - udsJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA"]                               = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_5_CA12_ready4comparison/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA"]         = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_5_CA12_ready4comparison/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK");   
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA");
 
  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (uds jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_CA_vs_AK_udsJets.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets and subjets - Explicit JTA, SVClustering - CA&AK - gluonJets QCD
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA"]                               = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_5_CA12_ready4comparison/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA"]         = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_5_CA12_ready4comparison/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_4_AK12_QCDeventWeightFalse/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK");   
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA");
 
  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (gluon jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_CA_vs_AK_gluonJets.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //------------------------------------------
  // Fat jets and subjets - Explicit JTA, SVClustering - CA&AK - W bkg (JetMinMass75GeV)
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA"]                               = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/CA/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/CA/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA"]         = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/CA/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/CA/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"]                               = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK"]         = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"); 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA");

  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD, Hadronic W)","","btagperfcomp_Pt200to400_FatJets_Subjets_CA_vs_AK_HadronicW.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets and subjets - Explicit JTA, SVClustering - CA&AK - Z bkg
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA"]                               = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/CA/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA"]         = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/CA/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"); 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA");

  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD, Hadronic Z)","","btagperfcomp_Pt200to400_FatJets_Subjets_CA_vs_AK_HadronicZ.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //------------------------------------------
  // Fat jets and subjets - Explicit JTA, SVClustering - CA&AK - Top bkg
  //------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA"]                               = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/CA/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA"]         = getEfficiencyCurve("ROOT_files/jobs_5_CA12_ready4comparison/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/CA/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"]                               = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK"]         = getEfficiencyCurve("ROOT_files/jobs_3_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - AK"); 
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA, SV Clustering) - CA");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA");

  //------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{R=1.2, 200<p_{T}<400 GeV/c}{90<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (QCD, Hadronic Top)","","btagperfcomp_Pt200to400_FatJets_Subjets_CA_vs_AK_HadronicTop.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

}


void makePlots_AK_IVFCSVonly()
{
  // for multiple plots on the same canvas

  // vectors storing the order of legend entries
  std::vector< std::string> orderingPt200To400;
  // maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt200To400;

  //==========================================
  // Post-BTV-13-001 setup
  //==========================================

  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Fat jet CSV (*)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),200,400);
  graphsPt200To400["Subjet CSV (*)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),200,400);
   
  orderingPt200To400.push_back("Fat jet IVFCSV");
  orderingPt200To400.push_back("Subjet IVFCSV");
  orderingPt200To400.push_back("Fat jet CSV (*)");
  orderingPt200To400.push_back("Subjet CSV (*)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Inclusive QCD)","*BTV-13-001 CA12","btagperfcomp_Pt200to400_FatJets_Subjets_AK12.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - bJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  //graphsPt200To400["Fat jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","bJets"),200,400);
  //graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","bJets"),200,400);
   
  //orderingPt200To400.push_back("Fat jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat jet IVFCSV");
  //orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (b jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_bJets.eps",0, 1, 0, 1,0);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - bJetsGSP QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  //graphsPt200To400["Fat jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","bJetsGSP"),200,400);
  //graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","bJetsGSP"),200,400);
   
  //orderingPt200To400.push_back("Fat jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat jet IVFCSV");
  //orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (g#rightarrowb#bar{b} splitting)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_bJetsGSP.eps",0, 1, 0, 1,0);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - cJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  //graphsPt200To400["Fat jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","cJets"),200,400);
  //graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","cJets"),200,400);
   
  //orderingPt200To400.push_back("Fat jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat jet IVFCSV");
  //orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (c jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_cJets.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - udsJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  //graphsPt200To400["Fat jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","udsJets"),200,400);
  //graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","udsJets"),200,400);
   
  //orderingPt200To400.push_back("Fat jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat jet IVFCSV");
  //orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (uds jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_udsJets.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - gluonJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  //graphsPt200To400["Fat jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","gluonJets"),200,400);
  //graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","gluonJets"),200,400);
   
  //orderingPt200To400.push_back("Fat jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat jet IVFCSV");
  //orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (gluon jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_gluonJets.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - udsgJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsgJets"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsgJets"),200,400);
  //graphsPt200To400["Fat jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","udsgJets"),200,400);
  //graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","udsgJets"),200,400);
   
  //orderingPt200To400.push_back("Fat jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat jet IVFCSV");
  //orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (udsg jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_udsgJets.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12  - with BTV-13-001 overlay - W bkg (JetMinMass75GeV)
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  //graphsPt200To400["Fat jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/BTV13001_fatjet/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BTV13001_fatjet/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),200,400);
  //graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/BTV13001_fatjet/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BTV13001_fatjet/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),200,400);

  //orderingPt200To400.push_back("Fat jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat jet IVFCSV");
  //orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Hadronic W)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_HadronicW_JetMinMass75GeV.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - Z bkg
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  //graphsPt200To400["Fat jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BTV13001_fatjet/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),200,400);
  //graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BTV13001_fatjet/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),200,400);

  //orderingPt200To400.push_back("Fat jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat jet IVFCSV");
  //orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Hadronic Z)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_HadronicZ.eps",0, 1, 1E-2, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12  - with BTV-13-001 overlay - Top bkg 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  //graphsPt200To400["Fat jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BTV13001_fatjet/TprimeToTHinc_M-1000_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),200,400);
  //graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BTV13001_fatjet/TprimeToTHinc_M-1000_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),200,400);

  //orderingPt200To400.push_back("Fat jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat jet IVFCSV");
  //orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Hadronic top)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_HadronicTop.eps",0, 1, 1E-2, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();
}

void makePlots_AK_MatchedK4()
{
  // for multiple plots on the same canvas

  // vectors storing the order of legend entries
  std::vector< std::string> orderingPt200To400;
  // maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt200To400;
  
  //==========================================
  // Post-BTV-13-001 setup
  //==========================================

  // === Fatjet, Subjet, matched AK4 comparison ===

  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - Matched AK4 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times1)"]                              = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetIVFCSV"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times2)"]                           = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetMinIVFCSV"),200,400); 

  orderingPt200To400.push_back("Fat jet IVFCSV");
  orderingPt200To400.push_back("Subjet IVFCSV");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times1)");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times2)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Inclusive QCD)","#DeltaR(AK4 jets,fat jet)<0.9","btagperfcomp_Pt200to400_FatJets_Subjets_MatchedAK4_AK12.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - Matched AK4 - bJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times1)"]                              = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetIVFCSV"),getHistName("Pruned","StdJetIVFCSV","bJets"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times2)"]                           = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetMinIVFCSV"),getHistName("Pruned","StdJetMinIVFCSV","bJets"),200,400);
   
  orderingPt200To400.push_back("Fat jet IVFCSV");
  orderingPt200To400.push_back("Subjet IVFCSV");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times1)");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times2)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (b jets)","#DeltaR(AK4 jets,fat jet)<0.9","btagperfcomp_Pt200to400_FatJets_Subjets_MatchedAK4_AK12_bJets.eps",0, 1, 0, 1,0);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - Matched AK4 - bJetsGSP QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times1)"]                              = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetIVFCSV"),getHistName("Pruned","StdJetIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times2)"]                           = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetMinIVFCSV"),getHistName("Pruned","StdJetMinIVFCSV","bJetsGSP"),200,400);
     
  orderingPt200To400.push_back("Fat jet IVFCSV");
  orderingPt200To400.push_back("Subjet IVFCSV");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times1)");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times2)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (g#rightarrowb#bar{b} splitting)","#DeltaR(AK4 jets,fat jet)<0.9","btagperfcomp_Pt200to400_FatJets_Subjets_MatchedAK4_AK12_bJetsGSP.eps",0, 1, 0, 1,0);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - Matched AK4 - cJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times1)"]                              = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetIVFCSV"),getHistName("Pruned","StdJetIVFCSV","cJets"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times2)"]                           = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetMinIVFCSV"),getHistName("Pruned","StdJetMinIVFCSV","cJets"),200,400);

  orderingPt200To400.push_back("Fat jet IVFCSV");
  orderingPt200To400.push_back("Subjet IVFCSV");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times1)");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times2)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (c jets)","#DeltaR(AK4 jets,fat jet)<0.9","btagperfcomp_Pt200to400_FatJets_Subjets_MatchedAK4_AK12_cJets.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - Matched AK4 - udsJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times1)"]                              = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetIVFCSV"),getHistName("Pruned","StdJetIVFCSV","udsJets"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times2)"]                           = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetMinIVFCSV"),getHistName("Pruned","StdJetMinIVFCSV","udsJets"),200,400);
   
  orderingPt200To400.push_back("Fat jet IVFCSV");
  orderingPt200To400.push_back("Subjet IVFCSV");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times1)");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times2)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (uds jets)","#DeltaR(AK4 jets,fat jet)<0.9","btagperfcomp_Pt200to400_FatJets_Subjets_MatchedAK4_AK12_udsJets.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - Matched AK4 - gluonJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times1)"]                              = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetIVFCSV"),getHistName("Pruned","StdJetIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times2)"]                           = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetMinIVFCSV"),getHistName("Pruned","StdJetMinIVFCSV","gluonJets"),200,400);
   
  orderingPt200To400.push_back("Fat jet IVFCSV");
  orderingPt200To400.push_back("Subjet IVFCSV");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times1)");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times2)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (gluon jets)","#DeltaR(AK4 jets,fat jet)<0.9","btagperfcomp_Pt200to400_FatJets_Subjets_MatchedAK4_AK12_gluonJets.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - Matched AK4 - udsgJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsgJets"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsgJets"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times1)"]                              = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetIVFCSV"),getHistName("Pruned","StdJetIVFCSV","udsgJets"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times2)"]                           = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetMinIVFCSV"),getHistName("Pruned","StdJetMinIVFCSV","udsgJets"),200,400);
   
  orderingPt200To400.push_back("Fat jet IVFCSV");
  orderingPt200To400.push_back("Subjet IVFCSV");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times1)");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times2)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (udsg jets)","#DeltaR(AK4 jets,fat jet)<0.9","btagperfcomp_Pt200to400_FatJets_Subjets_MatchedAK4_AK12_udsgJets.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12  - Matched AK4 - W bkg (JetMinMass75GeV)
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times1)"]                              = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetIVFCSV"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times2)"]                           = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetMinIVFCSV"),200,400);

  orderingPt200To400.push_back("Fat jet IVFCSV");
  orderingPt200To400.push_back("Subjet IVFCSV");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times1)");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times2)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Hadronic W)","#DeltaR(AK4 jets,fat jet)<0.9","btagperfcomp_Pt200to400_FatJets_Subjets_MatchedAK4_AK12_HadronicW_JetMinMass75GeV.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12 - Matched AK4 - Z bkg 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times1)"]                              = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetIVFCSV"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times2)"]                           = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetMinIVFCSV"),200,400);

  orderingPt200To400.push_back("Fat jet IVFCSV");
  orderingPt200To400.push_back("Subjet IVFCSV");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times1)");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times2)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Hadronic Z)","#DeltaR(AK4 jets,fat jet)<0.9","btagperfcomp_Pt200to400_FatJets_Subjets_MatchedAK4_AK12_HadronicZ.eps",0, 1, 1E-2, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - Explicit JTA, SVClustering - AK12  - Matched AK4 - Top bkg 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat jet IVFCSV"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times1)"]                              = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetIVFCSV"),200,400);
  graphsPt200To400["Matched AK4 IVFCSV (#times2)"]                           = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","StdJetMinIVFCSV"),200,400);

  orderingPt200To400.push_back("Fat jet IVFCSV");
  orderingPt200To400.push_back("Subjet IVFCSV");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times1)");
  orderingPt200To400.push_back("Matched AK4 IVFCSV (#times2)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Hadronic top)","#DeltaR(AK4 jets,fat jet)<0.9","btagperfcomp_Pt200to400_FatJets_Subjets_MatchedAK4_AK12_HadronicTop.eps",0, 1, 1E-2, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();
} 
 
void makePlots_AK_IVFCSV_BTV13001_kT()
{
  // for multiple plots on the same canvas

  // vectors storing the order of legend entries
  std::vector< std::string> orderingPt200To400;
  // maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt200To400;

  //==========================================
  // Post-BTV-13-001 setup
  //==========================================

  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned, kT) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]  = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Fat Jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),200,400);
  graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),200,400);
   
  orderingPt200To400.push_back("Fat Jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Inclusive QCD)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_Kt_BTV13001.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned, kT) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - bJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]  = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","bJets"),200,400);
  graphsPt200To400["Fat Jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","bJets"),200,400);
  graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","bJets"),200,400);
   
  orderingPt200To400.push_back("Fat Jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (b jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_Kt_BTV13001_bJets.eps",0, 1, 0, 1,0);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned, kT) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - bJetsGSP QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]  = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","bJetsGSP"),200,400);
  graphsPt200To400["Fat Jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","bJetsGSP"),200,400);
  graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","bJetsGSP"),200,400);
   
  orderingPt200To400.push_back("Fat Jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (g#rightarrowb#bar{b} splitting)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_Kt_BTV13001_bJetsGSP.eps",0, 1, 0, 1,0);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned, kT) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - cJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]  = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","cJets"),200,400);
  graphsPt200To400["Fat Jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","cJets"),200,400);
  graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","cJets"),200,400);
   
  orderingPt200To400.push_back("Fat Jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (c jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_Kt_BTV13001_cJets.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned, kT) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - udsJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]  = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","udsJets"),200,400);
  graphsPt200To400["Fat Jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","udsJets"),200,400);
  graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","udsJets"),200,400);
   
  orderingPt200To400.push_back("Fat Jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (uds jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_Kt_BTV13001_udsJets.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned, kT) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - gluonJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]  = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","gluonJets"),200,400);
  graphsPt200To400["Fat Jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","gluonJets"),200,400);
  graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","gluonJets"),200,400);
   
  orderingPt200To400.push_back("Fat Jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (gluon jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_Kt_BTV13001_gluonJets.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned, kT) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - udsgJets QCD 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),getHistName("Pruned","JetIVFCSV","udsgJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),getHistName("Pruned","SubJetMinIVFCSV","udsgJets"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]  = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files_AK12/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),getHistName("Kt","SubJetMinIVFCSV","udsgJets"),200,400);
  graphsPt200To400["Fat Jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),getHistName("Pruned","JetCSV","udsgJets"),200,400);
  graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),getHistName("Pruned","SubJetMinCSV","udsgJets"),200,400);
   
  orderingPt200To400.push_back("Fat Jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (udsg jets)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_Kt_BTV13001_udsgJets.eps",0, 1, 1E-4, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned, kT) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - W bkg (JetMinMass75GeV)
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]  = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Fat Jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/BTV13001_fatjet/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BTV13001_fatjet/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),200,400);
  graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files/jobs_6_forApproval/BTV13001_fatjet/JetMinMass75GeV/RadionToHH_4b_M-600_HiggsTagging_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BTV13001_fatjet/JetMinMass75GeV/TprimeTprimeToBWBWinc_M-800_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),200,400);

  orderingPt200To400.push_back("Fat Jet CSV (CA12, BTV-13-001)");  
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<150 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Hadronic W)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_Kt_BTV13001_HadronicW_JetMinMass75GeV.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned, kT) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - Z bkg
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]  = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Fat Jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BTV13001_fatjet/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),200,400);
  graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BTV13001_fatjet/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),200,400);

  orderingPt200To400.push_back("Fat Jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Hadronic Z)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_Kt_BTV13001_HadronicZ.eps",0, 1, 1E-2, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();


  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned, kT) - Explicit JTA, SVClustering - AK12 - with BTV-13-001 overlay - Top bkg 
  //--------------------------------------------------------------------------------------------
  graphsPt200To400["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","JetIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]  = getEfficiencyCurve("ROOT_files_AK12/RadionToHH_4b_M-600_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3_R1p2.root",getHistName("Kt","SubJetMinIVFCSV"),200,400);
  graphsPt200To400["Fat Jet CSV (CA12, BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BTV13001_fatjet/TprimeToTHinc_M-1000_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),200,400);
  graphsPt200To400["Subjet CSV (CA12, Pruned, BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files/jobs_6_forApproval/BTV13001_fatjet/TprimeToTHinc_M-1000_HiggsTagging_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),200,400);

  orderingPt200To400.push_back("Fat Jet CSV (CA12, BTV-13-001)");
  orderingPt200To400.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt200To400.push_back("Subjet CSV (CA12, Pruned, BTV-13-001)");
  orderingPt200To400.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt200To400.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");

  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{AK R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Hadronic top)","","btagperfcomp_Pt200to400_FatJets_Subjets_AK12_Kt_BTV13001_HadronicTop.eps",0, 1, 1E-2, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

}

void makePlots()
{
  writeExtraText = true;
  extraText  = "Simulation Preliminary";
  lumi_8TeV = "";
  relPosX = 0.15;

  // for multiple plots on the same canvas

  // vectors storing the order of legend entries
  std::vector< std::string> orderingPt200To400;
  // maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt200To400;

  //--------------------------------------------------------------------------------------------
  // Fat jets and subjets (pruned) - CA12 -  BTV-13-001 
  //--------------------------------------------------------------------------------------------
   graphsPt200To400["Fat jet CSV (BTV-13-001)"]                      = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","JetCSV"),200,400);
  graphsPt200To400["Subjet CSV (BTV-13-001)"]               = getEfficiencyCurve("ROOT_files_CA12/RadionToHH_4b_M-600_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root","ROOT_files_CA12/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3_R1p2.root",getHistName("Pruned","SubJetMinCSV"),200,400);
   
  orderingPt200To400.push_back("Fat jet CSV (BTV-13-001)");
  orderingPt200To400.push_back("Subjet CSV (BTV-13-001)");
 
  //-------------------------------------------------------------------------------------------
  
  plotEfficiencyCurves(graphsPt200To400,orderingPt200To400,"#splitline{CA R=1.2, 200<p_{T}<400 GeV/c}{75<m_{pruned}<135 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Inclusive QCD)","","btagperfcomp_Pt200to400_FatJets_Subjets_BTV13001.eps",0, 1, 1E-3, 1,1);
 
  graphsPt200To400.clear();
  orderingPt200To400.clear();

  //===========================

  //makePlots_AK_AllsubjetAllsetting();
  //makePlots_CA_AK_comparison();
  makePlots_AK_IVFCSVonly();
  makePlots_AK_MatchedK4();
  //makePlots_AK_IVFCSV_BTV13001_kT();
}
