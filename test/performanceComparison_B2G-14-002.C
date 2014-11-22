#include <iostream>
#include <map>
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
#include "Rtypes.h"
#include "exoStyle.C"
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

  file_S1->Close();
  file_B1->Close();

  return g_eff_1;
}


TGraph* getEfficiencyCurve(const string& fFileS1, const string& fFileB1, const string& fPlot, const double fXMin, const double fXMax)
{
  return getEfficiencyCurve(fFileS1, fFileB1, fPlot, fPlot, fXMin, fXMax);
}


TGraph* getEfficiencyCurve2D(const string& fFileS1, const string& fFileB1,const string& fPlot1,const string& fPlot2)
{
  //get files and histograms
  TFile *file_S1  = new TFile(fFileS1.c_str());
  TFile *file_B1 = new TFile(fFileB1.c_str());

  TH2 *h2_S_1 = (TH2*)file_S1->Get(fPlot1.c_str());
  TH2 *h2_B_1 = (TH2*)file_B1->Get(fPlot2.c_str());

  int nBinsX = h2_S_1->GetXaxis()->GetNbins();
  int overflowX = nBinsX+1;
  int nBinsY = h2_S_1->GetYaxis()->GetNbins();
  int overflowY = nBinsY+1;

  std::cout << "Producing optimal performance curve from 2D cuts..." << std::endl;

  //total jet count for denominator of efficiency calculation
  double denom_S_1 = h2_S_1->Integral(0,overflowX,0,overflowY);
  double denom_B_1 = h2_B_1->Integral(0,overflowX,0,overflowY);

  TH2F *h2_temp = new TH2F("h2_temp","h2_temp",1000,0.,1.,1000,0.,1.);

  int entryCounter = 0;

  for(int ix = 1; ix<=overflowX; ++ix)
  {
    for(int iy = 1; iy<=overflowY; ++iy)
    {
      if(entryCounter%((overflowX*overflowY)/10) == 0) std::cout << (float(entryCounter)/float(overflowX*overflowY)*100.) << "% completed" << std::endl;

      double num_S = h2_S_1->Integral(ix,overflowX,iy,overflowY);
      double num_B = h2_B_1->Integral(ix,overflowX,iy,overflowY);

      h2_temp->Fill((num_S/denom_S_1),(num_B/denom_B_1));

      ++entryCounter;
    }
  }

  float eff_x[1000] = {0.};
  float eff_y[1000] = {0.};
  int point = 0;

  for(int iy=1; iy<=h2_temp->GetYaxis()->GetNbins(); ++iy)
  {
    float eff = -1.;
    for(int ix=h2_temp->GetXaxis()->GetNbins(); ix>0; --ix)
    {
      if(h2_temp->GetBinContent(ix,iy)>0.5)
      {
        eff = h2_temp->GetXaxis()->GetBinCenter(ix);
        break;
      }
    }
    if(eff>0.)
    {
      eff_x[point] = eff;
      eff_y[point] = h2_temp->GetYaxis()->GetBinCenter(iy);
      //std::cout << "eff_x: " << eff << std::endl;
      //std::cout << "eff_y: " << h2_temp->GetYaxis()->GetBinCenter(iy) << std::endl;
      ++point;
    }
  }

  TGraph *g_eff_1 = new TGraph(point,eff_x,eff_y);

  std::cout << "Done!" << std::endl;

  //h2_temp->SaveAs("temp.root");

  return g_eff_1;
}


TGraph* getEfficiencyCurve2D(const string& fFileS1, const string& fFileB1, const string& fPlot)
{
  return getEfficiencyCurve2D(fFileS1, fFileB1, fPlot, fPlot);
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


std::string getHistName(const std::string & grooming = "Pruned", const std::string & algo = "CSV", const std::string & flavor = "", const std::string & postfix = "", const std::string & jet = "")
{
  return "jetAnalyzer"+jet+"FatJets_"+grooming+"Subjets" + (flavor != "" ? "_" + flavor : "") + (postfix != "" ? "_" + postfix : "") + "/h2_JetPt_"+algo+"_BosonMatched_JetMass";
}


std::string getHistName2D(const std::string & grooming, const std::string & algo, const std::string & ptbin, const std::string & flavor = "", const std::string & postfix = "")
{
  return "jetAnalyzerFatJets_"+grooming+"Subjets" + (flavor != "" ? "_" + flavor : "") + (postfix != "" ? "_" + postfix : "") + "/h2_SubJet1"+algo+"_SubJet2"+algo+"_BosonMatched_JetMass_Pt"+ptbin;
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

  TLegend *legend = new TLegend(.16,.75-float(ordering.size())*(0.05),.36,.75);
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

  legend->Draw();
  TLatex l1;
  l1.SetTextAlign(13);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.14+0.03,0.90, fTitle.c_str());

  //l1.SetTextAlign(12);
  //l1.SetTextSize(0.045);
  //l1.SetTextFont(62);
  //l1.DrawLatex(0.14,0.96, "CMS Simulation Preliminary, #sqrt{s} = 8 TeV");

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

  l1.SetTextFont(42);
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.48,0.20, fExtraInfo.c_str());

  if (fLogy) c->SetLogy();

  c->SaveAs(fOutputFile.c_str());

  graphs.clear();

  delete c;
  delete legend;
  delete bkg;
}

void makePlots()
{
  writeExtraText = false;       // if extra text
  extraText  = "Simulation Preliminary";  // default extra text is "Preliminary"
  relPosX    = 0.15;
  lumi_sqrtS = "8 TeV"; // default is "19.7 fb^{-1}"

  // for multiple plots on the same canvas

  // vectors storing the order of legend entries
  std::vector< std::string> orderingPt150To300;
  std::vector< std::string> orderingPt300To500;

  // maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt150To300;
  std::map< std::string,TGraph* > graphsPt300To500;

  //==========================================
  // BTV-13-001 setup
  //==========================================

  //-------------------------------------------------
  // Fat jets, subjets and standard jets -- H->bb vs Inclusive QCD
  //-------------------------------------------------

  graphsPt150To300["R=0.8 fat jet CSV"]       = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R0p8_BTV-13-001_PATTuple_v3.root", "QCDPythia6_HiggsTagging_R0p8_BTV-13-001_PATTuple_v3.root", getHistName("Pruned","JetCSV"),            150,300);
  graphsPt150To300["R=0.8 subjet CSV"]        = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R0p8_BTV-13-001_PATTuple_v3.root", "QCDPythia6_HiggsTagging_R0p8_BTV-13-001_PATTuple_v3.root", getHistName("Pruned","SubJetMinCSV"),      150,300);
  graphsPt150To300["R=1.5 fat jet CSV"]       = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "QCDPythia6_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","JetCSV"),          150,300);
  graphsPt150To300["R=1.5 subjet CSV"]        = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "QCDPythia6_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","SubJetMinCSV"),    150,300);
  graphsPt150To300["Matched AK5 CSV (#geq1)"] = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "QCDPythia6_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","StdJetIVFCSV"),    150,300);
  graphsPt150To300["Matched AK5 CSV (#geq2)"] = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "QCDPythia6_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","StdJetMinIVFCSV"), 150,300);

  orderingPt150To300.push_back("R=0.8 fat jet CSV");
  orderingPt150To300.push_back("R=0.8 subjet CSV");
  orderingPt150To300.push_back("R=1.5 fat jet CSV");
  orderingPt150To300.push_back("R=1.5 subjet CSV");
  orderingPt150To300.push_back("Matched AK5 CSV (#geq1)");
  orderingPt150To300.push_back("Matched AK5 CSV (#geq2)");
  //-------------------------------------------------
  graphsPt300To500["R=0.8 fat jet CSV"]       = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-1000_HiggsTagging_R0p8_BTV-13-001_PATTuple_v3.root", "QCDPythia6_HiggsTagging_R0p8_BTV-13-001_PATTuple_v3.root", getHistName("Pruned","JetCSV"),            300,500);
  graphsPt300To500["R=0.8 subjet CSV"]        = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-1000_HiggsTagging_R0p8_BTV-13-001_PATTuple_v3.root", "QCDPythia6_HiggsTagging_R0p8_BTV-13-001_PATTuple_v3.root", getHistName("Pruned","SubJetMinCSV"),      300,500);
  graphsPt300To500["R=1.5 fat jet CSV"]       = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-1000_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "QCDPythia6_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","JetCSV"),          300,500);
  graphsPt300To500["R=1.5 subjet CSV"]        = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-1000_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "QCDPythia6_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","SubJetMinCSV"),    300,500);
  graphsPt300To500["Matched AK5 CSV (#geq1)"] = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-1000_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "QCDPythia6_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","StdJetIVFCSV"),    300,500);
  graphsPt300To500["Matched AK5 CSV (#geq2)"] = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-1000_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "QCDPythia6_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","StdJetMinIVFCSV"), 300,500);

  orderingPt300To500.push_back("R=0.8 fat jet CSV");
  orderingPt300To500.push_back("R=0.8 subjet CSV");
  orderingPt300To500.push_back("R=1.5 fat jet CSV");
  orderingPt300To500.push_back("R=1.5 subjet CSV");
  orderingPt300To500.push_back("Matched AK5 CSV (#geq1)");
  orderingPt300To500.push_back("Matched AK5 CSV (#geq2)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt150To300,orderingPt150To300, "#splitline{CA, 150<p_{T}<300 GeV/c, |#eta|<2.4}{m>60 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Inclusive QCD)", "#DeltaR(AK5,CA15)<1.1", "btagperfcomp_Pt150to300_B2G-14-002_BTV-13-001.eps", 0, 1, 1E-3, 1, 1);
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, "#splitline{CA, 300<p_{T}<500 GeV/c, |#eta|<2.4}{m>60 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Inclusive QCD)", "#DeltaR(AK5,CA15)<1.1", "btagperfcomp_Pt300to500_B2G-14-002_BTV-13-001.eps", 0, 1, 1E-3, 1, 1);

  graphsPt150To300.clear();
  graphsPt300To500.clear();

  orderingPt150To300.clear();
  orderingPt300To500.clear();

  //-------------------------------------------------
  // Fat jets, subjets and standard jets -- H->bb vs Inclusive QCD
  //-------------------------------------------------

  graphsPt150To300["R=0.8 fat jet CSV"]         = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R0p8_BTV-13-001_PATTuple_v3.root", "TprimeTprimeToTHTHinc_M-600_HiggsTagging_TopBkg_R0p8_BTV-13-001_PATTuple_v3.root", getHistName("Pruned","JetCSV"),            150,300);
  graphsPt150To300["R=0.8 subjet CSV"]          = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R0p8_BTV-13-001_PATTuple_v3.root", "TprimeTprimeToTHTHinc_M-600_HiggsTagging_TopBkg_R0p8_BTV-13-001_PATTuple_v3.root", getHistName("Pruned","SubJetMinCSV"),      150,300);
  graphsPt150To300["R=1.5 fat jet CSV"]         = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "TprimeTprimeToTHTHinc_M-600_HiggsTagging_TopBkg_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","JetCSV"),          150,300);
  graphsPt150To300["R=1.5 subjet CSV"]          = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "TprimeTprimeToTHTHinc_M-600_HiggsTagging_TopBkg_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","SubJetMinCSV"),    150,300);
  graphsPt150To300["Matched AK5 CSV (#geq1)"] = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "TprimeTprimeToTHTHinc_M-600_HiggsTagging_TopBkg_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","StdJetIVFCSV"),    150,300);
  graphsPt150To300["Matched AK5 CSV (#geq2)"] = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-600_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "TprimeTprimeToTHTHinc_M-600_HiggsTagging_TopBkg_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","StdJetMinIVFCSV"), 150,300);

  orderingPt150To300.push_back("R=0.8 fat jet CSV");
  orderingPt150To300.push_back("R=0.8 subjet CSV");
  orderingPt150To300.push_back("R=1.5 fat jet CSV");
  orderingPt150To300.push_back("R=1.5 subjet CSV");
  orderingPt150To300.push_back("Matched AK5 CSV (#geq1)");
  orderingPt150To300.push_back("Matched AK5 CSV (#geq2)");
  //-------------------------------------------------
  graphsPt300To500["R=0.8 fat jet CSV"]         = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-1000_HiggsTagging_R0p8_BTV-13-001_PATTuple_v3.root", "TprimeTprimeToTHTHinc_M-1000_HiggsTagging_TopBkg_R0p8_BTV-13-001_PATTuple_v3.root", getHistName("Pruned","JetCSV"),            300,500);
  graphsPt300To500["R=0.8 subjet CSV"]          = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-1000_HiggsTagging_R0p8_BTV-13-001_PATTuple_v3.root", "TprimeTprimeToTHTHinc_M-1000_HiggsTagging_TopBkg_R0p8_BTV-13-001_PATTuple_v3.root", getHistName("Pruned","SubJetMinCSV"),      300,500);
  graphsPt300To500["R=1.5 fat jet CSV"]         = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-1000_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "TprimeTprimeToTHTHinc_M-1000_HiggsTagging_TopBkg_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","JetCSV"),          300,500);
  graphsPt300To500["R=1.5 subjet CSV"]          = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-1000_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "TprimeTprimeToTHTHinc_M-1000_HiggsTagging_TopBkg_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","SubJetMinCSV"),    300,500);
  graphsPt300To500["Matched AK5 CSV (#geq1)"] = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-1000_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "TprimeTprimeToTHTHinc_M-1000_HiggsTagging_TopBkg_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","StdJetIVFCSV"),    300,500);
  graphsPt300To500["Matched AK5 CSV (#geq2)"] = getEfficiencyCurve("TprimeTprimeToTHTHinc_M-1000_HiggsTagging_R1p5_BTV-13-001_PATTuple_v3.root", "TprimeTprimeToTHTHinc_M-1000_HiggsTagging_TopBkg_R1p5_BTV-13-001_PATTuple_v3.root", getHistName("Filtered","StdJetMinIVFCSV"), 300,500);

  orderingPt300To500.push_back("R=0.8 fat jet CSV");
  orderingPt300To500.push_back("R=0.8 subjet CSV");
  orderingPt300To500.push_back("R=1.5 fat jet CSV");
  orderingPt300To500.push_back("R=1.5 subjet CSV");
  orderingPt300To500.push_back("Matched AK5 CSV (#geq1)");
  orderingPt300To500.push_back("Matched AK5 CSV (#geq2)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt150To300,orderingPt150To300, "#splitline{CA, 150<p_{T}<300 GeV/c, |#eta|<2.4}{m>60 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Hadronic top)", "#DeltaR(AK5,CA15)<1.1", "btagperfcomp_Pt150to300_B2G-14-002_BTV-13-001_Hadronic_top.eps", 0, 1, 1E-2, 1, 1);
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, "#splitline{CA, 300<p_{T}<500 GeV/c, |#eta|<2.4}{m>60 GeV/c^{2}}", "Tagging efficiency (H#rightarrowb#bar{b})", "Misid. probability (Hadronic top)", "#DeltaR(AK5,CA15)<1.1", "btagperfcomp_Pt300to500_B2G-14-002_BTV-13-001_Hadronic_top.eps", 0, 1, 1E-2, 1, 1);

  graphsPt150To300.clear();
  graphsPt300To500.clear();

  orderingPt150To300.clear();
  orderingPt300To500.clear();
  //-------------------------------------------------
}
