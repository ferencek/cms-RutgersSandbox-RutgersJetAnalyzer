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


void makePlotsQCD(const string & dir = "ROOT_files", const string & algo = "AK", const double Ymin = 1e-3, const Int_t Logy = 1, const string & flavor = "", const string & ext = "eps")
{
  // for multiple plots on the same canvas

  // vectors storing the order of legend entries
  std::vector< std::string> orderingPt300To500;
  std::vector< std::string> orderingPt700ToInf;
  // maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt300To500;
  std::map< std::string,TGraph* > graphsPt700ToInf;
  
  //==========================================
  // Post-BTV-13-001 setup
  //==========================================
  // Fat jets
  //-------------------------------------------------
  graphsPt300To500["Fat Jet CSV"]              = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","JetCSV"),    getHistName("Pruned","JetCSV",flavor),    300,500);
  graphsPt300To500["Fat Jet IVFCSV"]           = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 300,500);
  graphsPt300To500["Fat Jet JP"]               = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","JetJP"),     getHistName("Pruned","JetJP",flavor),     300,500);
  graphsPt300To500["Fat Jet JBP"]              = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","JetJBP"),    getHistName("Pruned","JetJBP",flavor),    300,500);
  
  orderingPt300To500.push_back("Fat Jet CSV");
  orderingPt300To500.push_back("Fat Jet IVFCSV");
  orderingPt300To500.push_back("Fat Jet JP");
  orderingPt300To500.push_back("Fat Jet JBP");
  //-------------------------------------------------
  graphsPt700ToInf["Fat Jet CSV"]              = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","JetCSV"),    getHistName("Pruned","JetCSV",flavor),    700,1100);
  graphsPt700ToInf["Fat Jet IVFCSV"]           = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Fat Jet JP"]               = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","JetJP"),     getHistName("Pruned","JetJP",flavor),     700,1100);
  graphsPt700ToInf["Fat Jet JBP"]              = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","JetJBP"),    getHistName("Pruned","JetJBP",flavor),    700,1100);

  orderingPt700ToInf.push_back("Fat Jet CSV");
  orderingPt700ToInf.push_back("Fat Jet IVFCSV");
  orderingPt700ToInf.push_back("Fat Jet JP");
  orderingPt700ToInf.push_back("Fat Jet JBP");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_FatJets_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_FatJets_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //-------------------------------------------------
  // Fat jets: IVFCSV (Cone vs Explicit JTA)
  //-------------------------------------------------
  graphsPt300To500["Fat Jet IVFCSV"]                = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(),                          (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),                          getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 300,500);
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 300,500);
  
  orderingPt300To500.push_back("Fat Jet IVFCSV");
  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA)");
  //-------------------------------------------------
  graphsPt700ToInf["Fat Jet IVFCSV"]                = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(),                          (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),                          getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 700,1100);
  
  orderingPt700ToInf.push_back("Fat Jet IVFCSV");
  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA)");
  //------------------------------------------------- 
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_FatJets_IVFCSV_Cone_vs_ExplJTA_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_FatJets_IVFCSV_Cone_vs_ExplJTA_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //-------------------------------------------------
  // Fat jets: IVFCSV Explicit JTA comparison
  //-------------------------------------------------
  graphsPt300To500["Fat Jet IVFCSV"]                            = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(),                          (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),                          getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 300,500);
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA)"]             = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 300,500);
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA signal only)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),                          getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 300,500);
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA bkg only)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(),                          (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 300,500);
  
  orderingPt300To500.push_back("Fat Jet IVFCSV");
  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA signal only)");
  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA bkg only)");
  //-------------------------------------------------
  graphsPt700ToInf["Fat Jet IVFCSV"]                            = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(),                          (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),                          getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA)"]             = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA signal only)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),                          getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA bkg only)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(),                          (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"), getHistName("Pruned","JetIVFCSV",flavor), 700,1100);
  
  orderingPt700ToInf.push_back("Fat Jet IVFCSV");
  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA signal only)");
  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA bkg only)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_FatJets_IVFCSV_ExplJTA_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_FatJets_IVFCSV_ExplJTA_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //-------------------------------------------------
  // Subjets: CSV
  //-------------------------------------------------
  graphsPt300To500["Subjet CSV (Pruned)"]         = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinCSV"),         getHistName("Pruned","SubJetMinCSV",flavor),         300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Filtered","SubJetMinCSV"),       getHistName("Filtered","SubJetMinCSV",flavor),       300,500);
  graphsPt300To500["Subjet CSV (k_{T})"]          = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinCSV"),             getHistName("Kt","SubJetMinCSV",flavor),             300,500);
  graphsPt300To500["Subjet CSV (k_{T}+Filtered)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("KtBDRSFiltered","SubJetMinCSV"), getHistName("KtBDRSFiltered","SubJetMinCSV",flavor), 300,500);

  orderingPt300To500.push_back("Subjet CSV (Pruned)");
  orderingPt300To500.push_back("Subjet CSV (Filtered)");
  orderingPt300To500.push_back("Subjet CSV (k_{T})");
  orderingPt300To500.push_back("Subjet CSV (k_{T}+Filtered)");
  //-------------------------------------------------
  graphsPt700ToInf["Subjet CSV (Pruned)"]         = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinCSV"),         getHistName("Pruned","SubJetMinCSV",flavor),         700,1100);
  graphsPt700ToInf["Subjet CSV (Filtered)"]       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Filtered","SubJetMinCSV"),       getHistName("Filtered","SubJetMinCSV",flavor),       700,1100);
  graphsPt700ToInf["Subjet CSV (k_{T})"]          = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinCSV"),             getHistName("Kt","SubJetMinCSV",flavor),             700,1100);
  graphsPt700ToInf["Subjet CSV (k_{T}+Filtered)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("KtBDRSFiltered","SubJetMinCSV"), getHistName("KtBDRSFiltered","SubJetMinCSV",flavor), 700,1100);

  orderingPt700ToInf.push_back("Subjet CSV (Pruned)");
  orderingPt700ToInf.push_back("Subjet CSV (Filtered)");
  orderingPt700ToInf.push_back("Subjet CSV (k_{T})");
  orderingPt700ToInf.push_back("Subjet CSV (k_{T}+Filtered)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_Subjets_CSV_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_Subjets_CSV_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //-------------------------------------------------
  // Subjets: IVFCSV
  //-------------------------------------------------
  graphsPt300To500["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"),         getHistName("Pruned","SubJetMinIVFCSV",flavor),         300,500);
  graphsPt300To500["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Filtered","SubJetMinIVFCSV"),       getHistName("Filtered","SubJetMinIVFCSV",flavor),       300,500);
  graphsPt300To500["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinIVFCSV"),             getHistName("Kt","SubJetMinIVFCSV",flavor),             300,500);
  graphsPt300To500["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("KtBDRSFiltered","SubJetMinIVFCSV"), getHistName("KtBDRSFiltered","SubJetMinIVFCSV",flavor), 300,500);

  orderingPt300To500.push_back("Subjet IVFCSV (Pruned)");
  orderingPt300To500.push_back("Subjet IVFCSV (Filtered)");
  orderingPt300To500.push_back("Subjet IVFCSV (k_{T})");
  orderingPt300To500.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //-------------------------------------------------
  graphsPt700ToInf["Subjet IVFCSV (Pruned)"]         = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"),         getHistName("Pruned","SubJetMinIVFCSV",flavor),         700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Filtered)"]       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Filtered","SubJetMinIVFCSV"),       getHistName("Filtered","SubJetMinIVFCSV",flavor),       700,1100);
  graphsPt700ToInf["Subjet IVFCSV (k_{T})"]          = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinIVFCSV"),             getHistName("Kt","SubJetMinIVFCSV",flavor),             700,1100);
  graphsPt700ToInf["Subjet IVFCSV (k_{T}+Filtered)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("KtBDRSFiltered","SubJetMinIVFCSV"), getHistName("KtBDRSFiltered","SubJetMinIVFCSV",flavor), 700,1100);

  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Filtered)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T})");
  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T}+Filtered)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_Subjets_IVFCSV_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_Subjets_IVFCSV_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
  // Subjets: IVFCSV (with Explicit JTA and SV Clustering)
  //-------------------------------------------------
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]         = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"),         getHistName("Pruned","SubJetMinIVFCSV",flavor),         300,500);
  graphsPt300To500["Subjet IVFCSV (Filtered, Explicit JTA, SV Clustering)"]       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Filtered","SubJetMinIVFCSV"),       getHistName("Filtered","SubJetMinIVFCSV",flavor),       300,500);
  graphsPt300To500["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]          = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinIVFCSV"),             getHistName("Kt","SubJetMinIVFCSV",flavor),             300,500);
  graphsPt300To500["Subjet IVFCSV (k_{T}+Filtered, Explicit JTA, SV Clustering)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("KtBDRSFiltered","SubJetMinIVFCSV"), getHistName("KtBDRSFiltered","SubJetMinIVFCSV",flavor), 300,500);

  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (Filtered, Explicit JTA, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (k_{T}+Filtered, Explicit JTA, SV Clustering)");
  //-------------------------------------------------
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]         = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"),         getHistName("Pruned","SubJetMinIVFCSV",flavor),         700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Filtered, Explicit JTA, SV Clustering)"]       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Filtered","SubJetMinIVFCSV"),       getHistName("Filtered","SubJetMinIVFCSV",flavor),       700,1100);
  graphsPt700ToInf["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]          = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinIVFCSV"),             getHistName("Kt","SubJetMinIVFCSV",flavor),             700,1100);
  graphsPt700ToInf["Subjet IVFCSV (k_{T}+Filtered, Explicit JTA, SV Clustering)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("KtBDRSFiltered","SubJetMinIVFCSV"), getHistName("KtBDRSFiltered","SubJetMinIVFCSV",flavor), 700,1100);

  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Filtered, Explicit JTA, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T}+Filtered, Explicit JTA, SV Clustering)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_Subjets_ExplJTA_SV_IVFCSV_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_Subjets_ExplJTA_SV_IVFCSV_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //-------------------------------------------------
  // Subjets: JP
  //-------------------------------------------------
  graphsPt300To500["Subjet JP (Pruned)"]         = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinJP"),         getHistName("Pruned","SubJetMinJP",flavor),         300,500);
  graphsPt300To500["Subjet JP (Filtered)"]       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Filtered","SubJetMinJP"),       getHistName("Filtered","SubJetMinJP",flavor),       300,500);
  graphsPt300To500["Subjet JP (k_{T})"]          = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinJP"),             getHistName("Kt","SubJetMinJP",flavor),             300,500);
  graphsPt300To500["Subjet JP (k_{T}+Filtered)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("KtBDRSFiltered","SubJetMinJP"), getHistName("KtBDRSFiltered","SubJetMinJP",flavor), 300,500);

  orderingPt300To500.push_back("Subjet JP (Pruned)");
  orderingPt300To500.push_back("Subjet JP (Filtered)");
  orderingPt300To500.push_back("Subjet JP (k_{T})");
  orderingPt300To500.push_back("Subjet JP (k_{T}+Filtered)");
  //-------------------------------------------------
  graphsPt700ToInf["Subjet JP (Pruned)"]         = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinJP"),         getHistName("Pruned","SubJetMinJP",flavor),         700,1100);
  graphsPt700ToInf["Subjet JP (Filtered)"]       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Filtered","SubJetMinJP"),       getHistName("Filtered","SubJetMinJP",flavor),       700,1100);
  graphsPt700ToInf["Subjet JP (k_{T})"]          = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinJP"),             getHistName("Kt","SubJetMinJP",flavor),             700,1100);
  graphsPt700ToInf["Subjet JP (k_{T}+Filtered)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("KtBDRSFiltered","SubJetMinJP"), getHistName("KtBDRSFiltered","SubJetMinJP",flavor), 700,1100);

  orderingPt700ToInf.push_back("Subjet JP (Pruned)");
  orderingPt700ToInf.push_back("Subjet JP (Filtered)");
  orderingPt700ToInf.push_back("Subjet JP (k_{T})");
  orderingPt700ToInf.push_back("Subjet JP (k_{T}+Filtered)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_Subjets_JP_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_Subjets_JP_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //-------------------------------------------------
  // Subjets: JBP
  //-------------------------------------------------
  graphsPt300To500["Subjet JBP (Pruned)"]         = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinJBP"),         getHistName("Pruned","SubJetMinJBP",flavor),         300,500);
  graphsPt300To500["Subjet JBP (Filtered)"]       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Filtered","SubJetMinJBP"),       getHistName("Filtered","SubJetMinJBP",flavor),       300,500);
  graphsPt300To500["Subjet JBP (k_{T})"]          = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinJBP"),             getHistName("Kt","SubJetMinJBP",flavor),             300,500);
  graphsPt300To500["Subjet JBP (k_{T}+Filtered)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("KtBDRSFiltered","SubJetMinJBP"), getHistName("KtBDRSFiltered","SubJetMinJBP",flavor), 300,500);

  orderingPt300To500.push_back("Subjet JBP (Pruned)");
  orderingPt300To500.push_back("Subjet JBP (Filtered)");
  orderingPt300To500.push_back("Subjet JBP (k_{T})");
  orderingPt300To500.push_back("Subjet JBP (k_{T}+Filtered)");
  //-------------------------------------------------
  graphsPt700ToInf["Subjet JBP (Pruned)"]         = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),getHistName("Pruned","SubJetMinJBP"),         getHistName("Pruned","SubJetMinJBP",flavor),         700,1100);
  graphsPt700ToInf["Subjet JBP (Filtered)"]       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),getHistName("Filtered","SubJetMinJBP"),       getHistName("Filtered","SubJetMinJBP",flavor),       700,1100);
  graphsPt700ToInf["Subjet JBP (k_{T})"]          = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),getHistName("Kt","SubJetMinJBP"),             getHistName("Kt","SubJetMinJBP",flavor),             700,1100);
  graphsPt700ToInf["Subjet JBP (k_{T}+Filtered)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),getHistName("KtBDRSFiltered","SubJetMinJBP"), getHistName("KtBDRSFiltered","SubJetMinJBP",flavor), 700,1100);

  orderingPt700ToInf.push_back("Subjet JBP (Pruned)");
  orderingPt700ToInf.push_back("Subjet JBP (Filtered)");
  orderingPt700ToInf.push_back("Subjet JBP (k_{T})");
  orderingPt700ToInf.push_back("Subjet JBP (k_{T}+Filtered)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_Subjets_JBP_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_Subjets_JBP_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
  // Subjets: Pruned comparison
  //-------------------------------------------------
  graphsPt300To500["Subjet CSV (Pruned)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinCSV"),    getHistName("Pruned","SubJetMinCSV",flavor),    300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet JP (Pruned)"]     = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinJP"),     getHistName("Pruned","SubJetMinJP",flavor),     300,500);
  graphsPt300To500["Subjet JBP (Pruned)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinJBP"),    getHistName("Pruned","SubJetMinJBP",flavor),    300,500);
  
  orderingPt300To500.push_back("Subjet CSV (Pruned)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned)");
  orderingPt300To500.push_back("Subjet JP (Pruned)");
  orderingPt300To500.push_back("Subjet JBP (Pruned)");
  //-------------------------------------------------
  graphsPt700ToInf["Subjet CSV (Pruned)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinCSV"),    getHistName("Pruned","SubJetMinCSV",flavor),    700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet JP (Pruned)"]     = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinJP"),     getHistName("Pruned","SubJetMinJP",flavor),     700,1100);
  graphsPt700ToInf["Subjet JBP (Pruned)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinJBP"),    getHistName("Pruned","SubJetMinJBP",flavor),    700,1100);
  
  orderingPt700ToInf.push_back("Subjet CSV (Pruned)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned)");
  orderingPt700ToInf.push_back("Subjet JP (Pruned)");
  orderingPt700ToInf.push_back("Subjet JBP (Pruned)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(),"",("btagperfcomp_Pt300to500_Subjets_Pruned_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(),0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(),"",("btagperfcomp_Pt700toInf_Subjets_Pruned_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(),0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //-------------------------------------------------
  // Subjets: Pruned comparison (with Explicit JTA and SV Clustering)
  //-------------------------------------------------
  graphsPt300To500["Subjet CSV (Pruned, Explicit JTA, SV Clustering)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinCSV"),    getHistName("Pruned","SubJetMinCSV",flavor),    300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet JP (Pruned, Explicit JTA, SV Clustering)"]     = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinJP"),     getHistName("Pruned","SubJetMinJP",flavor),     300,500);
  graphsPt300To500["Subjet JBP (Pruned, Explicit JTA, SV Clustering)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinJBP"),    getHistName("Pruned","SubJetMinJBP",flavor),    300,500);
  
  orderingPt300To500.push_back("Subjet CSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt300To500.push_back("Subjet JP (Pruned, Explicit JTA, SV Clustering)");
  orderingPt300To500.push_back("Subjet JBP (Pruned, Explicit JTA, SV Clustering)");
  //-------------------------------------------------
  graphsPt700ToInf["Subjet CSV (Pruned, Explicit JTA, SV Clustering)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinCSV"),    getHistName("Pruned","SubJetMinCSV",flavor),    700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet JP (Pruned, Explicit JTA, SV Clustering)"]     = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinJP"),     getHistName("Pruned","SubJetMinJP",flavor),     700,1100);
  graphsPt700ToInf["Subjet JBP (Pruned, Explicit JTA, SV Clustering)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinJBP"),    getHistName("Pruned","SubJetMinJBP",flavor),    700,1100);
  
  orderingPt700ToInf.push_back("Subjet CSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet JP (Pruned, Explicit JTA, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet JBP (Pruned, Explicit JTA, SV Clustering)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(),"",("btagperfcomp_Pt300to500_Subjets_ExplJTA_SV_Pruned_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(),0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(),"",("btagperfcomp_Pt700toInf_Subjets_ExplJTA_SV_Pruned_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(),0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //-------------------------------------------------
  // Subjets: Pruned IVFCSV comparison
  //-------------------------------------------------
  graphsPt300To500["Subjet IVFCSV (Pruned)"]                                           = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(),                                     (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),                                     getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA)"]                             = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),                         (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),                         getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, SV Clustering)"]                            = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),                        (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),                        getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]              = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(),            (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(),            getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)");
  //-------------------------------------------------
  graphsPt700ToInf["Subjet IVFCSV (Pruned)"]                                           = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(),                                     (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),                                     getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA)"]                             = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),                         (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),                         getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, SV Clustering)"]                            = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),                        (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),                        getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]              = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(),            (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(),            getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, SV momentum)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_Subjets_Pruned_IVFCSV_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_Subjets_Pruned_IVFCSV_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
  // Subjets: Kt IVFCSV comparison
  //-------------------------------------------------
  graphsPt300To500["Subjet IVFCSV (k_{T})"]                                           = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(),                                     (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),                                     getHistName("Kt","SubJetMinIVFCSV"), getHistName("Kt","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (k_{T}, Explicit JTA)"]                             = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),                         (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),                         getHistName("Kt","SubJetMinIVFCSV"), getHistName("Kt","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (k_{T}, SV Clustering)"]                            = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),                        (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),                        getHistName("Kt","SubJetMinIVFCSV"), getHistName("Kt","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]              = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(),            (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(),            getHistName("Kt","SubJetMinIVFCSV"), getHistName("Kt","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering, SV momentum)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinIVFCSV"), getHistName("Kt","SubJetMinIVFCSV",flavor), 300,500);

  orderingPt300To500.push_back("Subjet IVFCSV (k_{T})");
  orderingPt300To500.push_back("Subjet IVFCSV (k_{T}, Explicit JTA)");
  orderingPt300To500.push_back("Subjet IVFCSV (k_{T}, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering, SV momentum)");
  //-------------------------------------------------
  graphsPt700ToInf["Subjet IVFCSV (k_{T})"]                                           = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(),                                     (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),                                     getHistName("Kt","SubJetMinIVFCSV"), getHistName("Kt","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (k_{T}, Explicit JTA)"]                             = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),                         (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),                         getHistName("Kt","SubJetMinIVFCSV"), getHistName("Kt","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (k_{T}, SV Clustering)"]                            = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),                        (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),                        getHistName("Kt","SubJetMinIVFCSV"), getHistName("Kt","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)"]              = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(),            (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(),            getHistName("Kt","SubJetMinIVFCSV"), getHistName("Kt","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering, SV momentum)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_SVMomentum_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinIVFCSV"), getHistName("Kt","SubJetMinIVFCSV",flavor), 700,1100);

  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T})");
  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T}, Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T}, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T}, Explicit JTA, SV Clustering, SV momentum)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_Subjets_Kt_IVFCSV_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_Subjets_Kt_IVFCSV_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
  // Subjets: Pruned and Kt IVFCSV comparison
  //-------------------------------------------------
  graphsPt300To500["Subjet IVFCSV (Pruned)"]               = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(),             (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),             getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (k_{T})"]                = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(),             (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),             getHistName("Kt","SubJetMinIVFCSV"),     getHistName("Kt","SubJetMinIVFCSV",flavor),     300,500);
  graphsPt300To500["Subjet IVFCSV (k_{T}, Explicit JTA)"]  = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinIVFCSV"),     getHistName("Kt","SubJetMinIVFCSV",flavor),     300,500);

  orderingPt300To500.push_back("Subjet IVFCSV (Pruned)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt300To500.push_back("Subjet IVFCSV (k_{T})");
  orderingPt300To500.push_back("Subjet IVFCSV (k_{T}, Explicit JTA)");

  //-------------------------------------------------
  graphsPt700ToInf["Subjet IVFCSV (Pruned)"]               = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(),             (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),             getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (k_{T})"]                = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(),             (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),             getHistName("Kt","SubJetMinIVFCSV"),     getHistName("Kt","SubJetMinIVFCSV",flavor),     700,1100);
  graphsPt700ToInf["Subjet IVFCSV (k_{T}, Explicit JTA)"]  = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(), getHistName("Kt","SubJetMinIVFCSV"),     getHistName("Kt","SubJetMinIVFCSV",flavor),     700,1100);

  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T})");
  orderingPt700ToInf.push_back("Subjet IVFCSV (k_{T}, Explicit JTA)");

  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_Subjets_Pruned_Kt_IVFCSV_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_Subjets_Pruned_Kt_IVFCSV_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //-------------------------------------------------
  // Subjets: Pruned IVFCSV Explicit JTA comparison
  //-------------------------------------------------
  graphsPt300To500["Subjet IVFCSV (Pruned, SV Clustering)"]                           = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),             (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),             getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)"]             = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),             getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),             (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);

  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)");
  //-------------------------------------------------
  graphsPt700ToInf["Subjet IVFCSV (Pruned, SV Clustering)"]                           = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),             (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),             getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)"]             = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),             getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(),             (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);

  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA signal only)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, SV Clustering, Explicit JTA bkg only)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_Subjets_Pruned_IVFCSV_ExplJTA_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_Subjets_Pruned_IVFCSV_ExplJTA_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
  // Subjets: Pruned IVFCSV SV comparison
  //-------------------------------------------------
  graphsPt300To500["Subjet IVFCSV (Pruned)"]                            = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(),              (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),              getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, SV Clustering)"]             = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, SV Clustering signal only)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),              getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, SV Clustering bkg only)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root").c_str(),              (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);

  orderingPt300To500.push_back("Subjet IVFCSV (Pruned)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, SV Clustering signal only)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, SV Clustering bkg only)");
  //-------------------------------------------------
  graphsPt700ToInf["Subjet IVFCSV (Pruned)"]                            = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(),              (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),              getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, SV Clustering)"]             = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, SV Clustering signal only)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_PATTuple_v3.root").c_str(),              getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, SV Clustering bkg only)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PATTuple_v3.root").c_str(),              (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);

  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, SV Clustering signal only)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, SV Clustering bkg only)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_Subjets_Pruned_IVFCSV_SV_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_Subjets_Pruned_IVFCSV_SV_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();
  
  //-------------------------------------------------
  // Subjets: Pruned IVFCSV, Explicit JTA, SV comparison
  //-------------------------------------------------
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA)"]                            = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),              (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),              getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),              getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),              (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);

  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)");
  //-------------------------------------------------
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA)"]                            = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),              (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),              getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),              getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)"]    = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_PATTuple_v3.root").c_str(),              (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);

  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering signal only)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering bkg only)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_Subjets_Pruned_IVFCSV_ExplJTA_SV_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_Subjets_Pruned_IVFCSV_ExplJTA_SV_comparison_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
}

void makePlotsQCDFatvsSub(const string & dir = "ROOT_files", const string & algo = "AK", const double Ymin = 1e-3, const Int_t Logy = 1, const string & flavor = "", const string & ext = "eps")
{
  // for multiple plots on the same canvas

  // vectors storing the order of legend entries
  std::vector< std::string> orderingPt300To500;
  std::vector< std::string> orderingPt700ToInf;
  // maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt300To500;
  std::map< std::string,TGraph* > graphsPt700ToInf;
  
  //==========================================
  // Post-BTV-13-001 setup
  //==========================================
  // Fat jets and subjets
  //-------------------------------------------------
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"),       getHistName("Pruned","JetIVFCSV",flavor),       300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);

  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //-------------------------------------------------
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"),       getHistName("Pruned","JetIVFCSV",flavor),       700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);

  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_FatJets_Subjets_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_FatJets_Subjets_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
}

void makePlotsQCDAKvsCA(const string & dir = "ROOT_files", const double Ymin = 1e-3, const Int_t Logy = 1, const string & flavor = "", const string & ext = "eps")
{
  // for multiple plots on the same canvas

  // vectors storing the order of legend entries
  std::vector< std::string> orderingPt300To500;
  std::vector< std::string> orderingPt700ToInf;
  // maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt300To500;
  std::map< std::string,TGraph* > graphsPt700ToInf;
  
  //==========================================
  // Post-BTV-13-001 setup
  //==========================================
  // Fat jets and subjets, AK8 vs CA8
  //-------------------------------------------------
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA) - AK"]                       = getEfficiencyCurve((dir + "_AK8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"),       getHistName("Pruned","JetIVFCSV",flavor),       300,500);
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA) - CA"]                       = getEfficiencyCurve((dir + "_CA8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_CA8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"),       getHistName("Pruned","JetIVFCSV",flavor),       300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK"] = getEfficiencyCurve((dir + "_AK8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA"] = getEfficiencyCurve((dir + "_CA8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_CA8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);

  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA) - AK");
  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA) - CA");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA");
  //-------------------------------------------------
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA) - AK"]                       = getEfficiencyCurve((dir + "_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"),       getHistName("Pruned","JetIVFCSV",flavor),       700,1100);
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA) - CA"]                       = getEfficiencyCurve((dir + "_CA8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_CA8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"),       getHistName("Pruned","JetIVFCSV",flavor),       700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK"] = getEfficiencyCurve((dir + "_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA"] = getEfficiencyCurve((dir + "_CA8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_CA8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);

  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA) - AK");
  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA) - CA");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - AK");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering) - CA");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, "#splitline{R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt300to500_FatJets_Subjets_AK_vs_CA" + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, "#splitline{R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}",     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "", ("btagperfcomp_Pt700toInf_FatJets_Subjets_AK_vs_CA" + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
}

void makePlotsQCDStdJets(const string & dir = "ROOT_files", const string & algo = "AK", const double Ymin = 1e-3, const Int_t Logy = 1, const string & flavor = "", const string & ext = "eps")
{
  // for multiple plots on the same canvas

  // vectors storing the order of legend entries
  std::vector< std::string> orderingPt300To500;
  std::vector< std::string> orderingPt700ToInf;
  // maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt300To500;
  std::map< std::string,TGraph* > graphsPt700ToInf;
  
  //==========================================
  // Post-BTV-13-001 setup
  //==========================================
  // Fat jets and subjets, AK8 vs CA8
  //-------------------------------------------------
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"),       getHistName("Pruned","JetIVFCSV",flavor),       300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 300,500);
  graphsPt300To500["Matched AK4 jet IVFCSV"]                              = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","StdJetIVFCSV"),    getHistName("Pruned","StdJetIVFCSV",flavor),    300,500);
  graphsPt300To500["2 matched AK4 jets IVFCSV"]                           = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","StdJetMinIVFCSV"), getHistName("Pruned","StdJetMinIVFCSV",flavor), 300,500);

  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt300To500.push_back("Matched AK4 jet IVFCSV");
  orderingPt300To500.push_back("2 matched AK4 jets IVFCSV");
  //-------------------------------------------------
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA)"]                       = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","JetIVFCSV"),       getHistName("Pruned","JetIVFCSV",flavor),       700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"] = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","SubJetMinIVFCSV"), getHistName("Pruned","SubJetMinIVFCSV",flavor), 700,1100);
  graphsPt700ToInf["Matched AK4 jet IVFCSV"]                              = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","StdJetIVFCSV"),    getHistName("Pruned","StdJetIVFCSV",flavor),    700,1000);
  graphsPt700ToInf["2 matched AK4 jets IVFCSV"]                           = getEfficiencyCurve((dir + "_" + algo + "8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), (dir + "_" + algo + "8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root").c_str(), getHistName("Pruned","StdJetMinIVFCSV"), getHistName("Pruned","StdJetMinIVFCSV",flavor), 700,1000);

  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  orderingPt700ToInf.push_back("Matched AK4 jet IVFCSV");
  orderingPt700ToInf.push_back("2 matched AK4 jets IVFCSV");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, ("#splitline{" + algo + " R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(), "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "#DeltaR(AK4 jets,fat jet)<0.5", ("btagperfcomp_Pt300to500_FatJets_Subjets_StdJets_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, ("#splitline{" + algo + " R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}").c_str(),     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", ("Misid. probability (QCD" + (flavor != "" ? ", " + flavor : "")  + ")").c_str(), "#DeltaR(AK4 jets,fat jet)<0.5", ("btagperfcomp_Pt700toInf_FatJets_Subjets_StdJets_" + algo + (flavor != "" ? "_" + flavor : "") + "." + ext).c_str(), 0, 1, Ymin, 1, Logy);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
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
  
  graphsPt300To500["Fat Jet CSV (BTV-13-001)"] = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_BTV-13-001_PATTuple_v3.root", "ROOT_files_CA8/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3.root", getHistName("Pruned","JetCSV"),       300,500);
  graphsPt300To500["Subjet CSV (BTV-13-001)"]  = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_BTV-13-001_PATTuple_v3.root", "ROOT_files_CA8/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3.root", getHistName("Pruned","SubJetMinCSV"), 300,500);
  
  orderingPt300To500.push_back("Fat Jet CSV (BTV-13-001)");
  orderingPt300To500.push_back("Subjet CSV (BTV-13-001)");
  //-------------------------------------------------
  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"] = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_BTV-13-001_PATTuple_v3.root", "ROOT_files_CA8/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3.root", getHistName("Pruned","JetCSV"),       700,1100);
  graphsPt700ToInf["Subjet CSV (BTV-13-001)"]  = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_BTV-13-001_PATTuple_v3.root", "ROOT_files_CA8/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3.root", getHistName("Pruned","SubJetMinCSV"), 700,1100);
  
  orderingPt700ToInf.push_back("Fat Jet CSV (BTV-13-001)");
  orderingPt700ToInf.push_back("Subjet CSV (BTV-13-001)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misid. probability (QCD)", "", "btagperfcomp_Pt300to500_BTV-13-001.eps", 0, 1, 1E-3, 1, 1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}",     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misid. probability (QCD)", "", "btagperfcomp_Pt700toInf_BTV-13-001.eps", 0, 1, 1E-3, 1, 1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //==========================================
  // Post-BTV-13-001 setup
  //==========================================

  // Inclusive QCD
  makePlotsQCD("ROOT_files", "CA");          // CA8
  makePlotsQCDFatvsSub("ROOT_files", "CA");  // CA8
  makePlotsQCD("ROOT_files", "AK");          // AK8
  makePlotsQCDFatvsSub("ROOT_files", "AK");  // AK8
  makePlotsQCDAKvsCA("ROOT_files");          // AK8 vs CA8
  makePlotsQCDStdJets("ROOT_files", "AK");   // AK8

  // b jets
  makePlotsQCD("ROOT_files", "AK", 0, 0, "bJets");          // AK8
  makePlotsQCDFatvsSub("ROOT_files", "AK", 0, 0, "bJets");  // AK8
  makePlotsQCDAKvsCA("ROOT_files", 0, 0, "bJets");          // AK8 vs CA8
  makePlotsQCDStdJets("ROOT_files", "AK", 0, 0, "bJets");   // AK8

  // b jets from gluon splitting
  makePlotsQCD("ROOT_files", "AK", 0, 0, "bJetsGSP");          // AK8
  makePlotsQCDFatvsSub("ROOT_files", "AK", 0, 0, "bJetsGSP");  // AK8
  makePlotsQCDAKvsCA("ROOT_files", 0, 0, "bJetsGSP");          // AK8 vs CA8
  makePlotsQCDStdJets("ROOT_files", "AK", 0, 0, "bJetsGSP");   // AK8

  // c jets
  makePlotsQCD("ROOT_files", "AK", 1e-3, 1, "cJets");          // AK8
  makePlotsQCDFatvsSub("ROOT_files", "AK", 1e-3, 1, "cJets");  // AK8
  makePlotsQCDAKvsCA("ROOT_files", 1e-3, 1, "cJets");          // AK8 vs CA8
  makePlotsQCDStdJets("ROOT_files", "AK", 1e-3, 1, "cJets");   // AK8

  // gluon jets
  makePlotsQCD("ROOT_files", "AK", 1e-4, 1, "gluonJets");          // AK8
  makePlotsQCDFatvsSub("ROOT_files", "AK", 1e-4, 1, "gluonJets");  // AK8
  makePlotsQCDAKvsCA("ROOT_files", 1e-4, 1, "gluonJets");          // AK8 vs CA8
  makePlotsQCDStdJets("ROOT_files", "AK", 1e-4, 1, "gluonJets");   // AK8

  // uds jets
  makePlotsQCD("ROOT_files", "AK", 1e-4, 1, "udsJets");          // AK8
  makePlotsQCDFatvsSub("ROOT_files", "AK", 1e-4, 1, "udsJets");  // AK8
  makePlotsQCDStdJets("ROOT_files", "AK", 1e-4, 1, "udsJets");   // AK8
  makePlotsQCDAKvsCA("ROOT_files", 1e-4, 1, "udsJets");          // AK8 vs CA8

  // udsg jets
  makePlotsQCDFatvsSub("ROOT_files", "AK", 1e-4, 1, "udsgJets");  // AK8
  makePlotsQCDStdJets("ROOT_files", "AK", 1e-4, 1, "udsgJets");   // AK8

  //-------------------------------------------------
  // Fat jets and subjets, CA8
  //-------------------------------------------------
  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]                                        = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_BTV-13-001_PATTuple_v3.root",               "ROOT_files_CA8/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3.root",               getHistName("Pruned","JetCSV"),          300,500);
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA)"]                                   = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_CA8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","JetIVFCSV"),       300,500);
  graphsPt300To500["Subjet CSV (Pruned, BTV-13-001)"]                                 = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_BTV-13-001_PATTuple_v3.root",               "ROOT_files_CA8/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3.root",               getHistName("Pruned","SubJetMinCSV"),    300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_CA8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","SubJetMinIVFCSV"), 300,500);
  //graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, optimal 2D)"] = getEfficiencyCurve2D("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_CA8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName2D("Pruned","IVFCSV","300to500"));

  orderingPt300To500.push_back("Fat Jet CSV (BTV-13-001)");
  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt300To500.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, optimal 2D)");
  //-------------------------------------------------
  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"]                                        = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_BTV-13-001_PATTuple_v3.root",               "ROOT_files_CA8/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3.root",               getHistName("Pruned","JetCSV"),          700,1100);
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA)"]                                   = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_CA8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","JetIVFCSV"),       700,1100);
  graphsPt700ToInf["Subjet CSV (Pruned, BTV-13-001)"]                                 = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_BTV-13-001_PATTuple_v3.root",               "ROOT_files_CA8/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3.root",               getHistName("Pruned","SubJetMinCSV"),    700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_CA8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","SubJetMinIVFCSV"), 700,1100);
  //graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, optimal 2D)"] = getEfficiencyCurve2D("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_CA8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName2D("Pruned","IVFCSV","700toInf"));

  orderingPt700ToInf.push_back("Fat Jet CSV (BTV-13-001)");
  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet CSV (Pruned, BTV-13-001)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, optimal 2D)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misid. probability (QCD)", "", "btagperfcomp_Pt300to500_FatJets_Subjets_CA_BTV-13-001.eps", 0, 1, 1E-3, 1, 1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}",     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misid. probability (QCD)", "", "btagperfcomp_Pt700toInf_FatJets_Subjets_CA_BTV-13-001.eps", 0, 1, 1E-3, 1, 1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
  // Fat jets and subjets, AK8
  //-------------------------------------------------
  graphsPt300To500["Fat Jet CSV (CA8, BTV-13-001)"]                                   = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_BTV-13-001_PATTuple_v3.root",               "ROOT_files_CA8/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3.root",               getHistName("Pruned","JetCSV"),          300,500);
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA)"]                                   = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","JetIVFCSV"),       300,500);
  graphsPt300To500["Subjet CSV (CA8 Pruned, BTV-13-001)"]                             = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_BTV-13-001_PATTuple_v3.root",               "ROOT_files_CA8/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3.root",               getHistName("Pruned","SubJetMinCSV"),    300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","SubJetMinIVFCSV"), 300,500);
  //graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, optimal 2D)"] = getEfficiencyCurve2D("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName2D("Pruned","IVFCSV","300to500"));

  orderingPt300To500.push_back("Fat Jet CSV (CA8, BTV-13-001)");
  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt300To500.push_back("Subjet CSV (CA8 Pruned, BTV-13-001)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, optimal 2D)");
  //-------------------------------------------------
  graphsPt700ToInf["Fat Jet CSV (CA8, BTV-13-001)"]                                   = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_BTV-13-001_PATTuple_v3.root",               "ROOT_files_CA8/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3.root",               getHistName("Pruned","JetCSV"),          700,1100);
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA)"]                                   = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","JetIVFCSV"),       700,1100);
  graphsPt700ToInf["Subjet CSV (CA8 Pruned, BTV-13-001)"]                             = getEfficiencyCurve("ROOT_files_CA8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_BTV-13-001_PATTuple_v3.root",               "ROOT_files_CA8/QCDPythia6_HiggsTagging_BTV-13-001_PATTuple_v3.root",               getHistName("Pruned","SubJetMinCSV"),    700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","SubJetMinIVFCSV"), 700,1100);
  //graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, optimal 2D)"] = getEfficiencyCurve2D("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/QCDPythia6_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName2D("Pruned","IVFCSV","700toInf"));

  orderingPt700ToInf.push_back("Fat Jet CSV (CA8, BTV-13-001)");
  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet CSV (CA8 Pruned, BTV-13-001)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering, optimal 2D)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, "#splitline{AK R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misid. probability (QCD)", "", "btagperfcomp_Pt300to500_FatJets_Subjets_AK_BTV-13-001.eps", 0, 1, 1E-3, 1, 1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, "#splitline{AK R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}",     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misid. probability (QCD)", "", "btagperfcomp_Pt700toInf_FatJets_Subjets_AK_BTV-13-001.eps", 0, 1, 1E-3, 1, 1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
  // Fat jets and subjets, AK8, boosted hadronic W
  //-------------------------------------------------
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA)"]                                   = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","JetIVFCSV"),       300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","SubJetMinIVFCSV"), 300,500);

  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //-------------------------------------------------
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA)"]                                   = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","JetIVFCSV"),       700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","SubJetMinIVFCSV"), 700,1100);

  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, "#splitline{AK R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misid. probability (Hadronic W)", "", "btagperfcomp_Pt300to500_FatJets_Subjets_AK_Hadronic_W.eps", 0, 1, 1E-3, 1, 1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, "#splitline{AK R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}",     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misid. probability (Hadronic W)", "", "btagperfcomp_Pt700toInf_FatJets_Subjets_AK_Hadronic_W.eps", 0, 1, 1E-3, 1, 1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
  // Fat jets and subjets, AK8, boosted hadronic Z
  //-------------------------------------------------
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA)"]                                   = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","JetIVFCSV"),       300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","SubJetMinIVFCSV"), 300,500);

  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //-------------------------------------------------
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA)"]                                   = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","JetIVFCSV"),       700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","SubJetMinIVFCSV"), 700,1100);

  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, "#splitline{AK R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misid. probability (Hadronic Z)", "", "btagperfcomp_Pt300to500_FatJets_Subjets_AK_Hadronic_Z.eps", 0, 1, 1E-2, 1, 1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, "#splitline{AK R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}",     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misid. probability (Hadronic Z)", "", "btagperfcomp_Pt700toInf_FatJets_Subjets_AK_Hadronic_Z.eps", 0, 1, 1E-2, 1, 1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
  // Fat jets and subjets, AK8, boosted hadronic top
  //-------------------------------------------------
  graphsPt300To500["Fat Jet IVFCSV (Explicit JTA)"]                                   = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","JetIVFCSV"),       300,500);
  graphsPt300To500["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/TprimeToTHinc_M-1000_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","SubJetMinIVFCSV"), 300,500);

  orderingPt300To500.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt300To500.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //-------------------------------------------------
  graphsPt700ToInf["Fat Jet IVFCSV (Explicit JTA)"]                                   = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/TprimeToTHinc_M-1700_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","JetIVFCSV"),       700,1100);
  graphsPt700ToInf["Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)"]             = getEfficiencyCurve("ROOT_files_AK8/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", "ROOT_files_AK8/TprimeToTHinc_M-1700_HiggsTagging_ExplicitJTA_SVClustering_PATTuple_v3.root", getHistName("Pruned","SubJetMinIVFCSV"), 700,1100);

  orderingPt700ToInf.push_back("Fat Jet IVFCSV (Explicit JTA)");
  orderingPt700ToInf.push_back("Subjet IVFCSV (Pruned, Explicit JTA, SV Clustering)");
  //-------------------------------------------------
  plotEfficiencyCurves(graphsPt300To500,orderingPt300To500, "#splitline{AK R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misid. probability (Hadronic top)", "", "btagperfcomp_Pt300to500_FatJets_Subjets_AK_Hadronic_top.eps", 0, 1, 1E-2, 1, 1);
  plotEfficiencyCurves(graphsPt700ToInf,orderingPt700ToInf, "#splitline{AK R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}",     "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misid. probability (Hadronic top)", "", "btagperfcomp_Pt700toInf_FatJets_Subjets_AK_Hadronic_top.eps", 0, 1, 1E-2, 1, 1);

  graphsPt300To500.clear();
  graphsPt700ToInf.clear();

  orderingPt300To500.clear();
  orderingPt700ToInf.clear();

  //-------------------------------------------------
}
