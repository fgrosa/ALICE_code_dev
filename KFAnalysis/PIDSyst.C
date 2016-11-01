#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <TInterpreter.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <THnSparse.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TLine.h>
#include <TList.h>
#include <TF1.h>
#include <TLatex.h>
#include <TArrayD.h>

#endif

const Int_t colors[] = {kRed,kBlue,kBlack,kOrange+7,kGreen+3,kMagenta,kCyan,kYellow+3,kOrange+3};

void PIDSyst(Double_t BR=0.0913) {
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.19);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.8,"y");
  gStyle->SetLegendBorderSize(0);

  TFile infile("HFPtSpectrum_combinedFD_cutset1.root","UPDATE");
  TH1F* hCross = (TH1F*)infile.Get("histoSigmaCorr");
  hCross->SetDirectory(0);
  hCross->SetMarkerStyle(20);
  hCross->SetMarkerSize(1.5);
  hCross->SetMarkerColor(colors[0]);
  hCross->SetLineColor(colors[0]);
  hCross->SetLineWidth(2);
  hCross->Scale(1./(1000000*BR));
  hCross->GetXaxis()->SetLabelFont(42);
  hCross->GetXaxis()->SetTitleFont(42);
  hCross->GetYaxis()->SetLabelFont(42);
  hCross->GetYaxis()->SetTitleFont(42);
  hCross->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}} (#mub c/GeV)");
  hCross->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hCross->GetYaxis()->SetTitleOffset(1.6);
  hCross->GetYaxis()->SetTitleSize(0.05);
  hCross->GetXaxis()->SetTitleSize(0.05);
  hCross->GetYaxis()->SetLabelSize(0.05);
  hCross->GetXaxis()->SetLabelSize(0.05);
  infile.Close();
  
  TFile infileNoPID("HFPtSpectrum_combinedFD_cutset1_noPID.root","UPDATE");
  TH1F* hCrossNoPID = (TH1F*)infileNoPID.Get("histoSigmaCorr");
  hCrossNoPID->SetDirectory(0);
  hCrossNoPID->SetMarkerStyle(20);
  hCrossNoPID->SetMarkerSize(1.5);
  hCrossNoPID->SetMarkerColor(colors[1]);
  hCrossNoPID->SetLineColor(colors[1]);
  hCrossNoPID->SetLineWidth(2);
  hCrossNoPID->Scale(1./(1000000*BR));
  hCrossNoPID->GetXaxis()->SetLabelFont(42);
  hCrossNoPID->GetXaxis()->SetTitleFont(42);
  hCrossNoPID->GetYaxis()->SetLabelFont(42);
  hCrossNoPID->GetYaxis()->SetTitleFont(42);
  infileNoPID.Close();

  const Int_t nPtBins=hCross->GetNbinsX();
  TArrayD* ptarray=(TArrayD*)hCross->GetXaxis()->GetXbins();
  Double_t *PtLims=(Double_t*)ptarray->GetArray();
  
  TH1F *hRatio = (TH1F*)hCross->Clone();
  hRatio->Divide(hCrossNoPID,hCross,1.,1.,"B");
  hRatio->SetDirectory(0);
  hRatio->SetTitle("");
  hRatio->GetXaxis()->SetTitleFont(42);
  hRatio->GetYaxis()->SetTitleFont(42);
  hRatio->GetXaxis()->SetLabelFont(42);
  hRatio->GetYaxis()->SetLabelFont(42);
  hRatio->GetXaxis()->SetTitleSize(0.05);
  hRatio->GetYaxis()->SetTitleSize(0.05);
  hRatio->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}}(no PID)/#frac{d#sigma}{d#it{p}_{T}}(PID)");
  hRatio->GetXaxis()->SetTitle(hCross->GetXaxis()->GetTitle());
 
  TLegend* l = new TLegend(0.4,0.7,0.85,0.89);
  l->SetTextSize(0.035);
  l->AddEntry(hCross,"PID","lpe");
  l->AddEntry(hCrossNoPID,"no PID","lpe");

  TCanvas *c = new TCanvas("c","",800,800);
  c->SetLogy();
  hCross->Draw();
  hCrossNoPID->Draw("same");
  l->Draw("same");

  TCanvas *cRatio = new TCanvas("cRatio","",800,800);
  hRatio->GetYaxis()->SetRangeUser(0.7,1.4);
  hRatio->Draw();
  
  c->SaveAs("CrossSection_KF_noPID.eps");
  cRatio->SaveAs("Ratio_KF_noPID.eps");

  delete c;
  delete cRatio;
}
