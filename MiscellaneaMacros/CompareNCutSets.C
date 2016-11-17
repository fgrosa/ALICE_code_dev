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
#include <TF1.h>
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

#endif

//_______________________________________________________________________________________________________
//function prototypes
void RawYieldsComparison(TString filename="RawYields");
void MassSpectraComparison(TString filename="RawYields");
//void EfficiencyComparison(TString filename="Efficiency", TString directoryname=".");
//void FPromptComparison(TString filename="fprompt_unbinned_sigmafree", TString directoryname=".");
//void CrossSectionComparison(TString filename1="DplusCrossSection", TString directoryname=".");
//_______________________________________________________________________________________________________
//global variables

const Int_t nCutSets=4;
const TString directoryname[nCutSets] = {"Cent3050","Cent3050","Cent3050","Cent3050"};
const TString cutsetsname[nCutSets] = {"_3050","_3050_d0cut","_3050_topocut","_3050_topod0cut"};
const TString legendnames[nCutSets] = {"w/o d_{0} and d_{0}-d_{0}exp cut","|d_{0}| cut","|norm max d_{0}-d_{0}exp| cut", "d_{0} and d_{0}-d_{0}exp cut"};
const Int_t colors[nCutSets] = {kRed,kBlue,kGreen+2,kBlack};
const Int_t styles[nCutSets] = {20,21,22,33};

const TString outpostpend = "_diffcuts";
const TString outprepend = "Cent3050/";

//_______________________________________________________________________________________________________
void RawYieldsComparison(TString filename) {
  
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.3);
  
  TH1F** hRawYields = new TH1F*[nCutSets];
  TH1F** hRawYieldsSigma = new TH1F*[nCutSets];
  TH1F** hRawYieldsMean = new TH1F*[nCutSets];
  TH1F** hRawYieldsSignificance = new TH1F*[nCutSets];
  TH1F** hRawYieldsSoverB = new TH1F*[nCutSets];
  TH1F** hRawYieldsSignal = new TH1F*[nCutSets];
  TH1F** hRawYieldsBkg = new TH1F*[nCutSets];
  TH1F** hRawYieldsChiSquare = new TH1F*[nCutSets];
  TH1F** hRawYieldsRelError = new TH1F*[nCutSets];
  
  TH1F** hRawYieldsRatio = new TH1F*[nCutSets-1];
  TH1F** hRawYieldsSigmaRatio = new TH1F*[nCutSets-1];
  TH1F** hRawYieldsMeanRatio = new TH1F*[nCutSets-1];
  TH1F** hRawYieldsSignificanceRatio = new TH1F*[nCutSets-1];
  TH1F** hRawYieldsSoverBRatio = new TH1F*[nCutSets-1];
  TH1F** hRawYieldsSignalRatio = new TH1F*[nCutSets-1];
  TH1F** hRawYieldsBkgRatio = new TH1F*[nCutSets-1];
  TH1F** hRawYieldsRelErrorRatio = new TH1F*[nCutSets-1];

  Double_t massD = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  TLine* line = new TLine(2,massD,100,massD);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->SetLineColor(kBlack);
  TLine* lineatone = new TLine(2,1,100,1);
  lineatone->SetLineWidth(2);
  lineatone->SetLineStyle(7);
  lineatone->SetLineColor(colors[0]);
  
  for(Int_t iCutSet=0; iCutSet<nCutSets; iCutSet++) {
    TFile* infile = TFile::Open(Form("%s/%s%s.root",directoryname[iCutSet].Data(),filename.Data(),cutsetsname[iCutSet].Data()),"READ");
    hRawYields[iCutSet] = (TH1F*)infile->Get("hRawYields");
    hRawYieldsSigma[iCutSet] = (TH1F*)infile->Get("hRawYieldsSigma");
    hRawYieldsMean[iCutSet] = (TH1F*)infile->Get("hRawYieldsMean");
    hRawYieldsSignificance[iCutSet] = (TH1F*)infile->Get("hRawYieldsSignificance");
    hRawYieldsSoverB[iCutSet] = (TH1F*)infile->Get("hRawYieldsSoverB");
    hRawYieldsSignal[iCutSet] = (TH1F*)infile->Get("hRawYieldsSignal");
    hRawYieldsBkg[iCutSet] = (TH1F*)infile->Get("hRawYieldsBkg");
    hRawYieldsChiSquare[iCutSet] = (TH1F*)infile->Get("hRawYieldsChiSquare");
    hRawYields[iCutSet]->SetDirectory(0);
    hRawYieldsSigma[iCutSet]->SetDirectory(0);
    hRawYieldsMean[iCutSet]->SetDirectory(0);
    hRawYieldsSignificance[iCutSet]->SetDirectory(0);
    hRawYieldsSoverB[iCutSet]->SetDirectory(0);
    hRawYieldsSignal[iCutSet]->SetDirectory(0);
    hRawYieldsBkg[iCutSet]->SetDirectory(0);
    hRawYieldsChiSquare[iCutSet]->SetDirectory(0);
    infile->Close();
    
    hRawYields[iCutSet]->GetYaxis()->SetTitleSize(0.07);
    hRawYields[iCutSet]->GetYaxis()->SetLabelSize(0.06);
    hRawYields[iCutSet]->GetYaxis()->SetTitleOffset(1.);
    hRawYields[iCutSet]->GetYaxis()->SetRangeUser(0.1,hRawYields[0]->GetMaximum()*1.5);
    hRawYields[iCutSet]->SetLineColor(colors[iCutSet]);
    hRawYields[iCutSet]->SetMarkerColor(colors[iCutSet]);
    hRawYields[iCutSet]->SetMarkerStyle(styles[iCutSet]);
    hRawYieldsSigma[iCutSet]->GetYaxis()->SetTitleSize(0.07);
    hRawYieldsSigma[iCutSet]->GetYaxis()->SetLabelSize(0.06);
    hRawYieldsSigma[iCutSet]->GetYaxis()->SetTitleOffset(1.);
    hRawYieldsSigma[iCutSet]->GetYaxis()->SetRangeUser(0.0021,hRawYieldsSigma[0]->GetMaximum()*2.);
    hRawYieldsSigma[iCutSet]->SetLineColor(colors[iCutSet]);
    hRawYieldsSigma[iCutSet]->SetMarkerColor(colors[iCutSet]);
    hRawYieldsSigma[iCutSet]->SetMarkerStyle(styles[iCutSet]);
    hRawYieldsMean[iCutSet]->GetYaxis()->SetTitleSize(0.07);
    hRawYieldsMean[iCutSet]->GetYaxis()->SetLabelSize(0.06);
    hRawYieldsMean[iCutSet]->GetYaxis()->SetTitleOffset(1.);
    hRawYieldsMean[iCutSet]->GetYaxis()->SetRangeUser(hRawYieldsMean[0]->GetMinimum()*0.9935,hRawYieldsMean[0]->GetMaximum()*1.014);
    hRawYieldsMean[iCutSet]->SetLineColor(colors[iCutSet]);
    hRawYieldsMean[iCutSet]->SetMarkerColor(colors[iCutSet]);
    hRawYieldsMean[iCutSet]->SetMarkerStyle(styles[iCutSet]);
    hRawYieldsSignificance[iCutSet]->GetYaxis()->SetTitleSize(0.07);
    hRawYieldsSignificance[iCutSet]->GetYaxis()->SetLabelSize(0.06);
    hRawYieldsSignificance[iCutSet]->GetYaxis()->SetTitleOffset(1.);
    hRawYieldsSignificance[iCutSet]->GetYaxis()->SetRangeUser(0.5,hRawYieldsSignificance[0]->GetMaximum()*1.8);
    hRawYieldsSignificance[iCutSet]->SetLineColor(colors[iCutSet]);
    hRawYieldsSignificance[iCutSet]->SetMarkerColor(colors[iCutSet]);
    hRawYieldsSignificance[iCutSet]->SetMarkerStyle(styles[iCutSet]);
    hRawYieldsSoverB[iCutSet]->GetYaxis()->SetTitleSize(0.07);
    hRawYieldsSoverB[iCutSet]->GetYaxis()->SetLabelSize(0.06);
    hRawYieldsSoverB[iCutSet]->GetYaxis()->SetTitleOffset(1.);
    hRawYieldsSoverB[iCutSet]->GetYaxis()->SetRangeUser(0.01,hRawYieldsSoverB[0]->GetMaximum()*2);
    hRawYieldsSoverB[iCutSet]->SetLineColor(colors[iCutSet]);
    hRawYieldsSoverB[iCutSet]->SetMarkerColor(colors[iCutSet]);
    hRawYieldsSoverB[iCutSet]->SetMarkerStyle(styles[iCutSet]);
    hRawYieldsSignal[iCutSet]->GetYaxis()->SetTitleSize(0.07);
    hRawYieldsSignal[iCutSet]->GetYaxis()->SetLabelSize(0.06);
    hRawYieldsSignal[iCutSet]->GetYaxis()->SetTitleOffset(1.);
    hRawYieldsSignal[iCutSet]->GetYaxis()->SetRangeUser(0.01,hRawYieldsSignal[0]->GetMaximum()*2);
    hRawYieldsSignal[iCutSet]->SetLineColor(colors[iCutSet]);
    hRawYieldsSignal[iCutSet]->SetMarkerColor(colors[iCutSet]);
    hRawYieldsSignal[iCutSet]->SetMarkerStyle(styles[iCutSet]);
    hRawYieldsBkg[iCutSet]->GetYaxis()->SetTitleSize(0.07);
    hRawYieldsBkg[iCutSet]->GetYaxis()->SetLabelSize(0.06);
    hRawYieldsBkg[iCutSet]->GetYaxis()->SetTitleOffset(1.);
    hRawYieldsBkg[iCutSet]->GetYaxis()->SetRangeUser(0.01,hRawYieldsBkg[0]->GetMaximum()*2);
    hRawYieldsBkg[iCutSet]->SetLineColor(colors[iCutSet]);
    hRawYieldsBkg[iCutSet]->SetMarkerColor(colors[iCutSet]);
    hRawYieldsBkg[iCutSet]->SetMarkerStyle(styles[iCutSet]);
    hRawYieldsChiSquare[iCutSet]->GetYaxis()->SetTitleSize(0.07);
    hRawYieldsChiSquare[iCutSet]->GetYaxis()->SetLabelSize(0.06);
    hRawYieldsChiSquare[iCutSet]->GetYaxis()->SetTitleOffset(1.);
    hRawYieldsChiSquare[iCutSet]->GetYaxis()->SetRangeUser(0.,3.5);
    hRawYieldsChiSquare[iCutSet]->SetLineColor(colors[iCutSet]);
    hRawYieldsChiSquare[iCutSet]->SetMarkerColor(colors[iCutSet]);
    hRawYieldsChiSquare[iCutSet]->SetMarkerStyle(styles[iCutSet]);
    hRawYieldsChiSquare[iCutSet]->GetYaxis()->SetTitleOffset(1.4);

    hRawYieldsRelError[iCutSet] = (TH1F*)hRawYields[iCutSet]->Clone();
    hRawYieldsRelError[iCutSet]->SetName(Form("hRawYieldsRelError_%d",iCutSet));
    for(Int_t iPt=0; iPt<hRawYieldsRelError[iCutSet]->GetNbinsX(); iPt++) {
      hRawYieldsRelError[iCutSet]->SetBinContent(iPt+1,hRawYields[iCutSet]->GetBinError(iPt+1)/hRawYields[iCutSet]->GetBinContent(iPt+1));
      hRawYieldsRelError[iCutSet]->SetBinError(iPt+1,1.e-15);
    }
    hRawYieldsRelError[iCutSet]->GetYaxis()->SetTitle("#sigma_{Y}(stat)/Y");
    hRawYieldsRelError[iCutSet]->GetYaxis()->SetRangeUser(-0.03,1.);

    if(iCutSet>0) {
      hRawYieldsRatio[iCutSet-1] = (TH1F*)hRawYields[iCutSet]->Clone();
      hRawYieldsRatio[iCutSet-1]->SetDirectory(0);
      hRawYieldsRatio[iCutSet-1]->Divide(hRawYields[iCutSet],hRawYields[0],1.,1.,"");
      for(Int_t iPt=0; iPt<hRawYieldsRatio[iCutSet-1]->GetNbinsX(); iPt++) {
        hRawYieldsRatio[iCutSet-1]->SetBinError(iPt+1,1.e-10);
      }
      hRawYieldsRatio[iCutSet-1]->GetYaxis()->SetTitle("Ratio");
      hRawYieldsRatio[iCutSet-1]->GetYaxis()->SetRangeUser(hRawYieldsRatio[0]->GetMinimum()*0.5,hRawYieldsRatio[0]->GetMaximum()*1.5);
      hRawYieldsRatio[iCutSet-1]->GetYaxis()->SetTitleSize(0.11);
      hRawYieldsRatio[iCutSet-1]->GetXaxis()->SetTitleSize(0.11);
      hRawYieldsRatio[iCutSet-1]->GetYaxis()->SetLabelSize(0.10);
      hRawYieldsRatio[iCutSet-1]->GetXaxis()->SetLabelSize(0.10);
      hRawYieldsRatio[iCutSet-1]->GetXaxis()->SetTitleOffset(1.2);
      hRawYieldsRatio[iCutSet-1]->GetYaxis()->SetTitleOffset(0.6);
      hRawYieldsRatio[iCutSet-1]->GetYaxis()->SetDecimals(2);
      
      hRawYieldsSigmaRatio[iCutSet-1] = (TH1F*)hRawYieldsSigma[iCutSet]->Clone();
      hRawYieldsSigmaRatio[iCutSet-1]->SetDirectory(0);
      hRawYieldsSigmaRatio[iCutSet-1]->Divide(hRawYieldsSigma[iCutSet],hRawYieldsSigma[0],1.,1.,"");
      for(Int_t iPt=0; iPt<hRawYieldsSigmaRatio[iCutSet-1]->GetNbinsX(); iPt++) {
        hRawYieldsSigmaRatio[iCutSet-1]->SetBinError(iPt+1,1.e-10);
      }
      hRawYieldsSigmaRatio[iCutSet-1]->GetYaxis()->SetTitle("Ratio");
      hRawYieldsSigmaRatio[iCutSet-1]->GetYaxis()->SetRangeUser(hRawYieldsSigmaRatio[0]->GetMinimum()*0.1,hRawYieldsSigmaRatio[0]->GetMaximum()*1.55);
      hRawYieldsSigmaRatio[iCutSet-1]->GetYaxis()->SetTitleSize(0.11);
      hRawYieldsSigmaRatio[iCutSet-1]->GetXaxis()->SetTitleSize(0.11);
      hRawYieldsSigmaRatio[iCutSet-1]->GetYaxis()->SetLabelSize(0.10);
      hRawYieldsSigmaRatio[iCutSet-1]->GetXaxis()->SetLabelSize(0.10);
      hRawYieldsSigmaRatio[iCutSet-1]->GetXaxis()->SetTitleOffset(1.2);
      hRawYieldsSigmaRatio[iCutSet-1]->GetYaxis()->SetTitleOffset(0.6);
      hRawYieldsSigmaRatio[iCutSet-1]->GetYaxis()->SetDecimals(2);
      
      hRawYieldsMeanRatio[iCutSet-1] = (TH1F*)hRawYieldsMean[iCutSet]->Clone();
      hRawYieldsMeanRatio[iCutSet-1]->SetDirectory(0);
      hRawYieldsMeanRatio[iCutSet-1]->Divide(hRawYieldsMean[iCutSet],hRawYieldsMean[0],1.,1.,"");
      for(Int_t iPt=0; iPt<hRawYieldsMeanRatio[iCutSet-1]->GetNbinsX(); iPt++) {
        hRawYieldsMeanRatio[iCutSet-1]->SetBinError(iPt+1,1.e-10);
      }
      hRawYieldsMeanRatio[iCutSet-1]->GetYaxis()->SetTitle("Ratio");
      hRawYieldsMeanRatio[iCutSet-1]->GetYaxis()->SetRangeUser(hRawYieldsMeanRatio[0]->GetMinimum()*0.992,hRawYieldsMeanRatio[0]->GetMaximum()*1.017);
      hRawYieldsMeanRatio[iCutSet-1]->GetYaxis()->SetTitleSize(0.11);
      hRawYieldsMeanRatio[iCutSet-1]->GetXaxis()->SetTitleSize(0.11);
      hRawYieldsMeanRatio[iCutSet-1]->GetYaxis()->SetLabelSize(0.10);
      hRawYieldsMeanRatio[iCutSet-1]->GetXaxis()->SetLabelSize(0.10);
      hRawYieldsMeanRatio[iCutSet-1]->GetXaxis()->SetTitleOffset(1.2);
      hRawYieldsMeanRatio[iCutSet-1]->GetYaxis()->SetTitleOffset(0.6);
      hRawYieldsMeanRatio[iCutSet-1]->GetYaxis()->SetDecimals(2);
      
      hRawYieldsSignificanceRatio[iCutSet-1] = (TH1F*)hRawYieldsSignificance[iCutSet]->Clone();
      hRawYieldsSignificanceRatio[iCutSet-1]->SetDirectory(0);
      hRawYieldsSignificanceRatio[iCutSet-1]->Divide(hRawYieldsSignificance[iCutSet],hRawYieldsSignificance[0],1.,1.,"");
      for(Int_t iPt=0; iPt<hRawYieldsSignificanceRatio[iCutSet-1]->GetNbinsX(); iPt++) {
        hRawYieldsSignificanceRatio[iCutSet-1]->SetBinError(iPt+1,1.e-10);
      }
      hRawYieldsSignificanceRatio[iCutSet-1]->GetYaxis()->SetTitle("Ratio");
      hRawYieldsSignificanceRatio[iCutSet-1]->GetYaxis()->SetRangeUser(hRawYieldsSignificanceRatio[0]->GetMinimum()*0.55,hRawYieldsSignificanceRatio[0]->GetMaximum()*1.45);
      hRawYieldsSignificanceRatio[iCutSet-1]->GetYaxis()->SetTitleSize(0.11);
      hRawYieldsSignificanceRatio[iCutSet-1]->GetXaxis()->SetTitleSize(0.11);
      hRawYieldsSignificanceRatio[iCutSet-1]->GetYaxis()->SetLabelSize(0.10);
      hRawYieldsSignificanceRatio[iCutSet-1]->GetXaxis()->SetLabelSize(0.10);
      hRawYieldsSignificanceRatio[iCutSet-1]->GetXaxis()->SetTitleOffset(1.2);
      hRawYieldsSignificanceRatio[iCutSet-1]->GetYaxis()->SetTitleOffset(0.6);
      hRawYieldsSignificanceRatio[iCutSet-1]->GetYaxis()->SetDecimals(2);
      
      hRawYieldsSoverBRatio[iCutSet-1] = (TH1F*)hRawYieldsSoverB[iCutSet]->Clone();
      hRawYieldsSoverBRatio[iCutSet-1]->SetDirectory(0);
      hRawYieldsSoverBRatio[iCutSet-1]->Divide(hRawYieldsSoverB[iCutSet],hRawYieldsSoverB[0],1.,1.,"");
      for(Int_t iPt=0; iPt<hRawYieldsRatio[iCutSet-1]->GetNbinsX(); iPt++) {
        hRawYieldsSoverBRatio[iCutSet-1]->SetBinError(iPt+1,1.e-10);
      }
      hRawYieldsSoverBRatio[iCutSet-1]->GetYaxis()->SetTitle("Ratio");
      hRawYieldsSoverBRatio[iCutSet-1]->GetYaxis()->SetRangeUser(hRawYieldsSoverBRatio[0]->GetMinimum()*0.5,hRawYieldsSoverBRatio[0]->GetMaximum()*1.5);
      hRawYieldsSoverBRatio[iCutSet-1]->GetYaxis()->SetTitleSize(0.11);
      hRawYieldsSoverBRatio[iCutSet-1]->GetXaxis()->SetTitleSize(0.11);
      hRawYieldsSoverBRatio[iCutSet-1]->GetYaxis()->SetLabelSize(0.10);
      hRawYieldsSoverBRatio[iCutSet-1]->GetXaxis()->SetLabelSize(0.10);
      hRawYieldsSoverBRatio[iCutSet-1]->GetXaxis()->SetTitleOffset(1.2);
      hRawYieldsSoverBRatio[iCutSet-1]->GetYaxis()->SetTitleOffset(0.6);
      hRawYieldsSoverBRatio[iCutSet-1]->GetYaxis()->SetDecimals(2);

      hRawYieldsSignalRatio[iCutSet-1] = (TH1F*)hRawYieldsSignal[iCutSet]->Clone();
      hRawYieldsSignalRatio[iCutSet-1]->SetDirectory(0);
      hRawYieldsSignalRatio[iCutSet-1]->Divide(hRawYieldsSignal[iCutSet],hRawYieldsSignal[0],1.,1.,"");
      for(Int_t iPt=0; iPt<hRawYieldsRatio[iCutSet-1]->GetNbinsX(); iPt++) {
        hRawYieldsSignalRatio[iCutSet-1]->SetBinError(iPt+1,1.e-10);
      }
      hRawYieldsSignalRatio[iCutSet-1]->GetYaxis()->SetTitle("Ratio");
      hRawYieldsSignalRatio[iCutSet-1]->GetYaxis()->SetRangeUser(hRawYieldsSignalRatio[0]->GetMinimum()*0.5,hRawYieldsSignalRatio[0]->GetMaximum()*1.5);
      hRawYieldsSignalRatio[iCutSet-1]->GetYaxis()->SetTitleSize(0.11);
      hRawYieldsSignalRatio[iCutSet-1]->GetXaxis()->SetTitleSize(0.11);
      hRawYieldsSignalRatio[iCutSet-1]->GetYaxis()->SetLabelSize(0.10);
      hRawYieldsSignalRatio[iCutSet-1]->GetXaxis()->SetLabelSize(0.10);
      hRawYieldsSignalRatio[iCutSet-1]->GetXaxis()->SetTitleOffset(1.2);
      hRawYieldsSignalRatio[iCutSet-1]->GetYaxis()->SetTitleOffset(0.6);
      hRawYieldsSignalRatio[iCutSet-1]->GetYaxis()->SetDecimals(2);
      
      hRawYieldsBkgRatio[iCutSet-1] = (TH1F*)hRawYieldsBkg[iCutSet]->Clone();
      hRawYieldsBkgRatio[iCutSet-1]->SetDirectory(0);
      hRawYieldsBkgRatio[iCutSet-1]->Divide(hRawYieldsBkg[iCutSet],hRawYieldsBkg[0],1.,1.,"");
      for(Int_t iPt=0; iPt<hRawYieldsRatio[iCutSet-1]->GetNbinsX(); iPt++) {
        hRawYieldsBkgRatio[iCutSet-1]->SetBinError(iPt+1,1.e-10);
      }
      hRawYieldsBkgRatio[iCutSet-1]->GetYaxis()->SetTitle("Ratio");
      hRawYieldsBkgRatio[iCutSet-1]->GetYaxis()->SetRangeUser(hRawYieldsBkgRatio[0]->GetMinimum()*0.5,hRawYieldsBkgRatio[0]->GetMaximum()*1.5);
      hRawYieldsBkgRatio[iCutSet-1]->GetYaxis()->SetTitleSize(0.11);
      hRawYieldsBkgRatio[iCutSet-1]->GetXaxis()->SetTitleSize(0.11);
      hRawYieldsBkgRatio[iCutSet-1]->GetYaxis()->SetLabelSize(0.10);
      hRawYieldsBkgRatio[iCutSet-1]->GetXaxis()->SetLabelSize(0.10);
      hRawYieldsBkgRatio[iCutSet-1]->GetXaxis()->SetTitleOffset(1.2);
      hRawYieldsBkgRatio[iCutSet-1]->GetYaxis()->SetTitleOffset(0.6);
      hRawYieldsBkgRatio[iCutSet-1]->GetYaxis()->SetDecimals(2);
      
      hRawYieldsRelErrorRatio[iCutSet-1] = (TH1F*)hRawYieldsRelError[iCutSet]->Clone();
      hRawYieldsRelErrorRatio[iCutSet-1]->SetDirectory(0);
      hRawYieldsRelErrorRatio[iCutSet-1]->Divide(hRawYieldsRelError[iCutSet],hRawYieldsRelError[0],1.,1.,"");
      for(Int_t iPt=0; iPt<hRawYieldsRatio[iCutSet-1]->GetNbinsX(); iPt++) {
        hRawYieldsRelErrorRatio[iCutSet-1]->SetBinError(iPt+1,1.e-10);
      }
      hRawYieldsRelErrorRatio[iCutSet-1]->GetYaxis()->SetTitle("Ratio");
      hRawYieldsRelErrorRatio[iCutSet-1]->GetYaxis()->SetRangeUser(hRawYieldsRelErrorRatio[0]->GetMinimum()*0.5,hRawYieldsRelErrorRatio[0]->GetMaximum()*1.5);
      hRawYieldsRelErrorRatio[iCutSet-1]->GetYaxis()->SetTitleSize(0.11);
      hRawYieldsRelErrorRatio[iCutSet-1]->GetXaxis()->SetTitleSize(0.11);
      hRawYieldsRelErrorRatio[iCutSet-1]->GetYaxis()->SetLabelSize(0.10);
      hRawYieldsRelErrorRatio[iCutSet-1]->GetXaxis()->SetLabelSize(0.10);
      hRawYieldsRelErrorRatio[iCutSet-1]->GetXaxis()->SetTitleOffset(1.2);
      hRawYieldsRelErrorRatio[iCutSet-1]->GetYaxis()->SetTitleOffset(0.6);
      hRawYieldsRelErrorRatio[iCutSet-1]->GetYaxis()->SetDecimals(2);
    }
  }
  
  TLegend* l = new TLegend(0.5,0.6,0.89,0.85);
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.06);
  TLegend* l2 = new TLegend(0.22,0.6,0.4,0.85);
  l2->SetFillStyle(0);
  l2->SetFillColor(0);
  l2->SetBorderSize(0);
  l2->SetTextSize(0.06);
  TLegend* l3 = new TLegend(0.22,0.6,0.4,0.85);
  l3->SetFillColor(0);
  l3->SetFillStyle(0);
  l3->SetBorderSize(0);
  l3->SetTextSize(0.06);
  TLegend* l4 = new TLegend(0.35,0.7,0.89,0.89);
  l4->SetFillColor(0);
  l4->SetFillStyle(0);
  l4->SetBorderSize(0);
  l4->SetTextSize(0.045);
  for(Int_t iCutSet=0; iCutSet<nCutSets; iCutSet++) {
    l->AddEntry(hRawYields[iCutSet],legendnames[iCutSet].Data(),"lpe");
    l2->AddEntry(hRawYields[iCutSet],legendnames[iCutSet].Data(),"lpe");
    l3->AddEntry(hRawYields[iCutSet],legendnames[iCutSet].Data(),"lpe");
    l4->AddEntry(hRawYields[iCutSet],legendnames[iCutSet].Data(),"lpe");
  }
  l3->AddEntry(line,"PDG value","l");

  TCanvas* cRawYields = new TCanvas("cRawYields","",10,10,1000,1000);
  cRawYields->cd();
  TPad *padRawYields1 = new TPad("padRawYields1","padRawYields1",0,0.35,1,1);
  padRawYields1->SetBottomMargin(0);
  padRawYields1->Draw();
  padRawYields1->cd();
  padRawYields1->SetLogx();
  hRawYields[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets; iCutSet++) {
    hRawYields[iCutSet]->Draw("Esame");
  }
  l->Draw("same");
  cRawYields->cd();
  TPad *padRawYields2 = new TPad("padRawYields2","padRawYields2",0.,0.,1,0.35);
  padRawYields2->SetTopMargin(0);
  padRawYields2->Draw();
  padRawYields2->cd();
  padRawYields2->SetLogx();
  hRawYieldsRatio[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets-1; iCutSet++) {
    hRawYieldsRatio[iCutSet]->Draw("Esame");
  }
  lineatone->Draw("same");
  
  TCanvas* cRawYieldsSigma = new TCanvas("cRawYieldsSigma","",10,10,1000,1000);
  cRawYieldsSigma->cd();
  TPad *padRawYieldsSigma1 = new TPad("padRawYieldsSigma1","padRawYieldsSigma1",0,0.35,1,1);
  padRawYieldsSigma1->SetBottomMargin(0);
  padRawYieldsSigma1->Draw();
  padRawYieldsSigma1->cd();
  padRawYieldsSigma1->SetLogx();
  hRawYieldsSigma[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets; iCutSet++) {
    hRawYieldsSigma[iCutSet]->Draw("Esame");
  }
  l2->Draw("same");
  cRawYieldsSigma->cd();
  TPad *padRawYieldsSigma2 = new TPad("padRawYieldsSigma2","padRawYieldsSigma2",0.,0.,1,0.35);
  padRawYieldsSigma2->SetTopMargin(0);
  padRawYieldsSigma2->Draw();
  padRawYieldsSigma2->cd();
  padRawYieldsSigma2->SetLogx();
  hRawYieldsSigmaRatio[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets-1; iCutSet++) {
    hRawYieldsSigmaRatio[iCutSet]->Draw("Esame");
  }
  lineatone->Draw("same");
  
  TCanvas* cRawYieldsMean = new TCanvas("cRawYieldsMean","",10,10,1000,1000);
  cRawYieldsMean->cd();
  TPad *padRawYieldsMean1 = new TPad("padRawYieldsMean1","padRawYieldsMean1",0,0.35,1,1);
  padRawYieldsMean1->SetBottomMargin(0);
  padRawYieldsMean1->Draw();
  padRawYieldsMean1->cd();
  padRawYieldsMean1->SetLogx();
  hRawYieldsMean[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets; iCutSet++) {
    hRawYieldsMean[iCutSet]->Draw("Esame");
  }
  l3->Draw("same");
  line->Draw("same");
  cRawYieldsMean->cd();
  TPad *padRawYieldsMean2 = new TPad("padRawYieldsMean2","padRawYieldsMean2",0.,0.,1,0.35);
  padRawYieldsMean2->SetTopMargin(0);
  padRawYieldsMean2->Draw();
  padRawYieldsMean2->cd();
  padRawYieldsMean2->SetLogx();
  hRawYieldsMeanRatio[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets-1; iCutSet++) {
    hRawYieldsMeanRatio[iCutSet]->Draw("Esame");
  }
  lineatone->Draw("same");
  
  TCanvas* cRawYieldsSignificance = new TCanvas("cRawYieldsSignificance","",10,10,1000,1000);
  cRawYieldsSignificance->cd();
  TPad *padRawYieldsSignificance1 = new TPad("padRawYieldsSignificance1","padRawYieldsSignificance1",0,0.35,1,1);
  padRawYieldsSignificance1->SetBottomMargin(0);
  padRawYieldsSignificance1->Draw();
  padRawYieldsSignificance1->cd();
  padRawYieldsSignificance1->SetLogx();
  hRawYieldsSignificance[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets; iCutSet++) {
    hRawYieldsSignificance[iCutSet]->Draw("Esame");
  }
  l->Draw("same");
  cRawYieldsSignificance->cd();
  TPad *padRawYieldsSignificance2 = new TPad("padRawYieldsSignificance2","padRawYieldsSignificance2",0.,0.,1,0.35);
  padRawYieldsSignificance2->SetTopMargin(0);
  padRawYieldsSignificance2->Draw();
  padRawYieldsSignificance2->cd();
  padRawYieldsSignificance2->SetLogx();
  hRawYieldsSignificanceRatio[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets-1; iCutSet++) {
    hRawYieldsSignificanceRatio[iCutSet]->Draw("Esame");
  }
  lineatone->Draw("same");

  TCanvas* cRawYieldsSoverB = new TCanvas("cRawYieldsSoverB","",10,10,1000,1000);
  cRawYieldsSoverB->cd();
  TPad *padRawYieldsSoverB1 = new TPad("padRawYieldsSoverB1","padRawYieldsSoverB1",0,0.35,1,1);
  padRawYieldsSoverB1->SetBottomMargin(0);
  padRawYieldsSoverB1->Draw();
  padRawYieldsSoverB1->cd();
  padRawYieldsSoverB1->SetLogx();
  hRawYieldsSoverB[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets; iCutSet++) {
    hRawYieldsSoverB[iCutSet]->Draw("Esame");
  }
  l2->Draw("same");
  cRawYieldsSoverB->cd();
  TPad *padRawYieldsSoverB2 = new TPad("padRawYieldsSoverB2","padRawYieldsSoverB2",0.,0.,1,0.35);
  padRawYieldsSoverB2->SetTopMargin(0);
  padRawYieldsSoverB2->Draw();
  padRawYieldsSoverB2->cd();
  padRawYieldsSoverB2->SetLogx();
  hRawYieldsSoverBRatio[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets-1; iCutSet++) {
    hRawYieldsSoverBRatio[iCutSet]->Draw("Esame");
  }
  lineatone->Draw("same");

  TCanvas* cRawYieldsSignal = new TCanvas("cRawYieldsSignal","",10,10,1000,1000);
  cRawYieldsSignal->cd();
  TPad *padRawYieldsSignal1 = new TPad("padRawYieldsSignal1","padRawYieldsSignal1",0,0.35,1,1);
  padRawYieldsSignal1->SetBottomMargin(0);
  padRawYieldsSignal1->Draw();
  padRawYieldsSignal1->cd();
  padRawYieldsSignal1->SetLogx();
  hRawYieldsSignal[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets; iCutSet++) {
    hRawYieldsSignal[iCutSet]->Draw("Esame");
  }
  l->Draw("same");
  cRawYieldsSignal->cd();
  TPad *padRawYieldsSignal2 = new TPad("padRawYieldsSignal2","padRawYieldsSignal2",0.,0.,1,0.35);
  padRawYieldsSignal2->SetTopMargin(0);
  padRawYieldsSignal2->Draw();
  padRawYieldsSignal2->cd();
  padRawYieldsSignal2->SetLogx();
  hRawYieldsSignalRatio[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets-1; iCutSet++) {
    hRawYieldsSignalRatio[iCutSet]->Draw("Esame");
  }
  lineatone->Draw("same");
  
  TCanvas* cRawYieldsBkg = new TCanvas("cRawYieldsBkg","",10,10,1000,1000);
  cRawYieldsBkg->cd();
  TPad *padRawYieldsBkg1 = new TPad("padRawYieldsBkg1","padRawYieldsBkg1",0,0.35,1,1);
  padRawYieldsBkg1->SetBottomMargin(0);
  padRawYieldsBkg1->Draw();
  padRawYieldsBkg1->cd();
  padRawYieldsBkg1->SetLogx();
  hRawYieldsBkg[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets; iCutSet++) {
    hRawYieldsBkg[iCutSet]->Draw("Esame");
  }
  l->Draw("same");
  cRawYieldsBkg->cd();
  TPad *padRawYieldsBkg2 = new TPad("padRawYieldsBkg2","padRawYieldsBkg2",0.,0.,1,0.35);
  padRawYieldsBkg2->SetTopMargin(0);
  padRawYieldsBkg2->Draw();
  padRawYieldsBkg2->cd();
  padRawYieldsBkg2->SetLogx();
  hRawYieldsBkgRatio[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets-1; iCutSet++) {
    hRawYieldsBkgRatio[iCutSet]->Draw("Esame");
  }
  lineatone->Draw("same");
  
  TCanvas* cChi = new TCanvas("cChi","",10,10,1000,1000);
  cChi->SetBottomMargin(0.12);
  cChi->SetLogx();
  hRawYieldsChiSquare[0]->Draw();
  for(Int_t iCutSet=1; iCutSet<nCutSets; iCutSet++) {
    hRawYieldsChiSquare[iCutSet]->Draw("same");
  }
  l4->Draw("same");

  TCanvas* cRawYieldsRelError = new TCanvas("cRawYieldsRelError","",10,10,1000,1000);
  cRawYieldsRelError->cd();
  TPad *padRawYieldsRelError1 = new TPad("padRawYieldsRelError1","padRawYieldsRelError1",0,0.35,1,1);
  padRawYieldsRelError1->SetBottomMargin(0);
  padRawYieldsRelError1->Draw();
  padRawYieldsRelError1->cd();
  padRawYieldsRelError1->SetLogx();
  hRawYieldsRelError[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets; iCutSet++) {
    hRawYieldsRelError[iCutSet]->Draw("Esame");
  }
  l->Draw("same");
  cRawYieldsRelError->cd();
  TPad *padRawYieldsRelError2 = new TPad("padRawYieldsRelError2","padRawYieldsRelError2",0.,0.,1,0.35);
  padRawYieldsRelError2->SetTopMargin(0);
  padRawYieldsRelError2->Draw();
  padRawYieldsRelError2->cd();
  padRawYieldsRelError2->SetLogx();
  hRawYieldsRelErrorRatio[0]->Draw("E");
  for(Int_t iCutSet=1; iCutSet<nCutSets-1; iCutSet++) {
    hRawYieldsRelErrorRatio[iCutSet]->Draw("Esame");
  }
  lineatone->Draw("same");
  
  cRawYields->SaveAs(Form("%sComparison_RawYields%s.pdf",outprepend.Data(),outpostpend.Data()));
  cRawYieldsSigma->SaveAs(Form("%sComparison_Sigma%s.pdf",outprepend.Data(),outpostpend.Data()));
  cRawYieldsMean->SaveAs(Form("%sComparison_Mean%s.pdf",outprepend.Data(),outpostpend.Data()));
  cRawYieldsSignificance->SaveAs(Form("%sComparison_Significance%s.pdf",outprepend.Data(),outpostpend.Data()));
  cRawYieldsSoverB->SaveAs(Form("%sComparison_SoverB%s.pdf",outprepend.Data(),outpostpend.Data()));
  cRawYieldsSignal->SaveAs(Form("%sComparison_Signal%s.pdf",outprepend.Data(),outpostpend.Data()));
  cRawYieldsBkg->SaveAs(Form("%sComparison_Bkg%s.pdf",outprepend.Data(),outpostpend.Data()));
  cRawYieldsRelError->SaveAs(Form("%sComparison_RawYieldRelErrors%s.pdf",outprepend.Data(),outpostpend.Data()));
  cChi->SaveAs(Form("%sComparison_ChiSquare%s.pdf",outprepend.Data(),outpostpend.Data()));
}

void MassSpectraComparison(TString filename) {
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.3);
  
  Int_t nPtBins;
  TH1F*** hMass = new TH1F**[nCutSets];
  TH1F*** hMassRatio = new TH1F**[nCutSets-1];
  for(Int_t iCutSet=0; iCutSet<nCutSets; iCutSet++) {
    TFile* infile = TFile::Open(Form("%s/%s%s.root",directoryname[iCutSet].Data(),filename.Data(),cutsetsname[iCutSet].Data()),"READ");
    TCanvas* c=(TCanvas*)infile->Get("cMass");
    TList* l=(TList*)c->GetListOfPrimitives();
    TH1F* hEv = (TH1F*)infile->Get("hEv");
    hEv->SetDirectory(0);
    TH1F* hRawYields = (TH1F*)infile->Get("hRawYields");
    hRawYields->SetDirectory(0);
    nPtBins = hRawYields->GetNbinsX();
    hMass[iCutSet] = new TH1F*[nPtBins];
    if(iCutSet>0) {
      hMassRatio[iCutSet-1] = new TH1F*[nPtBins];
    }
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      hMass[iCutSet][iPt] = (TH1F*)c->GetPad(iPt+1)->GetPrimitive("fhistoInvMass");
      hMass[iCutSet][iPt]->SetDirectory(0);
      hMass[iCutSet][iPt]->Sumw2();
      hMass[iCutSet][iPt]->Scale(1./hEv->GetBinContent(1));
      hMass[iCutSet][iPt]->SetLineColor(colors[iCutSet]);
      hMass[iCutSet][iPt]->SetMarkerColor(colors[iCutSet]);
      hMass[iCutSet][iPt]->SetMarkerStyle(styles[iCutSet]);
      hMass[iCutSet][iPt]->GetFunction("funcbkgFullRange")->SetBit(TF1::kNotDraw);
      hMass[iCutSet][iPt]->GetFunction("funcbkgRecalc")->SetBit(TF1::kNotDraw);
      hMass[iCutSet][iPt]->GetYaxis()->SetRangeUser(0.,hMass[iCutSet][iPt]->GetMaximum()*1.1);
      if(iCutSet>0) {
        hMassRatio[iCutSet-1][iPt] = (TH1F*)hMass[iCutSet][iPt]->Clone();
        hMassRatio[iCutSet-1][iPt]->SetDirectory(0);
        hMassRatio[iCutSet-1][iPt]->Divide(hMass[iCutSet][iPt],hMass[0][iPt],1.,1.);
        hMassRatio[iCutSet-1][iPt]->GetYaxis()->SetRangeUser(0.,2);
      }
    }
    infile->Close();
  }
  
  TCanvas *c2 = new TCanvas("c2","",10,10,2000,2000);
  c2->Divide(4,3);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    c2->cd(iPt+1);
    hMass[0][iPt]->Draw();
    for(Int_t iCutSet=1; iCutSet<nCutSets; iCutSet++) {
      hMass[iCutSet][iPt]->Draw("same");
    }
  }

  TCanvas *cRatio = new TCanvas("cRatio","",10,10,2000,2000);
  cRatio->Divide(4,3);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    cRatio->cd(iPt+1);
    hMassRatio[0][iPt]->Draw();
    for(Int_t iCutSet=1; iCutSet<nCutSets-1; iCutSet++) {
      hMassRatio[iCutSet][iPt]->Draw("same");
    }
  }
  
  c2->SaveAs(Form("%sComparison_MassSpectra%s.pdf",outprepend.Data(),outpostpend.Data()));
  cRatio->SaveAs(Form("%sComparison_MassSpectraRatio%s.pdf",outprepend.Data(),outpostpend.Data()));
  
}

/*
void EfficiencyComparison(TString filename) {
  
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.3);
  
  TFile file1(filename1,"READ");
  TH1F* hEffPrompt1 = (TH1F*)file1.Get("hEffPrompt");
  TH1F* hEffFD1 = (TH1F*)file1.Get("hEffFD");
  hEffPrompt1->SetDirectory(0);
  hEffFD1->SetDirectory(0);
  file1.Close();
  
  TFile file2(filename2,"READ");
  TH1F* hEffPrompt2 = (TH1F*)file2.Get("hEffPrompt");
  TH1F* hEffFD2 = (TH1F*)file2.Get("hEffFD");
  hEffPrompt2->SetDirectory(0);
  hEffFD2->SetDirectory(0);
  file2.Close();
  
  TLegend* l = new TLegend(0.55,0.25,0.85,0.5);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.06);
  l->AddEntry(hEffPrompt1,name1,"lpe");
  l->AddEntry(hEffPrompt2,name2,"lpe");
  
  TCanvas* cEffPrompt = new TCanvas("cEffPrompt","",10,10,1000,1000);
  cEffPrompt->cd();
  TPad *padPrompt1 = new TPad("padPrompt1","padPrompt1",0,0.35,1,1);
  padPrompt1->SetBottomMargin(0);
  padPrompt1->Draw();
  padPrompt1->SetLogy();
  padPrompt1->cd();
  cEffPrompt->SetLogy();
  hEffPrompt1->GetYaxis()->SetTitle("#epsilon_{prompt}");
  hEffPrompt1->GetYaxis()->SetTitleSize(0.08);
  hEffPrompt1->GetYaxis()->SetLabelSize(0.06);
  hEffPrompt1->GetYaxis()->SetTitleOffset(0.8);
  hEffPrompt1->SetLineColor(kBlue);
  hEffPrompt1->SetMarkerColor(kBlue);
  hEffPrompt2->SetLineColor(kRed);
  hEffPrompt2->SetMarkerColor(kRed);
  hEffPrompt2->SetMarkerStyle(21);
  hEffPrompt1->Draw();
  hEffPrompt2->Draw("same");
  l->Draw("same");
  cEffPrompt->cd();
  TPad *padPrompt2 = new TPad("padPrompt2","padPrompt2",0.,0.,1,0.35);
  padPrompt2->SetTopMargin(0);
  padPrompt2->Draw();
  padPrompt2->cd();
  TH1F* hEffPromptRatio = (TH1F*)hEffPrompt2->Clone();
  hEffPromptRatio->SetDirectory(0);
  hEffPromptRatio->Divide(hEffPrompt2,hEffPrompt1,1.,1.,"");
  for(Int_t iPt=0; iPt<hEffPromptRatio->GetNbinsX(); iPt++) {
    hEffPromptRatio->SetBinError(iPt+1,1.e-10);
  }
  hEffPromptRatio->GetYaxis()->SetTitle("Ratio");
  hEffPromptRatio->GetYaxis()->SetTitleSize(0.10);
  hEffPromptRatio->GetXaxis()->SetTitleSize(0.10);
  hEffPromptRatio->GetYaxis()->SetLabelSize(0.10);
  hEffPromptRatio->GetXaxis()->SetLabelSize(0.10);
  hEffPromptRatio->GetXaxis()->SetTitleOffset(1.2);
  hEffPromptRatio->GetYaxis()->SetTitleOffset(0.6);
  hEffPromptRatio->GetYaxis()->SetDecimals(2);
  hEffPromptRatio->Draw();

  
  TCanvas* cEffFD = new TCanvas("cEffFD","",10,10,1000,1000);
  cEffFD->cd();
  TPad *padFD1 = new TPad("padFD1","padFD1",0,0.35,1,1);
  padFD1->SetBottomMargin(0);
  padFD1->Draw();
  padFD1->SetLogy();
  padFD1->cd();
  cEffFD->SetLogy();
  hEffFD1->GetYaxis()->SetTitle("#epsilon_{feed-down}");
  hEffFD1->GetYaxis()->SetTitleSize(0.08);
  hEffFD1->GetYaxis()->SetLabelSize(0.06);
  hEffFD1->GetYaxis()->SetTitleOffset(0.8);
  hEffFD1->SetLineColor(kBlue);
  hEffFD1->SetMarkerColor(kBlue);
  hEffFD2->SetLineColor(kRed);
  hEffFD2->SetMarkerColor(kRed);
  hEffFD2->SetMarkerStyle(21);
  hEffFD1->Draw();
  hEffFD2->Draw("same");
  l->Draw("same");
  cEffFD->cd();
  TPad *padFD2 = new TPad("padFD2","padFD2",0.,0.,1,0.35);
  padFD2->SetTopMargin(0);
  padFD2->Draw();
  padFD2->cd();
  TH1F* hEffFDRatio = (TH1F*)hEffFD2->Clone();
  hEffFDRatio->SetDirectory(0);
  hEffFDRatio->Divide(hEffFD2,hEffFD1,1.,1.,"");
  for(Int_t iPt=0; iPt<hEffFDRatio->GetNbinsX(); iPt++) {
    hEffFDRatio->SetBinError(iPt+1,1.e-10);
  }
  hEffFDRatio->GetYaxis()->SetTitle("Ratio");
  hEffFDRatio->GetYaxis()->SetTitleSize(0.10);
  hEffFDRatio->GetXaxis()->SetTitleSize(0.10);
  hEffFDRatio->GetYaxis()->SetLabelSize(0.10);
  hEffFDRatio->GetXaxis()->SetLabelSize(0.10);
  hEffFDRatio->GetXaxis()->SetTitleOffset(1.2);
  hEffFDRatio->GetYaxis()->SetTitleOffset(0.6);
  hEffFDRatio->GetYaxis()->SetDecimals(2);
  hEffFDRatio->Draw();
  
  cEffPrompt->SaveAs("Comparison_EffPrompt.eps");
  cEffFD->SaveAs("Comparison_EffFD.eps");
  
}

void FPromptComparison(TString filename) {
  
  
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.3);
  
  TFile file1(filename1,"READ");
  TH1F* hFPrompt1 = (TH1F*)file1.Get("hFrac");
  hFPrompt1->SetDirectory(0);
  file1.Close();
  
  TFile file2(filename2,"READ");
  TH1F* hFPrompt2 = (TH1F*)file2.Get("hFrac");
  hFPrompt2->SetDirectory(0);
  file2.Close();
  
  
  TLegend* l = new TLegend(0.59,0.05,0.89,0.3);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.06);
  l->AddEntry(hFPrompt1,name1,"lpe");
  l->AddEntry(hFPrompt2,name2,"lpe");
  
  TCanvas* cFPrompt = new TCanvas("cFPrompt","",10,10,1000,1000);
  cFPrompt->cd();
  TPad *padPrompt1 = new TPad("padPrompt1","padPrompt1",0,0.35,1,1);
  padPrompt1->SetBottomMargin(0);
  padPrompt1->Draw();
  padPrompt1->cd();
  hFPrompt1->GetYaxis()->SetRangeUser(0.51,1.2);
  hFPrompt1->GetYaxis()->SetTitleSize(0.08);
  hFPrompt1->GetYaxis()->SetLabelSize(0.06);
  hFPrompt1->GetYaxis()->SetTitleOffset(0.8);
  hFPrompt1->SetLineColor(kBlue);
  hFPrompt1->SetMarkerColor(kBlue);
  hFPrompt2->SetLineColor(kRed);
  hFPrompt2->SetMarkerColor(kRed);
  hFPrompt2->SetMarkerStyle(21);
  hFPrompt1->Draw();
  hFPrompt2->Draw("same");
  l->Draw("same");
  cFPrompt->cd();
  TPad *padPrompt2 = new TPad("padPrompt2","padPrompt2",0.,0.,1,0.35);
  padPrompt2->SetTopMargin(0);
  padPrompt2->Draw();
  padPrompt2->cd();
  TH1F* hFPromptRatio = (TH1F*)hFPrompt2->Clone();
  hFPromptRatio->SetDirectory(0);
  hFPromptRatio->Divide(hFPrompt2,hFPrompt1,1.,1.,"");
  for(Int_t iPt=0; iPt<hFPromptRatio->GetNbinsX(); iPt++) {
    hFPromptRatio->SetBinError(iPt+1,1.e-10);
  }
  hFPromptRatio->GetYaxis()->SetTitle("Ratio");
  hFPromptRatio->GetYaxis()->SetTitleSize(0.10);
  hFPromptRatio->GetXaxis()->SetTitleSize(0.10);
  hFPromptRatio->GetYaxis()->SetLabelSize(0.10);
  hFPromptRatio->GetXaxis()->SetLabelSize(0.10);
  hFPromptRatio->GetXaxis()->SetTitleOffset(1.2);
  hFPromptRatio->GetYaxis()->SetTitleOffset(0.6);
  hFPromptRatio->GetYaxis()->SetDecimals(2);
  hFPromptRatio->Draw();
  
  cFPrompt->SaveAs("Comparison_Fprompt.eps");
  
}

void CrossSectionComparison(TString filename) {
  
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.3);
  
  TFile file1(filename1,"READ");
  TH1F* hCrossPrompt1 = (TH1F*)file1.Get("hPromptCrossSection");
  hCrossPrompt1->SetDirectory(0);
  file1.Close();
  
  TFile file2(filename2,"READ");
  TH1F* hCrossPrompt2 = (TH1F*)file2.Get("hPromptCrossSection");
  hCrossPrompt2->SetDirectory(0);
  file2.Close();
  
  TLegend* l = new TLegend(0.55,0.6,0.85,0.85);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.06);
  l->AddEntry(hCrossPrompt1,name1,"lpe");
  l->AddEntry(hCrossPrompt2,name2,"lpe");
  
  TCanvas* cCrossPrompt = new TCanvas("cCrossPrompt","",10,10,1000,1000);
  cCrossPrompt->cd();
  TPad *padPrompt1 = new TPad("padPrompt1","padPrompt1",0,0.35,1,1);
  padPrompt1->SetBottomMargin(0);
  padPrompt1->Draw();
  padPrompt1->cd();
  padPrompt1->SetLogy();
  hCrossPrompt1->GetYaxis()->SetRangeUser(hCrossPrompt1->GetMinimum()*0.7,hCrossPrompt1->GetMaximum());
  hCrossPrompt1->GetYaxis()->SetTitleSize(0.08);
  hCrossPrompt1->GetYaxis()->SetLabelSize(0.06);
  hCrossPrompt1->GetYaxis()->SetTitleOffset(0.8);
  hCrossPrompt1->SetLineColor(kBlue);
  hCrossPrompt1->SetMarkerColor(kBlue);
  hCrossPrompt2->SetLineColor(kRed);
  hCrossPrompt2->SetMarkerColor(kRed);
  hCrossPrompt2->SetMarkerStyle(21);
  hCrossPrompt1->Draw();
  hCrossPrompt2->Draw("same");
  l->Draw("same");
  cCrossPrompt->cd();
  TPad *padPrompt2 = new TPad("padPrompt2","padPrompt2",0.,0.,1,0.35);
  padPrompt2->SetTopMargin(0);
  padPrompt2->Draw();
  padPrompt2->cd();
  TH1F* hCrossPromptRatio = (TH1F*)hCrossPrompt2->Clone();
  hCrossPromptRatio->SetDirectory(0);
  hCrossPromptRatio->Divide(hCrossPrompt2,hCrossPrompt1,1.,1.,"");
  for(Int_t iPt=0; iPt<hCrossPromptRatio->GetNbinsX(); iPt++) {
    hCrossPromptRatio->SetBinError(iPt+1,1.e-10);
  }
  hCrossPromptRatio->GetYaxis()->SetTitle("Ratio");
  hCrossPromptRatio->GetYaxis()->SetTitleSize(0.10);
  hCrossPromptRatio->GetXaxis()->SetTitleSize(0.10);
  hCrossPromptRatio->GetYaxis()->SetLabelSize(0.10);
  hCrossPromptRatio->GetXaxis()->SetLabelSize(0.10);
  hCrossPromptRatio->GetXaxis()->SetTitleOffset(1.2);
  hCrossPromptRatio->GetYaxis()->SetTitleOffset(0.6);
  hCrossPromptRatio->GetYaxis()->SetDecimals(2);
  hCrossPromptRatio->Draw();
  
  cCrossPrompt->SaveAs("Comparison_PromptCrossSection.eps");
  
}*/

