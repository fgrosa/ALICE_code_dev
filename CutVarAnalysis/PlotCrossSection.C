#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TColor.h>
#include <TLine.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <TArrayD.h>
#include <TArrow.h>
#include <TLatex.h>
#include <vector>
#endif

const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

void PlotCrossSection() {

  TString system = "pPb";

  //________________________________________________________________________________________________________________
  //ptbins  
  const Int_t nPtBins = 7;
  const Int_t nPtLims = nPtBins+1;
  Double_t PtLims[nPtLims] = {2.,3.,4.,5.,6.,8.,12.,16.};

  //________________________________________________________________________________________________________________
  //systematic uncertainties

  Double_t CutVarSystPromptMinLow[nPtBins] = {0.09,0.04,0.04,0.06,0.06,0.04,0.04};
  Double_t CutVarSystPromptMinHigh[nPtBins] = {0.09,0.04,0.04,0.06,0.06,0.04,0.04};
  Double_t CutVarSystFDMinLow[nPtBins] = {0.26,0.28,0.10,0.20,0.11,0.11,0.12};
  Double_t CutVarSystFDMinHigh[nPtBins] = {0.26,0.28,0.10,0.20,0.11,0.11,0.12};
  Double_t CutVarSystPromptIncLow[nPtBins] = {0.08,0.07,0.04,0.05,0.04,0.04,0.05};
  Double_t CutVarSystPromptIncHigh[nPtBins] = {0.08,0.07,0.04,0.05,0.04,0.04,0.05};
  Double_t CutVarSystFDIncLow[nPtBins] = {0.25,0.23,0.09,0.21,0.10,0.11,0.12};
  Double_t CutVarSystFDIncHigh[nPtBins] = {0.25,0.23,0.09,0.21,0.10,0.11,0.12};

  Double_t RawSystPromptMinLow[nPtBins] = {0.12,0.08,0.08,0.08,0.09,0.09,0.15};
  Double_t RawSystPromptMinHigh[nPtBins] = {0.12,0.08,0.08,0.08,0.09,0.09,0.15};
  Double_t RawSystFDMinLow[nPtBins] = {0.26,0.20,0.20,0.16,0.10,0.11,0.20};
  Double_t RawSystFDMinHigh[nPtBins] = {0.26,0.20,0.20,0.16,0.10,0.11,0.20};
  Double_t RawSystPromptIncLow[nPtBins] = {0.12,0.08,0.08,0.08,0.10,0.12,0.15};
  Double_t RawSystPromptIncHigh[nPtBins] = {0.12,0.08,0.08,0.08,0.10,0.12,0.15};
  Double_t RawSystFDIncLow[nPtBins] = {0.30,0.20,0.20,0.18,0.11,0.11,0.20};
  Double_t RawSystFDIncHigh[nPtBins] = {0.30,0.20,0.20,0.18,0.11,0.11,0.20};

  Double_t EffSystPromptMinLow[nPtBins] = {0.01,0.01,0.01,0.01,0.01,0.02,0.02};
  Double_t EffSystPromptMinHigh[nPtBins] = {0.01,0.01,0.01,0.01,0.01,0.02,0.02};
  Double_t EffSystFDMinLow[nPtBins] = {0.01,0.01,0.01,0.01,0.01,0.05,0.05};
  Double_t EffSystFDMinHigh[nPtBins] = {0.01,0.01,0.01,0.01,0.01,0.05,0.05};  
  Double_t EffSystPromptIncLow[nPtBins] = {0.01,0.01,0.01,0.01,0.01,0.04,0.02};
  Double_t EffSystPromptIncHigh[nPtBins] = {0.01,0.01,0.01,0.01,0.01,0.04,0.02};
  Double_t EffSystFDIncLow[nPtBins] = {0.01,0.01,0.01,0.01,0.01,0.05,0.05};
  Double_t EffSystFDIncHigh[nPtBins] = {0.01,0.01,0.01,0.01,0.01,0.05,0.05};

  Double_t TrackSystPromptMinLow[nPtBins] = {0.09,0.09,0.09,0.09,0.09,0.09,0.09};
  Double_t TrackSystPromptMinHigh[nPtBins] = {0.09,0.09,0.09,0.09,0.09,0.09,0.09};
  Double_t TrackSystFDMinLow[nPtBins] = {0.09,0.09,0.09,0.09,0.09,0.09,0.09};
  Double_t TrackSystFDMinHigh[nPtBins] = {0.09,0.09,0.09,0.09,0.09,0.09,0.09};  
  Double_t TrackSystPromptIncLow[nPtBins] = {0.09,0.09,0.09,0.09,0.09,0.09,0.09};
  Double_t TrackSystPromptIncHigh[nPtBins] = {0.09,0.09,0.09,0.09,0.09,0.09,0.09};
  Double_t TrackSystFDIncLow[nPtBins] = {0.09,0.09,0.09,0.09,0.09,0.09,0.09};
  Double_t TrackSystFDIncHigh[nPtBins] = {0.09,0.09,0.09,0.09,0.09,0.09,0.09};
  
  //________________________________________________________________________________________________________________
  //input files
  
  //FONLL Cross sections for comparison
  TFile fileFONLLPrompt("../FONLL/FONLL5TeVpPbPrompt.root","READ");
  TGraphAsymmErrors* gFONLLPrompt =(TGraphAsymmErrors*)fileFONLLPrompt.Get("gFONLLPrompt");
  fileFONLLPrompt.Close();
  TFile fileFONLLFD("../FONLL/FONLL5TeVpPbFD.root","READ");
  TGraphAsymmErrors* gFONLLFD =(TGraphAsymmErrors*)fileFONLLFD.Get("gFONLLFD");
  fileFONLLFD.Close();
  
  //Published Cross sections for comparison (only prompt)
  TFile filePubPrompt("../Published/DplusCrossSec_method2_fd2_br1.root","READ");
  TH1F* hPubPrompt = (TH1F*)filePubPrompt.Get("hAAC");
  TGraphAsymmErrors* gPubSysPrompt =(TGraphAsymmErrors*)filePubPrompt.Get("gaaCsystTot");
  hPubPrompt->SetDirectory(0);
  filePubPrompt.Close();
  Double_t* puberrhigh = gPubSysPrompt->GetEYhigh();
  Double_t* puberrlow = gPubSysPrompt->GetEYlow();
  Double_t puberrhigh2[10];
  Double_t puberrlow2[10];
  
  for(Int_t iPt=0; iPt<10; iPt++) {
    puberrhigh2[iPt] = puberrhigh[iPt]/hPubPrompt->GetBinContent(iPt+1);
    puberrhigh2[iPt] = puberrhigh2[iPt]*puberrhigh2[iPt]-0.09*0.09;
    puberrhigh2[iPt] = TMath::Sqrt(puberrhigh2[iPt]);
    puberrlow2[iPt] = puberrlow[iPt]/hPubPrompt->GetBinContent(iPt+1);
    puberrlow2[iPt] = puberrlow2[iPt]*puberrlow2[iPt]-0.09*0.09;
    puberrlow2[iPt] = TMath::Sqrt(puberrlow2[iPt]);
  }
  
  Double_t puberrhighreb[7] = {puberrhigh2[1],puberrhigh2[2],puberrhigh2[3],puberrhigh2[4],(puberrhigh2[5]+puberrhigh2[6])/2,puberrhigh2[7],puberrhigh2[8]};
  Double_t puberrlowreb[7] = {puberrlow2[1],puberrlow2[2],puberrlow2[3],puberrlow2[4],(puberrlow2[5]+puberrlow2[6])/2,puberrlow2[7],puberrlow2[8]};

  //Cut variation files
  TFile cutvarfileMin("CrossSec_Min.root","READ");
  TH1F* hCrossPromptMin = (TH1F*)cutvarfileMin.Get("hCrossSecPrompt");
  TH1F* hCrossFDMin = (TH1F*)cutvarfileMin.Get("hCrossSecFD");
  hCrossPromptMin->SetDirectory(0);
  hCrossFDMin->SetDirectory(0);
  cutvarfileMin.Close();

  TFile cutvarfileInc("CrossSec_Inc.root","READ");
  TH1F* hCrossPromptInc = (TH1F*)cutvarfileInc.Get("hCrossSecPrompt");
  TH1F* hCrossFDInc = (TH1F*)cutvarfileInc.Get("hCrossSecFD");
  hCrossPromptInc->SetDirectory(0);
  hCrossFDInc->SetDirectory(0);
  cutvarfileInc.Close();

  //______________________________________________________________________________________
  //systematic uncertainties tgraphasymmerrors
  TGraphAsymmErrors* gRawSystPromptMin = new TGraphAsymmErrors(nPtBins);
  gRawSystPromptMin->SetLineColor(fillColors[0]);
  gRawSystPromptMin->SetFillStyle(20);
  gRawSystPromptMin->SetName("gRawSystPromptMin");
  TGraphAsymmErrors* gRawSystFDMin = new TGraphAsymmErrors(nPtBins);
  gRawSystFDMin->SetLineColor(fillColors[0]);
  gRawSystFDMin->SetFillStyle(20);
  gRawSystFDMin->SetName("gRawSystFDMin");
  TGraphAsymmErrors* gRawSystPromptInc = new TGraphAsymmErrors(nPtBins);
  gRawSystPromptInc->SetLineColor(fillColors[0]);
  gRawSystPromptInc->SetFillStyle(20);
  gRawSystPromptInc->SetName("gRawSystPromptInc");
  TGraphAsymmErrors* gRawSystFDInc = new TGraphAsymmErrors(nPtBins);
  gRawSystFDInc->SetLineColor(fillColors[0]);
  gRawSystFDInc->SetFillStyle(20);
  gRawSystFDInc->SetName("gRawSystFDInc");

  TGraphAsymmErrors* gCutSystPromptMin = new TGraphAsymmErrors(nPtBins);
  gCutSystPromptMin->SetLineColor(fillColors[1]);
  gCutSystPromptMin->SetFillStyle(20);
  gCutSystPromptMin->SetName("gCutSystPromptMin");
  TGraphAsymmErrors* gCutSystFDMin = new TGraphAsymmErrors(nPtBins);
  gCutSystFDMin->SetLineColor(fillColors[1]);
  gCutSystFDMin->SetFillStyle(20);
  gCutSystFDMin->SetName("gCutSystFDMin");
  TGraphAsymmErrors* gCutSystPromptInc = new TGraphAsymmErrors(nPtBins);
  gCutSystPromptInc->SetLineColor(fillColors[1]);
  gCutSystPromptInc->SetFillStyle(20);
  gCutSystPromptInc->SetName("gCutSystPromptInc");
  TGraphAsymmErrors* gCutSystFDInc = new TGraphAsymmErrors(nPtBins);
  gCutSystFDInc->SetLineColor(fillColors[1]);
  gCutSystFDInc->SetFillStyle(20);
  gCutSystFDInc->SetName("gCutSystFDInc");
  
  TGraphAsymmErrors* gEffSystPromptMin = new TGraphAsymmErrors(nPtBins);
  gEffSystPromptMin->SetLineColor(fillColors[2]);
  gEffSystPromptMin->SetFillStyle(20);
  gEffSystPromptMin->SetName("gEffSystPromptMin");
  TGraphAsymmErrors* gEffSystFDMin = new TGraphAsymmErrors(nPtBins);
  gEffSystFDMin->SetLineColor(fillColors[2]);
  gEffSystFDMin->SetFillStyle(20);
  gEffSystFDMin->SetName("gEffSystFDMin");
  TGraphAsymmErrors* gEffSystPromptInc = new TGraphAsymmErrors(nPtBins);
  gEffSystPromptInc->SetLineColor(fillColors[2]);
  gEffSystPromptInc->SetFillStyle(20);
  gEffSystPromptInc->SetName("gEffSystPromptInc");
  TGraphAsymmErrors* gEffSystFDInc = new TGraphAsymmErrors(nPtBins);
  gEffSystFDInc->SetLineColor(fillColors[2]);
  gEffSystFDInc->SetFillStyle(20);
  gEffSystFDInc->SetName("gEffSystFDInc");

  TGraphAsymmErrors* gTrackSystPromptMin = new TGraphAsymmErrors(nPtBins);
  gTrackSystPromptMin->SetLineColor(fillColors[2]);
  gTrackSystPromptMin->SetFillStyle(20);
  gTrackSystPromptMin->SetName("gTrackSystPromptMin");
  TGraphAsymmErrors* gTrackSystFDMin = new TGraphAsymmErrors(nPtBins);
  gTrackSystFDMin->SetLineColor(fillColors[2]);
  gTrackSystFDMin->SetFillStyle(20);
  gTrackSystFDMin->SetName("gTrackSystFDMin");
  TGraphAsymmErrors* gTrackSystPromptInc = new TGraphAsymmErrors(nPtBins);
  gTrackSystPromptInc->SetLineColor(fillColors[2]);
  gTrackSystPromptInc->SetFillStyle(20);
  gTrackSystPromptInc->SetName("gTrackSystPromptInc");
  TGraphAsymmErrors* gTrackSystFDInc = new TGraphAsymmErrors(nPtBins);
  gTrackSystFDInc->SetLineColor(fillColors[2]);
  gTrackSystFDInc->SetFillStyle(20);
  gTrackSystFDInc->SetName("gTrackSystFDInc");
  
  TGraphAsymmErrors* gSystPromptMin = new TGraphAsymmErrors(nPtBins);
  gSystPromptMin->SetLineColor(kBlack);
  gSystPromptMin->SetFillStyle(20);
  gSystPromptMin->SetName("gSystPromptMin");
  TGraphAsymmErrors* gSystFDMin = new TGraphAsymmErrors(nPtBins);
  gSystFDMin->SetLineColor(kBlack);
  gSystFDMin->SetFillStyle(20);
  gSystFDMin->SetName("gSystFDMin");
  TGraphAsymmErrors* gSystPromptInc = new TGraphAsymmErrors(nPtBins);
  gSystPromptInc->SetLineColor(kBlack);
  gSystPromptInc->SetFillStyle(20);
  gSystPromptInc->SetName("gSystPromptInc");
  TGraphAsymmErrors* gSystFDInc = new TGraphAsymmErrors(nPtBins);
  gSystFDInc->SetLineColor(kBlack);
  gSystFDInc->SetFillStyle(20);
  gSystFDInc->SetName("gSystFDInc");

  TGraphAsymmErrors* gCrossSystPromptMin = new TGraphAsymmErrors(nPtBins);
  gCrossSystPromptMin->SetLineColor(colors[2]);
  gCrossSystPromptMin->SetFillStyle(20);
  gCrossSystPromptMin->SetName("gCrossSystPromptMin");
  TGraphAsymmErrors* gCrossSystFDMin = new TGraphAsymmErrors(nPtBins);
  gCrossSystFDMin->SetLineColor(colors[2]);
  gCrossSystFDMin->SetFillStyle(20);
  gCrossSystFDMin->SetName("gCrossSystFDMin");
  TGraphAsymmErrors* gCrossSystPromptInc = new TGraphAsymmErrors(nPtBins);
  gCrossSystPromptInc->SetLineColor(colors[1]);
  gCrossSystPromptInc->SetFillStyle(20);
  gCrossSystPromptInc->SetName("gCrossSystPromptInc");
  TGraphAsymmErrors* gCrossSystFDInc = new TGraphAsymmErrors(nPtBins);
  gCrossSystFDInc->SetLineColor(colors[1]);
  gCrossSystFDInc->SetFillStyle(20);
  gCrossSystFDInc->SetName("gCrossSystFDInc");
   
  TGraphAsymmErrors* gSystRatioMin = new TGraphAsymmErrors(nPtBins);
  gSystRatioMin->SetLineColor(colors[1]);
  gSystRatioMin->SetFillStyle(20);
  gSystRatioMin->SetName("gSystRatioMin");
  TGraphAsymmErrors* gSystRatioInc = new TGraphAsymmErrors(nPtBins);
  gSystRatioInc->SetLineColor(colors[2]);
  gSystRatioInc->SetFillStyle(20);
  gSystRatioInc->SetName("gSystRatioInc");
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    Double_t ptcent = (PtLims[iPt]+PtLims[iPt+1])/2;

    gRawSystPromptMin->SetPoint(iPt,ptcent,0);
    gRawSystPromptMin->SetPointError(iPt,0.15,0.15,RawSystPromptMinLow[iPt],RawSystPromptMinHigh[iPt]);
    gRawSystFDMin->SetPoint(iPt,ptcent,0);
    gRawSystFDMin->SetPointError(iPt,0.15,0.15,RawSystFDMinLow[iPt],RawSystFDMinHigh[iPt]);
    gRawSystPromptInc->SetPoint(iPt,ptcent,0);
    gRawSystPromptInc->SetPointError(iPt,0.15,0.15,RawSystPromptIncLow[iPt],RawSystPromptIncHigh[iPt]);
    gRawSystFDInc->SetPoint(iPt,ptcent,0);
    gRawSystFDInc->SetPointError(iPt,0.15,0.15,RawSystFDIncLow[iPt],RawSystFDIncHigh[iPt]);

    gCutSystPromptMin->SetPoint(iPt,ptcent,0);
    gCutSystPromptMin->SetPointError(iPt,0.15,0.15,CutVarSystPromptMinLow[iPt],CutVarSystPromptMinHigh[iPt]);
    gCutSystFDMin->SetPoint(iPt,ptcent,0);
    gCutSystFDMin->SetPointError(iPt,0.15,0.15,CutVarSystFDMinLow[iPt],CutVarSystFDMinHigh[iPt]);
    gCutSystPromptInc->SetPoint(iPt,ptcent,0);
    gCutSystPromptInc->SetPointError(iPt,0.15,0.15,CutVarSystPromptIncLow[iPt],CutVarSystPromptIncHigh[iPt]);
    gCutSystFDInc->SetPoint(iPt,ptcent,0);
    gCutSystFDInc->SetPointError(iPt,0.15,0.15,CutVarSystFDIncLow[iPt],CutVarSystFDIncHigh[iPt]);

    gEffSystPromptMin->SetPoint(iPt,ptcent,0);
    gEffSystPromptMin->SetPointError(iPt,0.15,0.15,EffSystPromptMinLow[iPt],EffSystPromptMinHigh[iPt]);
    gEffSystFDMin->SetPoint(iPt,ptcent,0);
    gEffSystFDMin->SetPointError(iPt,0.15,0.15,EffSystFDMinLow[iPt],EffSystFDMinHigh[iPt]);
    gEffSystPromptInc->SetPoint(iPt,ptcent,0);
    gEffSystPromptInc->SetPointError(iPt,0.15,0.15,EffSystPromptIncLow[iPt],EffSystPromptIncHigh[iPt]);
    gEffSystFDInc->SetPoint(iPt,ptcent,0);
    gEffSystFDInc->SetPointError(iPt,0.15,0.15,EffSystFDIncLow[iPt],EffSystFDIncHigh[iPt]);

    gTrackSystPromptMin->SetPoint(iPt,ptcent,0);
    gTrackSystPromptMin->SetPointError(iPt,0.15,0.15,TrackSystPromptMinLow[iPt],TrackSystPromptMinHigh[iPt]);
    gTrackSystFDMin->SetPoint(iPt,ptcent,0);
    gTrackSystFDMin->SetPointError(iPt,0.15,0.15,TrackSystFDMinLow[iPt],TrackSystFDMinHigh[iPt]);
    gTrackSystPromptInc->SetPoint(iPt,ptcent,0);
    gTrackSystPromptInc->SetPointError(iPt,0.15,0.15,TrackSystPromptIncLow[iPt],TrackSystPromptIncHigh[iPt]);
    gTrackSystFDInc->SetPoint(iPt,ptcent,0);
    gTrackSystFDInc->SetPointError(iPt,0.15,0.15,TrackSystFDIncLow[iPt],TrackSystFDIncHigh[iPt]);

    Double_t SystPromptMinLow = TMath::Sqrt(RawSystPromptMinLow[iPt]*RawSystPromptMinLow[iPt]+CutVarSystPromptMinLow[iPt]*CutVarSystPromptMinLow[iPt]+EffSystPromptMinLow[iPt]*EffSystPromptMinLow[iPt]+TrackSystPromptMinLow[iPt]*TrackSystPromptMinLow[iPt]);
    Double_t SystPromptMinHigh = TMath::Sqrt(RawSystPromptMinHigh[iPt]*RawSystPromptMinHigh[iPt]+CutVarSystPromptMinHigh[iPt]*CutVarSystPromptMinHigh[iPt]+EffSystPromptMinHigh[iPt]*EffSystPromptMinHigh[iPt]+TrackSystPromptMinHigh[iPt]*TrackSystPromptMinHigh[iPt]);
    Double_t SystFDMinLow = TMath::Sqrt(RawSystFDMinLow[iPt]*RawSystFDMinLow[iPt]+CutVarSystFDMinLow[iPt]*CutVarSystFDMinLow[iPt]+EffSystFDMinLow[iPt]*EffSystFDMinLow[iPt]+TrackSystFDMinLow[iPt]*TrackSystFDMinLow[iPt]);
    Double_t SystFDMinHigh = TMath::Sqrt(RawSystFDMinHigh[iPt]*RawSystFDMinHigh[iPt]+CutVarSystFDMinHigh[iPt]*CutVarSystFDMinHigh[iPt]+EffSystFDMinHigh[iPt]*EffSystFDMinHigh[iPt]+TrackSystFDMinHigh[iPt]*TrackSystFDMinHigh[iPt]);
    Double_t SystPromptIncLow = TMath::Sqrt(RawSystPromptIncLow[iPt]*RawSystPromptIncLow[iPt]+CutVarSystPromptIncLow[iPt]*CutVarSystPromptIncLow[iPt]+EffSystPromptIncLow[iPt]*EffSystPromptIncLow[iPt]+TrackSystPromptIncLow[iPt]*TrackSystPromptIncLow[iPt]);
    Double_t SystPromptIncHigh = TMath::Sqrt(RawSystPromptIncHigh[iPt]*RawSystPromptIncHigh[iPt]+CutVarSystPromptIncHigh[iPt]*CutVarSystPromptIncHigh[iPt]+EffSystPromptIncHigh[iPt]*EffSystPromptIncHigh[iPt]+TrackSystPromptIncHigh[iPt]*TrackSystPromptIncHigh[iPt]);
    Double_t SystFDIncLow = TMath::Sqrt(RawSystFDIncLow[iPt]*RawSystFDIncLow[iPt]+CutVarSystFDIncLow[iPt]*CutVarSystFDIncLow[iPt]+EffSystFDIncLow[iPt]*EffSystFDIncLow[iPt]+TrackSystFDIncLow[iPt]*TrackSystFDIncLow[iPt]);
    Double_t SystFDIncHigh = TMath::Sqrt(RawSystFDIncHigh[iPt]*RawSystFDIncHigh[iPt]+CutVarSystFDIncHigh[iPt]*CutVarSystFDIncHigh[iPt]+EffSystFDIncHigh[iPt]*EffSystFDIncHigh[iPt]+TrackSystFDIncHigh[iPt]*TrackSystFDIncHigh[iPt]);

    gSystPromptMin->SetPoint(iPt,ptcent,0);
    gSystPromptMin->SetPointError(iPt,0.15,0.15,SystPromptMinLow,SystPromptMinHigh);
    gSystFDMin->SetPoint(iPt,ptcent,0);
    gSystFDMin->SetPointError(iPt,0.15,0.15,SystFDMinLow,SystFDMinHigh);
    gSystPromptInc->SetPoint(iPt,ptcent,0);
    gSystPromptInc->SetPointError(iPt,0.15,0.15,SystPromptIncLow,SystPromptIncHigh);
    gSystFDInc->SetPoint(iPt,ptcent,0);
    gSystFDInc->SetPointError(iPt,0.15,0.15,SystFDIncLow,SystFDIncHigh);

    gCrossSystPromptMin->SetPoint(iPt,ptcent,hCrossPromptMin->GetBinContent(iPt+1));
    gCrossSystPromptMin->SetPointError(iPt,0.15,0.15,SystPromptMinLow*hCrossPromptMin->GetBinContent(iPt+1),SystPromptMinHigh*hCrossPromptMin->GetBinContent(iPt+1));
    gCrossSystFDMin->SetPoint(iPt,ptcent,hCrossFDMin->GetBinContent(iPt+1));
    gCrossSystFDMin->SetPointError(iPt,0.15,0.15,SystFDMinLow*hCrossFDMin->GetBinContent(iPt+1),SystFDMinHigh*hCrossFDMin->GetBinContent(iPt+1));
    gCrossSystPromptInc->SetPoint(iPt,ptcent,hCrossPromptInc->GetBinContent(iPt+1));
    gCrossSystPromptInc->SetPointError(iPt,0.15,0.15,SystPromptIncLow*hCrossPromptInc->GetBinContent(iPt+1),SystPromptIncHigh*hCrossPromptInc->GetBinContent(iPt+1));
    gCrossSystFDInc->SetPoint(iPt,ptcent,hCrossFDInc->GetBinContent(iPt+1));
    gCrossSystFDInc->SetPointError(iPt,0.15,0.15,SystFDIncLow*hCrossFDInc->GetBinContent(iPt+1),SystFDIncHigh*hCrossFDInc->GetBinContent(iPt+1));
    
  }
  
  //______________________________________________________________________________________
  //calculate ratios between methods (only stat error)
  TH1F* hRatioMethodsPrompt = new TH1F("hRatioMethodsPrompt","",nPtBins,PtLims);
  hRatioMethodsPrompt->Divide(hCrossPromptMin,hCrossPromptInc,1.,1.);

  TH1F* hRatioMethodsFD = new TH1F("hRatioMethodsFD","",nPtBins,PtLims);
  hRatioMethodsFD->Divide(hCrossFDMin,hCrossFDInc,1.,1.);

  //______________________________________________________________________________________
  //calculate ratios w.r.t. the published cross section (only prompt) 
  TAxis* PtAxis = (TAxis*)hCrossPromptMin->GetXaxis();
  TArrayD* ptarray = (TArrayD*)PtAxis->GetXbins();

  const Int_t nPtBinsPub = hPubPrompt->GetNbinsX();
  const Int_t nPtLimsPub = nPtBins+1;
  TAxis* PtAxisPub = (TAxis*)hPubPrompt->GetXaxis();
  TArrayD* ptarrayPub = (TArrayD*)PtAxisPub->GetXbins();
  Double_t* PtLimsPub = (Double_t*)ptarrayPub->GetArray();

  Bool_t IsBinningEqual = kTRUE;

  if(nPtBins!=nPtBinsPub) 
    IsBinningEqual = kFALSE;

  Int_t Ptcounter=0;
  while(Ptcounter<nPtBins && IsBinningEqual==kTRUE) {
    if(PtLims[Ptcounter]!=PtLimsPub[Ptcounter])
      IsBinningEqual = kFALSE;
    Ptcounter++;
  }
  
  TH1F* hRatioPromptMin = new TH1F("hRatioPromptMin","",nPtBins,PtLims);
  TH1F* hRatioPromptInc = new TH1F("hRatioPromptInc","",nPtBins,PtLims);
  
  if(IsBinningEqual) {
    hRatioPromptMin->Divide(hCrossPromptMin,hPubPrompt,1.,1.);
    hRatioPromptInc->Divide(hCrossPromptInc,hPubPrompt,1.,1.);
  }
  else {
    TH1F* hPubPrompt2 = (TH1F*)hPubPrompt->Clone();
    for(Int_t iPt=0; iPt<nPtBinsPub; iPt++) {
      hPubPrompt2->SetBinContent(iPt+1,hPubPrompt2->GetBinContent(iPt+1)*(PtLimsPub[iPt+1]-PtLimsPub[iPt]));
      hPubPrompt2->SetBinError(iPt+1,hPubPrompt2->GetBinError(iPt+1)*(PtLimsPub[iPt+1]-PtLimsPub[iPt]));
    }
    TH1F* hPubPromptReb = (TH1F*)hPubPrompt2->Rebin(nPtBins,"hPubPromptReb",PtLims);
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      hPubPromptReb->SetBinContent(iPt+1,hPubPromptReb->GetBinContent(iPt+1)/(PtLims[iPt+1]-PtLims[iPt]));
      hPubPromptReb->SetBinError(iPt+1,hPubPromptReb->GetBinError(iPt+1)/(PtLims[iPt+1]-PtLims[iPt]));
      cout << hPubPromptReb->GetBinContent(iPt+1) << "    " << hPubPromptReb->GetBinError(iPt+1) << "   " <<hCrossPromptMin->GetBinContent(iPt+1) << "    " << hCrossPromptMin->GetBinError(iPt+1) <<endl;
    }
    hRatioPromptMin->Divide(hCrossPromptMin,hPubPromptReb,1.,1.);
    hRatioPromptInc->Divide(hCrossPromptInc,hPubPromptReb,1.,1.);    
  }

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    Double_t ptcent = (PtLims[iPt]+PtLims[iPt+1])/2;

    Double_t minerrnotrack = TMath::Sqrt(RawSystPromptMinLow[iPt]*RawSystPromptMinLow[iPt]+CutVarSystPromptMinLow[iPt]*CutVarSystPromptMinLow[iPt]+EffSystPromptMinLow[iPt]*EffSystPromptMinLow[iPt]);
    Double_t incerrnotrack = TMath::Sqrt(RawSystPromptIncLow[iPt]*RawSystPromptIncLow[iPt]+CutVarSystPromptIncLow[iPt]*CutVarSystPromptIncLow[iPt]+EffSystPromptIncLow[iPt]*EffSystPromptIncLow[iPt]);

    gSystRatioMin->SetPoint(iPt,ptcent,hRatioPromptMin->GetBinContent(iPt+1));
    gSystRatioMin->SetPointError(iPt,0.15,0.15,TMath::Sqrt(minerrnotrack*minerrnotrack+puberrlowreb[iPt]*puberrlowreb[iPt]),TMath::Sqrt(minerrnotrack*minerrnotrack+puberrhighreb[iPt]*puberrhighreb[iPt]));
    gSystRatioInc->SetPoint(iPt,ptcent,hRatioPromptInc->GetBinContent(iPt+1));
    gSystRatioInc->SetPointError(iPt,0.15,0.15,TMath::Sqrt(incerrnotrack*incerrnotrack+puberrlowreb[iPt]*puberrlowreb[iPt]),TMath::Sqrt(incerrnotrack*incerrnotrack+puberrhighreb[iPt]*puberrhighreb[iPt]));
  }
  
  //______________________________________________________________________________________
  //plot cross sections
//  gStyle->SetPadTickX(1);
//  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(1);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.04,"xyz");
  gStyle->SetTitleOffset(1.6,"y");
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);

  TCanvas *cCrossSectionPrompt = new TCanvas("cCrossSectionPrompt","",10,10,900,1000);
  cCrossSectionPrompt->Clear();
  TCanvas *cCrossSectionFD = new TCanvas("cCrossSectionFD","",10,10,800,800);
  cCrossSectionFD->Clear();
  TCanvas* cRatioMethodsPrompt = new TCanvas("cRatioMethodsPrompt","",10,10,800,600); 
  cRatioMethodsPrompt->Clear();
  TCanvas* cRatioMethodsFD = new TCanvas("cRatioMethodsFD","",10,10,800,600); 
  cRatioMethodsFD->Clear();
  TCanvas* cRatioPublishedPrompt = new TCanvas("cRatioPublishedPrompt","",10,10,800,600); 
  cRatioPublishedPrompt->Clear();
  
  TLine *line = new TLine(PtLims[0],1.,PtLims[nPtBins],1.);
  line->SetLineColor(colors[0]);
  line->SetLineStyle(7);
  line->SetLineWidth(2);
  
  TLegend *lFD = new TLegend(0.40,0.68,0.75,0.87);
  lFD->SetTextSizePixels(1000);
  lFD->SetTextSize(0.04);
  lFD->AddEntry(gFONLLFD,"FONLL X A","f");
  lFD->AddEntry(hCrossFDMin,"Analytic Minimisation","lpe");
  lFD->AddEntry(hCrossFDInc,"Incentre Minimisation","lpe");
    
  TLegend *lPrompt = new TLegend(0.40,0.62,0.75,0.87);
  lPrompt->SetTextSizePixels(1000);
  lPrompt->SetTextSize(0.045);
  lPrompt->AddEntry(gFONLLPrompt,"FONLL X A","f");
  lPrompt->AddEntry(hCrossPromptMin,"Analytic Minimisation","lpe");
  lPrompt->AddEntry(hCrossPromptInc,"Incentre Minimisation","lpe");  
  lPrompt->AddEntry(hPubPrompt,"Published Cross Section","lpe");  
  
  TLegend *lRatios = new TLegend(0.2,0.68,0.55,0.87);
  lRatios->SetTextSizePixels(1000);
  lRatios->SetTextSize(0.045);
  lRatios->AddEntry(hRatioPromptMin,"#frac{Minimisation}{Published}","lpe");
  lRatios->AddEntry(hRatioPromptInc,"#frac{Incentre}{Published}","lpe");

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.04);
  latex->SetTextFont(132);
  
  gFONLLPrompt->SetTitle("");
  gFONLLPrompt->SetFillStyle(3154);
  gFONLLPrompt->SetLineWidth(2);
  gFONLLPrompt->SetFillColor(fillColors[6]);
  gFONLLPrompt->SetLineColor(fillColors[6]);
  gFONLLPrompt->GetYaxis()->SetTitle(hCrossPromptMin->GetYaxis()->GetTitle());
  gFONLLPrompt->GetXaxis()->SetTitle(PtAxis->GetTitle());
  gFONLLPrompt->GetXaxis()->SetLimits(hPubPrompt->GetBinLowEdge(1)-hPubPrompt->GetBinWidth(1),hPubPrompt->GetBinLowEdge(hPubPrompt->GetNbinsX())+hPubPrompt->GetBinWidth(hPubPrompt->GetNbinsX())+hPubPrompt->GetBinWidth(1));
  gFONLLPrompt->GetYaxis()->SetRangeUser(8e-1,3e+4);

  gFONLLFD->SetTitle("");
  gFONLLFD->SetFillStyle(3154);
  gFONLLFD->SetFillColor(fillColors[6]);
  gFONLLFD->SetLineColor(fillColors[6]);
  gFONLLFD->SetLineWidth(2);
  gFONLLFD->GetXaxis()->SetLimits(hPubPrompt->GetBinLowEdge(1)-hPubPrompt->GetBinWidth(1),hPubPrompt->GetBinLowEdge(hPubPrompt->GetNbinsX())+hPubPrompt->GetBinWidth(hPubPrompt->GetNbinsX())+hPubPrompt->GetBinWidth(1));
  gFONLLFD->GetYaxis()->SetTitle(hCrossFDMin->GetYaxis()->GetTitle());
  gFONLLFD->GetXaxis()->SetTitle(PtAxis->GetTitle());
  gFONLLFD->GetYaxis()->SetRangeUser(1e-1,1e+4);

  hPubPrompt->SetLineColor(colors[0]);
  hPubPrompt->SetLineWidth(2);
  hPubPrompt->SetMarkerColor(colors[0]);
  hPubPrompt->SetMarkerSize(1);
  hPubPrompt->SetMarkerStyle(markers[0]);
  gPubSysPrompt->SetFillStyle(0);
  gPubSysPrompt->SetLineWidth(2);
  gPubSysPrompt->SetLineColor(colors[0]);
  gPubSysPrompt->SetLineWidth(2);

  hCrossPromptMin->SetStats(kFALSE);
  hCrossPromptMin->SetLineColor(colors[2]);
  hCrossPromptMin->SetLineWidth(2);
  hCrossPromptMin->SetMarkerStyle(markers[1]);
  hCrossPromptMin->SetMarkerSize(1);
  hCrossPromptMin->SetMarkerColor(colors[2]);
  hCrossFDMin->SetStats(kFALSE);
  hCrossFDMin->SetLineWidth(2);
  hCrossFDMin->SetLineColor(colors[2]);
  hCrossFDMin->SetMarkerStyle(markers[1]);
  hCrossFDMin->SetMarkerSize(1);
  hCrossFDMin->SetMarkerColor(colors[2]);
  
  hCrossPromptInc->SetStats(kFALSE);
  hCrossPromptInc->SetLineColor(colors[1]);
  hCrossPromptInc->SetLineWidth(2);
  hCrossPromptInc->SetMarkerStyle(markers[7]);
  hCrossPromptInc->SetMarkerSize(1);
  hCrossPromptInc->SetMarkerColor(colors[1]);
  hCrossFDInc->SetStats(kFALSE);
  hCrossFDInc->SetLineWidth(2);
  hCrossFDInc->SetLineColor(colors[1]);
  hCrossFDInc->SetMarkerStyle(markers[7]);
  hCrossFDInc->SetMarkerSize(1);
  hCrossFDInc->SetMarkerColor(colors[1]);
  hRatioPromptMin->SetLineWidth(2);
  hRatioPromptInc->SetLineWidth(2);
  
  hRatioMethodsPrompt->GetXaxis()->SetTitle(PtAxis->GetTitle());
  hRatioMethodsPrompt->GetYaxis()->SetTitle("#frac{Minimisation}{Incentre}");
  hRatioMethodsPrompt->SetStats(0);
  hRatioMethodsPrompt->SetMarkerStyle(markers[0]);
  hRatioMethodsPrompt->SetMarkerColor(colors[2]);
  hRatioMethodsPrompt->SetLineColor(colors[2]);
  hRatioMethodsPrompt->SetMarkerSize(1);
  
  hRatioMethodsFD->GetXaxis()->SetTitle(PtAxis->GetTitle());
  hRatioMethodsFD->GetYaxis()->SetTitle("#frac{Minimisation}{Incentre}");
  hRatioMethodsFD->SetStats(0);
  hRatioMethodsFD->SetMarkerStyle(markers[0]);
  hRatioMethodsFD->SetMarkerColor(colors[2]);
  hRatioMethodsFD->SetLineColor(colors[2]);
  hRatioMethodsFD->SetMarkerSize(1); 

  hRatioPromptMin->SetStats(0);
  hRatioPromptMin->SetMarkerStyle(markers[1]);
  hRatioPromptMin->SetMarkerSize(1);
  hRatioPromptMin->SetMarkerColor(colors[2]);
  hRatioPromptMin->SetLineColor(colors[2]);
  hRatioPromptMin->GetXaxis()->SetTitle(PtAxis->GetTitle());
  hRatioPromptMin->GetYaxis()->SetTitle("#frac{Cut Variation}{Published}");
  
  hRatioPromptInc->SetStats(0);
  hRatioPromptInc->SetMarkerStyle(markers[7]);
  hRatioPromptInc->SetMarkerSize(1);
  hRatioPromptInc->SetMarkerColor(colors[1]);
  hRatioPromptInc->SetLineColor(colors[1]);
  
  gCrossSystPromptMin->SetLineWidth(2);
  gCrossSystPromptInc->SetLineWidth(2);
  gCrossSystFDMin->SetLineWidth(2);
  gCrossSystFDInc->SetLineWidth(2);

  gSystRatioMin->SetLineColor(colors[2]);
  gSystRatioMin->SetLineWidth(2);
  gSystRatioInc->SetLineColor(colors[1]);
  gSystRatioInc->SetLineWidth(2);
  
  cCrossSectionPrompt->Clear();
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->SetLogy();
  pad1->cd();
  gFONLLPrompt->GetYaxis()->SetRangeUser(0.2,40000);
  gFONLLPrompt->GetXaxis()->SetTitleSize(0.05);
  gFONLLPrompt->GetYaxis()->SetTitleSize(0.05);
  gFONLLPrompt->Draw("a5");
  hPubPrompt->Draw("Esame");
  gPubSysPrompt->Draw("2");
  hCrossPromptMin->Draw("Esame");
  hCrossPromptInc->Draw("Esame");
  gCrossSystPromptMin->Draw("2");
  gCrossSystPromptInc->Draw("2");
  lPrompt->Draw("same");
  latex->DrawLatex(PtLims[0],hCrossPromptMin->GetMinimum()*0.4,"Prompt D^{+}");
  latex->DrawLatex(PtLims[0],hCrossPromptMin->GetMinimum()*0.18,"pPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->DrawLatex(PtLims[0],hCrossPromptMin->GetMinimum()*0.07,"#pm 2.1% BR unc. not shown");
  latex->DrawLatex(PtLims[0],hCrossPromptMin->GetMinimum()*0.03,"#pm 3.7% norm. unc. not shown");
  cCrossSectionPrompt->cd();
  TPad *pad2 = new TPad("pad2","pad2",0.,0.,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();
  Double_t PtLims2[12] = {0,1,2,3,4,5,6,8,12,16,24,25};
  TLine* line2 = (TLine*)line->Clone();
  line2->SetX1(0);
  line2->SetX2(25);
  TH1F* hRatioPromptMinCopy= new TH1F("hRatioPromptMinCopy","",11,PtLims2);
  hRatioPromptMinCopy->SetBinContent(1,0);
  hRatioPromptMinCopy->SetBinContent(2,0);
  hRatioPromptMinCopy->SetBinContent(10,0);
  hRatioPromptMinCopy->SetBinContent(11,0);
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    hRatioPromptMinCopy->SetBinContent(iPt+3,hRatioPromptMin->GetBinContent(iPt+1));
    hRatioPromptMinCopy->SetBinError(iPt+3,hRatioPromptMin->GetBinError(iPt+1));
  }
  hRatioPromptMinCopy->SetStats(0);
  hRatioPromptMinCopy->SetMarkerStyle(hRatioPromptMin->GetMarkerStyle());
  hRatioPromptMinCopy->SetMarkerSize(hRatioPromptMin->GetMarkerSize());
  hRatioPromptMinCopy->SetMarkerColor(hRatioPromptMin->GetMarkerColor());
  hRatioPromptMinCopy->SetLineColor(hRatioPromptMin->GetLineColor());
  hRatioPromptMinCopy->SetLineWidth(hRatioPromptMin->GetLineWidth());
  hRatioPromptMinCopy->GetYaxis()->SetTitle("#frac{Cut-variation}{Published}");
  hRatioPromptMinCopy->GetYaxis()->SetTitleOffset(0.8);
  hRatioPromptMinCopy->GetYaxis()->SetTitleSize(0.11);
  hRatioPromptMinCopy->GetYaxis()->SetLabelSize(0.1);
  hRatioPromptMinCopy->GetXaxis()->SetLabelSize(0.1);
  hRatioPromptMinCopy->GetYaxis()->SetRangeUser(0.5,1.19);
  hRatioPromptMinCopy->GetXaxis()->SetTitle(hRatioPromptMin->GetXaxis()->GetTitle());
  hRatioPromptMinCopy->GetXaxis()->SetTitleSize(0.12);
  hRatioPromptMinCopy->Draw("E");
  hRatioPromptInc->Draw("Esame");
  gSystRatioMin->Draw("2");
  gSystRatioInc->Draw("2");
  line2->Draw("same");
  
  cCrossSectionFD->Clear();
  cCrossSectionFD->SetLogy();
  gFONLLFD->GetXaxis()->SetTitleSize(0.05);
  gFONLLFD->GetYaxis()->SetTitleSize(0.05);
  gFONLLFD->GetYaxis()->SetRangeUser(0.03,9000);
  gFONLLFD->Draw("a5");
  hCrossFDMin->Draw("Esame");
  hCrossFDInc->Draw("Esame");
  gCrossSystFDMin->Draw("2");
  gCrossSystFDInc->Draw("2");
  lFD->Draw("same");
  latex->DrawLatex(PtLims[0],hCrossFDMin->GetMinimum()*0.28,"Feed-down D^{+}");
  latex->DrawLatex(PtLims[0],hCrossFDMin->GetMinimum()*0.13,"pPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->DrawLatex(PtLims[0],hCrossFDMin->GetMinimum()*0.055,"#pm 2.1% BR unc. not shown");
  latex->DrawLatex(PtLims[0],hCrossFDMin->GetMinimum()*0.025,"#pm 3.7% norm. unc. not shown");

  cRatioMethodsPrompt->Clear();
  hRatioMethodsPrompt->Draw("E");
  line->Draw("same");

  cRatioMethodsFD->Clear();
  hRatioMethodsFD->Draw("E");
  line->Draw("same");

  cRatioPublishedPrompt->Clear();
  hRatioPromptMin->GetYaxis()->SetRangeUser(0.5,1.5);
  hRatioPromptMin->Draw("E");
  hRatioPromptInc->Draw("Esame");
  gSystRatioMin->Draw("2");
  gSystRatioInc->Draw("2");
  line->Draw("same");
  lRatios->Draw("same");
      
  cCrossSectionFD->SaveAs(Form("CrossSectionFD%s.eps",system.Data()));
  cCrossSectionFD->SaveAs(Form("CrossSectionFD%s.root",system.Data()));
  cCrossSectionPrompt->SaveAs(Form("CrossSectionPrompt%s.eps",system.Data()));
  cCrossSectionPrompt->SaveAs(Form("CrossSectionPrompt%s.root",system.Data()));
  cRatioMethodsPrompt->SaveAs(Form("RatioMethodsPrompt%s.eps",system.Data()));
  cRatioMethodsPrompt->SaveAs(Form("RatioMethodsPrompt%s.root",system.Data()));
  cRatioMethodsFD->SaveAs(Form("RatioMethodsFD%s.eps",system.Data()));
  cRatioMethodsFD->SaveAs(Form("RatioMethodsFD%s.root",system.Data()));
  cRatioPublishedPrompt->SaveAs(Form("RatioPublishedPrompt%s.eps",system.Data()));
  cRatioPublishedPrompt->SaveAs(Form("RatioPublishedPrompt%s.root",system.Data()));

  TFile outfileMin("DplusCrossSection_CutVar_Min.root","RECREATE");
  hCrossPromptMin->Write();
  gCrossSystPromptMin->Write();
  hCrossFDMin->Write();
  gCrossSystFDMin->Write();
  outfileMin.Close();

  TFile systMin("DplusSystErr_CutVar_Min.root","RECREATE");
  gRawSystPromptMin->Write();
  gCutSystPromptMin->Write();
  gEffSystPromptMin->Write();
  gTrackSystPromptMin->Write();
  gSystPromptMin->Write();
  gRawSystFDMin->Write();
  gCutSystFDMin->Write();
  gEffSystFDMin->Write();
  gTrackSystFDMin->Write();
  gSystFDMin->Write();
  systMin.Close();
  
  TFile outfileInc("DplusCrossSection_CutVar_Inc.root","RECREATE");
  hCrossPromptInc->Write();
  gCrossSystPromptInc->Write();
  hCrossFDInc->Write();
  gCrossSystFDInc->Write();
  outfileInc.Close();

  TFile systInc("DplusSystErr_CutVar_Inc.root","RECREATE");
  gRawSystPromptInc->Write();
  gCutSystPromptInc->Write();
  gEffSystPromptInc->Write();
  gTrackSystPromptInc->Write();
  gSystPromptInc->Write();
  gRawSystFDInc->Write();
  gCutSystFDInc->Write();
  gEffSystFDInc->Write();
  gTrackSystFDInc->Write();
  gSystFDInc->Write();
  systInc.Close();
}
