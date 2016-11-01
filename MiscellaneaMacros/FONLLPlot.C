#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TMath.h>
#include <TH1.h>  
#include <TCanvas.h>
#include <TString.h>
#endif

void FONLLPlot(TString infilename = "D0DplusDstarPredictions_5023TeV_y05_all_050515_BDShapeCorrected.root", TString outfilename = "FONLL5TeVpPb") {
  
  //____________________________________________________________________________________________________
  
  const Int_t nPtBins = 10;
  const Int_t nPtLims = nPtBins+1;
  Double_t PtLims[nPtLims] = {1,2,3,4,5,6,7,8,12,16,24};
  
  //____________________________________________________________________________________________________

  TFile infile(infilename.Data(),"UPDATE");
  TH1F* hFONLLPrompt = (TH1F*)infile.Get("hDpluskpipipred_central");
  TH1F* hFONLLPromptLow = (TH1F*)infile.Get("hDpluskpipipred_min");
  TH1F* hFONLLPromptHigh = (TH1F*)infile.Get("hDpluskpipipred_max");
  TH1F* hFONLLFD = (TH1F*)infile.Get("hDpluskpipifromBpred_central_corr");
  TH1F* hFONLLFDLow = (TH1F*)infile.Get("hDpluskpipifromBpred_min_corr");
  TH1F* hFONLLFDHigh = (TH1F*)infile.Get("hDpluskpipifromBpred_max_corr");

  TH1F* hFONLLPromptReb = (TH1F*)hFONLLPrompt->Rebin(nPtBins,"hFONLLPromptReb",PtLims);
  TH1F* hFONLLPromptLowReb = (TH1F*)hFONLLPromptLow->Rebin(nPtBins,"hFONLLPromptLowReb",PtLims);
  TH1F* hFONLLPromptHighReb = (TH1F*)hFONLLPromptHigh->Rebin(nPtBins,"hFONLLPromptHighReb",PtLims);
  TH1F* hFONLLFDReb = (TH1F*)hFONLLFD->Rebin(nPtBins,"hFONLLFDReb",PtLims);
  TH1F* hFONLLFDLowReb = (TH1F*)hFONLLFDLow->Rebin(nPtBins,"hFONLLFDLowReb",PtLims);
  TH1F* hFONLLFDHighReb = (TH1F*)hFONLLFDHigh->Rebin(nPtBins,"hFONLLFDHighReb",PtLims);

  TGraphAsymmErrors* gFONLLPrompt = new TGraphAsymmErrors(nPtBins);
  TGraphAsymmErrors* gFONLLFD = new TGraphAsymmErrors(nPtBins);

  Double_t BR = 1.;//0.0913;
  Double_t A = 208./2;

  TH1F* hFONLLPromptCentral = new TH1F("hFONLLPromptCentral","",nPtBins,PtLims);
  TH1F* hFONLLFDCentral = new TH1F("hFONLLFDCentral","",nPtBins,PtLims);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    gFONLLFD->SetPoint(iPt,(PtLims[iPt]+PtLims[iPt+1])/2,hFONLLFDReb->GetBinContent(iPt+1)/hFONLLFDReb->GetBinWidth(iPt+1)/1000000/BR*A);
    gFONLLFD->SetPointError(iPt,hFONLLFDReb->GetBinWidth(iPt+1)/2,hFONLLFDReb->GetBinWidth(iPt+1)/2,(hFONLLFDReb->GetBinContent(iPt+1)-hFONLLFDLowReb->GetBinContent(iPt+1))/hFONLLFDReb->GetBinWidth(iPt+1)/1000000/BR*A,(hFONLLFDHighReb->GetBinContent(iPt+1)-hFONLLFDReb->GetBinContent(iPt+1))/hFONLLFDReb->GetBinWidth(iPt+1)/1000000/BR*A);
    
    gFONLLPrompt->SetPoint(iPt,(PtLims[iPt]+PtLims[iPt+1])/2,hFONLLPromptReb->GetBinContent(iPt+1)/hFONLLPromptReb->GetBinWidth(iPt+1)/1000000/BR*A);
    gFONLLPrompt->SetPointError(iPt,hFONLLPromptReb->GetBinWidth(iPt+1)/2,hFONLLPromptReb->GetBinWidth(iPt+1)/2,(hFONLLPromptReb->GetBinContent(iPt+1)-hFONLLPromptLowReb->GetBinContent(iPt+1))/hFONLLPromptReb->GetBinWidth(iPt+1)/1000000/BR*A,(hFONLLPromptHighReb->GetBinContent(iPt+1)-hFONLLPromptReb->GetBinContent(iPt+1))/hFONLLPromptReb->GetBinWidth(iPt+1)/1000000/BR*A);

    hFONLLPromptCentral->SetBinContent(iPt+1,hFONLLPromptReb->GetBinContent(iPt+1)/hFONLLPromptReb->GetBinWidth(iPt+1)/1000000/BR*A);
    hFONLLFDCentral->SetBinContent(iPt+1,hFONLLFDReb->GetBinContent(iPt+1)/hFONLLFDReb->GetBinWidth(iPt+1)/1000000/BR*A);
    hFONLLPromptCentral->SetBinError(iPt+1,1.e-10);
    hFONLLFDCentral->SetBinError(iPt+1,1.e-10);
  }
  
  gFONLLPrompt->SetName("gFONLLPrompt");
  gFONLLFD->SetName("gFONLLFD");
    
  TCanvas *cPrompt= new TCanvas("cPrompt","cPrompt",10,10,1200,800);
  cPrompt->Clear();
  cPrompt->SetLogy();
  gFONLLPrompt->SetFillColor(4);
  gFONLLPrompt->SetFillStyle(3001);
  gFONLLPrompt->Draw("A2");
  
  TFile outfilePrompt(Form("%sPrompt.root",outfilename.Data()),"RECREATE");  
  gFONLLPrompt->Write();
  hFONLLPromptCentral->Write();
  outfilePrompt.Close();

  TCanvas *cFD= new TCanvas("cFD","cFD",10,10,1200,800);
  cFD->Clear();
  cFD->SetLogy();
  gFONLLFD->SetFillColor(4);
  gFONLLFD->SetFillStyle(3001);
  gFONLLFD->Draw("A2");
  
  TFile outfileFD(Form("%sFD.root",outfilename.Data()),"RECREATE");  
  gFONLLFD->Write();
  hFONLLFDCentral->Write();
  outfileFD.Close();
}
