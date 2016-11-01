#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TInterpreter.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <THnSparse.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <TLine.h>
#include <TList.h>
#include <TArrayD.h>

#include "AliNormalizationCounter.h"
#endif

void PlotCrossSection();
void CrossSectionRawYieldsVariation();
TH1F* CalculateCrossSection(TH1F* hPtRawYields,TH1F* hPtEffXAcc,TH1F* hPtPromptFrac,Double_t Nev,Double_t BR=0.0913,Double_t sigma=2.09);
//sigma in barn

void PlotCrossSection() {
  
  //__________________________________________________________________________________
  //input files
  TString infilename = "/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root";
  TString dirdataname = "PWG3_D2H_InvMassDplus";
  TString listdataname = "coutputDplus_ImpParpPbData0100";

  TFile infiledata(infilename.Data(),"UPDATE");
  TDirectoryFile* dirdata=(TDirectoryFile*)infiledata.Get(dirdataname.Data()); 
  TList* listdata=(TList*)dirdata->Get(listdataname.Data());

  TString normobj=listdataname;
  normobj.ReplaceAll("coutputDplus","coutputDplusNorm");
  AliNormalizationCounter* nc=(AliNormalizationCounter*)dirdata->Get(normobj.Data());
  Double_t nEvents = nc->GetNEventsForNorm();
  
  TFile efffile("efficiency_and_acceptance.root","UPDATE");
  TH1F* hPtEffXAcc = (TH1F*)efffile.Get("hPtEffXAccPrompt");
  hPtEffXAcc->SetDirectory(0);
  efffile.Close();

  TFile fracfile("../fprompt.root","UPDATE");
  TH1F* hPtPromptFrac = (TH1F*)fracfile.Get("hPromptFraction");
  hPtPromptFrac->SetDirectory(0);
  TGraphAsymmErrors* gPtPromptFrac = (TGraphAsymmErrors*)fracfile.Get("gPromptFractionSyst");
  fracfile.Close();

  TFile fracuncfile("../fPrompt_syst_unc.root","UPDATE");
  TGraphAsymmErrors* gPtPromptFracSyst = (TGraphAsymmErrors*)fracuncfile.Get("gSyst");
  fracuncfile.Close();

  TFile rawYfile("rawyields.root","UPDATE");
  TH1F* hPtRawYields = (TH1F*)rawYfile.Get("hPtRawYields");
  hPtRawYields->SetDirectory(0);
  rawYfile.Close();

  TFile pubfile("DplusCrossSec_method2_fd2_br1.root","UPDATE");
  TH1F* hPubCrossSection = (TH1F*)pubfile.Get("hAAC");
  hPubCrossSection->SetDirectory(0);
  TGraphAsymmErrors* gPubCrossSection=(TGraphAsymmErrors*)pubfile.Get("gaaCsystTot");
  pubfile.Close();

  TFile FONLLpromptfile("FONLL5TeVpPbPrompt.root","UPDATE");
  TGraphAsymmErrors* gFONLLprompt=(TGraphAsymmErrors*)FONLLpromptfile.Get("gFONLLPrompt");
  TH1F* hFONLLprompt=(TH1F*)FONLLpromptfile.Get("hFONLLPromptCentral");
  hFONLLprompt->SetDirectory(0);
  FONLLpromptfile.Close();
  
  TFile FONLLFDfile("FONLL5TeVpPbFD.root","UPDATE");
  TGraphAsymmErrors* gFONLLFD=(TGraphAsymmErrors*)FONLLFDfile.Get("gFONLLFD");
  TH1F* hFONLLFD=(TH1F*)FONLLFDfile.Get("hFONLLFDCentral");
  hFONLLFD->SetDirectory(0);
  FONLLFDfile.Close();

  //____________________________________________________________________________________
  //prompt cross section calculation
  TH1F* hPtPromptCrossSection = (TH1F*)CalculateCrossSection(hPtRawYields,hPtEffXAcc,hPtPromptFrac,nEvents);

  //feed-down cross section calculation
  const Int_t nPtBins = hPtPromptCrossSection->GetNbinsX();
  const Int_t nPtLims = nPtBins+1;
  TAxis* PtAxis = (TAxis*)hPtPromptCrossSection->GetXaxis();
  TArrayD* ptarray = (TArrayD*)PtAxis->GetXbins();
  Double_t* PtLims = (Double_t*)ptarray->GetArray();

  TH1F* hPtFDFrac = new TH1F("hPtFDFrac","",nPtBins,PtLims); 
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    hPtFDFrac->SetBinContent(iPt+1,1-hPtPromptFrac->GetBinContent(iPt+1));
    hPtFDFrac->SetBinError(iPt+1,hPtPromptFrac->GetBinError(iPt+1));    
  }
  TH1F* hPtFDCrossSection = (TH1F*)CalculateCrossSection(hPtRawYields,hPtEffXAcc,hPtFDFrac,nEvents);

  //______________________________________________________________________________________
  //systematic uncertainty on the prompt cross section 
  Double_t trackerr[nPtBins]={0.09,0.09,0.09,0.09,0.09,0.09,0.09};
  Double_t efferr[nPtBins]={0.00,0.0,0.0,0.0,0.0,0.04,0.04};
  Double_t cutvarerr[nPtBins]={0.08,0.05,0.05,0.05,0.05,0.05,0.08};
  Double_t fracerr[nPtBins]={0.};
  Double_t systerr[nPtBins]={0.};
  TGraphAsymmErrors* gPtPromptCrossSection = new TGraphAsymmErrors(nPtBins);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    fracerr[iPt] = gPtPromptFracSyst->GetErrorYhigh(iPt);
    Double_t ptcent = (PtLims[iPt]+PtLims[iPt+1])/2;
    Double_t dsdpt = hPtPromptCrossSection->GetBinContent(iPt+1);
    Double_t syst = TMath::Sqrt(efferr[iPt]*efferr[iPt]+trackerr[iPt]*trackerr[iPt]+cutvarerr[iPt]*cutvarerr[iPt]+fracerr[iPt]*fracerr[iPt])*dsdpt;
    gPtPromptCrossSection->SetPoint(iPt,ptcent,dsdpt);
    gPtPromptCrossSection->SetPointError(iPt,0.15,0.15,syst,syst);
  }
  
  //______________________________________________________________________________________
  //calculate ratios w.r.t. the published cross section (only prompt) 
  const Int_t nPtBinsPub = hPubCrossSection->GetNbinsX();
  const Int_t nPtLimsPub = nPtBins+1;
  TAxis* PtAxisPub = (TAxis*)hPubCrossSection->GetXaxis();
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
  
  TH1F* hRatioCrossSection = new TH1F("hRatioCrossSection","",nPtBins,PtLims);
  
  if(IsBinningEqual) {
    hRatioCrossSection->Divide(hPtPromptCrossSection,hPubCrossSection,1.,1.);
  }
  else {
    TH1F* hPubCrossSection2 = (TH1F*)hPubCrossSection->Clone();
    for(Int_t iPt=0; iPt<nPtBinsPub; iPt++) {
      hPubCrossSection2->SetBinContent(iPt+1,hPubCrossSection2->GetBinContent(iPt+1)*(PtLimsPub[iPt+1]-PtLimsPub[iPt]));
      hPubCrossSection2->SetBinError(iPt+1,hPubCrossSection2->GetBinError(iPt+1)*(PtLimsPub[iPt+1]-PtLimsPub[iPt]));
    }
    TH1F* hPubCrossSectionReb = (TH1F*)hPubCrossSection2->Rebin(nPtBins,"hPubCrossSectionReb",PtLims);
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      hPubCrossSectionReb->SetBinContent(iPt+1,hPubCrossSectionReb->GetBinContent(iPt+1)/(PtLims[iPt+1]-PtLims[iPt]));
      hPubCrossSectionReb->SetBinError(iPt+1,hPubCrossSectionReb->GetBinError(iPt+1)/(PtLims[iPt+1]-PtLims[iPt]));
    }
    hRatioCrossSection->Divide(hPtPromptCrossSection,hPubCrossSectionReb,1.,1.);
  }
  
  //____________________________________________________________________________________
  //plot
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTickLength(0.02,"X");
  gStyle->SetTickLength(0.02,"Y"); 
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);
  
  TLegend *l = new TLegend(0.5,0.65,0.85,0.85);
  l->AddEntry(hPtPromptCrossSection,"Impact parameter fit","lpe");
  l->AddEntry(hPubCrossSection,"Published cross section","lpe");
  l->AddEntry(gFONLLprompt,"FONLL X A","f");
  l->SetTextSize(0.03);

  TLegend *l2 = new TLegend(0.5,0.65,0.85,0.85);
  l2->AddEntry(hPtFDCrossSection,"Impact parameter fit","lpe");
  l2->AddEntry(gFONLLFD,"FONLL X A","f");
  l2->SetTextSize(0.03);

  gFONLLprompt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  gFONLLprompt->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}} (#mub c/GeV)");
  gFONLLprompt->GetYaxis()->SetTitleOffset(1.5);
  gFONLLprompt->GetXaxis()->SetTitleOffset(1.2);

  gFONLLFD->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  gFONLLFD->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}} (#mub c/GeV)");
  gFONLLFD->GetYaxis()->SetTitleOffset(1.5);
  gFONLLFD->GetXaxis()->SetTitleOffset(1.2);
  
  hPubCrossSection->SetLineColor(kBlack);
  hPubCrossSection->SetMarkerColor(kBlack);
  hPubCrossSection->SetMarkerStyle(20);
  gPubCrossSection->SetLineColor(kBlack);
  hPubCrossSection->SetLineWidth(1.5);
  gPubCrossSection->SetLineWidth(1.5);

  hPtPromptCrossSection->SetLineWidth(1.5);
  hPtPromptCrossSection->SetMarkerStyle(20);
  hPtPromptCrossSection->SetMarkerColor(kBlue);
  hPtPromptCrossSection->SetLineColor(kBlue);
  gPtPromptCrossSection->SetLineWidth(1.5);
  gPtPromptCrossSection->SetFillStyle(20);
  gPtPromptCrossSection->SetLineColor(kBlue);

  hPtFDCrossSection->SetLineWidth(1.5);
  hPtFDCrossSection->SetMarkerStyle(20);
  hPtFDCrossSection->SetMarkerColor(kBlue);
  hPtFDCrossSection->SetLineColor(kBlue);

  gFONLLprompt->SetFillStyle(20);
  gFONLLprompt->SetLineColor(kRed);
  gFONLLprompt->GetYaxis()->SetRangeUser(hPubCrossSection->GetMinimum()*0.2,hPubCrossSection->GetMaximum()*2.);
  
  gFONLLFD->SetFillStyle(20);
  gFONLLFD->SetLineColor(kRed);
  gFONLLFD->GetYaxis()->SetRangeUser(hPtFDCrossSection->GetMinimum()*0.2,hPtFDCrossSection->GetMaximum()*2.);
   
  hFONLLprompt->SetLineColor(kRed);
  hFONLLprompt->SetMarkerColor(kRed);
  hFONLLFD->SetLineColor(kRed);
  hFONLLFD->SetMarkerColor(kRed);

  TCanvas *cPromptCross = new TCanvas("cPromptCross","",10,10,800,800);
  cPromptCross->Clear();
  cPromptCross->SetLogy();
  gFONLLprompt->GetXaxis()->SetRangeUser(PtLimsPub[0],PtLimsPub[nPtBinsPub]);
  gFONLLprompt->Draw("a2");
  hFONLLprompt->Draw("Esame");
  hPubCrossSection->Draw("Esame");
  gPubCrossSection->Draw("2same");
  hPtPromptCrossSection->Draw("Esame");
  gPtPromptCrossSection->Draw("2");
  l->Draw("same");

  TCanvas *cFDCross = new TCanvas("cFDCross","",10,10,800,800);
  cFDCross->Clear();
  cFDCross->SetLogy();
  gFONLLFD->GetXaxis()->SetRangeUser(PtLimsPub[0],PtLimsPub[nPtBinsPub]);
  gFONLLFD->Draw("a2");
  hFONLLFD->Draw("Esame");
  hPtFDCrossSection->Draw("Esame");
  l2->Draw("same");

  TLine* line = new TLine(PtLims[0],1.,PtLims[nPtBins],1.);
  line->SetLineColor(kBlack);
  line->SetLineStyle(6);

  hRatioCrossSection->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hRatioCrossSection->GetYaxis()->SetTitle("#frac{Impact parameter}{Published}");
  hRatioCrossSection->GetYaxis()->SetTitleOffset(1.5);
  hRatioCrossSection->GetXaxis()->SetTitleOffset(1.2);
  hRatioCrossSection->SetLineColor(kRed);
  hRatioCrossSection->SetMarkerStyle(20);
  hRatioCrossSection->SetMarkerColor(kRed);

  TCanvas *cRatio = new TCanvas("cRatio","",10,10,800,800);
  cRatio->Clear();
  hRatioCrossSection->Draw("E");
  line->Draw("same");
  
  cPromptCross->SaveAs("PromptCrossSection_ImpParFit.eps");
  cPromptCross->SaveAs("PromptCrossSection_ImpParFit.root");
  cFDCross->SaveAs("FDCrossSection_ImpParFit.eps");
  cFDCross->SaveAs("FDCrossSection_ImpParFit.root");
  cRatio->SaveAs("RatioPrompt_ImpParFit.eps");
  cRatio->SaveAs("RatioPrompt_ImpParFit.root");

  TFile outfile("DplusCrossSection_ImpParFit.root","RECREATE");
  hPtPromptCrossSection->SetName("hPromptCrossSec");
  gPtPromptCrossSection->SetName("gPromptCrossSec");
  hPtPromptCrossSection->Write();
  gPtPromptCrossSection->Write();
  hPtFDCrossSection->SetName("hFDCrossSec");
  hPtFDCrossSection->Write();
  outfile.Close();
  
}

void CrossSectionRawYieldsVariation() {

  //____________________________________________________________________________________
  //input files
  TString infilename = "/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root";
  TString dirdataname = "PWG3_D2H_InvMassDplus";
  TString listdataname = "coutputDplus_ImpParpPbData0100";

  TFile infiledata(infilename.Data(),"UPDATE");
  TDirectoryFile* dirdata=(TDirectoryFile*)infiledata.Get(dirdataname.Data()); 
  TList* listdata=(TList*)dirdata->Get(listdataname.Data());

  TString normobj=listdataname;
  normobj.ReplaceAll("coutputDplus","coutputDplusNorm");
  AliNormalizationCounter* nc=(AliNormalizationCounter*)dirdata->Get(normobj.Data());
  Double_t nEvents = nc->GetNEventsForNorm();
    
  TFile efffile("efficiency_and_acceptance.root","UPDATE");
  TH1F* hPtEffXAcc = (TH1F*)efffile.Get("hPtEffXAccPrompt");
  hPtEffXAcc->SetDirectory(0);
  efffile.Close();

  TFile fracfile("../fprompt.root","UPDATE");
  TH1F* hPtPromptFrac = (TH1F*)fracfile.Get("hPromptFraction");
  hPtPromptFrac->SetDirectory(0);
  TGraphAsymmErrors* gPtPromptFrac = (TGraphAsymmErrors*)fracfile.Get("gPromptFractionSyst");
  fracfile.Close();

  TFile fracSvarfile("../../systematics/SoverT/promptfraction_syst_SoverT.root","UPDATE");
  TH1F* hPtPromptFracLowS = (TH1F*)fracSvarfile.Get("hFrac3");
  TH1F* hPtPromptFracLowStat = (TH1F*)fracSvarfile.Get("hFrac4");
  TH1F* hPtPromptFracLowSyst = (TH1F*)fracSvarfile.Get("hFrac5");
  TH1F* hPtPromptFracHighS = (TH1F*)fracSvarfile.Get("hFrac0");
  TH1F* hPtPromptFracHighStat = (TH1F*)fracSvarfile.Get("hFrac1");
  TH1F* hPtPromptFracHighSyst = (TH1F*)fracSvarfile.Get("hFrac2");
  hPtPromptFracLowS->SetDirectory(0);
  hPtPromptFracLowStat->SetDirectory(0);
  hPtPromptFracLowSyst->SetDirectory(0);
  hPtPromptFracHighS->SetDirectory(0);
  hPtPromptFracHighStat->SetDirectory(0);
  hPtPromptFracHighSyst->SetDirectory(0);
  fracSvarfile.Close();

  TFile rawYfile("rawyields.root","UPDATE");
  TH1F* hPtRawYields = (TH1F*)rawYfile.Get("hPtRawYields");
  TH1F* hPtRawYieldsLowS = (TH1F*)rawYfile.Get("hPtRawYieldsLowS");
  TH1F* hPtRawYieldsLowStat = (TH1F*)rawYfile.Get("hPtRawYieldsLowStat");
  TH1F* hPtRawYieldsLowSyst = (TH1F*)rawYfile.Get("hPtRawYieldsLowSyst");
  TH1F* hPtRawYieldsHighS = (TH1F*)rawYfile.Get("hPtRawYieldsHighS");
  TH1F* hPtRawYieldsHighStat = (TH1F*)rawYfile.Get("hPtRawYieldsHighStat");
  TH1F* hPtRawYieldsHighSyst = (TH1F*)rawYfile.Get("hPtRawYieldsHighSyst");
  hPtRawYields->SetDirectory(0);
  hPtRawYieldsLowS->SetDirectory(0);
  hPtRawYieldsLowStat->SetDirectory(0);
  hPtRawYieldsLowSyst->SetDirectory(0);
  hPtRawYieldsHighS->SetDirectory(0);
  hPtRawYieldsHighStat->SetDirectory(0);
  hPtRawYieldsHighSyst->SetDirectory(0);
  rawYfile.Close();

  //____________________________________________________________________________________
  //pt bins 
  const Int_t nPtBins = hPtRawYields->GetNbinsX();
  const Int_t nPtLims = nPtBins+1;
  TArrayD* ptarray = (TArrayD*)hPtRawYields->GetXaxis()->GetXbins();
  Double_t* PtLims = (Double_t*)ptarray->GetArray();

  //____________________________________________________________________________________
  //calculate cross sections
  TH1F* hPtCrossSection = (TH1F*)CalculateCrossSection(hPtRawYields,hPtEffXAcc,hPtPromptFrac,nEvents);
  TH1F* hPtCrossSectionLowS = (TH1F*)CalculateCrossSection(hPtRawYieldsLowS,hPtEffXAcc,hPtPromptFracLowS,nEvents);
  TH1F* hPtCrossSectionLowStat = (TH1F*)CalculateCrossSection(hPtRawYieldsLowStat,hPtEffXAcc,hPtPromptFracLowStat,nEvents);
  TH1F* hPtCrossSectionLowSyst = (TH1F*)CalculateCrossSection(hPtRawYieldsLowSyst,hPtEffXAcc,hPtPromptFracLowSyst,nEvents);
  TH1F* hPtCrossSectionHighS = (TH1F*)CalculateCrossSection(hPtRawYieldsHighS,hPtEffXAcc,hPtPromptFracHighS,nEvents);
  TH1F* hPtCrossSectionHighStat = (TH1F*)CalculateCrossSection(hPtRawYieldsHighStat,hPtEffXAcc,hPtPromptFracHighStat,nEvents);
  TH1F* hPtCrossSectionHighSyst = (TH1F*)CalculateCrossSection(hPtRawYieldsHighSyst,hPtEffXAcc,hPtPromptFracHighSyst,nEvents);

  //____________________________________________________________________________________
  //calculate ratios
  TH1F* hPtRatioLowS = new TH1F("hPtRatioLowS","",nPtBins,PtLims);
  hPtRatioLowS->Divide(hPtCrossSectionLowS,hPtCrossSection,1.,1.);
  TH1F* hPtRatioLowStat = new TH1F("hPtRatioLowStat","",nPtBins,PtLims);
  hPtRatioLowStat->Divide(hPtCrossSectionLowStat,hPtCrossSection,1.,1.);
  TH1F* hPtRatioLowSyst = new TH1F("hPtRatioLowSyst","",nPtBins,PtLims);
  hPtRatioLowSyst->Divide(hPtCrossSectionLowSyst,hPtCrossSection,1.,1.);
  TH1F* hPtRatioHighS = new TH1F("hPtRatioHighS","",nPtBins,PtLims);
  hPtRatioHighS->Divide(hPtCrossSectionHighS,hPtCrossSection,1.,1.);
  TH1F* hPtRatioHighStat = new TH1F("hPtRatioHighStat","",nPtBins,PtLims);
  hPtRatioHighStat->Divide(hPtCrossSectionHighStat,hPtCrossSection,1.,1.);
  TH1F* hPtRatioHighSyst = new TH1F("hPtRatioHighSyst","",nPtBins,PtLims);
  hPtRatioHighSyst->Divide(hPtCrossSectionHighSyst,hPtCrossSection,1.,1.);

  //____________________________________________________________________________________
  //plot cross sections and ratios

  const Int_t colors[10] = {kRed,kBlue,kGreen+3,kBlack,kMagenta,kTeal-5,kOrange+7,kYellow-3,kCyan+3,kGreen};

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetHistLineWidth(2);
  gStyle->SetTickLength(0.02,"X");
  gStyle->SetTickLength(0.02,"Y"); 
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.14);

  hPtCrossSection->SetLineColor(colors[0]);
  hPtCrossSectionHighS->SetLineColor(colors[1]);
  hPtCrossSectionHighStat->SetLineColor(colors[2]);
  hPtCrossSectionHighSyst->SetLineColor(colors[3]);
  hPtCrossSectionLowS->SetLineColor(colors[4]);  
  hPtCrossSectionLowStat->SetLineColor(colors[5]);
  hPtCrossSectionLowSyst->SetLineColor(colors[6]);
  hPtRatioHighS->SetLineColor(colors[1]);
  hPtRatioHighStat->SetLineColor(colors[2]);
  hPtRatioHighSyst->SetLineColor(colors[3]);
  hPtRatioLowS->SetLineColor(colors[4]);
  hPtRatioLowStat->SetLineColor(colors[5]);
  hPtRatioLowSyst->SetLineColor(colors[6]);
  hPtCrossSection->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hPtCrossSection->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}} (#mub c/GeV)");
  hPtCrossSection->GetYaxis()->SetTitleOffset(1.5);
  hPtCrossSection->GetXaxis()->SetTitleOffset(1.2);
  hPtCrossSection->GetYaxis()->SetTitleSize(0.04);
  hPtCrossSection->GetXaxis()->SetTitleSize(0.04);
  hPtRatioLowS->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hPtRatioLowS->GetYaxis()->SetTitle("ratio of d#sigma/d#it{p}_{T} w.r.t the mean value");
  hPtRatioLowS->GetYaxis()->SetTitleOffset(1.5);
  hPtRatioLowS->GetXaxis()->SetTitleOffset(1.2);
  hPtRatioLowS->GetYaxis()->SetTitleSize(0.04);
  hPtRatioLowS->GetXaxis()->SetTitleSize(0.04);
  
  TLegend *l = new TLegend(0.58,0.58,0.88,0.88);
  l->SetTextSize(0.045);
  TLegend *l2 = (TLegend*)l->Clone();
  l->AddEntry(hPtCrossSection,"Reference value","l");
  l->AddEntry(hPtCrossSectionHighS,"S + #sigma_{S}","l");
  l->AddEntry(hPtCrossSectionHighStat,"S + #sigma_{S}(stat)","l");
  l->AddEntry(hPtCrossSectionHighSyst,"S + #sigma_{S}(syst)","l");
  l->AddEntry(hPtCrossSectionLowS,"S - #sigma_{S}","l");
  l->AddEntry(hPtCrossSectionLowStat,"S - #sigma_{S}(stat)","l");
  l->AddEntry(hPtCrossSectionLowSyst,"S - #sigma_{S}(syst)","l");
  l2->AddEntry(hPtCrossSectionHighS,"S + #sigma_{S}","l");
  l2->AddEntry(hPtCrossSectionHighStat,"S + #sigma_{S}(stat)","l");
  l2->AddEntry(hPtCrossSectionHighSyst,"S + #sigma_{S}(syst)","l");
  l2->AddEntry(hPtCrossSectionLowS,"S - #sigma_{S}","l");
  l2->AddEntry(hPtCrossSectionLowStat,"S - #sigma_{S}(stat)","l");
  l2->AddEntry(hPtCrossSectionLowSyst,"S - #sigma_{S}(syst)","l");

  TLine* line = new TLine(PtLims[0],1,PtLims[nPtBins],1);
  line->SetLineColor(colors[0]);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  
  TCanvas* cPromptCross = new TCanvas("cPromptCross","",10,10,1200,600);
  cPromptCross->Divide(2,1);
  cPromptCross->cd(1)->SetLogy();
  hPtCrossSection->GetXaxis()->SetTitleSize(0.04);
  hPtCrossSection->GetYaxis()->SetTitleSize(0.04);
  hPtCrossSection->Draw("hist");
  hPtCrossSectionLowS->Draw("hist same");
  hPtCrossSectionLowStat->Draw("hist same");
  hPtCrossSectionLowSyst->Draw("hist same");
  hPtCrossSectionHighS->Draw("hist same");
  hPtCrossSectionHighStat->Draw("hist same");
  hPtCrossSectionHighSyst->Draw("hist same");
  l->Draw("same");
  cPromptCross->cd(2);
  hPtRatioLowS->GetXaxis()->SetTitleSize(0.04);
  hPtRatioLowS->GetYaxis()->SetTitleSize(0.04);
  hPtRatioLowS->GetYaxis()->SetRangeUser(0.9,1.15);
  hPtRatioLowS->Draw("hist");
  hPtRatioLowStat->Draw("hist same");
  hPtRatioLowSyst->Draw("hist same");
  hPtRatioHighS->Draw("hist same");  
  hPtRatioHighStat->Draw("hist same");  
  hPtRatioHighSyst->Draw("hist same");  
  line->Draw("same");
  l2->Draw("same");

  TCanvas* cPromptCrossOnly = new TCanvas("cPromptCrossOnly","",10,10,800,800);
  cPromptCrossOnly->SetLogy();
  hPtCrossSection->GetXaxis()->SetTitleSize(0.05);
  hPtCrossSection->GetYaxis()->SetTitleSize(0.05);
  hPtCrossSection->GetXaxis()->SetLabelSize(0.05);
  hPtCrossSection->GetYaxis()->SetLabelSize(0.05);
  hPtCrossSection->Draw("hist");
  hPtCrossSectionLowS->Draw("hist same");
  hPtCrossSectionLowStat->Draw("hist same");
  hPtCrossSectionLowSyst->Draw("hist same");
  hPtCrossSectionHighS->Draw("hist same");
  hPtCrossSectionHighStat->Draw("hist same");
  hPtCrossSectionHighSyst->Draw("hist same");
  l->Draw("same");
  
  TCanvas* cRatio = new TCanvas("cRatio","",10,10,800,800);
  hPtRatioLowS->GetXaxis()->SetTitleSize(0.05);
  hPtRatioLowS->GetYaxis()->SetTitleSize(0.05);
  hPtRatioLowS->GetXaxis()->SetLabelSize(0.05);
  hPtRatioLowS->GetYaxis()->SetLabelSize(0.05);
  hPtRatioLowS->GetYaxis()->SetRangeUser(0.93,1.15);
  hPtRatioLowS->Draw("hist");
  hPtRatioLowStat->Draw("hist same");
  hPtRatioLowSyst->Draw("hist same");
  hPtRatioHighS->Draw("hist same");  
  hPtRatioHighStat->Draw("hist same");  
  hPtRatioHighSyst->Draw("hist same");  
  line->Draw("same");
  l2->Draw("same");
  
  cPromptCross->SaveAs("promptcrosssection_syst_SoverT.eps");
  cPromptCrossOnly->SaveAs("promptcrosssection_syst_SoverT_crossonly.eps");
  cRatio->SaveAs("promptcrosssection_syst_SoverT_ratioonly.eps");
  cPromptCross->SaveAs("promptcrosssection_syst_SoverT.root");
  
}

TH1F* CalculateCrossSection(TH1F* hPtRawYields,TH1F* hPtEffXAcc,TH1F* hPtPromptFrac,Double_t Nev,Double_t BR,Double_t sigma) {
  
  sigma = sigma*1000000;

  const Int_t nPtBins = hPtEffXAcc->GetNbinsX();
  const Int_t nPtLims = nPtBins+1;
  TArrayD* ptarray = (TArrayD*)hPtEffXAcc->GetXaxis()->GetXbins();
  Double_t* PtLims = (Double_t*)ptarray->GetArray();
  
  TH1F* hPtCrossSection = new TH1F("hPtCrossSection","",nPtBins,PtLims);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {

    Double_t crosssec = hPtRawYields->GetBinContent(iPt+1)*hPtPromptFrac->GetBinContent(iPt+1)*sigma/(2*(PtLims[iPt+1]-PtLims[iPt])*BR*Nev*hPtEffXAcc->GetBinContent(iPt+1));
    Double_t relerr = TMath::Sqrt((hPtRawYields->GetBinError(iPt+1)/hPtRawYields->GetBinContent(iPt+1)*hPtRawYields->GetBinError(iPt+1)/hPtRawYields->GetBinContent(iPt+1))+(hPtPromptFrac->GetBinError(iPt+1)/hPtPromptFrac->GetBinContent(iPt+1)*hPtPromptFrac->GetBinError(iPt+1)/hPtPromptFrac->GetBinContent(iPt+1))+(hPtEffXAcc->GetBinError(iPt+1)/hPtEffXAcc->GetBinContent(iPt+1)*hPtEffXAcc->GetBinError(iPt+1)/hPtEffXAcc->GetBinContent(iPt+1)));
    
    hPtCrossSection->SetBinContent(iPt+1,crosssec);
    hPtCrossSection->SetBinError(iPt+1,relerr*crosssec);
  }
  
  return hPtCrossSection;
  
}
