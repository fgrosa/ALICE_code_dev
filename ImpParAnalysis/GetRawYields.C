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

#endif

const Int_t colors[] = {kRed,kBlue,kGreen+3,kBlack,kMagenta,kCyan,kOrange+7,kYellow-3,kCyan+3,kGreen};

void GetRawYields(Int_t nSigma = 3) {

  //______________________________________________________________________________________
  //pT bins of the analysis

  const Int_t nPtBins = 7;
  const Int_t nPtLims = nPtBins+1;
  Double_t PtLims[nPtLims] = {2,3,4,5,6,8,12,16};

  //______________________________________________________________________________________
  //get raw yields from files
  TH1F* hPtRawYields = new TH1F("hPtRawYields","raw yields",nPtBins,PtLims);
  hPtRawYields->SetDirectory(0);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    TFile massfile(Form("Mass_%0.f-%0.f.root",PtLims[iPt],PtLims[iPt+1]),"UPDATE");
    TCanvas* c = (TCanvas*)massfile.Get("cMass");
    TF1* fSignal = (TF1*)c->GetPrimitive("funcmass");
    TF1* fBkg = (TF1*)c->GetPrimitive("funcbkgFullRange");
    TH1F* hMass = (TH1F*)c->GetPrimitive("fhistoInvMass");
    
    Double_t mean = fSignal->GetParameter(3);
    Double_t sigma = fSignal->GetParameter(4);
    Double_t ints=fSignal->Integral(mean-nSigma*sigma,mean+nSigma*sigma)/(Double_t)hMass->GetBinWidth(4);
    Double_t intb=fBkg->Integral(mean-nSigma*sigma,mean+nSigma*sigma)/(Double_t)hMass->GetBinWidth(2);
    Double_t signal = ints-intb;
    Double_t signalerr = fSignal->GetParError(fSignal->GetNpar()-3)/fSignal->GetParameter(fSignal->GetNpar()-3)*signal;
    hPtRawYields->SetBinContent(iPt+1,signal);
    hPtRawYields->SetBinError(iPt+1,signalerr);
    massfile.Close();
  }

  //______________________________________________________________________________________
  //raw yields variations
  Double_t relsysterr = 0.05;
  
  TH1F* hPtRawYieldsLowS = new TH1F("hPtRawYieldsLowS","",nPtBins,PtLims);
  TH1F* hPtRawYieldsLowStat = new TH1F("hPtRawYieldsLowStat","",nPtBins,PtLims);
  TH1F* hPtRawYieldsLowSyst = new TH1F("hPtRawYieldsLowSyst","",nPtBins,PtLims);
  TH1F* hPtRawYieldsHighS = new TH1F("hPtRawYieldsHighS","",nPtBins,PtLims);
  TH1F* hPtRawYieldsHighStat = new TH1F("hPtRawYieldsHighStat","",nPtBins,PtLims);
  TH1F* hPtRawYieldsHighSyst = new TH1F("hPtRawYieldsHighSyst","",nPtBins,PtLims);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    if(PtLims[iPt]>=12)
      relsysterr = 0.08;
    Double_t signal = hPtRawYields->GetBinContent(iPt+1);
    Double_t staterr = hPtRawYields->GetBinError(iPt+1);
    Double_t err = TMath::Sqrt(staterr*staterr+relsysterr*signal*relsysterr*signal);
    hPtRawYieldsLowS->SetBinContent(iPt+1,signal-err);
    hPtRawYieldsLowStat->SetBinContent(iPt+1,signal-staterr);
    hPtRawYieldsLowSyst->SetBinContent(iPt+1,signal-relsysterr*signal);
    hPtRawYieldsHighS->SetBinContent(iPt+1,signal+err);
    hPtRawYieldsHighStat->SetBinContent(iPt+1,signal+staterr);
    hPtRawYieldsHighSyst->SetBinContent(iPt+1,signal+relsysterr*signal);
  }

  TH1F* hPtRatioLowS = new TH1F("hPtRatioLowS","",nPtBins,PtLims);
  hPtRatioLowS->Divide(hPtRawYieldsLowS,hPtRawYields,1.,1.);
  TH1F* hPtRatioLowStat = new TH1F("hPtRatioLowStat","",nPtBins,PtLims);
  hPtRatioLowStat->Divide(hPtRawYieldsLowStat,hPtRawYields,1.,1.);
  TH1F* hPtRatioLowSyst = new TH1F("hPtRatioLowSyst","",nPtBins,PtLims);
  hPtRatioLowSyst->Divide(hPtRawYieldsLowSyst,hPtRawYields,1.,1.);
  TH1F* hPtRatioHighS = new TH1F("hPtRatioHighS","",nPtBins,PtLims);
  hPtRatioHighS->Divide(hPtRawYieldsHighS,hPtRawYields,1.,1.);
  TH1F* hPtRatioHighStat = new TH1F("hPtRatioHighStat","",nPtBins,PtLims);
  hPtRatioHighStat->Divide(hPtRawYieldsHighStat,hPtRawYields,1.,1.);
  TH1F* hPtRatioHighSyst = new TH1F("hPtRatioHighSyst","",nPtBins,PtLims);
  hPtRatioHighSyst->Divide(hPtRawYieldsHighSyst,hPtRawYields,1.,1.);
  
  //_______________________________________________________________________________________
  //plots
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.045,"xyz");
  gStyle->SetTitleOffset(1.3,"y");
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);

  hPtRawYields->SetLineColor(colors[0]);
  hPtRawYieldsLowS->SetLineColor(colors[4]);
  hPtRawYieldsLowStat->SetLineColor(colors[5]);
  hPtRawYieldsLowSyst->SetLineColor(colors[6]);
  hPtRawYieldsHighS->SetLineColor(colors[1]);
  hPtRawYieldsHighStat->SetLineColor(colors[2]);
  hPtRawYieldsHighSyst->SetLineColor(colors[3]);
  hPtRatioLowS->SetLineColor(colors[4]);
  hPtRatioLowStat->SetLineColor(colors[5]);
  hPtRatioLowSyst->SetLineColor(colors[6]);
  hPtRatioHighS->SetLineColor(colors[1]);
  hPtRatioHighStat->SetLineColor(colors[2]);
  hPtRatioHighSyst->SetLineColor(colors[3]);

  TLegend *l = new TLegend(0.55,0.6,0.85,0.85);
  l->AddEntry(hPtRawYields,"Reference value","l");
  l->AddEntry(hPtRawYieldsLowS,"S - #sigma_{S}","l");
  l->AddEntry(hPtRawYieldsLowStat,"S - #sigma_{S}(stat)","l");
  l->AddEntry(hPtRawYieldsLowSyst,"S - #sigma_{S}(syst)","l");
  l->AddEntry(hPtRawYieldsHighS,"S + #sigma_{S}","l");
  l->AddEntry(hPtRawYieldsHighStat,"S + #sigma_{S}(stat)","l");
  l->AddEntry(hPtRawYieldsHighSyst,"S + #sigma_{S}(syst)","l");
  l->SetTextSize(0.035);

  TLegend *l2 = new TLegend(0.55,0.15,0.89,0.4);
  l2->AddEntry(hPtRawYields,"Reference value","l");
  l2->AddEntry(hPtRawYieldsLowS,"S - #sigma_{S}","l");
  l2->AddEntry(hPtRawYieldsLowStat,"S - #sigma_{S}(stat)","l");
  l2->AddEntry(hPtRawYieldsLowSyst,"S - #sigma_{S}(syst)","l");
  l2->AddEntry(hPtRawYieldsHighS,"S + #sigma_{S}","l");
  l2->AddEntry(hPtRawYieldsHighStat,"S + #sigma_{S}(stat)","l");
  l2->AddEntry(hPtRawYieldsHighSyst,"S + #sigma_{S}(syst)","l");
  l2->SetTextSize(0.035);
  
  TLine *line = new TLine(PtLims[0],1,PtLims[nPtBins],1);
  line->SetLineColor(colors[0]);
  
  TCanvas* cRaw = new TCanvas("cRaw","",1200,600);
  cRaw->Divide(2,1);
  cRaw->cd(1);
  hPtRawYields->GetYaxis()->SetTitleOffset(1.3);
  hPtRawYields->GetYaxis()->SetRangeUser(hPtRawYields->GetMinimum()*0.4,hPtRawYields->GetMaximum()*1.2);
  hPtRawYields->GetYaxis()->SetTitle("Raw Yields");
  hPtRawYields->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hPtRawYields->Draw("hist");
  hPtRawYieldsLowS->Draw("hist same");
  hPtRawYieldsLowStat->Draw("hist same");
  hPtRawYieldsLowSyst->Draw("hist same");
  hPtRawYieldsHighS->Draw("hist same");
  hPtRawYieldsHighStat->Draw("hist same");
  hPtRawYieldsHighSyst->Draw("hist same");
  l->Draw("same");
  cRaw->cd(2);
  hPtRatioLowS->GetYaxis()->SetRangeUser(0.6,1.2);
  hPtRatioLowS->GetYaxis()->SetTitle("ratio of raw yield w.r.t. the central value");
  hPtRatioLowS->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hPtRatioLowS->Draw("hist");
  hPtRatioLowStat->Draw("hist same");
  hPtRatioLowSyst->Draw("hist same");
  hPtRatioHighS->Draw("hist same");
  hPtRatioHighStat->Draw("hist same");
  hPtRatioHighSyst->Draw("hist same");
  line->Draw("same");
  l2->Draw("same");

  TCanvas* cRawOnly = new TCanvas("cRawOnly","",800,800);
  hPtRawYields->Draw("hist");
  hPtRawYieldsLowS->Draw("hist same");
  hPtRawYieldsLowStat->Draw("hist same");
  hPtRawYieldsLowSyst->Draw("hist same");
  hPtRawYieldsHighS->Draw("hist same");
  hPtRawYieldsHighStat->Draw("hist same");
  hPtRawYieldsHighSyst->Draw("hist same");
  l->Draw("same");
  TCanvas* cRatio = new TCanvas("cRatio","",800,800);
  hPtRatioLowS->Draw("hist");
  hPtRatioLowStat->Draw("hist same");
  hPtRatioLowSyst->Draw("hist same");
  hPtRatioHighS->Draw("hist same");
  hPtRatioHighStat->Draw("hist same");
  hPtRatioHighSyst->Draw("hist same");
  line->Draw("same");
  l2->Draw("same");
  
  cRaw->SaveAs("RawYields.eps");
  cRawOnly->SaveAs("RawYieldsOnly.eps");
  cRatio->SaveAs("RawYieldsRatio.eps");
  cRaw->SaveAs("RawYields.root");

  delete cRawOnly;
  delete cRatio;
  
  TFile outfile("rawyields.root","RECREATE");
  hPtRawYields->Write();
  hPtRawYieldsLowS->Write();
  hPtRawYieldsLowStat->Write();
  hPtRawYieldsLowSyst->Write();
  hPtRawYieldsHighS->Write();
  hPtRawYieldsHighStat->Write();
  hPtRawYieldsHighSyst->Write();
  outfile.Close();
  
}
