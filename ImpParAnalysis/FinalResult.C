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
#include <TGraphAsymmErrors.h>
#include <TArrayD.h>
#include <TArrow.h>
#include <TLatex.h>
#include <vector>
#endif

void DrawResult(Double_t linewidth=2.);
Double_t GetLowerLimit(Double_t mean, Double_t sigmastat, Double_t sigmasyst, Double_t CL=0.95);

void DrawResult(Double_t linewidth) {

  //SYSTEMATIC UNCERTAINTIES______________________________________________________________________________
  const Int_t nPtBins = 7;
  const Int_t nPtLims = nPtBins+1;
  Double_t PtLims[nPtLims] = {2.,3.,4.,5.,6.,8.,12.,16.};
  
  Double_t RangeSystLow[nPtBins] = {0.02,0.02,0.02,0.02,0.02,0.02,0.02};
  Double_t RangeSystHigh[nPtBins] = {0.02,0.02,0.02,0.02,0.02,0.02,0.02};
  Double_t MassRangeSystLow[nPtBins] = {0.,0.,0.,0.,0.,0.,0.};
  Double_t MassRangeSystHigh[nPtBins] = {0.,0.,0.,0.,0.,0.,0.};
  Double_t SigmaFixedSystLow[nPtBins] = {0.,0.,0.,0.,0.,0.,0.};
  Double_t SigmaFixedSystHigh[nPtBins] = {0.,0.,0.,0.,0.,0.,0.};
  Double_t PrefitSystLow[nPtBins] = {0.04,0.01,0.01,0.01,0.01,0.01,0.01};
  Double_t PrefitSystHigh[nPtBins] = {0.04,0.01,0.01,0.01,0.01,0.01,0.01};
  Double_t SBSystLow[nPtBins] = {0.06,0.03,0.03,0.03,0.03,0.03,0.03};
  Double_t SBSystHigh[nPtBins] = {0.06,0.03,0.03,0.03,0.03,0.03,0.03};
  Double_t SignalSystLow[nPtBins] = {0.04,0.02,0.02,0.02,0.02,0.02,0.02};
  Double_t SignalSystHigh[nPtBins] = {0.04,0.02,0.02,0.02,0.02,0.02,0.02};
  Double_t PtShapeSystLow[nPtBins] = {0.,0.,0.,0.,0.,0.,0.};
  Double_t PtShapeSystHigh[nPtBins] = {0.,0.,0.,0.,0.,0.,0.};
  Double_t MCTestSystLow[nPtBins] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  Double_t MCTestSystHigh[nPtBins] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01};  
  //______________________________________________________________________________________________________
  
  TFile FONLLfile("HFPtSpectrum_combinedFD.root","UPDATE");
  TGraphAsymmErrors* gfPromptFONLL = (TGraphAsymmErrors*)FONLLfile.Get("gFcCorrConservative");
  TGraphAsymmErrors* gfPromptFONLLPoints = (TGraphAsymmErrors*)gfPromptFONLL->Clone();
  FONLLfile.Close();

  Int_t nFONLLBins = gfPromptFONLL->GetN();

  for(Int_t iPt=0; iPt<nFONLLBins; iPt++) {
    gfPromptFONLLPoints->SetPointEYlow(iPt,0);
    gfPromptFONLLPoints->SetPointEYhigh(iPt,0);
  }
  
  TFile fPromptfile("fraction_unbinned_sigmafree.root","UPDATE");
  TH1F *hPromptFraction = (TH1F*)fPromptfile.Get("hFrac");
  hPromptFraction->SetDirectory(0);
  fPromptfile.Close();
  
  TGraphAsymmErrors* gSystRange = new  TGraphAsymmErrors(nPtBins);
  gSystRange->SetName("gSystRange");
  gSystRange->SetTitle("Range Systematics");
  gSystRange->SetFillStyle(20);
  gSystRange->SetLineWidth(linewidth);
  gSystRange->SetLineColor(kRed);

  TGraphAsymmErrors* gSystMass = new  TGraphAsymmErrors(nPtBins);
  gSystMass->SetName("gSystMass");
  gSystMass->SetTitle("Mass range Systematics");
  gSystMass->SetFillStyle(20);
  gSystRange->SetLineWidth(linewidth);
  gSystMass->SetLineColor(kBlue);
  
  TGraphAsymmErrors* gSystSignal = new  TGraphAsymmErrors(nPtBins);
  gSystSignal->SetName("gSystSignal");
  gSystSignal->SetTitle("Signal Systematics");
  gSystSignal->SetFillStyle(20);
  gSystSignal->SetLineWidth(linewidth);
  gSystSignal->SetLineColor(kGreen+3);

  TGraphAsymmErrors* gSystPrefit = new  TGraphAsymmErrors(nPtBins);
  gSystPrefit->SetName("gSystPrefit");
  gSystPrefit->SetTitle("Prefit Systematics");
  gSystPrefit->SetFillStyle(20);
  gSystPrefit->SetLineWidth(linewidth);
  gSystPrefit->SetLineColor(kYellow);

  TGraphAsymmErrors* gSystSB = new  TGraphAsymmErrors(nPtBins);
  gSystSB->SetName("gSystSB");
  gSystSB->SetTitle("SB range Systematics");
  gSystSB->SetFillStyle(20);
  gSystSB->SetLineWidth(linewidth);
  gSystSB->SetLineColor(kMagenta);

  TGraphAsymmErrors* gSystSigma = new TGraphAsymmErrors(nPtBins);
  gSystSigma->SetName("gSystSigma");
  gSystSigma->SetTitle("Sigma Fix Systematics");
  gSystSigma->SetFillStyle(20);
  gSystSigma->SetLineWidth(linewidth);
  gSystSigma->SetLineColor(kCyan);

  TGraphAsymmErrors* gSystPtShape = new TGraphAsymmErrors(nPtBins);
  gSystPtShape->SetName("gSystPtShape");
  gSystPtShape->SetTitle("Pt Shape Systematics");
  gSystPtShape->SetFillStyle(20);
  gSystPtShape->SetLineWidth(linewidth);
  gSystPtShape->SetLineColor(kOrange+7);

  TGraphAsymmErrors* gSystMCTest = new TGraphAsymmErrors(nPtBins);
  gSystPtShape->SetName("gSystMCTest");
  gSystPtShape->SetTitle("MC closure test Systematics");
  gSystPtShape->SetFillStyle(20);
  gSystPtShape->SetLineWidth(linewidth);
  gSystPtShape->SetLineColor(kOrange+3);
  
  TGraphAsymmErrors* gSyst = new TGraphAsymmErrors(nPtBins);
  gSyst->SetName("gSyst");
  gSyst->SetTitle("Systematic Uncertainties");
  gSyst->SetLineWidth(linewidth);
  gSyst->SetMarkerColor(kBlack);
  gSyst->SetFillStyle(20);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {

    Double_t PtCentValue = (PtLims[iPt+1]+PtLims[iPt])/2; 
    
    gSystRange->SetPoint(iPt,PtCentValue,0);
    gSystRange->SetPointError(iPt,0.5,0.5,RangeSystLow[iPt],RangeSystHigh[iPt]);
    gSystMass->SetPoint(iPt,PtCentValue,0);
    gSystMass->SetPointError(iPt,0.5,0.5,MassRangeSystLow[iPt],MassRangeSystHigh[iPt]);
    gSystSigma->SetPoint(iPt,PtCentValue,0);
    gSystSigma->SetPointError(iPt,0.5,0.5,SigmaFixedSystLow[iPt],SigmaFixedSystHigh[iPt]);
    gSystSignal->SetPoint(iPt,PtCentValue,0);
    gSystSignal->SetPointError(iPt,0.5,0.5,SignalSystLow[iPt],SignalSystHigh[iPt]);
    gSystPrefit->SetPoint(iPt,PtCentValue,0);
    gSystPrefit->SetPointError(iPt,0.5,0.5,PrefitSystLow[iPt],PrefitSystHigh[iPt]);
    gSystSB->SetPoint(iPt,PtCentValue,0);
    gSystSB->SetPointError(iPt,0.5,0.5,SBSystLow[iPt],SBSystHigh[iPt]);
    gSystPtShape->SetPoint(iPt,PtCentValue,0);
    gSystPtShape->SetPointError(iPt,0.5,0.5,PtShapeSystLow[iPt],PtShapeSystHigh[iPt]);
    gSystPtShape->SetPoint(iPt,PtCentValue,0);
    gSystPtShape->SetPointError(iPt,0.5,0.5,MCTestSystLow[iPt],MCTestSystHigh[iPt]);
    
    Double_t systlow = TMath::Sqrt(RangeSystLow[iPt]*RangeSystLow[iPt]+MassRangeSystLow[iPt]*MassRangeSystLow[iPt]+SigmaFixedSystLow[iPt]*SigmaFixedSystLow[iPt]+SignalSystLow[iPt]*SignalSystLow[iPt]+PrefitSystLow[iPt]*PrefitSystLow[iPt]+SBSystLow[iPt]*SBSystLow[iPt]+PtShapeSystLow[iPt]*PtShapeSystLow[iPt]+MCTestSystLow[iPt]*MCTestSystLow[iPt]);
    Double_t systhigh = TMath::Sqrt(RangeSystHigh[iPt]*RangeSystHigh[iPt]+MassRangeSystHigh[iPt]*MassRangeSystHigh[iPt]+SigmaFixedSystHigh[iPt]*SigmaFixedSystHigh[iPt]+SignalSystHigh[iPt]*SignalSystHigh[iPt]+PrefitSystHigh[iPt]*PrefitSystHigh[iPt]+SBSystHigh[iPt]*SBSystHigh[iPt]+PtShapeSystHigh[iPt]*PtShapeSystHigh[iPt]+MCTestSystHigh[iPt]*MCTestSystHigh[iPt]);
    
    gSyst->SetPoint(iPt,PtCentValue,0);
    gSyst->SetPointError(iPt,0.5,0.5,systlow,systhigh);
  }
  
  TGraphAsymmErrors* gPromptFractionSyst = new TGraphAsymmErrors(nPtBins);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    Double_t Pt = hPromptFraction->GetBinCenter(iPt+1);
    Double_t frac = hPromptFraction->GetBinContent(iPt+1);
    gPromptFractionSyst->SetPoint(iPt,Pt,frac);
    gPromptFractionSyst->SetPointError(iPt,gSyst->GetErrorXlow(iPt),gSyst->GetErrorXhigh(iPt),gSyst->GetErrorYlow(iPt)*frac,gSyst->GetErrorYhigh(iPt)*frac);
    cout << hPromptFraction->GetBinContent(iPt+1) << endl;
  }

  Bool_t IsOverLimit=kFALSE;
  Bool_t OverLimit[nPtBins];
  Double_t LowerLimit[nPtBins];
  TGraphAsymmErrors *gPromptFractionLimits = new TGraphAsymmErrors(nPtBins);
  TGraphAsymmErrors *gPromptFractionLimits2 = new TGraphAsymmErrors(nPtBins);
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    OverLimit[iPt] = kFALSE;
    if(hPromptFraction->GetBinContent(iPt+1)>1) {
      IsOverLimit = kTRUE;
      OverLimit[iPt] = kTRUE;
      Double_t PtCentValue = (PtLims[iPt+1]+PtLims[iPt])/2;
      LowerLimit[iPt] = GetLowerLimit(hPromptFraction->GetBinContent(iPt+1),hPromptFraction->GetBinError(iPt+1),gPromptFractionSyst->GetErrorYlow(iPt),0.95);

      gPromptFractionSyst->RemovePoint(iPt);
      hPromptFraction->SetBinContent(iPt+1,0);
      hPromptFraction->SetBinError(iPt+1,0);
      gPromptFractionLimits->SetPoint(iPt+1,PtCentValue,LowerLimit[iPt]);
      gPromptFractionLimits2->SetPoint(iPt+1,PtCentValue,LowerLimit[iPt]);
      gPromptFractionLimits->SetPointError(iPt+1,0.,0.,0.,0.05);
      gPromptFractionLimits2->SetPointError(iPt+1,0.5,0.5,0.,0.);
    }
  }

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    if(!OverLimit[iPt]) {
      gPromptFractionLimits->RemovePoint(iPt);
      gPromptFractionLimits2->RemovePoint(iPt);
    }
  }
  
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTickLength(0.02,"X");
  gStyle->SetTickLength(0.02,"Y"); 
  
  TLegend *l = new TLegend(0.15,0.18,0.5,0.42);
  l->AddEntry(gfPromptFONLL,"FONLL based method", "f");
  l->AddEntry(hPromptFraction,"Impact Parameter method", "lpe");
  l->SetTextSize(0.035);
  TLatex* latex = new TLatex();
  latex->SetTextFont(42);
  latex->SetTextSize(0.045);

  TCanvas *cPromptFrac = new TCanvas("cPromptFrac","",800,600);
  gfPromptFONLL->SetTitle("");
  gfPromptFONLL->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  gfPromptFONLL->GetYaxis()->SetTitle("#it{f}_{prompt}");
  gfPromptFONLL->GetYaxis()->SetRangeUser(0.1,1.2);
  gfPromptFONLL->GetXaxis()->SetLabelFont(42);
  gfPromptFONLL->GetYaxis()->SetLabelFont(42);
  gfPromptFONLL->GetXaxis()->SetTitleFont(42);
  gfPromptFONLL->GetYaxis()->SetTitleFont(42);
  gfPromptFONLL->GetXaxis()->SetTitleSize(0.045);
  gfPromptFONLL->GetYaxis()->SetTitleSize(0.045);
  gfPromptFONLL->SetLineWidth(linewidth);
  gfPromptFONLL->Draw("A2");
  gfPromptFONLLPoints->SetLineColor(kRed);
  gfPromptFONLLPoints->SetLineStyle(9);
  gfPromptFONLLPoints->SetLineWidth(linewidth);
  gfPromptFONLLPoints->SetMarkerSize(1);
  gfPromptFONLLPoints->Draw("Esame");
  hPromptFraction->SetStats(0);
  hPromptFraction->SetMarkerStyle(20);
  hPromptFraction->SetMarkerColor(kBlack);
  hPromptFraction->SetLineColor(kBlack);
  hPromptFraction->SetLineWidth(linewidth);
  hPromptFraction->Draw("Esame");
  gPromptFractionSyst->SetFillStyle(20);
  gPromptFractionSyst->SetLineWidth(linewidth);
  gPromptFractionSyst->Draw("2");
  if(IsOverLimit) {
    TArrow* arrow = new TArrow(3.6,0.4,3.6,0.45);
    arrow->SetLineWidth(linewidth);
    arrow->SetArrowSize(0.008);
    TLine* line = new TLine(3.1,0.4,4.1,0.4);
    line->SetLineWidth(linewidth);
    l->SetY1(0.35);
    l->SetY2(0.13);
    gPromptFractionLimits->SetMarkerStyle(9);
    gPromptFractionLimits->SetLineWidth(linewidth);
    gPromptFractionLimits->SetLineColor(kBlack);
    gPromptFractionLimits2->SetMarkerSize(0);
    gPromptFractionLimits2->SetLineWidth(linewidth);
    gPromptFractionLimits2->SetLineColor(kBlack);
    gPromptFractionLimits->Draw("E |> same");
    gPromptFractionLimits2->Draw("EZ same");
    arrow->Draw("|> same");
    line->Draw("same");
    TLatex* latex2 = new TLatex();
    latex2->SetTextFont(42);
    latex2->SetTextSize(0.035);
    latex2->DrawLatex(4.5,0.45,"Lower Limit at 95% C.L.");
  }
  else {
    delete gPromptFractionLimits;
    delete gPromptFractionLimits2;
  }
  latex->DrawLatex(2,0.6,"D^{+} #rightarrow K^{-}#pi^{+}#pi^{+}");
  latex->DrawLatex(16,1.1,"p-Pb #sqrt{s_{NN}} = 5.02 TeV");
  l->Draw("same");
  
  cPromptFrac->SaveAs("prompt_fraction.eps");
  cPromptFrac->SaveAs("prompt_fraction.pdf");
  cPromptFrac->SaveAs("prompt_fraction.png");

  TFile outfile("fprompt.root","RECREATE");
  gfPromptFONLL->Write();
  hPromptFraction->SetName("hPromptFraction");
  hPromptFraction->Write();
  gPromptFractionSyst->SetName("gPromptFractionSyst");
  gPromptFractionSyst->Write();
  if(IsOverLimit) {
    gPromptFractionLimits->SetName("gPromptFractionLimits");
    gPromptFractionLimits->Write();
  }
  cPromptFrac->Write();
  outfile.Close();
  
  TFile outfilesyst("fPrompt_syst_unc.root","RECREATE");
  gSystRange->Write();
  gSystMass->Write();
  gSystSignal->Write();
  gSystPrefit->Write();
  gSystSB->Write();
  gSystSigma->Write();
  gSystPtShape->Write();
  gSystMCTest->Write();
  gSyst->Write();
  outfilesyst.Close();
}

Double_t GetLowerLimit(Double_t mean, Double_t sigmastat, Double_t sigmasyst, Double_t CL) {
  
  Double_t sigma = TMath::Sqrt(sigmastat*sigmastat+sigmasyst*mean*mean*sigmasyst);
  Double_t nSigma = TMath::Sqrt(2)*TMath::ErfInverse(2*CL-1); //takes all the right integral of the gaussian and fix the nSigma to the left value which gives the chosen CL
  
  return mean-nSigma*sigma;
}
