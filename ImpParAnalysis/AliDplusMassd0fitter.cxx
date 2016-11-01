#include "AliDplusMassd0fitter.h"

#include <Riostream.h>
#include <TMath.h>
#include <TAxis.h>
#include <TF1.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TVirtualFitter.h>
#include <TH2F.h>
#include <TArrayD.h>
#include <TMatrixD.h>
#include <vector>

#include "AliHFMassFitter.h"

/// \cond CLASSIMP
ClassImp(AliDplusMassd0fitter);
/// \endcond

AliDplusMassd0fitter::AliDplusMassd0fitter()
  :TObject()
  ,fMCPromptSparse(0x0)
  ,fMCTrueFDSparse(0x0)
  ,fMCRecoFDSparse(0x0)
  ,fMCBkgSparse(0x0)
  ,fDataSparse(0x0)
  ,fDataTree(0x0)
  ,fImpParHisto(0x0)
  ,fMassBranch("InvMass")
  ,fPtBranch("Pt")
  ,fImpParBranch("d0")
  ,fPtMin(0)
  ,fPtMax(0)
  ,fGenPromptFraction(0.9)
  ,fMassMean(0)
  ,fMassSigma(0)
  ,fNSigma(1.5)
  ,fFitOptions("REM")
  ,fFitOptionsPrefit("RELI")
  ,fReb(1)
  ,fNSigmaSBLow(3)
  ,fNSigmaSBHigh(15)
  ,fInitParamPrompt(0)
  ,fInitParamFD(0)
  ,fInitParamBkg(0)
  ,fSig(0)
  ,fSigErr(0)
  ,fSigErrStat(0)
  ,fSigRelErrSyst(0.05)
  ,fBkg(0)
  ,fBkgErrStat(0)
  ,fIntegral(0)
  ,fPromptFraction(0)
  ,fPromptFractionGauss(0)
  ,fPromptMean(0)
  ,fPromptSigma(0)
  ,fPromptSigmaMC(0)
  ,fPromptLambda(0)
  ,fFDFraction1(0)
  ,fFDMean(0)
  ,fFDLambda1(0)
  ,fFDLambda2(0)
  ,fBkgFractionGauss1(0)
  ,fBkgMean1(0)
  ,fBkgSigma1(0)
  ,fBkgLambda1(0)
  ,fBkgFractionGauss2(0)
  ,fBkgMean2(0)
  ,fBkgSigma2(0)
  ,fBkgLambda2(0)
  ,fBkgFracFunc1(0)
  ,fPromptFractionErr(0)
  ,fPromptFractionGaussErr(0)
  ,fPromptMeanErr(0)
  ,fPromptSigmaErr(0)
  ,fPromptSigmaMCErr(0)
  ,fPromptLambdaErr(0)
  ,fFDFraction1Err(0)
  ,fFDMeanErr(0)
  ,fFDLambda1Err(0)
  ,fFDLambda2Err(0)
  ,fBkgFractionGauss1Err(0)
  ,fBkgMean1Err(0)
  ,fBkgSigma1Err(0)
  ,fBkgLambda1Err(0)
  ,fBkgFractionGauss2Err(0)
  ,fBkgMean2Err(0)
  ,fBkgSigma2Err(0)
  ,fBkgLambda2Err(0)
  ,fBkgFracFunc1Err(0)
  ,fCovFracSigmaPrompt(0)
  ,fSigmaPromptFixed(kFALSE) 
  ,fSigmaPromptLim(0.2)
  ,fPromptFracMin(0.2) 
  ,fPromptFracMax(1.5)
  ,fChiSquare(0)
  ,fNDF(0)
  ,fMCTest(kFALSE)
  ,fPID(kTRUE)
  ,fGaussOnlyForPrompt(kFALSE)
  ,fBkgFromMC(kFALSE)
  ,fFDType(kConvolution)
  ,fBkgType(kDoubleGaussExpoSymm)
  ,fSBRegion(kBoth)
  ,fPtBWeight(kFALSE)
  ,fPtDWeight(kFALSE)
  ,fPtBWeightsHisto(0x0)
  ,fPtDWeightsHisto(0x0)
  ,fMassAxis(0)
  ,fPtAxis(1)
  ,fImpParAxis(2)
  ,fPIDAxis(3)
  ,fPIDbin(2)
  ,fVariableBinning(kFALSE)
  ,fnCountsMin(0)
  ,fSubBkg(kFALSE)
{
  /// Default constructor
}

//_________________________________________________________________________________
AliDplusMassd0fitter::~AliDplusMassd0fitter()
{
  if(fMCPromptSparse) delete fMCPromptSparse;
  if(fMCTrueFDSparse) delete fMCTrueFDSparse;
  if(fMCRecoFDSparse) delete fMCRecoFDSparse;
  if(fMCBkgSparse) delete fMCBkgSparse;
  if(fDataSparse) delete fDataSparse;
  if(fDataTree) delete fDataTree;
  if(fImpParHisto) delete fImpParHisto;
  if(fPtBWeightsHisto) delete fPtBWeightsHisto;
  if(fPtDWeightsHisto) delete fPtDWeightsHisto;
  
  /// Default destructor
}

//_________________________________________________________________________________
void AliDplusMassd0fitter::GetSignal(Int_t rebin, Int_t funcsig, Int_t funcbkg, Double_t min, Double_t max, Int_t sovert)
{
  if(!fDataSparse) {
    cerr << "Data THnSparse not set!!" << endl;
    return;
  }

  ResetAxes(fDataSparse);
  SetPtRange(fDataSparse);
  
  TH1F* hMass=(TH1F*)fDataSparse->Projection(fMassAxis);
  hMass->Rebin(rebin);
  
  AliHFMassFitter* fitter = new AliHFMassFitter(hMass,min,max,1,funcbkg,funcsig);  
  fitter->MassFitter(kFALSE);

  fMassMean = fitter->GetMean();
  fMassSigma = fitter->GetSigma();
  fitter->Signal(fNSigma,fSig,fSigErrStat);
  fitter->Background(fNSigma,fBkg,fBkgErrStat); //only used in the TOY MC to generate the entries (otherwise B = N - S)
  
  if(fPtMin < 1)
    fSigRelErrSyst = 0.1;
  if(fPtMin <= 2 || fPtMin >= 12)
    fSigRelErrSyst = 0.08;
  
  fSigErr = TMath::Sqrt(fSigErrStat*fSigErrStat+fSig*fSigRelErrSyst*fSig*fSigRelErrSyst);
  
  if(sovert==kUpperLimit) {
    fSig = fSig+fSigErr;
  }
  if(sovert==kLowerLimit) {
    fSig = fSig-fSigErr;
  }
  if(sovert==kUpperLimitOnlyStat) {
    fSig = fSig+fSigErrStat;
  }
  if(sovert==kLowerLimitOnlyStat) {
    fSig = fSig-fSigErrStat;
  }
  if(sovert==kUpperLimitOnlySyst) {
    fSig = fSig+fSig*fSigRelErrSyst;
  }
  if(sovert==kLowerLimitOnlySyst) {
    fSig = fSig-fSig*fSigRelErrSyst;
  }
  
  gStyle->SetOptFit(0);
  TCanvas *cMass = new TCanvas("cMass","cMass",900,900);  
  cMass->Clear();
  fitter->DrawHere(gPad);
  cMass->SaveAs(Form("Mass_%0.f-%0.f.pdf",fPtMin,fPtMax));
  cMass->SaveAs(Form("Mass_%0.f-%0.f.root",fPtMin,fPtMax));

  delete hMass;
  delete cMass;
  delete fitter;
}

//_________________________________________________________________________________
void AliDplusMassd0fitter::CheckSideBandsImpParDist(Int_t nSigmaWidth)
{
  ResetAxes(fDataSparse);
  
  TAxis* massaxis=(TAxis*)fDataSparse->GetAxis(fMassAxis);
  Int_t massbinmin;
  Int_t massbinmax;
  
  const Int_t nSideBands = 4;
  Int_t nSigma[nSideBands] = {4,6,8,10};
  Int_t color[nSideBands] = {1,2,3,4};
  TH1F** hImpParSideBands = new TH1F*[nSideBands];
  TH1F** hRatioSideBands = new TH1F*[nSideBands-1];

  TCanvas* cSideBands = new TCanvas("cSideBands","cSideBands",1200,800);
  cSideBands->Clear();
  
  TCanvas* cSideBandsRatio = new TCanvas("cSideBandsRatio","cSideBandsRatio",1200,800);
  cSideBandsRatio->Clear();

  TLegend* l = new TLegend(0.6,0.65,0.89,0.89);
  l->SetTextSize(0.04);
  TLegend* l2 = new TLegend(0.35,0.7,0.65,0.89);
  l->SetTextSize(0.04);
  
  for(Int_t iSide=0; iSide<nSideBands; iSide++) {
    massbinmin=massaxis->FindBin((fMassMean+nSigma[iSide]*fMassSigma)*1.0001);
    massbinmax=massaxis->FindBin((fMassMean+(nSigma[iSide]+nSigmaWidth)*fMassSigma)*0.9999);
    fDataSparse->GetAxis(fMassAxis)->SetRange(massbinmin,massbinmax);
    TH1F* hImpParSideBands1=(TH1F*)fDataSparse->Projection(fImpParAxis);
    hImpParSideBands1->Sumw2();
    hImpParSideBands1->Rebin(fReb);

    massbinmax=massaxis->FindBin((fMassMean-nSigma[iSide]*fMassSigma)*0.9999);
    massbinmin=massaxis->FindBin((fMassMean-(nSigma[iSide]+nSigmaWidth)*fMassSigma)*1.0001);
    fDataSparse->GetAxis(fMassAxis)->SetRange(massbinmin,massbinmax);
    TH1F* hImpParSideBands2=(TH1F*)fDataSparse->Projection(fImpParAxis);
    hImpParSideBands2->Sumw2();
    hImpParSideBands2->Rebin(fReb);

    hImpParSideBands[iSide]=new TH1F(Form("hImpParSideBands%d",iSide),"",hImpParSideBands1->GetNbinsX(),hImpParSideBands1->GetBinLowEdge(0),hImpParSideBands1->GetBinLowEdge(hImpParSideBands1->GetNbinsX())+hImpParSideBands1->GetBinWidth(1));
  
    for(Int_t iBin=0; iBin<hImpParSideBands1->GetNbinsX(); iBin++){
      hImpParSideBands[iSide]->SetBinContent(iBin+1,hImpParSideBands1->GetBinContent(iBin+1)+hImpParSideBands2->GetBinContent(iBin+1));
    }
    hImpParSideBands[iSide]->Sumw2();
    hImpParSideBands[iSide]->Scale(1./hImpParSideBands[iSide]->Integral());
    hImpParSideBands[iSide]->SetLineColor(color[iSide]);
    hImpParSideBands[iSide]->GetXaxis()->SetTitle("Imp Par XY (#mum)");
    hImpParSideBands[iSide]->GetYaxis()->SetTitle("Normalised Entries");
    
    l->AddEntry(hImpParSideBands[iSide], Form("n#sigma = %d, #DeltaM = %d#sigma",nSigma[iSide],nSigmaWidth),"l");

    cSideBands->cd()->SetLogy();
    if(iSide==0)
      hImpParSideBands[iSide]->Draw();
    else
      hImpParSideBands[iSide]->Draw("same");
    l->Draw("same");

    if(iSide>0) {
      hRatioSideBands[iSide-1]=new TH1F(Form("hRatioSideBands%d",iSide),"",hImpParSideBands[iSide]->GetNbinsX(),hImpParSideBands[iSide]->GetBinLowEdge(0),hImpParSideBands[iSide]->GetBinLowEdge(hImpParSideBands[iSide]->GetNbinsX())+hImpParSideBands[iSide]->GetBinWidth(1));
      hRatioSideBands[iSide-1]->Divide(hImpParSideBands[iSide],hImpParSideBands[0],1.,1.);
      hRatioSideBands[iSide-1]->SetLineColor(color[iSide]);
      hRatioSideBands[iSide-1]->GetXaxis()->SetTitle("Imp Par XY (#mum)");
      hRatioSideBands[iSide-1]->GetYaxis()->SetTitle("n#sigma/4#sigma");
    
      l2->AddEntry(hRatioSideBands[iSide-1], Form("n#sigma = %d, #DeltaM = %d#sigma",nSigma[iSide],nSigmaWidth),"lpe");

      cSideBandsRatio->cd();
      hRatioSideBands[iSide-1]->GetYaxis()->SetRangeUser(0.,hRatioSideBands[iSide-1]->GetMaximum()+3);
      if(iSide==1)
      hRatioSideBands[iSide-1]->Draw("E1");
    else
      hRatioSideBands[iSide-1]->Draw("E1same");
      l2->Draw("same");
    }
  }
 
  cSideBands->SaveAs(Form("SideBandsDist_%0.f-%0.f.pdf",fPtMin,fPtMax));
  cSideBands->SaveAs(Form("SideBandsDist_%0.f-%0.f.root",fPtMin,fPtMax));
  cSideBandsRatio->SaveAs(Form("SideBandsRatio_%0.f-%0.f.pdf",fPtMin,fPtMax));
  cSideBandsRatio->SaveAs(Form("SideBandsRatio_%0.f-%0.f.root",fPtMin,fPtMax));

  delete cSideBands;
  delete cSideBandsRatio;
  delete l;
  delete l2;
  
  for(Int_t iSide=0; iSide<nSideBands; iSide++) {
    delete hImpParSideBands[iSide];
    if(iSide!=nSideBands-1)
      delete hRatioSideBands[iSide];
  }
  delete[] hRatioSideBands;
  delete[] hImpParSideBands;
}

//__________________________________________________________________________________
void AliDplusMassd0fitter::PrefitOnPrompt(Double_t d0minPrompt, Double_t d0maxPrompt)
{
  gStyle->SetOptFit(1);
  gStyle->SetTextSize(0.045);

  if(!fMCPromptSparse) {
    cerr << "MC Prompt THnSparse not set!!" << endl;
    return;
  }

  ResetAxes(fMCPromptSparse);
  
  TH1F* hImpParPrompt=0x0;
  if(fPtDWeight) {
    hImpParPrompt = (TH1F*)GetPtreweightedHisto(fMCPromptSparse,kD);
  }
  else {
    SetPtRange(fMCPromptSparse);
    SetMassRange(fMCPromptSparse);
    hImpParPrompt = (TH1F*)fMCPromptSparse->Projection(fImpParAxis);
  }
  hImpParPrompt->Rebin(fReb);
  hImpParPrompt->SetMarkerStyle(21);
  hImpParPrompt->SetMarkerSize(1.);
  hImpParPrompt->SetLineWidth(2);
  
  TF1* ImpParPromptFunc = new TF1("ImpParPromptFunc",this,&AliDplusMassd0fitter::FunctionImpParPrompt,d0minPrompt,d0maxPrompt,5,"AliDplusMassd0fitter","FunctionImpParPrompt");
  ImpParPromptFunc->SetLineColor(kGreen+3);
  ImpParPromptFunc->SetParNames("fractionGauss","mean","sigma","tailPrompt","integral");
  ImpParPromptFunc->SetParameters(fInitParamPrompt);
  ImpParPromptFunc->FixParameter(4,hImpParPrompt->Integral()*hImpParPrompt->GetBinWidth(1));
  ImpParPromptFunc->SetParLimits(0,0.6,1.);
  ImpParPromptFunc->SetParLimits(2,0.000001,10000.);
  ImpParPromptFunc->SetParLimits(3,0.00001,10000.);
  ImpParPromptFunc->SetParLimits(1,-10.,10.);
  if(fGaussOnlyForPrompt)
    ImpParPromptFunc->FixParameter(0,1.);

  TCanvas *cPrompt = new TCanvas("cPrompt","cPrompt",900,900);
  cPrompt->SetLogy();
  hImpParPrompt->SetLineColor(kBlack);
  hImpParPrompt->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",fPtMin,fPtMax));
  hImpParPrompt->GetXaxis()->SetTitle("Imp Par XY (#mum)");
  hImpParPrompt->GetYaxis()->SetTitle(Form("Entries/(%0.f #mum)",hImpParPrompt->GetBinWidth(5)));
  hImpParPrompt->GetXaxis()->SetNdivisions(505);
  hImpParPrompt->SetTitle("");
  hImpParPrompt->Draw("E1");
  TVirtualFitter::SetDefaultFitter("Minuit2");
  hImpParPrompt->Fit("ImpParPromptFunc",fFitOptionsPrefit.Data());
  fPromptFractionGauss = ImpParPromptFunc->GetParameter(0);
  fPromptMean = ImpParPromptFunc->GetParameter(1);
  fPromptSigmaMC = ImpParPromptFunc->GetParameter(2);
  fPromptLambda = ImpParPromptFunc->GetParameter(3);
  fPromptFractionGaussErr = ImpParPromptFunc->GetParError(0);
  fPromptMeanErr = ImpParPromptFunc->GetParError(1);
  fPromptSigmaMCErr = ImpParPromptFunc->GetParError(2);
  fPromptLambdaErr = ImpParPromptFunc->GetParError(3);
  cPrompt->SaveAs(Form("ImpParPrompt_%0.f-%0.f.pdf",fPtMin,fPtMax));
  cPrompt->SaveAs(Form("ImpParPrompt_%0.f-%0.f.root",fPtMin,fPtMax));
  delete cPrompt;
  delete ImpParPromptFunc;
} 

//_________________________________________________________________________________
void AliDplusMassd0fitter::PrefitOnFD(Double_t d0minFD, Double_t d0maxFD)
{
  gStyle->SetOptFit(1);
  gStyle->SetTextSize(0.045);
  gStyle->SetTitleSize(0.05,"xy");
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.12);
 
  if(!fMCTrueFDSparse) {
    cerr << "MC true FD THnSparse not set!!" << endl;
    return;
  }

  if(!fMCRecoFDSparse) {
    cerr << "MC reco FD THnSparse not set!!" << endl;
    return;
  }

  ResetAxes(fMCTrueFDSparse);
  ResetAxes(fMCRecoFDSparse);

  TH1F* hImpParTrueFD=0x0;
  TH1F* hImpParRecoFD=0x0;
  if(fPtBWeight) {
    hImpParTrueFD=(TH1F*)GetPtreweightedHisto(fMCTrueFDSparse,kB);
    hImpParRecoFD=(TH1F*)GetPtreweightedHisto(fMCRecoFDSparse,kB);
  }
  else{
    SetPtRange(fMCTrueFDSparse);
    SetPtRange(fMCRecoFDSparse);
    hImpParTrueFD=(TH1F*)fMCTrueFDSparse->Projection(fImpParAxis);
    hImpParRecoFD=(TH1F*)fMCRecoFDSparse->Projection(fImpParAxis);  
  }
  hImpParTrueFD->Rebin(fReb);
  hImpParRecoFD->Rebin(fReb);
  hImpParTrueFD->SetMarkerStyle(21);
  hImpParTrueFD->SetMarkerSize(1.);
  hImpParTrueFD->SetLineWidth(2);
  hImpParRecoFD->SetMarkerStyle(21);
  hImpParRecoFD->SetMarkerSize(1.);
  hImpParRecoFD->SetLineWidth(2);
  
  TF1* ImpParRecoFDFunc=0x0;

  if(fFDType==kConvolution) {
    TF1* ImpParTrueFDFunc = new TF1("ImpParTrueFDFunc",this,&AliDplusMassd0fitter::FunctionTrueImpParFD,d0minFD,d0maxFD,5,"AliDplusMassd0fitter","FunctionTrueImpParFD");
    ImpParTrueFDFunc->SetLineColor(kBlue);
    ImpParTrueFDFunc->SetParNames("fracTailfromB","mean","lambda1","lambda2","integral");
    ImpParTrueFDFunc->SetParameters(fInitParamFD);
    ImpParTrueFDFunc->FixParameter(4,hImpParTrueFD->Integral()*hImpParTrueFD->GetBinWidth(1));
    ImpParTrueFDFunc->SetParLimits(0,0.01,1.);
    ImpParTrueFDFunc->SetParError(0,0.005);
    ImpParTrueFDFunc->SetParLimits(2,0.000001,1000.);
    ImpParTrueFDFunc->SetParError(2,0.5);
    ImpParTrueFDFunc->SetParLimits(3,0.000001,1000.);
    ImpParTrueFDFunc->SetParError(3,0.5);
    
    TCanvas *cTrueFD = new TCanvas("cTrueFD","cTrueFD",900,900);
    cTrueFD->SetLogy();
    hImpParTrueFD->SetLineColor(kBlack);
    hImpParTrueFD->SetTitle("");
    hImpParTrueFD->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",fPtMin,fPtMax));
    hImpParTrueFD->GetXaxis()->SetTitle("Imp Par XY (#mum)");
    hImpParTrueFD->GetYaxis()->SetTitle(Form("Entries/(%0.f #mum)",hImpParTrueFD->GetBinWidth(5)));
    hImpParTrueFD->Draw("E1");
    hImpParTrueFD->GetXaxis()->SetNdivisions(505);
    TVirtualFitter::SetDefaultFitter("Minuit2");
    hImpParTrueFD->Fit("ImpParTrueFDFunc",fFitOptionsPrefit.Data());
    fFDFraction1 = ImpParTrueFDFunc->GetParameter(0);
    fFDMean = ImpParTrueFDFunc->GetParameter(1);
    fFDLambda1 = ImpParTrueFDFunc->GetParameter(2);
    fFDLambda2 = ImpParTrueFDFunc->GetParameter(3);
    fFDFraction1Err = ImpParTrueFDFunc->GetParError(0);
    fFDMeanErr = ImpParTrueFDFunc->GetParError(1);
    fFDLambda1Err = ImpParTrueFDFunc->GetParError(2);
    fFDLambda2Err = ImpParTrueFDFunc->GetParError(3);
    cTrueFD->SaveAs(Form("ImpParTrueFD_%0.f-%0.f.pdf",fPtMin,fPtMax));
    cTrueFD->SaveAs(Form("ImpParTrueFD_%0.f-%0.f.root",fPtMin,fPtMax));
    delete cTrueFD;
    
    ImpParRecoFDFunc = new TF1("ImpParRecoFDFunc",this,&AliDplusMassd0fitter::FunctionRecoImpParFD,-1000,1000,1,"AliDplusMassd0fitter","FunctionRecoImpParFD");
    ImpParRecoFDFunc->SetLineColor(kBlue);
    ImpParRecoFDFunc->FixParameter(0,hImpParRecoFD->Integral()*hImpParRecoFD->GetBinWidth(1));
    TCanvas *cRecoFD = new TCanvas("cRecoFD","cRecoFD",900,900);
    cRecoFD->SetLogy();
    hImpParRecoFD->SetLineColor(kBlack);
    hImpParRecoFD->Draw("E1");
    hImpParRecoFD->SetTitle("");
    hImpParRecoFD->GetXaxis()->SetTitle("Imp Par XY (#mum)");
    hImpParRecoFD->GetYaxis()->SetTitle(Form("Entries/(%0.f #mum)",hImpParRecoFD->GetBinWidth(5)));
    hImpParRecoFD->GetXaxis()->SetNdivisions(505);
    ImpParRecoFDFunc->Draw("same");
    cRecoFD->SaveAs(Form("ImpParRecoFD_%0.f-%0.f.pdf",fPtMin,fPtMax));
    cRecoFD->SaveAs(Form("ImpParRecoFD_%0.f-%0.f.root",fPtMin,fPtMax));
    delete cRecoFD;
    delete ImpParRecoFDFunc;
    delete ImpParTrueFDFunc;
  }
  else if (fFDType==kGaussExpo) {
    ImpParRecoFDFunc = new TF1("ImpParRecoFDFunc",this,&AliDplusMassd0fitter::FunctionImpParPrompt,d0minFD,d0maxFD,5,"AliDplusMassd0fitter","FunctionRecoImpParPrompt");
    ImpParRecoFDFunc->SetLineColor(kBlue);
    ImpParRecoFDFunc->SetParNames("fractionGauss","mean","sigma","lambdaFD","integral");
    ImpParRecoFDFunc->SetParameters(fInitParamFD);
    ImpParRecoFDFunc->FixParameter(4,hImpParRecoFD->Integral()*hImpParRecoFD->GetBinWidth(1));
    ImpParRecoFDFunc->SetParLimits(0,0.01,1.);
    ImpParRecoFDFunc->SetParError(0,0.005);
    ImpParRecoFDFunc->SetParLimits(2,0.000001,1000.);
    ImpParRecoFDFunc->SetParError(2,0.5);
    ImpParRecoFDFunc->SetParLimits(3,0.000001,1000.);
    ImpParRecoFDFunc->SetParError(3,0.5);
    TVirtualFitter::SetDefaultFitter("Minuit2");
    hImpParRecoFD->Fit("ImpParRecoFDFunc","RELI");
    fFDFraction1 = ImpParRecoFDFunc->GetParameter(0);
    fFDMean = ImpParRecoFDFunc->GetParameter(1);
    fFDLambda1 = ImpParRecoFDFunc->GetParameter(2);
    fFDLambda2 = ImpParRecoFDFunc->GetParameter(3);
    fFDFraction1Err = ImpParRecoFDFunc->GetParError(0);
    fFDMeanErr = ImpParRecoFDFunc->GetParError(1);
    fFDLambda1Err = ImpParRecoFDFunc->GetParError(2);
    fFDLambda2Err = ImpParRecoFDFunc->GetParError(3); 
    TCanvas *cRecoFD = new TCanvas("cRecoFD","cRecoFD",900,900);
    cRecoFD->SetLogy();
    hImpParRecoFD->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",fPtMin,fPtMax));
    hImpParRecoFD->GetXaxis()->SetTitle("Imp Par XY (#mum)");
    hImpParRecoFD->GetYaxis()->SetTitle(Form("Entries/(%0.f #mum)",hImpParRecoFD->GetBinWidth(5)));
    hImpParRecoFD->GetXaxis()->SetNdivisions(505);
    hImpParRecoFD->Draw("E1");
    cRecoFD->SaveAs(Form("ImpParRecoFD_%0.f-%0.f.pdf",fPtMin,fPtMax));
    cRecoFD->SaveAs(Form("ImpParRecoFD_%0.f-%0.f.root",fPtMin,fPtMax));
    delete cRecoFD;
    delete ImpParRecoFDFunc;
  }
  else {
    cerr << "only kConvolution and kGaussExpo are supported for the prefit on MC FD distribution" << endl;
    return;
  }  
}

//___________________________________________________________________________________
void AliDplusMassd0fitter::PrefitOnBkg(Double_t d0minBkg, Double_t d0maxBkg, Bool_t printSB)
{
  gStyle->SetOptFit(1);
  gStyle->SetTextSize(0.045);
  gStyle->SetTitleSize(0.05,"xy");
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.12);

  if(fBkgFromMC && !fMCBkgSparse) {
    cerr << "MC background THnSparse not set!!" << endl;
    return;
  }
  
  if(!fBkgFromMC && !fDataSparse) {
    cerr << "Data THnSparse not set!!" << endl;
    return;
  }
  
  TH1F* hImpParBkg=0x0;  
  if(fBkgFromMC) {
    ResetAxes(fMCBkgSparse);
    SetPtRange(fMCBkgSparse);
    SetMassRange(fMCBkgSparse);
    hImpParBkg=(TH1F*)fMCBkgSparse->Projection(fImpParAxis);  
  }
  else {
    hImpParBkg=(TH1F*)GetSidebandsDist(printSB);  
  }
  
  hImpParBkg->Rebin(fReb);
  hImpParBkg->SetMarkerStyle(21);
  hImpParBkg->SetMarkerSize(1.);
  hImpParBkg->SetLineWidth(2);
  
  TF1* ImpParBkgFunc=0x0;
  if(fBkgType==kSingleGaussExpo) {
    ImpParBkgFunc = new TF1("ImpParBkgFunc",this,&AliDplusMassd0fitter::FunctionImpParPrompt,d0minBkg,d0maxBkg,5,"AliDplusMassd0fitter","FunctionImpParPrompt");
    ImpParBkgFunc->SetLineColor(kMagenta);
    ImpParBkgFunc->SetParNames("fractionGauss","mean","sigma","lambda","integral");
    ImpParBkgFunc->SetParameters(fInitParamBkg);
    ImpParBkgFunc->SetParLimits(0,0.,1.);
    ImpParBkgFunc->SetParLimits(1,0.,10000000.);
    ImpParBkgFunc->SetParLimits(2,0.,10000000.);
    ImpParBkgFunc->SetParLimits(3,0.,10000000.);
    ImpParBkgFunc->FixParameter(4,hImpParBkg->Integral()*hImpParBkg->GetBinWidth(1));
  }
  else if(fBkgType==kDoubleGaussExpo){
    ImpParBkgFunc = new TF1("ImpParBkgFunc",this,&AliDplusMassd0fitter::FunctionImpParBkg,d0minBkg,d0maxBkg,10,"AliDplusMassd0fitter","FunctionImpParBkg");
    ImpParBkgFunc->SetLineColor(kMagenta);
    ImpParBkgFunc->SetParNames("fractionGauss1","mean1","sigma1","lambda1","fractionGauss2","mean2","sigma2","lambda2","fractionFunc1","integral");
    ImpParBkgFunc->SetParLimits(0,0.,1.);
    ImpParBkgFunc->SetParLimits(1,-10000000.,0.);
    ImpParBkgFunc->SetParLimits(2,0.,10000000.);
    ImpParBkgFunc->SetParLimits(3,0.,10000000.);
    ImpParBkgFunc->SetParLimits(4,0.,1.);
    ImpParBkgFunc->SetParLimits(5,0.,10000000.);
    ImpParBkgFunc->SetParLimits(6,0.,10000000.);
    ImpParBkgFunc->SetParLimits(7,0.,10000000.);
    ImpParBkgFunc->SetParLimits(8,0.,1.);
    ImpParBkgFunc->SetParameters(fInitParamBkg);
    ImpParBkgFunc->FixParameter(9,hImpParBkg->Integral()*hImpParBkg->GetBinWidth(1));
  }
  else if(fBkgType==kDoubleGaussExpoSymm){
    ImpParBkgFunc = new TF1("ImpParBkgFunc",this,&AliDplusMassd0fitter::FunctionImpParBkg,d0minBkg,d0maxBkg,6,"AliDplusMassd0fitter","FunctionImpParBkg");
    ImpParBkgFunc->SetLineColor(kMagenta);
    ImpParBkgFunc->SetParNames("fractionGauss1","mean1","sigma1","lambda1","mean2","integral");
    ImpParBkgFunc->SetParLimits(0,0.,1.);
    ImpParBkgFunc->SetParLimits(1,-10000000.,0.);
    ImpParBkgFunc->SetParLimits(2,0.,10000000.);
    ImpParBkgFunc->SetParLimits(3,0.,10000000.);
    ImpParBkgFunc->SetParLimits(4,0.,10000000.);
    ImpParBkgFunc->SetParameters(fInitParamBkg);
    ImpParBkgFunc->FixParameter(5,hImpParBkg->Integral()*hImpParBkg->GetBinWidth(1));
  }
  else {
    cerr << "only kSingleGaussExpo, kDoubleGaussExpo and kDoubleGaussExpoSymm are supported for the prefit on bkg distribution" << endl;
    return;
  }
  
  TCanvas *cBkg = new TCanvas("cBkg","cBkg",900,900);
  cBkg->SetLogy();
  hImpParBkg->SetLineColor(kBlack);
  hImpParBkg->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",fPtMin,fPtMax));
  hImpParBkg->GetXaxis()->SetTitle("Imp Par XY (#mum)");
  hImpParBkg->GetYaxis()->SetTitle(Form("Entries/(%0.f #mum)",hImpParBkg->GetBinWidth(5)));
  hImpParBkg->GetXaxis()->SetNdivisions(505);
  hImpParBkg->Draw("E1");
  TVirtualFitter::SetDefaultFitter("Minuit2");
  hImpParBkg->Fit("ImpParBkgFunc",fFitOptionsPrefit.Data());
  fBkgFractionGauss1 = ImpParBkgFunc->GetParameter(0);
  fBkgMean1 = ImpParBkgFunc->GetParameter(1);
  fBkgSigma1 = ImpParBkgFunc->GetParameter(2);
  fBkgLambda1 = ImpParBkgFunc->GetParameter(3);
  fBkgFractionGauss1Err = ImpParBkgFunc->GetParError(0);
  fBkgMean1Err = ImpParBkgFunc->GetParError(1);
  fBkgSigma1Err = ImpParBkgFunc->GetParError(2);
  fBkgLambda1Err = ImpParBkgFunc->GetParError(3);
  if(fBkgType==kDoubleGaussExpo) { 
    fBkgFractionGauss2 = ImpParBkgFunc->GetParameter(4);
    fBkgMean2 = ImpParBkgFunc->GetParameter(5);
    fBkgSigma2 = ImpParBkgFunc->GetParameter(6);
    fBkgLambda2 = ImpParBkgFunc->GetParameter(7);
    fBkgFracFunc1 = ImpParBkgFunc->GetParameter(8);
    fBkgFractionGauss2Err = ImpParBkgFunc->GetParError(4);
    fBkgMean2Err = ImpParBkgFunc->GetParError(5);
    fBkgSigma2Err = ImpParBkgFunc->GetParError(6);
    fBkgLambda2Err = ImpParBkgFunc->GetParError(7);
    fBkgFracFunc1Err = ImpParBkgFunc->GetParError(8);
  }
  if(fBkgType==kDoubleGaussExpoSymm) {
    fBkgFractionGauss2 = ImpParBkgFunc->GetParameter(0);
    fBkgMean2 = ImpParBkgFunc->GetParameter(4);
    fBkgSigma2 = ImpParBkgFunc->GetParameter(2);
    fBkgLambda2 = ImpParBkgFunc->GetParameter(3);
    fBkgFractionGauss2Err = ImpParBkgFunc->GetParError(0);
    fBkgMean2Err = ImpParBkgFunc->GetParError(4);
    fBkgSigma2Err = ImpParBkgFunc->GetParError(2);
    fBkgLambda2Err = ImpParBkgFunc->GetParError(3);
  }
  
  cBkg->SaveAs(Form("ImpParBkg_%0.f-%0.f.pdf",fPtMin,fPtMax));
  cBkg->SaveAs(Form("ImpParBkg_%0.f-%0.f.root",fPtMin,fPtMax));
  delete cBkg;
  delete ImpParBkgFunc;
}

//_______________________________________________________________________________
void AliDplusMassd0fitter::PrefitStep(Double_t d0minPrompt, Double_t d0maxPrompt, Double_t d0minFD, Double_t d0maxFD, Double_t d0minBkg, Double_t d0maxBkg)
{
  PrefitOnPrompt(d0minPrompt,d0maxPrompt);
  PrefitOnFD(d0minFD,d0maxFD);
  PrefitOnBkg(d0minBkg,d0maxBkg);
}

//_______________________________________________________________________________
void AliDplusMassd0fitter::GenerateEntries(Bool_t background, Int_t gentype, Int_t treeOrhisto)
{
  if(!fMCPromptSparse) {
    cerr << "MC Prompt THnSparse not set!!" << endl;
    return;
  }

  if(!fMCTrueFDSparse) {
    cerr << "MC true FD THnSparse not set!!" << endl;
    return;
  }

  if(!fMCRecoFDSparse) {
    cerr << "MC reco FD THnSparse not set!!" << endl;
    return;
  }

  if(fBkgFromMC && !fMCBkgSparse) {
    cerr << "MC background THnSparse not set!!" << endl;
    return;
  }
  
  if(!fBkgFromMC && !fDataSparse) {
    cerr << "Data THnSparse not set!!" << endl;
    return;
  }

  fMCTest = kTRUE; 
  
  TH1F* hImpParPrompt=0x0;
  TH1F* hImpParRecoFD=0x0;
  TH1F* hImpParBkg=0x0;
  TF1* PromptFunc=0x0;
  TF1* RecoFDFunc=0x0;
  TF1* BkgFunc=0x0;
  
  if(gentype==kFromFunc) {
    if(fPromptSigmaMC==0){
      PrefitStep();
    }
    PromptFunc = new TF1("PromptFunc",this,&AliDplusMassd0fitter::FunctionImpParPrompt,-1000,1000,5,"AliDplusMassd0fitter","FunctionImpParPrompt");
    PromptFunc->SetParameters(fPromptFractionGauss,fPromptMean,fPromptSigmaMC,fPromptLambda,1.);
    if(fFDType==kConvolution) {
      RecoFDFunc = new TF1("RecoFDFunc",this,&AliDplusMassd0fitter::FunctionRecoImpParFD,-1000,1000,1,"AliDplusMassd0fitter","FunctionRecoImpParFD");
      RecoFDFunc->SetParameter(0,1.);
    }
    else if(fFDType==kGaussExpo) {
      RecoFDFunc = new TF1("RecoFDFunc",this,&AliDplusMassd0fitter::FunctionImpParPrompt,-1000,1000,5,"AliDplusMassd0fitter","FunctionRecoImpParPrompt");
      RecoFDFunc->SetParameters(fFDFraction1,fFDMean,fFDLambda1,fFDLambda2,1.);
    }
    else {
      cerr << "only kConvolution and kGaussExpo are supported for the prefit on MC FD distribution" << endl;
      return;
    }

    if(fBkgType==kSingleGaussExpo) {
      BkgFunc = new TF1("BkgFunc",this,&AliDplusMassd0fitter::FunctionImpParPrompt,-1000,1000,5,"AliDplusMassd0fitter","FunctionImpParPrompt");
      BkgFunc->SetParameters(fBkgFractionGauss1,fBkgMean1,fBkgSigma1,fBkgLambda1,1.);
    }
    else if(fBkgType==kDoubleGaussExpo || fBkgType==kDoubleGaussExpoSymm) {
      BkgFunc = new TF1("BkgFunc",this,&AliDplusMassd0fitter::FunctionImpParBkg,-1000,1000,9,"AliDplusMassd0fitter","FunctionImpParBkg");
      BkgFunc->SetParameters(fBkgFractionGauss1,fBkgMean1,fBkgSigma1,fBkgLambda1,fBkgFractionGauss2,fBkgMean2,fBkgSigma2,fBkgLambda2,1.);
    }
    else{
      cerr << "only kSingleGaussExpo, kDoubleGaussExpo and kDoubleGaussExpoSymm are supported for the prefit on SB distribution" << endl;
      return;
    }
  }
  else if(gentype==kFromHisto) {
    ResetAxes(fMCPromptSparse);
    ResetAxes(fMCRecoFDSparse);
    SetPtRange(fMCPromptSparse);
    SetPtRange(fMCRecoFDSparse);
    SetMassRange(fMCPromptSparse);
    SetMassRange(fMCRecoFDSparse);
    
    hImpParPrompt=(TH1F*)fMCPromptSparse->Projection(fImpParAxis);
    hImpParRecoFD=(TH1F*)fMCRecoFDSparse->Projection(fImpParAxis);  
    
    if(fBkgFromMC) {
      ResetAxes(fMCBkgSparse);
      SetPtRange(fMCBkgSparse);
      SetMassRange(fMCBkgSparse);
      hImpParBkg=(TH1F*)fMCBkgSparse->Projection(fImpParAxis);  
    }
    else {
      hImpParBkg=(TH1F*)GetSidebandsDist();
    }
  }
  
  Double_t RndImpPar;
  
  if(treeOrhisto==kTree) {
    fDataTree = new TTree(Form("ImpParTree_%0.f-%0.f",fPtMin,fPtMax),"impact parameter tree");
    fDataTree->Branch(fImpParBranch.Data(), &RndImpPar);
  }
  else if(treeOrhisto==kHisto) {
    if(fImpParHisto)
      fImpParHisto->Clear();
    fImpParHisto = new TH1F("fImpParHisto","",400,-1000,1000);
  }
  else {
    cerr << "Only tree or histo can be generated for the TOY MC!!!" << endl;
    return;
  }
    
  for(Int_t iGen=0; iGen<fGenPromptFraction*fSig; iGen++) {
    if(gentype==kFromFunc)
      RndImpPar = PromptFunc->GetRandom();
    else if(gentype==kFromHisto)
      RndImpPar = hImpParPrompt->GetRandom();
    else
      cout << "Generation of the entries can only be done from functions (0) or histos (1)" << endl;
    if(treeOrhisto==kTree) {
      fDataTree->Fill();
    }
    else {
      fImpParHisto->Fill(RndImpPar);
    }
  }

  for(Int_t iGen=0; iGen<(1-fGenPromptFraction)*fSig; iGen++) {
    if(gentype==kFromFunc)
      RndImpPar = RecoFDFunc->GetRandom();
    else if(gentype==kFromHisto)
      RndImpPar = hImpParRecoFD->GetRandom();
    else
      cout << "Generation of the entries can only be done from functions (0) or histos (1)" << endl;
    if(treeOrhisto==kTree) {
      fDataTree->Fill();
    }
    else {
      fImpParHisto->Fill(RndImpPar);
    }
  }
  
  if(!background)
    fBkg=0;
  
  for(Int_t iGen=0; iGen<fBkg; iGen++) {
    if(gentype==kFromFunc) {
      RndImpPar = BkgFunc->GetRandom();
    }
    else if(gentype==kFromHisto)
      RndImpPar = hImpParBkg->GetRandom();
    else
      cout << "Generation of the entries can only be done from functions (0) or histos (1)" <<endl;
    if(treeOrhisto==kTree) {
      fDataTree->Fill();
    }
    else {
      fImpParHisto->Fill(RndImpPar);
    }
  }
  
  delete hImpParPrompt;
  delete hImpParRecoFD;
  delete hImpParBkg;
  delete PromptFunc;
  delete RecoFDFunc;
  delete BkgFunc;
}

//_________________________________________________________________________________
void AliDplusMassd0fitter::FitTree(Double_t d0min, Double_t d0max, Bool_t print)
{
  if(fSig==0)
    GetSignal();
  if(fPromptSigmaMC==0)
    PrefitStep();

  if(!fDataTree) {
    cerr << "Data Tree not set!!" << endl;
    return;
  }

  if(fDataTree->GetBranch(fPtBranch.Data()) && fDataTree->GetBranch(fMassBranch.Data())) { ///real data tree
    Double_t massmin=fMassMean-fNSigma*fMassSigma;
    Double_t massmax=fMassMean+fNSigma*fMassSigma;
    fIntegral=fDataTree->GetEntries(Form("%s>%f && %s<%f && %s>%f && %s<%f",fPtBranch.Data(),fPtMin,fPtBranch.Data(),fPtMax,fMassBranch.Data(),massmin,fMassBranch.Data(),massmax));

  }
  else { ///tree generated with TOY MC
    fIntegral=fDataTree->GetEntries();
  }
  fBkg=fIntegral-fSig;
  fBkgErrStat=fSigErrStat;
  
  TF1* ImpParFitFunction = new TF1("ImpParFitFunction",this,&AliDplusMassd0fitter::FitFunction,d0min,d0max,2,"ImpParFitFunction","FitFunction");
  Int_t nfreepars=2;
  
  ImpParFitFunction->SetParNames("PromptFraction","SigmaPrompt");
  ImpParFitFunction->SetParameter(0,fGenPromptFraction);
  ImpParFitFunction->SetParLimits(0,fPromptFracMin,fPromptFracMax);
  if(fSigmaPromptFixed) {
    ImpParFitFunction->FixParameter(1,fPromptSigmaMC);
    nfreepars=1;
  }
  else {
    ImpParFitFunction->SetParameter(1,fPromptSigmaMC);
    Double_t sigmamin = fPromptSigmaMC-fSigmaPromptLim*fPromptSigmaMC;
    Double_t sigmamax = fPromptSigmaMC+fSigmaPromptLim*fPromptSigmaMC;
    ImpParFitFunction->SetParLimits(1,sigmamin,sigmamax);
  }

  if(fDataTree->GetBranch("Pt") && fDataTree->GetBranch("InvMass")) {
    Double_t massmin=fMassMean-fNSigma*fMassSigma;
    Double_t massmax=fMassMean+fNSigma*fMassSigma;
    fDataTree->UnbinnedFit("ImpParFitFunction",fImpParBranch.Data(),Form("%s>%f && %s<%f && %s>%f && %s<%f",fPtBranch.Data(),fPtMin,fPtBranch.Data(),fPtMax,fMassBranch.Data(),massmin,fMassBranch.Data(),massmax),fFitOptions.Data());
  }
  else {
    fDataTree->UnbinnedFit("ImpParFitFunction",fImpParBranch.Data(),"",fFitOptions.Data());
  }

  fPromptFraction = ImpParFitFunction->GetParameter(0);
  fPromptFractionErr = ImpParFitFunction->GetParError(0);
  fPromptSigma = ImpParFitFunction->GetParameter(1);
  fPromptSigmaErr = ImpParFitFunction->GetParError(1);
  
  if(nfreepars==2) {
    TVirtualFitter *fitter = TVirtualFitter::GetFitter();
    TMatrixD matrix(2,2,fitter->GetCovarianceMatrix());
    fCovFracSigmaPrompt = fitter->GetCovarianceMatrixElement(0,1);
  }
  else {
    fCovFracSigmaPrompt=0;
  }
  
  ///chi square -> as it was a binned fit (not very exact)
  Double_t massmin=fMassMean-fNSigma*fMassSigma;
  Double_t massmax=fMassMean+fNSigma*fMassSigma;
  fImpParHisto = new TH1F("fImpParHisto","",400,-1000,1000);
  if(fDataTree->GetBranch("Pt") && fDataTree->GetBranch("InvMass")) {
   fDataTree->Project(Form("%s>>fImpParHisto",fImpParBranch.Data()),Form("%s>%f && %s<%f && %s>%f && %s<%f",fPtBranch.Data(),fPtMin,fPtBranch.Data(),fPtMax,fMassBranch.Data(),massmin,fMassBranch.Data(),massmax));
   }
   else {
     fDataTree->Project("fImpParHisto",Form("%s",fImpParBranch.Data()));
   }
  fChiSquare = fImpParHisto->Chisquare(ImpParFitFunction);
  fNDF = 0;
  for(Int_t iBin=0; iBin<fImpParHisto->GetNbinsX(); iBin++) {
    if(fImpParHisto->GetBinContent(iBin+1)!=0) fNDF++;
  }
  fNDF -= nfreepars;
  
  if(print)
    DrawResult(kTRUE);
  
  delete ImpParFitFunction;
}

//___________________________________________________________________________
void AliDplusMassd0fitter::FitHisto(Double_t d0min, Double_t d0max, Bool_t print)
{
  if(fSig==0)
    GetSignal();
  if(fPromptSigmaMC==0)
    PrefitStep();
  
  if(!fDataSparse) {
    cerr << "Data THnSparse not set!!" << endl;
    return;
  }

  Double_t nbkg_temp=0;

  if(!fMCTest) {
    ResetAxes(fDataSparse);
    SetPtRange(fDataSparse);
    SetMassRange(fDataSparse);
    fImpParHisto = (TH1F*)fDataSparse->Projection(fImpParAxis);
  }
  
  fIntegral = fImpParHisto->GetEntries();
  fSig *= fImpParHisto->GetBinWidth(1); 
  fIntegral *= fImpParHisto->GetBinWidth(1);
  fBkg = fIntegral-fSig;
  fBkgErrStat=fSigErrStat;

  TF1* ImpParFitFunction = new TF1("ImpParFitFunction",this,&AliDplusMassd0fitter::FitFunction,d0min,d0max,2,"ImpParFitFunction","FitFunction");
  
  Int_t nfreepars=2;
  ImpParFitFunction->SetParNames("PromptFraction","SigmaPrompt");
  ImpParFitFunction->SetParameter(0,fGenPromptFraction);
  ImpParFitFunction->SetParLimits(0,fPromptFracMin,fPromptFracMax);
  if(fSigmaPromptFixed) {
    ImpParFitFunction->FixParameter(1,fPromptSigmaMC);
    nfreepars=1;
  }
  else {
    ImpParFitFunction->SetParameter(1,fPromptSigmaMC);
    Double_t sigmamin = fPromptSigmaMC-fSigmaPromptLim*fPromptSigmaMC;
    Double_t sigmamax = fPromptSigmaMC+fSigmaPromptLim*fPromptSigmaMC;
    ImpParFitFunction->SetParLimits(1,sigmamin,sigmamax);
  }

  if(fSubBkg && !fVariableBinning) {
    TH1F* hImpParSideBands=(TH1F*)GetSidebandsDist();

    Double_t bkgintegral = fImpParHisto->Integral()*fBkg/fIntegral;    
    hImpParSideBands->Sumw2();
    hImpParSideBands->Scale(bkgintegral/hImpParSideBands->Integral());
    hImpParSideBands->SetDirectory(0);

    fImpParHisto->Sumw2();
    TH1F* hImpParSub = (TH1F*)fImpParHisto->Clone();
    hImpParSub->SetDirectory(0);
    hImpParSub->Sumw2();
    
    fImpParHisto->Add(hImpParSub,hImpParSideBands,1.,-1.);

    //rebin -> normalisation must take into account the binning
    fSig /= fImpParHisto->GetBinWidth(1); 
    fBkg /= fImpParHisto->GetBinWidth(1);
    fIntegral /= fImpParHisto->GetBinWidth(1);

    fImpParHisto->Rebin(fReb);
    fBkg *= fImpParHisto->GetBinWidth(1);
    fIntegral = fSig;
    nbkg_temp = fBkg;
    fBkg = 0;
    fSig *= fImpParHisto->GetBinWidth(1); 
    fIntegral *= fImpParHisto->GetBinWidth(1);
  }
  else if(fVariableBinning){
    //rebin -> normalisation must take into account the binning
    fSig /= fImpParHisto->GetBinWidth(1); 
    fBkg /= fImpParHisto->GetBinWidth(1);
    fIntegral /= fImpParHisto->GetBinWidth(1);
    
    RebinVariableWidth();
    fImpParHisto->Sumw2();
    fImpParHisto->Scale(1.,"width"); ///rescales each bin for its bin width
    fFitOptions.ReplaceAll("M","I"); ///if variable binning in the fit options MUST be "WLI"
    if(!fFitOptions.Contains("I")) fFitOptions += "I"; ///if variable binning in the fit options MUST be "I" 
  }
  else {
    //rebin -> normalisation must take into account the binning
    fSig /= fImpParHisto->GetBinWidth(1); 
    fBkg /= fImpParHisto->GetBinWidth(1);
    fIntegral /= fImpParHisto->GetBinWidth(1);

    fImpParHisto->Rebin(fReb);
    fBkg *= fImpParHisto->GetBinWidth(1);
    fSig *= fImpParHisto->GetBinWidth(1); 
    fIntegral=fIntegral*fImpParHisto->GetBinWidth(1);
  }

  if(!fVariableBinning)
    TVirtualFitter::SetDefaultFitter("Minuit2");
  fImpParHisto->Fit("ImpParFitFunction",fFitOptions.Data());
  
  fPromptFraction = ImpParFitFunction->GetParameter(0);
  fPromptFractionErr = ImpParFitFunction->GetParError(0);
  fPromptSigma = ImpParFitFunction->GetParameter(1);
  fPromptSigmaErr = ImpParFitFunction->GetParError(1);

  if(nfreepars==2) {
    TVirtualFitter *fitter = TVirtualFitter::GetFitter();
    TMatrixD matrix(2,2,fitter->GetCovarianceMatrix());
    fCovFracSigmaPrompt = fitter->GetCovarianceMatrixElement(0,1);
  }
  else {
    fCovFracSigmaPrompt=0;
  }
  
  fNDF = ImpParFitFunction->GetNDF();
  fChiSquare = ImpParFitFunction->GetChisquare();
  
  if(print)
    DrawResult(kFALSE);

  if(fSubBkg)
    fBkg = nbkg_temp;

  if(!fVariableBinning) {
    fIntegral /= fImpParHisto->GetBinWidth(1);
    fBkg /= fImpParHisto->GetBinWidth(1);
    fSig /= fImpParHisto->GetBinWidth(1);
  }
  
  delete ImpParFitFunction;  
}

//_________________________________________________________________________
void AliDplusMassd0fitter::SetPID(Bool_t isPIDon)
{
  if(fDataSparse->GetNdimensions() >= fPIDAxis) {
    if(isPIDon) {
//      if(fDataSparse) fDataSparse->GetAxis(fPIDAxis)->SetRange(fPIDbin,fPIDbin);
      if(fMCPromptSparse) fMCPromptSparse->GetAxis(fPIDAxis)->SetRange(fPIDbin,fPIDbin);
      if(fMCTrueFDSparse) fMCTrueFDSparse->GetAxis(fPIDAxis)->SetRange(fPIDbin,fPIDbin);
      if(fMCRecoFDSparse) fMCRecoFDSparse->GetAxis(fPIDAxis)->SetRange(fPIDbin,fPIDbin);
      if(fMCPromptSparse) fMCPromptSparse->GetAxis(11)->SetRange(16.,25.);
      if(fMCTrueFDSparse) fMCTrueFDSparse->GetAxis(11)->SetRange(16.,25.);
      if(fMCRecoFDSparse) fMCRecoFDSparse->GetAxis(11)->SetRange(16.,25.);
      if(fMCBkgSparse) fMCBkgSparse->GetAxis(fPIDAxis)->SetRange(fPIDbin,fPIDbin);
    }
    else {
      //if(fDataSparse) fDataSparse->GetAxis(fPIDAxis)->SetRange(-1,-1);
      if(fMCPromptSparse) fMCPromptSparse->GetAxis(fPIDAxis)->SetRange(-1,-1);
      if(fMCTrueFDSparse) fMCTrueFDSparse->GetAxis(fPIDAxis)->SetRange(-1,-1);
      if(fMCRecoFDSparse) fMCRecoFDSparse->GetAxis(fPIDAxis)->SetRange(-1,-1);
      if(fMCBkgSparse) fMCBkgSparse->GetAxis(fPIDAxis)->SetRange(-1,-1);
    }
  }
  else {
    cerr << "PID axis = " << fPIDAxis << ", number of axes in the THnSparses = " << fDataSparse->GetNdimensions() << endl;
    cerr << "Apply the PID cut is not possible!" << endl;
  }
}

//_________________________________________________________________________________________________
void AliDplusMassd0fitter::GetPromptFractionWithIPCut(Double_t d0cut, Double_t &promptfracd0cut, Double_t &err, Bool_t genfrac)
{
  TF1* PromptFunc = new TF1("PromptFunc",this,&AliDplusMassd0fitter::FunctionImpParPrompt,-1000,1000,5,"AliDplusMassd0fitter","FunctionImpParPrompt");
  PromptFunc->SetParameters(fPromptFractionGauss,fPromptMean,fPromptSigma,fPromptLambda,1.);
  
  TF1* RecoFDFunc = 0x0;
  if(fFDType==kConvolution) {
    RecoFDFunc = new TF1("RecoFDFunc",this,&AliDplusMassd0fitter::FunctionRecoImpParFD,-1000,1000,1,"AliDplusMassd0fitter","FunctionRecoImpParFD");
    RecoFDFunc->SetParameter(0,1.);
  }
  else if(fFDType==kGaussExpo) {
    RecoFDFunc = new TF1("RecoFDFunc",this,&AliDplusMassd0fitter::FunctionImpParPrompt,-1000,1000,5,"AliDplusMassd0fitter","FunctionRecoImpParPrompt");
    RecoFDFunc->SetParameters(fFDFraction1,fFDMean,fFDLambda1,fFDLambda2,1.);
  }
  else {
    cerr << "only kConvolution and kGaussExpo are supported for the prefit on MC FD distribution" << endl;
    return;
  }
  
  Double_t PromptInt = PromptFunc->Integral(-d0cut,d0cut);
  Double_t RecoFDInt = RecoFDFunc->Integral(-d0cut,d0cut);
  
  delete PromptFunc;
  delete RecoFDFunc;
  
  if(!genfrac) {
    promptfracd0cut = PromptInt*fPromptFraction/(PromptInt*fPromptFraction+RecoFDInt*(1-fPromptFraction));
    
    //derivatives for error propagation
    Double_t dfdfprompt = (PromptInt*(fPromptFraction*(PromptInt-RecoFDInt)+RecoFDInt)-fPromptFraction*PromptInt*(PromptInt-RecoFDInt))/((PromptInt*fPromptFraction+RecoFDInt*(1-fPromptFraction))*(PromptInt*fPromptFraction+RecoFDInt*(1-fPromptFraction)));
    
    TF1* PropmtFuncSigmaDer = new TF1("PromptFuncSigmaDer",this,&AliDplusMassd0fitter::GaussSigmaDerivative,-1000,1000,3,"AliDplusMassd0fitter","GaussSigmaDerivative");
    PropmtFuncSigmaDer->SetParameters(fPromptFractionGauss,fPromptMean,fPromptSigma);
    Double_t PromptSigmaDerInt=PropmtFuncSigmaDer->Integral(-d0cut,d0cut);
    delete PropmtFuncSigmaDer;
    
    Double_t RecoFDSigmaDerInt=0;
    if(fFDType==kConvolution) {
      TF1* RecoFDFuncSigmaDer = new TF1("RecoFDFuncSigmaDer",this,&AliDplusMassd0fitter::RecoFDSigmaDerivative,-1000,1000,1,"AliDplusMassd0fitter","RecoFDSigmaDerivative");
      RecoFDFuncSigmaDer->SetParameter(0,1.);
      RecoFDSigmaDerInt=RecoFDFuncSigmaDer->Integral(-d0cut,d0cut);
      delete RecoFDFuncSigmaDer;
    }
    
    Double_t dfdsigmaprompt = (fPromptFraction*PromptSigmaDerInt*(PromptInt*fPromptFraction+RecoFDInt*(1-fPromptFraction))-fPromptFraction*PromptInt*(fPromptFraction*PromptSigmaDerInt+(1-fPromptFraction)*RecoFDSigmaDerInt))/((PromptInt*fPromptFraction+RecoFDInt*(1-fPromptFraction))*(PromptInt*fPromptFraction+RecoFDInt*(1-fPromptFraction)));
    
    //error on fprompt within [-d0cut,d0cut]
    err = TMath::Sqrt(dfdfprompt*fPromptFractionErr*dfdfprompt*fPromptFractionErr+
                      dfdsigmaprompt*fPromptSigmaErr*dfdsigmaprompt*fPromptSigmaErr+
                      2*dfdfprompt*dfdsigmaprompt*fCovFracSigmaPrompt);
  }
  else {
    promptfracd0cut = PromptInt*fGenPromptFraction/(PromptInt*fGenPromptFraction+RecoFDInt*(1-fGenPromptFraction));
    err = 0;
  }
  
}

//______________________________________________________________________________
Double_t AliDplusMassd0fitter::Gauss(Double_t d0,Double_t mean,Double_t sigma)
{
  return TMath::Gaus(d0,mean,sigma,kTRUE); //kTRUE means that is normalised
}

//______________________________________________________________________________
Double_t AliDplusMassd0fitter::ExpoDouble(Double_t d0, Double_t mean, Double_t lambda)
{
  return 1./(2.*lambda)*TMath::Exp(-TMath::Abs(d0-mean)/lambda);
}

//______________________________________________________________________________
Double_t AliDplusMassd0fitter::FunctionImpParPrompt(Double_t* x, Double_t* par)
{
  Double_t d0 = x[0];
  Double_t fractionG = par[0];
  Double_t meanD = par[1];
  Double_t sigma = par[2];
  Double_t lambdaD = par[3];
  Double_t norm = par[4];

  return norm*((1.-fractionG)*ExpoDouble(d0,meanD,lambdaD)+fractionG*Gauss(d0,meanD,sigma));
}

//______________________________________________________________________________
Double_t AliDplusMassd0fitter::FunctionTrueImpParFD(Double_t* x, Double_t* par)
{    
  Double_t d0 = x[0];
  Double_t fraction1 = par[0];
  Double_t meanB = par[1];
  Double_t lambda1 = par[2];
  Double_t lambda2 = par[3];
  Double_t norm = par[4];
  
  return norm*(fraction1*ExpoDouble(d0,meanB,lambda1)+(1.-fraction1)*ExpoDouble(d0,meanB,lambda2));
}

//______________________________________________________________________________
Double_t AliDplusMassd0fitter::FunctionRecoImpParFD(Double_t *x, Double_t *par)
{
  Double_t norm = par[0];
  
  return norm*Convolution(x[0],-1000,1000,1000); 
}

//______________________________________________________________________________
Double_t AliDplusMassd0fitter::FunctionImpParBkg(Double_t* x, Double_t* par)
{
  Double_t fracfunc1=0.5;
  Double_t fractiongauss1 = par[0];
  Double_t mean1 = par[1];
  Double_t sigma1 = par[2];
  Double_t lambda1 = par[3];
  Double_t fractiongauss2=0;
  Double_t mean2=0;
  Double_t sigma2=0;
  Double_t lambda2=0;
  Double_t norm=0;

  Double_t pars1[5] = {fractiongauss1,mean1,sigma1,lambda1,1.};
  Double_t pars2[5];

  if(fBkgType==kDoubleGaussExpoSymm) {
    mean2=par[4];
    norm=par[5];
    pars2[0]=fractiongauss1;
    pars2[1]=mean2;
    pars2[2]=sigma1;
    pars2[3]=lambda1;
  }
  else {
    fractiongauss2=par[4];
    mean2=par[5];
    sigma2=par[6];
    lambda2=par[7];
    fracfunc1=par[8];
    norm=par[9];
    pars2[0]=fractiongauss2;
    pars2[1]=mean2;
    pars2[2]=sigma2;
    pars2[3]=lambda2;
  }
  pars2[4]=1.;
  
  return norm*(fracfunc1*FunctionImpParPrompt(x,pars1)+(1-fracfunc1)*FunctionImpParPrompt(x,pars2));
}

//_________________________________________________________________________
Double_t AliDplusMassd0fitter::Convolution(Double_t x, Double_t xmin, Double_t xmax, Int_t nstpdf)
{
  Double_t parprompt[5] = {fPromptFractionGauss,fPromptMean,fPromptSigma,fPromptLambda,1.};
  Double_t parFD[5] = {fFDFraction1,fFDMean,fFDLambda1,fFDLambda2,1.};
  
  Double_t d0true[1];
  Double_t diffd0[1];

  Double_t sum=0.;
  Double_t dx=(xmax-xmin)/(Double_t)nstpdf;
    for(Int_t j=0;j<nstpdf;j++){
      d0true[0]=xmin+dx*j;
      diffd0[0]=x-d0true[0];
      Double_t molt = FunctionTrueImpParFD(d0true,parFD)*FunctionImpParPrompt(diffd0,parprompt);
      sum = sum+molt*dx;
    }
    
    return sum;
}

//_________________________________________________________________________
Double_t AliDplusMassd0fitter::FitFunction(Double_t* x, Double_t *par)
{
  Double_t fractionPrompt = par[0];
  fPromptSigma = par[1];
  Double_t promptlambda = fPromptLambda/fPromptSigmaMC*par[1];

  Double_t parprompt[5] = {fPromptFractionGauss,fPromptMean,fPromptSigma,promptlambda,1.};
  Double_t parFD[5];
  Double_t parBkg[10];
  parBkg[0]=fBkgFractionGauss1;
  parBkg[1]=fBkgMean1;
  parBkg[2]=fBkgSigma1;
  parBkg[3]=fBkgLambda1;
  
  if(fBkgType==kSingleGaussExpo && fFDType==kConvolution) {
    parFD[0]=1.;
    parBkg[4]=1.;

    return fIntegral*((fSig/fIntegral)*(fractionPrompt*FunctionImpParPrompt(x,parprompt)+(1-fractionPrompt)*FunctionRecoImpParFD(x,parFD))+(1-fSig/fIntegral)*FunctionImpParPrompt(x,parBkg));
  }
  else if(fBkgType==kDoubleGaussExpo && fFDType==kConvolution) {
    parFD[0]=1.;
    parBkg[4]=fBkgFractionGauss2;
    parBkg[5]=fBkgMean2;
    parBkg[6]=fBkgSigma2;
    parBkg[7]=fBkgLambda2;
    parBkg[8]=fBkgFracFunc1;
    parBkg[9]=1.;

    return fIntegral*((fSig/fIntegral)*(fractionPrompt*FunctionImpParPrompt(x,parprompt)+(1-fractionPrompt)*FunctionRecoImpParFD(x,parFD))+(1-fSig/fIntegral)*FunctionImpParBkg(x,parBkg));
  }
  else if(fBkgType==kDoubleGaussExpoSymm && fFDType==kConvolution) {
    parFD[0]=1.;
    parBkg[4]=fBkgMean2;
    parBkg[5]=1.;

    return fIntegral*((fSig/fIntegral)*(fractionPrompt*FunctionImpParPrompt(x,parprompt)+(1-fractionPrompt)*FunctionRecoImpParFD(x,parFD))+(1-fSig/fIntegral)*FunctionImpParBkg(x,parBkg));
  }
  else {
    parFD[0]=fFDFraction1;
    parFD[1]=fFDMean;
    parFD[2]=fFDLambda1;
    parFD[3]=fFDLambda2;
    parFD[4]=1.;
    parBkg[4]=fBkgFractionGauss2;
    parBkg[5]=fBkgMean2;
    parBkg[6]=fBkgSigma2;
    parBkg[7]=fBkgLambda2;
    parBkg[8]=fBkgFracFunc1;
    parBkg[9]=1.;

    return fIntegral*((fSig/fIntegral)*(fractionPrompt*FunctionImpParPrompt(x,parprompt)+(1-fractionPrompt)*FunctionImpParPrompt(x,parFD))+(1-fSig/fIntegral)*FunctionImpParBkg(x,parBkg));  
  }
}

//_________________________________________________________________________________________________
Double_t AliDplusMassd0fitter::GaussSigmaDerivative(Double_t *d0,Double_t *pars)
{
  Double_t d00 = d0[0];
  Double_t fgaus = pars[0];
  Double_t mean = pars[1];
  Double_t sigma = pars[2];
  
  return fgaus*(-1/(2*sigma)+(d00-mean)*(d00-mean)/(sigma*sigma*sigma))*TMath::Gaus(d00,mean,sigma,kTRUE);
}

//_________________________________________________________________________
Double_t AliDplusMassd0fitter::ConvolutionSigmaDerivative(Double_t x, Double_t xmin, Double_t xmax, Int_t nstpdf)
{
  Double_t pargausder[3] = {fPromptFractionGauss,fPromptMean,fPromptSigma};
  Double_t parFD[5] = {fFDFraction1,fFDMean,fFDLambda1,fFDLambda2,1.};
  
  Double_t d0true[1];
  Double_t diffd0[1];
  
  Double_t sum=0.;
  Double_t dx=(xmax-xmin)/(Double_t)nstpdf;
  for(Int_t j=0;j<nstpdf;j++){
    d0true[0]=xmin+dx*j;
    diffd0[0]=x-d0true[0];
    Double_t molt = FunctionTrueImpParFD(d0true,parFD)*GaussSigmaDerivative(diffd0,pargausder);
    sum = sum+molt*dx;
  }
  
  return sum;
}

//_________________________________________________________________________
Double_t AliDplusMassd0fitter::RecoFDSigmaDerivative(Double_t *d0, Double_t *pars)
{
  Double_t d00=d0[0];
  Double_t norm=pars[0];
  
  return norm*ConvolutionSigmaDerivative(d00,-1000,1000,1000);
}

//_________________________________________________________________________
TH1F* AliDplusMassd0fitter::GetPtreweightedHisto(THnSparseF* sparse, Int_t meson)
{
  ResetAxes(sparse);
  
  TH2F* hPtImpPar = (TH2F*)sparse->Projection(fPtAxis,fImpParAxis);
  hPtImpPar->SetDirectory(0);

  TAxis* ptaxisweights = 0x0;
  Double_t binwidthweights = 0;
  
  if(meson==kB) {
    ptaxisweights = (TAxis*)fPtBWeightsHisto->GetXaxis();
    binwidthweights = fPtBWeightsHisto->GetBinWidth(1);
  }
  else {
    ptaxisweights = (TAxis*)fPtDWeightsHisto->GetXaxis();
    binwidthweights = fPtDWeightsHisto->GetBinWidth(1);
  }  
  
  TAxis* ptaxis = (TAxis*)hPtImpPar->GetYaxis();
  Double_t ptbinwidth = ptaxis->GetBinWidth(1);
  
  TAxis* impparaxis = (TAxis*)hPtImpPar->GetXaxis();
  const Int_t nImpParBins = impparaxis->GetNbins();
  Double_t ImpParMin = impparaxis->GetBinLowEdge(1);
  Double_t ImpParMax = impparaxis->GetBinLowEdge(nImpParBins)+impparaxis->GetBinWidth(nImpParBins);  
  
  Double_t binwidthratio = ptbinwidth/binwidthweights;
  if(binwidthratio>1) {
    Int_t rebin = (Int_t)binwidthratio;
    if(meson==kB) {
      fPtBWeightsHisto->Rebin(rebin);
      ptaxisweights = (TAxis*)fPtBWeightsHisto->GetXaxis();
    }
    else {
      fPtDWeightsHisto->Rebin(rebin);
      ptaxisweights = (TAxis*)fPtDWeightsHisto->GetXaxis();
    }
  }
  if(binwidthratio<1) {
    Int_t rebin = (Int_t)(1/binwidthratio);
    hPtImpPar->RebinY(rebin);
    ptaxis = (TAxis*)hPtImpPar->GetYaxis();
  }
  
  Int_t PtBinMin = ptaxisweights->FindBin(fPtMin*1.0001);
  Int_t PtBinMax = ptaxisweights->FindBin(fPtMax*0.9999);
  
  const Int_t nPtBins = PtBinMax-PtBinMin+1;
  
  TH1F** hImpPar = new TH1F*[nPtBins];
  Int_t PtBins[nPtBins];
  for(Int_t iPt=0; iPt<nPtBins; iPt++){
    PtBins[iPt] = PtBinMin+iPt;
  }
    
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    sparse->GetAxis(fPtAxis)->SetRange(PtBins[iPt],PtBins[iPt]);
    hImpPar[iPt] = (TH1F*)sparse->Projection(fImpParAxis);
    hImpPar[iPt]->SetDirectory(0);
  }

  TH1F* hImpParReweighted = new TH1F("hImpParReweighted","",nImpParBins,ImpParMin,ImpParMax);
  
  for(Int_t iImpPar=0; iImpPar<nImpParBins; iImpPar++) {
    Double_t ImpPar = 0;
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      if(meson==kB)
        ImpPar = ImpPar+hImpPar[iPt]->GetBinContent(iImpPar+1)*fPtBWeightsHisto->GetBinContent(PtBins[iPt]); 
      else
        ImpPar = ImpPar+hImpPar[iPt]->GetBinContent(iImpPar+1)*fPtDWeightsHisto->GetBinContent(PtBins[iPt]); 
    }
    hImpParReweighted->SetBinContent(iImpPar+1,ImpPar);
  }
  
  return hImpParReweighted;
}

//___________________________________________________________________________________
void AliDplusMassd0fitter::DrawResult(Bool_t isFromTree)
{
  gStyle->SetOptFit(1);
  gStyle->SetTextSize(0.045);
  gStyle->SetTitleSize(0.05,"xy");
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.12);
  
  if(isFromTree) {
    fIntegral=fIntegral*fMCPromptSparse->GetAxis(fImpParAxis)->GetBinWidth(1);
    fBkg=fBkg*fMCPromptSparse->GetAxis(fImpParAxis)->GetBinWidth(1);
    fSig=fSig*fMCPromptSparse->GetAxis(fImpParAxis)->GetBinWidth(1);
  }
  
  TF1* ImpParPromptFunc = new TF1("ImpParPromptFunc",this,&AliDplusMassd0fitter::FunctionImpParPrompt,-1000,1000,5,"ImpParPromptFunc","FunctionImpParPrompt");
  ImpParPromptFunc->SetParameters(fPromptFractionGauss,fPromptMean,fPromptSigma,fPromptLambda,fPromptFraction*fSig);
  ImpParPromptFunc->SetLineColor(kGreen+3);
  
  TF1* ImpParBkgFunc = 0x0;
  TF1* ImpParTotFunc = 0x0;
  TF1* ImpParRecoFDFunc = 0x0;
  
  if(fFDType==kConvolution && fBkgType==kSingleGaussExpo) {
    ImpParRecoFDFunc = new TF1("ImpParRecoFDFunc",this,&AliDplusMassd0fitter::FunctionRecoImpParFD,-1000,1000,1,"ImpParRecoFDFunc","FunctionRecoImpParFD");
    ImpParRecoFDFunc->SetParameter(0,fSig*(1-fPromptFraction));
    
    ImpParBkgFunc = new TF1("ImpParBkgFunc",this,&AliDplusMassd0fitter::FunctionImpParPrompt,-1000,1000,5,"ImpParBkgFunc","FunctionImpParPrompt");
    ImpParBkgFunc->SetParameters(fBkgFractionGauss1,fBkgMean1,fBkgSigma1,fBkgLambda1,fBkg);
    
    ImpParTotFunc = new TF1("ImpParTotFunc",this,&AliDplusMassd0fitter::FitFunction,-1000,1000,2,"ImpParTotFunc","FitFunction");
  }
  else if(fFDType==kConvolution && fBkgType==kDoubleGaussExpo) {
    ImpParRecoFDFunc = new TF1("ImpParRecoFDFunc",this,&AliDplusMassd0fitter::FunctionRecoImpParFD,-1000,1000,1,"ImpParRecoFDFunc","FunctionRecoImpParFD");
    ImpParRecoFDFunc->SetParameter(0,fSig*(1-fPromptFraction));
    
    ImpParBkgFunc = new TF1("ImpParBkgFunc",this,&AliDplusMassd0fitter::FunctionImpParBkg,-1000,1000,10,"ImpParBgkFunc","FunctionImpParBkg");
    ImpParBkgFunc->SetParameters(fBkgFractionGauss1,fBkgMean1,fBkgSigma1,fBkgLambda1,fBkgFractionGauss2,fBkgMean2,fBkgSigma2,fBkgLambda2,fBkgFracFunc1,fBkg);
    ImpParTotFunc = new TF1("ImpParTotFunc",this,&AliDplusMassd0fitter::FitFunction,-1000,1000,2,"ImpParTotFunc","FitFunction");
  }
  else if(fFDType==kConvolution && fBkgType==kDoubleGaussExpoSymm) {
    ImpParRecoFDFunc = new TF1("ImpParRecoFDFunc",this,&AliDplusMassd0fitter::FunctionRecoImpParFD,-1000,1000,1,"ImpParRecoFDFunc","FunctionRecoImpParFD");
    ImpParRecoFDFunc->SetParameter(0,fSig*(1-fPromptFraction));
    
    ImpParBkgFunc = new TF1("ImpParBkgFunc",this,&AliDplusMassd0fitter::FunctionImpParBkg,-1000,1000,6,"ImpParBgkFunc","FunctionImpParBkg");
    ImpParBkgFunc->SetParameters(fBkgFractionGauss1,fBkgMean1,fBkgSigma1,fBkgLambda1,fBkgMean2,fBkg);
    ImpParTotFunc = new TF1("ImpParTotFunc",this,&AliDplusMassd0fitter::FitFunction,-1000,1000,2,"ImpParTotFunc","FitFunction");
  }
  else {
    ImpParRecoFDFunc = new TF1("ImpParRecoFDFunc",this,&AliDplusMassd0fitter::FunctionImpParPrompt,-1000,1000,5,"ImpParRecoFDFunc","FunctionRecoImpParPrompt");
    ImpParRecoFDFunc->SetParameters(fFDFraction1,fFDMean,fFDLambda1,fFDLambda2,fSig*(1-fPromptFraction));
    
    ImpParBkgFunc = new TF1("ImpParBkgFunc",this,&AliDplusMassd0fitter::FunctionImpParBkg,-1000,1000,9,"ImpParBkgFunc","FunctionImpParBkg");
    ImpParBkgFunc->SetParameters(fBkgFractionGauss1,fBkgMean1,fBkgSigma1,fBkgLambda1,fBkgFractionGauss2,fBkgMean2,fBkgSigma2,fBkgLambda2,fBkg);
    ImpParTotFunc = new TF1("ImpParTotFunc",this,&AliDplusMassd0fitter::FitFunction,-1000,1000,2,"ImpParTotFunc","FitFunction");
  }        
  
  ImpParTotFunc->SetParameters(fPromptFraction,fPromptSigma);

  ImpParPromptFunc->SetNpx(500);
  ImpParRecoFDFunc->SetNpx(500);
  ImpParBkgFunc->SetNpx(500);
  ImpParTotFunc->SetNpx(500);
  
  ImpParPromptFunc->SetLineColor(kGreen+3);
  ImpParRecoFDFunc->SetLineColor(kBlue);
  ImpParBkgFunc->SetLineColor(kMagenta);
  ImpParTotFunc->SetLineColor(kRed);

  ImpParPromptFunc->SetLineWidth(2);
  ImpParRecoFDFunc->SetLineWidth(2);
  ImpParBkgFunc->SetLineWidth(2);
  ImpParTotFunc->SetLineWidth(2);

  ImpParPromptFunc->SetLineStyle(7);
  ImpParRecoFDFunc->SetLineStyle(10);
  ImpParBkgFunc->SetLineStyle(9);
  
  TCanvas *cFit = new TCanvas("cFit","",900,900);
  cFit->SetLogy();

  if(isFromTree) {
    fImpParHisto = new TH1F("fImpParHisto","",400,-1000,1000);
    if(fMCTest){  
      fDataTree->Draw(Form("%s>>fImpParHisto",fImpParBranch.Data()));      
    }
    else { 
      Double_t massmin=fMassMean-fNSigma*fMassSigma;
      Double_t massmax=fMassMean+fNSigma*fMassSigma;
      fDataTree->Draw(Form("%s>>fImpParHisto",fImpParBranch.Data()),Form("%s>%f && %s<%f && %s>%f && %s<%f",fPtBranch.Data(),fPtMin,fPtBranch.Data(),fPtMax,fMassBranch.Data(),massmin,fMassBranch.Data(),massmax));
    }
  }

  Double_t d0min=ImpParRecoFDFunc->GetMaximum()/50;
  Double_t d0max=fImpParHisto->GetMaximum()*1.5;
  if(d0min<0 || d0min>fImpParHisto->GetMinimum())
    d0min=fImpParHisto->GetMinimum()/5;
  if(d0min<=0)
    d0min=1.e-1;
    
  fImpParHisto->GetYaxis()->SetRangeUser(d0min,d0max);
  fImpParHisto->SetLineColor(kBlack);
  fImpParHisto->SetLineWidth(2);
  fImpParHisto->SetMarkerStyle(21);
  fImpParHisto->SetMarkerSize(1.);
  if(isFromTree){
    fImpParHisto->SetStats(0);
    fImpParHisto->SetDrawOption("E1");
  }
  else
    fImpParHisto->Draw("E1");
  
  fImpParHisto->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",fPtMin,fPtMax));
  fImpParHisto->GetXaxis()->SetTitle("Imp Par XY (#mum)");
  if(fVariableBinning)
    fImpParHisto->GetYaxis()->SetTitle("Entries/bin width");
  else
    fImpParHisto->GetYaxis()->SetTitle(Form("Entries/(%0.f #mum)",fImpParHisto->GetBinWidth(5)));    
  fImpParHisto->GetXaxis()->SetNdivisions(505);
  ImpParTotFunc->Draw("same");
  ImpParPromptFunc->Draw("same");
  ImpParRecoFDFunc->Draw("same");
  if(!fSubBkg)
    ImpParBkgFunc->Draw("same");
  
  TFile *outfile = 0x0;
  TString filename;
  
  if(isFromTree)
    filename= "FitUnbinned";
  else if(!isFromTree && fSubBkg)
    filename="FitBinnedBkgSub";
  else
    filename="FitBinned";

  if(fMCTest) {
    cFit->SaveAs(Form("%s_%0.f-%0.f_genfrac_%0.2f.pdf",filename.Data(),fPtMin,fPtMax,fGenPromptFraction));
    outfile = new TFile(Form("%s_%0.f-%0.f_genfrac_%0.2f.root",filename.Data(),fPtMin,fPtMax,fGenPromptFraction),"RECREATE");
  }
  else {
    cFit->SaveAs(Form("%s_%0.f-%0.f.pdf",filename.Data(),fPtMin,fPtMax));
    outfile = new TFile(Form("%s_%0.f-%0.f.root",filename.Data(),fPtMin,fPtMax),"RECREATE");
  }
  
  cFit->Write();
  fImpParHisto->Write();
  ImpParPromptFunc->Write();
  ImpParRecoFDFunc->Write();
  ImpParBkgFunc->Write();
  ImpParTotFunc->Write();
  outfile->Close();

  delete ImpParPromptFunc;
  delete ImpParRecoFDFunc;
  delete ImpParBkgFunc;
  delete ImpParTotFunc;
  delete cFit;
  
  if(isFromTree) {
    fBkg /= fMCPromptSparse->GetAxis(fImpParAxis)->GetBinWidth(1);
    fSig /= fMCPromptSparse->GetAxis(fImpParAxis)->GetBinWidth(1);
    fIntegral /= fMCPromptSparse->GetAxis(fImpParAxis)->GetBinWidth(1); 
  }
}

//___________________________________________________________________________________
void AliDplusMassd0fitter::RebinVariableWidth()
{
  Int_t bincounter=0;
  const Int_t nbins = fImpParHisto->GetNbinsX()/2;
  vector<Double_t> d0limsleft;
  vector<Double_t> d0limsright;
  d0limsleft.push_back(fImpParHisto->GetBinLowEdge(1));
  d0limsright.push_back(fImpParHisto->GetBinLowEdge(nbins*2+1));

  while(bincounter<nbins) {
    Int_t countscounter = 0;
    while(countscounter<fnCountsMin) {
      countscounter += fImpParHisto->GetBinContent(bincounter+1);
      bincounter++;
      if(bincounter>=nbins-1) break;
    }
    d0limsleft.push_back(fImpParHisto->GetBinLowEdge(bincounter)+fImpParHisto->GetBinWidth(bincounter));
  }
  bincounter=2*nbins;
  while(bincounter>nbins+1) {
    Int_t countscounter = 0;
    while(countscounter<fnCountsMin) {
      countscounter += fImpParHisto->GetBinContent(bincounter);
      bincounter--;
      if(bincounter<=nbins+2) break;
    }
    d0limsright.push_back(fImpParHisto->GetBinLowEdge(bincounter+1));
  }

  const Int_t nvarbins = d0limsleft.size()+d0limsright.size();
  Double_t d0lims[nvarbins];
  for(UInt_t iVarBin=0; iVarBin<d0limsleft.size(); iVarBin++) 
    d0lims[iVarBin] = d0limsleft[iVarBin];
  for(UInt_t iVarBin=0; iVarBin<d0limsright.size(); iVarBin++) 
    d0lims[d0limsleft.size()+iVarBin] = d0limsright[d0limsright.size()-(iVarBin+1)];

  TH1F* hClone = (TH1F*)fImpParHisto->Clone();
  fImpParHisto = (TH1F*)hClone->Rebin(nvarbins-1,"fImpParHistoVarBin",d0lims);

  delete hClone;
}

//__________________________________________________________________________________
void AliDplusMassd0fitter::ResetAxes(THnSparseF* sparse)
{
  for(Int_t iAxis=0; iAxis<sparse->GetNdimensions(); iAxis++) {
    if(iAxis!=fPIDAxis) {
      TAxis* ax=(TAxis*)sparse->GetAxis(iAxis);
      ax->SetRange(-1,-1);
    }
  }
}

//___________________________________________________________________________________
TH1F* AliDplusMassd0fitter::GetSidebandsDist(Bool_t printSB)
{
  ResetAxes(fDataSparse);
  SetPtRange(fDataSparse);

  TAxis* massax = (TAxis*)fDataSparse->GetAxis(fMassAxis);
  Int_t massbinminR = massax->FindBin((fMassMean+fNSigmaSBLow*fMassSigma)*1.001);
  Int_t massbinmaxR = massax->FindBin((fMassMean+fNSigmaSBHigh*fMassSigma)*0.999);
  Int_t massbinminL = massax->FindBin((fMassMean-fNSigmaSBHigh*fMassSigma)*1.001);
  Int_t massbinmaxL = massax->FindBin((fMassMean-fNSigmaSBLow*fMassSigma)*0.999);

  TH1F* hLeftDist = 0x0;
  TH1F* hRightDist = 0x0;
  
  if(fSBRegion==kLeft || fSBRegion==kBoth) {
    massax->SetRange(massbinminL,massbinmaxL);
    hLeftDist=(TH1F*)fDataSparse->Projection(fImpParAxis); 
    hLeftDist->SetName("hLeftDist");
    hLeftDist->SetDirectory(0);
    hLeftDist->Sumw2();
  }
  if(fSBRegion==kRight || fSBRegion==kBoth) {
    massax->SetRange(massbinminR,massbinmaxR);
    hRightDist=(TH1F*)fDataSparse->Projection(fImpParAxis); 
    hRightDist->SetName("hRightDist");
    hRightDist->SetDirectory(0);
    hRightDist->Sumw2();
 }
  ResetAxes(fDataSparse);

  Double_t nLeft=hLeftDist->GetEntries();
  Double_t nRight=hRightDist->GetEntries();
  Double_t nTot=nLeft+nRight;

  TFile outfile(Form("Sidebands_Pt_%0.f-%0.f.root",fPtMin,fPtMax),"RECREATE");
  
  if(fSBRegion==kLeft) {
    if(printSB) { 
      outfile.cd();
      hLeftDist->Write();
      outfile.Close();
    }
    return hLeftDist;
  }
  else if(fSBRegion==kRight) {
    if(printSB) { 
      outfile.cd();
      hRightDist->Write();
      outfile.Close();
    }
    return hRightDist;
  }
  else if(fSBRegion==kBoth) {
    TH1F* hSumDist=(TH1F*)hLeftDist->Clone();
    hSumDist->SetDirectory(0);
    hSumDist->Add(hLeftDist,hRightDist,(nTot/(2*nLeft)),(nTot/(2*nRight)));
    hSumDist->SetName("hSumDist");
    if(printSB) { 
      outfile.cd();
      hLeftDist->Write();
      hRightDist->Write();
      hSumDist->Write();
      outfile.Close();
    }
    return hSumDist;
  }
  else {
    cerr << "You can keep the sidebands distribution only from left, right or both sides of the mass peak!" << endl;
    return 0x0;
  }
}

//_______________________________________________________________________________________
void AliDplusMassd0fitter::SetPtRange(THnSparseF* sparse)
{
  TAxis* ptax = (TAxis*)sparse->GetAxis(fPtAxis);
  Int_t ptbinmin=ptax->FindBin(fPtMin*1.001);
  Int_t ptbinmax=ptax->FindBin(fPtMax*0.999);
  ptax->SetRange(ptbinmin,ptbinmax);
}

//_______________________________________________________________________________________
void AliDplusMassd0fitter::SetMassRange(THnSparseF* sparse)
{
  TAxis* massax = (TAxis*)sparse->GetAxis(fMassAxis);
  Int_t massbinmin=massax->FindBin((fMassMean-fNSigma*fMassSigma)*1.001);
  Int_t massbinmax=massax->FindBin((fMassMean+fNSigma*fMassSigma)*0.999);
  massax->SetRange(massbinmin,massbinmax);
}
