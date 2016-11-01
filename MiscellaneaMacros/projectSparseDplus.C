#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TFile.h>
#include <THnSparse.h>
#include <TString.h>
#include <TList.h>
#include <Riostream.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TH1.h>
#include <TDirectoryFile.h>

#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFMassFitter.h"
#include "AliNormalizationCounter.h"

#endif

///////////////////////////////////////////
// main functions:                       //
// - plotEfficiency()                    //
// - plotInvMassSpectra()                //
///////////////////////////////////////////

Double_t projectAccStep(THnSparse *hSparse, Int_t ipt);
TH1D* projectRecoStep(THnSparse *hSparse, Int_t ipt);
TH1F* RebinHisto(TH1F* hOrig, Int_t reb, Int_t firstUse=-1);

///////////////////////////////////////////
/////         Source files           //////
///////////////////////////////////////////

TString fileeff = "/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/MC/AnalysisResultspPbMC.root";
TString direff  = "PWG3_D2H_InvMassDplus";
TString listeff = "coutputDplus_ImpParpPbMC0100";

TString filedata = "/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root";
TString dirdata  = "PWG3_D2H_InvMassDplus";
TString listdata = "coutputDplus_ImpParpPbData0100";

///////////////////////////////////////////
/////     Set here the cut values    //////
///////////////////////////////////////////

const Int_t nptbins = 10;
Double_t ptlims[nptbins+1] = {1.,2.,3.,4.,5.,6.,7.,8.,12.,16.,24.};
Double_t DLxy[nptbins]     = {-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.};
Double_t NDLxy[nptbins]    = {-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.};
Double_t CosPxy[nptbins]   = {-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.};
//Double_t MinDD0[nptbins]   = {-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.};
//Double_t MaxDD0[nptbins]   = {1.,1.,1.,1.,1.,1.,1.,1.};
//Double_t MinDD0[nptbins]   = {-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5};
//Double_t MaxDD0[nptbins]   = {1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5};
//Double_t MinDD0[nptbins]   = {-2.,-2.,-2.,-2.,-2.,-2.,-2.,-2.};
//Double_t MaxDD0[nptbins]   = {2.,2.,2.,2.,2.,2.,2.,2.};
Double_t MinDD0[nptbins]   = {-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.};
Double_t MaxDD0[nptbins]   = {999.,999.,999.,999.,999.,999.,999.,999.,999.,999.};
Bool_t usePID[nptbins]   = {1,1,1,1,1,1,1,1,1,1};
Double_t Y = 0.9;

///////////////////////////////////////////
/////  Settings for IM spectra fits  //////
///////////////////////////////////////////

enum {kD0toKpi, kDplusKpipi, kDStarD0pi, kDsKKpi};
enum {kBoth, kParticleOnly, kAntiParticleOnly};
enum {kExpo=0, kLinear, kPol2};
enum {kGaus=0, kDoubleGaus};

// Common variables: to be configured by the user
Int_t rebin[nptbins]={5,5,5,5,5,5,5,6,6,7};
Double_t lowLim[nptbins]={1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7};
Double_t upLim[nptbins]={2.05,2.05,2.05,2.05,2.05,2.05,2.05,2.05,2.05,2.05};
Int_t firstUsedBin[nptbins]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

Double_t sigmaMC[nptbins] = {0.0082,0.010,0.012,0.0136}; //pass2

Int_t typeb=kExpo;
Int_t types=kGaus;
Int_t optPartAntiPart=kBoth;
Int_t factor4refl=0;
Double_t nEventsForNorm=0.;

//___________________________________________________________

void plotEfficiency(){
  
  TFile *infileEff = TFile::Open(fileeff.Data());
  TDirectoryFile *dirEff = (TDirectoryFile*)infileEff->Get(direff.Data());
  TList *listEff = (TList*)dirEff->Get(listeff.Data());
  listEff->ls();

  THnSparse* hSparseAcc_C = (THnSparse*)listEff->FindObject("hMCAccPrompt");
  THnSparse* hSparseAcc_B = (THnSparse*)listEff->FindObject("hMCAccBFeed");
  THnSparse* hSparseReco_C = (THnSparse*)listEff->FindObject("hMassPtImpParPrompt");
  THnSparse* hSparseReco_B = (THnSparse*)listEff->FindObject("hMassPtImpParBfeed");
  
  TH1D* hAccStep_C  = new TH1D("hAccStep_C","",nptbins,ptlims);
  TH1D* hAccStep_B  = new TH1D("hAccStep_B","",nptbins,ptlims);
  TH1D* hRecoStep_C = new TH1D("hRecoStep_C","",nptbins,ptlims);
  TH1D* hRecoStep_B = new TH1D("hRecoStep_B","",nptbins,ptlims);
  TH1D* hEff_C = new TH1D("hEff_C","",nptbins,ptlims);
  TH1D* hEff_B = new TH1D("hEff_B","",nptbins,ptlims);

  TH1D* hMass;
  
  for(Int_t ipt=0; ipt<nptbins; ipt++) {
    Double_t accStep_C  = projectAccStep(hSparseAcc_C,ipt);
    Double_t accStep_B  = projectAccStep(hSparseAcc_B,ipt);
    hMass = (TH1D*)projectRecoStep(hSparseReco_C,ipt);
    hMass->SetName(Form("hMassPromptPtBin%d",ipt));
    Double_t recoStep_C = hMass->GetEntries();
    hMass = (TH1D*)projectRecoStep(hSparseReco_B,ipt);
    hMass->SetName(Form("hMassBfeedPtBin%d",ipt));
    Double_t recoStep_B = hMass->GetEntries();
    
    hAccStep_C->SetBinContent(ipt+1,accStep_C);
    hAccStep_B->SetBinContent(ipt+1,accStep_B);
    hRecoStep_C->SetBinContent(ipt+1,recoStep_C);
    hRecoStep_B->SetBinContent(ipt+1,recoStep_B);
    hAccStep_C->SetBinError(ipt+1,TMath::Sqrt(accStep_C));
    hAccStep_B->SetBinError(ipt+1,TMath::Sqrt(accStep_B));
    hRecoStep_C->SetBinError(ipt+1,TMath::Sqrt(recoStep_C));
    hRecoStep_B->SetBinError(ipt+1,TMath::Sqrt(recoStep_B));
  }
   
  TFile fout("efficiencies.root","recreate");
  fout.cd();
  hAccStep_C->Write();
  hAccStep_B->Write();
  hRecoStep_C->Write();
  hRecoStep_B->Write();
  hRecoStep_C->Divide(hRecoStep_C,hAccStep_C,1,1,"B");
  hRecoStep_B->Divide(hRecoStep_B,hAccStep_B,1,1,"B");
  for (Int_t ii=1; ii<=nptbins; ii++) {
    hEff_C->SetBinContent(ii,hRecoStep_C->GetBinContent(ii));
    hEff_B->SetBinContent(ii,hRecoStep_B->GetBinContent(ii));
    hEff_C->SetBinError(ii,hRecoStep_C->GetBinError(ii));
    hEff_B->SetBinError(ii,hRecoStep_B->GetBinError(ii));
  }
  hEff_C->SetStats(0);
  hEff_B->SetStats(0);
  hEff_C->Write();
  hEff_B->Write();
  fout.Close();
  
  
  hRecoStep_C->Divide(hAccStep_C);
  hRecoStep_B->Divide(hAccStep_B);
  
  
  TCanvas *c = new TCanvas("c","",600,400);
  hEff_C->SetLineColor(kRed);
  hEff_B->SetLineColor(kBlue);
  hEff_C->SetLineWidth(2);
  hEff_B->SetLineWidth(2);
  hEff_C->SetMaximum(0.2);
  hEff_C->Draw();
  hEff_B->Draw("same");
  TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(hEff_C,"prompt");
  leg->AddEntry(hEff_B,"feeddw");
  leg->Draw("same");
   
  delete hMass;
  
}


//___________________________________________________________

void plotInvMassSpectraMC(){
 
  TFile *infileEff = TFile::Open(fileeff.Data());
  TDirectoryFile *dirEff = (TDirectoryFile*)infileEff->Get(direff.Data());
  TList *listEff = (TList*)dirEff->Get(listeff.Data());
  THnSparse* hSparse = (THnSparse*)listEff->FindObject("hMassPtImpParPrompt");
 
  TH1D* hMass         = new TH1D("hMass","hMass",nptbins,ptlims);
  TH1D* hSigma        = new TH1D("hSigma","hSigma",nptbins,ptlims);

  Double_t massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  
  TH1D *hmass[nptbins];
  AliHFMassFitter** fitter = new AliHFMassFitter*[nptbins];
  
  TCanvas *c = new TCanvas("c","InvMass",1500,900);
  c->Divide(4,3);
  for(Int_t ipt=0; ipt<nptbins; ipt++) {
    
    hmass[ipt] = (TH1D*)projectRecoStep(hSparse,ipt);
    hmass[ipt]->SetName(Form("hInvMass_PtBin%d",ipt));
    c->cd(ipt+1);
    TH1F* hRebinned = RebinHisto((TH1F*)hmass[ipt],1,-1);
    hRebinned->Fit("gaus");
    TF1* fg=(TF1*)hRebinned->GetListOfFunctions()->FindObject("gaus");
    hMass->SetBinContent(ipt+1,fg->GetParameter(1));
    hMass->SetBinError(ipt+1,fg->GetParError(1));    
    hSigma->SetBinContent(ipt+1,fg->GetParameter(2));
    hSigma->SetBinError(ipt+1,fg->GetParError(2));    
  }

  TCanvas *cpar=new TCanvas("cpar","Fit params",1200,600);
  cpar->Divide(2,1);
  cpar->cd(1);
  hMass->SetMarkerStyle(20);
  hMass->Draw("PE");
  hMass->GetXaxis()->SetTitle("Pt (GeV/c)");
  hMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
  cpar->cd(2);
  hSigma->SetMarkerStyle(20);
  hSigma->Draw("PE");
  hSigma->GetXaxis()->SetTitle("Pt (GeV/c)");
  hSigma->GetXaxis()->SetTitle("Sigma (GeV/c^{2})");

  TFile* outf=new TFile("GaussFitMC.root","recreate");
  outf->cd();
  hMass->Write();
  hSigma->Write();
  outf->Close();
}

//___________________________________________________________

void plotInvMassSpectra(){
  
  TFile *infileData = TFile::Open(filedata.Data());
  TDirectoryFile *dirData = (TDirectoryFile*)infileData->Get(dirdata.Data());
  TList *listData = (TList*)dirData->Get(listdata.Data());
  TString listcuts=listdata.Data();
  listcuts.ReplaceAll("coutputDplus","coutputDplusCuts");
  TList *listCuts = (TList*)dirData->Get(listcuts.Data());
  THnSparse* hSparse = (THnSparse*)listData->FindObject("hMassPtImpParAll");
  TH1F** hmass1d=new TH1F*[nptbins];
  for(Int_t i=0;i<nptbins;i++) hmass1d[i]=0x0;
  Bool_t retCode=LoadDplusHistos(listCuts,listData,hmass1d);

  TH1D* hSignal       = new TH1D("hSignal","hSignal",nptbins,ptlims);
  TH1D* hRelErrSig    = new TH1D("hRelErrSig","hRelErrSig",nptbins,ptlims);
  TH1D* hInvSignif    = new TH1D("hInvSignif","hInvSignif",nptbins,ptlims);
  TH1D* hBackground   = new TH1D("hBackground","hBackground",nptbins,ptlims);
  TH1D* hSignificance = new TH1D("hSignificance","hSignificance",nptbins,ptlims);
  TH1D* hMass         = new TH1D("hMass","hMass",nptbins,ptlims);
  TH1D* hSigma        = new TH1D("hSigma","hSigma",nptbins,ptlims);

  TH1D* hCntSig1=new TH1D("hCntSig1","hCntSig1",nptbins,ptlims);
  TH1D* hCntSig2=new TH1D("hCntSig2","hCntSig2",nptbins,ptlims);
  TH1D* hNDiffCntSig1=new TH1D("hNDiffCntSig1","hNDiffCntSig1",nptbins,ptlims);
  TH1D* hNDiffCntSig2=new TH1D("hNDiffCntSig2","hNDiffCntSig2",nptbins,ptlims);
  
  Double_t massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  
  TH1D *hmass[nptbins];
  AliHFMassFitter** fitter = new AliHFMassFitter*[nptbins];
  
  TCanvas *c = new TCanvas("c","InvMass",1500,900);
  c->Divide(4,3);

  // for(Int_t ipt=0; ipt<nptbins; ipt++) {
  //   c->cd(ipt+1);
  //   hmass[ipt] = (TH1D*)projectRecoStep(hSparse,ipt);
  //   hmass[ipt]->Draw("e");
  //   hmass1d[ipt]->SetLineColor(2);
  //   hmass1d[ipt]->Draw("esames");
  // }
  // return;

  Double_t sig,errsig,s,errs,b,errb;

  TFile MCfile("GaussFitMC.root","UPDATE");
  TH1F* hSigmaMC = (TH1F*)MCfile.Get("hSigma");
  hSigmaMC->SetDirectory(0);
  MCfile.Close();
  
  for(Int_t ipt=0; ipt<nptbins; ipt++) {
    
    hmass[ipt] = (TH1D*)projectRecoStep(hSparse,ipt);
    hmass[ipt]->SetName(Form("hInvMass_PtBin%d",ipt));
    c->cd(ipt+1);
    TH1F* hRebinned = RebinHisto((TH1F*)hmass[ipt],rebin[ipt],firstUsedBin[ipt]);
    fitter[ipt]=new AliHFMassFitter(hRebinned,lowLim[ipt],upLim[ipt],1,typeb,types);
    fitter[ipt]->SetReflectionSigmaFactor(factor4refl);
    fitter[ipt]->SetInitialGaussianMean(massD);
    fitter[ipt]->SetInitialGaussianSigma(hSigmaMC->GetBinContent(ipt+1));
    fitter[ipt]->SetUseLikelihoodFit();
    
    Bool_t out=fitter[ipt]->MassFitter(0);
    if(!out) {
      fitter[ipt]->GetHistoClone()->Draw();
      continue;
    }
    fitter[ipt]->DrawHere(gPad);
    fitter[ipt]->Signal(3,s,errs);
    fitter[ipt]->Background(3,b,errb);
    fitter[ipt]->Significance(3,sig,errsig);
    Double_t mass  = fitter[ipt]->GetMean();
    Double_t sigma = fitter[ipt]->GetSigma();
    Double_t ry=fitter[ipt]->GetRawYield();
    Double_t ery=fitter[ipt]->GetRawYieldError();
    Double_t cntSig1=0.;
    Double_t cntSig2=0.;
    Double_t cntErr=0.;
    Double_t massRangeForCounting=3.5*sigma;
    Double_t minBinSum=hmass[ipt]->FindBin(massD-massRangeForCounting);
    Double_t maxBinSum=hmass[ipt]->FindBin(massD+massRangeForCounting);
    TF1* fB1=fitter[ipt]->GetBackgroundFullRangeFunc();
    TF1* fB2=fitter[ipt]->GetBackgroundRecalcFunc();
     for(Int_t iMB=minBinSum; iMB<=maxBinSum; iMB++){
      Double_t bkg1=fB1 ? fB1->Eval(hmass[ipt]->GetBinCenter(iMB))/rebin[ipt] : 0;
      Double_t bkg2=fB2 ? fB2->Eval(hmass[ipt]->GetBinCenter(iMB))/rebin[ipt] : 0;
      cntSig1+=(hmass[ipt]->GetBinContent(iMB)-bkg1);
      cntSig2+=(hmass[ipt]->GetBinContent(iMB)-bkg2);
      cntErr+=(hmass[ipt]->GetBinContent(iMB));
    }
    hCntSig1->SetBinContent(ipt+1,cntSig1);
    hCntSig1->SetBinError(ipt+1,TMath::Sqrt(cntErr));
    hNDiffCntSig1->SetBinContent(ipt+1,(s-cntSig1)/s);
    hNDiffCntSig1->SetBinError(ipt+1,TMath::Sqrt(cntErr)/s);
    hCntSig2->SetBinContent(ipt+1,cntSig2);
    hNDiffCntSig2->SetBinContent(ipt+1,(s-cntSig2)/s);
    hNDiffCntSig2->SetBinError(ipt+1,TMath::Sqrt(cntErr)/s);
    hCntSig2->SetBinError(ipt+1,TMath::Sqrt(cntErr));
    
    hSignal->SetBinContent(ipt+1,s);
    hSignal->SetBinError(ipt+1,errs);
    hRelErrSig->SetBinContent(ipt+1,errs/s);
    hInvSignif->SetBinContent(ipt+1,1/sig);
    hInvSignif->SetBinError(ipt+1,errsig/(sig*sig));
    hBackground->SetBinContent(ipt+1,b);
    hBackground->SetBinError(ipt+1,errb);
    hSignificance->SetBinContent(ipt+1,sig);
    hSignificance->SetBinError(ipt+1,errsig);
    hMass->SetBinContent(ipt+1,mass);
    hMass->SetBinError(ipt+1,fitter[ipt]->GetMeanUncertainty());
    hSigma->SetBinContent(ipt+1,sigma);
    hSigma->SetBinError(ipt+1,fitter[ipt]->GetSigmaUncertainty());
  }
  TString normobj=listdata;
  normobj.ReplaceAll("coutputDplus","coutputDplusNorm");
  AliNormalizationCounter* nc=(AliNormalizationCounter*)dirData->Get(normobj.Data());
  printf("Events for normalization = %f\n",nc->GetNEventsForNorm());
  nEventsForNorm+=nc->GetNEventsForNorm();
  printf("Events for norm = %f\n",nEventsForNorm);
  TH1F* hNEvents=new TH1F("hNEvents","",1,0.,1.);
  hNEvents->SetBinContent(1,nEventsForNorm);
  
  
  TFile* outf=new TFile("RawYield.root","recreate");
  outf->cd();
  for(Int_t ii=0;ii<nptbins;ii++)hmass[ii]->Write();
  hNEvents->Write();
  hMass->Write();
  hSigma->Write();
  hSignal->Write();
  hRelErrSig->Write();
  hInvSignif->Write();
  hBackground->Write();
  hSignificance->Write();
  outf->Close();

  TCanvas *cpar=new TCanvas("cpar","Fit params",1200,600);
  cpar->Divide(2,1);
  cpar->cd(1);
  hMass->SetMarkerStyle(20);
  hMass->Draw("PE");
  hMass->GetXaxis()->SetTitle("Pt (GeV/c)");
  hMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
  cpar->cd(2);
  hSigma->SetMarkerStyle(20);
  hSigma->Draw("PE");
  hSigma->GetXaxis()->SetTitle("Pt (GeV/c)");
  hSigma->GetXaxis()->SetTitle("Sigma (GeV/c^{2})");

  TCanvas* csig=new TCanvas("csig","Results",1200,600);
  csig->Divide(3,1);
  csig->cd(1);
  hSignal->SetMarkerStyle(20);
  hSignal->SetMarkerColor(4);
  hSignal->SetLineColor(4);
  hSignal->GetXaxis()->SetTitle("Pt (GeV/c)");
  hSignal->GetYaxis()->SetTitle("Signal");
  hSignal->Draw("P");
  hCntSig1->SetMarkerStyle(26);
  hCntSig1->SetMarkerColor(2);
  hCntSig1->SetLineColor(2);
  hCntSig1->Draw("PSAME");
  hCntSig2->SetMarkerStyle(29);
  hCntSig2->SetMarkerColor(kGray+1);
  hCntSig2->SetLineColor(kGray+1);
  hCntSig2->Draw("PSAME");
  TLegend* leg=new TLegend(0.4,0.7,0.89,0.89);
  leg->SetFillColor(0);
  TLegendEntry* ent=leg->AddEntry(hSignal,"From Fit","PL");
  ent->SetTextColor(hSignal->GetMarkerColor());
  ent=leg->AddEntry(hCntSig1,"From Counting1","PL");
  ent->SetTextColor(hCntSig1->GetMarkerColor());
  ent=leg->AddEntry(hCntSig2,"From Counting2","PL");
  ent->SetTextColor(hCntSig2->GetMarkerColor());
  leg->Draw();
  csig->cd(2);
  hBackground->SetMarkerStyle(20);
  hBackground->Draw("P");
  hBackground->GetXaxis()->SetTitle("Pt (GeV/c)");
  hBackground->GetYaxis()->SetTitle("Background");
  csig->cd(3);
  hSignificance->SetMarkerStyle(20);
  hSignificance->Draw("P");
  hSignificance->GetXaxis()->SetTitle("Pt (GeV/c)");
  hSignificance->GetYaxis()->SetTitle("Significance");

  TCanvas* cDiffS=new TCanvas("cDiffS","Difference",1200,600);
  cDiffS->Divide(2,1);
  cDiffS->cd(1);
  hRelErrSig->SetMarkerStyle(20); //fullcircle
  hRelErrSig->SetTitleOffset(1.2);  
  hRelErrSig->SetTitle("Rel Error from Fit;p_{t} (GeV/c);Signal Relative Error");
  hRelErrSig->Draw("P");
  hInvSignif->SetMarkerStyle(21); //fullsquare
  hInvSignif->SetMarkerColor(kMagenta+1);
  hInvSignif->SetLineColor(kMagenta+1);
  hInvSignif->Draw("PSAMES");
  TLegend* leg2=new TLegend(0.4,0.7,0.89,0.89);
  leg2->SetFillColor(0);
  TLegendEntry* ent2=leg2->AddEntry(hRelErrSig,"From Fit","P");
  ent2->SetTextColor(hRelErrSig->GetMarkerColor());
  ent2=leg2->AddEntry(hInvSignif,"1/Significance","PL");
  ent2->SetTextColor(hInvSignif->GetMarkerColor());
  leg2->Draw();

  cDiffS->cd(2);
  hNDiffCntSig1->SetMarkerStyle(26);
  hNDiffCntSig1->SetMarkerColor(2);
  hNDiffCntSig1->SetLineColor(2);
  hNDiffCntSig1->SetTitle("Cmp Fit-Count;p_{t} (GeV/c);(S_{fit}-S_{count})/S_{fit}");
  hNDiffCntSig1->Draw("P");
  hNDiffCntSig2->SetMarkerStyle(29);
  hNDiffCntSig2->SetMarkerColor(kGray+1);
  hNDiffCntSig2->SetLineColor(kGray+1);
  hNDiffCntSig2->Draw("PSAME");
  TLegend* leg1=new TLegend(0.4,0.7,0.89,0.89);
  leg1->SetFillColor(0);
  TLegendEntry* ent1=leg1->AddEntry(hNDiffCntSig1,"From Counting1","PL");
  ent1->SetTextColor(hNDiffCntSig1->GetMarkerColor());
  ent1=leg1->AddEntry(hNDiffCntSig2,"From Counting2","PL");
  ent1->SetTextColor(hNDiffCntSig2->GetMarkerColor());
  leg1->Draw();

}


//___________________________________________________________

Double_t projectAccStep(THnSparse *hSparse, Int_t ipt){
  
  TAxis* axPt=(TAxis*)hSparse->GetAxis(0);
  TAxis* axY=(TAxis*)hSparse->GetAxis(1);
  
  Int_t first=0;
  Int_t last=0;
  Double_t binWidth;
  
  //Pt
  binWidth=axPt->GetBinWidth(1);
  first=axPt->FindBin(ptlims[ipt]+binWidth/2);
  last=axPt->FindBin(ptlims[ipt+1]+binWidth/2)-1;
  axPt->SetRange(first,last);
  cout<<"pt  "<<axPt->GetBinLowEdge(first)<<"    "<<axPt->GetBinUpEdge(last)<<endl;
  
  //Y
  binWidth=axY->GetBinWidth(1);
  first=axY->FindBin(-Y-binWidth/2)+1;//GetFirst();
  last=axY->FindBin(Y+binWidth/2)-1;//GetLast();
  axY->SetRange(0,axY->GetNbins()+1);
  cout<<"Y  "<<axY->GetBinLowEdge(0)<<"    "<<axY->GetBinUpEdge(axY->GetNbins()+1)<<endl;
  
  TH1D* hproj = (TH1D*)hSparse->Projection(0);
  Double_t entries = hproj->GetEntries();
  delete hproj;
  
  return entries;
}
//___________________________________________________________

TH1D *projectRecoStep(THnSparse *hSparse, Int_t ipt){
  Int_t nAxes=hSparse->GetNdimensions();
  Int_t nbins;
  for(Int_t iax=0; iax<nAxes; iax++){
    printf("Axis %d: %s\n",iax,hSparse->GetAxis(iax)->GetTitle());
    nbins = hSparse->GetAxis(iax)->GetNbins();
    hSparse->GetAxis(iax)->SetRange(0,nbins+1);
  }


  TAxis* axIMDs=(TAxis*)hSparse->GetAxis(0);
  TAxis* axPt=(TAxis*)hSparse->GetAxis(1);
  TAxis* axDip=(TAxis*)hSparse->GetAxis(2);
  TAxis* axCosPxy=(TAxis*)hSparse->GetAxis(7);
  TAxis* axDLxy=(TAxis*)hSparse->GetAxis(9);
  TAxis* axNDLxy=(TAxis*)hSparse->GetAxis(10);
  TAxis* axDeltaD0=(TAxis*)hSparse->GetAxis(11);
  TAxis* axPID=(TAxis*)hSparse->GetAxis(3);

  
  Int_t first=0;
  Int_t last=0;
  Double_t binWidth;
  
  
  
  //Mass
  axIMDs->SetRange(1,axIMDs->GetNbins());

  //Pt
  binWidth=axPt->GetBinWidth(1);
  first=axPt->FindBin(ptlims[ipt]+binWidth/2);
  last=axPt->FindBin(ptlims[ipt+1]+binWidth/2)-1;
  axPt->SetRange(first,last);
  cout<<"pt  "<<axPt->GetBinLowEdge(first)<<"    "<<axPt->GetBinUpEdge(last)<<endl;
  
  //DLxy
  binWidth=axDLxy->GetBinWidth(1);
  first=axDLxy->FindBin(DLxy[ipt]+binWidth/2);
  last=axDLxy->GetNbins()+1;//  GetLast();
  axDLxy->SetRange(first,last);
  cout<<"dec l xy "<< axDLxy->GetBinLowEdge(first)<<"    "<<axDLxy->GetBinUpEdge(last)<<endl;
   
  //NDLxy
  binWidth=axNDLxy->GetBinWidth(1);
  first=axNDLxy->FindBin(NDLxy[ipt]+binWidth/2);
  last=axNDLxy->GetNbins()+1;//  GetLast();
  axNDLxy->SetRange(first,last);
  cout<<"NDLxy  "<< axNDLxy->GetBinLowEdge(first)<<"    "<<axNDLxy->GetBinUpEdge(last)<<endl;
  
  //cosPxy
  binWidth=axCosPxy->GetBinWidth(1);
  first=axCosPxy->FindBin(CosPxy[ipt]+binWidth/2);
  last=axCosPxy->GetNbins()+1;//  GetLast();
  axCosPxy->SetRange(first,last);
  cout<<"cosPxy  "<< axCosPxy->GetBinLowEdge(first)<<"    "<<axCosPxy->GetBinUpEdge(last)<<endl;

  //d0meas-d0exp
  binWidth=axDeltaD0->GetBinWidth(1);
  first=axDeltaD0->FindBin(MinDD0[ipt]+binWidth/2);
  last=axDeltaD0->FindBin(MaxDD0[ipt]-binWidth/2);
  axDeltaD0->SetRange(first,last);
  cout<<"deltaD0  "<< axDeltaD0->GetBinLowEdge(first)<<"    "<<axDeltaD0->GetBinUpEdge(last)<<endl;
  
  if(usePID[ipt]){
    first=2;
    last=2;
    axPID->SetRange(first,last);
    cout<<"pid  "<< axPID->GetBinLowEdge(first)<<"    "<<axPID->GetBinUpEdge(last)<<endl;
  }

  TH1D* hproj = (TH1D*)hSparse->Projection(0);
  hproj->SetTitle(Form("%.1f<pt<%.1f",ptlims[ipt],ptlims[ipt+1]));
  return hproj;
}


//___________________________________________________________

TH1F* RebinHisto(TH1F* hOrig, Int_t reb, Int_t firstUse){
  // Rebin histogram, from bin firstUse to lastUse
  // Use all bins if firstUse=-1
  
  Int_t nBinOrig=hOrig->GetNbinsX();
  Int_t firstBinOrig=1;
  Int_t lastBinOrig=nBinOrig;
  Int_t nBinOrigUsed=nBinOrig;
  Int_t nBinFinal=nBinOrig/reb;
  if(firstUse>=1){
    firstBinOrig=firstUse;
    nBinFinal=(nBinOrig-firstUse+1)/reb;
    nBinOrigUsed=nBinFinal*reb;
    lastBinOrig=firstBinOrig+nBinOrigUsed-1;
  }else{
    Int_t exc=nBinOrigUsed%reb;
    if(exc!=0){
      nBinOrigUsed-=exc;
      firstBinOrig+=exc/2;
      lastBinOrig=firstBinOrig+nBinOrigUsed-1;
    }
  }
  
  printf("Rebin from %d bins to %d bins -- Used bins=%d in range %d-%d\n",nBinOrig,nBinFinal,nBinOrigUsed,firstBinOrig,lastBinOrig);
  Float_t lowLim=hOrig->GetXaxis()->GetBinLowEdge(firstBinOrig);
  Float_t hiLim=hOrig->GetXaxis()->GetBinUpEdge(lastBinOrig);
  TH1F* hRebin=new TH1F(Form("%s-rebin",hOrig->GetName()),hOrig->GetTitle(),nBinFinal,lowLim,hiLim);
  Int_t lastSummed=firstBinOrig-1;
  for(Int_t iBin=1;iBin<=nBinFinal; iBin++){
    Float_t sum=0.;
    for(Int_t iOrigBin=0;iOrigBin<reb;iOrigBin++){
      sum+=hOrig->GetBinContent(lastSummed+1);
      lastSummed++;
    }
    hRebin->SetBinContent(iBin,sum);
  }
  return hRebin;
}


Bool_t LoadDplusHistos(TList* listcut, TList* hlist, TH1F** hMass){

  AliRDHFCutsDplustoKpipi* dcuts=(AliRDHFCutsDplustoKpipi*)listcut->At(0);

  Int_t nPtBinsCuts=dcuts->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBinsCuts);
  Float_t *ptlimsCuts=dcuts->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;

  Int_t iFinBin=0;
  for(Int_t i=0;i<nPtBinsCuts;i++){
    if(ptlimsCuts[i]>=ptlims[iFinBin+1]) iFinBin+=1; 
    if(iFinBin>nptbins) break;
    if(ptlimsCuts[i]>=ptlims[iFinBin] && 
       ptlimsCuts[i+1]<=ptlims[iFinBin+1]){
      TString histoName;
      if(optPartAntiPart==kBoth) histoName.Form("hMassPt%dTC",i);
      else if(optPartAntiPart==kParticleOnly) histoName.Form("hMassPt%dTCPlus",i);
      else if(optPartAntiPart==kAntiParticleOnly) histoName.Form("hMassPt%dTCMinus",i);
      TH1F* htemp=(TH1F*)hlist->FindObject(histoName.Data());
      if(!htemp){
	printf("ERROR: Histogram %s not found\n",histoName.Data());
	return kFALSE;
      }
      if(!hMass[iFinBin]){
	hMass[iFinBin]=new TH1F(*htemp);
      }else{
	hMass[iFinBin]->Add(htemp);
      }
    }
  }

  return kTRUE;

}
