#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TColor.h>
#include <TLine.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TArrayD.h>
#include <TArrow.h>
#include <TLatex.h>
#include <vector>

#include "AliHFMassFitter.h"
#endif

void DplusMass(Int_t massreb = 4, Int_t funsig = 0, Int_t funbkg = 0, Bool_t applyPID = kTRUE) {
  
  //________________________________________________________________________________________________

  const Int_t nPtBins = 10;
  const Int_t nPtLims = nPtBins+1;
  Double_t PtLims[nPtLims] = {1,2,3,4,5,6,7,8,12,16,24};

  //________________________________________________________________________________________________
  
  TFile *infile = new TFile("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root","read");
  TDirectoryFile *DplusDir = (TDirectoryFile*)infile->Get("PWG3_D2H_InvMassDplus");
  TList *DplusList = (TList*)DplusDir->Get("coutputDplus_ImpParpPbData0100");
  THnSparseF *hMassPtImpParAll = (THnSparseF*)DplusList->FindObject("hMassPtImpParAll");

  //________________________________________________________________________________________________
  
  TH1F* hMass = (TH1F*)hMassPtImpParAll->Projection(0);
  TH1F* hPt = (TH1F*)hMassPtImpParAll->Projection(1);
  TH1F* hImpPar = (TH1F*)hMassPtImpParAll->Projection(2);
  TH1F* hCosPXY = (TH1F*)hMassPtImpParAll->Projection(3); 
  TH1F* hDecLXY = (TH1F*)hMassPtImpParAll->Projection(4); 
  TH1F* hNormDecLXY = (TH1F*)hMassPtImpParAll->Projection(5); 
  TH1F* hD0d0Exp = (TH1F*)hMassPtImpParAll->Projection(6); 
  TH1F* hPID = (TH1F*)hMassPtImpParAll->Projection(7); 

  TCanvas* cVars = new TCanvas("cVars","",1200,900);
  cVars->Divide(4,2);
  cVars->cd(1);
  hMass->Draw();
  cVars->cd(2);
  hPt->Draw();
  cVars->cd(3);
  hImpPar->Draw();
  cVars->cd(4);
  hCosPXY->Draw();
  cVars->cd(5);
  hDecLXY->Draw();
  cVars->cd(6);
  hNormDecLXY->Draw();
  cVars->cd(7);
  hD0d0Exp->Draw();
  cVars->cd(8);
  hPID->Draw();

  cVars->SaveAs("Variables.eps");
  cVars->SaveAs("Variables.root");  

  delete cVars;
  
  //__________________________________________________________________________________________________
  
  TH1F** hMassPt = new TH1F*[nPtBins];
  AliHFMassFitter** Fitter = new AliHFMassFitter*[nPtBins];
  TCanvas* cMassFit = new TCanvas("cMassFit","",1200,900);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {    
    TAxis* PtAxis = (TAxis*)hMassPtImpParAll->GetAxis(1);
    Double_t ptbinmin = PtAxis->FindBin(PtLims[iPt]*1.001);
    Double_t ptbinmax = PtAxis->FindBin(PtLims[iPt+1]*0.999);
    PtAxis->SetRange(ptbinmin,ptbinmax);

    TAxis* PIDAxis = hMassPtImpParAll->GetAxis(3);
    if(applyPID) {
      PIDAxis->SetRange(2,2);//apply PID cutting on PID axis
    }
    
    hMassPt[iPt]=(TH1F*)hMassPtImpParAll->Projection(0);
    hMassPt[iPt]->Rebin(massreb);
    Double_t binwidth = hMassPt[iPt]->GetBinWidth(1)*1000; //bin width in MeV/c
    Double_t minmass = hMassPt[iPt]->GetBinLowEdge(3);
    Double_t maxmass = hMassPt[iPt]->GetBinLowEdge(hMassPt[iPt]->GetNbinsX());
    hMassPt[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",PtLims[iPt],PtLims[iPt+1]));
    hMassPt[iPt]->GetYaxis()->SetTitle(Form("Entries/(%0.f MeV/c)",binwidth));
    
    Fitter[iPt] = new AliHFMassFitter(hMassPt[iPt],minmass,maxmass,funbkg,funsig);
    Fitter[iPt]->MassFitter(0);
    
    cMassFit->Clear();
    Fitter[iPt]->DrawHere(gPad);
    
    cMassFit->SaveAs(Form("MassFit_pt_%0.f-%0.f.eps",PtLims[iPt],PtLims[iPt+1]));
    cMassFit->SaveAs(Form("MassFit_pt_%0.f-%0.f.root",PtLims[iPt],PtLims[iPt+1]));
    cMassFit->Modified();
    cMassFit->Update();
    
    //reset axes
    PtAxis->SetRange(-1,-1);
    PIDAxis->SetRange(-1,-1);
  }

  delete cMassFit;
  
}
