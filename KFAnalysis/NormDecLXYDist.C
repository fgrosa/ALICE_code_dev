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

#include "AliHFMassFitter.h"

#endif
void NormDecLXYDist() {

  TFile standardfile("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/MC/AnalysisResultspPbMC.root","READ");
  TDirectoryFile* dirMC=(TDirectoryFile*)standardfile.Get("PWG3_D2H_InvMassDplus"); 
  TList* listMC=(TList*)dirMC->Get("coutputDplus_CutVarpPbMC0100");
  THnSparseF* PromptSparse =(THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
  THnSparseF* RecoFDSparse =(THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
  standardfile.Close();
  PromptSparse->GetAxis(10)->SetRange(1,20);
  RecoFDSparse->GetAxis(10)->SetRange(1,20);

  TH1F* hPrompt = (TH1F*)PromptSparse->Projection(10);
  hPrompt->SetDirectory(0);
  hPrompt->Scale(1./hPrompt->Integral());
  TH1F* hFD = (TH1F*)RecoFDSparse->Projection(10);
  hFD->SetDirectory(0);
  hFD->Scale(1./hFD->Integral());
  
  TString MCfilename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pPb_LHC13/MC/fgrosa_Dplus_KF_pPbMC.root";
  TString MClistname = "coutputDplusKF";
  TFile infile(MCfilename.Data(),"UPDATE");
  TList* list = (TList*)infile.Get(MClistname.Data());  

  THnSparseF *PromptSparseKF = (THnSparseF*)list->FindObject("fSparsePrompt");
  THnSparseF *FDSparseKF = (THnSparseF*)list->FindObject("fSparseFD");

  TH1F* hKFPrompt = (TH1F*)PromptSparseKF->Projection(6);
  hKFPrompt->Rebin(2);
  hKFPrompt->SetDirectory(0);
  hKFPrompt->SetLineColor(kRed);
  hKFPrompt->Scale(1./hKFPrompt->Integral());
  TH1F* hKFFD = (TH1F*)FDSparseKF->Projection(6);
  hKFFD->Rebin(2);
  hKFFD->SetDirectory(0);
  hKFFD->SetLineColor(kRed);
  hKFFD->Scale(1./hKFFD->Integral());

  TCanvas* cPrompt = new TCanvas("cPrompt","",1200,900);
  hPrompt->GetYaxis()->SetRangeUser(0.,0.2);
  hPrompt->Draw();
  hKFPrompt->Draw("same");
  TCanvas* cFD = new TCanvas("cFD","",1200,900);
  hFD->GetYaxis()->SetRangeUser(0.,0.1);
  hFD->Draw();
  hKFFD->Draw("same");
}
