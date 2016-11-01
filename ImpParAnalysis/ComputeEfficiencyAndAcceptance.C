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

void ComputeEfficiencyAndAcceptance(Int_t PIDaxis=3) {

  //______________________________________________________________________________________
  //pT bins of the analysis

  Int_t nPtBins = 7;
  const Int_t nPtLims = nPtBins+1;
  Double_t PtLims[nPtLims] = {2,3,4,5,6,8,12,16};

  //______________________________________________________________________________________
  //input file
  
  TFile infileMC("../../files/AnalysisResultspPbMC.root","READ");
  TDirectoryFile* dirMC=(TDirectoryFile*)infileMC.Get("PWG3_D2H_InvMassDplus"); 
  TList* listMC=(TList*)dirMC->Get("coutputDplus_ImpParpPbMC0100");
  THnSparseF* RecoPromptSparse = (THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
  THnSparseF* RecoFDSparse = (THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
  THnSparseF* GenAccPromptSparse = (THnSparseF*)listMC->FindObject("hMCAccPrompt");
  THnSparseF* GenAccFDSparse = (THnSparseF*)listMC->FindObject("hMCAccBFeed");
  infileMC.Close();

  TFile accfile("Acceptance_Toy_DplusKpipi_yfidPtDep_etaDau09_ptDau100_FONLL5ptshape.root","READ");
  TH1F* hPtGenLimAcc = (TH1F*)accfile.Get("hPtGenLimAcc");
  TH1F* hPtGenAcc = (TH1F*)accfile.Get("hPtGenAcc");
  hPtGenLimAcc->SetDirectory(0);
  hPtGenAcc->SetDirectory(0);
  hPtGenLimAcc->Sumw2();
  hPtGenAcc->Sumw2();
  accfile.Close();
  
  //______________________________________________________________________________________
  //efficiency calculation
  
  for(Int_t iAxis=0; iAxis<RecoPromptSparse->GetNdimensions(); iAxis++) {
    RecoPromptSparse->GetAxis(iAxis)->SetRange(-1,-1);
    RecoFDSparse->GetAxis(iAxis)->SetRange(-1,-1);
  }

  RecoPromptSparse->GetAxis(PIDaxis)->SetRange(2,2);
  RecoFDSparse->GetAxis(PIDaxis)->SetRange(2,2);
  
  TH1F* hPtRecoPrompt = (TH1F*)RecoPromptSparse->Projection(1);
  TH1F* hPtRecoFD = (TH1F*)RecoFDSparse->Projection(1);
  TH1F* hPtGenAccPrompt = (TH1F*)GenAccPromptSparse->Projection(0);
  TH1F* hPtGenAccFD = (TH1F*)GenAccFDSparse->Projection(0);
  hPtRecoPrompt->SetDirectory(0);
  hPtRecoFD->SetDirectory(0);
  hPtGenAccPrompt->SetDirectory(0);
  hPtGenAccFD->SetDirectory(0);
  hPtRecoPrompt->Sumw2();
  hPtRecoFD->Sumw2();
  hPtGenAccPrompt->Sumw2();
  hPtGenAccFD->Sumw2();
  
  TH1F* hPtRecoPromptReb = (TH1F*)hPtRecoPrompt->Rebin(nPtBins,"hPtRecoPromptReb",PtLims);
  TH1F* hPtRecoFDReb = (TH1F*)hPtRecoFD->Rebin(nPtBins,"hPtRecoFDReb",PtLims);
  TH1F* hPtGenAccPromptReb = (TH1F*)hPtGenAccPrompt->Rebin(nPtBins,"hPtGenAccPromptReb",PtLims);
  TH1F* hPtGenAccFDReb = (TH1F*)hPtGenAccFD->Rebin(nPtBins,"hPtGenAccFDReb",PtLims);
  hPtRecoPromptReb->SetDirectory(0);
  hPtRecoFDReb->SetDirectory(0);
  hPtGenAccPromptReb->SetDirectory(0);
  hPtGenAccFDReb->SetDirectory(0);
  hPtRecoPromptReb->Sumw2();
  hPtRecoFDReb->Sumw2();
  hPtGenAccPromptReb->Sumw2();
  hPtGenAccFDReb->Sumw2();
  
  TH1F* hPtEffPrompt = new TH1F("hPtEffPrompt","Prompt Efficiency",nPtBins,PtLims);
  hPtEffPrompt->Divide(hPtRecoPromptReb,hPtGenAccPromptReb,1.,1.,"B");
  TH1F* hPtEffFD = new TH1F("hPtEffFD","FD Efficiency",nPtBins,PtLims);
  hPtEffFD->Divide(hPtRecoFDReb,hPtGenAccFDReb,1.,1.,"B");
  hPtEffPrompt->SetDirectory(0);
  hPtEffFD->SetDirectory(0);

  //________________________________________________________________________________________
  //acceptance calculation
  TH1F* hPtGenLimAccReb = (TH1F*)hPtGenLimAcc->Rebin(nPtBins,"hPtGenLimAccReb",PtLims);
  TH1F* hPtGenAccReb = (TH1F*)hPtGenAcc->Rebin(nPtBins,"hPtGenAccReb",PtLims);
  hPtGenLimAccReb->SetDirectory(0);
  hPtGenAccReb->SetDirectory(0);
  hPtGenLimAccReb->Sumw2();
  hPtGenAccReb->Sumw2();

  TH1F* hPtAcc = new TH1F("hPtAcc","acceptance",nPtBins,PtLims);
  hPtAcc->Divide(hPtGenAccReb,hPtGenLimAccReb,1.,1.,"B");
  hPtAcc->SetDirectory(0);
    
  //_________________________________________________________________________________________
  //eff x acc
  TH1F* hPtEffXAccPrompt = new TH1F("hPtEffXAccPrompt","Prompt efficiency x acceptance",nPtBins,PtLims);
  hPtEffXAccPrompt->Multiply(hPtEffPrompt,hPtAcc,1.,1.);
  hPtEffXAccPrompt->SetDirectory(0);
  
  //_________________________________________________________________________________________
  //plot 
  TCanvas* cEff = new TCanvas("cEff","",1200,900);
  cEff->Clear();
  hPtEffFD->Draw("E");
  hPtEffPrompt->Draw("Esame");

  TCanvas* cAcc = new TCanvas("cAcc","",1200,900);
  cAcc->Clear();
  hPtAcc->Draw("E");

  TCanvas* cEffXAcc = new TCanvas("cEffXAcc","",1200,900);
  cEffXAcc->Clear();
  hPtEffXAccPrompt->Draw("E");
  
  TFile outfile("efficiency_and_acceptance.root","RECREATE");
  hPtEffPrompt->Write();
  hPtEffFD->Write();
  hPtAcc->Write();
  hPtEffXAccPrompt->Write();
  outfile.Close();
}
