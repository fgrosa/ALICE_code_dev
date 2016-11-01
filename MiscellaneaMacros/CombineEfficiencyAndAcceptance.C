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
#include <TColor.h>

#endif

//macro for combining the efficiency and acceptance of D+ mesons
//author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino

Int_t CombineEfficiencyAndAcceptance(TString efffilename="Efficiency.root",
                                     TString accfilename="$HOME/ALICE_WORK/Files/Acceptance/Acceptance_Toy_DplusKpipi_yfidPtDep_etaDau09_ptDau100_FONLL5ptshape.root",
                                     TString outfilename="Efficiency_times_Acceptance.root") {
  
  
  TFile *efffile = TFile::Open(efffilename,"READ");
  TH1F* hEffPrompt = 0x0;
  TH1F* hEffFD = 0x0;
  if(!efffile) return 1;
  else {
    hEffPrompt = (TH1F*)efffile->Get("hEffPrompt");
    hEffFD = (TH1F*)efffile->Get("hEffFD");
    if(!hEffPrompt || !hEffFD) {
      cerr << "The efficiency histos in " << efffilename << " are missing! Exit." << endl;
      return 2;
    }
    else {
      hEffPrompt->SetDirectory(0);
      hEffFD->SetDirectory(0);
    }
  }
  efffile->Close();
  
  TFile *accfile = TFile::Open(accfilename,"READ");
  TH1F* hGenLimAcc = 0x0;
  TH1F* hGenAcc = 0x0;
  if(!accfile) return 1;
  else {
    hGenLimAcc = (TH1F*)accfile->Get("hPtGenLimAcc");
    hGenAcc = (TH1F*)accfile->Get("hPtGenAcc");
    if(!hGenLimAcc || !hGenAcc) {
      cerr << "The acceptance histos in " << accfilename << " are missing! Exit." << endl;
      return 2;
    }
    else {
      hGenLimAcc->SetDirectory(0);
      hGenAcc->SetDirectory(0);
      hGenLimAcc->Sumw2();
      hGenAcc->Sumw2();
    }
  }
  accfile->Close();
  
  //________________________________________________________________________________________
  //get pt bins
  TAxis* ptaxis = (TAxis*)hEffPrompt->GetXaxis();
  const Int_t nPtBins = ptaxis->GetNbins();
  TArrayD* ptarray = (TArrayD*)ptaxis->GetXbins();
  Double_t* PtLims = (Double_t*)ptarray->GetArray();
  
  //________________________________________________________________________________________
  //acceptance calculation
  TH1F* hGenLimAccReb = (TH1F*)hGenLimAcc->Rebin(nPtBins,"hGenLimAccReb",PtLims);
  TH1F* hGenAccReb = (TH1F*)hGenAcc->Rebin(nPtBins,"hGenAccReb",PtLims);
  hGenLimAccReb->SetDirectory(0);
  hGenAccReb->SetDirectory(0);
  hGenLimAccReb->Sumw2();
  hGenAccReb->Sumw2();
  
  TH1F* hAcc = new TH1F("hAcc","",nPtBins,PtLims);
  hAcc->Divide(hGenAccReb,hGenLimAccReb,1.,1.,"B");
  hAcc->SetDirectory(0);
  
  //_________________________________________________________________________________________
  //eff x acc
  TH1F* hEffAccPrompt = (TH1F*)hEffPrompt->Clone();
  hEffAccPrompt->Multiply(hEffPrompt,hAcc,1.,1.);
  hEffAccPrompt->GetYaxis()->SetTitle("efficiency x acceptance");
  hEffAccPrompt->SetDirectory(0);
  TH1F* hEffAccFD = (TH1F*)hEffFD->Clone();
  hEffAccFD->Multiply(hEffFD,hAcc,1.,1.);
  hEffAccFD->GetYaxis()->SetTitle("efficiency x acceptance");
  hEffAccFD->GetYaxis()->SetTitleOffset(1.4);
  hEffAccFD->SetDirectory(0);
  hEffAccPrompt->SetName("hEffAccPrompt");
  hEffAccFD->SetName("hEffAccFD");
  
  //_________________________________________________________________________________________
  //draw
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.12);
  TCanvas *cEffAcc = new TCanvas("cEffAcc","",10,10,600,600);
  cEffAcc->SetLogy();
  hEffAccFD->Draw();
  hEffAccPrompt->Draw("same");
  
  //_________________________________________________________________________________________
  //output file
  TFile outfile(outfilename,"RECREATE");
  hEffAccPrompt->Write();
  hEffAccFD->Write();
  outfile.Close();
  
  return 0;
  
}
