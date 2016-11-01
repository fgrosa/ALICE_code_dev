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

//macro for the cross section calculation of D+ mesons
//author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino

//****************************************************//
//                                                    //
//    Main Function: CalculateDplusCrossSection       //
//                                                    //
//****************************************************//

//________________________________________________________________________________________________________________
Int_t CalculateDplusCrossSection(TString rawfilename="RawYields.root",
                                 TString efffilename="Efficiency_times_Acceptance.root",
                                 TString fpromptfilename="fprompt_unbinned_sigmafree.root",
                                 Double_t sigma=2.09/*barn*/, Double_t BR=0.0913,
                                 TString outfilename="DplusCrossSection.root") {

  //________________________________________________________________________________________________________________
  //input files
  TFile *rawfile = TFile::Open(rawfilename,"READ");
  TH1F* hRawYields = 0x0;
  TH1F* hEv = 0x0;
  if(!rawfile) return 1;
  else {
    hRawYields = (TH1F*)rawfile->Get("hRawYields");
    hEv = (TH1F*)rawfile->Get("hEv");
    if(hRawYields || hEv) {
      hRawYields->SetDirectory(0);
      hEv->SetDirectory(0);
    }
    else {
      cerr << "Histo hRawYields or hEv does not exists. Exit" << endl;
      return 2;
    }
  }
  rawfile->Close();
  
  TFile *efffile = TFile::Open(efffilename,"READ");
  TH1F* hEffPrompt = 0x0;
  if(!efffile) return 1;
  else {
    hEffPrompt = (TH1F*)efffile->Get("hEffAccPrompt");
    if(hEffPrompt) {
      hEffPrompt->SetDirectory(0);
    }
    else {
      cerr << "Histo hEffAccPrompt does not exists. Exit" << endl;
      return 2;
    }
  }
  efffile->Close();
  
  TFile *fpromptfile = TFile::Open(fpromptfilename,"READ");
  TH1F* hFprompt = 0x0;
  if(!fpromptfile) return 1;
  else {
    hFprompt = (TH1F*)fpromptfile->Get("hFrac");
    if(hFprompt) {
      hFprompt->SetDirectory(0);
    }
    else {
      cerr << "Histo hFrac does not exists. Exit" << endl;
      return 2;
    }
  }
  fpromptfile->Close();
  
  //________________________________________________________________________________________________________________
  //cross section calculation
  TH1F* hPromptCrossSection = (TH1F*)hFprompt->Clone();
  for(Int_t iPt=0; iPt<hPromptCrossSection->GetNbinsX(); iPt++) {
    Double_t cross = hRawYields->GetBinContent(iPt+2)*hFprompt->GetBinContent(iPt+1)*(sigma*1000000)/(2*hEffPrompt->GetBinContent(iPt+2)*hRawYields->GetBinWidth(iPt+2)*BR*hEv->GetBinContent(1));
    Double_t crosserr = TMath::Sqrt(hRawYields->GetBinError(iPt+2)/hRawYields->GetBinContent(iPt+2)*hRawYields->GetBinError(iPt+2)/hRawYields->GetBinContent(iPt+2)+hEffPrompt->GetBinError(iPt+2)/hEffPrompt->GetBinContent(iPt+2)*hEffPrompt->GetBinError(iPt+2)/hEffPrompt->GetBinContent(iPt+2)+hFprompt->GetBinError(iPt+1)/hFprompt->GetBinContent(iPt+1)*hFprompt->GetBinError(iPt+1)/hFprompt->GetBinContent(iPt+1))*cross;
    
    hPromptCrossSection->SetBinContent(iPt+1,cross);
    hPromptCrossSection->SetBinError(iPt+1,crosserr);
  }
  
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.18);
  
  TCanvas* cCrossSection = new TCanvas("cCrossSection","",10,10,800,800);
  cCrossSection->SetLogy();
  hPromptCrossSection->SetName("hPromptCrossSection");
  hPromptCrossSection->GetYaxis()->SetRangeUser(hPromptCrossSection->GetBinContent(hPromptCrossSection->GetNbinsX())*0.5,hPromptCrossSection->GetBinContent(1)*1.5);
  hPromptCrossSection->GetYaxis()->SetTitle("#frac{d#sigma_{prompt}}{d#it{p}_{T}} (#mub c/GeV)");
  hPromptCrossSection->GetYaxis()->SetTitleOffset(1.6);
  hPromptCrossSection->Draw();

  //________________________________________________________________________________________________________________
  //output file
  TFile outfile(outfilename,"RECREATE");
  hPromptCrossSection->Write();
  outfile.Close();
  cout << "\nFile " << outfilename << " saved.\n" << endl;
  
  return 0;
  
}
