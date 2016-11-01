#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <THnSparse.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TColor.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include <TLine.h>
#include <TDatabasePDG.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>

#include "AliDplusCharmFractionIPfitter.h"

#endif

//Macro that performs the impact parameter fits for the measurement of the fraction of prompt D+
//Author: Fabrizio Grosa, INFN Turin grosa@to.infn.it

//***************************************//
//                                       //
//      Main Function: FitImpPar()       //
//                                       //
//***************************************//

//_____________________________________________________________________________________________
//GLOBAL VARIABLES

//PtBins of the analysis
const Int_t nPtBins = 9;
const Int_t nPtLims = nPtBins+1;
const Double_t PtLims[nPtLims] = {2,3,4,5,6,8,10,12,16,24};
//IP range limits
Double_t d0limFD[nPtBins] = {200,200,300,300,300,300,300,400};
Double_t d0lim[nPtBins] = {200,200,300,300,300,300,300,400};
Double_t d0limPrompt[nPtBins] = {300,300,300,300,300,300,300,300};

//input file names
const TString infileMCname="$HOME/ALICE_WORK/Files/Trains/Run1/LHC13/AnalysisResultspPbMC.root";
const TString dirMCname="PWG3_D2H_InvMassDplus";
const TString listMCname="coutputDplus_ImpParpPbMC0100";
const TString infileDataname="$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/AnalysisResults_3050_central_imppar.root";//"$HOME/ALICE_WORK/Files/Trains/Run1/LHC13/AnalysisResultspPbData.root";
const TString dirDataname=dirMCname;
const TString listDataname="coutputDplus_3050_CentralCuts_kINT73050";//"coutputDplus_ImpParpPbData0100";

enum {kBinned,kUnbinned};

//_____________________________________________________________________________________________
//FUNCTION PROTOTYPES
Int_t FitImpPar(Int_t method=kUnbinned,
                TString fitoption="RLEM0",
                Bool_t isSigmaFixed=kFALSE,
                Double_t nSigmas=1.5,
                Bool_t isVarBinning=kFALSE,
                Bool_t isBkgSub=kFALSE,
                Bool_t PIDcut=kTRUE,
                Int_t sovert=AliDplusCharmFractionIPfitter::kCentralValue,
                Int_t SBregion=AliDplusCharmFractionIPfitter::kBoth,
                Int_t FDtype=AliDplusCharmFractionIPfitter::kConvolution,
                Int_t bkgtype=AliDplusCharmFractionIPfitter::kDoubleGaussExpo,
                Int_t impparreb=4, Int_t massreb=5, Int_t SBlow=4, Int_t SBhigh=15,
                Double_t d0cut=80);

Int_t LoadMCSparses(THnSparseF *&promptsparse, THnSparseF *&trueFDsparse, THnSparseF *&recoFDsparse, THnSparseF *&bkgsparse);
Int_t LoadDataSparsesAndTree(Int_t method, THnSparseF *&datasparse, TNtuple *&tree);
void SetStyle();

//_____________________________________________________________________________________________
//IMPACT PARAMETER FIT FUNCTION
Int_t FitImpPar(Int_t method,TString fitoption,Bool_t isSigmaFixed,Double_t nSigmas,Bool_t isVarBinning,Bool_t isBkgSub,Bool_t PIDcut,Int_t sovert,Int_t SBregion,Int_t FDtype,Int_t bkgtype,Int_t impparreb, Int_t massreb, Int_t SBlow, Int_t SBhigh,Double_t d0cut) {

  //___________________________________________________________________________________________
  //input files
  
  THnSparseF* hMassPtImpPrompt=0x0;
  THnSparseF* hMassPtImpRecoFD=0x0;
  THnSparseF* hMassPtImpTrueFD=0x0;
  THnSparseF* hMassPtImpParBkg=0x0;
  Int_t loadMC=LoadMCSparses(hMassPtImpPrompt,hMassPtImpRecoFD,hMassPtImpTrueFD,hMassPtImpParBkg);
  if(loadMC>0) {return 1;}

  THnSparseF* hMassPtImpParAll=0x0;
  TNtuple* dataTree=0x0;
  Int_t loadData=LoadDataSparsesAndTree(method,hMassPtImpParAll,dataTree);
  if(loadData>0) {return 2;}

  //____________________________________________________________________________________________
  //analysis

  TH1F* hFrac = new TH1F("hFrac","",nPtBins,PtLims);
  TH1F* hFracd0Cut = new TH1F("hFrac","",nPtBins,PtLims);
  TH1F* hSigma = new TH1F("hSigma","",nPtBins,PtLims);
  hFrac->SetDirectory(0);
  hFracd0Cut->SetDirectory(0);
  hSigma->SetDirectory(0);
  
  AliDplusCharmFractionIPfitter *ImpParFitter = new AliDplusCharmFractionIPfitter();
  ImpParFitter->SetDataSparse(hMassPtImpParAll);
  ImpParFitter->SetMCPromptSparse(hMassPtImpPrompt);
  ImpParFitter->SetMCTrueFDSparse(hMassPtImpTrueFD);
  ImpParFitter->SetMCRecoFDSparse(hMassPtImpRecoFD);
  if(method==kUnbinned && dataTree) ImpParFitter->SetDataTree(dataTree);
  
  ImpParFitter->SetSideBandsRegion(SBregion);
  ImpParFitter->SetBkgFunction(bkgtype);
  ImpParFitter->SetFDFunction(FDtype);
  ImpParFitter->SetPID(PIDcut);
  
  Double_t initparprompt[5] = {0.9,0.,35,300,1};
  Double_t initparFD[5] = {0.5,0.,100,10,1};
  Double_t initparBkg[10] = {0.3,-30,50,50,0.5,50,50,50,0.5,1};

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {

    ImpParFitter->SetRebinImpParHistos(impparreb);
    ImpParFitter->SetPtLims(PtLims[iPt],PtLims[iPt+1]);
    ImpParFitter->SetNSigmas(nSigmas);
    ImpParFitter->SetNSigmaSBLimits(SBlow,SBhigh);
    ImpParFitter->GetSignal(massreb,0,0,1.68,2.05,sovert);
    ImpParFitter->SetFitOptions(fitoption);
    ImpParFitter->FixSigmaPromptFromMC(isSigmaFixed);
    ImpParFitter->SetRelativeLimitSigmaPrompt(0.3);
    ImpParFitter->SetLimitsPromptFrac(0.2,1.5);
    ImpParFitter->SetInitialParameters(initparprompt,initparFD,initparBkg);
    ImpParFitter->PrefitStep(-d0limPrompt[iPt],d0limPrompt[iPt],-d0limFD[iPt],d0limFD[iPt],-1000,1000);
    if(method==kUnbinned)
      ImpParFitter->FitTree(-d0lim[iPt],d0lim[iPt],kTRUE);
    else if(method==kBinned) {
      ImpParFitter->SetVariableBinningHisto(isVarBinning,5);
      ImpParFitter->SetBkgSubtraction(isBkgSub);
      ImpParFitter->FitHisto(-d0lim[iPt],d0lim[iPt],kTRUE);
    }
    else
      cerr << "only binned or unbinned fits are supported" << endl;

    Double_t promptfracd0cut;
    Double_t promptfracd0cuterr;
    ImpParFitter->GetPromptFractionWithIPCut(d0cut,promptfracd0cut,promptfracd0cuterr);
    Double_t promptfrac = ImpParFitter->GetPromptFraction();
    Double_t promptfracerr = ImpParFitter->GetPromptFractionErr();
    Double_t promptsigma = ImpParFitter->GetPromptSigma();
    Double_t promptsigmaerr = ImpParFitter->GetPromptSigmaErr();
    if(isSigmaFixed) {
      promptsigma = ImpParFitter->GetPromptSigmaMC();
      promptsigmaerr = ImpParFitter->GetPromptSigmaMCErr();
    }
    hFrac->SetBinContent(iPt+1,promptfrac);
    hFrac->SetBinError(iPt+1,promptfracerr);
    hFracd0Cut->SetBinContent(iPt+1,promptfracd0cut);
    hFracd0Cut->SetBinError(iPt+1,promptfracd0cuterr);
    hSigma->SetBinContent(iPt+1,promptsigma);
    hSigma->SetBinError(iPt+1,promptsigmaerr);
  }

  //Set draw style
  SetStyle();
  
  TCanvas* cFrac = new TCanvas("cFrac","",800,600);
  cFrac->Clear();
  hFrac->SetLineWidth(2);
  hFrac->SetMarkerSize(1.5);
  hFrac->SetMarkerStyle(20);
  hFrac->SetLineColor(kBlack);
  hFrac->SetMarkerColor(kBlack);
  hFrac->GetYaxis()->SetRangeUser(0,1.2);
  hFrac->GetYaxis()->SetTitle("#it{f}_{prompt}");
  hFrac->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hFracd0Cut->GetYaxis()->SetTitle("#it{f}_{prompt}");
  hFracd0Cut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hFrac->Draw("E1");
  hFracd0Cut->SetLineWidth(2);
  hFracd0Cut->SetMarkerSize(1.5);
  hFracd0Cut->SetMarkerStyle(20);
  hFracd0Cut->SetLineColor(kRed);
  hFracd0Cut->SetMarkerColor(kRed);
  hFracd0Cut->GetYaxis()->SetRangeUser(0,1.2);
  hFracd0Cut->Draw("E1same");
  TLegend *leg = new TLegend(0.6,0.3,0.89,0.5);
  leg->SetTextSize(0.05);
  leg->AddEntry(hFrac,"w/o d_{0} cut","lpe");
  leg->AddEntry(hFracd0Cut,Form("|d_{0}| < %0.f #mum",d0cut),"lpe");
  leg->Draw("same");
  
  TCanvas* cSigma = new TCanvas("cSigma","",800,600);
  cSigma->Clear();
  hSigma->SetLineWidth(2);
  hSigma->SetMarkerSize(1.5);
  hSigma->SetMarkerStyle(20);
  hSigma->SetLineColor(kBlack);
  hSigma->SetMarkerColor(kBlack);
  hSigma->GetYaxis()->SetTitle("#sigma_{prompt} (#mum)");
  hSigma->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hSigma->Draw("E1");

  TString sigmaname="sigmafree";
  if(isSigmaFixed) {sigmaname="sigmafixed";}
  TString methodname="unbinned";
  if(method==kBinned) {methodname="binned";}
  
  TFile outfile(Form("fprompt_%s_%s.root",methodname.Data(),sigmaname.Data()),"RECREATE");
  hFrac->Write();
  hSigma->Write();
  outfile.Close();
  TFile outfile2(Form("fprompt_%s_%s_d0cut.root",methodname.Data(),sigmaname.Data()),"RECREATE");
  hFracd0Cut->Write();
  hSigma->Write();
  outfile2.Close();
  
  cFrac->SaveAs(Form("fprompt_%s_%s.pdf",methodname.Data(),sigmaname.Data()));
  cSigma->SaveAs(Form("sigmaprompt_%s_%s.pdf",methodname.Data(),sigmaname.Data()));
  
  return 0;
}

//_____________________________________________________________________________________________
//LOAD MC SPARSES FUNCTION
Int_t LoadMCSparses(THnSparseF *&promptsparse, THnSparseF *&trueFDsparse, THnSparseF *&recoFDsparse, THnSparseF *&bkgsparse) {
  
  cout<<"Opening MC file " <<infileMCname<< "..." << endl;
  TFile* infileMC = TFile::Open(infileMCname.Data(),"READ");
  TDirectoryFile* dirMC=0x0;
  TList* listMC=0x0;
  THnSparseF* hMassPtImpPrompt=0x0;
  THnSparseF* hMassPtImpRecoFD=0x0;
  THnSparseF* hMassPtImpTrueFD=0x0;
  
  if(infileMC) {dirMC=(TDirectoryFile*)infileMC->Get(dirMCname.Data()); cout << "MC file opened!" << endl;}
  else {cerr << "Error: File " << infileMCname << " not found. Exit." << endl; return 1;}
  if(dirMC) listMC=(TList*)dirMC->Get(listMCname.Data());
  else {cerr << "Error: Wrong TDirectoryFile name " << dirMCname << ". Exit." << endl; return 2;}
  if(listMC) {
    promptsparse=(THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
    recoFDsparse=(THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
    trueFDsparse=(THnSparseF*)listMC->FindObject("hMassPtImpParTrueBfeed");
    bkgsparse=(THnSparseF*)listMC->FindObject("hMassPtImpParBkg");
    cout << "Sparses got!" << endl;
  }
  else {cerr << "Error: Wrong TList name " << listMCname << ". Exit." << endl; return 3;}
  if(!promptsparse) {cerr << "Error: No MC sparse for prompt D+. Check if the name hMassPtImpParPrompt is right!" << endl; return 4;}
  if(!recoFDsparse) {cerr << "Error: No MC sparse for feed-down D+. Check if the name hMassPtImpParBfeed is right!" << endl; return 5;}
  if(!trueFDsparse) {cerr << "Error: No MC sparse for true feed-down D+. Check if the name hMassPtImpParTrueBfeed is right!" << endl; return 6;}
  if(!bkgsparse) {cout << "Warning: No MC sparse for the background. Check if the name hMassPtImpParBkg is right!" << endl;}
  
  infileMC->Close();
  cout<<"MC file closed."<< endl;

  return 0;
}

//_____________________________________________________________________________________________
//LOAD DATA SPARSE AND TREE FUNCTION
Int_t LoadDataSparsesAndTree(Int_t method, THnSparseF *&datasparse, TNtuple *&tree) {

  cout<<"Opening data file " <<infileMCname<< "..." << endl;
  TFile* infileData = TFile::Open(infileDataname.Data(),"READ");
  TDirectoryFile* dirData=0x0;
  TList* listData=0x0;
  THnSparseF* hMassPtImpParAll=0x0;
  TNtuple* dataTree=0x0;
  
  if(infileData) {dirData=(TDirectoryFile*)infileData->Get(dirDataname.Data()); cout << "Data file opened!" << endl;}
  else {cerr << "Error: File " << infileDataname << " not found. Exit." << endl; return 1;}
  if(dirData) listData=(TList*)dirData->Get(listDataname.Data());
  else {cerr << "Error: Wrong TDirectoryFile name " << dirDataname << ". Exit." << endl; return 2;}
  if(listData) {
    datasparse=(THnSparseF*)listData->FindObject("hMassPtImpParAll");
    cout << "Sparse got!" << endl;
    if(method==kUnbinned) {tree=(TNtuple*)dirData->Get("fNtupleDplus"); cout << "TNtuple got!" << endl;}
  }
  else {cerr << "Error: Wrong TList name " << listDataname << ". Exit." << endl; return 3;}
  if(!datasparse && method==kUnbinned) {cout << "Warning: No data sparse. Check if the name hMassPtImpParAll is right!" << endl;}
  if(!datasparse && method==kBinned) {cout << "Error: No data sparse. Check if the name hMassPtImpParAll is right!" << endl; return 4;}
  if(!tree && method==kUnbinned) {cout << "Error: No data tree. Check if the name fNtupleDplus is right!" << endl; return 5;}
  
  infileData->Close();
  cout<<"Data file closed."<< endl;

  return 0;
}

//_____________________________________________________________________________________________
//DRAW STYLE FUNCTION
void SetStyle() {

  cout << "Setting drawing style!" << endl;
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);
}
