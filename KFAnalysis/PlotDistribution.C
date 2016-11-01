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

enum {kDzero,kDplus};

TString axesfile = "axes";
TString cutsfile = "cutset";
Int_t PIDaxnum=13;

void PlotDistribution(TString filename="ScatterPlot",Int_t axis1=11, Int_t axis2=12,TString varname1="#chi^{2}/ndf", TString varname2="Norm. max d0-d0exp",Int_t meson=kDplus, Double_t min1=0, Double_t max1=50, Double_t min2=-10, Double_t max2=10);
void PlotDist2D(THnSparseF* AllSparse, THnSparseF* PromptSparse, THnSparseF* FDSparse,Int_t axis1,Int_t axis2,TString varname1, TString varname2, Double_t min1, Double_t max1,Double_t min2, Double_t max2, TString filename);
void PlotDist1D(THnSparseF* AllSparse, THnSparseF* PromptSparse, THnSparseF* FDSparse,Int_t axis,TString varname, Double_t min1, Double_t max1, TString filename);
void ResetAxes(THnSparseF* sparse);
void ReadAxesNum(TString Filename, vector<string> &axesanmes, vector<int> &axesno);
void ReadSet(TString Filename, vector<string> &varnames, vector<double> &cutset);
void UsePID(THnSparseF* sparse);
void SetPtRange(THnSparseF* sparse, Int_t iPt, Int_t iAxis, Double_t &ptmin, Double_t &ptmax);
void ApplyTopologicalCuts(THnSparseF* sparse, Int_t iPt, Bool_t chicut);

void PlotDistribution(TString filename, Int_t axis1, Int_t axis2,TString varname1, TString varname2,Int_t meson, Double_t min1, Double_t max1, Double_t min2, Double_t max2) {

  TString mesonname;
  if(meson==kDzero) mesonname = "Dzero";
  else mesonname = "Dplus";
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetLegendBorderSize(0);

  //INPUT FILES________________________________________________________________
  TString MCfilename;
  TString MClistname;

  if(meson==kDzero) {
    MCfilename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pPb_LHC13/MC/fgrosa_Dzero_KF_pPbMC.root";
    //MCfilename = "/home/fabrizio/fgrosa_Dzero_KF_MC.root";
    MClistname = "coutputDzeroKF";
  }
  else {
    MCfilename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pPb_LHC13/MC/fgrosa_Dplus_KF_pPbMC.root";
    //MCfilename = "/home/fabrizio/fgrosa_Dplus_KF_MC.root";
    MClistname = "coutputDplusKF";
  }
  TFile infile(MCfilename.Data(),"UPDATE");
  TList* list = (TList*)infile.Get(MClistname.Data());  

  THnSparseF *AllSparse = (THnSparseF*)list->FindObject("fSparseAll");
  THnSparseF *PromptSparse = (THnSparseF*)list->FindObject("fSparsePrompt");
  THnSparseF *FDSparse = (THnSparseF*)list->FindObject("fSparseFD");

  infile.Close();
    
  //PLOT DISTRIBUTION_______________________________________________________
  if(axis1>=0 && axis2>=0)
    PlotDist2D(AllSparse,PromptSparse,FDSparse,axis1,axis2,varname1,varname2,min1,max1,min2,max2,filename);
  else 
    PlotDist1D(AllSparse,PromptSparse,FDSparse,axis1,varname1,min1,max1,filename);
}

void PlotDist1D(THnSparseF* AllSparse, THnSparseF* PromptSparse, THnSparseF* FDSparse, Int_t axis,TString varname, Double_t min1, Double_t max1, TString filename) {

  Double_t ptmin, ptmax;
  
  ResetAxes(AllSparse);
  ResetAxes(PromptSparse);
  ResetAxes(FDSparse);
  ApplyTopologicalCuts(AllSparse,0,kTRUE);
  ApplyTopologicalCuts(PromptSparse,0,kTRUE);
  ApplyTopologicalCuts(FDSparse,0,kTRUE);
  SetPtRange(AllSparse,0,1,ptmin, ptmax);
  SetPtRange(PromptSparse,0,1,ptmin, ptmax);
  SetPtRange(FDSparse,0,1,ptmin, ptmax);
  
  TH1F* hAll = (TH1F*)AllSparse->Projection(axis);
  TH1F* hPrompt = (TH1F*)PromptSparse->Projection(axis);
  TH1F* hFD = (TH1F*)FDSparse->Projection(axis);

  hAll->GetXaxis()->SetRangeUser(min1,max1);
  hPrompt->GetXaxis()->SetRangeUser(min1,max1);
  hFD->GetXaxis()->SetRangeUser(min1,max1);
  
  TCanvas* cPrompt = new TCanvas("cPrompt","Prompt",800,800);
  cPrompt->SetLogy();
  hPrompt->Draw();
  hPrompt->GetXaxis()->SetTitle(varname.Data());
  hPrompt->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptmin,ptmax));
  TCanvas* cFD = new TCanvas("cFD","Feed-down",800,800);
  cFD->SetLogy();
  hFD->Draw();
  hFD->GetXaxis()->SetTitle(varname.Data());
  hFD->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptmin,ptmax));
  TCanvas* cAll = new TCanvas("cAll","All",800,800);
  cAll->SetLogy();
  hAll->Draw();
  hAll->GetXaxis()->SetTitle(varname.Data());
  hAll->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptmin,ptmax));

  cPrompt->SaveAs(Form("%sPrompt_%d_Pt_%0.f-%0.f.eps",filename.Data(),axis,ptmin,ptmax));
  cFD->SaveAs(Form("%sFD_%d_Pt_%0.f-%0.f.eps",filename.Data(),axis,ptmin,ptmax));
  cAll->SaveAs(Form("%sAll_%d_Pt_%0.f-%0.f.eps",filename.Data(),axis,ptmin,ptmax));

}

void PlotDist2D(THnSparseF* AllSparse, THnSparseF* PromptSparse, THnSparseF* FDSparse, Int_t axis1, Int_t axis2,TString varname1, TString varname2, Double_t min1, Double_t max1, Double_t min2, Double_t max2, TString filename) {

  Double_t ptmin, ptmax;

  ResetAxes(AllSparse);
  ResetAxes(PromptSparse);
  ResetAxes(FDSparse);
  ApplyTopologicalCuts(AllSparse,0,kTRUE);
  ApplyTopologicalCuts(PromptSparse,0,kTRUE);
  ApplyTopologicalCuts(FDSparse,0,kTRUE);
  SetPtRange(AllSparse,0,1,ptmin,ptmax);
  SetPtRange(PromptSparse,0,1,ptmin,ptmax);
  SetPtRange(FDSparse,0,1,ptmin,ptmax);
  
  TH2F* h2DAll = (TH2F*)AllSparse->Projection(axis1,axis2);
  TH2F* h2DPrompt = (TH2F*)PromptSparse->Projection(axis1,axis2);
  TH2F* h2DFD = (TH2F*)FDSparse->Projection(axis1,axis2);

  h2DAll->GetXaxis()->SetRangeUser(min1,max1);
  h2DAll->GetYaxis()->SetRangeUser(min2,max2);
  h2DAll->GetXaxis()->SetTitle(varname1.Data());
  h2DAll->GetYaxis()->SetTitle(varname2.Data());
  h2DPrompt->GetXaxis()->SetRangeUser(min1,max1);
  h2DPrompt->GetYaxis()->SetRangeUser(min2,max2);
  h2DPrompt->GetXaxis()->SetTitle(varname1.Data());
  h2DPrompt->GetYaxis()->SetTitle(varname2.Data());
  h2DFD->GetXaxis()->SetRangeUser(min1,max1);
  h2DFD->GetYaxis()->SetRangeUser(min2,max2);
  h2DFD->GetXaxis()->SetTitle(varname1.Data());
  h2DFD->GetYaxis()->SetTitle(varname2.Data());
 
  TCanvas* cPrompt = new TCanvas("cPrompt","Prompt",800,800);
  cPrompt->SetLogz();
  h2DPrompt->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptmin,ptmax));
  h2DPrompt->Draw("colz");
  TCanvas* cFD = new TCanvas("cFD","Feed-down",800,800);
  cFD->SetLogz();
  h2DFD->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptmin,ptmax));
  h2DFD->Draw("colz");
  TCanvas* cAll = new TCanvas("cAll","All",800,800);
  cAll->SetLogz();
  h2DAll->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptmin,ptmax));
  h2DAll->Draw("colz");

  cPrompt->SaveAs(Form("%sPrompt_%d-%d_Pt_%0.f-%0.f.eps",filename.Data(),axis1,axis2,ptmin,ptmax));
  cFD->SaveAs(Form("%sFD_%d-%d_Pt_%0.f-%0.f.eps",filename.Data(),axis1,axis2,ptmin,ptmax));
  cAll->SaveAs(Form("%sAll_%d-%d_Pt_%0.f-%0.f.eps",filename.Data(),axis1,axis2,ptmin,ptmax));
  
}

void UsePID(THnSparseF* sparse) {
  TAxis* PIDax = (TAxis*)sparse->GetAxis(PIDaxnum);
  Double_t binmin = PIDax->FindBin(1*1.0001);
  Double_t binmax = PIDax->FindBin(4*0.9999);
  sparse->GetAxis(PIDaxnum)->SetRange(binmin,binmax);
}

void ApplyTopologicalCuts(THnSparseF* sparse, Int_t iPt, Bool_t chicut) {
  vector<string> axesnames;
  vector<int> axesno;
  vector<string> cutvarnames;
  vector<double> cutset;
  ReadAxesNum(axesfile,axesnames,axesno);
  ReadSet(cutsfile,cutvarnames,cutset);
  const Int_t nVars = axesno.size();
  Int_t chiax = 0;
  for(Int_t iVar=0; iVar<nVars; iVar++) {
    chiax = axesno[iVar];
    if(axesnames[iVar]=="chi") break;
  }

  for(Int_t iVar=1; iVar<nVars; iVar++) {
    if(!chicut && axesno[iVar]==chiax) {continue;} ///if chicut=kFALSE don't apply cut on chi square
    TAxis* ax = (TAxis*)sparse->GetAxis(axesno[iVar]);
    Double_t binmin = ax->FindBin(cutset[(2*iVar)+iPt*nVars]+0.0001);
    Double_t binmax = ax->FindBin(cutset[(2*iVar+1)+iPt*nVars]-0.0001);
    sparse->GetAxis(axesno[iVar])->SetRange(binmin,binmax);
  }
}

void ResetAxes(THnSparseF* sparse) {
  Int_t nAxes = sparse->GetNdimensions();
  for(Int_t iAxis=0; iAxis<nAxes; iAxis++) {
    sparse->GetAxis(iAxis)->SetRange(-1,-1);
  }
}

void SetPtRange(THnSparseF* sparse, Int_t iPt, Int_t iAxis, Double_t &ptmin, Double_t &ptmax) {
 vector<string> axesnames;
  vector<int> axesno;
  vector<string> cutvarnames;
  vector<double> cutset;
  ReadAxesNum(axesfile,axesnames,axesno);
  ReadSet(cutsfile,cutvarnames,cutset);
  const Int_t nVars = axesno.size();

  ptmin=cutset[iPt*2*nVars];
  ptmax=cutset[1+iPt*2*nVars];
  
  TAxis* ptax = (TAxis*)sparse->GetAxis(iAxis);
  Int_t binmin = ptax->FindBin(ptmin+0.0001);
  Int_t binmax = ptax->FindBin(ptmax-0.0001);
  sparse->GetAxis(iAxis)->SetRange(binmin,binmax);
}

void ReadAxesNum(TString FileName, vector<string> &axesnames, vector<int> &axesno) {  
  ifstream inAxes(FileName.Data());
  if(!inAxes){
    cerr<<"Il file "<<FileName.Data() <<" non esiste "<<endl;
    return;
  }

  string axisname;
  string axestring;

  getline(inAxes,axestring);
  stringstream ss(axestring);
  
  while(ss >> axisname)
    axesnames.push_back(axisname);

  Int_t num;
  while(inAxes>>num)
    axesno.push_back(num);
  
  inAxes.close();
}

void ReadSet(TString FileName, vector<string> &varnames, vector<double> &cutset) {
  ifstream inSet(FileName.Data());
  if(!inSet){
    cerr<<"Il file "<<FileName.Data() <<" non esiste "<<endl;
    return;
  }

  string varname;
  string varstring;

  getline(inSet,varstring);
  stringstream ss(varstring);
  
  while(ss >> varname)
    varnames.push_back(varname);

  Double_t num;
  while(inSet>>num)
    cutset.push_back(num);
  
  inSet.close();
}

