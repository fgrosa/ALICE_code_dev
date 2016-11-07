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
#include <AliNormalizationCounter.h>

#include "AliHFMassFitter.h"
#include "AliHFMassFitterVAR.h"

#endif

//macro for the standard analysis of D+ mesons
//author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino

//**************************************//
//                                      //
//    Main Function: AnalyseDplus       //
//                                      //
//**************************************//

//________________________________________________________________________________________________________________
//global variables
const TString cent = "3050";
const TString axesfile = "axes";
const TString cutsfile = Form("Cent%s/cutset_topocut2.txt",cent.Data());

const Int_t nDataFiles = 2;
const TString datafilename[nDataFiles] = {Form("$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/HR_bunch4/AnalysisResults_%s.root",cent.Data())
                                         ,Form("$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/HR_bunch1-3/AnalysisResults_%s.root",cent.Data())};
                                          //,Form("$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/HR_bunch1-3/AnalysisResults_%s_part2.root",cent.Data())
                                          //,Form("$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/HR_bunch1-3/AnalysisResults_%s_part3.root",cent.Data())};

//const TString datafilename[nDataFiles] = {Form("$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/HR_bunch4/AnalysisResults_%s.root",cent.Data())};
const TString datadirname = "PWG3_D2H_InvMassDplus";
const TString datalistname = Form("coutputDplus%s_kINT7%s",cent.Data(),cent.Data());

const TString MCfilename = "";
const TString MCdirname = "PWG3_D2H_InvMassDplus";
const TString MClistname = "";

const Int_t rebin[11] = {5,5,5,5,5,5,5,8,8,6,6};
const Int_t sgnfcn = AliHFMassFitter::kGaus;
const Int_t bkgfcn = AliHFMassFitter::kExpo;

//________________________________________________________________________________________________________________
//functions prototypes
Int_t AnalyseDplus(Bool_t fixmean=kFALSE, Bool_t fixsigma=kFALSE, TString outfilerawname=Form("Cent%s/RawYields_%s_d0cut_bunch1-3.root",cent.Data(),cent.Data()), TString outfileeffname=Form("Cent%s/Efficiency.root",cent.Data()));
Int_t ReadAxes(TString FileName, vector<string> &axesanmes, vector<int> &axesno);
Int_t ReadSet(TString FileName, vector<string> &varnames, vector<double> &cutset);
THnSparseF* GetSparse(TString filename, TString dirname, TString listname, TString sparsename);
void ApplyCuts(THnSparseF* sparse, vector<int> axesno, vector<string> axesnames, vector<double> cutset, Int_t iPt);
void ResetAxes(THnSparseF* sparse);
void GetPtLims(Double_t ptlims[], vector<double> cutset, const Int_t nCutVars);
Double_t GetNevents(const TString filename);
TH1F* GetEfficiency(TString effhistoname, TString recosparsename, TString gensparsename, vector<int> axesno, vector<string> axesnames, vector<double> cutset, Int_t nPtBins, Double_t *PtLims, Int_t color=kBlue);
void SetStyle();
void DivideCanvas(TCanvas* c, const Int_t nPtBins);

//________________________________________________________________________________________________________________
Int_t AnalyseDplus(Bool_t fixmean, Bool_t fixsigma, TString outfilerawname, TString outfileeffname) {

  //get thnsparse from files
  THnSparseF** datasparse = new THnSparseF*[nDataFiles];
  
  for(Int_t iDataFile=0; iDataFile<nDataFiles; iDataFile++) {
    datasparse[iDataFile] = (THnSparseF*)GetSparse(datafilename[iDataFile],datadirname,datalistname,"hMassPtCutVarsAll");
    if(!datasparse)
      return 1;
  }
  
  //cuts set
  vector<string> axesnames;
  vector<int> axesno;
  Int_t readaxes = ReadAxes(axesfile,axesnames,axesno);
  if(readaxes==1) return 2;

  vector<string> varnames;
  vector<double> cutset;
  Int_t readset = ReadSet(cutsfile,varnames,cutset);
  if(readset==1) return 3;

  const Int_t nCutVars = axesnames.size(); //pt included
  const Int_t nPtBins = cutset.size()/(2*nCutVars); //2 values for each cut variable (min,max)
  const Int_t nPtLims = nPtBins+1;
  Double_t PtLims[nPtLims];
  GetPtLims(PtLims,cutset,nCutVars);

  //setting drawing style
  SetStyle();

  //invariant-mass distribution fits
  TH1F** hMassPart = new TH1F*[nDataFiles];
  TH1F** hMass = new TH1F*[nPtBins];
  AliHFMassFitter** fitter = new AliHFMassFitter*[nPtBins];
  TH1F* hRawYields = new TH1F("hRawYields","",nPtBins,PtLims);
  TH1F* hRawYieldsSigma = new TH1F("hRawYieldsSigma","",nPtBins,PtLims);
  TH1F* hRawYieldsMean = new TH1F("hRawYieldsMean","",nPtBins,PtLims);
  TH1F* hRawYieldsSignificance = new TH1F("hRawYieldsSignificance","",nPtBins,PtLims);
  TH1F* hRawYieldsSoverB = new TH1F("hRawYieldsSoverB","",nPtBins,PtLims);
  TH1F* hRawYieldsSignal = new TH1F("hRawYieldsSignal","",nPtBins,PtLims);
  TH1F* hRawYieldsBkg = new TH1F("hRawYieldsBkg","",nPtBins,PtLims);
  TH1F* hRawYieldsChiSquare = new TH1F("hRawYieldsChiSquare","",nPtBins,PtLims);
  hRawYields->GetYaxis()->SetTitle("raw yield");
  hRawYields->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hRawYieldsSigma->GetYaxis()->SetTitle("width (GeV/c^{2})");
  hRawYieldsSigma->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hRawYieldsMean->GetYaxis()->SetTitle("mean (GeV/c^{2})");
  hRawYieldsMean->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hRawYieldsSignificance->GetYaxis()->SetTitle("significance (3#sigma)");
  hRawYieldsSignificance->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hRawYieldsSoverB->GetYaxis()->SetTitle("S/B (3#sigma)");
  hRawYieldsSoverB->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hRawYieldsSignal->GetYaxis()->SetTitle("Signal (3#sigma)");
  hRawYieldsSignal->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hRawYieldsBkg->GetYaxis()->SetTitle("Background (3#sigma)");
  hRawYieldsBkg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hRawYieldsChiSquare->GetYaxis()->SetTitle("#chi^{2}");
  hRawYieldsChiSquare->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hRawYields->SetMarkerSize(1.);
  hRawYields->SetMarkerStyle(20);
  hRawYields->SetLineWidth(2);
  hRawYields->SetMarkerColor(kBlack);
  hRawYields->SetLineColor(kBlack);
  hRawYieldsSigma->SetMarkerSize(1.);
  hRawYieldsSigma->SetMarkerStyle(20);
  hRawYieldsSigma->SetLineWidth(2);
  hRawYieldsSigma->SetMarkerColor(kBlack);
  hRawYieldsSigma->SetLineColor(kBlack);
  hRawYieldsMean->SetMarkerSize(1.);
  hRawYieldsMean->SetMarkerStyle(20);
  hRawYieldsMean->SetLineWidth(2);
  hRawYieldsMean->SetMarkerColor(kBlack);
  hRawYieldsMean->SetLineColor(kBlack);
  hRawYieldsSignificance->SetMarkerSize(1.);
  hRawYieldsSignificance->SetMarkerStyle(20);
  hRawYieldsSignificance->SetLineWidth(2);
  hRawYieldsSignificance->SetMarkerColor(kBlack);
  hRawYieldsSignificance->SetLineColor(kBlack);
  hRawYieldsSoverB->SetMarkerSize(1.);
  hRawYieldsSoverB->SetMarkerStyle(20);
  hRawYieldsSoverB->SetLineWidth(2);
  hRawYieldsSoverB->SetMarkerColor(kBlack);
  hRawYieldsSoverB->SetLineColor(kBlack);
  hRawYieldsSignal->SetMarkerSize(1.);
  hRawYieldsSignal->SetMarkerStyle(20);
  hRawYieldsSignal->SetLineWidth(2);
  hRawYieldsSignal->SetMarkerColor(kBlack);
  hRawYieldsSignal->SetLineColor(kBlack);
  hRawYieldsBkg->SetMarkerSize(1.);
  hRawYieldsBkg->SetMarkerStyle(20);
  hRawYieldsBkg->SetLineWidth(2);
  hRawYieldsBkg->SetMarkerColor(kBlack);
  hRawYieldsBkg->SetLineColor(kBlack);
  hRawYieldsChiSquare->SetMarkerSize(1.);
  hRawYieldsChiSquare->SetMarkerStyle(20);
  hRawYieldsChiSquare->SetLineWidth(2);
  hRawYieldsChiSquare->SetMarkerColor(kBlack);
  hRawYieldsChiSquare->SetLineColor(kBlack);

  TCanvas* cMass = new TCanvas("cMass","",10,10,1920,1080);
  DivideCanvas(cMass, nPtBins);

  Double_t massD = TDatabasePDG::Instance()->GetParticle(411)->Mass();

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    for(Int_t iDataFile=0; iDataFile<nDataFiles; iDataFile++) {
      ApplyCuts(datasparse[iDataFile],axesno,axesnames,cutset,iPt);
      hMassPart[iDataFile] = (TH1F*)datasparse[iDataFile]->Projection(0);
    }
    hMass[iPt]=(TH1F*)hMassPart[0]->Clone();
    if(nDataFiles>1) {
      for(Int_t iBin=0; iBin<hMass[iPt]->GetNbinsX(); iBin++) {
        Double_t nentries=0;
        for(Int_t iDataFile=0; iDataFile<nDataFiles; iDataFile++) {
          nentries += hMassPart[iDataFile]->GetBinContent(iBin+1);
        }
        hMass[iPt]->SetBinContent(iBin+1,nentries);
      }
    }
    
    hMass[iPt]->Rebin(rebin[iPt]);
    hMass[iPt]->GetYaxis()->SetTitle(Form("Entries / %0.f MeV/c^{2}",hMass[iPt]->GetBinWidth(5)*1000));
    hMass[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    fitter[iPt] = new AliHFMassFitter(hMass[iPt],1.70,2.03,1,bkgfcn,sgnfcn);
    if(fixmean)
      fitter[iPt]->SetFixGaussianMean(massD);
    else
      fitter[iPt]->SetInitialGaussianMean(massD);
    if(fixsigma)
      fitter[iPt]->SetFixGaussianSigma(0.015);
    else
      fitter[iPt]->SetInitialGaussianSigma(0.015);
    fitter[iPt]->MassFitter(kFALSE);
    hRawYields->SetBinContent(iPt+1,fitter[iPt]->GetRawYield());
    hRawYields->SetBinError(iPt+1,fitter[iPt]->GetRawYieldError());
    hRawYieldsSigma->SetBinContent(iPt+1,fitter[iPt]->GetSigma());
    hRawYieldsSigma->SetBinError(iPt+1,fitter[iPt]->GetSigmaUncertainty());
    hRawYieldsMean->SetBinContent(iPt+1,fitter[iPt]->GetMean());
    hRawYieldsMean->SetBinError(iPt+1,fitter[iPt]->GetMeanUncertainty());
    Double_t signif;
    Double_t signiferr;
    Double_t sgn;
    Double_t sgnerr;
    Double_t bkg;
    Double_t bkgerr;
    fitter[iPt]->Significance(3,signif,signiferr);
    fitter[iPt]->Signal(3,sgn,sgnerr);
    fitter[iPt]->Background(3,bkg,bkgerr);
    hRawYieldsSignificance->SetBinContent(iPt+1,signif);
    hRawYieldsSignificance->SetBinError(iPt+1,signiferr);
    hRawYieldsSoverB->SetBinContent(iPt+1,sgn/bkg);
    hRawYieldsSoverB->SetBinError(iPt+1,sgn/bkg*TMath::Sqrt(sgnerr/sgn*sgnerr/sgn+bkgerr/bkg*bkgerr/bkg));
    hRawYieldsSignal->SetBinContent(iPt+1,sgn);
    hRawYieldsSignal->SetBinError(iPt+1,sgnerr);
    hRawYieldsBkg->SetBinContent(iPt+1,bkg);
    hRawYieldsBkg->SetBinError(iPt+1,bkgerr);
    hRawYieldsChiSquare->SetBinContent(iPt+1,fitter[iPt]->GetReducedChiSquare());
    hRawYieldsChiSquare->SetBinError(iPt+1,0);
    if(nPtBins>1)
      cMass->cd(iPt+1);
    fitter[iPt]->DrawHere(gPad);
    cMass->Modified();
    cMass->Update();
    for(Int_t iDataFile=0; iDataFile<nDataFiles; iDataFile++) {
      ResetAxes(datasparse[iDataFile]);
    }
  }

  //efficiency
  /*
  TH1F* hEffPrompt=(TH1F*)GetEfficiency("hEffPrompt","hMassPtImpParPrompt","hMCAccPrompt",axesno,axesnames,cutset,nPtBins,PtLims,kBlue);
  TH1F* hEffFD=(TH1F*)GetEfficiency("hEffFD","hMassPtImpParBfeed","hMCAccBFeed",axesno,axesnames,cutset,nPtBins,PtLims,kRed);
  TCanvas *cEff=new TCanvas("cEff","efficiency",10,10,600,600);
  cEff->SetLogy();
  Double_t min;
  Double_t max;
  if(hEffPrompt->GetMinimum()<hEffFD->GetMinimum()) min = hEffPrompt->GetMinimum()*0.8;
  else min = hEffFD->GetMinimum()*0.8;
  if(hEffPrompt->GetMaximum()>hEffFD->GetMaximum()) max = hEffPrompt->GetMaximum()*1.5;
  else max = hEffFD->GetMaximum()*1.5;

  hEffFD->GetYaxis()->SetRangeUser(min,max);
  hEffFD->Draw("E1");
  hEffPrompt->Draw("E1same");
  */

  //number of events (corrected)
  TH1F* hEv = new TH1F("hEv","Number of events",1,0.5,1.5);
  Double_t nevents=0;
  for(Int_t iDataFile=0; iDataFile<nDataFiles; iDataFile++) {
    nevents += GetNevents(datafilename[iDataFile]);
  }
  hEv->SetBinContent(1,nevents);
  hEv->SetMarkerSize(1.);
  hEv->SetMarkerStyle(20);
  hEv->SetLineWidth(2);
  hEv->SetMarkerColor(kBlack);
  hEv->SetLineColor(kBlack);

  //output files
  TFile outfileraw(outfilerawname.Data(),"RECREATE");
  cMass->Write();
  hRawYields->Write();
  hRawYieldsSigma->Write();
  hRawYieldsMean->Write();
  hRawYieldsSignificance->Write();
  hRawYieldsSoverB->Write();
  hRawYieldsSignal->Write();
  hRawYieldsBkg->Write();
  hRawYieldsChiSquare->Write();
  hEv->Write();
  outfileraw.Close();
  cout << "\n" <<outfilerawname << " saved." <<endl;
  outfilerawname.ReplaceAll("root","pdf");
  cMass->SaveAs(outfilerawname.Data());
/*  TFile outfileeff(outfileeffname.Data(),"RECREATE");
  cEff->Write();
  hEffPrompt->Write();
  hEffFD->Write();
  outfileeff.Close();
  cout <<outfileeffname << " saved.\n" <<endl;
*/

  return 0;
}

//________________________________________________________________________________________________________________
Int_t ReadAxes(TString FileName, vector<string> &axesnames, vector<int> &axesno) {
  ifstream inAxes(FileName.Data());
  if(!inAxes){
    cerr<<"The file "<<FileName.Data() <<" does not exist. "<<endl;
    return 1;
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

  return 0;
}

//________________________________________________________________________________________________________________
Int_t ReadSet(TString FileName, vector<string> &varnames, vector<double> &cutset) {
  ifstream inSet(FileName.Data());
  if(!inSet){
    cerr<<"The file "<<FileName.Data() <<" does not exist. "<<endl;
    return 1;
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

  return 0;
}

//__________________________________________________________________________________________________________________
THnSparseF* GetSparse(TString filename, TString dirname, TString listname, TString sparsename) {
  cout << "opening file " << filename << "..." <<endl;

  TFile* infile = TFile::Open(filename.Data(),"READ");
  TDirectoryFile* dir = 0x0;
  TList* list = 0x0;
  THnSparseF* sparse = 0x0;

  if(infile) {
    dir = (TDirectoryFile*)infile->Get(dirname);
    if(dir) {
      list = (TList*)dir->Get(listname);
      if(list) {
        sparse = (THnSparseF*)list->FindObject(sparsename);
        if(!sparse)
          cout << "Error: Wrong THnSparse name! Exit." << endl;
      }
      else
        cout << "Error: Wrong TList name! Exit." << endl;
    }
    else
      cout << "Error: Wrong TDirectoryFile name! Exit." << endl;
  }

  if(infile && dir && list && sparse) {
    cout << "file opened and THnSparse " <<sparsename<< " got!" <<endl;
    cout << "file closed."<<endl;
  }

  return sparse;
}

//__________________________________________________________________________________________________________________
void ApplyCuts(THnSparseF* sparse, vector<int> axesno, vector<string> axesnames, vector<double> cutset, Int_t iPt) {
  cout << endl;
  for(UInt_t iAxis=0; iAxis<axesno.size(); iAxis++) {
    TAxis* ax = (TAxis*)sparse->GetAxis(axesno[iAxis]);
    Double_t binmin = ax->FindBin(cutset[2*iPt*axesno.size() + 2*iAxis]*1.00001);
    Double_t binmax = ax->FindBin(cutset[2*iPt*axesno.size() + 2*iAxis+1]*0.9999);
    cout << "Cutting " << axesnames[iAxis]  << " axis " << " from bin " <<binmin << " to bin " << binmax <<endl;
    ax->SetRange(binmin,binmax);
  }
  cout << endl;
}

//__________________________________________________________________________________________________________________
void ResetAxes(THnSparseF* sparse) {
  for(Int_t iAxis=0; iAxis<sparse->GetNdimensions(); iAxis++) {
    TAxis* ax = (TAxis*)sparse->GetAxis(iAxis);
    ax->SetRange(-1,-1);
  }
}

//__________________________________________________________________________________________________________________
void GetPtLims(Double_t ptlims[], vector<double> cutset, const Int_t nCutVars) {
  const Int_t nPtBins = cutset.size()/(2*nCutVars);
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ptlims[iPt] = cutset[iPt*2*nCutVars];
  }
  ptlims[nPtBins] = cutset[(nPtBins-1)*2*nCutVars+1];
}

//__________________________________________________________________________________________________________________
Double_t GetNevents(const TString filename) {

  TFile* infile = TFile::Open(filename.Data(),"READ");
  TDirectoryFile* dir = (TDirectoryFile*)infile->Get("PWG3_D2H_InvMassDplus");
  TString normobj=datalistname;
  normobj.ReplaceAll("coutputDplus","coutputDplusNorm");

  if(!infile || !dir)
    return 0;

  AliNormalizationCounter* nc=(AliNormalizationCounter*)dir->Get(normobj.Data());
  Double_t nEvents = nc->GetNEventsForNorm();

  infile->Close();

  return nEvents;
}

//__________________________________________________________________________________________________________________
TH1F* GetEfficiency(TString effhistoname, TString recosparsename, TString gensparsename, vector<int> axesno, vector<string> axesnames, vector<double> cutset, Int_t nPtBins, Double_t *PtLims, Int_t color) {

  //get thnsparse from MC file
  THnSparseF* recosparse = (THnSparseF*)GetSparse(MCfilename,MCdirname,MClistname,recosparsename);
  THnSparseF* gensparse = (THnSparseF*)GetSparse(MCfilename,MCdirname,MClistname,gensparsename);
  if(!recosparse || !gensparse)
    return 0x0;

  TH1F* hEff = new TH1F(effhistoname,"efficiency",nPtBins,PtLims);
  TH1F* hNumPt = new TH1F("hNumPt","",nPtBins,PtLims);
  TH1F* hDenPt = new TH1F("hDenPt","",nPtBins,PtLims);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ApplyCuts(recosparse,axesno,axesnames,cutset,iPt);
    TH1F* hRecoPt = (TH1F*)recosparse->Projection(1);
    TAxis* ptgenaxis = (TAxis*)gensparse->GetAxis(0);
    Int_t binmin = ptgenaxis->FindBin(PtLims[iPt]*1.0001);
    Int_t binmax = ptgenaxis->FindBin(PtLims[iPt+1]*0.9999);
    ptgenaxis->SetRange(binmin,binmax);
    TH1F* hGenPt = (TH1F*)gensparse->Projection(0);
    hNumPt->SetBinContent(iPt+1,hRecoPt->GetEntries());
    hNumPt->SetBinError(iPt+1,TMath::Sqrt(hRecoPt->GetEntries()));
    hDenPt->SetBinContent(iPt+1,hGenPt->GetEntries());
    hDenPt->SetBinError(iPt+1,TMath::Sqrt(hGenPt->GetEntries()));
    ResetAxes(recosparse);
    ResetAxes(gensparse);
  }
  hEff->Divide(hNumPt,hDenPt,1.,1.,"B");
  hEff->SetMarkerSize(1.);
  hEff->SetLineWidth(2);
  hEff->SetMarkerStyle(20);
  hEff->SetMarkerColor(color);
  hEff->SetLineColor(color);
  hEff->GetYaxis()->SetTitle("efficiency");
  hEff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");

  delete hNumPt;
  delete hDenPt;

  return hEff;

}

//__________________________________________________________________________________________________________________
void SetStyle() {
  cout << "Setting drawing style!" << endl;
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetOptStat(0);
}

//__________________________________________________________________________________________________________________
void DivideCanvas(TCanvas* c, const Int_t nPtBins) {

  if(nPtBins<2)
    c->Divide(1,1);
  if(nPtBins==2 || nPtBins==3)
    c->Divide(nPtBins,1);
  else if(nPtBins==4 || nPtBins==6 || nPtBins==8)
    c->Divide(nPtBins/2,2);
  else if(nPtBins==5 || nPtBins==7)
    c->Divide((nPtBins+1)/2,2);
  else if(nPtBins==9 || nPtBins==12 || nPtBins==15)
    c->Divide(nPtBins/3,3);
  else if(nPtBins==10 || nPtBins==11)
    c->Divide(4,3);
  else if(nPtBins==13 || nPtBins==14)
    c->Divide(5,3);
  else if(nPtBins>15 && nPtBins%2==0)
    c->Divide(nPtBins/2,2);
  else
    c->Divide((nPtBins+1)/2,2);
}
