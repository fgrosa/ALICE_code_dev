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

#include "AliHFMassFitter.h"

#include "AliHFCutVarFDsubAxis.h"
#include "AliHFCutVarFDsubCut.h"
#include "AliHFCutVarFDsubCutSet.h"
#include "AliHFCutVarFDsubEfficiency.h"
#include "AliHFCutVarFDsubMassFitter.h"
#include "AliHFCutVarFDsubMinimiser.h"
#include "AliHFCutVarFDsubAnalysisManager.h"
#include "AliHFCutVarFDsubAnalysisManagerDplus.h"

#endif

//macro for the standard analysis using the Cut Variation code for the D+
//author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino

//******************************************************//
//                                                      //
//    Main Functions: GetEfficienciesAndRawYields       //
//                                                      //
//******************************************************//

//________________________________________________________________________________________________________________
//Input files
TString datafile = "/home/fabrizio/ALICE/Files/Trains/LHC13/AnalysisResultspPbData.root";
TString datadir = "PWG3_D2H_InvMassDplus";
TString datalist = "coutputDplus_ImpParpPbData0100";

TString MCfile = "/home/fabrizio/ALICE/Files/Trains/LHC13/AnalysisResultspPbMC.root";
TString MCdir = "PWG3_D2H_InvMassDplus";
TString MClist = "coutputDplus_ImpParpPbMC0100";

TString axesfile = "axes";
TString cutsfile = "cutset";

//________________________________________________________________________________________________________________
//Pt bins
const Int_t nSets=1;
const Int_t nPtBins=10;
const Int_t nPtLims=nPtBins+1;
Double_t PtLims[nPtLims] = {1.,2.,3.,4.,5.,6.,7.,8.,12.,16.,24.};

//________________________________________________________________________________________________________________
//Parameters for invariant mass fits
Int_t rebin = 5;
Double_t massmin = 1.7;
Double_t massmax = 2.03;
Int_t funsig = AliHFMassFitter::kGaus;
Int_t funbkg = AliHFMassFitter::kExpo;

//________________________________________________________________________________________________________________
//functions prototypes
Int_t GetEfficienciesAndRawYields(Bool_t PID=kTRUE, Bool_t DrawMassFits=kTRUE);
void ReadAxes(TString FileName, vector<string> &axesanmes, vector<int> &axesno);
void ReadSet(TString FileName, vector<string> &varnames, vector<double> &cutset);

//________________________________________________________________________________________________________________
Int_t GetEfficienciesAndRawYields(Bool_t PID, Bool_t DrawMassFits) {

  //______________________________________________________________________________________________________________
  //cuts set
  vector<string> axesnames;
  vector<int> axesno;
  ReadAxes(axesfile,axesnames,axesno);
  const Int_t nCutVars=axesnames.size(); //pt included

  UInt_t dataAxesNo[nCutVars];
  UInt_t MCgenAxesNo[nCutVars];
  UInt_t MCcutAxesNo[nCutVars];
  TString AxesNames[nCutVars];
  Bool_t isCutSymm[nCutVars];

  for(Int_t iCutVar=0; iCutVar<nCutVars; iCutVar++) {
    AxesNames[iCutVar] = axesnames[iCutVar];
    dataAxesNo[iCutVar] = axesno[iCutVar];
    MCgenAxesNo[iCutVar] = axesno[iCutVar+nCutVars];
    MCcutAxesNo[iCutVar] = axesno[iCutVar+2*nCutVars];
    isCutSymm[iCutVar] = kFALSE;
  }

  vector<string> varnames;
  vector<double> cutset;
  ReadSet(cutsfile,varnames,cutset);
  Int_t nPtBins= cutset.size()/(2*nCutVars);
  const Int_t nptbins = nPtBins;

  Double_t*** cutlowset = new Double_t**[nSets]; //first: set, second: pt bin, third: cut variable
  Double_t*** cuthighset = new Double_t**[nSets];

  for(Int_t iSet=0; iSet<nSets; iSet++) {
    cutlowset[iSet] = new Double_t*[nPtBins];
    cuthighset[iSet] = new Double_t*[nPtBins];
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      cutlowset[iSet][iPt] = new Double_t[nCutVars];
      cuthighset[iSet][iPt] = new Double_t[nCutVars];
    }
  }

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    for(Int_t iCutVar=0; iCutVar<nCutVars; iCutVar++) {
      cutlowset[0][iPt][iCutVar] = cutset[(2*iCutVar)+(iPt*nCutVars*2)];
      cuthighset[0][iPt][iCutVar] = cutset[(2*iCutVar+1)+(iPt*nCutVars*2)];
     }
  }

  Double_t massD = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  Double_t **means = new Double_t*[nSets];
  Double_t **sigmas = new Double_t*[nSets];

  for(Int_t iSet=0; iSet<nSets; iSet++) {
    means[iSet] = new Double_t[nPtBins];
    sigmas[iSet] = new Double_t[nPtBins];

    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      means[iSet][iPt] = -massD;
      sigmas[iSet][iPt] = -0.008;
    }
  }

  //______________________________________________________________________________________________________________
  //get efficiencies and raw yields
  AliHFCutVarFDsubAnalysisManagerDplus *AnalysisManagerDplus = new AliHFCutVarFDsubAnalysisManagerDplus();
  AnalysisManagerDplus->SetPID(PID,3);
  Int_t loadTH = AnalysisManagerDplus->GetTHnSparses(MCfile,datafile,MCdir,datadir,MClist,datalist,kFALSE);

  if(loadTH>0)
    return 1;

  AnalysisManagerDplus->GetAxes(dataAxesNo,MCgenAxesNo,MCcutAxesNo,AxesNames,nCutVars,isCutSymm);
  AnalysisManagerDplus->GetCuts(cutlowset,cuthighset,means,sigmas,rebin,funsig,funbkg,massmin,massmax,nSets,nPtBins,nCutVars);

  AnalysisManagerDplus->GetXaxisInformation();
  AnalysisManagerDplus->GetEfficiencies(".","Efficiency.root");
  AnalysisManagerDplus->GetRawYields(DrawMassFits,".","Rawyields.root");

  //______________________________________________________________________________________________________________
  //delete dynamically allocated arrays
  for(Int_t iSet=0; iSet<nSets; iSet++) {
    delete means[iSet];
    delete sigmas[iSet];
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      delete cutlowset[iSet][iPt];
      delete cuthighset[iSet][iPt];
    }
  }
  delete[] means;
  delete[] sigmas;
  delete[] cutlowset;
  delete[] cuthighset;
  delete AnalysisManagerDplus;

  return 0;
}

void ReadAxes(TString FileName, vector<string> &axesnames, vector<int> &axesno) {
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
