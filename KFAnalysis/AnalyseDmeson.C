#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TSystem.h>
#include <TROOT.h>
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

#include "GetEfficiencyAndRawYields.C"
#include "CreateEfficFile.C"
#include "HFPtSpectrum.C"
#include "CombineFeedDownMCSubtractionMethodsUncertainties.C"

#endif

const Int_t nCutSets=1;
TString cutsfile[nCutSets] = {"cutsettopo"};
//TString cutsfile[nCutSets] = {"cutset1","cutset2","cutset3","cutset4","cutset5"};
//TString cutsfile[nCutSets] = {"cutset6","cutset7","cutset8","cutset9","cutset10"};
//TString cutsfile[nCutSets] = {"cutset11","cutset12","cutset13","cutset14","cutset15"};
//TString cutsfile[nCutSets] = {"cutset16","cutset17","cutset18","cutset19","cutset20"};
//TString cutsfile[nCutSets] = {"cutset21","cutset22","cutset23","cutset24","cutset25"};
//TString cutsfile[nCutSets] = {"cutset26","cutset27","cutset28","cutset29","cutset30"};
//TString cutsfile[nCutSets] = {"cutset31","cutset32","cutset33","cutset34","cutset35","cutset36"};

void AnalyseDmeson(Int_t meson=kDplus) {
  TString axesfile = "axes";
  TString outfile[nCutSets];
  TString fcfile[nCutSets];
  TString Nbfile[nCutSets];
  TString rawfile[nCutSets];
  TString rawfileMC[nCutSets];
  TString efffile[nCutSets];
  TString mesonname = "Dplus";
  if(meson==kDzero)
    mesonname = "Dzero";
  
  for(Int_t iSet=0; iSet<nCutSets; iSet++) {
    outfile[iSet] = Form("result/HFPtSpectrum_combinedFD_%s.root",cutsfile[iSet].Data());
    rawfile[iSet] = Form("result/rawyields_%s_%s.root",mesonname.Data(),cutsfile[iSet].Data());
    rawfileMC[iSet] = Form("result/rawyieldsMC_%s_%s.root",mesonname.Data(),cutsfile[iSet].Data());
    efffile[iSet] = Form("result/efficiency_%s_%s.root",mesonname.Data(),cutsfile[iSet].Data());
    fcfile[iSet] = Form("result/HFPtSpectrum_fc_%s.root",cutsfile[iSet].Data());
    Nbfile[iSet] = Form("result/HFPtSpectrum_Nb_%s.root",cutsfile[iSet].Data());  

    GetEfficiencyAndRawYields(kDplus,axesfile,cutsfile[iSet],rawfile[iSet],rawfileMC[iSet],efffile[iSet],kTRUE);
    CreateEfficFile(efffile[iSet]);
    HFPtSpectrum(1,fcfile[iSet],Form("result/AccEff_%s_%s.root",mesonname.Data(),cutsfile[iSet].Data()),rawfile[iSet]);
    HFPtSpectrum(2,Nbfile[iSet],Form("result/AccEff_%s_%s.root",mesonname.Data(),cutsfile[iSet].Data()),rawfile[iSet]);
    CombineFeedDownMCSubtractionMethodsUncertainties(fcfile[iSet],Nbfile[iSet],outfile[iSet]);
  }
}
