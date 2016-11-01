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
#include <TAxis.h>
#include <TGaxis.h>

#endif

enum {kDzero,kDplus};
const Int_t colors[] = {kRed,kBlue,kGreen+3,kBlack,kMagenta,kCyan,kOrange+7,kYellow+3,kGreen};

const Int_t nPtBins = 13;
const Int_t nPtLims = nPtBins+1;
Double_t PtLims[nPtLims] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16};

void PromptDQuality(Int_t meson=kDplus);
TH1F* GetHistoVsPt(THnSparseF* sparse, Int_t axis, Int_t scalefactor=1,TString outfilename="outfile.eps", Bool_t noFIT=kFALSE, Bool_t limits=kFALSE, Bool_t print=kFALSE, Bool_t singlegauss=kTRUE);
void SetPtRange(THnSparseF* sparse, Double_t ptmin, Double_t ptmax);
void ResetAxes(THnSparseF* sparse);
Double_t DoubleGauss(Double_t *x, Double_t *pars);

void PromptDQuality(Int_t meson) {
  TString mesonname;
  if(meson==kDzero) mesonname = "Dzero";
  else mesonname = "Dplus";
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHistLineWidth(2);
  gStyle->SetFuncWidth(2);
  
  TGaxis::SetMaxDigits(4);
  
  //INPUT FILES________________________________________________________________
  TString filename;
  TString listname;
  
  if(meson==kDzero) {
    //filename = "../fgrosa_Dzero_KF.root";
    filename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pPb_LHC13/MC/fgrosa_Dzero_KF_pPbMC_DLXY.root";
    listname = "coutputDzeroKF";
  }
  else {
    //filename = "../fgrosa_Dplus_KF.root";
    filename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pPb_LHC13/MC/fgrosa_Dplus_KF_pPbMC_DLXY.root";
    listname = "coutputDplusKF";
  }
  
  TFile infile(filename.Data(),"UPDATE");
  TList* list = (TList*)infile.Get(listname.Data());

  //D-MESON QA_________________________________________________________________
  THnSparseF* SparseResKF = (THnSparseF*)list->FindObject("fHistMres");
  THnSparseF* SparsePullsKF = (THnSparseF*)list->FindObject("fHistMpulls");
  THnSparseF* SparseResKFtopo = (THnSparseF*)list->FindObject("fHistMresTopo");
  THnSparseF* SparsePullsKFtopo = (THnSparseF*)list->FindObject("fHistMpullsTopo");
  THnSparseF* SparseResESD = (THnSparseF*)list->FindObject("fHistAliVertMres");
  THnSparseF* SparsePullsESD = (THnSparseF*)list->FindObject("fHistAliVertMpulls");
  TH2F* hChi2VsPtKF = (TH2F*)list->FindObject("fHistMchiS");
  TH2F* hChi2VsPtKFtopo = (TH2F*)list->FindObject("fHistMchiSTopo");
  TH2F* hProbVsPtKF = (TH2F*)list->FindObject("fHistMprob");
  TH2F* hProbVsPtKFtopo = (TH2F*)list->FindObject("fHistMprobTopo");

  infile.Close();
  TString outfilename = Form("Dmeson_%s_checks.root",mesonname.Data());
  TFile outfile(outfilename.Data(),"RECREATE");
  outfile.Close();

  TH1F* hChi2KF = (TH1F*)hChi2VsPtKF->ProjectionX();
  hChi2KF->SetDirectory(0);
  hChi2KF->SetLineColor(colors[0]);
  hChi2KF->SetLineWidth(2);
  hChi2KF->GetXaxis()->SetTitle("#chi^{2}/ndf");
  TH1F* hChi2KFtopo = (TH1F*)hChi2VsPtKFtopo->ProjectionX();
  hChi2KFtopo->SetDirectory(0);
  hChi2KFtopo->SetLineColor(colors[1]);
  hChi2KFtopo->SetLineWidth(2);

  TH1F* hProbKF = (TH1F*)hProbVsPtKF->ProjectionX();
  hProbKF->SetDirectory(0);
  hProbKF->SetLineColor(colors[0]);
  hProbKF->SetLineWidth(2);
  TH1F* hProbKFtopo = (TH1F*)hProbVsPtKFtopo->ProjectionX();
  hProbKFtopo->SetDirectory(0);
  hProbKFtopo->SetLineColor(colors[1]);
  hProbKFtopo->SetLineWidth(2);

  TH1F* hResXKF = GetHistoVsPt(SparseResKF,0,10000,outfilename,kFALSE,kTRUE);
  TH1F* hResYKF = GetHistoVsPt(SparseResKF,1,10000,outfilename,kFALSE,kTRUE);
  TH1F* hResZKF = GetHistoVsPt(SparseResKF,2,10000,outfilename,kFALSE,kTRUE,kTRUE);
  TH1F* hResPVXKF = GetHistoVsPt(SparseResKF,9,10000,outfilename,kFALSE,kTRUE);
  TH1F* hResPVYKF = GetHistoVsPt(SparseResKF,10,10000,outfilename,kFALSE,kTRUE);
  TH1F* hResPVZKF = GetHistoVsPt(SparseResKF,11,10000,outfilename,kFALSE,kTRUE);
  TH1F* hResMKF = GetHistoVsPt(SparseResKF,7,1,outfilename,kFALSE);
  TH1F* hResPtKF = GetHistoVsPt(SparseResKF,8,1,outfilename,kFALSE);
  TH1F* hResImpParKF = GetHistoVsPt(SparseResKF,14,10000,outfilename,kFALSE,kTRUE);
  hResXKF->SetName("hResXKF");
  hResYKF->SetName("hResYKF");
  hResZKF->SetName("hResZKF");
  hResPVXKF->SetName("hResPVXKF");
  hResPVYKF->SetName("hResPVYKF");
  hResPVZKF->SetName("hResPVZKF");
  hResMKF->SetName("hResMKF");
  hResPtKF->SetName("hResPtKF");
  hResImpParKF->SetName("hResImpParKF");   
  hResXKF->SetLineColor(colors[0]);
  hResYKF->SetLineColor(colors[0]);
  hResZKF->SetLineColor(colors[0]);
  hResPVXKF->SetLineColor(colors[0]);
  hResPVYKF->SetLineColor(colors[0]);
  hResPVZKF->SetLineColor(colors[0]);
  hResMKF->SetLineColor(colors[0]);
  hResPtKF->SetLineColor(colors[0]);
  hResImpParKF->SetLineColor(colors[0]);
  hResXKF->SetMarkerColor(colors[0]);
  hResYKF->SetMarkerColor(colors[0]);
  hResZKF->SetMarkerColor(colors[0]);
  hResPVXKF->SetMarkerColor(colors[0]);
  hResPVYKF->SetMarkerColor(colors[0]);
  hResPVZKF->SetMarkerColor(colors[0]);
  hResMKF->SetMarkerColor(colors[0]);
  hResPtKF->SetMarkerColor(colors[0]);
  hResImpParKF->SetMarkerColor(colors[0]);
  hResXKF->SetMarkerStyle(20);
  hResYKF->SetMarkerStyle(20);
  hResZKF->SetMarkerStyle(20);
  hResPVXKF->SetMarkerStyle(20);
  hResPVYKF->SetMarkerStyle(20);
  hResPVZKF->SetMarkerStyle(20);
  hResMKF->SetMarkerStyle(20);
  hResPtKF->SetMarkerStyle(20);
  hResImpParKF->SetMarkerStyle(21);
  hResXKF->SetMarkerSize(1.5);
  hResYKF->SetMarkerSize(1.5);
  hResZKF->SetMarkerSize(1.5);
  hResPVXKF->SetMarkerSize(1.5);
  hResPVYKF->SetMarkerSize(1.5);
  hResPVZKF->SetMarkerSize(1.5);
  hResMKF->SetMarkerSize(1.5);
  hResPtKF->SetMarkerSize(1.5);
  hResImpParKF->SetMarkerSize(1.5);
  
  TH1F* hResXKFtopo = GetHistoVsPt(SparseResKFtopo,0,10000,outfilename,kFALSE,kTRUE);
  TH1F* hResYKFtopo = GetHistoVsPt(SparseResKFtopo,1,10000,outfilename,kFALSE,kTRUE);
  TH1F* hResZKFtopo = GetHistoVsPt(SparseResKFtopo,2,10000,outfilename,kFALSE,kTRUE);
  TH1F* hResPVXKFtopo = GetHistoVsPt(SparseResKFtopo,9,10000,outfilename,kFALSE);
  TH1F* hResPVYKFtopo = GetHistoVsPt(SparseResKFtopo,10,10000,outfilename,kFALSE);
  TH1F* hResPVZKFtopo = GetHistoVsPt(SparseResKFtopo,11,10000,outfilename,kFALSE);
  TH1F* hResMKFtopo = GetHistoVsPt(SparseResKFtopo,7,1,outfilename,kFALSE);
  TH1F* hResPtKFtopo = GetHistoVsPt(SparseResKFtopo,8,1,outfilename,kFALSE);
  TH1F* hResDecLKFtopo = GetHistoVsPt(SparseResKFtopo,13,10000,outfilename,kFALSE);
  TH1F* hResImpParKFtopo = GetHistoVsPt(SparseResKFtopo,14,10000,outfilename,kFALSE,kTRUE);
  hResXKFtopo->SetName("hResXKFtopo");
  hResYKFtopo->SetName("hResYKFtopo");
  hResZKFtopo->SetName("hResZKFtopo");
  hResPVXKFtopo->SetName("hResPVXKFtopo");
  hResPVYKFtopo->SetName("hResPVYKFtopo");
  hResPVZKFtopo->SetName("hResPVZKFtopo");
  hResMKFtopo->SetName("hResMKFtopo");
  hResPtKFtopo->SetName("hResPtKFtopo");
  hResDecLKFtopo->SetName("hResDecLKFtopo");   
  hResImpParKFtopo->SetName("hResImpParKFtopo");   
  hResXKFtopo->SetLineColor(colors[1]);
  hResYKFtopo->SetLineColor(colors[1]);
  hResZKFtopo->SetLineColor(colors[1]);
  hResPVXKFtopo->SetLineColor(colors[1]);
  hResPVYKFtopo->SetLineColor(colors[1]);
  hResPVZKFtopo->SetLineColor(colors[1]);
  hResMKFtopo->SetLineColor(colors[1]);
  hResPtKFtopo->SetLineColor(colors[1]);
  hResDecLKFtopo->SetLineColor(colors[1]);
  hResImpParKFtopo->SetLineColor(colors[1]);
  hResXKFtopo->SetMarkerColor(colors[1]);
  hResYKFtopo->SetMarkerColor(colors[1]);
  hResZKFtopo->SetMarkerColor(colors[1]);
  hResPVXKFtopo->SetMarkerColor(colors[1]);
  hResPVYKFtopo->SetMarkerColor(colors[1]);
  hResPVZKFtopo->SetMarkerColor(colors[1]);
  hResMKFtopo->SetMarkerColor(colors[1]);
  hResPtKFtopo->SetMarkerColor(colors[1]);
  hResDecLKFtopo->SetMarkerColor(colors[1]);
  hResImpParKFtopo->SetMarkerColor(colors[1]);
  hResXKFtopo->SetMarkerStyle(20);
  hResYKFtopo->SetMarkerStyle(20);
  hResZKFtopo->SetMarkerStyle(20);
  hResPVXKFtopo->SetMarkerStyle(20);
  hResPVYKFtopo->SetMarkerStyle(20);
  hResPVZKFtopo->SetMarkerStyle(20);
  hResMKFtopo->SetMarkerStyle(20);
  hResPtKFtopo->SetMarkerStyle(20);
  hResDecLKFtopo->SetMarkerStyle(20);
  hResImpParKFtopo->SetMarkerStyle(20);
  hResXKFtopo->SetMarkerSize(1.5);
  hResYKFtopo->SetMarkerSize(1.5);
  hResZKFtopo->SetMarkerSize(1.5);
  hResPVXKFtopo->SetMarkerSize(1.5);
  hResPVYKFtopo->SetMarkerSize(1.5);
  hResPVZKFtopo->SetMarkerSize(1.5);
  hResMKFtopo->SetMarkerSize(1.5);
  hResPtKFtopo->SetMarkerSize(1.5);
  hResDecLKFtopo->SetMarkerSize(1.5);
  hResImpParKFtopo->SetMarkerSize(1.5);
  
  TH1F* hResXESD = GetHistoVsPt(SparseResESD,0,10000,outfilename,kFALSE,kTRUE);
  TH1F* hResYESD = GetHistoVsPt(SparseResESD,1,10000,outfilename,kFALSE,kTRUE);
  TH1F* hResZESD = GetHistoVsPt(SparseResESD,2,10000,outfilename,kFALSE,kTRUE);
  TH1F* hResMESD = GetHistoVsPt(SparseResESD,7,1,outfilename,kFALSE);
  TH1F* hResPtESD = GetHistoVsPt(SparseResESD,8,1,outfilename,kFALSE);
  TH1F* hResDecLESD = GetHistoVsPt(SparseResESD,10,10000,outfilename,kFALSE);
  TH1F* hResImpParESD = GetHistoVsPt(SparseResESD,11,10000,outfilename,kFALSE,kTRUE);
  hResXESD->SetName("hResXESD");
  hResYESD->SetName("hResYESD");
  hResZESD->SetName("hResZESD");
  hResMESD->SetName("hResMESD");
  hResPtESD->SetName("hResPtESD");
  hResDecLESD->SetName("hResDecLESD");
  hResImpParESD->SetName("hResImpParESD");
  hResXESD->SetLineColor(colors[2]);
  hResYESD->SetLineColor(colors[2]);
  hResZESD->SetLineColor(colors[2]);
  hResMESD->SetLineColor(colors[2]);
  hResPtESD->SetLineColor(colors[2]);
  hResDecLESD->SetLineColor(colors[2]);
  hResImpParESD->SetLineColor(colors[2]);
  hResXESD->SetMarkerColor(colors[2]);
  hResYESD->SetMarkerColor(colors[2]);
  hResZESD->SetMarkerColor(colors[2]);
  hResMESD->SetMarkerColor(colors[2]);
  hResPtESD->SetMarkerColor(colors[2]);
  hResDecLESD->SetMarkerColor(colors[2]);
  hResImpParESD->SetMarkerColor(colors[2]);
  hResXESD->SetMarkerStyle(20);
  hResYESD->SetMarkerStyle(20);
  hResZESD->SetMarkerStyle(20);
  hResMESD->SetMarkerStyle(20);
  hResPtESD->SetMarkerStyle(20);
  hResDecLESD->SetMarkerStyle(22);
  hResImpParESD->SetMarkerStyle(22);
  hResXESD->SetMarkerSize(1.5);
  hResYESD->SetMarkerSize(1.5);
  hResZESD->SetMarkerSize(1.5);
  hResMESD->SetMarkerSize(1.5);
  hResPtESD->SetMarkerSize(1.5);
  hResDecLESD->SetMarkerSize(1.5);
  hResImpParESD->SetMarkerSize(1.5);
  
  TH1F* hPullsXKF = GetHistoVsPt(SparsePullsKF,0,1,outfilename,kFALSE);
  TH1F* hPullsYKF = GetHistoVsPt(SparsePullsKF,1,1,outfilename,kFALSE);
  TH1F* hPullsZKF = GetHistoVsPt(SparsePullsKF,2,1,outfilename,kFALSE,kTRUE,kTRUE);
  TH1F* hPullsPVXKF = GetHistoVsPt(SparsePullsKF,9,1,outfilename,kFALSE);
  TH1F* hPullsPVYKF = GetHistoVsPt(SparsePullsKF,10,1,outfilename,kFALSE);
  TH1F* hPullsPVZKF = GetHistoVsPt(SparsePullsKF,11,1,outfilename,kFALSE);
  TH1F* hPullsMKF = GetHistoVsPt(SparsePullsKF,7,1,outfilename,kFALSE);
  TH1F* hPullsPtKF = GetHistoVsPt(SparsePullsKF,8,1,outfilename,kFALSE);
  TH1F* hPullsImpParKF = GetHistoVsPt(SparsePullsKF,14,1,outfilename,kFALSE);
  hPullsXKF->SetName("hPullsXKF");
  hPullsYKF->SetName("hPullsYKF");
  hPullsZKF->SetName("hPullsZKF");
  hPullsPVXKF->SetName("hPullsPVXKF");
  hPullsPVYKF->SetName("hPullsPVYKF");
  hPullsPVZKF->SetName("hPullsPVZKF");
  hPullsMKF->SetName("hPullsMKF");
  hPullsPtKF->SetName("hPullsPtKF");
  hPullsImpParKF->SetName("hPullsImpParKF");   
  hPullsXKF->SetLineColor(colors[0]);
  hPullsYKF->SetLineColor(colors[0]);
  hPullsZKF->SetLineColor(colors[0]);
  hPullsPVXKF->SetLineColor(colors[0]);
  hPullsPVYKF->SetLineColor(colors[0]);
  hPullsPVZKF->SetLineColor(colors[0]);
  hPullsMKF->SetLineColor(colors[0]);
  hPullsPtKF->SetLineColor(colors[0]);
  hPullsImpParKF->SetLineColor(colors[0]);
  hPullsXKF->SetMarkerColor(colors[0]);
  hPullsYKF->SetMarkerColor(colors[0]);
  hPullsZKF->SetMarkerColor(colors[0]);
  hPullsPVXKF->SetMarkerColor(colors[0]);
  hPullsPVYKF->SetMarkerColor(colors[0]);
  hPullsPVZKF->SetMarkerColor(colors[0]);
  hPullsMKF->SetMarkerColor(colors[0]);
  hPullsPtKF->SetMarkerColor(colors[0]);
  hPullsImpParKF->SetMarkerColor(colors[0]);
  hPullsXKF->SetMarkerStyle(20);
  hPullsYKF->SetMarkerStyle(20);
  hPullsZKF->SetMarkerStyle(20);
  hPullsPVXKF->SetMarkerStyle(20);
  hPullsPVYKF->SetMarkerStyle(20);
  hPullsPVZKF->SetMarkerStyle(20);
  hPullsMKF->SetMarkerStyle(20);
  hPullsPtKF->SetMarkerStyle(20);
  hPullsImpParKF->SetMarkerStyle(21);
  hPullsXKF->SetMarkerSize(1.5);
  hPullsYKF->SetMarkerSize(1.5);
  hPullsZKF->SetMarkerSize(1.5);
  hPullsPVXKF->SetMarkerSize(1.5);
  hPullsPVYKF->SetMarkerSize(1.5);
  hPullsPVZKF->SetMarkerSize(1.5);
  hPullsMKF->SetMarkerSize(1.5);
  hPullsPtKF->SetMarkerSize(1.5);
  hPullsImpParKF->SetMarkerSize(1.5);

  TH1F* hPullsXKFtopo = GetHistoVsPt(SparsePullsKFtopo,0,1,outfilename,kFALSE);
  TH1F* hPullsYKFtopo = GetHistoVsPt(SparsePullsKFtopo,1,1,outfilename,kFALSE);
  TH1F* hPullsZKFtopo = GetHistoVsPt(SparsePullsKFtopo,2,1,outfilename,kFALSE);
  TH1F* hPullsPVXKFtopo = GetHistoVsPt(SparsePullsKFtopo,9,1,outfilename,kFALSE);
  TH1F* hPullsPVYKFtopo = GetHistoVsPt(SparsePullsKFtopo,10,1,outfilename,kFALSE);
  TH1F* hPullsPVZKFtopo = GetHistoVsPt(SparsePullsKFtopo,11,1,outfilename,kFALSE);
  TH1F* hPullsMKFtopo = GetHistoVsPt(SparsePullsKFtopo,7,1,outfilename,kFALSE);
  TH1F* hPullsPtKFtopo = GetHistoVsPt(SparsePullsKFtopo,8,1,outfilename,kFALSE);
  TH1F* hPullsDecLKFtopo = GetHistoVsPt(SparsePullsKFtopo,13,1,outfilename,kFALSE);
  TH1F* hPullsImpParKFtopo = GetHistoVsPt(SparsePullsKFtopo,14,1,outfilename,kFALSE);
  hPullsXKFtopo->SetName("hPullsXKFtopo");
  hPullsYKFtopo->SetName("hPullsYKFtopo");
  hPullsZKFtopo->SetName("hPullsZKFtopo");
  hPullsPVXKFtopo->SetName("hPullsPVXKFtopo");
  hPullsPVYKFtopo->SetName("hPullsPVYKFtopo");
  hPullsPVZKFtopo->SetName("hPullsPVZKFtopo");
  hPullsMKFtopo->SetName("hPullsMKFtopo");
  hPullsPtKFtopo->SetName("hPullsPtKFtopo");
  hPullsDecLKFtopo->SetName("hPullsDecLKFtopo");   
  hPullsImpParKFtopo->SetName("hPullsImpParKFtopo");   
  hPullsXKFtopo->SetLineColor(colors[1]);
  hPullsYKFtopo->SetLineColor(colors[1]);
  hPullsZKFtopo->SetLineColor(colors[1]);
  hPullsPVXKFtopo->SetLineColor(colors[1]);
  hPullsPVYKFtopo->SetLineColor(colors[1]);
  hPullsPVZKFtopo->SetLineColor(colors[1]);
  hPullsMKFtopo->SetLineColor(colors[1]);
  hPullsPtKFtopo->SetLineColor(colors[1]);
  hPullsDecLKFtopo->SetLineColor(colors[1]);
  hPullsImpParKFtopo->SetLineColor(colors[1]);
  hPullsXKFtopo->SetMarkerColor(colors[1]);
  hPullsYKFtopo->SetMarkerColor(colors[1]);
  hPullsZKFtopo->SetMarkerColor(colors[1]);
  hPullsPVXKFtopo->SetMarkerColor(colors[1]);
  hPullsPVYKFtopo->SetMarkerColor(colors[1]);
  hPullsPVZKFtopo->SetMarkerColor(colors[1]);
  hPullsMKFtopo->SetMarkerColor(colors[1]);
  hPullsPtKFtopo->SetMarkerColor(colors[1]);
  hPullsDecLKFtopo->SetMarkerColor(colors[1]);
  hPullsImpParKFtopo->SetMarkerColor(colors[1]);
  hPullsXKFtopo->SetMarkerStyle(20);
  hPullsYKFtopo->SetMarkerStyle(20);
  hPullsZKFtopo->SetMarkerStyle(20);
  hPullsPVXKFtopo->SetMarkerStyle(20);
  hPullsPVYKFtopo->SetMarkerStyle(20);
  hPullsPVZKFtopo->SetMarkerStyle(20);
  hPullsMKFtopo->SetMarkerStyle(20);
  hPullsPtKFtopo->SetMarkerStyle(20);
  hPullsDecLKFtopo->SetMarkerStyle(20);
  hPullsImpParKFtopo->SetMarkerStyle(20);
  hPullsXKFtopo->SetMarkerSize(1.5);
  hPullsYKFtopo->SetMarkerSize(1.5);
  hPullsZKFtopo->SetMarkerSize(1.5);
  hPullsPVXKFtopo->SetMarkerSize(1.5);
  hPullsPVYKFtopo->SetMarkerSize(1.5);
  hPullsPVZKFtopo->SetMarkerSize(1.5);
  hPullsMKFtopo->SetMarkerSize(1.5);
  hPullsPtKFtopo->SetMarkerSize(1.5);
  hPullsDecLKFtopo->SetMarkerSize(1.5);
  hPullsImpParKFtopo->SetMarkerSize(1.5);
  
  TH1F* hPullsXESD = GetHistoVsPt(SparsePullsESD,0,1,outfilename,kFALSE);
  TH1F* hPullsYESD = GetHistoVsPt(SparsePullsESD,1,1,outfilename,kFALSE);
  TH1F* hPullsZESD = GetHistoVsPt(SparsePullsESD,2,1,outfilename,kFALSE);
  TH1F* hPullsDecLESD = GetHistoVsPt(SparsePullsESD,7,1,outfilename,kFALSE);
  hPullsXESD->SetName("hPullsXESD");
  hPullsYESD->SetName("hPullsYESD");
  hPullsZESD->SetName("hPullsZESD");
  hPullsDecLESD->SetName("hPullsDecLESD");
  hPullsXESD->SetLineColor(colors[2]);
  hPullsYESD->SetLineColor(colors[2]);
  hPullsZESD->SetLineColor(colors[2]);
  hPullsDecLESD->SetLineColor(colors[2]);
  hPullsXESD->SetMarkerColor(colors[2]);
  hPullsYESD->SetMarkerColor(colors[2]);
  hPullsZESD->SetMarkerColor(colors[2]);
  hPullsDecLESD->SetMarkerColor(colors[2]);
  hPullsXESD->SetMarkerStyle(20);
  hPullsYESD->SetMarkerStyle(20);
  hPullsZESD->SetMarkerStyle(20);
  hPullsDecLESD->SetMarkerStyle(22);
  hPullsXESD->SetMarkerSize(1.5);
  hPullsYESD->SetMarkerSize(1.5);
  hPullsZESD->SetMarkerSize(1.5);
  hPullsDecLESD->SetMarkerSize(1.5);

  TLegend *lchi = new TLegend(0.4,0.65,0.8,0.85);
  lchi->SetTextSize(0.045);
  lchi->AddEntry(hChi2KF,"KF","l");
  lchi->AddEntry(hChi2KFtopo,"KF topo constraint","l");

  TLegend *l2 = new TLegend(0.4,0.65,0.85,0.85);
  l2->SetTextSize(0.045);
  l2->AddEntry(hResImpParKF,"KF","lpe");
  l2->AddEntry(hResImpParKFtopo,"KF topo constraint","lpe");

  TLegend *l3 = new TLegend(0.4,0.65,0.85,0.85);
  l3->SetTextSize(0.045);
  l3->AddEntry(hResImpParKF,"KF","lpe");
  l3->AddEntry(hResImpParKFtopo,"KF topo constraint","lpe");
  l3->AddEntry(hResImpParESD,"AliVertexer","lpe");

  TLegend *l4 = new TLegend(0.2,0.65,0.5,0.85);
  l4->SetTextSize(0.045);
  l4->AddEntry(hResImpParKF,"KF","lpe");
  l4->AddEntry(hResImpParKFtopo,"KF topo constraint","lpe");
  l4->AddEntry(hResImpParESD,"4-vector","lpe");

  TLegend *ld = new TLegend(0.2,0.65,0.5,0.85);
  ld->SetTextSize(0.045);
  ld->AddEntry(hResDecLKFtopo,"KF topo constraint","lpe");
  ld->AddEntry(hResDecLESD,"AliVertexer","lpe");

  TLatex* latKF = new TLatex();
  latKF->SetTextFont(132);
  latKF->SetTextSize(0.05);
  latKF->SetTextColor(colors[0]);
  TLatex* latKFtopo = new TLatex();
  latKFtopo->SetTextFont(132);
  latKFtopo->SetTextSize(0.05);
  latKFtopo->SetTextColor(colors[1]);
  
  TCanvas* cChi2 = new TCanvas("cChi2","",800,800);
  cChi2->SetLogy();
  hChi2KF->SetTitle("");
  hChi2KF->GetYaxis()->SetTitle("Entries");
  hChi2KF->GetXaxis()->SetTitleSize(0.05);
  hChi2KF->GetXaxis()->SetLabelSize(0.05);
  hChi2KF->Draw();
  hChi2KFtopo->Draw("same");
  latKF->DrawLatex(8,hChi2KF->GetMaximum()*0.5,"KF");
  latKF->DrawLatex(8,hChi2KF->GetMaximum()*0.2,Form("<#chi^{2}/ndf> = %0.3f #pm %0.3f",hChi2KF->GetMean(),hChi2KF->GetMeanError()));
  latKFtopo->DrawLatex(8,hChi2KF->GetMaximum()*0.05,"KF topo constraint");
  latKFtopo->DrawLatex(8,hChi2KF->GetMaximum()*0.02,Form("<#chi^{2}/ndf> = %0.3f #pm %0.3f",hChi2KFtopo->GetMean(),hChi2KFtopo->GetMeanError()));
  
  TCanvas* cProb = new TCanvas("cProb","",800,800);
  cProb->SetLogy();
  hProbKF->SetTitle("");
  hProbKF->GetYaxis()->SetTitle("Entries");
  hProbKF->GetXaxis()->SetTitle("Probability");
  hProbKF->Draw();
  hProbKF->GetXaxis()->SetTitleSize(0.05);
  hProbKF->GetXaxis()->SetLabelSize(0.05);
  hProbKF->GetXaxis()->SetNdivisions(505);
  hProbKFtopo->Draw("same");
  latKF->DrawLatex(0.2,hProbKF->GetMaximum()*1.25,"KF");
  latKF->DrawLatex(0.2,hProbKF->GetMaximum()*0.95,Form("<prob> = %0.3f #pm %0.3f",hProbKF->GetMean(),hProbKF->GetMeanError()));
  latKFtopo->DrawLatex(0.2,hProbKF->GetMaximum()*0.65,"KF topo constraint");
  latKFtopo->DrawLatex(0.2,hProbKF->GetMaximum()*0.5,Form("<prob> = %0.3f #pm %0.3f",hProbKFtopo->GetMean(),hProbKFtopo->GetMeanError()));

  TCanvas *cResSVX = new TCanvas("cResSVX","",800,800);
  hResXKF->GetYaxis()->SetRangeUser(0.0,180);
  hResXKF->Draw();
  hResXKFtopo->Draw("same");
  hResXESD->Draw("same");
  l3->Draw("same");
  TCanvas *cResSVY = new TCanvas("cResSVY","",800,800);
  hResYKF->GetYaxis()->SetRangeUser(0.0,180);
  hResYKF->Draw();
  hResYKFtopo->Draw("same");
  hResYESD->Draw("same");
  l3->Draw("same");
  TCanvas *cResSVZ = new TCanvas("cResSVZ","",800,800);
  hResZKF->GetYaxis()->SetRangeUser(0.0,180);
  hResZKF->Draw();
  hResZKFtopo->Draw("same");
  hResZESD->Draw("same");
  l3->Draw("same");
  TCanvas *cResSVXZ = new TCanvas("cResSVXZ","",800,800);
  TH1F* hResSVXKFCopy=(TH1F*)hResXKF->Clone();
  TH1F* hResSVZKFCopy=(TH1F*)hResZKF->Clone();
  TH1F* hResSVXESDCopy=(TH1F*)hResXESD->Clone();
  TH1F* hResSVZESDCopy=(TH1F*)hResZESD->Clone();
  TH1F* hResSVXKFtopoCopy=(TH1F*)hResXKFtopo->Clone();
  TH1F* hResSVZKFtopoCopy=(TH1F*)hResZKFtopo->Clone();
  TLegend* lXZSV = new TLegend(0.35,0.6,0.89,0.89);
  lXZSV->SetTextSize(0.045);
  lXZSV->AddEntry(hResSVXKFCopy,"KF x","lpe");
  lXZSV->AddEntry(hResSVZKFCopy,"KF z","lpe");
  lXZSV->AddEntry(hResSVXKFtopoCopy,"KF topo constraint x","lpe");
  lXZSV->AddEntry(hResSVZKFtopoCopy,"KF topo constraint z","lpe");
  lXZSV->AddEntry(hResSVXESDCopy,"AliVertexer x","lpe");
  lXZSV->AddEntry(hResSVZESDCopy,"AliVertexer z","lpe");
  hResSVZKFCopy->SetMarkerStyle(26);
  hResSVZKFtopoCopy->SetMarkerStyle(26);
  hResSVZESDCopy->SetMarkerStyle(26);
  hResSVXKFCopy->GetYaxis()->SetTitle("decay vtx #sigma(pos_{reco}-pos_{MC}) (#mum)");
  hResSVXKFCopy->GetYaxis()->SetRangeUser(50,180);
  hResSVXKFCopy->Draw();
  hResSVXKFtopoCopy->Draw("same");
  hResSVZKFCopy->Draw("same");
  hResSVZKFtopoCopy->Draw("same");
  hResSVXESDCopy->Draw("same");
  hResSVZESDCopy->Draw("same");
  lXZSV->Draw("same");
  TCanvas *cResPVX = new TCanvas("cResPVX","",800,800);
  hResPVXKF->GetYaxis()->SetRangeUser(0.,150);
  hResPVXKF->Draw();
  hResPVXKFtopo->Draw("same");
  l2->Draw("same");
  TCanvas *cResPVY = new TCanvas("cResPVY","",800,800);
  hResPVYKF->GetYaxis()->SetRangeUser(0.,150);
  hResPVYKF->Draw();
  hResPVYKFtopo->Draw("same");
  l2->Draw("same");
  TCanvas *cResPVZ = new TCanvas("cResPVZ","",800,800);
  hResPVZKF->GetYaxis()->SetRangeUser(0.,150);
  hResPVZKF->Draw();
  hResPVZKFtopo->Draw("same");
  l2->Draw("same");
  TCanvas *cResPVXZ = new TCanvas("cResPVXZ","",800,800);
  TH1F* hResPVXKFCopy=(TH1F*)hResPVXKF->Clone();
  TH1F* hResPVZKFCopy=(TH1F*)hResPVZKF->Clone();
  TH1F* hResPVXKFtopoCopy=(TH1F*)hResPVXKFtopo->Clone();
  TH1F* hResPVZKFtopoCopy=(TH1F*)hResPVZKFtopo->Clone();
  TLegend* lXZPV = new TLegend(0.35,0.6,0.89,0.89);
  lXZPV->SetTextSize(0.045);
  lXZPV->AddEntry(hResPVXKFCopy,"KF x","lpe");
  lXZPV->AddEntry(hResPVZKFCopy,"KF z","lpe");
  lXZPV->AddEntry(hResPVXKFtopoCopy,"KF topo constraint x","lpe");
  lXZPV->AddEntry(hResPVZKFtopoCopy,"KF topo constraint z","lpe");
  hResPVZKFCopy->SetMarkerStyle(26);
  hResPVZKFtopoCopy->SetMarkerStyle(26);
  hResPVXKFCopy->GetYaxis()->SetTitle("prod vtx #sigma(pos_{reco}-pos_{MC}) (#mum)");
  hResPVXKFCopy->Draw();
  hResPVXKFtopoCopy->Draw("same");
  hResPVZKFCopy->Draw("same");
  hResPVZKFtopoCopy->Draw("same");
  lXZPV->Draw("same");
  TCanvas *cResSVM = new TCanvas("cResSVM","",800,800);
  hResMKF->GetYaxis()->SetRangeUser(60,180);
  hResMKF->Draw();
  hResMKFtopo->Draw("same");
  hResMESD->Draw("same");
  l4->Draw("same");
  TCanvas *cResSVPt = new TCanvas("cResSVPt","",800,800);
  hResPtKF->GetYaxis()->SetRangeUser(0.001,0.2);
  hResPtKF->Draw();
  hResPtKFtopo->Draw("same");
  hResPtESD->Draw("same");
  l4->Draw("same");
  TCanvas *cResDecL = new TCanvas("cResDecL","",800,800);
  hResDecLKFtopo->GetYaxis()->SetRangeUser(100,240);
  hResDecLKFtopo->Draw();
  hResDecLESD->Draw("same");
  ld->Draw("same");
  TCanvas *cResImpPar = new TCanvas("cResImpPar","",800,800);
  hResImpParKFtopo->GetYaxis()->SetRangeUser(0.,140);
  hResImpParKFtopo->Draw();
  hResImpParKF->Draw("same");
  hResImpParESD->Draw("same");
  l3->Draw("same");
  
  TCanvas *cPullsSVX = new TCanvas("cPullsSVX","",800,800);
  hPullsXKF->GetYaxis()->SetRangeUser(0.,2.);
  hPullsXKF->Draw();
  hPullsXKFtopo->Draw("same");
  hPullsXESD->Draw("same");
  l3->Draw("same");
  TCanvas *cPullsSVY = new TCanvas("cPullsSVY","",800,800);
  hPullsYKF->GetYaxis()->SetRangeUser(0.,2.);
  hPullsYKF->Draw();
  hPullsYKFtopo->Draw("same");
  hPullsYESD->Draw("same");
  l3->Draw("same");
  TCanvas *cPullsSVZ = new TCanvas("cPullsSVZ","",800,800);
  hPullsZKF->GetYaxis()->SetRangeUser(0.,2.5);
  hPullsZKF->Draw();
  hPullsZKFtopo->Draw("same");
  hPullsZESD->Draw("same");
  l3->Draw("same");
  TCanvas *cPullsSVXZ = new TCanvas("cPullsSVXZ","",800,800);
  TH1F* hPullsSVXKFCopy=(TH1F*)hPullsXKF->Clone();
  TH1F* hPullsSVZKFCopy=(TH1F*)hPullsZKF->Clone();
  TH1F* hPullsSVXESDCopy=(TH1F*)hPullsXESD->Clone();
  TH1F* hPullsSVZESDCopy=(TH1F*)hPullsZESD->Clone();
  TH1F* hPullsSVXKFtopoCopy=(TH1F*)hPullsXKFtopo->Clone();
  TH1F* hPullsSVZKFtopoCopy=(TH1F*)hPullsZKFtopo->Clone();
  hPullsSVZKFCopy->SetMarkerStyle(26);
  hPullsSVZESDCopy->SetMarkerStyle(26);
  hPullsSVZKFtopoCopy->SetMarkerStyle(26);
  hPullsSVXKFCopy->GetYaxis()->SetTitle("decay vtx #sigma((pos_{reco}-pos_{MC})/#sigma_{pos_{reco}})");
  hPullsSVXKFCopy->GetYaxis()->SetRangeUser(0,2.5);
  hPullsSVXKFCopy->Draw();
  hPullsSVXKFtopoCopy->Draw("same");
  hPullsSVZKFCopy->Draw("same");
  hPullsSVZKFtopoCopy->Draw("same");
  hPullsSVXESDCopy->Draw("same");
  hPullsSVZESDCopy->Draw("same");
  TLegend* lXZSV2 = (TLegend*)lXZSV->Clone();
  lXZSV2->SetX1(0.2);
  lXZSV2->SetX2(0.74);
  lXZSV2->Draw("same");
  TCanvas *cPullsPVX = new TCanvas("cPullsPVX","",800,800);
  hPullsPVXKF->GetYaxis()->SetRangeUser(0.5,1.8);
  hPullsPVXKF->Draw();
  hPullsPVXKFtopo->Draw("same");
  l2->Draw("same");
  TCanvas *cPullsPVY = new TCanvas("cPullsPVY","",800,800);
  hPullsPVYKF->GetYaxis()->SetRangeUser(0.5,1.5);
  hPullsPVYKF->Draw();
  hPullsPVYKFtopo->Draw("same");
  l2->Draw("same");
  TCanvas *cPullsPVZ = new TCanvas("cPullsPVZ","",800,800);
  hPullsPVZKF->GetYaxis()->SetRangeUser(0.9,1.3);
  hPullsPVZKF->Draw();
  hPullsPVZKFtopo->Draw("same");
  l2->Draw("same");
  TCanvas *cPullsPVXZ = new TCanvas("cPullsPVXZ","",800,800);
  TH1F* hPullsPVXKFCopy=(TH1F*)hPullsPVXKF->Clone();
  TH1F* hPullsPVZKFCopy=(TH1F*)hPullsPVZKF->Clone();
  TH1F* hPullsPVXKFtopoCopy=(TH1F*)hPullsPVXKFtopo->Clone();
  TH1F* hPullsPVZKFtopoCopy=(TH1F*)hPullsPVZKFtopo->Clone();
  hPullsPVZKFCopy->SetMarkerStyle(26);
  hPullsPVZKFtopoCopy->SetMarkerStyle(26);
  hPullsPVXKFCopy->GetYaxis()->SetTitle("prod vtx #sigma((pos_{reco}-pos_{MC})/#sigma_{pos_{reco}})");
  hPullsPVXKFCopy->Draw();
  hPullsPVXKFtopoCopy->Draw("same");
  hPullsPVZKFCopy->Draw("same");
  hPullsPVZKFtopoCopy->Draw("same");
  lXZPV->Draw("same");
  TCanvas *cPullsSVM = new TCanvas("cPullsSVM","",800,800);
  hPullsMKF->GetYaxis()->SetRangeUser(0.5,1.5);
  hPullsMKF->Draw();
  hPullsMKFtopo->Draw("same");
  l2->Draw("same");
  TCanvas *cPullsSVPt = new TCanvas("cPullsSVPt","",800,800);
  hPullsPtKF->GetYaxis()->SetRangeUser(0.5,1.5);
  hPullsPtKF->Draw();
  hPullsPtKFtopo->Draw("same");
  l2->Draw("same");
  TCanvas *cPullsDecL = new TCanvas("cPullsDecL","",800,800);
  hPullsDecLKFtopo->GetYaxis()->SetRangeUser(0.5,2.);
  hPullsDecLKFtopo->Draw();
  hPullsDecLESD->Draw("same");
  ld->Draw("same");
  TCanvas *cPullsImpPar = new TCanvas("cPullsImpPar","",800,800);
  hPullsImpParKF->GetYaxis()->SetRangeUser(0.,1.8);
  hPullsImpParKF->Draw();
  hPullsImpParKFtopo->Draw("same");
  l2->Draw("same");

  cResSVX->SaveAs("ResSVX.eps");
  cResSVXZ->SaveAs("ResSVXZ.eps");
  cResSVY->SaveAs("ResSVY.eps");
  cResSVZ->SaveAs("ResSVZ.eps");
  cResPVX->SaveAs("ResPVX.eps");
  cResPVY->SaveAs("ResPVY.eps");
  cResPVZ->SaveAs("ResPVZ.eps");
  cResPVXZ->SaveAs("ResPVXZ.eps");
  cResSVM->SaveAs("ResM.eps");
  cResSVPt->SaveAs("ResPt.eps");
  cResDecL->SaveAs("ResDecL.eps");
  cResImpPar->SaveAs("ResImpPar.eps");
  cPullsSVX->SaveAs("PullsSVX.eps");
  cPullsSVXZ->SaveAs("PullsSVXZ.eps");
  cPullsSVY->SaveAs("PullsSVY.eps");
  cPullsSVZ->SaveAs("PullsSVZ.eps");
  cPullsPVX->SaveAs("PullsPVX.eps");
  cPullsPVXZ->SaveAs("PullsPVXZ.eps");
  cPullsPVY->SaveAs("PullsPVY.eps");
  cPullsPVZ->SaveAs("PullsPVZ.eps");
  cPullsSVM->SaveAs("PullsM.eps");
  cPullsSVPt->SaveAs("PullsPt.eps");
  cPullsDecL->SaveAs("PullsDecL.eps");
  cPullsImpPar->SaveAs("PullsImpPar.eps");
  cChi2->SaveAs("ChiSquare.eps");
  cProb->SaveAs("Probability.eps");
  
  TFile outfile2(outfilename.Data(),"UPDATE");
  outfile2.cd();
  hResXKF->Write();
  hResYKF->Write();
  hResZKF->Write();
  hResPVXKF->Write();
  hResPVYKF->Write();
  hResPVZKF->Write();
  hResMKF->Write();
  hResPtKF->Write();
  hResImpParKF->Write();
  hResXKFtopo->Write();
  hResYKFtopo->Write();
  hResZKFtopo->Write();
  hResPVXKFtopo->Write();
  hResPVYKFtopo->Write();
  hResPVZKFtopo->Write();
  hResMKFtopo->Write();
  hResPtKFtopo->Write();
  hResDecLKFtopo->Write();
  hResImpParKFtopo->Write();
  hResXESD->Write();
  hResYESD->Write();
  hResZESD->Write();
  hResMESD->Write();
  hResPtESD->Write();
  hResDecLESD->Write();
  hResImpParESD->Write();
  cResSVX->Write();
  cResSVY->Write();
  cResSVZ->Write();
  cResPVX->Write();
  cResPVY->Write();
  cResPVZ->Write();
  cResSVM->Write();
  cResSVPt->Write();
  cResDecL->Write();
  cResImpPar->Write();
  hPullsXKF->Write();
  hPullsYKF->Write();
  hPullsZKF->Write();
  hPullsPVXKF->Write();
  hPullsPVYKF->Write();
  hPullsPVZKF->Write();
  hPullsMKF->Write();
  hPullsPtKF->Write();
  hPullsXKFtopo->Write();
  hPullsYKFtopo->Write();
  hPullsZKFtopo->Write();
  hPullsPVXKFtopo->Write();
  hPullsPVYKFtopo->Write();
  hPullsPVZKFtopo->Write();
  hPullsMKFtopo->Write();
  hPullsPtKFtopo->Write();
  hPullsDecLKFtopo->Write();
  hPullsXESD->Write();
  hPullsYESD->Write();
  hPullsZESD->Write();
  hPullsImpParKF->Write();
  hPullsImpParKFtopo->Write();
  cPullsSVX->Write();
  cPullsSVY->Write();
  cPullsSVZ->Write();
  cPullsPVX->Write();
  cPullsPVY->Write();
  cPullsPVZ->Write();
  cPullsSVM->Write();
  cPullsSVPt->Write();
  cPullsDecL->Write();
  cPullsImpPar->Write();
  cProb->Write();
  cChi2->Write();
  outfile2.Close();

  delete cResSVX;
  delete cResSVY;
  delete cResSVZ;
  delete cResPVX;
  delete cResPVY;
  delete cResPVZ;
  delete cResSVM;
  delete cResSVPt;
  delete cResDecL;
  delete cResImpPar;
  delete cPullsSVX;
  delete cPullsSVY;
  delete cPullsSVZ;
  delete cPullsPVX;
  delete cPullsPVY;
  delete cPullsPVZ;
  delete cPullsSVM;
  delete cPullsSVPt;
  delete cPullsDecL;
  delete cPullsImpPar;
  delete cChi2;
  delete cProb;

}

TH1F* GetHistoVsPt(THnSparseF* sparse, Int_t axis, Int_t unityfactor, TString outfilename, Bool_t noFIT, Bool_t limits, Bool_t print, Bool_t singlegauss) {

  TH1F* hVsPt = new TH1F("hVsPt","",nPtBins,PtLims);
  hVsPt->SetDirectory(0);

  TFile outfile(outfilename.Data(),"UPDATE");
  TCanvas *c = new TCanvas();
  TF1* f = 0x0;
  TF1* f2 = 0x0;

  if(!singlegauss) {
    f = new TF1("f",DoubleGauss,-0.06,0.06,5);
    f->SetParameter(1,0.);
    f->SetParameter(2,0.001);
    f->SetParameter(3,0.1);
    f->SetParLimits(2,0,1000000.);
    f->SetParLimits(3,0,1000000.);
    f->SetParLimits(4,0.,1.);
    f->SetParameter(4,0.7);
 }
  else {
    f = new TF1("f","gaus",-200,200);
    f->SetParameter(1,0.);
    f->SetParameter(2,1.);
  }
  TString fitoption="";
  if(!singlegauss)
    fitoption="LEM";
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ResetAxes(sparse);
    SetPtRange(sparse,PtLims[iPt],PtLims[iPt+1]);

    TH1F* h = (TH1F*)sparse->Projection(axis);
    h->Sumw2();
    Int_t nbins = h->GetNbinsX();
    Double_t min = h->GetBinLowEdge(1)*unityfactor;
    Double_t max = (h->GetBinLowEdge(nbins)+h->GetBinWidth(nbins))*unityfactor;
    TH1F* hCopy = new TH1F("hCopy","",nbins,min,max);
    hCopy->SetName(h->GetName());
    for(Int_t iBin=0; iBin<nbins; iBin++) {
      hCopy->SetBinContent(iBin+1,h->GetBinContent(iBin+1));
      hCopy->SetBinError(iBin+1,h->GetBinError(iBin+1));        
    }
    hCopy->SetFillColor(kBlue);
    hCopy->SetFillStyle(3004);
    hCopy->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",PtLims[iPt],PtLims[iPt+1]));
    TString xtitle = h->GetXaxis()->GetTitle();
    xtitle.ReplaceAll("#sigma(","(");
    xtitle.ReplaceAll("cm","#mum");
    hCopy->GetXaxis()->SetTitle(xtitle.Data());
    hCopy->GetYaxis()->SetTitle("Entries");
    hCopy->GetXaxis()->SetNdivisions(508);
  
    if(!noFIT) {
      if(!singlegauss) f->FixParameter(0,h->GetEntries()*hCopy->GetBinWidth(2));
      if(hCopy->GetEntries()!=0) {
        c->Clear();
        c->Update();
        hCopy->Fit("f",fitoption.Data());
        if(limits)
          hCopy->Fit("f",fitoption.Data(),"",-200,200);
        Double_t chi = f->GetChisquare()/f->GetNDF();
        Bool_t badfit=kFALSE;
        if(!singlegauss && (f->GetParError(2)>f->GetParameter(2) || f->GetParError(3)>f->GetParameter(3))) {
          f2 = new TF1("f2","gaus",-200,200);
          f2->SetParameter(1,0.);
          f2->SetParameter(2,1.);
          hCopy->Fit("f2");
          if(limits)
            hCopy->Fit("f2","","",-200,200);
         badfit=kTRUE;
        }
        if(!singlegauss && chi>3) {
          f2 = new TF1("f2","gaus",-200,200);
          f2->SetParameter(1,0.);
          f2->SetParameter(2,1.);
          hCopy->Fit("f2");
          Double_t chi2 = f2->GetChisquare()/f2->GetNDF();
          if(chi2<chi)
            badfit=kTRUE;
          else
            hCopy->Fit("f",fitoption.Data());
        }
        outfile.cd();
        hCopy->Write();
        Double_t res = f->GetParameter(2);
        Double_t reserr = f->GetParError(2);
        if(!singlegauss && f->GetParameter(2)>f->GetParameter(3)) {
          res = f->GetParameter(3);
          reserr = f->GetParError(3);
      }
        if(badfit) {
          res = f2->GetParameter(2);
          reserr = f2->GetParError(2);
        }
        hVsPt->SetBinContent(iPt+1,res); 
        hVsPt->SetBinError(iPt+1,reserr); 
      }
    }
    else {
      hVsPt->SetBinContent(iPt+1,hCopy->GetRMS()*unityfactor);
      hVsPt->SetBinError(iPt+1,hCopy->GetRMSError()*unityfactor); 
      outfile.cd();
      hCopy->Write();
    }
    if(print && iPt==2) {
      TCanvas* cDist = new TCanvas("cDist","",1000,800);
      gPad->SetRightMargin(0.12);
      hCopy->Draw("hist");
      if(!singlegauss && (f->GetParError(2)>f->GetParameter(2) || f->GetParError(3)>f->GetParameter(3))) 
        f2->Draw("same");
      else
        f->Draw("same");

      cDist->SaveAs("Dist.eps");
      delete cDist;
    }
  }
  
  delete c;
  outfile.Close();
  
  hVsPt->GetXaxis()->SetTitle(sparse->GetAxis(sparse->GetNdimensions()-1)->GetTitle());
  TString ytitle = sparse->GetAxis(axis)->GetTitle();
  ytitle.ReplaceAll("cm","#mum");
  cout << ytitle << endl;
  hVsPt->GetYaxis()->SetTitle(ytitle.Data());
  
  return hVsPt;
}

void SetPtRange(THnSparseF* sparse, Double_t ptmin, Double_t ptmax) {
  TAxis* ptax = (TAxis*)sparse->GetAxis(sparse->GetNdimensions()-1);
  Double_t binmin = ptax->FindBin(ptmin*1.0001);
  Double_t binmax = ptax->FindBin(ptmax*0.9999);
  sparse->GetAxis(sparse->GetNdimensions()-1)->SetRange(binmin,binmax);
}

void ResetAxes(THnSparseF* sparse) {
  Int_t nAxes = sparse->GetNdimensions();
  for(Int_t iAxis=0; iAxis<nAxes; iAxis++) {
    sparse->GetAxis(iAxis)->SetRange(-1,-1);
  }
}

Double_t DoubleGauss(Double_t *x, Double_t *pars) {

  Double_t xx = x[0];
  Double_t norm = pars[0];
  Double_t mean = pars[1];
  Double_t sigma1 = pars[2];
  Double_t sigma2 = pars[3];
  Double_t frac1= pars[4];

  return norm*(frac1*TMath::Gaus(xx,mean,sigma1,kTRUE)+(1-frac1)*TMath::Gaus(xx,mean,sigma2,kTRUE));
}
