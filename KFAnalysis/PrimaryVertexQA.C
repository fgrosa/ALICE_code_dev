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
#include <TGaxis.h>

#endif

enum {kDzero,kDplus};
enum {kRMS,kMean};
const Int_t colors[] = {kBlack,kBlue,kRed,kGreen+3,kMagenta,kCyan,kOrange+7,kYellow,kGreen};
Int_t Nmin=5;
Int_t Nmax=55;
Int_t nbins=10;

void PrimaryVertexQA(Int_t meson=kDplus);
void ResetAxes(THnSparseF* sparse);
TH1F* GetHistoVsNtracks(THnSparseF* sparse, Int_t axis, TString outfilename, Int_t quantity=kRMS);

void PrimaryVertexQA(Int_t meson) {

  TString mesonname;
  if(meson==kDzero) mesonname = "Dzero";
  else mesonname = "Dplus";

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetLegendBorderSize(0);

  TGaxis::SetMaxDigits(2);
    
  //INPUT FILES________________________________________________________________
  TString filename;
  TString listname;
  
  if(meson==kDzero) {
    filename = "$HOME/ALICE/Files/TrainsGSI/Run1/LHC13/fgrosa_Dplus_KF_pPbMC_DLXY.root";
    listname = "coutputDzeroKF";
  }
  else {
    filename = "$HOME/ALICE/Files/TrainsGSI/Run1/LHC13/fgrosa_Dplus_KF_pPbMC_DLXY.root";
    listname = "coutputDplusKF";
  }
  
  TFile infile(filename.Data(),"UPDATE");
  TList* list = (TList*)infile.Get(listname.Data());
  
  //PRIMARY VERTEX QA__________________________________________________________
  THnSparseF* hPVHFres = (THnSparseF*)list->FindObject("fHistPVHFres"); 
  THnSparseF* hPVHFpulls = (THnSparseF*)list->FindObject("fHistPVHFpulls"); 
  TH1F* hPVHFchiS = (TH1F*)list->FindObject("fHistPVHFchiS");
  TH1F* hPVHFprob = (TH1F*)list->FindObject("fHistPVHFprob");
  hPVHFchiS->SetDirectory(0);
  hPVHFprob->SetDirectory(0);
  hPVHFchiS->SetLineColor(colors[0]);
  hPVHFprob->SetLineColor(colors[0]);
  hPVHFchiS->GetYaxis()->SetTitle("Entries"); 
  hPVHFprob->GetYaxis()->SetTitle("Entries");
  
  THnSparseF* hPVHFremDaures = (THnSparseF*)list->FindObject("fHistPVHFremDaures"); 
  THnSparseF* hPVHFremDaupulls = (THnSparseF*)list->FindObject("fHistPVHFremDaupulls"); 
  TH1F* hPVHFremDauchiS = (TH1F*)list->FindObject("fHistPVHFremDauchiS");
  TH1F* hPVHFremDauprob = (TH1F*)list->FindObject("fHistPVHFremDauprob");
  hPVHFremDauchiS->SetDirectory(0);
  hPVHFremDauprob->SetDirectory(0);
  hPVHFremDauchiS->SetLineColor(colors[1]);
  hPVHFremDauprob->SetLineColor(colors[1]);
  hPVHFremDauchiS->GetYaxis()->SetTitle("Entries"); 
  hPVHFremDauprob->GetYaxis()->SetTitle("Entries");
 
  THnSparseF* hPVESDHFres = (THnSparseF*)list->FindObject("fHistPVESDHFres"); 
  THnSparseF* hPVESDHFpulls = (THnSparseF*)list->FindObject("fHistPVESDHFpulls"); 
  TH1F* hPVESDHFchiS = (TH1F*)list->FindObject("fHistPVESDHFchiS");
  TH1F* hPVESDHFprob = (TH1F*)list->FindObject("fHistPVESDHFprob");
  hPVESDHFchiS->SetDirectory(0);
  hPVESDHFprob->SetDirectory(0);
  hPVESDHFchiS->SetLineColor(colors[2]);
  hPVESDHFprob->SetLineColor(colors[2]);
  hPVESDHFchiS->GetYaxis()->SetTitle("Entries"); 
  hPVESDHFprob->GetYaxis()->SetTitle("Entries");
  
  THnSparseF* hPVESDHFremDaures = (THnSparseF*)list->FindObject("fHistPVESDHFremDaures"); 
  THnSparseF* hPVESDHFremDaupulls = (THnSparseF*)list->FindObject("fHistPVESDHFremDaupulls"); 
  TH1F* hPVESDHFremDauchiS = (TH1F*)list->FindObject("fHistPVESDHFremDauchiS");
  TH1F* hPVESDHFremDauprob = (TH1F*)list->FindObject("fHistPVESDHFremDauprob");
  hPVESDHFremDauchiS->SetDirectory(0);
  hPVESDHFremDauprob->SetDirectory(0);
  hPVESDHFremDauchiS->SetLineColor(colors[3]);
  hPVESDHFremDauprob->SetLineColor(colors[3]);
  hPVESDHFremDauchiS->GetYaxis()->SetTitle("Entries"); 
  hPVESDHFremDauprob->GetYaxis()->SetTitle("Entries");

  TString outfilename = Form("PVHF_%s_checks.root",mesonname.Data());
  TFile outfile(outfilename.Data(),"RECREATE");
  outfile.Close();
  
  TH1F* hPVHFresXVsN = (TH1F*)GetHistoVsNtracks(hPVHFres,0,outfilename);
  TH1F* hPVHFpullsXVsN = (TH1F*)GetHistoVsNtracks(hPVHFpulls,0,outfilename);
  TH1F* hPVHFresYVsN = (TH1F*)GetHistoVsNtracks(hPVHFres,1,outfilename);
  TH1F* hPVHFpullsYVsN = (TH1F*)GetHistoVsNtracks(hPVHFpulls,1,outfilename);
  TH1F* hPVHFresZVsN = (TH1F*)GetHistoVsNtracks(hPVHFres,2,outfilename);
  TH1F* hPVHFpullsZVsN = (TH1F*)GetHistoVsNtracks(hPVHFpulls,2,outfilename);
  hPVHFresXVsN->SetName("hPVHFresXVsN");
  hPVHFpullsXVsN->SetName("hPVHFpullsXVsN");
  hPVHFresYVsN->SetName("hPVHFresYVsN");
  hPVHFpullsYVsN->SetName("hPVHFpullsYVsN");
  hPVHFresZVsN->SetName("hPVHFresZVsN");
  hPVHFpullsZVsN->SetName("hPVHFpullsZVsN");
  TH1F* hPVHFremDauresXVsN = (TH1F*)GetHistoVsNtracks(hPVHFremDaures,0,outfilename);
  TH1F* hPVHFremDaupullsXVsN = (TH1F*)GetHistoVsNtracks(hPVHFremDaupulls,0,outfilename);
  TH1F* hPVHFremDauresYVsN = (TH1F*)GetHistoVsNtracks(hPVHFremDaures,1,outfilename);
  TH1F* hPVHFremDaupullsYVsN = (TH1F*)GetHistoVsNtracks(hPVHFremDaupulls,1,outfilename);
  TH1F* hPVHFremDauresZVsN = (TH1F*)GetHistoVsNtracks(hPVHFremDaures,2,outfilename);
  TH1F* hPVHFremDaupullsZVsN = (TH1F*)GetHistoVsNtracks(hPVHFremDaupulls,2,outfilename);
  hPVHFremDauresXVsN->SetName("hPVHFremDauresXVsN");
  hPVHFremDaupullsXVsN->SetName("hPVHFremDaupullsXVsN");
  hPVHFremDauresYVsN->SetName("hPVHFremDauresYVsN");
  hPVHFremDaupullsYVsN->SetName("hPVHFremDaupullsYVsN");
  hPVHFremDauresZVsN->SetName("hPVHFremDauresZVsN");
  hPVHFremDaupullsZVsN->SetName("hPVHFremDaupullsZVsN");
  TH1F* hPVESDHFresXVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFres,0,outfilename);
  TH1F* hPVESDHFpullsXVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFpulls,0,outfilename);
  TH1F* hPVESDHFresYVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFres,1,outfilename);
  TH1F* hPVESDHFpullsYVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFpulls,1,outfilename);
  TH1F* hPVESDHFresZVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFres,2,outfilename);
  TH1F* hPVESDHFpullsZVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFpulls,2,outfilename);
  hPVESDHFresXVsN->SetName("hPVESDHFresXVsN");
  hPVESDHFpullsXVsN->SetName("hPVESDHFpullsXVsN");
  hPVESDHFresYVsN->SetName("hPVESDHFresYVsN");
  hPVESDHFpullsYVsN->SetName("hPVESDHFpullsYVsN");
  hPVESDHFresZVsN->SetName("hPVESDHFresZVsN");
  hPVESDHFpullsZVsN->SetName("hPVESDHFpullsZVsN");
  TH1F* hPVESDHFremDauresXVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFremDaures,0,outfilename);
  TH1F* hPVESDHFremDaupullsXVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFremDaupulls,0,outfilename);
  TH1F* hPVESDHFremDauresYVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFremDaures,1,outfilename);
  TH1F* hPVESDHFremDaupullsYVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFremDaupulls,1,outfilename);
  TH1F* hPVESDHFremDauresZVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFremDaures,2,outfilename);
  TH1F* hPVESDHFremDaupullsZVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFremDaupulls,2,outfilename);
  hPVESDHFremDauresXVsN->SetName("hPVESDHFremDauresXVsN");
  hPVESDHFremDaupullsXVsN->SetName("hPVESDHFremDaupullsXVsN");
  hPVESDHFremDauresYVsN->SetName("hPVESDHFremDauresYVsN");
  hPVESDHFremDaupullsYVsN->SetName("hPVESDHFremDaupullsYVsN");
  hPVESDHFremDauresZVsN->SetName("hPVESDHFremDauresZVsN");
  hPVESDHFremDaupullsZVsN->SetName("hPVESDHFremDaupullsZVsN");
  
  TH1F* hPVHFmeanXVsN = (TH1F*)GetHistoVsNtracks(hPVHFres,0,outfilename,kMean);
  TH1F* hPVHFmeanYVsN = (TH1F*)GetHistoVsNtracks(hPVHFres,1,outfilename,kMean);
  TH1F* hPVHFmeanZVsN = (TH1F*)GetHistoVsNtracks(hPVHFres,2,outfilename,kMean);
  hPVHFmeanXVsN->SetName("hPVHFmeanXVsN");
  hPVHFmeanYVsN->SetName("hPVHFmeanYVsN");
  hPVHFmeanZVsN->SetName("hPVHFmeanZVsN");
  TH1F* hPVHFremDaumeanXVsN = (TH1F*)GetHistoVsNtracks(hPVHFremDaures,0,outfilename,kMean);
  TH1F* hPVHFremDaumeanYVsN = (TH1F*)GetHistoVsNtracks(hPVHFremDaures,1,outfilename,kMean);
  TH1F* hPVHFremDaumeanZVsN = (TH1F*)GetHistoVsNtracks(hPVHFremDaures,2,outfilename,kMean);
  hPVHFremDaumeanXVsN->SetName("hPVHFremDaumeanXVsN");
  hPVHFremDaumeanYVsN->SetName("hPVHFremDaumeanYVsN");
  hPVHFremDaumeanZVsN->SetName("hPVHFremDaumeanZVsN");
  TH1F* hPVESDHFmeanXVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFres,0,outfilename,kMean);
  TH1F* hPVESDHFmeanYVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFres,1,outfilename,kMean);
  TH1F* hPVESDHFmeanZVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFres,2,outfilename,kMean);
  hPVESDHFmeanXVsN->SetName("hPVESDHFmeanXVsN");
  hPVESDHFmeanYVsN->SetName("hPVESDHFmeanYVsN");
  hPVESDHFmeanZVsN->SetName("hPVESDHFmeanZVsN");
  TH1F* hPVESDHFremDaumeanXVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFremDaures,0,outfilename,kMean);
  TH1F* hPVESDHFremDaumeanYVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFremDaures,1,outfilename,kMean);
  TH1F* hPVESDHFremDaumeanZVsN = (TH1F*)GetHistoVsNtracks(hPVESDHFremDaures,2,outfilename,kMean);
  hPVESDHFremDaumeanXVsN->SetName("hPVESDHFremDaumeanXVsN");
  hPVESDHFremDaumeanYVsN->SetName("hPVESDHFremDaumeanYVsN");
  hPVESDHFremDaumeanZVsN->SetName("hPVESDHFremDaumeanZVsN");
  
  hPVHFresXVsN->SetLineColor(colors[0]);
  hPVHFpullsXVsN->SetLineColor(colors[0]);
  hPVHFresYVsN->SetLineColor(colors[0]);
  hPVHFpullsYVsN->SetLineColor(colors[0]);
  hPVHFresZVsN->SetLineColor(colors[0]);
  hPVHFpullsZVsN->SetLineColor(colors[0]);
  hPVHFremDauresXVsN->SetLineColor(colors[1]);
  hPVHFremDaupullsXVsN->SetLineColor(colors[1]);
  hPVHFremDauresYVsN->SetLineColor(colors[1]);
  hPVHFremDaupullsYVsN->SetLineColor(colors[1]);
  hPVHFremDauresZVsN->SetLineColor(colors[1]);
  hPVHFremDaupullsZVsN->SetLineColor(colors[1]);
  hPVESDHFresXVsN->SetLineColor(colors[2]);
  hPVESDHFpullsXVsN->SetLineColor(colors[2]);
  hPVESDHFresYVsN->SetLineColor(colors[2]);
  hPVESDHFpullsYVsN->SetLineColor(colors[2]);
  hPVESDHFresZVsN->SetLineColor(colors[2]);
  hPVESDHFpullsZVsN->SetLineColor(colors[2]);
  hPVESDHFremDauresXVsN->SetLineColor(colors[3]);
  hPVESDHFremDaupullsXVsN->SetLineColor(colors[3]);
  hPVESDHFremDauresYVsN->SetLineColor(colors[3]);
  hPVESDHFremDaupullsYVsN->SetLineColor(colors[3]);
  hPVESDHFremDauresZVsN->SetLineColor(colors[3]);
  hPVESDHFremDaupullsZVsN->SetLineColor(colors[3]);
  hPVHFresXVsN->SetMarkerColor(colors[0]);
  hPVHFpullsXVsN->SetMarkerColor(colors[0]);
  hPVHFresYVsN->SetMarkerColor(colors[0]);
  hPVHFpullsYVsN->SetMarkerColor(colors[0]);
  hPVHFresZVsN->SetMarkerColor(colors[0]);
  hPVHFpullsZVsN->SetMarkerColor(colors[0]);
  hPVHFremDauresXVsN->SetMarkerColor(colors[1]);
  hPVHFremDaupullsXVsN->SetMarkerColor(colors[1]);
  hPVHFremDauresYVsN->SetMarkerColor(colors[1]);
  hPVHFremDaupullsYVsN->SetMarkerColor(colors[1]);
  hPVHFremDauresZVsN->SetMarkerColor(colors[1]);
  hPVHFremDaupullsZVsN->SetMarkerColor(colors[1]);
  hPVESDHFresXVsN->SetMarkerColor(colors[2]);
  hPVESDHFpullsXVsN->SetMarkerColor(colors[2]);
  hPVESDHFresYVsN->SetMarkerColor(colors[2]);
  hPVESDHFpullsYVsN->SetMarkerColor(colors[2]);
  hPVESDHFresZVsN->SetMarkerColor(colors[2]);
  hPVESDHFpullsZVsN->SetMarkerColor(colors[2]);
  hPVESDHFremDauresXVsN->SetMarkerColor(colors[3]);
  hPVESDHFremDaupullsXVsN->SetMarkerColor(colors[3]);
  hPVESDHFremDauresYVsN->SetMarkerColor(colors[3]);
  hPVESDHFremDaupullsYVsN->SetMarkerColor(colors[3]);
  hPVESDHFremDauresZVsN->SetMarkerColor(colors[3]);
  hPVESDHFremDaupullsZVsN->SetMarkerColor(colors[3]);
  hPVHFresXVsN->SetMarkerStyle(20);
  hPVHFpullsXVsN->SetMarkerStyle(20);
  hPVHFresYVsN->SetMarkerStyle(20);
  hPVHFpullsYVsN->SetMarkerStyle(20);
  hPVHFresZVsN->SetMarkerStyle(20);
  hPVHFpullsZVsN->SetMarkerStyle(20);
  hPVHFremDauresXVsN->SetMarkerStyle(20);
  hPVHFremDaupullsXVsN->SetMarkerStyle(20);
  hPVHFremDauresYVsN->SetMarkerStyle(20);
  hPVHFremDaupullsYVsN->SetMarkerStyle(20);
  hPVHFremDauresZVsN->SetMarkerStyle(20);
  hPVHFremDaupullsZVsN->SetMarkerStyle(20);
  hPVESDHFresXVsN->SetMarkerStyle(20);
  hPVESDHFpullsXVsN->SetMarkerStyle(20);
  hPVESDHFresYVsN->SetMarkerStyle(20);
  hPVESDHFpullsYVsN->SetMarkerStyle(20);
  hPVESDHFresZVsN->SetMarkerStyle(20);
  hPVESDHFpullsZVsN->SetMarkerStyle(20);
  hPVESDHFremDauresXVsN->SetMarkerStyle(20);
  hPVESDHFremDaupullsXVsN->SetMarkerStyle(20);
  hPVESDHFremDauresYVsN->SetMarkerStyle(20);
  hPVESDHFremDaupullsYVsN->SetMarkerStyle(20);
  hPVESDHFremDauresZVsN->SetMarkerStyle(20);
  hPVESDHFremDaupullsZVsN->SetMarkerStyle(20);
  hPVHFresXVsN->SetMarkerSize(1.5);
  hPVHFpullsXVsN->SetMarkerSize(1.5);
  hPVHFresYVsN->SetMarkerSize(1.5);
  hPVHFpullsYVsN->SetMarkerSize(1.5);
  hPVHFresZVsN->SetMarkerSize(1.5);
  hPVHFpullsZVsN->SetMarkerSize(1.5);
  hPVHFremDauresXVsN->SetMarkerSize(1.5);
  hPVHFremDaupullsXVsN->SetMarkerSize(1.5);
  hPVHFremDauresYVsN->SetMarkerSize(1.5);
  hPVHFremDaupullsYVsN->SetMarkerSize(1.5);
  hPVHFremDauresZVsN->SetMarkerSize(1.5);
  hPVHFremDaupullsZVsN->SetMarkerSize(1.5);
  hPVESDHFresXVsN->SetMarkerSize(1.5);
  hPVESDHFpullsXVsN->SetMarkerSize(1.5);
  hPVESDHFresYVsN->SetMarkerSize(1.5);
  hPVESDHFpullsYVsN->SetMarkerSize(1.5);
  hPVESDHFresZVsN->SetMarkerSize(1.5);
  hPVESDHFpullsZVsN->SetMarkerSize(1.5);
  hPVESDHFremDauresXVsN->SetMarkerSize(1.5);
  hPVESDHFremDaupullsXVsN->SetMarkerSize(1.5);
  hPVESDHFremDauresYVsN->SetMarkerSize(1.5);
  hPVESDHFremDaupullsYVsN->SetMarkerSize(1.5);
  hPVESDHFremDauresZVsN->SetMarkerSize(1.5);
  hPVESDHFremDaupullsZVsN->SetMarkerSize(1.5);
  
  hPVHFmeanXVsN->SetLineColor(colors[0]);
  hPVHFmeanYVsN->SetLineColor(colors[0]);
  hPVHFmeanZVsN->SetLineColor(colors[0]);
  hPVHFremDaumeanXVsN->SetLineColor(colors[1]);
  hPVHFremDaumeanYVsN->SetLineColor(colors[1]);
  hPVHFremDaumeanZVsN->SetLineColor(colors[1]);
  hPVESDHFmeanXVsN->SetLineColor(colors[2]);
  hPVESDHFmeanYVsN->SetLineColor(colors[2]);
  hPVESDHFmeanZVsN->SetLineColor(colors[2]);
  hPVESDHFremDaumeanXVsN->SetLineColor(colors[3]);
  hPVESDHFremDaumeanYVsN->SetLineColor(colors[3]);
  hPVESDHFremDaumeanZVsN->SetLineColor(colors[3]);
  hPVHFmeanXVsN->SetMarkerColor(colors[0]);
  hPVHFmeanYVsN->SetMarkerColor(colors[0]);
  hPVHFmeanZVsN->SetMarkerColor(colors[0]);
  hPVHFremDaumeanXVsN->SetMarkerColor(colors[1]);
  hPVHFremDaumeanYVsN->SetMarkerColor(colors[1]);
  hPVHFremDaumeanZVsN->SetMarkerColor(colors[1]);
  hPVESDHFmeanXVsN->SetMarkerColor(colors[2]);
  hPVESDHFmeanYVsN->SetMarkerColor(colors[2]);
  hPVESDHFmeanZVsN->SetMarkerColor(colors[2]);
  hPVESDHFremDaumeanXVsN->SetMarkerColor(colors[3]);
  hPVESDHFremDaumeanYVsN->SetMarkerColor(colors[3]);
  hPVESDHFremDaumeanZVsN->SetMarkerColor(colors[3]);
  hPVHFmeanXVsN->SetMarkerStyle(20);
  hPVHFmeanYVsN->SetMarkerStyle(20);
  hPVHFmeanZVsN->SetMarkerStyle(20);
  hPVHFremDaumeanXVsN->SetMarkerStyle(20);
  hPVHFremDaumeanYVsN->SetMarkerStyle(20);
  hPVHFremDaumeanZVsN->SetMarkerStyle(20);
  hPVESDHFmeanXVsN->SetMarkerStyle(20);
  hPVESDHFmeanYVsN->SetMarkerStyle(20);
  hPVESDHFmeanZVsN->SetMarkerStyle(20);
  hPVESDHFremDaumeanXVsN->SetMarkerStyle(20);
  hPVESDHFremDaumeanYVsN->SetMarkerStyle(20);
  hPVESDHFremDaumeanZVsN->SetMarkerStyle(20);
  hPVHFmeanXVsN->SetMarkerSize(1.5);
  hPVHFmeanYVsN->SetMarkerSize(1.5);
  hPVHFmeanZVsN->SetMarkerSize(1.5);
  hPVHFremDaumeanXVsN->SetMarkerSize(1.5);
  hPVHFremDaumeanYVsN->SetMarkerSize(1.5);
  hPVHFremDaumeanZVsN->SetMarkerSize(1.5);
  hPVESDHFmeanXVsN->SetMarkerSize(1.5);
  hPVESDHFmeanYVsN->SetMarkerSize(1.5);
  hPVESDHFmeanZVsN->SetMarkerSize(1.5);
  hPVESDHFremDaumeanXVsN->SetMarkerSize(1.5);
  hPVESDHFremDaumeanYVsN->SetMarkerSize(1.5);
  hPVESDHFremDaumeanZVsN->SetMarkerSize(1.5);
  
  hPVHFresXVsN->GetYaxis()->SetRangeUser(0.,0.008);
  hPVHFresYVsN->GetYaxis()->SetRangeUser(0.,0.008);
  hPVHFresZVsN->GetYaxis()->SetRangeUser(0.,0.016);
  hPVHFpullsXVsN->GetYaxis()->SetRangeUser(0.5,2.);
  hPVHFpullsYVsN->GetYaxis()->SetRangeUser(0.5,2.);
  hPVHFpullsZVsN->GetYaxis()->SetRangeUser(0.5,2.);
    
  TLegend* leg = new TLegend(0.25,0.7,0.85,0.89);
  leg->SetTextSize(0.045);
  leg->AddEntry(hPVESDHFresXVsN,"AliVertexer","l");
  leg->AddEntry(hPVESDHFremDauresYVsN,"AliVertexer rem dau","l");
  leg->AddEntry(hPVHFremDauresYVsN,"KF rem dau & add mother","l");

  TLegend* leg2 = new TLegend(0.25,0.7,0.85,0.89);
  leg2->SetTextSize(0.045);
  leg2->AddEntry(hPVESDHFresXVsN,"AliVertexer","lpe");
  leg2->AddEntry(hPVESDHFremDauresYVsN,"AliVertexer rem dau","lpe");
  leg2->AddEntry(hPVHFremDauresYVsN,"KF rem dau & add mother","lpe");

  TLatex* latKF = new TLatex();
  latKF->SetTextFont(132);
  latKF->SetTextSize(0.035);
  latKF->SetTextColor(colors[0]);
  TLatex* latKFremDau = new TLatex();
  latKFremDau->SetTextFont(132);
  latKFremDau->SetTextSize(0.035);
  latKFremDau->SetTextColor(colors[1]);
  TLatex* latESD = new TLatex();
  latESD->SetTextFont(132);
  latESD->SetTextSize(0.035);
  latESD->SetTextColor(colors[2]);
  TLatex* latESDremDau = new TLatex();
  latESDremDau->SetTextFont(132);
  latESDremDau->SetTextSize(0.035);
  latESDremDau->SetTextColor(colors[3]);
  
  TCanvas* cPVHFchiS = new TCanvas("cPVHFchiS","cPVHFchiS",800,800);
  cPVHFchiS->SetLogy();
  hPVHFchiS->Draw();
  hPVHFremDauchiS->Draw("same");
  hPVESDHFchiS->Draw("same");
  hPVESDHFremDauchiS->Draw("same");
  leg->Draw("same");
  latKF->DrawLatex(12,hPVHFchiS->GetMaximum()*0.05,Form("mean = %0.2f #pm %0.2f",hPVHFchiS->GetMean(),hPVHFchiS->GetMeanError()));
  latKFremDau->DrawLatex(12,hPVHFchiS->GetMaximum()*0.025,Form("mean = %0.2f #pm %0.2f",hPVHFremDauchiS->GetMean(),hPVHFremDauchiS->GetMeanError()));
  latESD->DrawLatex(12,hPVHFchiS->GetMaximum()*0.010,Form("mean = %0.2f #pm %0.2f",hPVESDHFchiS->GetMean(),hPVESDHFchiS->GetMeanError()));
  latESDremDau->DrawLatex(12,hPVHFchiS->GetMaximum()*0.005,Form("mean = %0.2f #pm %0.2f",hPVESDHFremDauchiS->GetMean(),hPVESDHFremDauchiS->GetMeanError()));
  
  TCanvas* cPVHFprob = new TCanvas("cPVHFprob","cPVHFprob",800,800);
  cPVHFprob->SetLogy();
  hPVHFprob->GetYaxis()->SetRangeUser(100,hPVHFprob->GetMaximum()*5.);
  hPVHFprob->Draw();
  hPVHFremDauprob->Draw("same");
  hPVESDHFprob->Draw("same");
  hPVESDHFremDauprob->Draw("same");
  leg->Draw("same");
  latKF->DrawLatex(0.5,hPVHFremDauprob->GetMaximum()*0.5,Form("mean = %0.2f #pm %0.2f",hPVHFprob->GetMean(),hPVHFprob->GetMeanError()));
  latKFremDau->DrawLatex(0.5,hPVHFremDauprob->GetMaximum()*0.4,Form("mean = %0.2f #pm %0.2f",hPVHFremDauprob->GetMean(),hPVHFremDauprob->GetMeanError()));
  latESD->DrawLatex(0.5,hPVHFremDauprob->GetMaximum()*0.32,Form("mean = %0.2f #pm %0.2f",hPVESDHFprob->GetMean(),hPVESDHFprob->GetMeanError()));
  latESDremDau->DrawLatex(0.5,hPVHFremDauprob->GetMaximum()*0.24,Form("mean = %0.2f #pm %0.2f",hPVESDHFremDauprob->GetMean(),hPVESDHFremDauprob->GetMeanError()));
  
  TCanvas* cPVHFresX = new TCanvas("cPVHFresX","cPVHFresX",800,800);
  hPVHFresXVsN->Draw();
  hPVHFremDauresXVsN->Draw("same");
  hPVESDHFresXVsN->Draw("same");
  hPVESDHFremDauresXVsN->Draw("same");
  leg2->Draw("same");
  TCanvas* cPVHFpullsX = new TCanvas("cPVHFpullsX","cPVHFpullsX",800,800);
  hPVHFpullsXVsN->Draw();
  hPVHFremDaupullsXVsN->Draw("same");
  hPVESDHFpullsXVsN->Draw("same");
  hPVESDHFremDaupullsXVsN->Draw("same");
  leg2->Draw("same");
  TCanvas* cPVHFresY = new TCanvas("cPVHFresY","cPVHFresY",800,800);
  hPVHFresYVsN->Draw();
  hPVHFremDauresYVsN->Draw("same");
  hPVESDHFresYVsN->Draw("same");
  hPVESDHFremDauresYVsN->Draw("same");
  leg2->Draw("same");
  TCanvas* cPVHFpullsY = new TCanvas("cPVHFpullsY","cPVHFpullsY",800,800);
  hPVHFpullsYVsN->Draw();
  hPVHFremDaupullsYVsN->Draw("same");
  hPVESDHFpullsYVsN->Draw("same");
  hPVESDHFremDaupullsYVsN->Draw("same");
  leg2->Draw("same");
  TCanvas* cPVHFresZ = new TCanvas("cPVHFresZ","cPVHFresZ",800,800);
  hPVHFresZVsN->Draw();
  hPVHFremDauresZVsN->Draw("same");
  hPVESDHFresZVsN->Draw("same");
  hPVESDHFremDauresZVsN->Draw("same");
  leg2->Draw("same");
  TCanvas* cPVHFpullsZ = new TCanvas("cPVHFpullsZ","cPVHFpullsZ",800,800);
  hPVHFpullsZVsN->Draw();
  hPVHFremDaupullsZVsN->Draw("same");
  hPVESDHFpullsZVsN->Draw("same");
  hPVESDHFremDaupullsZVsN->Draw("same");
  leg2->Draw("same");

  TCanvas* cPVHFresXZ = new TCanvas("cPVHFresXZ","cPVHFresXZ",800,800);
  TH1F* hPVHFresXVsNCopy = (TH1F*)hPVHFresXVsN->Clone();
  TH1F* hPVHFresZVsNCopy = (TH1F*)hPVHFresZVsN->Clone();
  TH1F* hPVHFremDauresXVsNCopy = (TH1F*)hPVHFremDauresXVsN->Clone();
  TH1F* hPVHFremDauresZVsNCopy = (TH1F*)hPVHFremDauresZVsN->Clone();
  TH1F* hPVESDHFresXVsNCopy = (TH1F*)hPVESDHFresXVsN->Clone();
  TH1F* hPVESDHFresZVsNCopy = (TH1F*)hPVESDHFresZVsN->Clone();
  TH1F* hPVESDHFremDauresXVsNCopy = (TH1F*)hPVESDHFremDauresXVsN->Clone();
  TH1F* hPVESDHFremDauresZVsNCopy = (TH1F*)hPVESDHFremDauresZVsN->Clone();
  TLegend* lXZ = new TLegend(0.25,0.55,0.89,0.89);
  lXZ->SetTextSize(0.04);
  lXZ->SetFillStyle(0);
  lXZ->AddEntry(hPVESDHFresXVsNCopy,"AliVertexer x","lpe");
  lXZ->AddEntry(hPVESDHFremDauresXVsNCopy,"AliVertexer rem dau x","lpe");
  lXZ->AddEntry(hPVHFremDauresXVsNCopy,"KF rem dau & add mother x","lpe");
  lXZ->AddEntry(hPVESDHFresZVsNCopy,"AliVertexer z","lpe");
  lXZ->AddEntry(hPVESDHFremDauresZVsNCopy,"AliVertexer rem dau z","lpe");
  lXZ->AddEntry(hPVHFremDauresZVsNCopy,"KF rem dau & add mother z","lpe");
  hPVHFresZVsNCopy->SetMarkerStyle(26);
  hPVHFremDauresZVsNCopy->SetMarkerStyle(26);
  hPVESDHFresZVsNCopy->SetMarkerStyle(26);
  hPVESDHFremDauresZVsNCopy->SetMarkerStyle(26);
  hPVHFresXVsNCopy->GetYaxis()->SetRangeUser(0.001,0.015);
  hPVHFresXVsNCopy->GetYaxis()->SetTitle("#sigma(pos_{reco}-pos_{MC}) (#mum)");
  hPVHFresXVsNCopy->Draw();
  hPVHFremDauresXVsNCopy->Draw("same");
  hPVESDHFresXVsNCopy->Draw("same");
  hPVESDHFremDauresXVsNCopy->Draw("same");
  hPVHFresZVsNCopy->Draw("same");
  hPVHFremDauresZVsNCopy->Draw("same");
  hPVESDHFresZVsNCopy->Draw("same");
  hPVESDHFremDauresZVsNCopy->Draw("same");
  lXZ->Draw("same");
  
  TCanvas* cPVHFmeanX = new TCanvas("cPVHFmeanX","cPVHFmeanX",800,800);
  hPVHFmeanXVsN->GetYaxis()->SetRangeUser(-0.0002,0.0005);
  hPVHFmeanXVsN->Draw();
  hPVHFremDaumeanXVsN->Draw("same");
  hPVESDHFmeanXVsN->Draw("same");
  hPVESDHFremDaumeanXVsN->Draw("same");
  leg2->Draw("same");
  TCanvas* cPVHFmeanY = new TCanvas("cPVHFmeanY","cPVHFmeanY",800,800);
  hPVHFmeanYVsN->GetYaxis()->SetRangeUser(-0.0002,0.0005);
  hPVHFmeanYVsN->Draw();
  hPVHFremDaumeanYVsN->Draw("same");
  hPVESDHFmeanYVsN->Draw("same");
  hPVESDHFremDaumeanYVsN->Draw("same");
  leg2->Draw("same");
  TCanvas* cPVHFmeanZ = new TCanvas("cPVHFmeanZ","cPVHFmeanZ",800,800);
  hPVHFmeanZVsN->GetYaxis()->SetRangeUser(-0.0005,0.0002);
  hPVHFmeanZVsN->Draw();
  hPVHFremDaumeanZVsN->Draw("same");
  hPVESDHFmeanZVsN->Draw("same");
  hPVESDHFremDaumeanZVsN->Draw("same");
  leg2->Draw("same");

  TCanvas* cPVHFmeanXZ = new TCanvas("cPVHFmeanXZ","cPVHFmeanXZ",800,800);
  TH1F* hPVHFmeanXVsNCopy = (TH1F*)hPVHFmeanXVsN->Clone();
  TH1F* hPVHFmeanZVsNCopy = (TH1F*)hPVHFmeanZVsN->Clone();
  TH1F* hPVHFremDaumeanXVsNCopy = (TH1F*)hPVHFremDaumeanXVsN->Clone();
  TH1F* hPVHFremDaumeanZVsNCopy = (TH1F*)hPVHFremDaumeanZVsN->Clone();
  TH1F* hPVESDHFmeanXVsNCopy = (TH1F*)hPVESDHFmeanXVsN->Clone();
  TH1F* hPVESDHFmeanZVsNCopy = (TH1F*)hPVESDHFmeanZVsN->Clone();
  TH1F* hPVESDHFremDaumeanXVsNCopy = (TH1F*)hPVESDHFremDaumeanXVsN->Clone();
  TH1F* hPVESDHFremDaumeanZVsNCopy = (TH1F*)hPVESDHFremDaumeanZVsN->Clone();
  hPVHFmeanZVsNCopy->SetMarkerStyle(26);
  hPVHFremDaumeanZVsNCopy->SetMarkerStyle(26);
  hPVESDHFmeanZVsNCopy->SetMarkerStyle(26);
  hPVESDHFremDaumeanZVsNCopy->SetMarkerStyle(26);
  hPVHFmeanXVsNCopy->GetYaxis()->SetRangeUser(-0.0008,0.0015);
  hPVHFmeanXVsNCopy->GetYaxis()->SetTitle("<pos_{reco}-pos_{MC}> (#mum)");
  hPVHFmeanXVsNCopy->Draw();
  hPVHFremDaumeanXVsNCopy->Draw("same");
  hPVESDHFmeanXVsNCopy->Draw("same");
  hPVESDHFremDaumeanXVsNCopy->Draw("same");
  hPVHFmeanZVsNCopy->Draw("same");
  hPVHFremDaumeanZVsNCopy->Draw("same");
  hPVESDHFmeanZVsNCopy->Draw("same");
  hPVESDHFremDaumeanZVsNCopy->Draw("same");
  lXZ->Draw("same");
  
  cPVHFresX->SaveAs("cPVresX.eps");
  cPVHFresY->SaveAs("cPVresY.eps");
  cPVHFresZ->SaveAs("cPVresZ.eps");
  cPVHFmeanX->SaveAs("cPVmeanX.eps");
  cPVHFmeanY->SaveAs("cPVmeanY.eps");
  cPVHFmeanZ->SaveAs("cPVmeanZ.eps");
  cPVHFpullsX->SaveAs("cPVpullsX.eps");
  cPVHFpullsY->SaveAs("cPVpullsY.eps");
  cPVHFpullsZ->SaveAs("cPVpullsZ.eps");
  cPVHFchiS->SaveAs("cPVchiS.eps");
  cPVHFprob->SaveAs("cPVprob.eps");
  cPVHFresXZ->SaveAs("cPVresXZ.eps");
  cPVHFmeanXZ->SaveAs("cPVmeanXZ.eps");

  TFile outfilePVHF(outfilename.Data(),"UPDATE");
  hPVHFchiS->Write();
  hPVHFprob->Write();
  hPVHFresXVsN->Write();
  hPVHFmeanXVsN->Write();
  hPVHFpullsXVsN->Write();
  hPVHFresYVsN->Write();
  hPVHFmeanYVsN->Write();
  hPVHFpullsYVsN->Write();
  hPVHFresZVsN->Write();
  hPVHFmeanZVsN->Write();
  hPVHFpullsZVsN->Write();
  hPVHFremDauchiS->Write();
  hPVHFremDauprob->Write();
  hPVHFremDauresXVsN->Write();
  hPVHFremDaumeanXVsN->Write();
  hPVHFremDaupullsXVsN->Write();
  hPVHFremDauresYVsN->Write();
  hPVHFremDaumeanYVsN->Write();
  hPVHFremDaupullsYVsN->Write();
  hPVHFremDauresZVsN->Write();
  hPVHFremDaumeanZVsN->Write();
  hPVHFremDaupullsZVsN->Write();  
  hPVESDHFchiS->Write();
  hPVESDHFprob->Write();
  hPVESDHFresXVsN->Write();
  hPVESDHFmeanXVsN->Write();
  hPVESDHFpullsXVsN->Write();
  hPVESDHFresYVsN->Write();
  hPVESDHFmeanYVsN->Write();
  hPVESDHFpullsYVsN->Write();
  hPVESDHFresZVsN->Write();
  hPVESDHFmeanZVsN->Write();
  hPVESDHFpullsZVsN->Write();
  hPVESDHFremDauchiS->Write();
  hPVESDHFremDauprob->Write();
  hPVESDHFremDauresXVsN->Write();
  hPVESDHFremDaumeanXVsN->Write();
  hPVESDHFremDaupullsXVsN->Write();
  hPVESDHFremDauresYVsN->Write();
  hPVESDHFremDaumeanYVsN->Write();
  hPVESDHFremDaupullsYVsN->Write();
  hPVESDHFremDauresZVsN->Write();
  hPVESDHFremDaumeanZVsN->Write();
  hPVESDHFremDaupullsZVsN->Write();
  cPVHFchiS->Write();
  cPVHFprob->Write();
  cPVHFresX->Write();
  cPVHFmeanX->Write();
  cPVHFpullsX->Write();
  cPVHFresY->Write();
  cPVHFmeanY->Write();
  cPVHFpullsY->Write();
  cPVHFresZ->Write();
  cPVHFmeanZ->Write();
  cPVHFpullsZ->Write();
  outfilePVHF.Close();

  delete cPVHFresX;
  delete cPVHFresY;
  delete cPVHFresZ;
  delete cPVHFmeanX;
  delete cPVHFmeanY;
  delete cPVHFmeanZ;
  delete cPVHFpullsX;
  delete cPVHFpullsY;
  delete cPVHFpullsZ;
  delete cPVHFprob;
  delete cPVHFchiS;
  delete cPVHFresXZ;
  delete cPVHFmeanXZ;
  
}

void ResetAxes(THnSparseF* sparse) {
  Int_t nAxes = sparse->GetNdimensions();
  for(Int_t iAxis=0; iAxis<nAxes; iAxis++) {
    sparse->GetAxis(iAxis)->SetRange(-1,-1);
  }
}

TH1F* GetHistoVsNtracks(THnSparseF* sparse, Int_t axis, TString outfilename, Int_t quantity) {
  Int_t par=2;
  if(quantity==kMean)
    par=1;
  
  TAxis* Naxis = (TAxis*)sparse->GetAxis(3);
  
  while((Nmax-Nmin)%nbins!=0) {
    cout << "the number of bins in the ThnSparse cannot be divided by "<< nbins << endl;
    nbins--;
    cout << "decreasing nbins to " << nbins << endl;
  } 
  
  Int_t step = (Nmax-Nmin)/nbins;

  TH1F* hVsN = new TH1F("hVsN","",nbins,Nmin,Nmax);
  TFile outfile(outfilename.Data(),"UPDATE");
  TCanvas *c = new TCanvas(); 
  for(Int_t iBin=0; iBin<nbins; iBin++) {
    ResetAxes(sparse);
    Int_t binmin = Naxis->FindBin(step*(iBin+1)*1.001);    
    Int_t binmax = Naxis->FindBin(step*(iBin+2)*0.999);
    Naxis->SetRange(binmin,binmax);
    TH1F* h = (TH1F*)sparse->Projection(axis);
    c->Clear();
    c->Update();
    if(h->GetEntries()!=0) {
      h->Fit("gaus");
      outfile.cd();
      h->Write();
      TF1* f=(TF1*)h->GetListOfFunctions()->FindObject("gaus");
      cout << f->GetParameter(par) << endl;
      hVsN->SetBinContent(iBin+1,f->GetParameter(par));
      hVsN->SetBinError(iBin+1,f->GetParError(par));
    }
    else {
      hVsN->SetBinContent(iBin+1,0);
    }
  }

  delete c;

  outfile.Close();
  hVsN->SetDirectory(0);
  TString ytitle = sparse->GetAxis(axis)->GetTitle();
  if(quantity==kMean) {
    ytitle.ReplaceAll("#sigma(","<");
    ytitle.ReplaceAll("C})","C}>");
    ytitle.ReplaceAll("reco})","reco}>");
  }
  hVsN->GetXaxis()->SetTitle(sparse->GetAxis(3)->GetTitle());
  hVsN->GetYaxis()->SetTitle(ytitle.Data());
  hVsN->SetStats(0);
  
  return hVsN;

}
