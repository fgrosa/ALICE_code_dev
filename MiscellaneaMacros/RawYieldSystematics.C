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

//*********************************************//
//                                             //
//    Main Function: RawYieldSystematics       //
//                                             //
//*********************************************//

//________________________________________________________________________________________________________________
//global variables
const TString cent = "3050";
const TString axesfile = "axes.txt";
const TString cutsfile = Form("Cent%s/cutset_topocut.txt",cent.Data());

const Int_t nDataFiles = 3;
const TString datafilename[nDataFiles] =
  {Form("$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/HR_bunch4/AnalysisResults_%s.root",cent.Data())
  ,Form("$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/HR_bunch1-3/AnalysisResults_%s.root",cent.Data())
  //,Form("$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/HR_bunch1-3/AnalysisResults_%s_part2.root",cent.Data())
  //,Form("$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/HR_bunch1-3/AnalysisResults_%s_part3.root",cent.Data())
  ,Form("$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/LR/AnalysisResults_%s.root",cent.Data())};

//const TString datafilename[nDataFiles] = {Form("$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/HR_bunch4/AnalysisResults_%s.root",cent.Data())};
const TString datadirname = "PWG3_D2H_InvMassDplus";
const TString datalistname = Form("coutputDplus%s_kINT7%s",cent.Data(),cent.Data());
const TString sparsename = "hMassPtCutVarsAll";

const TString reffilename = Form("Cent%s/RawYields_%s_topocut.root",cent.Data(),cent.Data());
const TString reffileMCname = Form("Cent%s/RawYieldsMC_%s.root",cent.Data(),cent.Data());

enum {kFree,kFixed};

const Int_t nMins = 4;
const Int_t nMaxs = 4;
const Int_t nReb = 5;
const Int_t nBkgFunc = 3;
const Double_t mins[nMins] = {1.68,1.69,1.70,1.71};
const Double_t maxs[nMaxs] = {2.01,2.02,2.03,2.04};
const Int_t rebin[nReb] = {2,4,5,6,8};
const Int_t sgnfcn = AliHFMassFitter::kGaus;
const Int_t bkgfcn[nBkgFunc] = {AliHFMassFitter::kExpo,AliHFMassFitter::kLin,AliHFMassFitter::kPol2};
const Double_t maxchisquare = 2.;
const Double_t nSigmaBinCount = 3.5;

//________________________________________________________________________________________________________________
//functions prototypes
Int_t RawYieldSystematics(TString outfilerawname=Form("Cent%s/RawYieldsSyst_%s_topocut.root",cent.Data(),cent.Data()));
Int_t ReadAxes(TString FileName, vector<string> &axesanmes, vector<int> &axesno);
Int_t ReadSet(TString FileName, vector<string> &varnames, vector<double> &cutset);
THnSparseF* GetSparse(TString filename, TString dirname, TString listname, TString sparsename);
Int_t LoadRefFiles(TString reffilename, TString reffileMCname, TH1F *&hRawYieldRef, TH1F *&hSigmaRef, TH1F *&hMeanRef, TH1F *&hSigmaMC);
void ApplyCuts(THnSparseF* sparse, vector<int> axesno, vector<string> axesnames, vector<double> cutset, Int_t iPt);
void ResetAxes(THnSparseF* sparse);
void GetPtLims(Double_t ptlims[], vector<double> cutset, const Int_t nCutVars);
void SetStyle();
void DivideCanvas(TCanvas* c, const Int_t nPtBins);

//________________________________________________________________________________________________________________
Int_t RawYieldSystematics(TString outfilerawname) {

  //get thnsparse from files
  THnSparseF** datasparse = new THnSparseF*[nDataFiles];
  
  for(Int_t iDataFile=0; iDataFile<nDataFiles; iDataFile++) {
    datasparse[iDataFile] = (THnSparseF*)GetSparse(datafilename[iDataFile],datadirname,datalistname,sparsename);
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
  TH1F* hMass = 0x0;
  TH1F* hMassToFit = 0x0;
  AliHFMassFitter* fitter = 0x0;
  TH1F** hRawYield = new TH1F*[nPtBins];
  TH1F** hSigma = new TH1F*[nPtBins];
  TH1F** hMean = new TH1F*[nPtBins];
  TH1F** hChiSquare = new TH1F*[nPtBins];
  TH1F** hRawYieldVsTrial = new TH1F*[nPtBins];
  TH1F** hSigmaVsTrial = new TH1F*[nPtBins];
  TH1F** hMeanVsTrial = new TH1F*[nPtBins];
  TH1F** hChiSquareVsTrial = new TH1F*[nPtBins];
  
  TLine** lRawRef = new TLine*[nPtBins];
  TLine** lSigmaRef = new TLine*[nPtBins];
  TLine** lMeanRef = new TLine*[nPtBins];
  TLine** lRawRefVsTrial = new TLine*[nPtBins];
  TLine** lSigmaRefVsTrial = new TLine*[nPtBins];
  TLine** lMeanRefVsTrial = new TLine*[nPtBins];
  
  TCanvas* cRaw = new TCanvas("cRaw","",10,10,1920,1080);
  DivideCanvas(cRaw, nPtBins);
  TCanvas* cSigma = new TCanvas("cSigma","",10,10,1920,1080);
  DivideCanvas(cSigma, nPtBins);
  TCanvas* cMean = new TCanvas("cMean","",10,10,1920,1080);
  DivideCanvas(cMean, nPtBins);
  TCanvas* cChi = new TCanvas("cChi","",10,10,1920,1080);
  DivideCanvas(cChi, nPtBins);
  TCanvas* cRawVsTrial = new TCanvas("cRawVsTrial","",10,10,1920,1080);
  DivideCanvas(cRawVsTrial, nPtBins);
  TCanvas* cSigmaVsTrial = new TCanvas("cSigmaVsTrial","",10,10,1920,1080);
  DivideCanvas(cSigmaVsTrial, nPtBins);
  TCanvas* cMeanVsTrial = new TCanvas("cMeanVsTrial","",10,10,1920,1080);
  DivideCanvas(cMeanVsTrial, nPtBins);
  TCanvas* cChiVsTrial = new TCanvas("cChiVsTrial","",10,10,1920,1080);
  DivideCanvas(cChiVsTrial, nPtBins);
  
  Double_t massD = TDatabasePDG::Instance()->GetParticle(411)->Mass();

  const Int_t nbins=50;
  const Int_t ntrials = nMins*nMaxs*2*2*nBkgFunc*nReb;
  TH1F* hRawYieldRef=0x0;
  TH1F* hSigmaRef=0x0;
  TH1F* hMeanRef=0x0;
  TH1F* hSigmaMC=0x0;
  Int_t loadref = LoadRefFiles(reffilename,reffileMCname,hRawYieldRef,hSigmaRef,hMeanRef,hSigmaMC);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    Double_t rawmin = 1.;
    Double_t rawmax = -1.;
    if(loadref!=1 && loadref!=2 && loadref!=4) {
      rawmin = hRawYieldRef->GetBinContent(iPt+1)*(1-0.5);
      rawmax = hRawYieldRef->GetBinContent(iPt+1)*(1+0.5);
      lRawRefVsTrial[iPt] = new TLine(-0.5,hRawYieldRef->GetBinContent(iPt+1),ntrials-0.5,hRawYieldRef->GetBinContent(iPt+1));
      lRawRefVsTrial[iPt]->SetLineColor(kRed);
      lRawRefVsTrial[iPt]->SetLineWidth(2);
    }
    if(PtLims[iPt]>=36) rawmin=0.;
    Double_t sigmamin = 1.;
    Double_t sigmamax = -1.;
    if(loadref!=1 && loadref!=2 && loadref!=5) {
      sigmamin = hSigmaRef->GetBinContent(iPt+1)*(1-0.5);
      sigmamax = hSigmaRef->GetBinContent(iPt+1)*(1+0.5);
      lSigmaRefVsTrial[iPt] = new TLine(-0.5,hSigmaRef->GetBinContent(iPt+1),ntrials-0.5,hSigmaRef->GetBinContent(iPt+1));
      lSigmaRefVsTrial[iPt]->SetLineColor(kRed);
      lSigmaRefVsTrial[iPt]->SetLineWidth(2);
    }
    Double_t meanmin = massD*(1-0.005);
    Double_t meanmax = massD*(1+0.005);
    if(loadref!=1 && loadref!=2 && loadref!=6) {
      lMeanRefVsTrial[iPt] = new TLine(-0.5,hMeanRef->GetBinContent(iPt+1),ntrials-0.5,hMeanRef->GetBinContent(iPt+1));
      lMeanRefVsTrial[iPt]->SetLineColor(kRed);
      lMeanRefVsTrial[iPt]->SetLineWidth(2);
    }
    
    hRawYield[iPt] = new TH1F(Form("hRawYield_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",nbins,rawmin,rawmax);
    hSigma[iPt] = new TH1F(Form("hMean_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",nbins,sigmamin,sigmamax);
    hMean[iPt] = new TH1F(Form("hSigma_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",nbins,meanmin,meanmax);
    hChiSquare[iPt] = new TH1F(Form("hChiSquare_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",nbins,0.,maxchisquare);
    hRawYield[iPt]->GetXaxis()->SetTitle("raw yield");
    hRawYield[iPt]->GetYaxis()->SetTitle("Entries");
    hRawYield[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hRawYield[iPt]->SetFillStyle(3004);
    hRawYield[iPt]->SetLineWidth(2);
    hRawYield[iPt]->SetLineColor(kBlue+1);
    hRawYield[iPt]->SetFillColor(kBlue+1);
    hSigma[iPt]->GetXaxis()->SetTitle("width (GeV/c^{2})");
    hSigma[iPt]->GetYaxis()->SetTitle("Entries");
    hSigma[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hSigma[iPt]->SetFillStyle(3004);
    hSigma[iPt]->SetLineWidth(2);
    hSigma[iPt]->SetLineColor(kBlue+1);
    hSigma[iPt]->SetFillColor(kBlue+1);
    hMean[iPt]->GetXaxis()->SetTitle("mean (GeV/c^{2})");
    hMean[iPt]->GetYaxis()->SetTitle("Entries");
    hMean[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hMean[iPt]->SetFillStyle(3004);
    hMean[iPt]->SetLineWidth(2);
    hMean[iPt]->SetLineColor(kBlue+1);
    hMean[iPt]->SetFillColor(kBlue+1);
    hChiSquare[iPt]->GetXaxis()->SetTitle("#chi^{2}");
    hChiSquare[iPt]->GetYaxis()->SetTitle("Entries");
    hChiSquare[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hChiSquare[iPt]->SetFillStyle(3004);
    hChiSquare[iPt]->SetLineWidth(2);
    hChiSquare[iPt]->SetLineColor(kBlue+1);
    hChiSquare[iPt]->SetFillColor(kBlue+1);
   
    hRawYieldVsTrial[iPt] = new TH1F(Form("hRawYieldVsTrial_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",ntrials,-0.5,ntrials-0.5);
    hSigmaVsTrial[iPt] = new TH1F(Form("hMeanVsTrial_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",ntrials,-0.5,ntrials-0.5);
    hMeanVsTrial[iPt] = new TH1F(Form("hSigmaVsTrial_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",ntrials,-0.5,ntrials-0.5);
    hChiSquareVsTrial[iPt] = new TH1F(Form("hChiSquareVsTrial_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",ntrials,-0.5,ntrials-0.5);
    hRawYieldVsTrial[iPt]->GetYaxis()->SetRangeUser(rawmin,rawmax);
    hRawYieldVsTrial[iPt]->GetYaxis()->SetTitle("raw yield");
    hRawYieldVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hRawYieldVsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hRawYieldVsTrial[iPt]->SetMarkerStyle(20);
    hRawYieldVsTrial[iPt]->SetLineWidth(2);
    hRawYieldVsTrial[iPt]->SetMarkerSize(0.5);
    hRawYieldVsTrial[iPt]->SetLineColor(kBlue+1);
    hRawYieldVsTrial[iPt]->SetMarkerColor(kBlue+1);
    hSigmaVsTrial[iPt]->GetYaxis()->SetRangeUser(sigmamin,sigmamax);
    hSigmaVsTrial[iPt]->GetYaxis()->SetTitle("width (GeV/c^{2})");
    hSigmaVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hSigmaVsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hSigmaVsTrial[iPt]->SetMarkerStyle(20);
    hSigmaVsTrial[iPt]->SetLineWidth(2);
    hSigmaVsTrial[iPt]->SetMarkerSize(0.5);
    hSigmaVsTrial[iPt]->SetLineColor(kBlue+1);
    hSigmaVsTrial[iPt]->SetMarkerColor(kBlue+1);
    hMeanVsTrial[iPt]->GetYaxis()->SetRangeUser(meanmin,meanmax);
    hMeanVsTrial[iPt]->GetYaxis()->SetTitle("mean (GeV/c^{2})");
    hMeanVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hMeanVsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hMeanVsTrial[iPt]->SetMarkerStyle(20);
    hMeanVsTrial[iPt]->SetLineWidth(2);
    hMeanVsTrial[iPt]->SetMarkerSize(0.5);
    hMeanVsTrial[iPt]->SetLineColor(kBlue+1);
    hMeanVsTrial[iPt]->SetMarkerColor(kBlue+1);
    hChiSquareVsTrial[iPt]->GetYaxis()->SetRangeUser(0,maxchisquare);
    hChiSquareVsTrial[iPt]->GetYaxis()->SetTitle("#chi^{2}");
    hChiSquareVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hChiSquareVsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hChiSquareVsTrial[iPt]->SetMarkerStyle(20);
    hChiSquareVsTrial[iPt]->SetLineWidth(2);
    hChiSquareVsTrial[iPt]->SetMarkerSize(0.5);
    hChiSquareVsTrial[iPt]->SetLineColor(kBlue+1);
    hChiSquareVsTrial[iPt]->SetMarkerColor(kBlue+1);
    
    for(Int_t iDataFile=0; iDataFile<nDataFiles; iDataFile++) {
      ApplyCuts(datasparse[iDataFile],axesno,axesnames,cutset,iPt);
      hMassPart[iDataFile] = (TH1F*)datasparse[iDataFile]->Projection(0);
    }
    hMass=(TH1F*)hMassPart[0]->Clone();
    if(nDataFiles>1) {
      for(Int_t iBin=0; iBin<hMass->GetNbinsX(); iBin++) {
        Double_t nentries=0;
        for(Int_t iDataFile=0; iDataFile<nDataFiles; iDataFile++) {
          nentries += hMassPart[iDataFile]->GetBinContent(iBin+1);
        }
        hMass->SetBinContent(iBin+1,nentries);
      }
    }
    
    Int_t iTrial=0;
    
    for(Int_t iSigma=kFree; iSigma<=kFixed; iSigma++) {
      for(Int_t iMean=kFree; iMean<=kFixed; iMean++) {
        for(Int_t iReb=0; iReb<nReb; iReb++) {
          for(Int_t iMin=0; iMin<nMins; iMin++) {
            for(Int_t iMax=0; iMax<nMaxs; iMax++) {
              for(Int_t iBkgFunc=0; iBkgFunc<nBkgFunc; iBkgFunc++) {
        
                hMassToFit=(TH1F*)hMass->Clone();
                hMassToFit->Rebin(rebin[iReb]);
                fitter = new AliHFMassFitter(hMassToFit,mins[iMin],maxs[iMax],1,bkgfcn[iBkgFunc],sgnfcn);
                if(iMean==kFixed) {fitter->SetFixGaussianMean(massD);}
                else {fitter->SetInitialGaussianMean(massD);}
                
                if(loadref!=2 && loadref!=3 && loadref!=7) {
                  if(iSigma==kFixed) {fitter->SetFixGaussianSigma(hSigmaMC->GetBinContent(iPt+1));}
                  else {fitter->SetInitialGaussianSigma(hSigmaMC->GetBinContent(iPt+1));}
                }
                else {
                  if(iSigma==kFixed) {continue;}
                  else {fitter->SetInitialGaussianSigma(0.012);}
                }
                fitter->MassFitter(kFALSE);
    
                Double_t chi=fitter->GetReducedChiSquare();
                if(chi<maxchisquare) {
                  Double_t rawyield=fitter->GetRawYield();
                  Double_t sigma=fitter->GetSigma();
                  Double_t mean=fitter->GetMean();
                  Double_t rawyielderr=fitter->GetRawYieldError();
                  Double_t sigmaerr=fitter->GetSigmaUncertainty();
                  Double_t meanerr=fitter->GetMeanUncertainty();
                  hRawYield[iPt]->Fill(rawyield);
                  hSigma[iPt]->Fill(sigma);
                  hMean[iPt]->Fill(mean);
                  hChiSquare[iPt]->Fill(chi);
                  hRawYieldVsTrial[iPt]->SetBinContent(iTrial+1,rawyield);
                  hSigmaVsTrial[iPt]->SetBinContent(iTrial+1,sigma);
                  hMeanVsTrial[iPt]->SetBinContent(iTrial+1,mean);
                  hChiSquareVsTrial[iPt]->SetBinContent(iTrial+1,chi);
                  hRawYieldVsTrial[iPt]->SetBinError(iTrial+1,rawyielderr);
                  hSigmaVsTrial[iPt]->SetBinError(iTrial+1,sigmaerr);
                  hMeanVsTrial[iPt]->SetBinError(iTrial+1,meanerr);
                }
                
                iTrial++;
                
                delete hMassToFit;
                hMassToFit=0x0;
              }
            }
          }
        }
      }
    }
    
    if(nPtBins>1) {cRaw->cd(iPt+1);}
    else {cRaw->cd();}
    hRawYield[iPt]->Draw();
    if(loadref!=1 && loadref!=2 && loadref!=4) {
      lRawRef[iPt] = new TLine(hRawYieldRef->GetBinContent(iPt+1),0,hRawYieldRef->GetBinContent(iPt+1),hRawYield[iPt]->GetMaximum());
      lRawRef[iPt]->SetLineColor(kRed);
      lRawRef[iPt]->SetLineWidth(2);
      lRawRef[iPt]->Draw("same");
    }
    if(nPtBins>1) {cSigma->cd(iPt+1);}
    else {cSigma->cd();}
    hSigma[iPt]->Draw();
    if(loadref!=1 && loadref!=2 && loadref!=5) {
      lSigmaRef[iPt] = new TLine(hSigmaRef->GetBinContent(iPt+1),0,hSigmaRef->GetBinContent(iPt+1),hSigma[iPt]->GetMaximum());
      lSigmaRef[iPt]->SetLineColor(kRed);
      lSigmaRef[iPt]->SetLineWidth(2);
      lSigmaRef[iPt]->Draw("same");
    }
    if(nPtBins>1) {cMean->cd(iPt+1);}
    else {cMean->cd();}
    hMean[iPt]->Draw();
    if(loadref!=1 && loadref!=2 && loadref!=6) {
      lMeanRef[iPt] = new TLine(hMeanRef->GetBinContent(iPt+1),0,hMeanRef->GetBinContent(iPt+1),hMean[iPt]->GetMaximum());
      lMeanRef[iPt]->SetLineColor(kRed);
      lMeanRef[iPt]->SetLineWidth(2);
      lMeanRef[iPt]->Draw("same");
    }
    if(nPtBins>1) {cChi->cd(iPt+1);}
    else {cChi->cd();}
    hChiSquare[iPt]->Draw();
    
    if(nPtBins>1) {cRawVsTrial->cd(iPt+1);}
    else {cRawVsTrial->cd();}
    hRawYieldVsTrial[iPt]->Draw();
    if(lRawRefVsTrial[iPt]) {lRawRefVsTrial[iPt]->Draw("same");}
    if(nPtBins>1) {cSigmaVsTrial->cd(iPt+1);}
    else {cSigmaVsTrial->cd();}
    hSigmaVsTrial[iPt]->Draw();
    if(lSigmaRefVsTrial[iPt]) {lSigmaRefVsTrial[iPt]->Draw("same");}
    if(nPtBins>1) {cMeanVsTrial->cd(iPt+1);}
    else {cMeanVsTrial->cd();}
    hMeanVsTrial[iPt]->Draw();
    if(lMeanRefVsTrial[iPt]) {lMeanRefVsTrial[iPt]->Draw("same");}
    if(nPtBins>1) {cChiVsTrial->cd(iPt+1);}
    else {cChiVsTrial->cd();}
    hChiSquareVsTrial[iPt]->Draw();
    
    for(Int_t iDataFile=0; iDataFile<nDataFiles; iDataFile++) {ResetAxes(datasparse[iDataFile]);}
  }
  
  //output files
  TFile outfileraw(outfilerawname.Data(),"RECREATE");
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    hRawYield[iPt]->Write();
    hSigma[iPt]->Write();
    hMean[iPt]->Write();
    hChiSquare[iPt]->Write();
  }
  outfileraw.Close();
  cout << "\n" <<outfilerawname << " saved." <<endl;
  outfilerawname.ReplaceAll(".root","_RawYield.pdf");
  cRaw->SaveAs(outfilerawname.Data());
  outfilerawname.ReplaceAll("_RawYield.pdf","_Sigma.pdf");
  cSigma->SaveAs(outfilerawname.Data());
  outfilerawname.ReplaceAll("_Sigma.pdf","_Mean.pdf");
  cMean->SaveAs(outfilerawname.Data());
  outfilerawname.ReplaceAll("_Mean.pdf","_ChiSquare.pdf");
  cChi->SaveAs(outfilerawname.Data());
  
  delete fitter;
  delete hMass;
  for(Int_t iDataFile=0; iDataFile<nDataFiles; iDataFile++) {if(hMassPart[iDataFile]) delete hMassPart[iDataFile];}
  delete[] hMassPart;

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
Int_t LoadRefFiles(TString reffilename, TString reffileMCname, TH1F *&hRawYieldRef, TH1F *&hSigmaRef, TH1F *&hMeanRef, TH1F *&hSigmaMC) {

  TFile* reffile = TFile::Open(reffilename.Data(),"READ");
  if(reffile) {
    hRawYieldRef=(TH1F*)reffile->Get("hRawYields");
    hSigmaRef=(TH1F*)reffile->Get("hRawYieldsSigma");
    hMeanRef=(TH1F*)reffile->Get("hRawYieldsMean");
    if(hRawYieldRef) {hRawYieldRef->SetDirectory(0);}
    if(hSigmaRef) {hSigmaRef->SetDirectory(0);}
    if(hMeanRef) {hMeanRef->SetDirectory(0);}
    reffile->Close();
  }
  
  TFile* reffileMC = TFile::Open(reffileMCname.Data(),"READ");
  if(reffileMC) {
    hSigmaMC=(TH1F*)reffileMC->Get("hRawYieldsSigma");
    if(hSigmaMC) {hSigmaMC->SetDirectory(0);}
    reffileMC->Close();
  }
  
  if(!reffile && !reffileMC) return 1;
  if(!reffile) return 2;
  if(!reffileMC) return 3;
  if(!hRawYieldRef) return 4;
  if(!hSigmaRef) return 5;
  if(!hMeanRef) return 6;
  if(!hSigmaMC) return 7;
  return 0;
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
void SetStyle() {
  cout << "Setting drawing style!" << endl;
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetOptStat(1111);
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
