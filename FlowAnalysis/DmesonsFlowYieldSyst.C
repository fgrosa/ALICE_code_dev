#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TString.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TASImage.h>
#include <TPad.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <TLatex.h>
#include <TLine.h>

#include "AliHFMassFitter.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliVertexingHFUtils.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliEventPlaneResolutionHandler.h"
#include "AliVertexingHFUtils.h"

#endif

//methods for the extraction yield systematics for the v2 computed with the Event Plane method
//Author: Fabrizio Grosa, INFN Turin grosa@to.infn.it

//*************************************************//
//                                                 //
//      Main Function: DmesonsFlowYieldSyst()      //
//                                                 //
//*************************************************//

//_________________________________________________________________
//GLOBAL VARIABLES TO BE SET
//input file
const TString filename="$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/AnalysisResults_3050_flow.root";
const TString suffix="TPC";//"_3050_CentralCuts_TPC";
const TString partname="Dplus";
const Int_t minCent=30;
const Int_t maxCent=50;

//EP resolution
//kTPCFullEta, kTPCPosEta,kVZERO,kVZEROA,kVZEROC
const Int_t evPlane=AliEventPlaneResolutionHandler::kTPCFullEta;
//resolution flag fromAliEventPlaneResolutionHandler:
//kTwoRandSub,kTwoChargeSub,kTwoEtaSub,kThreeSub,kThreeSubTPCGap
const Bool_t useAliHandlerForRes=kFALSE;
Int_t evPlaneRes=AliEventPlaneResolutionHandler::kTwoEtaSub;
const Bool_t useNcollWeight=kFALSE;

// pt and phi binning
const Int_t nptbinsnew=8;
const Double_t ptbinsnew[nptbinsnew+1]={2.,3.,4.,5.,6.,8.,10.,12.,16.};
const Int_t nphibins=4;
const Double_t phibinslim[nphibins+1]={0,TMath::Pi()/4,TMath::Pi()/2,3*TMath::Pi()/4,TMath::Pi()};

// mass fit configuration
const Int_t nReb=4;
const Int_t rebin[nptbinsnew]={2,4,5,6};
const Int_t nBkgFcn=3;
const Int_t typeb[nBkgFcn]={AliHFMassFitter::kExpo,AliHFMassFitter::kLin,AliHFMassFitter::kPol2};
const Int_t nMins=4;
const Double_t minMassForFit[nMins]={1.67,1.68,1.69,1.7};
const Int_t nMaxs=4;
const Double_t maxMassForFit[nMaxs]={2.05,2.04,2.03,2.02};
const Double_t nSigmaForCounting=3.5;
const Bool_t fixAlsoMass=kFALSE;

//not to be set
Int_t minPtBin[nptbinsnew]={-1,-1,-1,-1};
Int_t maxPtBin[nptbinsnew]={-1,-1,-1,-1};
const Double_t effInOverEffOut=1.03;
Double_t massD;

const Int_t colors[] = {kRed+1,kBlack,kBlue+1,kGreen+2,kOrange+7,kBlue-7};
const Int_t markers[] = {kFullSquare,kFullCircle,kFullTriangleUp,kFullDiamond,kOpenSquare,kOpenCircle,kOpenTriangleUp,kOpenDiamond};

//_________________________________________________________________
//METHODS PROTOTYPES
void DmesonsFlowYieldSyst(Bool_t inoutanis=kTRUE);
TList *LoadMassHistos(TList *inputlist,Bool_t inoutanis);
TList* LoadResolutionHistos(TList *inputlist);
Int_t FindPtBin(Int_t nbins, Double_t* array,Double_t value);
void FillSignalGraph(TList *histlist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2, Bool_t inoutanis, Int_t bkgfunc, Int_t minfit, Int_t maxfit, Int_t rebin);
TGraphAsymmErrors* Computev2(TGraphAsymmErrors **gSignal, Double_t resol, Float_t *averagePt, Bool_t inoutanis, TGraphAsymmErrors *gRelSystEff);
Double_t GetEventPlaneResolution(Double_t& error, TH1F* hsubev1, TH1F* hsubev2, TH1F* hsubev3);
Bool_t DefinePtBins(AliRDHFCuts *cutobj);
Int_t LoadRefGraphs(TString reffilename, TH1F **hRawYieldRef, TH1F **hRawYieldfsRef, TH1F **hRawYieldBC1Ref, TH1F **hRawYieldBC2Ref, TGraphAsymmErrors *&gv2Ref, TGraphAsymmErrors *&gv2fsRef, TGraphAsymmErrors *&gv2BC1Ref, TGraphAsymmErrors *&gv2BC2Ref, Bool_t inoutanis);
void DivideCanvas(TCanvas* c, Int_t nPtBins);
void SetStyle(Int_t optfit=0);

//_________________________________________________________________
//METHODS IMPLEMENTATION
void DmesonsFlowYieldSyst(Bool_t inoutanis){
  
  TString dirname=Form("PWGHF_D2H_HFv2_%s%s",partname.Data(),suffix.Data());
  TString listname=Form("coutputv2%s%s",partname.Data(),suffix.Data());
  
  AliRDHFCuts *cutsobj=0x0;
  //Load input data from AliAnalysisTaskSEHFv2
  TFile *f = TFile::Open(filename.Data());
  if(!f){
    printf("file %s not found, please check file name\n",filename.Data());return;
  }
  TDirectoryFile* dir=(TDirectoryFile*)f->Get(dirname.Data());
  if(!dir){
    printf("Directory %s not found, please check dir name\n",dirname.Data());return;
  }
  if(partname.Contains("Dzero")) {
    cutsobj=((AliRDHFCutsD0toKpi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
    massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  }
  if(partname.Contains("Dplus")){
    cutsobj=((AliRDHFCutsDplustoKpipi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
    massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  }
  if(partname.Contains("Dstar")) {
    cutsobj=((AliRDHFCutsDStartoKpipi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
    massD=(TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass());
  }
  if(partname.Contains("Ds")) {
    cutsobj=((AliRDHFCutsDstoKKpi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
    massD=(TDatabasePDG::Instance()->GetParticle(431)->Mass());
  }
  
  TList *list =(TList*)dir->Get(listname.Data());
  if(!list){
    printf("list %s not found in file, please check list name\n",listname.Data());return;
  }
  if(!cutsobj){
    printf("cut object not found in file, please check keylist number\n");return;
  }
  //Define new pt bins
  if(!DefinePtBins(cutsobj)){
    printf("cut not define pt bins\n");return;
  }
  
  //Load mass histograms corresponding to the required centrality, pt range and phi binning
  printf("Load mass histos \n");
  TList *histlist=LoadMassHistos(list,inoutanis);
  TString aniss="";
  if(inoutanis)aniss+="anis";
  histlist->SaveAs(Form("v2Histograms_%d_%d_%s_%s.root",minCent,maxCent,aniss.Data(),suffix.Data()),"RECREATE");
  
  Int_t nphi=nphibins;
  if(inoutanis)nphi=2;
  
  printf("average pt for pt bin \n");
  //average pt for pt bin
  AliVertexingHFUtils *utils=new AliVertexingHFUtils();
  Int_t minCentTimesTen=minCent*10;
  Int_t maxCentTimesTen=maxCent*10;
  TH2F* hmasspt=(TH2F*)list->FindObject(Form("hMPtCandcentr%d_%d",minCentTimesTen,minCentTimesTen+25));
  for(Int_t icent=minCentTimesTen+25;icent<maxCentTimesTen;icent=icent+25)hmasspt->Add((TH2F*)list->FindObject(Form("hMPtCandcentr%d_%d",icent,icent+25)));
  Float_t averagePt[nptbinsnew];
  Float_t errorPt[nptbinsnew];
  for(Int_t ipt=0;ipt<nptbinsnew;ipt++){
    Int_t binMin=hmasspt->FindBin(ptbinsnew[ipt]);
    Int_t binMax=hmasspt->FindBin(ptbinsnew[ipt+1]-0.001);
    if(TMath::Abs(hmasspt->GetXaxis()->GetBinLowEdge(binMin)-ptbinsnew[ipt])>0.001 ||
       TMath::Abs(hmasspt->GetXaxis()->GetBinUpEdge(binMax)-ptbinsnew[ipt+1])>0.001){
      printf("Error in pt bin limits for projection!\n");
      return;
    }
    TH1F *histtofit = (TH1F*)hmasspt->ProjectionY("_py",binMin,binMax);
    Int_t nMassBins=histtofit->GetNbinsX();
    Double_t hmin=histtofit->GetBinLowEdge(2); // need wide range for <pt>
    Double_t hmax=histtofit->GetBinLowEdge(nMassBins-2); // need wide range for <pt>
    AliHFMassFitter fitter(histtofit,hmin,hmax,1);
    fitter.MassFitter(kFALSE);
    Double_t massFromFit=fitter.GetMean();
    Double_t sigmaFromFit=fitter.GetSigma();
    TF1* funcB2=fitter.GetBackgroundRecalcFunc();
    utils->AveragePt(averagePt[ipt],errorPt[ipt],ptbinsnew[ipt],ptbinsnew[ipt+1],hmasspt,massFromFit,sigmaFromFit,funcB2,2.5,4.5,0.,3.,1);
  }
  printf("Average pt\n");
  for(Int_t ipt=0;ipt<nptbinsnew;ipt++) printf("%f +- %f\n",averagePt[ipt],errorPt[ipt]);
  
  printf("Fill TGraphs for signal \n");
  //Fill TGraphs for signal
  TGraphAsymmErrors *gSignal[nptbinsnew];
  TGraphAsymmErrors *gSignalfs[nptbinsnew];
  TGraphAsymmErrors *gSignalBC1[nptbinsnew];
  TGraphAsymmErrors *gSignalBC2[nptbinsnew];
  for(Int_t i=0;i<nptbinsnew;i++){
    gSignal[i]=new TGraphAsymmErrors(nphi);
    gSignal[i]->SetName(Form("gasigpt%d",i));
    gSignal[i]->SetTitle(Form("Signal %.1f<#it{p}_{T}<%.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gSignal[i]->SetMarkerStyle(25);
    gSignalfs[i]=new TGraphAsymmErrors(nphi);
    gSignalfs[i]->SetName(Form("gasigfspt%d",i));
    gSignalfs[i]->SetTitle(Form("Signal (fixed sigma) %.1f<#it{p}_{T}<%.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gSignalfs[i]->SetMarkerStyle(21);
    gSignalBC1[i]=new TGraphAsymmErrors(nphi);
    gSignalBC1[i]->SetName(Form("gasigBC1pt%d",i));
    gSignalBC1[i]->SetTitle(Form("Signal (BC1) %.1f<#it{p}_{T}<%.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gSignalBC2[i]=new TGraphAsymmErrors(nphi);
    gSignalBC2[i]->SetName(Form("gasigBC2pt%d",i));
    gSignalBC2[i]->SetTitle(Form("Signal (BC2) %.1f<#it{p}_{T}<%.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
  }
  
  //EP resolution
  Double_t resol=-1.;
  Double_t errorres=-1.;
  
  if(useAliHandlerForRes) {
    AliEventPlaneResolutionHandler* epres=new AliEventPlaneResolutionHandler();
    epres->SetEventPlane(evPlane);
    epres->SetResolutionOption(evPlaneRes);
    if(useNcollWeight)
      epres->SetUseNcollWeights();
    resol=epres->GetEventPlaneResolution(minCent,maxCent);
    delete epres;
  }
  else {
    TList* resolhist=LoadResolutionHistos(list);
    TString suffixcentr=Form("centr%d_%d",minCent,maxCent);
    TH1F* hevplresos[3];
    TString namereso[3]={"Reso","Reso2","Reso3"};
    Int_t nSubRes=1;
    if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
       evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
    for(Int_t ires=0;ires<nSubRes;ires++){
      hevplresos[ires]=(TH1F*)resolhist->FindObject(Form("hEvPlane%s%s",namereso[ires].Data(),suffixcentr.Data()));
    }
    resol=GetEventPlaneResolution(errorres,hevplresos[0],hevplresos[1],hevplresos[2]);
  }
  
  printf("Event plane resolution %f\n",resol);
  printf("Compute v2 \n");
  //compute v2
  
  //load reference graphs
  TGraphAsymmErrors* gv2Ref=0x0;
  TGraphAsymmErrors* gv2fsRef=0x0;
  TGraphAsymmErrors* gv2BC1Ref=0x0;
  TGraphAsymmErrors* gv2BC2Ref=0x0;
  TH1F** hRawYieldRef = new TH1F*[nphi];
  TH1F** hRawYieldfsRef = new TH1F*[nphi];
  TH1F** hRawYieldBC1Ref = new TH1F*[nphi];
  TH1F** hRawYieldBC2Ref = new TH1F*[nphi];
  TString reffilename = Form("v2Output_%d_%d_%s_%s.root",minCent,maxCent,aniss.Data(),suffix.Data());
  Int_t loadref=LoadRefGraphs(reffilename, hRawYieldRef, hRawYieldfsRef, hRawYieldBC1Ref, hRawYieldBC2Ref, gv2Ref, gv2fsRef, gv2BC1Ref, gv2BC2Ref, inoutanis);
  
  TCanvas *cv2 =new TCanvas("cv2","v2 - systematic on yield extraction",1920,1080);
  DivideCanvas(cv2,nptbinsnew);
  TCanvas *cv2VsTrial =new TCanvas("cv2VsTrial","v2 vs. Trial- systematic on yield extraction",1920,1080);
  DivideCanvas(cv2VsTrial,nptbinsnew);

  TH1F** hv2 = new TH1F*[nptbinsnew];
  TH1F** hv2fs = new TH1F*[nptbinsnew];
  TH1F** hv2BC1 = new TH1F*[nptbinsnew];
  TH1F** hv2BC2 = new TH1F*[nptbinsnew];
  TH1F** hv2VsTrial = new TH1F*[nptbinsnew];
  TH1F** hv2fsVsTrial = new TH1F*[nptbinsnew];
  TH1F** hv2BC1VsTrial = new TH1F*[nptbinsnew];
  TH1F** hv2BC2VsTrial = new TH1F*[nptbinsnew];

  TCanvas** cRawYield=new TCanvas*[nphi];
  
  TH1F*** hRawYield=new TH1F**[nphi];
  TH1F*** hRawYieldfs=new TH1F**[nphi];
  TH1F*** hRawYieldBC1=new TH1F**[nphi];
  TH1F*** hRawYieldBC2=new TH1F**[nphi];
  for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
    hRawYield[iPhi]=new TH1F*[nptbinsnew];
    hRawYieldfs[iPhi]=new TH1F*[nptbinsnew];
    hRawYieldBC1[iPhi]=new TH1F*[nptbinsnew];
    hRawYieldBC2[iPhi]=new TH1F*[nptbinsnew];
    cRawYield[iPhi]=new TCanvas(Form("cRawYield_phi%d",iPhi),"Y - systematic on yield extraction",1920,1080);
    DivideCanvas(cRawYield[iPhi],nptbinsnew);
  }
  
  const Int_t nbins=50;
  
  const Int_t nTrials=nReb*nMins*nMaxs*nBkgFcn;
  
  for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
    hv2[iPt] = new TH1F(Form("hv2_%d",iPt),"",nbins,-0.3,0.7);
    hv2fs[iPt] = new TH1F(Form("hv2fs_%d",iPt),"",nbins,-0.3,0.7);
    hv2BC1[iPt] = new TH1F(Form("hv2BC1_%d",iPt),"",nbins,-0.3,0.7);
    hv2BC2[iPt] = new TH1F(Form("hv2BC2_%d",iPt),"",nbins,-0.3,0.7);
    
    hv2VsTrial[iPt] = new TH1F(Form("hv2VsTrial_%d",iPt),"",nTrials,-0.5,nTrials-0.5);
    hv2fsVsTrial[iPt] = new TH1F(Form("hv2fsVsTrial_%d",iPt),"",nTrials,-0.5,nTrials-0.5);
    hv2BC1VsTrial[iPt] = new TH1F(Form("hv2BC1VsTrial_%d",iPt),"",nTrials,-0.5,nTrials-0.5);
    hv2BC2VsTrial[iPt] = new TH1F(Form("hv2BC2VsTrial_%d",iPt),"",nTrials,-0.5,nTrials-0.5);
    hv2VsTrial[iPt]->SetStats(kFALSE);
    hv2fsVsTrial[iPt]->SetStats(kFALSE);
    hv2BC1VsTrial[iPt]->SetStats(kFALSE);
    hv2BC2VsTrial[iPt]->SetStats(kFALSE);
    
    hv2[iPt]->GetXaxis()->SetTitle("v_{2} {EP}");
    hv2[iPt]->GetYaxis()->SetTitle("Entries");
    hv2[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2fs[iPt]->GetXaxis()->SetTitle("v_{2} {EP}");
    hv2fs[iPt]->GetYaxis()->SetTitle("Entries");
    hv2fs[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2BC1[iPt]->GetXaxis()->SetTitle("v_{2} {EP}");
    hv2BC1[iPt]->GetYaxis()->SetTitle("Entries");
    hv2BC1[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2BC2[iPt]->GetXaxis()->SetTitle("v_{2} {EP}");
    hv2BC2[iPt]->GetYaxis()->SetTitle("Entries");
    hv2BC2[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2[iPt]->SetLineColor(colors[0]);
    hv2fs[iPt]->SetLineColor(colors[1]);
    hv2BC1[iPt]->SetLineColor(colors[2]);
    hv2BC2[iPt]->SetLineColor(colors[3]);
    hv2[iPt]->SetFillColor(colors[0]);
    hv2fs[iPt]->SetFillColor(colors[1]);
    hv2BC1[iPt]->SetFillColor(colors[2]);
    hv2BC2[iPt]->SetFillColor(colors[3]);
    hv2[iPt]->SetFillStyle(3004);
    hv2fs[iPt]->SetFillStyle(3004);
    hv2BC1[iPt]->SetFillStyle(3004);
    hv2BC2[iPt]->SetFillStyle(3004);
 
    hv2VsTrial[iPt]->GetYaxis()->SetTitle("v_{2} {EP}");
    hv2VsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hv2VsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2fsVsTrial[iPt]->GetYaxis()->SetTitle("v_{2} {EP}");
    hv2fsVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hv2fsVsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2BC1VsTrial[iPt]->GetYaxis()->SetTitle("v_{2} {EP}");
    hv2BC1VsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hv2BC1VsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2BC2VsTrial[iPt]->GetYaxis()->SetTitle("v_{2} {EP}");
    hv2BC2VsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hv2BC2VsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2VsTrial[iPt]->SetLineColor(colors[0]);
    hv2fsVsTrial[iPt]->SetLineColor(colors[1]);
    hv2BC1VsTrial[iPt]->SetLineColor(colors[2]);
    hv2BC2VsTrial[iPt]->SetLineColor(colors[3]);
    hv2VsTrial[iPt]->SetMarkerColor(colors[0]);
    hv2fsVsTrial[iPt]->SetMarkerColor(colors[1]);
    hv2BC1VsTrial[iPt]->SetMarkerColor(colors[2]);
    hv2BC2VsTrial[iPt]->SetMarkerColor(colors[3]);
    hv2VsTrial[iPt]->SetMarkerSize(0.5);
    hv2fsVsTrial[iPt]->SetMarkerSize(0.5);
    hv2BC1VsTrial[iPt]->SetMarkerSize(0.5);
    hv2BC2VsTrial[iPt]->SetMarkerSize(0.5);
    
    for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
      
      Int_t minraw=hRawYieldRef[iPhi]->GetBinContent(iPt+1)*(1-0.5);
      Int_t maxraw=hRawYieldRef[iPhi]->GetBinContent(iPt+1)*(1+0.5);
      
      hRawYield[iPhi][iPt] = new TH1F(Form("hRawYield_%d_%d",iPhi,iPt),"",nbins,minraw,maxraw);
      hRawYieldfs[iPhi][iPt] = new TH1F(Form("hRawYieldfs_%d_%d",iPhi,iPt),"",nbins,minraw,maxraw);
      hRawYieldBC1[iPhi][iPt] = new TH1F(Form("hRawYieldBC1_%d_%d",iPhi,iPt),"",nbins,minraw,maxraw);
      hRawYieldBC2[iPhi][iPt] = new TH1F(Form("hRawYieldBC2_%d_%d",iPhi,iPt),"",nbins,minraw,maxraw);
      hRawYield[iPhi][iPt]->GetXaxis()->SetTitle("raw yield");
      hRawYield[iPhi][iPt]->GetYaxis()->SetTitle("Entries");
      hRawYield[iPhi][iPt]->SetTitle(Form("#phi%d, %0.f < #it{p}_{T} < %0.f GeV/c",iPhi,ptbinsnew[iPt],ptbinsnew[iPt+1]));
      hRawYieldfs[iPhi][iPt]->GetXaxis()->SetTitle("raw yield");
      hRawYieldfs[iPhi][iPt]->GetYaxis()->SetTitle("Entries");
      hRawYieldfs[iPhi][iPt]->SetTitle(Form("#phi%d, %0.f < #it{p}_{T} < %0.f GeV/c",iPhi,ptbinsnew[iPt],ptbinsnew[iPt+1]));
      hRawYieldBC1[iPhi][iPt]->GetXaxis()->SetTitle("raw yield");
      hRawYieldBC1[iPhi][iPt]->GetYaxis()->SetTitle("Entries");
      hRawYieldBC1[iPhi][iPt]->SetTitle(Form("#phi%d, %0.f < #it{p}_{T} < %0.f GeV/c",iPhi,ptbinsnew[iPt],ptbinsnew[iPt+1]));
      hRawYieldBC2[iPhi][iPt]->GetXaxis()->SetTitle("raw yield");
      hRawYieldBC2[iPhi][iPt]->GetYaxis()->SetTitle("Entries");
      hRawYieldBC2[iPhi][iPt]->SetTitle(Form("#phi%d, %0.f < #it{p}_{T} < %0.f GeV/c",iPhi,ptbinsnew[iPt],ptbinsnew[iPt+1]));
      hRawYield[iPhi][iPt]->SetLineColor(colors[0]);
      hRawYield[iPhi][iPt]->SetFillColor(colors[0]);
      hRawYield[iPhi][iPt]->SetFillStyle(3004);
      hRawYieldfs[iPhi][iPt]->SetLineColor(colors[1]);
      hRawYieldfs[iPhi][iPt]->SetFillColor(colors[1]);
      hRawYieldfs[iPhi][iPt]->SetFillStyle(3004);
      hRawYieldBC1[iPhi][iPt]->SetLineColor(colors[2]);
      hRawYieldBC1[iPhi][iPt]->SetFillColor(colors[2]);
      hRawYieldBC1[iPhi][iPt]->SetFillStyle(3004);
      hRawYieldBC2[iPhi][iPt]->SetLineColor(colors[3]);
      hRawYieldBC2[iPhi][iPt]->SetFillColor(colors[3]);
      hRawYieldBC2[iPhi][iPt]->SetFillStyle(3004);
    }
  }
  
  Int_t iTrial=0;
  for(Int_t iReb=0; iReb<nReb; iReb++) {
    for(Int_t iMin=0; iMin<nMins; iMin++) {
      for(Int_t iMax=0; iMax<nMaxs; iMax++) {
        for(Int_t iBkgFcn=0; iBkgFcn<nBkgFcn; iBkgFcn++) {
      
          FillSignalGraph(histlist,gSignal,gSignalfs,gSignalBC1,gSignalBC2,inoutanis,typeb[iBkgFcn],minMassForFit[iMin],maxMassForFit[iMax],rebin[iReb]);

          TGraphAsymmErrors *gv2=Computev2(gSignal,resol,averagePt,inoutanis,0x0);
          TGraphAsymmErrors *gv2fs=Computev2(gSignalfs,resol,averagePt,inoutanis,0x0);
          TGraphAsymmErrors *gv2BC1=Computev2(gSignalBC1,resol,averagePt,inoutanis,0x0);
          TGraphAsymmErrors *gv2BC2=Computev2(gSignalBC2,resol,averagePt,inoutanis,0x0);

          for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
            Double_t v2, v2err, pt;
            gv2->GetPoint(iPt,pt,v2);
            v2err=gv2->GetErrorY(iPt);
            hv2[iPt]->Fill(v2);
            hv2VsTrial[iPt]->SetBinContent(iTrial+1,v2);
            hv2VsTrial[iPt]->SetBinError(iTrial+1,v2err);
            gv2fs->GetPoint(iPt,pt,v2);
            v2err=gv2fs->GetErrorY(iPt);
            hv2fs[iPt]->Fill(v2);
            hv2fsVsTrial[iPt]->SetBinContent(iTrial+1,v2);
            hv2fsVsTrial[iPt]->SetBinError(iTrial+1,v2err);
            gv2BC1->GetPoint(iPt,pt,v2);
            v2err=gv2BC1->GetErrorY(iPt);
            hv2BC1[iPt]->Fill(v2);
            hv2BC1VsTrial[iPt]->SetBinContent(iTrial+1,v2);
            hv2BC1VsTrial[iPt]->SetBinError(iTrial+1,v2err);
            gv2BC2->GetPoint(iPt,pt,v2);
            v2err=gv2BC2->GetErrorY(iPt);
            hv2BC2VsTrial[iPt]->SetBinContent(iTrial+1,v2);
            hv2BC2VsTrial[iPt]->SetBinError(iTrial+1,v2err);
            hv2BC2[iPt]->Fill(v2);
            
            for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
              Double_t Y,phi;
              gSignal[iPt]->GetPoint(iPhi,phi,Y);
              hRawYield[iPhi][iPt]->Fill(Y);
              gSignalfs[iPt]->GetPoint(iPhi,phi,Y);
              hRawYieldfs[iPhi][iPt]->Fill(Y);
              gSignalBC1[iPt]->GetPoint(iPhi,phi,Y);
              hRawYieldBC1[iPhi][iPt]->Fill(Y);
              gSignalBC2[iPt]->GetPoint(iPhi,phi,Y);
              hRawYieldBC2[iPhi][iPt]->Fill(Y);
            }
          }
          iTrial++;
        }
      }
    }
  }

  SetStyle();
  
  //Prepare output file
  TFile *fout=new TFile(Form("v2RawYieldSyst_%d_%d_%s_%s.root",minCent,maxCent,aniss.Data(),suffix.Data()),"RECREATE");
  
  TPaveStats **pv2 = new TPaveStats*[nptbinsnew];
  TPaveStats **pv2fs = new TPaveStats*[nptbinsnew];
  TPaveStats **pv2BC1 = new TPaveStats*[nptbinsnew];
  TPaveStats **pv2BC2 = new TPaveStats*[nptbinsnew];

  TPaveStats ***pRawYield = new TPaveStats**[nphi];
  TPaveStats ***pRawYieldfs = new TPaveStats**[nphi];
  TPaveStats ***pRawYieldBC1 = new TPaveStats**[nphi];
  TPaveStats ***pRawYieldBC2 = new TPaveStats**[nphi];
  
  for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
    pRawYield[iPhi] = new TPaveStats*[nptbinsnew];
    pRawYieldfs[iPhi] = new TPaveStats*[nptbinsnew];
    pRawYieldBC1[iPhi] = new TPaveStats*[nptbinsnew];
    pRawYieldBC2[iPhi] = new TPaveStats*[nptbinsnew];
  }
  
  for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
    cv2->cd(iPt+1);
    hv2[iPt]->Draw();
    hv2fs[iPt]->Draw("sames");
    hv2BC1[iPt]->Draw("sames");
    hv2BC2[iPt]->Draw("sames");
    cv2->cd(iPt+1)->Update();
    pv2[iPt] = (TPaveStats*)hv2[iPt]->FindObject("stats");
    pv2[iPt]->SetTextColor(colors[0]);
    pv2[iPt]->SetY1NDC(0.74);
    pv2[iPt]->SetY2NDC(0.89);
    pv2fs[iPt] = (TPaveStats*)hv2fs[iPt]->FindObject("stats");
    pv2fs[iPt]->SetTextColor(colors[1]);
    pv2fs[iPt]->SetY1NDC(0.59);
    pv2fs[iPt]->SetY2NDC(0.74);
    pv2BC1[iPt] = (TPaveStats*)hv2BC1[iPt]->FindObject("stats");
    pv2BC1[iPt]->SetTextColor(colors[2]);
    pv2BC1[iPt]->SetY1NDC(0.44);
    pv2BC1[iPt]->SetY2NDC(0.59);
    pv2BC2[iPt] = (TPaveStats*)hv2BC2[iPt]->FindObject("stats");
    pv2BC2[iPt]->SetTextColor(colors[3]);
    pv2BC2[iPt]->SetY1NDC(0.29);
    pv2BC2[iPt]->SetY2NDC(0.44);
    cv2->cd(iPt+1)->Modified();
    
    cv2VsTrial->cd(iPt+1);
    hv2VsTrial[iPt]->GetYaxis()->SetRangeUser(-0.3,0.7);
    hv2fsVsTrial[iPt]->GetYaxis()->SetRangeUser(-0.3,0.7);
    hv2BC1VsTrial[iPt]->GetYaxis()->SetRangeUser(-0.3,0.7);
    hv2BC2VsTrial[iPt]->GetYaxis()->SetRangeUser(-0.3,0.7);
    hv2VsTrial[iPt]->Draw();
    hv2fsVsTrial[iPt]->Draw("same");
    hv2BC1VsTrial[iPt]->Draw("same");
    hv2BC2VsTrial[iPt]->Draw("same");
    
    fout->cd();
    hv2[iPt]->Write();
    hv2fs[iPt]->Write();
    hv2BC1[iPt]->Write();
    hv2BC2[iPt]->Write();
    hv2VsTrial[iPt]->Write();
    hv2fsVsTrial[iPt]->Write();
    hv2BC1VsTrial[iPt]->Write();
    hv2BC2VsTrial[iPt]->Write();
    
    for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
      cRawYield[iPhi]->cd(iPt+1);
      hRawYield[iPhi][iPt]->Draw();
      hRawYieldfs[iPhi][iPt]->Draw("sames");
      hRawYieldBC1[iPhi][iPt]->Draw("sames");
      hRawYieldBC2[iPhi][iPt]->Draw("sames");
      cRawYield[iPhi]->cd(iPt+1)->Update();
      pRawYield[iPhi][iPt] = (TPaveStats*)hRawYield[iPhi][iPt]->FindObject("stats");
      pRawYield[iPhi][iPt]->SetTextColor(colors[0]);
      pRawYield[iPhi][iPt]->SetY1NDC(0.74);
      pRawYield[iPhi][iPt]->SetY2NDC(0.89);
      pRawYieldfs[iPhi][iPt] = (TPaveStats*)hRawYieldfs[iPhi][iPt]->FindObject("stats");
      pRawYieldfs[iPhi][iPt]->SetTextColor(colors[1]);
      pRawYieldfs[iPhi][iPt]->SetY1NDC(0.59);
      pRawYieldfs[iPhi][iPt]->SetY2NDC(0.74);
      pRawYieldBC1[iPhi][iPt] = (TPaveStats*)hRawYieldBC1[iPhi][iPt]->FindObject("stats");
      pRawYieldBC1[iPhi][iPt]->SetTextColor(colors[2]);
      pRawYieldBC1[iPhi][iPt]->SetY1NDC(0.44);
      pRawYieldBC1[iPhi][iPt]->SetY2NDC(0.59);
      pRawYieldBC2[iPhi][iPt] = (TPaveStats*)hRawYieldBC2[iPhi][iPt]->FindObject("stats");
      pRawYieldBC2[iPhi][iPt]->SetTextColor(colors[3]);
      pRawYieldBC2[iPhi][iPt]->SetY1NDC(0.29);
      pRawYieldBC2[iPhi][iPt]->SetY2NDC(0.44);
      cRawYield[iPhi]->cd(iPt+1)->Modified();
      
      fout->cd();
      hRawYield[iPhi][iPt]->Write();
      hRawYieldfs[iPhi][iPt]->Write();
      hRawYieldBC1[iPhi][iPt]->Write();
      hRawYieldBC2[iPhi][iPt]->Write();
    }
  }
  fout->Close();
  
  cv2->SaveAs(Form("v2RawYieldSyst_%d_%d_%s_%s.pdf",minCent,maxCent,aniss.Data(),suffix.Data()));
  for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
    cRawYield[iPhi]->SaveAs(Form("RawYieldSyst_phi%d_%d_%d_%s_%s.pdf",iPhi,minCent,maxCent,aniss.Data(),suffix.Data()));
  }
}

//______________________________________________________________
Double_t GetEventPlaneResolution(Double_t& error, TH1F* hsubev1, TH1F* hsubev2, TH1F* hsubev3){
  Double_t resolFull=1.;
  if(evPlaneRes==AliEventPlaneResolutionHandler::kTwoRandSub ||
     evPlaneRes==AliEventPlaneResolutionHandler::kTwoChargeSub){
    resolFull=AliVertexingHFUtils::GetFullEvResol(hsubev1);
    error = TMath::Abs(resolFull-AliVertexingHFUtils::GetFullEvResolLowLim(hsubev1));
  }else if(evPlaneRes==AliEventPlaneResolutionHandler::kTwoEtaSub){
    if(evPlane==AliEventPlaneResolutionHandler::kTPCFullEta){
      resolFull=AliVertexingHFUtils::GetFullEvResol(hsubev1);
      error = TMath::Abs(resolFull-AliVertexingHFUtils::GetFullEvResolLowLim(hsubev1));
    }else if(evPlane==AliEventPlaneResolutionHandler::kTPCPosEta){
      resolFull=AliVertexingHFUtils::GetSubEvResol(hsubev1);
      error = TMath::Abs(resolFull-AliVertexingHFUtils::GetSubEvResolLowLim(hsubev1));      
    }
  }else{
    Double_t resolSub[3];
    Double_t errors[3];
    TH1F* hevplresos[3];
    hevplresos[0]=hsubev1;
    hevplresos[1]=hsubev2;
    hevplresos[2]=hsubev3;
    for(Int_t ires=0;ires<3;ires++){
      resolSub[ires]=hevplresos[ires]->GetMean();
      errors[ires]=hevplresos[ires]->GetMeanError();
    }
    Double_t lowlim[3];for(Int_t ie=0;ie<3;ie++)lowlim[ie]=TMath::Abs(resolSub[ie]-errors[ie]);
    if(evPlane==AliEventPlaneResolutionHandler::kVZEROC ||
       evPlane==AliEventPlaneResolutionHandler::kVZERO){
      resolFull=TMath::Sqrt(resolSub[1]*resolSub[2]/resolSub[0]);
      error=resolFull-TMath::Sqrt(lowlim[2]*lowlim[1]/lowlim[0]);
    }
    else if(evPlane==AliEventPlaneResolutionHandler::kVZEROA){
      resolFull=TMath::Sqrt(resolSub[0]*resolSub[2]/resolSub[1]);
      error=resolFull-TMath::Sqrt(lowlim[2]*lowlim[0]/lowlim[1]);
    }
    else if(evPlane==AliEventPlaneResolutionHandler::kTPCFullEta ||
	    evPlane==AliEventPlaneResolutionHandler::kTPCPosEta){
      resolFull=TMath::Sqrt(resolSub[0]*resolSub[1]/resolSub[2]);
      error=resolFull-TMath::Sqrt(lowlim[0]*lowlim[1]/lowlim[2]);
    }
  }
  return resolFull;
}

//____________________________________________________________________
TList* LoadResolutionHistos(TList *inputlist){

  TList *outlist = new TList();
  outlist->SetName("eventplanehistlist");

  const Int_t nBins=20;
  Double_t ncoll[nBins]={1790.77,1578.44,1394.82,1236.17
    ,1095.08,969.86,859.571,759.959,669.648,589.588,516.039
    ,451.409,392.853,340.493,294.426,252.385,215.484,183.284
    ,155.101,130.963};

  Int_t minCentTimesTen=minCent*10;
  Int_t maxCentTimesTen=maxCent*10;
  TGraphErrors* gResolVsCent=new TGraphErrors(0);
  Int_t iPt=0;
  Int_t nSubRes=1;
  if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
     evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
  TString namereso[3]={"Reso","Reso2","Reso3"};
  TString suffixcentr=Form("centr%d_%d",minCent,maxCent);
  TH2F* hevpls=(TH2F*)inputlist->FindObject(Form("hEvPlanecentr%d_%d",minCentTimesTen,minCentTimesTen+25));
  hevpls->SetName(Form("hEvPlane%s",suffixcentr.Data()));
  hevpls->SetTitle(Form("Event Plane angle %s",suffixcentr.Data()));
  TH1F* hevplresos[3];
  Int_t ncBin=minCentTimesTen/25;
  
  for(Int_t ires=0;ires<nSubRes;ires++){
    hevplresos[ires]=(TH1F*)inputlist->FindObject(Form("hEvPlane%scentr%d_%d",namereso[ires].Data(),minCentTimesTen,minCentTimesTen+25));
    if(hevplresos[ires]){
      hevplresos[ires]->SetName(Form("hEvPlane%s%s",namereso[ires].Data(),suffixcentr.Data()));
      hevplresos[ires]->SetTitle(Form("Event Plane Resolution %s%s",namereso[ires].Data(),suffixcentr.Data()));
      if(useNcollWeight){
	printf("Centr %d Bin %d  Ncoll %f\n",minCentTimesTen,ncBin,ncoll[ncBin]);
	hevplresos[ires]->Scale(ncoll[ncBin]);
      }
    }
  }
  Double_t error;
  Double_t lowestRes=1;
  Double_t highestRes=0;
  Double_t resolBin=GetEventPlaneResolution(error,hevplresos[0],hevplresos[1],hevplresos[2]);
  if(resolBin<lowestRes) lowestRes=resolBin;
  if(resolBin>highestRes) highestRes=resolBin;

  Double_t binHalfWid=25./20.;
  Double_t binCentr=(Double_t)minCentTimesTen/10.+binHalfWid;
  gResolVsCent->SetPoint(iPt,binCentr,resolBin);
  gResolVsCent->SetPointError(iPt,binHalfWid,error);
  ++iPt;
  
  for(Int_t icentr=minCentTimesTen+25;icentr<maxCentTimesTen;icentr=icentr+25){
    TH2F* h=(TH2F*)inputlist->FindObject(Form("hEvPlanecentr%d_%d",icentr,icentr+25));
    if(h)hevpls->Add(h);
    else cout<<"skipping ev plane "<<icentr<<"_"<<icentr+5<<endl;
    TH1F* htmpresos[3];
    for(Int_t ires=0;ires<nSubRes;ires++){
      htmpresos[ires]=(TH1F*)inputlist->FindObject(Form("hEvPlane%scentr%d_%d",namereso[ires].Data(),icentr,icentr+25));
      if(!htmpresos[ires])cout<<"skipping ev pl reso "<<icentr<<"_"<<icentr+25<<endl;
    }
    resolBin=GetEventPlaneResolution(error,htmpresos[0],htmpresos[1],htmpresos[2]);
    if(resolBin<lowestRes) lowestRes=resolBin;
    if(resolBin>highestRes) highestRes=resolBin;
    binCentr=(Double_t)icentr/10.+binHalfWid;
    gResolVsCent->SetPoint(iPt,binCentr,resolBin);
    gResolVsCent->SetPointError(iPt,binHalfWid,error);
    ++iPt;
    ncBin=icentr/25;
    for(Int_t ires=0;ires<nSubRes;ires++){
      if(htmpresos[ires]){
	if(useNcollWeight){
	  printf("Centr %d Bin %d  Ncoll %f\n",icentr,ncBin,ncoll[ncBin]);
	  htmpresos[ires]->Scale(ncoll[ncBin]);
	}
	hevplresos[ires]->Add(htmpresos[ires]);
      }
    }
  }
  outlist->Add(hevpls->Clone());
  for(Int_t ires=0;ires<nSubRes;ires++){
    if(hevplresos[ires]) outlist->Add(hevplresos[ires]->Clone());
  }
  gResolVsCent->SetName("gResolVsCent");
  gResolVsCent->SetTitle("Resolution vs. Centrality");
  outlist->Add(gResolVsCent->Clone());
  return outlist;
}

//__________________________________________________________
TList *LoadMassHistos(TList *inputlist,Bool_t inoutanis){
  // printf("Start load histos\n");
  //  const Int_t nptbins=cutobj->GetNPtBins();
  TList *outlist = new TList();
  outlist->SetName("azimuthalhistoslist");
  
  Int_t nphi=nphibins;
  if(inoutanis)nphi=2;
  Int_t minCentTimesTen=minCent*10;
  Int_t maxCentTimesTen=maxCent*10;
  
  //Create 2D histogram in final pt bins
  for(Int_t iFinalPtBin=0; iFinalPtBin<nptbinsnew; iFinalPtBin++){
    for(Int_t iphi=0;iphi<nphi;iphi++){
      TH1F *hMass=0x0;//=new TH1F();
      for(Int_t iPtBin=minPtBin[iFinalPtBin]; iPtBin<=maxPtBin[iFinalPtBin];iPtBin++){
	for(Int_t iHisC=minCentTimesTen; iHisC<=maxCentTimesTen-25; iHisC+=25){    
	  TString hisname=Form("hMdeltaphi_pt%dcentr%d_%d",iPtBin,iHisC,iHisC+25);
	  TH2F* htmp=(TH2F*)inputlist->FindObject(hisname.Data());
	  printf("---> Histogram: %s\n",htmp->GetName());
	  Int_t startX=htmp->FindBin(phibinslim[iphi]);
	  Int_t endX=htmp->FindBin(phibinslim[iphi+1]-0.0001); // -0.0001 to be sure that the upper limit of the bin is properly set
	  TH1F *h1tmp;
	  if(inoutanis){
	    if(iphi==0){
	      Int_t firstBin=htmp->FindBin(0);
	      Int_t lastBin=htmp->FindBin(TMath::Pi()/4.-0.0001); // -0.0001 to be sure that the upper limit of the bin is pi/4
	      h1tmp=(TH1F*)htmp->ProjectionY(Form("hMass%d_phi0",iPtBin),firstBin,lastBin);
	      printf("In-plane, Range: bins %d-%d -> phi %f - %f\n",firstBin,lastBin,htmp->GetXaxis()->GetBinLowEdge(firstBin),htmp->GetXaxis()->GetBinUpEdge(lastBin));
	      firstBin=htmp->FindBin(3.*TMath::Pi()/4.);
	      lastBin=htmp->FindBin(TMath::Pi()-0.0001); // -0.0001 to be sure that the upper limit of the bin is pi
	      h1tmp->Add((TH1F*)htmp->ProjectionY(Form("hMass%d",iPtBin),firstBin,lastBin));
	      printf("In-plane, Range: bins %d-%d -> phi %f - %f\n",firstBin,lastBin,htmp->GetXaxis()->GetBinLowEdge(firstBin),htmp->GetXaxis()->GetBinUpEdge(lastBin));
	    }else{
	      Int_t firstBin=htmp->FindBin(TMath::Pi()/4.);
	      Int_t lastBin=htmp->FindBin(3.*TMath::Pi()/4.-0.0001); // -0.0001 to be sure that the upper limit of the bin is pi/4
	      h1tmp=(TH1F*)htmp->ProjectionY(Form("hMass%d_phi1",iPtBin),firstBin,lastBin);
	      printf("Out-of-plane, Range: bins %d-%d -> phi %f - %f\n",firstBin,lastBin,htmp->GetXaxis()->GetBinLowEdge(firstBin),htmp->GetXaxis()->GetBinUpEdge(lastBin));
	    }
	  }else{
	    h1tmp=(TH1F*)htmp->ProjectionY(Form("hMass%d_phi%d",iPtBin,iphi),startX,endX);
	  }
	  if(hMass==0)hMass=(TH1F*)h1tmp->Clone();
	  else hMass->Add((TH1F*)h1tmp->Clone());
	}
      }
      hMass->SetTitle(Form("hMass_pt%d_phi%d",iFinalPtBin,iphi));
      hMass->SetName(Form("hMass_pt%d_phi%d",iFinalPtBin,iphi));
      outlist->Add(hMass->Clone());
      delete hMass;
      hMass=0x0;
    }
  }
  return outlist;
}
//______________________________________________________________
Bool_t DefinePtBins(AliRDHFCuts *cutobj){
  Int_t nPtBinsCuts=cutobj->GetNPtBins();
  Float_t *ptlimsCuts=(Float_t*)cutobj->GetPtBinLimits();
  for(Int_t iPtCuts=0; iPtCuts<nPtBinsCuts; iPtCuts++){
    for(Int_t iFinalPtBin=0; iFinalPtBin<nptbinsnew; iFinalPtBin++){  
      if(TMath::Abs(ptlimsCuts[iPtCuts]-ptbinsnew[iFinalPtBin])<0.0001){ 
        minPtBin[iFinalPtBin]=iPtCuts;
        if(iFinalPtBin>0) maxPtBin[iFinalPtBin-1]=iPtCuts-1;
      }
    }
    if(TMath::Abs(ptlimsCuts[iPtCuts]-ptbinsnew[nptbinsnew])<0.0001) maxPtBin[nptbinsnew-1]=iPtCuts-1;
  }
  if(TMath::Abs(ptbinsnew[nptbinsnew]-ptlimsCuts[nPtBinsCuts])<0.0001) maxPtBin[nptbinsnew-1]=nPtBinsCuts-1;
  for(Int_t iFinalPtBin=0; iFinalPtBin<nptbinsnew; iFinalPtBin++){
    printf("Pt bins to be merged: %d %d\n",minPtBin[iFinalPtBin],maxPtBin[iFinalPtBin]);
    if(minPtBin[iFinalPtBin]<0 || maxPtBin[iFinalPtBin]<0) return kFALSE;
  }

  return kTRUE;
}
//______________________________________________________________
Int_t GetPadNumber(Int_t ix,Int_t iy){
  return (iy)*nptbinsnew+ix+1;
}
//______________________________________________________________
void FillSignalGraph(TList *histlist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2, Bool_t inoutanis, Int_t bkgfunc, Int_t minfit, Int_t maxfit, Int_t rebin){

  Int_t nphi=nphibins;
  if(inoutanis)nphi=2;
  
  //Canvases for drawing histograms
  Int_t nMassBins;
  for(Int_t ipt=0;ipt<nptbinsnew;ipt++){
    TH1F *histtofitfullsigma=(TH1F*)histlist->FindObject(Form("hMass_pt%d_phi0",ipt))->Clone();
    for(Int_t iphi=0;iphi<nphi;iphi++){
      Int_t ipad=GetPadNumber(ipt,iphi);
      Double_t signal=0,esignal=0;
      TH1F *histtofit=(TH1F*)histlist->FindObject(Form("hMass_pt%d_phi%d",ipt,iphi))->Clone();
      if(iphi>0)histtofitfullsigma->Add((TH1F*)histtofit->Clone());
      if(!histtofit){
        gSignal[ipt]->SetPoint(iphi,iphi,signal);
        gSignal[ipt]->SetPointError(iphi,0,0,esignal,esignal);
        return;
      }
      histtofit->SetTitle(Form("%.0f < #it{p}_{T} < %.0f, #phi%d",ptbinsnew[ipt],ptbinsnew[ipt+1],iphi));
      nMassBins=histtofit->GetNbinsX();
      histtofit->Rebin(rebin);
      AliHFMassFitter fitter(histtofit,minfit,maxfit,1,bkgfunc);
      fitter.SetInitialGaussianMean(massD);
      fitter.SetInitialGaussianSigma(0.012);
      Bool_t ok=fitter.MassFitter(kFALSE);
      Double_t sigmaforcounting=0;
      Double_t meanforcounting=0;
      if(ok){
        fitter.Signal(3,signal,esignal);
        sigmaforcounting=fitter.GetSigma();
        meanforcounting=fitter.GetMean();
      }
      gSignal[ipt]->SetPoint(iphi,iphi,signal);
      gSignal[ipt]->SetPointError(iphi,0,0,esignal,esignal);
      TF1* fB1=fitter.GetBackgroundFullRangeFunc();
      TF1* fB2=fitter.GetBackgroundRecalcFunc();
      Double_t minBinSum=histtofit->FindBin(meanforcounting-nSigmaForCounting*sigmaforcounting);
      Double_t maxBinSum=histtofit->FindBin(meanforcounting+nSigmaForCounting*sigmaforcounting);
      Double_t cntSig1=0.;
      Double_t cntSig2=0.;
      Double_t cntErr=0.;
      for(Int_t iMB=minBinSum; iMB<=maxBinSum; iMB++){
        Double_t bkg1=fB1 ? fB1->Eval(histtofit->GetBinCenter(iMB)) : 0;
        Double_t bkg2=fB2 ? fB2->Eval(histtofit->GetBinCenter(iMB)) : 0;
        cntSig1+=(histtofit->GetBinContent(iMB)-bkg1);
        cntSig2+=(histtofit->GetBinContent(iMB)-bkg2);
        cntErr+=(histtofit->GetBinContent(iMB));
      }
      cntErr=TMath::Sqrt(cntErr);
      gSignalBC1[ipt]->SetPoint(iphi,iphi,cntSig1);
      gSignalBC1[ipt]->SetPointError(iphi,0,0,cntErr,cntErr);
      gSignalBC2[ipt]->SetPoint(iphi,iphi,cntSig2);
      gSignalBC2[ipt]->SetPointError(iphi,0,0,cntErr,cntErr);
    }
    //fit for fixed sigma
    histtofitfullsigma->SetTitle(Form("%.0f < #it{p}_{T} < %.0f GeV/c",ptbinsnew[ipt],ptbinsnew[ipt+1]));
    histtofitfullsigma->GetXaxis()->SetTitle("M_{K#pi#pi} (GeV/c^{2})");
    histtofitfullsigma->GetXaxis()->SetTitleSize(0.05);
    nMassBins=histtofitfullsigma->GetNbinsX();
    histtofitfullsigma->Rebin(rebin);
    AliHFMassFitter fitter(histtofitfullsigma,minfit,maxfit,1,bkgfunc);
    fitter.SetInitialGaussianMean(massD);
    Bool_t ok=fitter.MassFitter(kFALSE);
    Double_t sigma=fitter.GetSigma();
    Double_t massFromFit=fitter.GetMean();
    for(Int_t iphi=0;iphi<nphi;iphi++){
      Int_t ipad=GetPadNumber(ipt,iphi);
      TH1F *histtofit=(TH1F*)histlist->FindObject(Form("hMass_pt%d_phi%d",ipt,iphi))->Clone();
      histtofit->SetTitle(Form("%.1f<#it{p}_{T}<%.1f, #phi%d",ptbinsnew[ipt],ptbinsnew[ipt+1],iphi));
      nMassBins=histtofit->GetNbinsX();
      histtofit->Rebin(rebin);
      AliHFMassFitter fitter2(histtofit,minfit,maxfit,1,bkgfunc);
      fitter2.SetInitialGaussianMean(massD);
      fitter2.SetFixGaussianSigma(sigma);
      if(fixAlsoMass) fitter2.SetFixGaussianMean(massFromFit);
      Bool_t ok2=fitter2.MassFitter(kFALSE);
      Double_t signal=0,esignal=0;
      if(ok2){
        fitter2.Signal(3,signal,esignal);
      }
      gSignalfs[ipt]->SetPoint(iphi,iphi,signal);
      gSignalfs[ipt]->SetPointError(iphi,0,0,esignal,esignal);
    }
  }//end loop on pt bin
}

TGraphAsymmErrors* Computev2(TGraphAsymmErrors **gSignal, Double_t resol, Float_t *averagePt, Bool_t inoutanis, TGraphAsymmErrors *gRelSystEff) {
  
  TGraphAsymmErrors* gv2 = new TGraphAsymmErrors(nptbinsnew);
  
  if(inoutanis) {
    for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
      Double_t *y=gSignal[iPt]->GetY();
      Double_t nIn=y[0];
      Double_t nOut=y[1];
      Double_t enIn=gSignal[iPt]->GetErrorY(0);
      Double_t enOut=gSignal[iPt]->GetErrorY(1);
      Double_t anis=(nIn-nOut)/(nIn+nOut);
      Double_t eAnis=2./((nIn+nOut)*(nIn+nOut))*TMath::Sqrt(nIn*nIn*enOut*enOut+nOut*nOut*enIn*enIn);
      Double_t v2=anis*TMath::Pi()/4./resol;
      Double_t ev2=eAnis*TMath::Pi()/4./resol;
      gv2->SetPoint(iPt,averagePt[iPt],v2);
      gv2->SetPointError(iPt,averagePt[iPt]-ptbinsnew[iPt],ptbinsnew[iPt+1]-averagePt[iPt],ev2,ev2);
      if(gRelSystEff) {
        //systematic uncertainty for in-out efficiency
        Double_t anis1=(nIn-nOut*effInOverEffOut)/(nIn+nOut*effInOverEffOut);
        Double_t anis2=(nIn*effInOverEffOut-nOut)/(nIn*effInOverEffOut+nOut);
        Double_t systEffUp=0.,systEffDown=0.;
        if(anis1>anis && anis1>anis2) systEffUp=anis1/anis;
        if(anis2>anis && anis2>anis1) systEffUp=anis2/anis;
        if(anis1<anis && anis1<anis2) systEffDown=anis1/anis;
        if(anis2<anis && anis2<anis1) systEffDown=anis2/anis;
        cout << Form(" Bin %d <pt>=%.3f  v2=%f+-%f systEff=%f %f\n",iPt,averagePt[iPt],v2,ev2,systEffUp*v2,systEffDown*v2)<<endl;
        gRelSystEff->SetPoint(iPt,averagePt[iPt],v2);
        gRelSystEff->SetPointError(iPt,0.4,0.4,v2*(1-systEffDown),v2*(systEffUp-1));
      }
    }
    return gv2;
  }
  else {
    TF1 *flowFunc = new TF1("flow","[0]*(1.+2.*[1]*TMath::Cos(2.*x))");
    for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
      //v2 from fit to Deltaphi distribution
      gSignal[iPt]->Fit(flowFunc);
      Double_t v2 = flowFunc->GetParameter(1)/resol;
      Double_t ev2=flowFunc->GetParError(1)/resol;
      gv2->SetPoint(iPt,averagePt[iPt],v2);
      gv2->SetPointError(iPt,averagePt[iPt]-ptbinsnew[iPt],ptbinsnew[iPt+1]-averagePt[iPt],ev2,ev2);
    }
    return gv2;
  }
}

//___________________________________________________________
Int_t FindPtBin(Int_t nbins, Double_t* array,Double_t value){
  for (Int_t i=0;i<nbins;i++){
    if(value>=array[i] && value<array[i+1]){
      return i;
    }
  }
  cout<<value<< " out of range "<<array[0]<<", "<<array[nbins]<<endl;
  return -1;
}

//___________________________________________________________
Int_t LoadRefGraphs(TString reffilename, TH1F **hRawYieldRef, TH1F **hRawYieldfsRef, TH1F **hRawYieldBC1Ref, TH1F **hRawYieldBC2Ref, TGraphAsymmErrors *&gv2Ref, TGraphAsymmErrors *&gv2fsRef, TGraphAsymmErrors *&gv2BC1Ref, TGraphAsymmErrors *&gv2BC2Ref, Bool_t inoutanis) {
  
  Int_t nphi=nphibins;
  if(inoutanis) nphi=2;
  
  TFile* reffile = TFile::Open(reffilename.Data(),"READ");
  if(reffile) {
    gv2Ref=(TGraphAsymmErrors*)reffile->Get("gav2");
    gv2fsRef=(TGraphAsymmErrors*)reffile->Get("gav2fs");
    gv2BC1Ref=(TGraphAsymmErrors*)reffile->Get("gav2BC1");
    gv2BC2Ref=(TGraphAsymmErrors*)reffile->Get("gav2BC2");
  
    if(!gv2Ref) {return 2;}
    if(!gv2fsRef) {return 3;}
    if(!gv2BC1Ref) {return 4;}
    if(!gv2BC2Ref) {return 5;}
    
    for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
      hRawYieldRef[iPhi] = new TH1F(Form("hRawYieldRef_phi%d",iPhi),"",nptbinsnew,ptbinsnew);
      hRawYieldfsRef[iPhi]= new TH1F(Form("hRawYieldfsRef_phi%d",iPhi),"",nptbinsnew,ptbinsnew);
      hRawYieldBC1Ref[iPhi]= new TH1F(Form("hRawYieldBC1Ref_phi%d",iPhi),"",nptbinsnew,ptbinsnew);
      hRawYieldBC2Ref[iPhi]= new TH1F(Form("hRawYieldBC2Ref_phi%d",iPhi),"",nptbinsnew,ptbinsnew);
      hRawYieldRef[iPhi]->SetDirectory(0);
      hRawYieldfsRef[iPhi]->SetDirectory(0);
      hRawYieldBC1Ref[iPhi]->SetDirectory(0);
      hRawYieldBC2Ref[iPhi]->SetDirectory(0);
    }
    
    for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
      TGraphAsymmErrors* gtmp = (TGraphAsymmErrors*)reffile->Get(Form("gasigpt%d",iPt));
      TGraphAsymmErrors* gtmpfs = (TGraphAsymmErrors*)reffile->Get(Form("gasigfspt%d",iPt));
      TGraphAsymmErrors* gtmpBC1 = (TGraphAsymmErrors*)reffile->Get(Form("gasigBC1pt%d",iPt));
      TGraphAsymmErrors* gtmpBC2 = (TGraphAsymmErrors*)reffile->Get(Form("gasigBC2pt%d",iPt));
      
      if(gtmp && gtmpfs && gtmpBC1 && gtmpBC2) {
        for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
          Double_t phi,Y, errY;
          gtmp->GetPoint(iPhi,phi,Y);
          hRawYieldRef[iPhi]->SetBinContent(iPt+1,Y);
          gtmpfs->GetPoint(iPhi,phi,Y);
          hRawYieldfsRef[iPhi]->SetBinContent(iPt+1,Y);
          gtmpBC1->GetPoint(iPhi,phi,Y);
          hRawYieldBC1Ref[iPhi]->SetBinContent(iPt+1,Y);
          gtmpBC2->GetPoint(iPhi,phi,Y);
          hRawYieldBC2Ref[iPhi]->SetBinContent(iPt+1,Y);
        }
      }
    }
    reffile->Close();
  }
  else {return 1;}
  
  return 0;
}

//___________________________________________________________
void DivideCanvas(TCanvas* c, Int_t nPtBins) {
  
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

//___________________________________________________________
void SetStyle(Int_t optfit) {
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(optfit);
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTextFont(42);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetStatFont(42);
  gStyle->SetStatY(0.89);
  gStyle->SetStatX(0.89);
  gStyle->SetTitleFont(42,"xyzg");
  gStyle->SetHistLineWidth(2);
  gStyle->SetLegendBorderSize(0);
}
