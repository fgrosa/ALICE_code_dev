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

#endif

enum {kDzero,kDplus};
enum {kLeft,kRight};

Int_t PIDaxnum=3;
Int_t minPID=0;
Int_t maxPID=4;

Int_t rebin[10] = {6,6,5,5,5,5,5,5,5,3};

Double_t minfit = 1.72;
Double_t maxfit = 2.02;
Int_t sgnfunc = AliHFMassFitter::kGaus;
Int_t bkgfunc = AliHFMassFitter::kExpo;

void GetEfficiencyAndRawYields(Int_t meson=kDplus,
                               TString axesfile="axes",
                               TString cutsfile="cutset1",
                               TString rawyieldfile="cutvar/rawyields_Dplus_cutset1",
                               TString rawyieldfileMC="cutvar/rawyieldsMC_Dplus_cutset1",
                               TString efffile="cutvar/efficiency_Dplus_cutset1",
                               Bool_t ptreweights=kFALSE,
                               Bool_t drawdist=kTRUE);

void PlotDistVsPt(THnSparseF* RecoPromptSparse,THnSparseF* RecoFDSparse, THnSparseF* DataSparse, Int_t axis,TString mesonname,TString axesfile,TString cutsfile);
TH1F* GetEfficiency(THnSparseF* RecoSparse, THnSparseF* MCGenAccSparse, TString histoname, TString axesfile, TString cutsfile, Int_t meson, Bool_t ptreweights, Bool_t chicut=kTRUE);
void FitRawYields(THnSparseF* sparse, Int_t meson, TString axesfile, TString cutsfile, Bool_t isMC, Double_t Nev,TString rawyieldfile, TString rawyieldfileMC);

void UsePID(THnSparseF* sparse);
void ApplyTopologicalCuts(THnSparseF* sparse, Int_t iPt, TString axesfile, TString cutsfile, Bool_t chicut=kTRUE);
void SetPtRange(THnSparseF* sparse, Int_t iPt, Int_t iAxis, TString axesfile, TString cutsfile);
void ResetAxes(THnSparseF* sparse);
void ReadAxesNum(TString Filename, vector<string> &axesanmes, vector<int> &axesno);
void ReadSet(TString Filename, vector<string> &varnames, vector<double> &cutset);
Int_t GetNPtBins(TString axesfile, TString cutsfile);
Double_t* GetPtLims(TString axesfile, TString cutsfile);
Double_t PtWeightsFromFONLL5overLHC13d3(Double_t *x, Double_t *pars);
Double_t PtBWeightsFromFONLL5overLHC13d3(Double_t *x, Double_t *pars);
void PrintCuts(TString axesfile, TString cutsfile);
void SetSideBandsRegion(THnSparseF* sparse, Int_t side, Int_t iPt);

void GetEfficiencyAndRawYields(Int_t meson,TString axesfile,TString cutsfile,TString rawyieldfile,TString rawyieldfileMC,TString efffile,Bool_t ptreweights, Bool_t drawdist) {
  
  TString mesonname;
  if(meson==kDzero) mesonname = "Dzero";
  else mesonname = "Dplus";
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetTitleSize(0.045,"xyz");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetLegendBorderSize(0);
    
  //INPUT FILES________________________________________________________________
  TString MCfilename;
  TString Datafilename;
  TString listname;
  
  if(meson==kDzero) {
    MCfilename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pPb_LHC13/MC/fgrosa_Dzero_KF_pPbMC.root";
    Datafilename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pPb_LHC13/Data/fgrosa_Dzero_KF_pPb.root";
    listname = "coutputDzeroKF";
  }
  else {
    MCfilename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pPb_LHC13/MC/fgrosa_Dplus_KF_pPbMC.root";
    Datafilename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pPb_LHC13/Data/fgrosa_Dplus_KF_pPb.root";
    listname = "coutputDplusKF";
  }

  //PRINT CUTS VALUES__________________________________________________________
  PrintCuts(axesfile,cutsfile);
  
  TFile infileMC(MCfilename.Data(),"UPDATE");
  TList* listMC = (TList*)infileMC.Get(listname.Data());  
  THnSparseF *RecoPromptSparse = (THnSparseF*)listMC->FindObject("fSparsePrompt");
  THnSparseF *MCPromptSparse = (THnSparseF*)listMC->FindObject("fMCGenAccPrompt");
  THnSparseF *RecoFDSparse = (THnSparseF*)listMC->FindObject("fSparseFD");
  THnSparseF *MCFDSparse = (THnSparseF*)listMC->FindObject("fMCGenAccFD");
  infileMC.Close();

  TFile infileData(Datafilename.Data(),"UPDATE");
  TList* listData = (TList*)infileData.Get(listname.Data());
  THnSparseF *AllSparse = (THnSparseF*)listData->FindObject("fSparseAll");
  TH1F *hEv = (TH1F*)listData->FindObject("fHistNEvents");

  //CHI SQUARE DISTRIBUTIONS___________________________________________________
  if(drawdist) {
    PlotDistVsPt(RecoPromptSparse,RecoFDSparse,AllSparse,12,mesonname,axesfile,cutsfile);
    PlotDistVsPt(RecoPromptSparse,RecoFDSparse,AllSparse,13,mesonname,axesfile,cutsfile);
  }
  
  //EFFICIENCY_________________________________________________________________
  TH1F* hEffPromptChiCut = (TH1F*)GetEfficiency(RecoPromptSparse,MCPromptSparse,"hEffD",axesfile,cutsfile,meson,ptreweights,kTRUE);
  TH1F* hEffFDChiCut = (TH1F*)GetEfficiency(RecoFDSparse,MCFDSparse,"hEffB",axesfile,cutsfile,meson,ptreweights,kTRUE);
  hEffPromptChiCut->SetLineColor(kBlue);
  hEffFDChiCut->SetLineColor(kRed);
  hEffPromptChiCut->SetMarkerColor(kBlue);
  hEffFDChiCut->SetMarkerColor(kRed);
  hEffPromptChiCut->SetMarkerStyle(20);
  hEffFDChiCut->SetMarkerStyle(20);
  
  TLegend *l = new TLegend(0.5,0.2,0.89,0.45);
  l->SetTextSize(0.04);
  l->SetTextFont(42);
  l->AddEntry(hEffPromptChiCut,"Prompt","ple");
  l->AddEntry(hEffFDChiCut,"Feed-down","ple");

  TCanvas* cEff = new TCanvas("cEff","Efficiency",800,800);
  cEff->SetLogy();
  hEffPromptChiCut->GetYaxis()->SetRangeUser(hEffPromptChiCut->GetMinimum()*0.8, hEffPromptChiCut->GetMaximum()*2);
  hEffPromptChiCut->GetXaxis()->SetTitleFont(42);
  hEffPromptChiCut->GetYaxis()->SetTitleFont(42);
  hEffPromptChiCut->GetXaxis()->SetLabelFont(42);
  hEffPromptChiCut->GetYaxis()->SetLabelFont(42);
  hEffPromptChiCut->Draw("E");
  hEffFDChiCut->Draw("Esame");
  l->Draw("same");

  TString effplotfile=efffile;
  effplotfile.ReplaceAll("root","eps");
  cEff->SaveAs(effplotfile.Data());

  TFile outfile(efffile,"RECREATE");
  hEffPromptChiCut->Write();
  hEffFDChiCut->Write();
  outfile.Close();

  delete hEffPromptChiCut;
  delete hEffFDChiCut;
  
  //NUMBER OF EVENTS___________________________________________________________
  TH1F *hNEvents = new TH1F("hNEvents","",1,0,1);
  Int_t nRecoVertMax10 = hEv->GetBinContent(2); //corresponds to the number of analysed event
  Int_t nRecoVertMin10 = hEv->GetBinContent(5); //events rejected because zvtx > 10 cm
  Int_t nRecoVert = hEv->GetBinContent(3); //events with vtx reconstructed
  Int_t nNoVert = hEv->GetBinContent(4); //events without vtx reconstructed

  Double_t norm = nRecoVertMax10+nNoVert-nNoVert*(nRecoVertMin10/nRecoVert);
  
  //MC FIT FOR WIDTH COMPARISON________________________________________________
  FitRawYields(RecoPromptSparse,meson,axesfile,cutsfile,kTRUE,norm,rawyieldfile,rawyieldfileMC);

  //SIGNAL EXTRACTION__________________________________________________________
  FitRawYields(AllSparse,meson,axesfile,cutsfile,kFALSE,norm,rawyieldfile,rawyieldfileMC);
  
  hNEvents->SetBinContent(1,norm);
  TFile outfileraw(rawyieldfile,"UPDATE");
  hNEvents->Write();
  outfileraw.Close();

  delete hNEvents;
}

TH1F* GetEfficiency(THnSparseF* RecoSparse, THnSparseF* MCGenAccSparse, TString histoname, TString axesfile, TString cutsfile, Int_t meson, Bool_t ptreweights, Bool_t chicut) {

  const Int_t nPtBins = GetNPtBins(axesfile,cutsfile);
  Double_t *PtLims = GetPtLims(axesfile,cutsfile);
  TH1F* hEff = new TH1F(histoname.Data(),"",nPtBins,PtLims);
  hEff->SetDirectory(0);
  hEff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hEff->GetYaxis()->SetTitle("Efficiency");

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ResetAxes(RecoSparse);
    ResetAxes(MCGenAccSparse);
    SetPtRange(RecoSparse,iPt,1,axesfile,cutsfile);
    SetPtRange(MCGenAccSparse,iPt,0,axesfile,cutsfile);    
    UsePID(RecoSparse);
    ApplyTopologicalCuts(RecoSparse,iPt,axesfile,cutsfile,chicut);
    TH1F* hReco = (TH1F*)RecoSparse->Projection(1);
    TH1F* hMC = (TH1F*)MCGenAccSparse->Projection(0);
    
    if(ptreweights) {
      TF1* fFuncWeight = 0x0;
      TString name = RecoSparse->GetName();
      if(name.Contains("Prompt")) {
        fFuncWeight = new TF1("fFuncWeight",PtWeightsFromFONLL5overLHC13d3,0,40,9);
        fFuncWeight->SetParameters(2.94999e+00,3.47032e+00,2.81278e+00,2.5,1.93370e-02,3.86865e+00,-1.54113e-01,8.86944e-02,2.56267e-02);
      }
      else {
        fFuncWeight = new TF1("fFuncWeight",PtBWeightsFromFONLL5overLHC13d3,0,40,8);
        fFuncWeight->SetParameters(0.983543,2.72857,1.49631,-82.4067,51.359,5.05355,0.00136537,2.13765e-05);
      }
      
      for(Int_t iPt=0; iPt<hReco->GetNbinsX(); iPt++) {
        hReco->SetBinContent(iPt+1,hReco->GetBinContent(iPt+1)*fFuncWeight->Eval(hReco->GetXaxis()->GetBinCenter(iPt+1)));
        hMC->SetBinContent(iPt+1,hMC->GetBinContent(iPt+1)*fFuncWeight->Eval(hMC->GetXaxis()->GetBinCenter(iPt+1)));
      }
      delete fFuncWeight;
    }

    Double_t counts[2] = {hMC->GetEntries(),hReco->GetEntries()};

    if(meson==kDzero && (maxPID>3 || minPID<2)) {
      TAxis* PIDax = (TAxis*)RecoSparse->GetAxis(PIDaxnum);
      PIDax->SetRange(2,3);
      TH1F *h23 = (TH1F*)RecoSparse->Projection(1);
      Double_t count23 = h23->Integral();
      PIDax->SetRange(1,1);
      TH1F *h1 = (TH1F*)RecoSparse->Projection(1);
      Double_t count1 = h1->Integral()/2;
      PIDax->SetRange(4,4);
      TH1F *h4 = (TH1F*)RecoSparse->Projection(1);
      Double_t count4 = h4->Integral()/2;
    }
  
    TH1F *hTempNum=new TH1F("hTempNum","hTempNum",1,0,1);
    hTempNum->SetBinContent(1,counts[1]);
    hTempNum->SetBinError(1,TMath::Sqrt(counts[1]));
    
    TH1F *hTempDeNum=new TH1F("hTempDeNum","hTempDeNum",1,0,1);
    hTempDeNum->SetBinContent(1,counts[0]);
    hTempDeNum->SetBinError(1,TMath::Sqrt(counts[0]));

    TH1F* hEffPt = new TH1F("hEffPt","",1,0,1);
    hEffPt->Divide(hTempNum,hTempDeNum,1.,1.,"B");
    hEffPt->SetDirectory(0);
    hEffPt->SetStats(0);

    hEff->SetBinContent(iPt+1,hEffPt->GetBinContent(1));
    hEff->SetBinError(iPt+1,hEffPt->GetBinError(1));

    delete hTempNum;
    delete hTempDeNum;
    delete hEffPt;
  }
  
  return hEff;
}

void PlotDistVsPt(THnSparseF* RecoPromptSparse,THnSparseF* RecoFDSparse,THnSparseF* AllSparse, Int_t axis,TString  mesonname,TString axesfile,TString cutsfile) {

  vector<string> axesnames;
  vector<int> axesno;
  ReadAxesNum(axesfile,axesnames,axesno);
  Int_t iAxis=0;  
  while(axesno[iAxis]!=axis) {
    iAxis++;
  }
  TString varname = axesnames[iAxis];

  const Int_t nPtBins = GetNPtBins(axesfile,cutsfile);
  Double_t *PtLims = GetPtLims(axesfile,cutsfile);

  TCanvas **c= new TCanvas*[nPtBins];
  TLegend* l = new TLegend(0.6,0.7,0.89,0.89);
  l->SetTextSize(0.04);
  TFile outfile(Form("%s_%s.root",varname.Data(),mesonname.Data()),"RECREATE");
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ResetAxes(RecoPromptSparse);
    ResetAxes(RecoFDSparse);
    SetPtRange(RecoPromptSparse,iPt,1,axesfile,cutsfile);
    SetPtRange(RecoFDSparse,iPt,1,axesfile,cutsfile);
    UsePID(RecoPromptSparse);
    UsePID(RecoFDSparse);    
    TH1F* hPrompt = (TH1F*)RecoPromptSparse->Projection(axis);
    hPrompt->Scale(1./hPrompt->Integral());
    TH1F* hFD = (TH1F*)RecoFDSparse->Projection(axis);    
    hFD->Scale(1./hFD->Integral());

    ResetAxes(AllSparse);
    SetPtRange(AllSparse,iPt,1,axesfile,cutsfile);
    UsePID(AllSparse);
    SetSideBandsRegion(AllSparse,kLeft,iPt);
    TH1F* hSBLeft = (TH1F*)AllSparse->Projection(axis);
    ResetAxes(AllSparse);
    SetPtRange(AllSparse,iPt,1,axesfile,cutsfile);
    UsePID(AllSparse);
    SetSideBandsRegion(AllSparse,kRight,iPt);
    TH1F* hSBRight = (TH1F*)AllSparse->Projection(axis);
    TH1F* hBkg = (TH1F*)hSBLeft->Clone();
    hBkg->Add(hSBLeft,hSBRight,1.,1.);
    hBkg->Scale(1./hBkg->Integral());
      
    hPrompt->SetDirectory(0);
    hFD->SetDirectory(0);
    hBkg->SetDirectory(0);
    hPrompt->SetStats(0);
    hPrompt->SetLineColor(kBlue);
    hFD->SetLineColor(kRed);
    hBkg->SetLineColor(kGreen+3);
    hPrompt->GetYaxis()->SetTitle("Normalised entries");
    hFD->GetYaxis()->SetTitle("Normalised entries");
    hBkg->GetYaxis()->SetTitle("Normalised entries");
    hPrompt->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hFD->SetTitle("");
    hBkg->SetTitle("");
    c[iPt] = new TCanvas(Form("%s_pT%d",varname.Data(),iPt),"",800,800);
    c[iPt]->SetLogy();
    if(iPt==0) {
      l->AddEntry(hPrompt,"Prompt","l");
      l->AddEntry(hFD,"Feed-down","l");
      l->AddEntry(hBkg,"Background","l");
    }
    hPrompt->Draw();
    hFD->Draw("same");
    hBkg->Draw("same");
    l->Draw("same");
    Double_t y1 = hPrompt->GetMaximum()*0.12;
    Double_t y2 = hPrompt->GetMaximum()*0.06;
    if(iPt>7) {
      y1 = hPrompt->GetMaximum()*0.22;
      y2 = hPrompt->GetMaximum()*0.14;
    }
    if(iPt==9) {
     y1 = hPrompt->GetMaximum()*0.28;
     y2 = hPrompt->GetMaximum()*0.18;
    }
    outfile.cd();
    c[iPt]->Write();
  }
  outfile.Close();
  for(Int_t iPt=0; iPt<nPtBins; iPt++)
    delete c[iPt];
  delete[] c;
  delete l;
}

void FitRawYields(THnSparseF* sparse, Int_t meson, TString axesfile, TString cutsfile, Bool_t isMC, Double_t norm, TString rawyieldfile, TString rawyieldfileMC) {

  TString mesonname;
  if(meson==kDzero) mesonname = "Dzero";
  else mesonname = "Dplus";

  const Int_t nPtBins = GetNPtBins(axesfile,cutsfile);
  Double_t *PtLims = GetPtLims(axesfile,cutsfile);

  TCanvas** cMass = new TCanvas*[nPtBins];
  TCanvas** cRaw = new TCanvas*[nPtBins];
  TLegend* l = new TLegend(0.6,0.7,0.89,0.89);
  l->SetTextSize(0.04);

  TFile *outfile;
  if(isMC)
    outfile = TFile::Open(rawyieldfileMC,"RECREATE"); 
  else
    outfile = TFile::Open(rawyieldfile,"RECREATE"); 

  TH1F* hSignal = new TH1F("hSignal","",nPtBins,PtLims);
  TH1F* hSignalPerEvent = new TH1F("hSignalPerEvent","",nPtBins,PtLims);
  TH1F* hSigma = new TH1F("hSigma","",nPtBins,PtLims);
  TH1F* hSignificance = new TH1F("hSignificance","",nPtBins,PtLims);
  TH1F* hMean = new TH1F("hMean","",nPtBins,PtLims);
  hSignal->SetDirectory(0);
  hSignalPerEvent->SetDirectory(0);
  hSigma->SetDirectory(0);
  hMean->SetDirectory(0);
  hSignificance->SetDirectory(0);

  //open file for MC sigmas
  TH1F *hSigmaMC = 0x0;
  if(!isMC) {
    TFile *infileMC = TFile::Open(rawyieldfileMC.Data(),"READ");
    if(!infileMC)
      cerr << "File " << rawyieldfileMC << " not found!" << endl;
    else {
      hSigmaMC = (TH1F*)infileMC->Get("hSigma");
      hSigmaMC->SetDirectory(0);
      infileMC->Close();
    }
  }
  
  //PDG mass value
  Int_t pdgcode;
  if(meson==kDplus)
    pdgcode=411;
  else
    pdgcode=421;
  Double_t massD = TDatabasePDG::Instance()->GetParticle(pdgcode)->Mass();

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ResetAxes(sparse);
    SetPtRange(sparse,iPt,1,axesfile,cutsfile);
    UsePID(sparse);
    if(!isMC)
      ApplyTopologicalCuts(sparse,iPt,axesfile,cutsfile,kFALSE);
    TH1F* hMass = (TH1F*)sparse->Projection(0);
    hMass->SetDirectory(0);
    hMass->SetStats(0);
    hMass->SetName(Form("hMass_Pt%d",iPt));
    ResetAxes(sparse);
    SetPtRange(sparse,iPt,1,axesfile,cutsfile);
    UsePID(sparse);
    if(!isMC)
      ApplyTopologicalCuts(sparse,iPt,axesfile,cutsfile,kTRUE);
    TH1F* hMassChi2Cut = (TH1F*)sparse->Projection(0);
    hMassChi2Cut->SetDirectory(0);
    hMassChi2Cut->SetName(Form("hMassChi2Cut_Pt%d",iPt));
    hMassChi2Cut->SetLineColor(kRed);
    cMass[iPt] = new TCanvas(Form("cMass_Pt%d",iPt),"",1200,900);
    if(iPt==0) {
      l->AddEntry(hMass,"w/o #chi^{2} cut","l");
      l->AddEntry(hMassChi2Cut,"with #chi^{2} cut","l");
    }
    hMass->GetYaxis()->SetTitle(Form("Entries/(%0.f MeV/c^{2})",hMass->GetBinWidth(3)));
    hMass->Draw();
    hMass->GetYaxis()->SetRangeUser(hMassChi2Cut->GetMinimum()*0.7,hMass->GetMaximum()*1.1);
    hMassChi2Cut->Draw("same");
    l->Draw("same");
    //fit
    TH1F* hMassForFit = (TH1F*)hMassChi2Cut->Clone();
    hMassForFit->SetLineColor(kBlue);
    hMassForFit->SetStats(0);
    hMassForFit->Rebin(rebin[iPt]);
    hMassForFit->GetYaxis()->SetRangeUser(hMassChi2Cut->GetMinimum()*0.6*rebin[iPt],hMassChi2Cut->GetMaximum()*1.1*rebin[iPt]);
    if(!isMC) {
      AliHFMassFitter *fitter = new AliHFMassFitter(hMassForFit,minfit,maxfit,1,bkgfunc,sgnfunc);
      //set sigma (from MC) and mean (from PDG)
      if(hSigmaMC)
        fitter->SetFixGaussianSigma(hSigmaMC->GetBinContent(iPt+1));
      fitter->SetInitialGaussianMean(massD);
      fitter->SetUseLikelihoodFit();
      fitter->MassFitter(0);
      Double_t signal;
      Double_t signalerr;
      fitter->Signal(3,signal,signalerr);
      Double_t sigma = fitter->GetSigma();
      Double_t sigmaerr = fitter->GetSigmaUncertainty();
      Double_t mean = fitter->GetMean();
      Double_t meanerr = fitter->GetMeanUncertainty();
      Double_t significance;
      Double_t signiferr;
      fitter->Significance(3,significance,signiferr);
      hSignal->SetBinContent(iPt+1,signal);
      hSignal->SetBinError(iPt+1,signalerr);
      hSignalPerEvent->SetBinContent(iPt+1,signal/norm);
      hSignalPerEvent->SetBinError(iPt+1,signalerr/norm);
      hSigma->SetBinContent(iPt+1,sigma);
      hSigma->SetBinError(iPt+1,sigmaerr);
      hMean->SetBinContent(iPt+1,mean);
      hMean->SetBinError(iPt+1,meanerr);
      hSignificance->SetBinContent(iPt+1,significance);
      hSignificance->SetBinError(iPt+1,signiferr);      
      cRaw[iPt] = new TCanvas(Form("cRaw_Pt%d",iPt),"",1200,900);
      fitter->DrawHere(gPad);      
      //save
      outfile->cd();
      cMass[iPt]->Write();
      cRaw[iPt]->Write();
    }
    else {
      cRaw[iPt] = new TCanvas(Form("cMassPrompt_Pt%d",iPt),"",1200,900);
      hMassForFit->Fit("gaus");
      TF1 *f=(TF1*)hMassForFit->GetListOfFunctions()->FindObject("gaus");
      hSigma->SetBinContent(iPt+1,f->GetParameter(2));
      hSigma->SetBinError(iPt+1,f->GetParError(2));
      hMean->SetBinContent(iPt+1,f->GetParameter(1));
      hMean->SetBinError(iPt+1,f->GetParError(1));
      outfile->cd();
      cRaw[iPt]->Write();
    } 
  }
  
  if(!isMC) {
    hSignal->Write();
    hSignalPerEvent->Write();
    hSignificance->Write();
  }
  hSigma->Write();
  hMean->Write();
  outfile->Close();
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    delete cMass[iPt];
    delete cRaw[iPt];
  }
  delete[] cMass;
  delete[] cRaw;
  delete hSignal;
  delete hSignificance;
  delete hMean;
  delete hSignalPerEvent;
}

void UsePID(THnSparseF* sparse) {
  TAxis* PIDax = (TAxis*)sparse->GetAxis(PIDaxnum);
  PIDax->SetRange(minPID,maxPID);
}

void ApplyTopologicalCuts(THnSparseF* sparse, Int_t iPt,TString axesfile,TString cutsfile,Bool_t chicut) {
  vector<string> axesnames;
  vector<int> axesno;
  vector<string> cutvarnames;
  vector<double> cutset;
  ReadAxesNum(axesfile,axesnames,axesno);
  ReadSet(cutsfile,cutvarnames,cutset);
  const Int_t nVars = axesno.size();
  Int_t KFchiax = 0;
  for(Int_t iVar=0; iVar<nVars; iVar++) {
    KFchiax = axesno[iVar];
    if(axesnames[iVar]=="KFchi") break;
  }
  Int_t PVchiax = 0;
  for(Int_t iVar=0; iVar<nVars; iVar++) {
    PVchiax = axesno[iVar];
    if(axesnames[iVar]=="PVchi") break;
  }
  
  for(Int_t iVar=1; iVar<nVars; iVar++) {
    if(!chicut && (axesno[iVar]==KFchiax || axesno[iVar]==PVchiax)) {continue;} ///if chicut=kFALSE don't apply cut on chi square
    TAxis* ax = (TAxis*)sparse->GetAxis(axesno[iVar]);
    Int_t binmin = ax->FindBin(cutset[(2*iVar)+iPt*2*nVars]+0.0001);
    Int_t binmax = ax->FindBin(cutset[(2*iVar+1)+iPt*2*nVars]-0.0001);
    sparse->GetAxis(axesno[iVar])->SetRange(binmin,binmax);
  }
}

void SetPtRange(THnSparseF* sparse, Int_t iPt, Int_t iAxis,TString axesfile,TString cutsfile) {
 vector<string> axesnames;
  vector<int> axesno;
  vector<string> cutvarnames;
  vector<double> cutset;
  ReadAxesNum(axesfile,axesnames,axesno);
  ReadSet(cutsfile,cutvarnames,cutset);
  const Int_t nVars = axesno.size();

  TAxis* ptax = (TAxis*)sparse->GetAxis(iAxis);
  Int_t binmin = ptax->FindBin(cutset[iPt*2*nVars]+0.0001);
  Int_t binmax = ptax->FindBin(cutset[1+iPt*2*nVars]-0.0001);
  sparse->GetAxis(iAxis)->SetRange(binmin,binmax);
}

void ResetAxes(THnSparseF* sparse) {
  Int_t nAxes = sparse->GetNdimensions();
  for(Int_t iAxis=0; iAxis<nAxes; iAxis++) {
    sparse->GetAxis(iAxis)->SetRange(-1,-1);
  }
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

Int_t GetNPtBins(TString axesfile,TString cutsfile) {
  
  vector<string> axesnames;
  vector<int> axesno;
  vector<string> cutvarnames;
  vector<double> cutset;
  ReadAxesNum(axesfile,axesnames,axesno);
  ReadSet(cutsfile,cutvarnames,cutset);

  return cutset.size()/(axesno.size()*2);
}

Double_t* GetPtLims(TString axesfile,TString cutsfile) {
  
  vector<string> axesnames;
  vector<int> axesno;
  vector<string> cutvarnames;
  vector<double> cutset;
  ReadAxesNum(axesfile,axesnames,axesno);
  ReadSet(cutsfile,cutvarnames,cutset);
  
  const Int_t nVars = axesno.size();
  const Int_t nPtBins = cutset.size()/(2*nVars);
  const Int_t nPtLims = nPtBins+1;
  Double_t *PtLims = (Double_t*)calloc(nPtBins,sizeof(Double_t)+1);
 
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    PtLims[iPt] = cutset[iPt*2*nVars];
  }
  PtLims[nPtBins] = cutset[2*nVars*(nPtBins-1)+1];

  return PtLims;
}

Double_t PtWeightsFromFONLL5overLHC13d3(Double_t *x, Double_t *pars){
  // weight function from the ratio of the LHC13d3 MC
  // and FONLL calculations for pp data
  Double_t pt = x[0];
  Double_t weight = (pars[0]*pt)/TMath::Power(pars[2],(1+TMath::Power(pars[3],pt/pars[1])))+pars[4]*TMath::Exp(pars[5]+pars[6]*pt)+pars[7]*TMath::Exp(pars[8]*pt);
  
  return weight;
}

//_________________________________________________________________________
Double_t PtBWeightsFromFONLL5overLHC13d3(Double_t *x, Double_t *pars){
  // weight function from the ratio of the LHC13d3 MC
  // and FONLL calculations for pp data rescaled for A
  Double_t pt = x[0];
  Double_t weight = pars[0]*TMath::Gaus(pt,pars[1],pars[2],kTRUE)+pars[3]+pars[4]*TMath::Log(pars[5]+pars[6]*pt+pars[7]*pt*pt);
  
  return weight;
}

void PrintCuts(TString axesfile,TString cutsfile) {

  vector<string> axesnames;
  vector<int> axesno;
  vector<string> cutvarnames;
  vector<double> cutset;
  ReadAxesNum(axesfile,axesnames,axesno);
  ReadSet(cutsfile,cutvarnames,cutset);

  const Int_t nVars = axesno.size();
  const Int_t nPtBins = cutset.size()/(2*nVars);
  
  cout <<"\n\n********************************** TOPOLOGICAL CUTS **********************************\n" << endl;
  for(Int_t iVar=1; iVar<nVars; iVar++) {
    cout << "* " << axesnames[iVar] << "\npT (GeV/c)    ";
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      cout <<" | " <<cutset[iPt*2*nVars] <<  "-" << cutset[iPt*2*nVars+1] << " |" << setw(5);
    }
    cout << endl;
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      cout << setw(5)<<cutset[iPt*2*nVars+(2*iVar)] <<  "-" << cutset[iPt*2*nVars+(2*iVar+1)] << "  ";
    }
    cout << "\n_______________________________"<<endl;
  }
  cout <<"\n**************************************************************************************\n\n" << endl;
}

void SetSideBandsRegion(THnSparseF* sparse, Int_t side, Int_t iPt) {

  Double_t approxsigma = 0.008+iPt*0.004/3;
  Double_t leftlimit;
  Double_t rightlimit;
  if(side==kLeft) {
    leftlimit = 1.69;
    rightlimit = 1.869-4*approxsigma;
  }
  else {
    leftlimit = 1.869+4*approxsigma;
    rightlimit = 2.05;
  }
    
  TAxis* massax = (TAxis*)sparse->GetAxis(0);
  Double_t binmin = massax->FindBin(leftlimit*1.0001);
  Double_t binmax = massax->FindBin(rightlimit*0.999);
  massax->SetRange(binmin,binmax);
  
}
