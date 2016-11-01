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

//_____________________________________________________________________________________________
//GLOBAL VARIABLES
//PtBins of the analysis
const Int_t nPtBins = 7;
const Int_t nPtLims = nPtBins+1;
const Double_t PtLims[nPtLims] = {2,3,4,5,6,8,12,16};
//range limits
const Double_t d0limFD[nPtBins] = {240,260,280,300,350,400,450};
const Double_t d0lim[nPtBins];
const Double_t d0limPrompt[nPtBins]={200,200,200,200,200,180,180};
//initial parameters
const Double_t initparprompt[nPtBins] = {0.9,0.,35,300,1};
const Double_t initparFD[nPtBins] = {0.5,0.,100,10,1};
const Double_t initparBkg[nPtBins] = {0.5,-30,50,100,50};
//colors
const Int_t colors[] = {kRed,kBlue,kGreen+3,kBlack,kMagenta,kCyan,kYellow-3,kOrange+7,kCyan+3,kGreen};
//line
TLine *line= new TLine(PtLims[0],1,PtLims[nPtBins],1);

//_____________________________________________________________________________________________
//functions for the systematic evaluation

void RangeSystematics();
void SigmaSystematics();
void MassRangeSystematics();
void SoverTSystematics();
void SideBandsSystematics();
void SideBandsDistanceSystematics();
void PrefitSystematics();
void BinningSystematics();
void PtReweightSystematics();

//_____________________________________________________________________________________________

void RangeSystematics() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);
  //___________________________________________________________________________________________
  //input files

  TFile infileMC("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/MC/AnalysisResultspPbMC.root","READ");
  TDirectoryFile* dirMC=(TDirectoryFile*)infileMC.Get("PWG3_D2H_InvMassDplus");
  TList* listMC=(TList*)dirMC->Get("coutputDplus_ImpParpPbMC0100");
  THnSparseF* hMassPtImpPrompt=(THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
  THnSparseF* hMassPtImpRecoFD=(THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
  THnSparseF* hMassPtImpTrueFD=(THnSparseF*)listMC->FindObject("hMassPtImpParTrueBfeed");
  infileMC.Close();

  TFile infileData("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root","READ");
  TDirectoryFile* dirData=(TDirectoryFile*)infileData.Get("PWG3_D2H_InvMassDplus");
  TList* listData=(TList*)dirData->Get("coutputDplus_ImpParpPbData0100");
  THnSparseF* hMassPtImpParAll=(THnSparseF*)listData->FindObject("hMassPtImpParAll");
  TNtuple* dataTree=(TNtuple*)dirData->Get("fNtupleDplus");

  //____________________________________________________________________________________________
  //analysis
  Int_t nRanges = 9;
  TH1F** hFrac = new TH1F*[nRanges];
  TH1F** hRatio = new TH1F*[nRanges];
  TLegend *l = new TLegend(0.25,0.55,0.8,0.84);
  l->SetTextSize(0.045);
  TLegend *l2 = new TLegend(0.25,0.55,0.8,0.84);
  l2->SetTextSize(0.045);

  for(Int_t iRange=0; iRange<nRanges; iRange++) {
    hFrac[iRange] = new TH1F(Form("hFrac%d",iRange),"",nPtBins,PtLims);
    hRatio[iRange] = new TH1F(Form("hRatio%d",iRange),"",nPtBins,PtLims);

    hFrac[iRange]->SetLineColor(colors[iRange+1]);
    hRatio[iRange]->SetLineColor(colors[iRange+1]);
    hFrac[iRange]->SetMarkerColor(colors[iRange+1]);
    hRatio[iRange]->SetMarkerColor(colors[iRange+1]);
    hFrac[iRange]->SetLineWidth(2);
    hRatio[iRange]->SetLineWidth(2);
    hFrac[iRange]->SetMarkerSize(1.5);
    hRatio[iRange]->SetMarkerSize(1.5);
    hFrac[iRange]->SetMarkerStyle(20+iRange);
    hRatio[iRange]->SetMarkerStyle(20+iRange);
    hFrac[iRange]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hFrac[iRange]->GetYaxis()->SetTitle("#it{f}_{prompt}");
    hFrac[iRange]->GetYaxis()->SetRangeUser(0.6,1.5);
    hRatio[iRange]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hRatio[iRange]->GetYaxis()->SetTitle("ratio of #it{f}_{prompt} w.r.t. the central value ");
    hRatio[iRange]->GetYaxis()->SetRangeUser(0.6,1.5);
    hRatio[iRange]->SetStats(0);
    hFrac[iRange]->SetStats(0);
    hFrac[iRange]->SetDirectory(0);
    hRatio[iRange]->SetDirectory(0);
    Double_t limit = 200+100*iRange;

    l->AddEntry(hFrac[iRange],Form("Imp Par Range [-%0.f,%0.f] #mum",limit,limit),"l");
    l2->AddEntry(hRatio[iRange],Form("Imp Par Range [-%0.f,%0.f] #mum",limit,limit),"p");
  }

  AliDplusCharmFractionIPfitter *ImpParFitter = new AliDplusCharmFractionIPfitter();
  ImpParFitter->SetDataSparse(hMassPtImpParAll);
  ImpParFitter->SetMCPromptSparse(hMassPtImpPrompt);
  ImpParFitter->SetMCTrueFDSparse(hMassPtImpTrueFD);
  ImpParFitter->SetMCRecoFDSparse(hMassPtImpRecoFD);
  ImpParFitter->SetDataTree(dataTree);

  ImpParFitter->SetSideBandsRegion(AliDplusCharmFractionIPfitter::kBoth);
  ImpParFitter->SetFDFunction(AliDplusCharmFractionIPfitter::kConvolution);
  ImpParFitter->SetBkgFunction(AliDplusCharmFractionIPfitter::kDoubleGaussExpoSymm);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ImpParFitter->SetRebinImpParHistos(2);
    ImpParFitter->SetPtLims(PtLims[iPt],PtLims[iPt+1]);
    ImpParFitter->SetNSigmas(2);
    ImpParFitter->FixSigmaPromptFromMC(kFALSE);
    ImpParFitter->SetRelativeLimitSigmaPrompt(0.2);
    ImpParFitter->SetLimitsPromptFrac(0.2,1.5);
    ImpParFitter->SetInitialParameters(initparprompt,initparFD,initparBkg);
    ImpParFitter->SetMassRegion(4,0,0,AliDplusCharmFractionIPfitter::kCentralValue);
    ImpParFitter->SetNSigmaSBLimits(3,15);

    for(Int_t iRange=1; iRange<=nRanges; iRange++) {
      Double_t limit = iRange*100+100;
      ImpParFitter->PrefitStep(-limit,limit,-limit,limit,-limit,limit);
      ImpParFitter->FitTree(-limit,limit,kFALSE);

      Double_t promptfrac = ImpParFitter->GetPromptFraction();
      Double_t promptfracerr = ImpParFitter->GetPromptFractionErr();

      if(iPt>3 && iPt<6) {
        if(iRange>2) {
          hFrac[iRange-1]->SetBinContent(iPt+1,promptfrac);
          hFrac[iRange-1]->SetBinError(iPt+1,promptfracerr);
        }
      }
      else if(iPt==6) {
        if(iRange>3) {
          hFrac[iRange-1]->SetBinContent(iPt+1,promptfrac);
          hFrac[iRange-1]->SetBinError(iPt+1,promptfracerr);
        }
      }
      else {
        hFrac[iRange-1]->SetBinContent(iPt+1,promptfrac);
        hFrac[iRange-1]->SetBinError(iPt+1,promptfracerr);
      }
    }
  }

  TFile resultfile("fraction_unbinned_sigmafree.root","UPDATE");
  TH1F* hPromptFrac=(TH1F*)resultfile.Get("hFrac");
  hPromptFrac->SetDirectory(0);

  for(Int_t iRange=0; iRange<nRanges; iRange++) {
    hRatio[iRange]->Divide(hFrac[iRange],hPromptFrac,1.,1.);
    for(Int_t iPt=0; iPt<hRatio[iRange]->GetNbinsX(); iPt++) {
      hRatio[iRange]->SetBinError(iPt+1,1.e-10);
    }
  }

  TCanvas *cFrac = new TCanvas("cFrac","",800,800);
  hFrac[0]->Draw("E1");
  l->Draw("same");

  TCanvas *cRatio = new TCanvas("cRatio","",800,800);
  hRatio[0]->Draw("E1");
  l2->Draw("same");
  line->SetLineColor(colors[0]);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw("same");

  for(Int_t iRange=1; iRange<nRanges; iRange++) {
    cFrac->cd();
    hFrac[iRange]->Draw("E1same");
    cRatio->cd();
    hRatio[iRange]->Draw("E1same");
  }

  cFrac->SaveAs("range/promptfraction_syst_range_onlyfrac.eps");
  cRatio->SaveAs("range/promptfraction_syst_range_onlyratio.eps");

  TFile outfile("range/promptfraction_syst_range.root","RECREATE");
  for(Int_t iRange=0; iRange<nRanges; iRange++) {
    hFrac[iRange]->Write();
    hRatio[iRange]->Write();
  }
  outfile.Close();
}

void SigmaSystematics() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);

  //___________________________________________________________________________________________
  //range limits
  for(Int_t iBin=0; iBin<nPtBins; iBin++) {
    d0lim[iBin] = d0limFD[iBin];
  }
  //___________________________________________________________________________________________
  //input files

  TFile infileMC("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/MC/AnalysisResultspPbMC.root","READ");
  TDirectoryFile* dirMC=(TDirectoryFile*)infileMC.Get("PWG3_D2H_InvMassDplus");
  TList* listMC=(TList*)dirMC->Get("coutputDplus_ImpParpPbMC0100");
  THnSparseF* hMassPtImpPrompt=(THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
  THnSparseF* hMassPtImpRecoFD=(THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
  THnSparseF* hMassPtImpTrueFD=(THnSparseF*)listMC->FindObject("hMassPtImpParTrueBfeed");
  infileMC.Close();

  TFile infileData("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root","READ");
  TDirectoryFile* dirData=(TDirectoryFile*)infileData.Get("PWG3_D2H_InvMassDplus");
  TList* listData=(TList*)dirData->Get("coutputDplus_ImpParpPbData0100");
  THnSparseF* hMassPtImpParAll=(THnSparseF*)listData->FindObject("hMassPtImpParAll");
  TNtuple* dataTree=(TNtuple*)dirData->Get("fNtupleDplus");

  //infileData.Close();

  //____________________________________________________________________________________________
  //analysis

  enum sigma {kFree, kFixed};

  TH1F** hFrac = new TH1F*[2];
  TH1F** hSigma = new TH1F*[2];("hSigma","",nPtBins,PtLims);
  hFrac[kFree] = new TH1F("hFracFree","",nPtBins,PtLims);
  hFrac[kFixed] = new TH1F("hFracFixed","",nPtBins,PtLims);
  hSigma[kFree] = new TH1F("hSigmaFree","",nPtBins,PtLims);
  hSigma[kFixed] = new TH1F("hSigmaFixed","",nPtBins,PtLims);

  for(Int_t iSigma=kFree; iSigma<=kFixed; iSigma++) {
    hFrac[iSigma]->SetDirectory(0);
    hSigma[iSigma]->SetDirectory(0);
    hFrac[iSigma]->SetLineColor(colors[iSigma]);
    hSigma[iSigma]->SetLineColor(colors[iSigma]);
    hFrac[iSigma]->SetLineWidth(2);
    hSigma[iSigma]->SetLineWidth(2);
    hFrac[iSigma]->SetMarkerColor(colors[iSigma]);
    hSigma[iSigma]->SetMarkerColor(colors[iSigma]);
    hFrac[iSigma]->SetMarkerSize(1.5);
    hSigma[iSigma]->SetMarkerSize(1.5);
    hFrac[iSigma]->SetMarkerStyle(20+iSigma);
    hSigma[iSigma]->SetMarkerStyle(20+iSigma);
    hFrac[iSigma]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hFrac[iSigma]->GetYaxis()->SetTitle("#it{f}_{prompt}");
    hSigma[iSigma]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hSigma[iSigma]->GetYaxis()->SetTitle("#sigma_{prompt} (#mum)");
  }

  AliDplusCharmFractionIPfitter *ImpParFitter = new AliDplusCharmFractionIPfitter();
  ImpParFitter->SetDataSparse(hMassPtImpParAll);
  ImpParFitter->SetMCPromptSparse(hMassPtImpPrompt);
  ImpParFitter->SetMCTrueFDSparse(hMassPtImpTrueFD);
  ImpParFitter->SetMCRecoFDSparse(hMassPtImpRecoFD);
  ImpParFitter->SetDataTree(dataTree);

  ImpParFitter->SetSideBandsRegion(AliDplusCharmFractionIPfitter::kBoth);
  ImpParFitter->SetFDFunction(AliDplusCharmFractionIPfitter::kConvolution);
  ImpParFitter->SetBkgFunction(AliDplusCharmFractionIPfitter::kDoubleGaussExpoSymm);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ImpParFitter->SetRebinImpParHistos(2);
    ImpParFitter->SetPtLims(PtLims[iPt],PtLims[iPt+1]);
    ImpParFitter->SetNSigmas(2);
    ImpParFitter->GetSignal(4,0,0,1.68,2.05,AliDplusCharmFractionIPfitter::kCentralValue);
    ImpParFitter->PrefitStep(-d0limPrompt[iPt],d0limPrompt[iPt],-d0limFD[iPt],d0limFD[iPt],-1000,1000);
    ImpParFitter->SetRelativeLimitSigmaPrompt(0.2);
    ImpParFitter->SetLimitsPromptFrac(0.2,1.5);
    ImpParFitter->SetInitialParameters(initparprompt,initparFD,initparBkg);

    for(Int_t iSigma=kFree; iSigma<=kFixed; iSigma++) {
      if(iSigma==kFree)
        ImpParFitter->FixSigmaPromptFromMC(kFALSE);
      else
        ImpParFitter->FixSigmaPromptFromMC(kTRUE);

      ImpParFitter->FitTree(-d0limFD[iPt],d0limFD[iPt],kFALSE);

      Double_t promptfrac = ImpParFitter->GetPromptFraction();
      Double_t promptfracerr = ImpParFitter->GetPromptFractionErr();
      Double_t promptsigma = ImpParFitter->GetPromptSigma();
      Double_t promptsigmaerr = ImpParFitter->GetPromptSigmaErr();
      if(iSigma==kFixed) {
        promptsigma = ImpParFitter->GetPromptSigmaMC();
        promptsigmaerr = ImpParFitter->GetPromptSigmaMCErr();
      }
      hFrac[iSigma]->SetBinContent(iPt+1,promptfrac);
      hFrac[iSigma]->SetBinError(iPt+1,promptfracerr);
      hSigma[iSigma]->SetBinContent(iPt+1,promptsigma);
      hSigma[iSigma]->SetBinError(iPt+1,promptsigmaerr);
    }
  }

  TH1F* hRatio = new TH1F("hRatio","",nPtBins,PtLims);
  hRatio->Divide(hFrac[kFixed],hFrac[kFree],1.,1.);
  hRatio->SetDirectory(0);
  hRatio->SetLineColor(colors[kFixed]);
  hRatio->SetMarkerColor(colors[kFixed]);
  hRatio->SetLineWidth(2);
  hRatio->SetMarkerSize(1.5);
  hRatio->SetMarkerStyle(20);
  hRatio->GetYaxis()->SetTitle("#it{f}_{prompt}^{#sigma fixed}/#it{f}_{prompt}^{#sigma free}");
  hRatio->GetXaxis()->SetTitle(hFrac[0]->GetXaxis()->GetTitle());

  for(Int_t iPt=0; iPt<hRatio->GetNbinsX(); iPt++) {
    hRatio->SetBinError(iPt+1,1.e-10);
  }

  TLegend *lFrac = new TLegend(0.4,0.2,0.8,0.49);
  TLegend *lSigma = new TLegend(0.5,0.6,0.89,0.89);
  lFrac->AddEntry(hFrac[kFree], "#sigma_{prompt} free", "lpe");
  lFrac->AddEntry(hFrac[kFixed], "#sigma_{prompt} fixed", "lpe");
  lSigma->AddEntry(hSigma[kFree], "#sigma_{prompt} from fit on data", "lpe");
  lSigma->AddEntry(hSigma[kFixed], "#sigma_{prompt} from MC prefit", "lpe");

  TCanvas* cFrac = new TCanvas("cFrac","",800,800);
  hFrac[kFree]->GetYaxis()->SetRangeUser(0,1.4);
  hFrac[kFree]->Draw("E1");
  hFrac[kFixed]->Draw("E1same");
  lFrac->Draw("same");

  TCanvas* cRatio = new TCanvas("cRatio","",800,800);
  hRatio->Draw("E1");
  line->SetLineColor(colors[0]);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw("same");

  TCanvas* cSigma = new TCanvas("cSigma","",800,600);
  cSigma->Clear();
  hSigma[kFree]->GetYaxis()->SetRangeUser(20,60);
  hSigma[kFree]->Draw("E1");
  hSigma[kFixed]->Draw("E1same");
  lSigma->Draw("same");

  TFile outfile("sigma/promptfraction_syst_sigma.root","RECREATE");
  for(Int_t iSigma=kFree; iSigma<=kFixed; iSigma++) {
    hFrac[iSigma]->Write();
    hSigma[iSigma]->Write();
  }
  hRatio->Write();
  outfile.Close();

  cFrac->SaveAs("sigma/promptfraction_syst_sigma_onlyfrac.eps");
  cFrac->SaveAs("sigma/promptfraction_syst_sigma_onlyratio.eps");
  cSigma->SaveAs("sigma/sigmaprompt.eps");
}

void MassRangeSystematics() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);

  //___________________________________________________________________________________________
  //range limits
  for(Int_t iBin=0; iBin<nPtBins; iBin++) {
    d0lim[iBin] = d0limFD[iBin];
  }
  //___________________________________________________________________________________________
  //input files

  TFile infileMC("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/MC/AnalysisResultspPbMC.root","READ");
  TDirectoryFile* dirMC=(TDirectoryFile*)infileMC.Get("PWG3_D2H_InvMassDplus");
  TList* listMC=(TList*)dirMC->Get("coutputDplus_ImpParpPbMC0100");
  THnSparseF* hMassPtImpPrompt=(THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
  THnSparseF* hMassPtImpRecoFD=(THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
  THnSparseF* hMassPtImpTrueFD=(THnSparseF*)listMC->FindObject("hMassPtImpParTrueBfeed");
  infileMC.Close();

  TFile infileData("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root","READ");
  TDirectoryFile* dirData=(TDirectoryFile*)infileData.Get("PWG3_D2H_InvMassDplus");
  TList* listData=(TList*)dirData->Get("coutputDplus_ImpParpPbData0100");
  THnSparseF* hMassPtImpParAll=(THnSparseF*)listData->FindObject("hMassPtImpParAll");
  TNtuple* dataTree=(TNtuple*)dirData->Get("fNtupleDplus");

  //infileData.Close();
  //____________________________________________________________________________________________
  //analysis

  const Int_t nSigma = 4;
  Double_t nSigmas[nSigma] = {1.,1.5,2.5,3.};

  TH1F** hFrac = new TH1F*[nSigma];
  TH1F** hRatio = new TH1F*[nSigma];
  for(Int_t iSigma=0; iSigma<nSigma; iSigma++) {
    hFrac[iSigma] = new TH1F(Form("hFrac%d",iSigma),"",nPtBins,PtLims);
    hFrac[iSigma]->SetDirectory(0);
    hFrac[iSigma]->SetLineColor(colors[iSigma+1]);
    hFrac[iSigma]->SetMarkerColor(colors[iSigma+1]);
    hFrac[iSigma]->SetLineWidth(2);
    hFrac[iSigma]->SetMarkerSize(1.5);
    hFrac[iSigma]->SetMarkerStyle(20+iSigma);
    hRatio[iSigma] = new TH1F(Form("hRatio%d",iSigma),"",nPtBins,PtLims);
    hRatio[iSigma]->SetDirectory(0);
    hRatio[iSigma]->SetLineColor(colors[iSigma+1]);
    hRatio[iSigma]->SetMarkerColor(colors[iSigma+1]);
    hRatio[iSigma]->SetLineWidth(2);
    hRatio[iSigma]->SetMarkerSize(1.5);
    hRatio[iSigma]->SetMarkerStyle(20+iSigma);
    hFrac[iSigma]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hFrac[iSigma]->GetYaxis()->SetTitle("#it{f}_{prompt}");
    hRatio[iSigma]->GetYaxis()->SetTitle("ratio of #it{f}_{prompt} w.r.t. the central value");
    hRatio[iSigma]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  }

  AliDplusCharmFractionIPfitter *ImpParFitter = new AliDplusCharmFractionIPfitter();
  ImpParFitter->SetDataSparse(hMassPtImpParAll);
  ImpParFitter->SetMCPromptSparse(hMassPtImpPrompt);
  ImpParFitter->SetMCTrueFDSparse(hMassPtImpTrueFD);
  ImpParFitter->SetMCRecoFDSparse(hMassPtImpRecoFD);
  ImpParFitter->SetDataTree(dataTree);

  ImpParFitter->SetSideBandsRegion(AliDplusCharmFractionIPfitter::kBoth);
  ImpParFitter->SetFDFunction(AliDplusCharmFractionIPfitter::kConvolution);
  ImpParFitter->SetBkgFunction(AliDplusCharmFractionIPfitter::kDoubleGaussExpoSymm);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ImpParFitter->SetPtLims(PtLims[iPt],PtLims[iPt+1]);

    for(Int_t iSigma=0; iSigma<nSigma; iSigma++) {
    ImpParFitter->SetRelativeLimitSigmaPrompt(0.2);
    ImpParFitter->SetLimitsPromptFrac(0.2,1.5);
    ImpParFitter->SetInitialParameters(initparprompt,initparFD,initparBkg);
    ImpParFitter->FixSigmaPromptFromMC(kFALSE);
    ImpParFitter->SetRebinImpParHistos(2);
    ImpParFitter->SetNSigmas(nSigmas[iSigma]);
    ImpParFitter->GetSignal(4,0,0,1.68,2.05,AliDplusCharmFractionIPfitter::kCentralValue);
    ImpParFitter->PrefitStep(-d0limPrompt[iPt],d0limPrompt[iPt],-d0limFD[iPt],d0limFD[iPt],-1000,1000);
    ImpParFitter->FitTree(-d0limFD[iPt],d0limFD[iPt],kFALSE);

      Double_t promptfrac = ImpParFitter->GetPromptFraction();
      Double_t promptfracerr = ImpParFitter->GetPromptFractionErr();
      hFrac[iSigma]->SetBinContent(iPt+1,promptfrac);
      hFrac[iSigma]->SetBinError(iPt+1,promptfracerr);
    }
  }

  TFile resultfile("fraction_unbinned_sigmafree.root","UPDATE");
  TH1F* hPromptFrac=(TH1F*)resultfile.Get("hFrac");
  hPromptFrac->SetDirectory(0);
  resultfile.Close();

  TLegend *l = new TLegend(0.5,0.2,0.89,0.4);
  l->SetTextSize(0.045);
  TLegend *l2 = new TLegend(0.5,0.2,0.89,0.4);
  l2->SetTextSize(0.045);

  for(Int_t iSigma=0; iSigma<nSigma; iSigma++) {
    hRatio[iSigma]->Divide(hFrac[iSigma],hPromptFrac,1.,1.);
    for(Int_t iPt=0; iPt<hRatio[0]->GetNbinsX(); iPt++) {
      hRatio[iSigma]->SetBinError(iPt+1,1.e-10);
    }
    l->AddEntry(hFrac[iSigma], Form("n#sigma = %0.1f",nSigmas[iSigma]), "lpe");
    l2->AddEntry(hRatio[iSigma], Form("n#sigma = %0.1f",nSigmas[iSigma]), "lpe");
  }

  TCanvas* cFrac = new TCanvas("cFrac","",800,800);
  hFrac[0]->Draw("E1");
  hFrac[0]->GetYaxis()->SetRangeUser(0,1.4);
  for(Int_t iSigma=1; iSigma<nSigma; iSigma++) {
    hFrac[iSigma]->Draw("E1same");
  }
  l->Draw("same");

  TCanvas* cRatio = new TCanvas("cRatio","",800,800);
  hRatio[0]->Draw("E1");
  hRatio[0]->GetYaxis()->SetRangeUser(0.75,1.15);
  for(Int_t iSigma=1; iSigma<nSigma; iSigma++) {
    hRatio[iSigma]->Draw("E1same");
  }
  l2->Draw("same");
  line->SetLineColor(colors[0]);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw("same");

  cFrac->SaveAs("massrange/promptfraction_syst_MassRange_onlyfrac.eps");
  cRatio->SaveAs("massrange/promptfraction_syst_MassRange_onlyratio.eps");

  TFile outfile("massrange/promptfraction_syst_MassRange.root","RECREATE");
  for(Int_t iSigma=0; iSigma<nSigma; iSigma++) {
    hFrac[iSigma]->Write();
    hRatio[iSigma]->Write();
  }
  outfile.Close();
}

void SoverTSystematics() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);

  //___________________________________________________________________________________________
  //range limits
  for(Int_t iBin=0; iBin<nPtBins; iBin++) {
    d0lim[iBin] = d0limFD[iBin];
  }
  //___________________________________________________________________________________________
  //input files

  TFile infileMC("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/MC/AnalysisResultspPbMC.root","READ");
  TDirectoryFile* dirMC=(TDirectoryFile*)infileMC.Get("PWG3_D2H_InvMassDplus");
  TList* listMC=(TList*)dirMC->Get("coutputDplus_ImpParpPbMC0100");
  THnSparseF* hMassPtImpPrompt=(THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
  THnSparseF* hMassPtImpRecoFD=(THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
  THnSparseF* hMassPtImpTrueFD=(THnSparseF*)listMC->FindObject("hMassPtImpParTrueBfeed");
  infileMC.Close();

  TFile infileData("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root","READ");
  TDirectoryFile* dirData=(TDirectoryFile*)infileData.Get("PWG3_D2H_InvMassDplus");
  TList* listData=(TList*)dirData->Get("coutputDplus_ImpParpPbData0100");
  THnSparseF* hMassPtImpParAll=(THnSparseF*)listData->FindObject("hMassPtImpParAll");
  TNtuple* dataTree=(TNtuple*)dirData->Get("fNtupleDplus");

  //infileData.Close();
  //____________________________________________________________________________________________
  //analysis
  const Int_t nSoverT = 6;

  Int_t sovert[nSoverT] = {AliDplusCharmFractionIPfitter::kUpperLimit,
                           AliDplusCharmFractionIPfitter::kUpperLimitOnlyStat,
                           AliDplusCharmFractionIPfitter::kUpperLimitOnlySyst,
                           AliDplusCharmFractionIPfitter::kLowerLimit,
                           AliDplusCharmFractionIPfitter::kLowerLimitOnlyStat,
                           AliDplusCharmFractionIPfitter::kLowerLimitOnlySyst};

  TH1F** hFrac = new TH1F*[nSoverT];
  TH1F** hRatio = new TH1F*[nSoverT];
  for(Int_t iSoverT=0; iSoverT<nSoverT; iSoverT++) {
    hFrac[iSoverT] = new TH1F(Form("hFrac%d",iSoverT),"",nPtBins,PtLims);
    hFrac[iSoverT]->SetDirectory(0);
    hFrac[iSoverT]->SetLineColor(colors[iSoverT+1]);
    hFrac[iSoverT]->SetLineWidth(2);
    hFrac[iSoverT]->SetMarkerStyle(20+iSoverT);
    hFrac[iSoverT]->SetMarkerSize(1.5);
    hFrac[iSoverT]->SetMarkerColor(colors[iSoverT+1]);
    hRatio[iSoverT] = new TH1F(Form("hRatio%d",iSoverT),"",nPtBins,PtLims);
    hRatio[iSoverT]->SetDirectory(0);
    hRatio[iSoverT]->SetLineColor(colors[iSoverT+1]);
    hRatio[iSoverT]->SetLineWidth(2);
    hRatio[iSoverT]->SetMarkerStyle(20+iSoverT);
    hRatio[iSoverT]->SetMarkerSize(1.5);
    hRatio[iSoverT]->SetMarkerColor(colors[iSoverT+1]);
    hFrac[iSoverT]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hFrac[iSoverT]->GetYaxis()->SetTitle("#it{f}_{prompt}");
    hRatio[iSoverT]->GetYaxis()->SetTitle("ratio of #it{f}_{prompt} w.r.t. the central value");
    hRatio[iSoverT]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  }

  AliDplusCharmFractionIPfitter *ImpParFitter = new AliDplusCharmFractionIPfitter();
  ImpParFitter->SetDataSparse(hMassPtImpParAll);
  ImpParFitter->SetMCPromptSparse(hMassPtImpPrompt);
  ImpParFitter->SetMCTrueFDSparse(hMassPtImpTrueFD);
  ImpParFitter->SetMCRecoFDSparse(hMassPtImpRecoFD);
  ImpParFitter->SetDataTree(dataTree);

  ImpParFitter->SetSideBandsRegion(AliDplusCharmFractionIPfitter::kBoth);
  ImpParFitter->SetFDFunction(AliDplusCharmFractionIPfitter::kConvolution);
  ImpParFitter->SetBkgFunction(AliDplusCharmFractionIPfitter::kDoubleGaussExpoSymm);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ImpParFitter->SetRelativeLimitSigmaPrompt(0.2);
    ImpParFitter->SetLimitsPromptFrac(0.2,1.5);
    ImpParFitter->SetInitialParameters(initparprompt,initparFD,initparBkg);
    ImpParFitter->FixSigmaPromptFromMC(kFALSE);
    ImpParFitter->SetRebinImpParHistos(2);
    ImpParFitter->SetPtLims(PtLims[iPt],PtLims[iPt+1]);
    ImpParFitter->SetNSigmas(2);

    for(Int_t iSoverT=0; iSoverT<nSoverT; iSoverT++) {
      ImpParFitter->GetSignal(4,0,0,1.68,2.05,sovert[iSoverT]);
      ImpParFitter->PrefitStep(-d0limPrompt[iPt],d0limPrompt[iPt],-d0limFD[iPt],d0limFD[iPt],-1000,1000);
      ImpParFitter->FitTree(-d0limFD[iPt],d0limFD[iPt],kFALSE);

      Double_t promptfrac = ImpParFitter->GetPromptFraction();
      Double_t promptfracerr = ImpParFitter->GetPromptFractionErr();
      hFrac[iSoverT]->SetBinContent(iPt+1,promptfrac);
      hFrac[iSoverT]->SetBinError(iPt+1,promptfracerr);
    }
  }

  TFile resultfile("fraction_unbinned_sigmafree.root","UPDATE");
  TH1F* hPromptFrac=(TH1F*)resultfile.Get("hFrac");
  hPromptFrac->SetDirectory(0);
  resultfile.Close();

  TLegend *l = new TLegend(0.4,0.15,0.89,0.4);
  l->SetTextSize(0.045);

  for(Int_t iSoverT=0; iSoverT<nSoverT; iSoverT++) {
    hRatio[iSoverT]->Divide(hFrac[iSoverT],hPromptFrac,1.,1.);
    for(Int_t iPt=0; iPt<hRatio[0]->GetNbinsX(); iPt++) {
      hRatio[iSoverT]->SetBinError(iPt+1,1.e-10);
    }
  }
  l->AddEntry(hFrac[0], "S + #sigma_{S}", "l");
  l->AddEntry(hFrac[1], "S + #sigma_{S}(stat)", "l");
  l->AddEntry(hFrac[2], "S + #sigma_{S}(syst)", "l");
  l->AddEntry(hFrac[3], "S - #sigma_{S}", "l");
  l->AddEntry(hFrac[4], "S - #sigma_{S}(stat)", "l");
  l->AddEntry(hFrac[5], "S - #sigma_{S}(syst)", "l");

  TCanvas* cFrac = new TCanvas("cFrac","",800,800);
  hFrac[0]->Draw("E1");
  hFrac[0]->GetYaxis()->SetRangeUser(0,1.4);
  for(Int_t iSoverT=1; iSoverT<nSoverT; iSoverT++) {
    hFrac[iSoverT]->Draw("E1same");
  }
  l->Draw("same");

  TCanvas* cRatio = new TCanvas("cRatio","",800,800);
  hRatio[0]->Draw("E1");
  hRatio[0]->GetYaxis()->SetRangeUser(0.6,1.2);
  for(Int_t iSoverT=1; iSoverT<nSoverT; iSoverT++) {
    hRatio[iSoverT]->Draw("E1same");
  }
  l->Draw("same");
  line->SetLineColor(colors[0]);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw("same");

  cFrac->SaveAs("SoverT/promptfraction_syst_SoverT_onlyfrac.eps");
  cRatio->SaveAs("SoverT/promptfraction_syst_SoverT_onlyratio.eps");

  TFile outfile("SoverT/promptfraction_syst_SoverT.root","RECREATE");
  for(Int_t iSoverT=0; iSoverT<nSoverT; iSoverT++) {
    hFrac[iSoverT]->Write();
    hRatio[iSoverT]->Write();
  }
  outfile.Close();

}

void SideBandsSystematics() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);

  //___________________________________________________________________________________________
  //range limits
  for(Int_t iBin=0; iBin<nPtBins; iBin++) {
    d0lim[iBin] = d0limFD[iBin];
  }
  //___________________________________________________________________________________________
  //input files

  TFile infileMC("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/MC/AnalysisResultspPbMC.root","READ");
  TDirectoryFile* dirMC=(TDirectoryFile*)infileMC.Get("PWG3_D2H_InvMassDplus");
  TList* listMC=(TList*)dirMC->Get("coutputDplus_ImpParpPbMC0100");
  THnSparseF* hMassPtImpPrompt=(THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
  THnSparseF* hMassPtImpRecoFD=(THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
  THnSparseF* hMassPtImpTrueFD=(THnSparseF*)listMC->FindObject("hMassPtImpParTrueBfeed");
  infileMC.Close();

  TFile infileData("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root","READ");
  TDirectoryFile* dirData=(TDirectoryFile*)infileData.Get("PWG3_D2H_InvMassDplus");
  TList* listData=(TList*)dirData->Get("coutputDplus_ImpParpPbData0100");
  THnSparseF* hMassPtImpParAll=(THnSparseF*)listData->FindObject("hMassPtImpParAll");
  TNtuple* dataTree=(TNtuple*)dirData->Get("fNtupleDplus");

  //infileData.Close();
  //____________________________________________________________________________________________
  //analysis

  Int_t SB[2] = {AliDplusCharmFractionIPfitter::kLeft, AliDplusCharmFractionIPfitter::kRight};

  TH1F** hFrac = new TH1F*[2];
  TH1F** hRatio = new TH1F*[2];
  for(Int_t iSB=0; iSB<2; iSB++) {
    hFrac[iSB] = new TH1F(Form("hFrac%d",iSB),"",nPtBins,PtLims);
    hFrac[iSB]->SetDirectory(0);
    hFrac[iSB]->SetLineColor(colors[iSB+1]);
    hFrac[iSB]->SetLineWidth(2);
    hFrac[iSB]->SetMarkerStyle(20+iSB);
    hFrac[iSB]->SetMarkerSize(1.5);
    hFrac[iSB]->SetMarkerColor(colors[iSB+1]);
    hRatio[iSB] = new TH1F(Form("hRatio%d",iSB),"",nPtBins,PtLims);
    hRatio[iSB]->SetDirectory(0);
    hRatio[iSB]->SetLineColor(colors[iSB+1]);
    hRatio[iSB]->SetLineWidth(2);
    hRatio[iSB]->SetMarkerStyle(20+iSB);
    hRatio[iSB]->SetMarkerSize(1.5);
    hRatio[iSB]->SetMarkerColor(colors[iSB+1]);
    hFrac[iSB]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hFrac[iSB]->GetYaxis()->SetTitle("#it{f}_{prompt}");
    hRatio[iSB]->GetYaxis()->SetTitle("ratio of #it{f}_{prompt} w.r.t. the central value");
    hRatio[iSB]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  }

  AliDplusCharmFractionIPfitter *ImpParFitter = new AliDplusCharmFractionIPfitter();
  ImpParFitter->SetDataSparse(hMassPtImpParAll);
  ImpParFitter->SetMCPromptSparse(hMassPtImpPrompt);
  ImpParFitter->SetMCTrueFDSparse(hMassPtImpTrueFD);
  ImpParFitter->SetMCRecoFDSparse(hMassPtImpRecoFD);
  ImpParFitter->SetDataTree(dataTree);
  ImpParFitter->SetFDFunction(AliDplusCharmFractionIPfitter::kConvolution);
  ImpParFitter->SetBkgFunction(AliDplusCharmFractionIPfitter::kDoubleGaussExpoSymm);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ImpParFitter->SetRebinImpParHistos(2);
    ImpParFitter->SetPtLims(PtLims[iPt],PtLims[iPt+1]);
    ImpParFitter->SetNSigmas(2.);

    for(Int_t iSB=0; iSB<2; iSB++) {
      ImpParFitter->SetRelativeLimitSigmaPrompt(0.2);
      ImpParFitter->SetLimitsPromptFrac(0.2,1.5);
      ImpParFitter->SetInitialParameters(initparprompt,initparFD,initparBkg);
      ImpParFitter->FixSigmaPromptFromMC(kFALSE);
      ImpParFitter->SetSideBandsRegion(SB[iSB]);
      ImpParFitter->GetSignal(4,0,0,1.68,2.05,AliDplusCharmFractionIPfitter::kCentralValue);
      ImpParFitter->PrefitStep(-d0limPrompt[iPt],d0limPrompt[iPt],-d0limFD[iPt],d0limFD[iPt],-1000,1000);
      ImpParFitter->FitTree(-d0limFD[iPt],d0limFD[iPt],kFALSE);

      Double_t promptfrac = ImpParFitter->GetPromptFraction();
      Double_t promptfracerr = ImpParFitter->GetPromptFractionErr();
      hFrac[iSB]->SetBinContent(iPt+1,promptfrac);
      hFrac[iSB]->SetBinError(iPt+1,promptfracerr);
    }
  }

  TFile resultfile("fraction_unbinned_sigmafree.root","UPDATE");
  TH1F* hPromptFrac=(TH1F*)resultfile.Get("hFrac");
  hPromptFrac->SetDirectory(0);

  TLegend *l = new TLegend(0.3,0.2,0.89,0.45);
  l->SetTextSize(0.045);

  for(Int_t iSB=0; iSB<2; iSB++) {
    hRatio[iSB]->Divide(hFrac[iSB],hPromptFrac,1.,1.);
    for(Int_t iPt=0; iPt<hRatio[0]->GetNbinsX(); iPt++) {
      hRatio[iSB]->SetBinError(iPt+1,1.e-10);
    }
  }
  l->AddEntry(hFrac[0], "left SB distribution", "lpe");
  l->AddEntry(hFrac[1], "right SB distribution", "lpe");

  TCanvas* cFrac = new TCanvas("cFrac","",1200,600);
  hFrac[0]->Draw("E1");
  hFrac[0]->GetYaxis()->SetRangeUser(0,1.5);
  hFrac[1]->Draw("E1same");
  l->Draw("same");

  TCanvas* cRatio = new TCanvas("cRatio","",1200,600);
  hRatio[0]->Draw("E1");
  hRatio[0]->GetYaxis()->SetRangeUser(0.6,1.2);
  hRatio[1]->Draw("E1same");
  l->Draw("same");
  line->SetLineColor(colors[0]);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw("same");

  cFrac->SaveAs("sidebands/promptfraction_syst_SB_onlyfrac.eps");
  cRatio->SaveAs("sidebands/promptfraction_syst_SB_onlyratio.eps");

  TFile outfile("sidebands/promptfraction_syst_SB.root","RECREATE");
  for(Int_t iSB=0; iSB<2; iSB++) {
    hFrac[iSB]->Write();
    hRatio[iSB]->Write();
  }
  outfile.Close();

}

void SideBandsDistanceSystematics() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);

  //___________________________________________________________________________________________
  //range limits
  for(Int_t iBin=0; iBin<nPtBins; iBin++) {
    d0lim[iBin] = d0limFD[iBin];
  }
  //___________________________________________________________________________________________
  //input files

  TFile infileMC("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/MC/AnalysisResultspPbMC.root","READ");
  TDirectoryFile* dirMC=(TDirectoryFile*)infileMC.Get("PWG3_D2H_InvMassDplus");
  TList* listMC=(TList*)dirMC->Get("coutputDplus_ImpParpPbMC0100");
  THnSparseF* hMassPtImpPrompt=(THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
  THnSparseF* hMassPtImpRecoFD=(THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
  THnSparseF* hMassPtImpTrueFD=(THnSparseF*)listMC->FindObject("hMassPtImpParTrueBfeed");
  infileMC.Close();

  TFile infileData("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root","READ");
  TDirectoryFile* dirData=(TDirectoryFile*)infileData.Get("PWG3_D2H_InvMassDplus");
  TList* listData=(TList*)dirData->Get("coutputDplus_ImpParpPbData0100");
  THnSparseF* hMassPtImpParAll=(THnSparseF*)listData->FindObject("hMassPtImpParAll");
  TNtuple* dataTree=(TNtuple*)dirData->Get("fNtupleDplus");

  //infileData.Close();
  //____________________________________________________________________________________________
  //analysis

  const Int_t nSB = 3;
  Int_t SB[nSB] = {4,5,6};
  Int_t Delta = 10;

  TH1F** hFrac = new TH1F*[nSB];
  TH1F** hRatio = new TH1F*[nSB];
  for(Int_t iSB=0; iSB<nSB; iSB++) {
    hFrac[iSB] = new TH1F(Form("hFrac%d",iSB),"",nPtBins,PtLims);
    hFrac[iSB]->SetDirectory(0);
    hFrac[iSB]->SetLineColor(colors[iSB+3]);
    hFrac[iSB]->SetLineWidth(2);
    hFrac[iSB]->SetMarkerStyle(20+iSB);
    hFrac[iSB]->SetMarkerSize(1.5);
    hFrac[iSB]->SetMarkerColor(colors[iSB+3]);
    hRatio[iSB] = new TH1F(Form("hRatio%d",iSB),"",nPtBins,PtLims);
    hRatio[iSB]->SetDirectory(0);
    hRatio[iSB]->SetLineColor(colors[iSB+3]);
    hRatio[iSB]->SetLineWidth(2);
    hRatio[iSB]->SetMarkerStyle(20+iSB);
    hRatio[iSB]->SetMarkerSize(1.5);
    hRatio[iSB]->SetMarkerColor(colors[iSB+3]);
    hFrac[iSB]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hFrac[iSB]->GetYaxis()->SetTitle("#it{f}_{prompt}");
    hRatio[iSB]->GetYaxis()->SetTitle("ratio of #it{f}_{prompt} w.r.t the central value");
    hRatio[iSB]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  }

  AliDplusCharmFractionIPfitter *ImpParFitter = new AliDplusCharmFractionIPfitter();
  ImpParFitter->SetDataSparse(hMassPtImpParAll);
  ImpParFitter->SetMCPromptSparse(hMassPtImpPrompt);
  ImpParFitter->SetMCTrueFDSparse(hMassPtImpTrueFD);
  ImpParFitter->SetMCRecoFDSparse(hMassPtImpRecoFD);
  ImpParFitter->SetDataTree(dataTree);

  ImpParFitter->SetFDFunction(AliDplusCharmFractionIPfitter::kConvolution);
  ImpParFitter->SetBkgFunction(AliDplusCharmFractionIPfitter::kDoubleGaussExpoSymm);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ImpParFitter->SetRebinImpParHistos(2);
    ImpParFitter->SetPtLims(PtLims[iPt],PtLims[iPt+1]);
    ImpParFitter->SetNSigmas(2.);

    for(Int_t iSB=0; iSB<nSB; iSB++) {
      ImpParFitter->SetRelativeLimitSigmaPrompt(0.2);
      ImpParFitter->SetLimitsPromptFrac(0.2,1.5);
      ImpParFitter->SetInitialParameters(initparprompt,initparFD,initparBkg);
      ImpParFitter->FixSigmaPromptFromMC(kFALSE);
      ImpParFitter->SetSideBandsRegion(AliDplusCharmFractionIPfitter::kBoth);
      ImpParFitter->GetSignal(4,0,0,1.68,2.05,AliDplusCharmFractionIPfitter::kCentralValue);
      ImpParFitter->SetNSigmaSBLimits(SB[iSB],SB[iSB]+Delta);
      ImpParFitter->PrefitStep(-d0limPrompt[iPt],d0limPrompt[iPt],-d0limFD[iPt],d0limFD[iPt],-1000,1000);
      ImpParFitter->FitTree(-d0limFD[iPt],d0limFD[iPt],kFALSE);

      Double_t promptfrac = ImpParFitter->GetPromptFraction();
      Double_t promptfracerr = ImpParFitter->GetPromptFractionErr();
      hFrac[iSB]->SetBinContent(iPt+1,promptfrac);
      hFrac[iSB]->SetBinError(iPt+1,promptfracerr);
    }
  }

  TFile resultfile("fraction_unbinned_sigmafree.root","UPDATE");
  TH1F* hPromptFrac=(TH1F*)resultfile.Get("hFrac");
  hPromptFrac->SetDirectory(0);
  resultfile.Close();

  TLegend *l = new TLegend(0.4,0.2,0.89,0.45);
  l->SetTextSize(0.045);

  for(Int_t iSB=0; iSB<nSB; iSB++) {
    hRatio[iSB]->Divide(hFrac[iSB],hPromptFrac,1.,1.);
    for(Int_t iPt=0; iPt<hRatio[iSB]->GetNbinsX(); iPt++) {
      hRatio[iSB]->SetBinError(iPt+1,1.e-10);
    }
    l->AddEntry(hFrac[iSB], Form("|M-M_{PEAK}| > %d#sigma",SB[iSB]),"lpe");
  }

  TCanvas* cFrac = new TCanvas("cFrac","",800,800);
  hFrac[0]->Draw("E1");
  hFrac[0]->GetYaxis()->SetRangeUser(0,1.4);
  for(Int_t iSB=1; iSB<nSB; iSB++) {
    hFrac[iSB]->Draw("E1same");
  }
  l->Draw("same");

  TCanvas* cRatio = new TCanvas("cRatio","",800,800);
  hRatio[0]->Draw("E1");
  hRatio[0]->GetYaxis()->SetRangeUser(0.6,1.2);
  for(Int_t iSB=1; iSB<nSB; iSB++) {
    hRatio[iSB]->Draw("E1same");
  }
  l->Draw("same");
  line->SetLineColor(colors[0]);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw("same");

  cFrac->SaveAs("sidebands/promptfraction_syst_SBdist_onlyfrac.eps");
  cRatio->SaveAs("sidebands/promptfraction_syst_SBdist_onlyratio.eps");

  TFile outfile("sidebands/promptfraction_syst_SBdist.root","RECREATE");
  for(Int_t iSB=0; iSB<nSB; iSB++) {
    hFrac[iSB]->Write();
    hRatio[iSB]->Write();
  }
  outfile.Close();

}

void PrefitSystematics() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);

  //___________________________________________________________________________________________
  //range limits
  for(Int_t iBin=0; iBin<nPtBins; iBin++) {
    d0lim[iBin] = d0limFD[iBin];
  }
  //___________________________________________________________________________________________
  //input files

  TFile infileMC("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/MC/AnalysisResultspPbMC.root","READ");
  TDirectoryFile* dirMC=(TDirectoryFile*)infileMC.Get("PWG3_D2H_InvMassDplus");
  TList* listMC=(TList*)dirMC->Get("coutputDplus_ImpParpPbMC0100");
  THnSparseF* hMassPtImpPrompt=(THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
  THnSparseF* hMassPtImpRecoFD=(THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
  THnSparseF* hMassPtImpTrueFD=(THnSparseF*)listMC->FindObject("hMassPtImpParTrueBfeed");
  infileMC.Close();

  TFile infileData("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root","READ");
  TDirectoryFile* dirData=(TDirectoryFile*)infileData.Get("PWG3_D2H_InvMassDplus");
  TList* listData=(TList*)dirData->Get("coutputDplus_ImpParpPbData0100");
  THnSparseF* hMassPtImpParAll=(THnSparseF*)listData->FindObject("hMassPtImpParAll");
  TNtuple* dataTree=(TNtuple*)dirData->Get("fNtupleDplus");

  //infileData.Close();
  //____________________________________________________________________________________________
  //analysis

  TH1F* hFrac = new TH1F("hFrac","",nPtBins,PtLims);
  TH1F* hRatio = new TH1F("hRatio","",nPtBins,PtLims);
  hFrac->SetDirectory(0);
  hFrac->SetLineColor(colors[1]);
  hFrac->SetLineWidth(2);
  hFrac->SetMarkerStyle(20);
  hFrac->SetMarkerSize(1.5);
  hFrac->SetMarkerColor(colors[1]);
  hRatio->SetDirectory(0);
  hRatio->SetLineColor(colors[1]);
  hRatio->SetLineWidth(2);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(1.5);
  hRatio->SetMarkerColor(colors[1]);
  hFrac->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hFrac->GetYaxis()->SetTitle("#it{f}_{prompt}");
  hRatio->GetYaxis()->SetTitle("#frac{asymmetric}{symmetric}");
  hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");

  AliDplusCharmFractionIPfitter *ImpParFitter = new AliDplusCharmFractionIPfitter();
  ImpParFitter->SetDataSparse(hMassPtImpParAll);
  ImpParFitter->SetMCPromptSparse(hMassPtImpPrompt);
  ImpParFitter->SetMCTrueFDSparse(hMassPtImpTrueFD);
  ImpParFitter->SetMCRecoFDSparse(hMassPtImpRecoFD);
  ImpParFitter->SetDataTree(dataTree);

  ImpParFitter->SetFDFunction(AliDplusCharmFractionIPfitter::kConvolution);
  ImpParFitter->SetBkgFunction(AliDplusCharmFractionIPfitter::kDoubleGaussExpo);
  ImpParFitter->SetSideBandsRegion(AliDplusCharmFractionIPfitter::kBoth);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ImpParFitter->SetRelativeLimitSigmaPrompt(0.2);
    ImpParFitter->SetLimitsPromptFrac(0.2,1.5);
    ImpParFitter->SetInitialParameters(initparprompt,initparFD,initparBkg);
    ImpParFitter->FixSigmaPromptFromMC(kFALSE);
    ImpParFitter->SetRebinImpParHistos(2);
    ImpParFitter->SetPtLims(PtLims[iPt],PtLims[iPt+1]);
    ImpParFitter->SetNSigmas(2.);
    ImpParFitter->GetSignal(4,0,0,1.68,2.05,AliDplusCharmFractionIPfitter::kCentralValue);
    ImpParFitter->PrefitStep(-d0limPrompt[iPt],d0limPrompt[iPt],-d0limFD[iPt],d0limFD[iPt],-1000,1000);
    ImpParFitter->FitTree(-d0limFD[iPt],d0limFD[iPt],kFALSE);

    Double_t promptfrac = ImpParFitter->GetPromptFraction();
    Double_t promptfracerr = ImpParFitter->GetPromptFractionErr();
    hFrac->SetBinContent(iPt+1,promptfrac);
    hFrac->SetBinError(iPt+1,promptfracerr);
  }

  TFile resultfile("fraction_unbinned_sigmafree.root","UPDATE");
  TH1F* hPromptFrac=(TH1F*)resultfile.Get("hFrac");
  hPromptFrac->SetDirectory(0);
  hPromptFrac->SetLineColor(colors[0]);
  resultfile.Close();

  hRatio->Divide(hFrac,hPromptFrac,1.,1.);
  for(Int_t iPt=0; iPt<hRatio->GetNbinsX(); iPt++) {
    hRatio->SetBinError(iPt+1,1.e-10);
  }

  TLegend *l = new TLegend(0.4,0.2,0.89,0.45);
  l->SetTextSize(0.045);
  l->AddEntry(hPromptFrac, "SB distr symmetric", "lpe");
  l->AddEntry(hFrac, "SB distr asymmetric", "lpe");

  TCanvas* cFrac = new TCanvas("cFrac","",800,800);
  hFrac->GetYaxis()->SetRangeUser(0.2,1.4);
  hFrac->Draw("E1");
  hPromptFrac->Draw("E1same");
  l->Draw("same");

  TCanvas* cRatio = new TCanvas("cRatio","",800,800);
  hRatio->GetYaxis()->SetRangeUser(0.9,1.1);
  hRatio->Draw("E1");
  line->SetLineColor(colors[0]);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw("same");

  cFrac->SaveAs("prefit/promptfraction_syst_prefit_onlyfrac.eps");
  cRatio->SaveAs("prefit/promptfraction_syst_prefit_onlyratio.eps");

  TFile outfile("prefit/promptfraction_syst_prefit.root","RECREATE");
  hFrac->Write();
  hRatio->Write();
  outfile.Close();

}

void BinningSystematics() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);

  //___________________________________________________________________________________________
  //range limits
  for(Int_t iBin=0; iBin<nPtBins; iBin++) {
    d0lim[iBin] = d0limFD[iBin];
  }
  //___________________________________________________________________________________________
  //input files

  TFile infileMC("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/MC/AnalysisResultspPbMC.root","READ");
  TDirectoryFile* dirMC=(TDirectoryFile*)infileMC.Get("PWG3_D2H_InvMassDplus");
  TList* listMC=(TList*)dirMC->Get("coutputDplus_ImpParpPbMC0100");
  THnSparseF* hMassPtImpPrompt=(THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
  THnSparseF* hMassPtImpRecoFD=(THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
  THnSparseF* hMassPtImpTrueFD=(THnSparseF*)listMC->FindObject("hMassPtImpParTrueBfeed");
  infileMC.Close();

  TFile infileData("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root","READ");
  TDirectoryFile* dirData=(TDirectoryFile*)infileData.Get("PWG3_D2H_InvMassDplus");
  TList* listData=(TList*)dirData->Get("coutputDplus_ImpParpPbData0100");
  THnSparseF* hMassPtImpParAll=(THnSparseF*)listData->FindObject("hMassPtImpParAll");
  TNtuple* dataTree=(TNtuple*)dirData->Get("fNtupleDplus");

  //infileData.Close();
  //____________________________________________________________________________________________
  //analysis
  Int_t nBinSize = 3;

  Int_t binsize[nBinSize] = {1,2,4};

  TH1F** hFrac = new TH1F*[nBinSize];
  TH1F** hRatio = new TH1F*[nBinSize];
  for(Int_t iBinSize=0; iBinSize<nBinSize; iBinSize++) {
    hFrac[iBinSize] = new TH1F(Form("hFrac%d",iBinSize),"",nPtBins,PtLims);
    hFrac[iBinSize]->SetDirectory(0);
    hFrac[iBinSize]->SetLineColor(colors[iBinSize+1]);
    hFrac[iBinSize]->SetLineWidth(2);
    hFrac[iBinSize]->SetMarkerColor(colors[iBinSize+1]);
    hFrac[iBinSize]->SetMarkerSize(1.5);
    hFrac[iBinSize]->SetMarkerStyle(20+iBinSize);
    hRatio[iBinSize] = new TH1F(Form("hRatio%d",iBinSize),"",nPtBins,PtLims);
    hRatio[iBinSize]->SetDirectory(0);
    hRatio[iBinSize]->SetLineColor(colors[iBinSize+1]);
    hRatio[iBinSize]->SetLineWidth(2);
    hRatio[iBinSize]->SetMarkerColor(colors[iBinSize+1]);
    hRatio[iBinSize]->SetMarkerSize(1.5);
    hRatio[iBinSize]->SetMarkerStyle(20+iBinSize);
    hFrac[iBinSize]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hFrac[iBinSize]->GetYaxis()->SetTitle("#it{f}_{prompt}");
    hRatio[iBinSize]->GetYaxis()->SetTitle("ratio of #it{f}_{prompt} w.r.t the central value");
    hRatio[iBinSize]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  }

  AliDplusCharmFractionIPfitter *ImpParFitter = new AliDplusCharmFractionIPfitter();
  ImpParFitter->SetDataSparse(hMassPtImpParAll);
  ImpParFitter->SetMCPromptSparse(hMassPtImpPrompt);
  ImpParFitter->SetMCTrueFDSparse(hMassPtImpTrueFD);
  ImpParFitter->SetMCRecoFDSparse(hMassPtImpRecoFD);
  ImpParFitter->SetDataTree(dataTree);

  ImpParFitter->SetSideBandsRegion(AliDplusCharmFractionIPfitter::kBoth);
  ImpParFitter->SetFDFunction(AliDplusCharmFractionIPfitter::kConvolution);
  ImpParFitter->SetBkgFunction(AliDplusCharmFractionIPfitter::kDoubleGaussExpoSymm);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ImpParFitter->SetPtLims(PtLims[iPt],PtLims[iPt+1]);
    ImpParFitter->SetNSigmas(2.);
    ImpParFitter->GetSignal(4,0,0,1.68,2.05,AliDplusCharmFractionIPfitter::kCentralValue);

    for(Int_t iBinSize=0; iBinSize<nBinSize; iBinSize++) {
      ImpParFitter->SetRebinImpParHistos(binsize[iBinSize]);
      ImpParFitter->PrefitStep(-d0limPrompt[iPt],d0limPrompt[iPt],-d0limFD[iPt],d0limFD[iPt],-1000,1000);
      ImpParFitter->SetFitOptions("RLEM0");
      ImpParFitter->SetRelativeLimitSigmaPrompt(0.2);
      ImpParFitter->SetLimitsPromptFrac(0.2,1.5);
      ImpParFitter->SetInitialParameters(initparprompt,initparFD,initparBkg);
      ImpParFitter->FixSigmaPromptFromMC(kFALSE);
      ImpParFitter->FitHisto(-d0limFD[iPt],d0limFD[iPt],kFALSE);      

      Double_t promptfrac = ImpParFitter->GetPromptFraction();
      Double_t promptfracerr = ImpParFitter->GetPromptFractionErr();
      hFrac[iBinSize]->SetBinContent(iPt+1,promptfrac);
      hFrac[iBinSize]->SetBinError(iPt+1,promptfracerr);
    }
  }

  TFile resultfile("fraction_unbinned_sigmafree.root","UPDATE");
  TH1F* hPromptFrac=(TH1F*)resultfile.Get("hFrac");
  hPromptFrac->SetDirectory(0);

  TLegend *l = new TLegend(0.6,0.2,0.89,0.45);

  for(Int_t iBinSize=0; iBinSize<nBinSize; iBinSize++) {
    hRatio[iBinSize]->Divide(hFrac[iBinSize],hPromptFrac,1.,1.);
    for(Int_t iPt=0; iPt<hRatio[iBinSize]->GetNbinsX(); iPt++) {
      hRatio[iBinSize]->SetBinError(iPt+1,1.e-10);
    }
    Double_t size = 2000./400*binsize[iBinSize];
    l->AddEntry(hFrac[iBinSize], Form("bin size = %0.f #mum",size), "lpe");
  }

  TCanvas* cFrac = new TCanvas("cFrac","",800,800);
  hFrac[0]->Draw("E1");
  hFrac[0]->GetYaxis()->SetRangeUser(0,1.4);
  hFrac[1]->Draw("E1same");
  hFrac[2]->Draw("E1same");
  l->Draw("same");

  TCanvas* cRatio = new TCanvas("cRatio","",800,800);
  hRatio[0]->Draw("E1");
  hRatio[0]->GetYaxis()->SetRangeUser(0,1.8);
  hRatio[1]->Draw("E1same");
  hRatio[2]->Draw("E1same");
  l->Draw("same");
  line->SetLineColor(colors[0]);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw("same");

  cFrac->SaveAs("binning/promptfraction_syst_binning_onlyfrac.eps");
  cRatio->SaveAs("binning/promptfraction_syst_binning_onlyratio.eps");

  TFile outfile("binning/promptfraction_syst_binning.root","RECREATE");
  for(Int_t iBinSize=0; iBinSize<nBinSize; iBinSize++) {
    hFrac[iBinSize]->Write();
    hRatio[iBinSize]->Write();
  }
  outfile.Close();

}

void PtReweightSystematics() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);

  //___________________________________________________________________________________________
  //range limits
  for(Int_t iBin=0; iBin<nPtBins; iBin++) {
    d0lim[iBin] = d0limFD[iBin];
  }
  //___________________________________________________________________________________________
  //input files

  TFile infilePtB("/home/fabrizio/tesi/PtReweights/pPb_LHC13/PtB/PtB.root","READ");
  TH1F* hPtBWeights = (TH1F*)infilePtB.Get("hWeights");
  hPtBWeights->SetDirectory(0);
  infilePtB.Close();

  TFile infilePtD("/home/fabrizio/tesi/PtReweights/pPb_LHC13/PtD/PtD.root","READ");
  TH1F* hPtDWeights = (TH1F*)infilePtD.Get("hWeights");
  hPtDWeights->SetDirectory(0);
  infilePtD.Close();

  TFile infileMC("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/MC/AnalysisResultspPbMC.root","READ");
  TDirectoryFile* dirMC=(TDirectoryFile*)infileMC.Get("PWG3_D2H_InvMassDplus");
  TList* listMC=(TList*)dirMC->Get("coutputDplus_ImpParpPbMC0100");
  THnSparseF* hMassPtImpPrompt=(THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
  THnSparseF* hMassPtImpRecoFD=(THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
  THnSparseF* hMassPtImpTrueFD=(THnSparseF*)listMC->FindObject("hMassPtImpParTrueBfeed");
  infileMC.Close();

  TFile infileData("/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root","READ");
  TDirectoryFile* dirData=(TDirectoryFile*)infileData.Get("PWG3_D2H_InvMassDplus");
  TList* listData=(TList*)dirData->Get("coutputDplus_ImpParpPbData0100");
  THnSparseF* hMassPtImpParAll=(THnSparseF*)listData->FindObject("hMassPtImpParAll");
  TNtuple* dataTree=(TNtuple*)dirData->Get("fNtupleDplus");

  //infileData.Close();

  //____________________________________________________________________________________________
  //analysis

  AliDplusCharmFractionIPfitter *ImpParFitter = new AliDplusCharmFractionIPfitter();
  ImpParFitter->SetDataSparse(hMassPtImpParAll);
  ImpParFitter->SetMCPromptSparse(hMassPtImpPrompt);
  ImpParFitter->SetMCTrueFDSparse(hMassPtImpTrueFD);
  ImpParFitter->SetMCRecoFDSparse(hMassPtImpRecoFD);
  ImpParFitter->SetDataTree(dataTree);

  ImpParFitter->SetSideBandsRegion(AliDplusCharmFractionIPfitter::kBoth);
  ImpParFitter->SetFDFunction(AliDplusCharmFractionIPfitter::kConvolution);
  ImpParFitter->SetBkgFunction(AliDplusCharmFractionIPfitter::kDoubleGaussExpoSymm);

  Int_t nReweights = 3;

  TH1F** hFrac = new TH1F*[nReweights];
  TH1F** hRatio = new TH1F*[nReweights];

  for(Int_t iRew=0; iRew<nReweights; iRew++) {

    if(iRew==0) {
      ImpParFitter->SetPtDWeightsHisto(hPtDWeights);
      ImpParFitter->SetPtBWeights(kFALSE);
    }
    else if(iRew==1) {
      ImpParFitter->SetPtBWeightsHisto(hPtBWeights);
      ImpParFitter->SetPtDWeights(kFALSE);
    }
    else {
      ImpParFitter->SetPtDWeightsHisto(hPtDWeights);
      ImpParFitter->SetPtBWeightsHisto(hPtBWeights);
    }

    hFrac[iRew] = new TH1F(Form("hFrac%d",iRew),"",nPtBins,PtLims);
    hFrac[iRew]->SetDirectory(0);
    hFrac[iRew]->SetLineColor(colors[iRew+1]);
    hFrac[iRew]->SetLineWidth(2);
    hFrac[iRew]->SetMarkerColor(colors[iRew+1]);
    hFrac[iRew]->SetMarkerSize(1.5);
    hFrac[iRew]->SetMarkerStyle(20+iRew);
    hRatio[iRew] = new TH1F(Form("hRatio%d",iRew),"",nPtBins,PtLims);
    hRatio[iRew]->SetDirectory(0);
    hRatio[iRew]->SetLineColor(colors[iRew+1]);
    hRatio[iRew]->SetLineWidth(2);
    hRatio[iRew]->SetMarkerColor(colors[iRew+1]);
    hRatio[iRew]->SetMarkerSize(1.5);
    hRatio[iRew]->SetMarkerStyle(20+iRew);
    hFrac[iRew]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hFrac[iRew]->GetYaxis()->SetTitle("#it{f}_{prompt}");
    hRatio[iRew]->GetYaxis()->SetTitle("ratio of #it{f}_{prompt} w.r.t. the central value");
    hRatio[iRew]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");

    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      ImpParFitter->SetRelativeLimitSigmaPrompt(0.2);
      ImpParFitter->SetLimitsPromptFrac(0.2,1.5);
      ImpParFitter->SetInitialParameters(initparprompt,initparFD,initparBkg);
      ImpParFitter->FixSigmaPromptFromMC(kFALSE);
      ImpParFitter->SetRebinImpParHistos(2);
      ImpParFitter->SetPtLims(PtLims[iPt],PtLims[iPt+1]);
      ImpParFitter->SetNSigmas(2.);
      ImpParFitter->SetNSigmaSBLimits(3,15);
      ImpParFitter->GetSignal(4,0,0,1.68,2.05,AliDplusCharmFractionIPfitter::kCentralValue);
      ImpParFitter->CheckSideBandsImpParDist(5);
      ImpParFitter->PrefitStep(-d0limPrompt[iPt],d0limPrompt[iPt],-d0limFD[iPt],d0limFD[iPt],-1000,1000);
      ImpParFitter->FitTree(-d0lim[iPt],d0lim[iPt],kFALSE);

      Double_t promptfrac = ImpParFitter->GetPromptFraction();
      Double_t promptfracerr = ImpParFitter->GetPromptFractionErr();
      hFrac[iRew]->SetBinContent(iPt+1,promptfrac);
      hFrac[iRew]->SetBinError(iPt+1,promptfracerr);
    }
  }

  TFile resultfile("fraction_unbinned_sigmafree.root","UPDATE");
  TH1F* hPromptFrac=(TH1F*)resultfile.Get("hFrac");
  hPromptFrac->SetDirectory(0);

  TLegend *l = new TLegend(0.35,0.15,0.89,0.4);
  l->SetTextSize(0.045);

  for(Int_t iRew=0; iRew<nReweights; iRew++) {
    hRatio[iRew]->Divide(hFrac[iRew],hPromptFrac,1.,1.);
    for(Int_t iPt=0; iPt<hRatio[iRew]->GetNbinsX(); iPt++) {
      hRatio[iRew]->SetBinError(iPt+1,1.e-10);
    }
  }

  l->AddEntry(hFrac[0], "#it{p}_{T}^{D} reweighted", "lpe");
  l->AddEntry(hFrac[1], "#it{p}_{T}^{B} reweighted", "lpe");
  l->AddEntry(hFrac[2], "#it{p}_{T}^{D} and #it{p}_{T}^{B} reweighted", "lpe");

  TCanvas* cFrac = new TCanvas("cFrac","",800,800);
  hFrac[0]->Draw("E1");
  hFrac[0]->GetYaxis()->SetRangeUser(0,1.4);
  hFrac[1]->Draw("E1same");
  hFrac[2]->Draw("E1same");
  l->Draw("same");

  TCanvas* cRatio = new TCanvas("cRatio","",800,800);
  hRatio[0]->Draw("E1");
  hRatio[0]->GetYaxis()->SetRangeUser(0.99,1.005);
  hRatio[1]->Draw("E1same");
  hRatio[2]->Draw("E1same");
  l->Draw("same");

  cFrac->SaveAs("ptreweight/promptfraction_syst_pTreweight_onlyfrac.eps");
  cRatio->SaveAs("ptreweight/promptfraction_syst_pTreweight_onlyratio.eps");

  TFile outfile("ptreweight/promptfraction_syst_pTreweight.root","RECREATE");
  for(Int_t iRew=0; iRew<nReweights; iRew++) {
    hFrac[iRew]->Write();
    hRatio[iRew]->Write();
  }
  outfile.Close();
}
