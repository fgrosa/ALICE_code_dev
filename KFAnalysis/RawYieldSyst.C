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
#include <TPaveText.h>
#include <TPad.h>
#include <TGaxis.h>

#include "AliHFMassFitter.h"

#endif

enum {kDzero,kDplus};
Int_t minPID=2;
Int_t maxPID=2;
Int_t PIDaxnum=3;

void RawYieldSyst(Int_t meson=kDplus,
                  TString reffilename="../result/rawyields_Dplus_cutset1.root",
                  TString reffileMCname="../result/rawyieldsMC_Dplus_cutset1.root",
                  TString axesfile="../axes",
                  TString cutsfile="../cutset1");

void UsePID(THnSparseF* sparse);
void ApplyTopologicalCuts(THnSparseF* sparse, Int_t iPt, TString axesfile, TString cutsfile);
void SetPtRange(THnSparseF* sparse, Int_t iPt, Int_t iAxis, TString axesfile, TString cutsfile);
void ResetAxes(THnSparseF* sparse);
void ReadAxesNum(TString Filename, vector<string> &axesanmes, vector<int> &axesno);
void ReadSet(TString Filename, vector<string> &varnames, vector<double> &cutset);
Int_t GetNPtBins(TString axesfile, TString cutsfile);
Double_t* GetPtLims(TString axesfile, TString cutsfile);
Double_t GetMax(TH1F* h);
Double_t GetMin(TH1F* h);
void PrintCuts(TString axesfile, TString cutsfile);

const Int_t colors[] = {kGreen+3,kBlue+2,kBlack,kRed,kOrange+7};

void RawYieldSyst(Int_t meson, TString reffilename, TString reffileMCname, TString axesfile, TString cutsfile) {

  Int_t pdgcode=411;
  if(meson==kDzero)
    pdgcode=421;
  Double_t massD = TDatabasePDG::Instance()->GetParticle(pdgcode)->Mass();

  TString mesonname;
  if(meson==kDzero) mesonname = "Dzero";
  else mesonname = "Dplus";
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.6,"y");
  gStyle->SetLegendBorderSize(0);
  TGaxis::SetMaxDigits(3);
  
  //INPUT FILES________________________________________________________________
  TString Datafilename;
  TString listname;
  
  if(meson==kDzero) {
    Datafilename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pPb_LHC13/Data/fgrosa_Dzero_KF_pPb.root";
    listname = "coutputDzeroKF";
  }
  else {
    Datafilename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pPb_LHC13/Data/fgrosa_Dplus_KF_pPb.root";
    listname = "coutputDplusKF";
  }


  TFile reffile(reffilename,"UPDATE");
  TH1F* hRawRef = (TH1F*)reffile.Get("hSignal");
  hRawRef->SetDirectory(0);
  reffile.Close();
  
  TFile reffileMC(reffileMCname,"UPDATE");
  TH1F* hSigmaMC = (TH1F*)reffileMC.Get("hSigma");
  hSigmaMC->SetDirectory(0);
  reffileMC.Close();
  
  //PRINT CUTS VALUES__________________________________________________________
  PrintCuts(axesfile,cutsfile);
  
  TFile infileData(Datafilename.Data(),"UPDATE");
  TList* listData = (TList*)infileData.Get(listname.Data());
  THnSparseF *AllSparse = (THnSparseF*)listData->FindObject("fSparseAll");

  const Int_t nPtBins=GetNPtBins(axesfile,cutsfile);
  const Int_t nPtLims=nPtBins+1;
  Double_t *PtLims = (Double_t*)GetPtLims(axesfile,cutsfile);
  
  TH1F** hMass = new TH1F*[nPtBins];
  TH1F** hRawYieldVsTrial = new TH1F*[nPtBins];
  TH1F** hRawYieldDist = new TH1F*[nPtBins];
  TH1F** hCntVsTrial = new TH1F*[nPtBins];
  TH1F** hCntDist = new TH1F*[nPtBins];
  TH1F** hSigmaVsTrial = new TH1F*[nPtBins];
  TH1F** hChiVsTrial = new TH1F*[nPtBins];
  TH1F** hMeanVsTrial = new TH1F*[nPtBins];
  
  const Int_t nBinVar=4;
  Int_t rebin[nBinVar] = {3,4,5,6};
  const Int_t nMins=6;
  Double_t mins[nMins] = {1.68,1.69,1.70,1.71,1.72,1.73};
  const Int_t nMaxs=6;
  Double_t maxs[nMaxs] = {2.05,2.04,2.03,2.02,2.01,2.00};
  const Int_t nBkgFunc=3;
  Int_t bkgfunc[nBkgFunc] = {AliHFMassFitter::kExpo,AliHFMassFitter::kLin,AliHFMassFitter::kPol2}; 
  
  Int_t nBins=40;
  Int_t nTrials = nBinVar*nMins*nMaxs*2*nBkgFunc;

  Double_t perc[nPtBins] = {0.5,0.3,0.2,0.2,0.15,0.2,0.2,0.3,0.4,0.6};
  Double_t percw[nPtBins] = {0.3,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.3};
  Double_t smin[nPtBins];
  Double_t smax[nPtBins];
  Double_t wmin[nPtBins];
  Double_t wmax[nPtBins];
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    smin[iPt]=hRawRef->GetBinContent(iPt+1)-perc[iPt]*hRawRef->GetBinContent(iPt+1);
    smax[iPt]=hRawRef->GetBinContent(iPt+1)+perc[iPt]*hRawRef->GetBinContent(iPt+1);
    wmin[iPt]=hSigmaMC->GetBinContent(iPt+1)-percw[iPt]*hSigmaMC->GetBinContent(iPt+1);
    wmax[iPt]=hSigmaMC->GetBinContent(iPt+1)+percw[iPt]*hSigmaMC->GetBinContent(iPt+1);
    
    hRawYieldDist[iPt] = new TH1F(Form("hRawYieldDist_Pt%d",iPt),"",nBins,smin[iPt],smax[iPt]);
    hRawYieldVsTrial[iPt] = new TH1F(Form("hRawYieldVsTrial_Pt%d",iPt),"",nTrials,0.5,nTrials+0.5);
    hCntDist[iPt] = new TH1F(Form("hCntDist_Pt%d",iPt),"",nBins,smin[iPt],smax[iPt]);
    hCntVsTrial[iPt] = new TH1F(Form("hCntVsTrial_Pt%d",iPt),"",nTrials,0.5,nTrials+0.5);
    hSigmaVsTrial[iPt] = new TH1F(Form("hSigmaVsTrial_Pt%d",iPt),"",nTrials,0.5,nTrials+0.5);
    hMeanVsTrial[iPt] = new TH1F(Form("hMeanVsTrial_Pt%d",iPt),"",nTrials,0.5,nTrials+0.5);
    hChiVsTrial[iPt] = new TH1F(Form("hChiVsTrial_Pt%d",iPt),"",nTrials,0.5,nTrials+0.5);

    TString ptrange = Form("  %0.f < #it{p}_{T} < %0.f GeV/c",PtLims[iPt],PtLims[iPt+1]);
    
    hRawYieldDist[iPt]->GetXaxis()->SetTitle("Raw Yield");
    hRawYieldDist[iPt]->GetYaxis()->SetTitle("Entires");
    hRawYieldDist[iPt]->SetTitle(ptrange.Data());
    hRawYieldVsTrial[iPt]->GetYaxis()->SetTitle("Raw Yield");
    hRawYieldVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hRawYieldVsTrial[iPt]->SetTitle(ptrange.Data());
    hCntDist[iPt]->GetXaxis()->SetTitle("Raw Yield");
    hCntDist[iPt]->GetYaxis()->SetTitle("Entires");
    hCntDist[iPt]->SetTitle(ptrange.Data());
    hCntVsTrial[iPt]->GetYaxis()->SetTitle("Raw Yield");
    hCntVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hCntVsTrial[iPt]->SetTitle(ptrange.Data());
    hSigmaVsTrial[iPt]->GetYaxis()->SetTitle("Sigma (GeV/c^{2})");
    hSigmaVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hSigmaVsTrial[iPt]->SetTitle(ptrange.Data());
    hMeanVsTrial[iPt]->GetYaxis()->SetTitle("Mean (GeV/c^{2})");
    hMeanVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hMeanVsTrial[iPt]->SetTitle(ptrange.Data());
    hChiVsTrial[iPt]->GetYaxis()->SetTitle("#chi^{2}/ndf");
    hChiVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hChiVsTrial[iPt]->SetTitle(ptrange.Data());

    hRawYieldVsTrial[iPt]->GetYaxis()->SetRangeUser(smin[iPt]*1.1,smax[iPt]*1.35); 
    hRawYieldVsTrial[iPt]->GetYaxis()->SetNdivisions(508);
    hRawYieldVsTrial[iPt]->GetXaxis()->SetNdivisions(508); 
    hCntVsTrial[iPt]->GetYaxis()->SetRangeUser(smin[iPt]*0.65,smax[iPt]*1.3); 
    hSigmaVsTrial[iPt]->GetYaxis()->SetRangeUser(wmin[iPt],wmax[iPt]);
    hMeanVsTrial[iPt]->GetYaxis()->SetRangeUser(massD-0.005*massD,massD+0.005*massD);
    hChiVsTrial[iPt]->GetYaxis()->SetRangeUser(0.,2.5);
      
    ResetAxes(AllSparse);
    SetPtRange(AllSparse,iPt,1,axesfile,cutsfile);
    UsePID(AllSparse);
    ApplyTopologicalCuts(AllSparse,iPt,axesfile,cutsfile);
    hMass[iPt]=(TH1F*)AllSparse->Projection(0);
    hMass[iPt]->SetDirectory(0);

    Int_t iTrial=0;
    
    for(Int_t iSigma=0; iSigma<2; iSigma++) {
      for(Int_t iMean=0; iMean<1; iMean++) {
        for(Int_t iBinning=0; iBinning<nBinVar; iBinning++) {
          for(Int_t iMin=0; iMin<nMins; iMin++) {
            for(Int_t iMax=0; iMax<nMaxs; iMax++) {
              for(Int_t iBkg=0; iBkg<nBkgFunc; iBkg++) {
                if((iTrial==173 || iTrial==23) && iPt==9) continue;
                cout <<iPt << "  " << iTrial << "  " << iMean << "  " << iSigma << "  " << iBinning << "  " << iMin << "  " << iMax << "  " << iBkg << endl;
                TH1F* h=(TH1F*)hMass[iPt]->Clone();
                h->Rebin(rebin[iBinning]);
                AliHFMassFitter* fitter= new AliHFMassFitter(h,mins[iMin],maxs[iMax],1,bkgfunc[iBkg],AliHFMassFitter::kGaus);
                fitter->SetUseLikelihoodFit();
                if(iMean==0)
                  fitter->SetInitialGaussianMean(massD);
                else
                  fitter->SetFixGaussianMean(massD);
                if(iSigma==0)
                  fitter->SetInitialGaussianSigma(hSigmaMC->GetBinContent(iPt+1));
                else
                  fitter->SetFixGaussianSigma(hSigmaMC->GetBinContent(iPt+1));
                fitter->MassFitter(kFALSE);
                
                Double_t rawyield;
                Double_t rawyielderr;
                Double_t bkg;
                Double_t bkgerr;
                Double_t signif;
                Double_t signiferr;
                
                fitter->Signal(3,rawyield,rawyielderr);
                fitter->Background(3,bkg,bkgerr);
                fitter->Significance(3,signif,signiferr);
                Double_t sigma=fitter->GetSigma();
                Double_t sigmaerr=fitter->GetSigmaUncertainty();
                Double_t mean=fitter->GetMean();
                Double_t meanerr=fitter->GetMeanUncertainty();
                Double_t chiS=fitter->GetReducedChiSquare();
                
                Double_t ry=fitter->GetRawYield();
                Double_t ery=fitter->GetRawYieldError();
                Double_t cntSig1=0.;
                Double_t cntSig2=0.;
                Double_t cntErr=0.;
                Double_t massRangeForCounting=3.5*sigma;
                Double_t minBinSum=h->FindBin(mean-massRangeForCounting);
                Double_t maxBinSum=h->FindBin(mean+massRangeForCounting);
                TF1* fB1=fitter->GetBackgroundFullRangeFunc();
                TF1* fB2=fitter->GetBackgroundRecalcFunc();
                
                for(Int_t iMB=minBinSum; iMB<=maxBinSum; iMB++){
                  Double_t bkg1=fB1 ? fB1->Eval(h->GetBinCenter(iMB)) : 0;
                  Double_t bkg2=fB2 ? fB2->Eval(h->GetBinCenter(iMB)) : 0;
                  cntSig1+=(h->GetBinContent(iMB)-bkg1);
                  cntSig2+=(h->GetBinContent(iMB)-bkg2);
                  cntErr+=(h->GetBinContent(iMB));  
                }
                
                hCntDist[iPt]->Fill(cntSig1);
                hCntVsTrial[iPt]->SetBinContent(iTrial+1,cntSig1);
                hCntVsTrial[iPt]->SetBinError(iTrial+1,TMath::Sqrt(cntSig1));
                
                hRawYieldDist[iPt]->Fill(rawyield);
                hRawYieldVsTrial[iPt]->SetBinContent(iTrial+1,rawyield);
                hRawYieldVsTrial[iPt]->SetBinError(iTrial+1,rawyielderr);
                hSigmaVsTrial[iPt]->SetBinContent(iTrial+1,sigma);
                hSigmaVsTrial[iPt]->SetBinError(iTrial+1,sigmaerr);
                hMeanVsTrial[iPt]->SetBinContent(iTrial+1,mean);
                hMeanVsTrial[iPt]->SetBinError(iTrial+1,meanerr);
                hChiVsTrial[iPt]->SetBinContent(iTrial+1,chiS);
                iTrial++;
              }
            }
          }
        }
      }
    }
  }

  TCanvas** cRawDist = new TCanvas*[nPtBins];
  TCanvas** cRawVsTrial = new TCanvas*[nPtBins];
  TCanvas** cSigmaVsTrial = new TCanvas*[nPtBins];
  TCanvas** cMeanVsTrial = new TCanvas*[nPtBins];
  TCanvas** cChiVsTrial = new TCanvas*[nPtBins];
  TLine** line = new TLine*[nPtBins];
  TLine** line2 = new TLine*[nPtBins];
  TPaveText** lat = new TPaveText*[nPtBins];
  TPaveText** latBin = new TPaveText*[nPtBins];


  TLegend* l = new TLegend(0.5,0.38,0.8,0.43);
  l->SetTextSize(0.045);
  l->SetTextFont(132);
  l->SetFillStyle(0);
  TLegend* l2 = 0x0;
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    line[iPt] = new TLine(hRawRef->GetBinContent(iPt+1),0,hRawRef->GetBinContent(iPt+1),hRawYieldDist[iPt]->GetMaximum()*0.8);
    line[iPt]->SetLineColor(colors[3]);
    line[iPt]->SetLineWidth(2);
    line2[iPt] = new TLine(0,hRawRef->GetBinContent(iPt+1),nTrials,hRawRef->GetBinContent(iPt+1));
    line2[iPt]->SetLineColor(colors[3]);
    line2[iPt]->SetLineWidth(2);
    if(iPt==0) {
      l->AddEntry(line[iPt],"Reference value","l");
      l2=(TLegend*)l->Clone();
      l2->SetY1(0.38);
      l2->SetY2(0.43);
    }
      
    hCntDist[iPt]->SetLineColor(colors[0]);
    hCntVsTrial[iPt]->SetLineColor(colors[0]);
    hCntVsTrial[iPt]->SetMarkerColor(colors[2]);
    hCntDist[iPt]->SetFillColor(colors[0]);
    hCntDist[iPt]->SetFillStyle(3004);
    hRawYieldDist[iPt]->SetFillColor(colors[1]);
    hRawYieldDist[iPt]->SetLineColor(colors[1]);
    hRawYieldVsTrial[iPt]->SetLineColor(colors[1]);
    hRawYieldVsTrial[iPt]->SetMarkerColor(colors[2]);
    hRawYieldDist[iPt]->SetFillStyle(3004);
    
    cRawDist[iPt] = new TCanvas(Form("cRawDist_Pt%d",iPt),"",800,800);
    hCntDist[iPt]->GetYaxis()->SetRangeUser(0,hCntDist[iPt]->GetMaximum()*1.3);
    hCntDist[iPt]->Draw();
    hRawYieldDist[iPt]->Draw("same");
    line[iPt]->Draw("same");
    cRawVsTrial[iPt] = new TCanvas(Form("cRawVsTrial_Pt%d",iPt),"",800,800);
    hCntVsTrial[iPt]->Draw();
    hRawYieldVsTrial[iPt]->Draw("same");
    Double_t min = GetMin(hRawYieldDist[iPt]);
    Double_t max = GetMax(hRawYieldDist[iPt]);
    Double_t ext = (max-min)/TMath::Sqrt(12);
    Double_t extperc = ext/hRawYieldDist[iPt]->GetMean()*100;
    Double_t minBin = GetMin(hCntDist[iPt]);
    Double_t maxBin = GetMax(hCntDist[iPt]);
    Double_t extBin = (maxBin-minBin)/TMath::Sqrt(12);
    Double_t extpercBin = extBin/hCntDist[iPt]->GetMean()*100;
    
    if(iPt<1 || iPt>=nPtBins-2)
      l->Draw("same");
    else
      l2->Draw("same");

    lat[iPt] = new TPaveText(0.2,0.18,0.6,0.39,"NDC");
    lat[iPt]->SetBorderSize(0);
    lat[iPt]->SetFillStyle(0);
    lat[iPt]->SetTextFont(132);
    lat[iPt]->SetTextSize(0.045);
    lat[iPt]->SetTextColor(colors[1]);
    lat[iPt]->SetTextAlign(11);
    latBin[iPt] = new TPaveText(0.2,0.64,0.6,0.85,"NDC");
    latBin[iPt]->SetBorderSize(0);
    latBin[iPt]->SetFillStyle(0);
    latBin[iPt]->SetTextFont(132);
    latBin[iPt]->SetTextSize(0.045);
    latBin[iPt]->SetTextColor(colors[0]);
    latBin[iPt]->SetTextAlign(11);

    lat[iPt]->AddText("Fit");
    lat[iPt]->AddText(Form("mean = %0.1f",hRawYieldDist[iPt]->GetMean()));
    lat[iPt]->AddText(Form("RMS = %0.1f (%0.1f%)",hRawYieldDist[iPt]->GetRMS(),hRawYieldDist[iPt]->GetRMS()*100/hRawYieldDist[iPt]->GetMean()));
    lat[iPt]->AddText(Form("(x_{max}-x_{min})/#sqrt{12} = %0.1f (%0.1f%)",ext,extperc));
    
    latBin[iPt]->AddText("Bin counting");
    latBin[iPt]->AddText(Form("mean = %0.1f",hCntDist[iPt]->GetMean()));
    latBin[iPt]->AddText(Form("RMS = %0.1f (%0.1f%)",hCntDist[iPt]->GetRMS(),hCntDist[iPt]->GetRMS()*100/hCntDist[iPt]->GetMean()));
    latBin[iPt]->AddText(Form("(x_{max}-x_{min})/#sqrt{12} = %0.1f (%0.1f%)",extBin,extpercBin));

    lat[iPt]->Draw("same");
    latBin[iPt]->Draw("same");
    
    line2[iPt]->Draw("same");
    cSigmaVsTrial[iPt] = new TCanvas(Form("cSigmaVsTrial_Pt%d",iPt),"",800,800);
    hSigmaVsTrial[iPt]->Draw();
    cMeanVsTrial[iPt] = new TCanvas(Form("cMeanVsTrial_Pt%d",iPt),"",800,800);
    hMeanVsTrial[iPt]->Draw();
    cChiVsTrial[iPt] = new TCanvas(Form("cChiVsTrial_Pt%d",iPt),"",800,800);
    hChiVsTrial[iPt]->Draw();
  }
  
  TFile outfile(Form("rawyieldsyst_%s.root",mesonname.Data()),"RECREATE");
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    hRawYieldDist[iPt]->Write();
    hRawYieldVsTrial[iPt]->Write();
    hCntDist[iPt]->Write();
    hCntVsTrial[iPt]->Write();
    hSigmaVsTrial[iPt]->Write();
    hMeanVsTrial[iPt]->Write();
    hChiVsTrial[iPt]->Write();
    cRawDist[iPt]->Write();
    cRawVsTrial[iPt]->Write();
    cSigmaVsTrial[iPt]->Write();
    cMeanVsTrial[iPt]->Write();
    cChiVsTrial[iPt]->Write();
    cRawDist[iPt]->SaveAs(Form("RawDist_Pt%d.eps",iPt));
    cRawVsTrial[iPt]->SaveAs(Form("RawVsTrial_Pt%d.eps",iPt));
    cSigmaVsTrial[iPt]->SaveAs(Form("SigmaVsTrial_Pt%d.eps",iPt));
    cMeanVsTrial[iPt]->SaveAs(Form("MeanVsTrial_Pt%d.eps",iPt));
    cChiVsTrial[iPt]->SaveAs(Form("ChiVsTrial_Pt%d.eps",iPt));
  }
  outfile.Close();

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    delete cRawDist[iPt];
    delete cRawVsTrial[iPt];
    delete cSigmaVsTrial[iPt];
    delete cMeanVsTrial[iPt];
    delete cChiVsTrial[iPt];
  }
}

void UsePID(THnSparseF* sparse) {
  TAxis* PIDax = (TAxis*)sparse->GetAxis(PIDaxnum);
  PIDax->SetRange(minPID,maxPID);
}

void ApplyTopologicalCuts(THnSparseF* sparse, Int_t iPt,TString axesfile,TString cutsfile) {
  vector<string> axesnames;
  vector<int> axesno;
  vector<string> cutvarnames;
  vector<double> cutset;
  ReadAxesNum(axesfile,axesnames,axesno);
  ReadSet(cutsfile,cutvarnames,cutset);
  const Int_t nVars = axesno.size();

  for(Int_t iVar=1; iVar<nVars; iVar++) {
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

Double_t GetMax(TH1F* h) {

  Double_t binmax;
  Bool_t foundMax = kFALSE;
  Int_t bincounterMax=h->GetNbinsX();
  
  while(!foundMax) {
    if(h->GetBinContent(bincounterMax)!=0)
      foundMax=kTRUE;
    bincounterMax--;
  }
  
  Double_t max = h->GetBinCenter(bincounterMax+1);

  return max;
}

Double_t GetMin(TH1F* h) {

  Double_t binmin;
  Bool_t foundMin = kFALSE;
  Int_t bincounterMin=0;
  
  while(!foundMin) {
    if(h->GetBinContent(bincounterMin+1)!=0)
      foundMin=kTRUE;
    bincounterMin++;
  }
  
  Double_t min = h->GetBinCenter(bincounterMin);

  return min;
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
