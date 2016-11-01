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

#endif

enum {kDzero,kDplus};
const Int_t colors[] = {kRed,kBlue,kGreen+3,kBlack,kMagenta,kCyan,kOrange+7,kYellow,kGreen};

const Int_t nPtBins = 11;
const Int_t nPtLims = nPtBins+1;
Double_t PtLims[nPtLims] = {0,1,2,3,4,5,6,7,8,12,16,24};

void DaughtersQuality(Int_t meson=kDplus);
TH1F* GetHistoVsPt(THnSparseF* sparse, Int_t axis, TString outfilename, Bool_t noFIT=kFALSE, Bool_t singlegauss=kTRUE);
void SetPtRange(THnSparseF* sparse, Double_t ptmin, Double_t ptmax);
void ResetAxes(THnSparseF* sparse);
Double_t DoubleGauss(Double_t *x, Double_t *pars);

void DaughtersQuality(Int_t meson) {

  TString mesonname;
  if(meson==kDzero) mesonname = "Dzero";
  else mesonname = "Dplus";
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleSize(0.045,"xyz");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.0,"y");
  gStyle->SetLegendBorderSize(0);
  
  //INPUT FILES________________________________________________________________
  TString filename;
  TString listname;
  
  if(meson==kDzero) {
    //filename = "../fgrosa_Dzero_KF.root";
    filename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pp_LHC10/MC/fgrosa_Dzero_KF_ppMC.root";
    listname = "coutputDzeroKF";
  }
  else {
    //filename = "../fgrosa_Dplus_KF.root";
    filename = "/home/fabrizio/tesi/GSI/KF/task_GSI/Trains/pp_LHC10/MC/fgrosa_Dplus_KF_ppMC.root";
    listname = "coutputDplusKF";
  }
  
  TFile infile(filename.Data(),"UPDATE");
  TList* list = (TList*)infile.Get(listname.Data());

  //Daughters QA_________________________________________________________________
  THnSparseF* SparseRes = (THnSparseF*)list->FindObject("fHistDres");
  THnSparseF* SparsePulls = (THnSparseF*)list->FindObject("fHistDpulls");

  infile.Close();
  TString outfilename = Form("Daughters_%s_checks.root",mesonname.Data());
  TFile outfile(outfilename.Data(),"RECREATE");
  outfile.Close();

  TH1F* hResX = GetHistoVsPt(SparseRes,0,outfilename,kTRUE);
  TH1F* hResY = GetHistoVsPt(SparseRes,1,outfilename,kTRUE);
  TH1F* hResZ = GetHistoVsPt(SparseRes,2,outfilename,kTRUE);
  TH1F* hResPx = GetHistoVsPt(SparseRes,3,outfilename,kTRUE);
  TH1F* hResPy = GetHistoVsPt(SparseRes,4,outfilename,kTRUE);
  TH1F* hResPz = GetHistoVsPt(SparseRes,5,outfilename,kTRUE);
  TH1F* hResPt = GetHistoVsPt(SparseRes,8,outfilename,kTRUE);
  hResX->SetName("hResX");
  hResY->SetName("hResY");
  hResZ->SetName("hResZ");
  hResPx->SetName("hResPx");
  hResPy->SetName("hResPy");
  hResPz->SetName("hResPz");
  hResPt->SetName("hResPz");
  hResX->SetMarkerStyle(20);
  hResY->SetMarkerStyle(20);
  hResZ->SetMarkerStyle(20);
  hResPx->SetMarkerStyle(20);
  hResPy->SetMarkerStyle(20);
  hResPz->SetMarkerStyle(20);
  hResX->SetMarkerColor(colors[0]);
  hResY->SetMarkerColor(colors[1]);
  hResZ->SetMarkerColor(colors[2]);
  hResPx->SetMarkerColor(colors[0]);
  hResPy->SetMarkerColor(colors[1]);
  hResPz->SetMarkerColor(colors[2]);
  hResX->SetLineColor(colors[0]);
  hResY->SetLineColor(colors[1]);
  hResZ->SetLineColor(colors[2]);
  hResPx->SetLineColor(colors[0]);
  hResPy->SetLineColor(colors[1]);
  hResPz->SetLineColor(colors[2]);
  
  TH1F* hPullsX = GetHistoVsPt(SparsePulls,0,outfilename);
  TH1F* hPullsY = GetHistoVsPt(SparsePulls,1,outfilename);
  TH1F* hPullsZ = GetHistoVsPt(SparsePulls,2,outfilename);
  TH1F* hPullsPx = GetHistoVsPt(SparsePulls,3,outfilename);
  TH1F* hPullsPy = GetHistoVsPt(SparsePulls,4,outfilename);
  TH1F* hPullsPz = GetHistoVsPt(SparsePulls,5,outfilename);
  TH1F* hPullsPt = GetHistoVsPt(SparsePulls,8,outfilename);
  hPullsX->SetName("hPullsX");
  hPullsY->SetName("hPullsY");
  hPullsZ->SetName("hPullsZ");
  hPullsPx->SetName("hPullsPx");
  hPullsPy->SetName("hPullsPy");
  hPullsPz->SetName("hPullsPz");
  hPullsPt->SetName("hPullsPz");
  hPullsX->SetMarkerStyle(20);
  hPullsY->SetMarkerStyle(20);
  hPullsZ->SetMarkerStyle(20);
  hPullsPx->SetMarkerStyle(20);
  hPullsPy->SetMarkerStyle(20);
  hPullsPz->SetMarkerStyle(20);
  hPullsX->SetMarkerColor(colors[0]);
  hPullsY->SetMarkerColor(colors[1]);
  hPullsZ->SetMarkerColor(colors[2]);
  hPullsPx->SetMarkerColor(colors[0]);
  hPullsPy->SetMarkerColor(colors[1]);
  hPullsPz->SetMarkerColor(colors[2]);
  hPullsX->SetLineColor(colors[0]);
  hPullsY->SetLineColor(colors[1]);
  hPullsZ->SetLineColor(colors[2]);
  hPullsPx->SetLineColor(colors[0]);
  hPullsPy->SetLineColor(colors[1]);
  hPullsPz->SetLineColor(colors[2]);

  TLegend *l = new TLegend(0.6,0.7,0.89,0.89);
  l->SetTextSize(0.04);
  l->AddEntry(hResX,"decay vertex X","lpe");
  l->AddEntry(hResY,"decay vertex Y","lpe");
  l->AddEntry(hResZ,"decay vertex Z","lpe");
  TLegend *lp = new TLegend(0.2,0.7,0.5,0.89);
  lp->SetTextSize(0.04);
  lp->AddEntry(hResX,"p_{X}","lpe");
  lp->AddEntry(hResY,"p_{Y}","lpe");
  lp->AddEntry(hResZ,"p_{Z}","lpe");
  
  TCanvas *cResPos = new TCanvas("cResPos","",1200,900);
  hResX->GetYaxis()->SetTitle("#sigma(pos_{MC}-pos_{reco}) (cm)");
  hResX->GetYaxis()->SetRangeUser(0.012,0.03);
  hResX->Draw();
  hResY->Draw("same");
  hResZ->Draw("same");
  l->Draw("same");
  
  TCanvas *cResP = new TCanvas("cResP","",1200,900);
  hResPx->GetYaxis()->SetTitle("#sigma(p_{i}^{MC}-p_{i}^{reco}) (GeV/c)");
  hResPx->GetYaxis()->SetRangeUser(0.004,0.05);
  hResPx->Draw();
  hResPy->Draw("same");
  hResPz->Draw("same");
  lp->Draw("same");  

  TCanvas *cPullsPos = new TCanvas("cPullsPos","",1200,900);
  hPullsX->GetYaxis()->SetTitle("#sigma((pos_{MC}-pos_{reco})/#sigma_{pos_{reco}})");
  hPullsX->GetYaxis()->SetRangeUser(0.,5.);
  hPullsX->Draw();
  hPullsY->Draw("same");
  hPullsZ->Draw("same");
  l->Draw("same");
  
  TCanvas *cPullsP = new TCanvas("cPullsP","",1200,900);
  hPullsPx->GetYaxis()->SetTitle("#sigma((p_{i}^{MC}-p_{i}^{reco})/#sigma_{p_{i}^{reco}})");
  hPullsPx->GetYaxis()->SetRangeUser(0.,2.);
  hPullsPx->Draw();
  hPullsPy->Draw("same");
  hPullsPz->Draw("same");
  lp->Draw("same");

  cResPos->SaveAs("DauResPos.eps");
  cResP->SaveAs("DauResP.eps");
  cPullsPos->SaveAs("DauPullsPos.eps");
  cPullsP->SaveAs("DauPullsP.eps");
  
  TFile outfile2(outfilename.Data(),"UPDATE");
  hResX->Write();
  hResY->Write();
  hResZ->Write();
  hResPx->Write();
  hResPy->Write();
  hResPz->Write();
  hPullsX->Write();
  hPullsY->Write();
  hPullsZ->Write();
  hPullsPx->Write();
  hPullsPy->Write();
  hPullsPz->Write();
  cResPos->Write();
  cResP->Write();
  cPullsPos->Write();
  cPullsP->Write();
  outfile2.Close();

  delete cResPos;
  delete cResP;
  delete cPullsPos;
  delete cPullsP;
 
}

TH1F* GetHistoVsPt(THnSparseF* sparse, Int_t axis, TString outfilename, Bool_t noFIT, Bool_t singlegauss) {

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
    f = new TF1("f","gaus",-10,10);
    f->SetParameter(1,0.);
    f->SetParameter(2,1.);
  }
  TString fitoption="";
  if(!singlegauss)
    fitoption="RLEM";
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ResetAxes(sparse);
    SetPtRange(sparse,PtLims[iPt],PtLims[iPt+1]);

    TH1F* h = (TH1F*)sparse->Projection(axis);

    if(!noFIT) {
      if(!singlegauss) f->FixParameter(0,h->GetEntries()*h->GetBinWidth(2));
      if(h->GetEntries()!=0) {
        c->Clear();
        c->Update();
        h->Fit("f",fitoption.Data());
        Double_t chi = f->GetChisquare()/f->GetNDF();
        Bool_t badfit=kFALSE;
        if(!singlegauss && (f->GetParError(2)>f->GetParameter(2) || f->GetParError(3)>f->GetParameter(3))) {
          f2 = new TF1("f2","gaus",-10,10);
          f2->SetParameter(1,0.);
          f2->SetParameter(2,1.);
          h->Fit("f2");
          badfit=kTRUE;
        }
        if(!singlegauss && chi>3) {
          f2 = new TF1("f2","gaus",-10,10);
          f2->SetParameter(1,0.);
          f2->SetParameter(2,1.);
          h->Fit("f2");
          Double_t chi2 = f2->GetChisquare()/f2->GetNDF();
          if(chi2<chi)
            badfit=kTRUE;
          else
            h->Fit("f",fitoption.Data());
        }
        outfile.cd();
        h->Write();
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
      hVsPt->SetBinContent(iPt+1,h->GetRMS());
      hVsPt->SetBinError(iPt+1,h->GetRMSError()); 
    }
  }
  
  delete c;
  outfile.Close();
  
  hVsPt->GetXaxis()->SetTitle(sparse->GetAxis(sparse->GetNdimensions()-1)->GetTitle());
  hVsPt->GetYaxis()->SetTitle(sparse->GetAxis(axis)->GetTitle());
  
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
