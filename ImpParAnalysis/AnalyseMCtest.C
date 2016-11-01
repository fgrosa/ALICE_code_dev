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
#include <TRandom3.h>
#include <TLine.h>
#include <TDatabasePDG.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>

#endif

//______________________________________________________________________________________________________
//global variables
const Int_t nFracs = 1;
TString FracValues[nFracs] = {"0.95"};
const Int_t colors[] = {kBlack};
const Int_t markers[] = {20};
TString infilename = Form("PromptFraction_TOYMC_%s_bkg_unbinned_sigmafree.root",FracValues[0].Data()); 

//______________________________________________________________________________________________________
//function prototypes
void DrawPromptFractionResult(Bool_t d0cutapplied=kTRUE);
void DrawSigmaPromptResult();

//______________________________________________________________________________________________________
void DrawPromptFractionResult(Bool_t d0cutapplied) {
  
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.12);

  TH1F** hBias = new TH1F*[nFracs];
  TH1F** hRMS = new TH1F*[nFracs];
  TH1F** hPulls = new TH1F*[nFracs];
  
  for(Int_t iGenFrac=0; iGenFrac<nFracs; iGenFrac++) {
    infilename.ReplaceAll("SigmaPrompt","PromptFraction");
    if(FracValues[iGenFrac]=="#it{f}_{N_{b}}")
      FracValues[iGenFrac] = "Nbfrac";
    if(iGenFrac>0)
      infilename.ReplaceAll(FracValues[iGenFrac-1].Data(),FracValues[iGenFrac].Data());
    TFile infile(infilename.Data(),"READ");
    hBias[iGenFrac]=(TH1F*)infile.Get("hBiasVsPt");
    hRMS[iGenFrac]=(TH1F*)infile.Get("hResVsPt");
    hPulls[iGenFrac]=(TH1F*)infile.Get("hPullsVsPt");
    hBias[iGenFrac]->SetDirectory(0);
    hRMS[iGenFrac]->SetDirectory(0);
    hPulls[iGenFrac]->SetDirectory(0);
    infile.Close();
    if(d0cutapplied) {
      TString title = hBias[iGenFrac]->GetYaxis()->GetTitle();
      hBias[iGenFrac]->GetYaxis()->SetTitle(title.ReplaceAll(">","> [-d_{0}^{cut},+d_{0}^{cut}]"));
      title = hRMS[iGenFrac]->GetYaxis()->GetTitle();
      hRMS[iGenFrac]->GetYaxis()->SetTitle(title.ReplaceAll(")",") [-d_{0}^{cut},+d_{0}^{cut}]"));
      title = hPulls[iGenFrac]->GetYaxis()->GetTitle();
      hPulls[iGenFrac]->GetYaxis()->SetTitle(title.ReplaceAll("#sigma_{#it{f}_{prompt}^{meas}})","#sigma_{#it{f}_{prompt}^{meas}}) [-d_{0}^{cut},+d_{0}^{cut}]"));
    }
  }
  
  TLine* line = new TLine(2,0,16,0);
  line->SetLineWidth(2);
  line->SetLineColor(kBlue);
  line->SetLineStyle(7);
  
  TLegend *l = new TLegend(0.5,0.6,0.87,0.89);
  for(Int_t iGenFrac=0; iGenFrac<nFracs; iGenFrac++) {
    if(FracValues[iGenFrac]=="Nbfrac")
      FracValues[iGenFrac] = "#it{f}_{N_{b}}";
    
    l->AddEntry(hBias[iGenFrac],Form("#it{f}^{gen}_{prompt} = %s",FracValues[iGenFrac].Data()),"lpe");
  }
  l->SetTextSize(0.045);
  l->SetBorderSize(0);
  l->SetFillStyle(0);

  for(Int_t iGenFrac=0; iGenFrac<nFracs; iGenFrac++) {
    hBias[iGenFrac]->SetStats(0);
    hBias[iGenFrac]->SetLineColor(colors[iGenFrac]);
    hBias[iGenFrac]->SetMarkerColor(colors[iGenFrac]);
    hBias[iGenFrac]->SetMarkerStyle(markers[iGenFrac]);
    hBias[iGenFrac]->GetYaxis()->SetRangeUser(-0.08,0.2);
    hBias[iGenFrac]->GetYaxis()->SetTitleOffset(1.7);
    hRMS[iGenFrac]->SetStats(0);
    hRMS[iGenFrac]->SetLineColor(colors[iGenFrac]);
    hRMS[iGenFrac]->SetMarkerColor(colors[iGenFrac]);
    hRMS[iGenFrac]->SetMarkerStyle(markers[iGenFrac]);
    hRMS[iGenFrac]->GetYaxis()->SetRangeUser(0.,0.2);
    hRMS[iGenFrac]->GetYaxis()->SetTitleOffset(1.7);
    hPulls[iGenFrac]->SetStats(0);
    hPulls[iGenFrac]->SetLineColor(colors[iGenFrac]);
    hPulls[iGenFrac]->SetMarkerColor(colors[iGenFrac]);
    hPulls[iGenFrac]->SetMarkerStyle(markers[iGenFrac]);
    hPulls[iGenFrac]->GetYaxis()->SetRangeUser(0.,2.);
    hPulls[iGenFrac]->GetYaxis()->SetTitleOffset(1.7);
  } 

  TCanvas* cBias= new TCanvas("cBias","cBias",800,800);
  cBias->Clear();
  hBias[0]->Draw("E1");
  if(nFracs>1) {
    for(Int_t iGenFrac=1; iGenFrac<nFracs; iGenFrac++) {
      hBias[iGenFrac]->Draw("E1same");
    }
  }
  line->Draw("same");
  l->Draw("same");

  TCanvas* cRMS= new TCanvas("cRMS","cRMS",800,800);
  cRMS->Clear();
  hRMS[0]->Draw("E1");
  if(nFracs>1) {
    for(Int_t iGenFrac=1; iGenFrac<nFracs; iGenFrac++) {
      hRMS[iGenFrac]->Draw("E1same");
    }
  }
  l->Draw("same");

  TCanvas* cPulls= new TCanvas("cPulls","cPulls",800,800);
  cPulls->Clear();
  hPulls[0]->Draw("E1");
  if(nFracs>1) {
    for(Int_t iGenFrac=1; iGenFrac<nFracs; iGenFrac++) {
      hPulls[iGenFrac]->Draw("E1same");
    }
  }
  l->Draw("same");

  TString bkg="bkg";
  if(!infilename.Contains("_bkg_")) bkg="nobkg";
  TString sigmapar="sigmafree";
  if(!infilename.Contains("sigmafree")) sigmapar="sigmafixed";
  TString method="unbinned";
  if(!infilename.Contains("unbinned")) method="binned";
  
  cBias->SaveAs(Form("PromptFrac_Bias_%s_%s_%s.pdf",bkg.Data(),method.Data(),sigmapar.Data()));
  cRMS->SaveAs(Form("PromptFrac_RMS_%s_%s_%s.pdf",bkg.Data(),method.Data(),sigmapar.Data()));
  cPulls->SaveAs(Form("PromptFrac_Pulls_%s_%s_%s.pdf",bkg.Data(),method.Data(),sigmapar.Data()));

  for(Int_t iGenFrac=0; iGenFrac<nFracs; iGenFrac++) {
    delete hBias[iGenFrac];
    delete hRMS[iGenFrac];
    delete hPulls[iGenFrac];
  }
  delete[] hBias;
  delete[] hPulls;
  delete[] hRMS;
  delete cBias;
  delete cRMS;
  delete cPulls;
  
  for(Int_t iGenFrac=0; iGenFrac<nFracs; iGenFrac++) {
    if(FracValues[iGenFrac]=="#it{f}_{N_{b}}")
      FracValues[iGenFrac] = "Nbfrac";
  } 
  infilename.ReplaceAll(FracValues[nFracs-1].Data(),FracValues[0].Data());
}

//______________________________________________________________________________________________________
void DrawSigmaPromptResult() {
  
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.12);

  TH1F** hSigma = new TH1F*[nFracs];
  TH1F* hSigmaTrue = 0x0;
  TH1F** hSigmaRMS = new TH1F*[nFracs];
  TH1F** hSigmaPulls = new TH1F*[nFracs];
  
  for(Int_t iGenFrac=0; iGenFrac<nFracs; iGenFrac++) {
    infilename.ReplaceAll("PromptFraction","SigmaPrompt");
    if(iGenFrac>0)
      infilename.ReplaceAll(FracValues[iGenFrac-1].Data(),FracValues[iGenFrac].Data());
    TFile infile(infilename.Data(),"READ");
    hSigma[iGenFrac]=(TH1F*)infile.Get("hSigmaVsPt");
    hSigmaRMS[iGenFrac]=(TH1F*)infile.Get("hSigmaResVsPt");
    hSigmaPulls[iGenFrac]=(TH1F*)infile.Get("hSigmaPullsVsPt");
    hSigma[iGenFrac]->SetDirectory(0);
    hSigmaRMS[iGenFrac]->SetDirectory(0);
    hSigmaPulls[iGenFrac]->SetDirectory(0);
    if(iGenFrac==0) {
      hSigmaTrue=(TH1F*)infile.Get("hSigmaTrueVsPt");
      hSigmaTrue->SetDirectory(0);
    }
    infile.Close();
  }
  
  TLegend *l = new TLegend(0.5,0.6,0.87,0.89);
  for(Int_t iGenFrac=0; iGenFrac<nFracs; iGenFrac++) {
    if(FracValues[iGenFrac]=="Nbfrac")
      FracValues[iGenFrac] = "#it{f}_{N_{b}}";
    l->AddEntry(hSigma[iGenFrac],Form("#it{f}^{gen}_{prompt} = %s",FracValues[iGenFrac].Data()),"lpe");
  }
  l->SetTextSize(0.045);
  l->SetBorderSize(0);
  l->SetFillStyle(0);

  TLegend *l2 = new TLegend(0.5,0.5,0.87,0.89);
  l2->SetTextSize(0.045);
  l2->SetBorderSize(0);
  l2->SetFillStyle(0);
  for(Int_t iGenFrac=0; iGenFrac<nFracs; iGenFrac++) {
    l2->AddEntry(hSigma[iGenFrac],Form("#it{f}^{gen}_{prompt} = %s",FracValues[iGenFrac].Data()),"lpe");
  }
  l2->AddEntry(hSigmaTrue,"MC prefit","lpe");
  
  for(Int_t iGenFrac=0; iGenFrac<nFracs; iGenFrac++) {
    hSigma[iGenFrac]->SetStats(0);
    hSigma[iGenFrac]->SetLineColor(colors[iGenFrac]);
    hSigma[iGenFrac]->SetMarkerColor(colors[iGenFrac]);
    hSigma[iGenFrac]->SetMarkerStyle(markers[iGenFrac]);
    hSigma[iGenFrac]->GetYaxis()->SetRangeUser(20.,55.);
    hSigma[iGenFrac]->GetYaxis()->SetTitleOffset(1.7);
    hSigmaRMS[iGenFrac]->SetStats(0);
    hSigmaRMS[iGenFrac]->SetLineColor(colors[iGenFrac]);
    hSigmaRMS[iGenFrac]->SetMarkerColor(colors[iGenFrac]);
    hSigmaRMS[iGenFrac]->SetMarkerStyle(markers[iGenFrac]);
    hSigmaRMS[iGenFrac]->GetYaxis()->SetRangeUser(0.,5);
    hSigmaRMS[iGenFrac]->GetYaxis()->SetTitleOffset(1.7);
    hSigmaPulls[iGenFrac]->SetStats(0);
    hSigmaPulls[iGenFrac]->SetLineColor(colors[iGenFrac]);
    hSigmaPulls[iGenFrac]->SetMarkerColor(colors[iGenFrac]);
    hSigmaPulls[iGenFrac]->SetMarkerStyle(markers[iGenFrac]);
    hSigmaPulls[iGenFrac]->GetYaxis()->SetRangeUser(0.,2.);
    hSigmaPulls[iGenFrac]->GetYaxis()->SetTitleOffset(1.7);
  }

  hSigmaTrue->SetStats(0);
  hSigmaTrue->SetMarkerStyle(23);
  hSigmaTrue->GetYaxis()->SetTitle("< #sigma_{prompt} > (#mum)");
  hSigmaTrue->GetYaxis()->SetTitleOffset(1.7);
  hSigmaTrue->GetYaxis()->SetRangeUser(20.,55.);

  TCanvas* cSigma= new TCanvas("cSigma","cSigma",800,800);
  cSigma->Clear();
  hSigmaTrue->Draw("E1");
  for(Int_t iGenFrac=0; iGenFrac<nFracs; iGenFrac++) {
    hSigma[iGenFrac]->Draw("E1same");
  }
  l2->Draw("same");

  TCanvas* cSigmaRMS= new TCanvas("cSigmaRMS","cSigmaRMS",800,800);
  cSigmaRMS->Clear();
  hSigmaRMS[0]->Draw("E1");
  if(nFracs>1) {
    for(Int_t iGenFrac=1; iGenFrac<nFracs; iGenFrac++) {
      hSigmaRMS[iGenFrac]->Draw("E1same");
    }
  }
  l->Draw("same");

  TCanvas* cSigmaPulls= new TCanvas("cSigmaPulls","cSigmaPulls",800,800);
  cSigmaPulls->Clear();
  hSigmaPulls[0]->Draw("E1");
  if(nFracs>1) {
    for(Int_t iGenFrac=1; iGenFrac<nFracs; iGenFrac++) {
      hSigmaPulls[iGenFrac]->Draw("E1same");
    }
  }
  l->Draw("same");

  TString bkg="bkg";
  if(!infilename.Contains("_bkg_")) bkg="nobkg";
  TString sigmapar="sigmafree";
  if(!infilename.Contains("sigmafree")) sigmapar="sigmafixed";
  TString method="unbinned";
  if(!infilename.Contains("unbinned")) method="binned";
  
  cSigma->SaveAs(Form("SigmaPrompt_%s_%s_%s.pdf",bkg.Data(),method.Data(),sigmapar.Data()));
  cSigmaRMS->SaveAs(Form("SigmaPrompt_RMS_%s_%s_%s.pdf",bkg.Data(),method.Data(),sigmapar.Data()));
  cSigmaPulls->SaveAs(Form("SigmaPrompt_Pulls_%s_%s_%s.pdf",bkg.Data(),method.Data(),sigmapar.Data()));

  for(Int_t iGenFrac=0; iGenFrac<nFracs; iGenFrac++) {
    delete hSigma[iGenFrac];
    delete hSigmaRMS[iGenFrac];
    delete hSigmaPulls[iGenFrac];
  }
  delete hSigmaTrue;
  delete[] hSigma;
  delete[] hSigmaPulls;
  delete[] hSigmaRMS;
  delete cSigma;
  delete cSigmaRMS;
  delete cSigmaPulls;

 
  for(Int_t iGenFrac=0; iGenFrac<nFracs; iGenFrac++) {
    if(FracValues[iGenFrac]=="#it{f}_{N_{b}}")
      FracValues[iGenFrac] = "Nbfrac";
  }
  infilename.ReplaceAll(FracValues[nFracs-1].Data(),FracValues[0].Data());
}

