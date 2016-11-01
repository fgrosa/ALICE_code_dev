#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <vector>
#include <TInterpreter.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <THnSparse.h>
#include <TPad.h>
#include <TString.h>
#include <TLine.h>
#include <TList.h>
#include <TLatex.h>

#endif

const Int_t nPtBins = 7;
const Int_t nPtLims = nPtBins+1;
Double_t PtLims[nPtLims] = {2,3,4,5,6,8,12,16};

void PlotCxyDist();
void PlotLxyDist();
void PlotNLxyDist();
void PlotD0D0expDist();
void PlotDist(Int_t iPt=3, Double_t PtMin=5, Double_t PtMax=6, Int_t iVar=4, Int_t nProj = 6, TString VarName="d0d0exp", TString labelname="Max norm d_{0}-d_{0}^{exp}", Double_t legX1 = 0.45, Double_t legX2 = 0.85, Double_t legY1 = 0.7, Double_t legY2 = 0.89);

void PlotCxyDist() {
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    PlotDist(iPt,PtLims[iPt],PtLims[iPt+1],1,7,"cxy","Cos(#theta_{P}^{XY})",0.2,0.55);
  }
}

void PlotLxyDist() {
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    PlotDist(iPt,PtLims[iPt],PtLims[iPt+1],2,9,"Lxy","Decay length XY (cm)");
  }
}

void PlotNLxyDist() {
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    PlotDist(iPt,PtLims[iPt],PtLims[iPt+1],3,10,"NLxy","Normalised decay length XY");
  }
}

void PlotD0D0expDist() {
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    PlotDist(iPt,PtLims[iPt],PtLims[iPt+1],4,11,"d0d0exp","Max. norm. d_{0}-d_{0}^{exp}",0.55);
  }
}

void PlotDist(Int_t iPt, Double_t PtMin, Double_t PtMax, Int_t iVar, Int_t nProj, TString VarName, TString labelname, Double_t legX1, Double_t legX2, Double_t legY1, Double_t legY2) {
    
  gStyle->SetOptStat(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.16);  

  TFile infile(Form("CutVarDist_%d.root",iPt),"READ");
  TCanvas* c=(TCanvas*)infile.Get("cDist");
  TList *list = (TList*)c->GetListOfPrimitives();
  TH1F* hPrompt=(TH1F*)c->GetPad(iVar)->GetPrimitive(Form("hMassPtImpParPrompt_proj_%d",nProj));
  TH1F* hFD=(TH1F*)c->GetPad(iVar)->GetPrimitive(Form("hMassPtImpParBfeed_proj_%d",nProj));
  hPrompt->SetDirectory(0);
  hFD->SetDirectory(0);
  infile.Close();
  cout << "ciao" << endl;
  
  TLegend* l1 = new TLegend(legX1,legY1,legX2,legY2);
  l1->AddEntry(hPrompt,"Prompt","l");
  l1->AddEntry(hFD,"Feed-down","l");
  l1->SetTextSize(0.05);
  l1->SetFillStyle(0);
  
  TLatex* lat = new TLatex();
  lat->SetTextFont(132);
  lat->SetTextSize(0.055);
  
  TCanvas* cDist = new TCanvas("cDist","",800,800);
  cDist->Clear();
  cDist->SetLogy();
  hPrompt->GetYaxis()->SetTitleSize(0.05);
  hPrompt->GetXaxis()->SetTitleSize(0.05);
  hPrompt->GetYaxis()->SetTitleOffset(1.4);
  hPrompt->GetXaxis()->SetTitleOffset(1.2);
  hPrompt->GetYaxis()->SetLabelSize(0.05);
  hPrompt->SetLineWidth(2);
  hFD->SetLineWidth(2);  
  hFD->SetLineStyle(9);
  Double_t min = hPrompt->GetMinimum();
  Double_t max = hPrompt->GetMaximum();
  if(hFD->GetMinimum()<hPrompt->GetMinimum()) min = hFD->GetMinimum();
  if(hFD->GetMaximum()>hPrompt->GetMaximum()) max = hFD->GetMaximum();
  if(min==0) min=1.e-4;
  if(max>1) max=1;
  Double_t xlat = hPrompt->GetBinLowEdge(8);
  Double_t ylat = min*3;
  hPrompt->GetXaxis()->SetTitle(labelname.Data());
  if(VarName=="NLxy" || (VarName=="Lxy" && PtMin>5)) {
    if(min==1.e-4) min=1.e-5;
    hPrompt->GetYaxis()->SetRangeUser(min*10,max);
    ylat = min*30;
  }
  if(VarName=="Lxy" && PtMin<=5) {
    if(min==0) min=1.e-5;
    hPrompt->GetYaxis()->SetRangeUser(min,max);
    xlat = 0.5;
    ylat=max/40;
  }
  if(VarName=="cxy") {
    if(max==1) max=1.5;
    hPrompt->GetYaxis()->SetRangeUser(min,max);
    xlat=0.972;
    ylat=max/40;
  }
  if(VarName=="d0d0exp") {
    if(max==1) max=3.5;
    hPrompt->GetYaxis()->SetRangeUser(min,max);
    xlat=hFD->GetBinLowEdge(3);
    ylat=max/4;
  }
  hPrompt->GetXaxis()->SetNdivisions(508);
  hPrompt->GetXaxis()->SetLabelSize(0.05);
  hPrompt->SetTitle("");
  hPrompt->Draw();
  hFD->Draw("same");
  l1->Draw("same");      
  lat->DrawLatex(xlat,ylat,Form("%0.f < #it{p}_{T} < %0.f GeV/c",PtMin,PtMax));
  cDist->SaveAs(Form("Dist_%s_%0.f-%0.f.eps",VarName.Data(),PtMin,PtMax));
  
  delete cDist;

}
