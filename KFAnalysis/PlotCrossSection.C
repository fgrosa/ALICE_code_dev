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
#include <TArrayD.h>

#endif

const Int_t colors[] = {kRed,kBlack,kBlue,kOrange+7,kGreen+3,kMagenta,kCyan,kYellow+3,kOrange+3};

void PlotCrossSection(Int_t iSet=1,Double_t BR=0.0913) {
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetLegendBorderSize(0);
  
  TFile pubfile("/home/fabrizio/tesi/Risultati_pubblicati/DplusCrossSec_method2_fd2_br1.root","UPDATE");
  TH1F* hPub = (TH1F*)pubfile.Get("hAAC");
  hPub->SetDirectory(0);
  hPub->SetTitle("");
  hPub->SetLineColor(colors[1]);
  hPub->SetMarkerColor(colors[1]);
  hPub->SetMarkerStyle(20);
  hPub->SetLineWidth(2);
  hPub->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}} (#mub c/GeV)");
  hPub->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hPub->GetYaxis()->SetTitleSize(0.05);
  hPub->GetXaxis()->SetTitleSize(0.05);
  hPub->GetYaxis()->SetLabelSize(0.05);
  hPub->GetXaxis()->SetLabelSize(0.05);
  hPub->GetYaxis()->SetTitleOffset(1.5);
  hPub->GetXaxis()->SetTitleOffset(1.);
  TGraphAsymmErrors* gPub = (TGraphAsymmErrors*)pubfile.Get("gaaCsystTot");
  gPub->SetLineColor(colors[1]);
  gPub->SetLineWidth(2);
  pubfile.Close();

  TFile infile(Form("../PID/HFPtSpectrum_combinedFD_cutset%d.root",iSet),"UPDATE");
  TH1F* hCross = (TH1F*)infile.Get("histoSigmaCorr");
  hCross->SetDirectory(0);
  hCross->SetMarkerStyle(33);
  hCross->SetMarkerColor(colors[0]);
  hCross->SetLineColor(colors[0]);
  hCross->Scale(1./(1000000*BR));
  hCross->GetXaxis()->SetLabelFont(42);
  hCross->GetXaxis()->SetTitleFont(42);
  hCross->GetYaxis()->SetLabelFont(42);
  hCross->GetYaxis()->SetTitleFont(42);
  hCross->SetLineWidth(2);
  TGraphAsymmErrors* gFracKF = (TGraphAsymmErrors*)infile.Get("gFcCorrConservative");
  gFracKF->SetMarkerColor(colors[0]);
  gFracKF->SetLineWidth(2);
  gFracKF->SetMarkerSize(1.5);
  gFracKF->SetMarkerStyle(33);
  TGraphAsymmErrors *gFracSyst = (TGraphAsymmErrors*)infile.Get("gSigmaCorrConservativePC");
  gFracSyst->SetName("gFracSyst");
  gFracSyst->SetLineColor(colors[7]);
  gFracSyst->GetXaxis()->SetLabelFont(42);
  gFracSyst->GetXaxis()->SetTitleFont(42);
  gFracSyst->GetYaxis()->SetLabelFont(42);
  gFracSyst->GetYaxis()->SetTitleFont(42);
  infile.Close();
  
  TFile infileStd("HFPtSpectrum_combinedFD_standard.root","UPDATE");
  TGraphAsymmErrors* gFracPub = (TGraphAsymmErrors*)infileStd.Get("gFcCorrConservative");
  gFracPub->SetLineColor(colors[1]);
  gFracPub->SetMarkerColor(colors[1]);
  gFracPub->GetXaxis()->SetTitleFont(42);
  gFracPub->GetYaxis()->SetTitleFont(42);
  gFracPub->GetXaxis()->SetLabelFont(42);
  gFracPub->GetYaxis()->SetLabelFont(42);
  gFracPub->GetYaxis()->SetTitle("#it{f}_{prompt}");
  gFracPub->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  gFracPub->GetYaxis()->SetTitleSize(0.05);
  gFracPub->GetXaxis()->SetTitleSize(0.05);
  gFracPub->GetYaxis()->SetLabelSize(0.05);
  gFracPub->GetXaxis()->SetLabelSize(0.05);
  gFracPub->GetYaxis()->SetTitleOffset(1.3);
  gFracPub->GetXaxis()->SetTitleOffset(1.);
  gFracPub->SetLineWidth(2);
  gFracPub->SetMarkerSize(1.5);
  gFracPub->SetMarkerStyle(20);
  gFracPub->SetTitle("");
  TGraphAsymmErrors *gFracSystPub = (TGraphAsymmErrors*)infileStd.Get("gSigmaCorrConservativePC");
  infileStd.Close();

  const Int_t nPtBins=hCross->GetNbinsX();
  TArrayD* ptarray=(TArrayD*)hCross->GetXaxis()->GetXbins();
  Double_t *PtLims=(Double_t*)ptarray->GetArray();
  
  //systematics
  Double_t cutvarsystlow[nPtBins] = {0.11,0.06,0.04,0.04,0.04,0.04,0.05,0.05,0.05,0.06};
  Double_t cutvarsysthigh[nPtBins] = {0.11,0.06,0.04,0.04,0.04,0.04,0.05,0.05,0.05,0.06};
  Double_t rawsystlow[nPtBins] = {0.12,0.07,0.04,0.04,0.04,0.04,0.04,0.04,0.08,0.10};
  Double_t rawsysthigh[nPtBins] = {0.12,0.07,0.04,0.04,0.04,0.04,0.04,0.04,0.08,0.10};
  Double_t effsystlow[nPtBins] = {0.01,0.,0.,0.,0.,0.,0.,0.,0.01,0.02};
  Double_t effsysthigh[nPtBins] = {0.01,0.,0.,0.,0.,0.,0.,0.,0.01,0.02};
  Double_t tracksystlow[nPtBins] = {0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09};
  Double_t tracksysthigh[nPtBins] = {0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09};
  Double_t PIDsystlow[nPtBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t PIDsysthigh[nPtBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t errsystlow[nPtBins];
  Double_t errsysthigh[nPtBins];
  
  TGraphAsymmErrors *gCutVarSyst = new TGraphAsymmErrors(nPtBins);
  gCutVarSyst->SetLineColor(colors[2]);
  gCutVarSyst->SetName("gCutVarSyst");
  TGraphAsymmErrors *gRawSyst = new TGraphAsymmErrors(nPtBins);
  gRawSyst->SetLineColor(colors[3]);
  gRawSyst->SetName("gRawSyst");
  TGraphAsymmErrors *gEffSyst = new TGraphAsymmErrors(nPtBins);
  gEffSyst->SetLineColor(colors[4]);
  gEffSyst->SetName("gEffSyst");
  TGraphAsymmErrors *gTrackSyst = new TGraphAsymmErrors(nPtBins);
  gTrackSyst->SetLineColor(colors[5]);
  gTrackSyst->SetName("gTrackSyst");
  TGraphAsymmErrors *gPIDSyst = new TGraphAsymmErrors(nPtBins);
  gPIDSyst->SetLineColor(colors[6]);
  gPIDSyst->SetName("gPIDSyst");
  TGraphAsymmErrors *gTotSyst = new TGraphAsymmErrors(nPtBins);
  gTotSyst->SetLineColor(colors[0]);
  gTotSyst->SetName("gTotSyst");
  TGraphAsymmErrors *gCrossSyst = new TGraphAsymmErrors(nPtBins);
  gCrossSyst->SetLineColor(colors[0]);
  gCrossSyst->SetFillStyle(20);
  gCrossSyst->SetName("gCrossSyst");
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    Double_t ptcent = (PtLims[iPt+1]+PtLims[iPt])/2;
    Double_t step = (PtLims[iPt+1]-PtLims[iPt])/2;
    gCutVarSyst->SetPoint(iPt,ptcent,0);
    gCutVarSyst->SetPointError(iPt,step,step,cutvarsystlow[iPt],cutvarsysthigh[iPt]);
    gRawSyst->SetPoint(iPt,ptcent,0);
    gRawSyst->SetPointError(iPt,step,step,rawsystlow[iPt],rawsysthigh[iPt]);
    gEffSyst->SetPoint(iPt,ptcent,0);
    gEffSyst->SetPointError(iPt,step,step,effsystlow[iPt],effsysthigh[iPt]);
    gTrackSyst->SetPoint(iPt,ptcent,0);
    gTrackSyst->SetPointError(iPt,step,step,tracksystlow[iPt],tracksysthigh[iPt]);
    gPIDSyst->SetPoint(iPt,ptcent,0);
    gPIDSyst->SetPointError(iPt,step,step,PIDsystlow[iPt],PIDsysthigh[iPt]);
    errsystlow[iPt] = TMath::Sqrt(gCutVarSyst->GetErrorYlow(iPt)*gCutVarSyst->GetErrorYlow(iPt)+
                                  gRawSyst->GetErrorYlow(iPt)*gRawSyst->GetErrorYlow(iPt)+
                                  gEffSyst->GetErrorYlow(iPt)*gEffSyst->GetErrorYlow(iPt)+
                                  gTrackSyst->GetErrorYlow(iPt)*gTrackSyst->GetErrorYlow(iPt)+
                                  gPIDSyst->GetErrorYlow(iPt)*gPIDSyst->GetErrorYlow(iPt)+
                                  gFracSyst->GetErrorYlow(iPt+1)*gFracSyst->GetErrorYlow(iPt+1));
    errsysthigh[iPt] = TMath::Sqrt(gCutVarSyst->GetErrorYhigh(iPt)*gCutVarSyst->GetErrorYhigh(iPt)+
                                  gRawSyst->GetErrorYhigh(iPt)*gRawSyst->GetErrorYhigh(iPt)+
                                  gEffSyst->GetErrorYhigh(iPt)*gEffSyst->GetErrorYhigh(iPt)+
                                  gTrackSyst->GetErrorYhigh(iPt)*gTrackSyst->GetErrorYhigh(iPt)+
                                  gPIDSyst->GetErrorYhigh(iPt)*gPIDSyst->GetErrorYhigh(iPt)+
                                  gFracSyst->GetErrorYhigh(iPt+1)*gFracSyst->GetErrorYhigh(iPt+1));
    gTotSyst->SetPoint(iPt,ptcent,0);
    gTotSyst->SetPointError(iPt,step,step,errsystlow[iPt],errsysthigh[iPt]);
    gCrossSyst->SetPoint(iPt,ptcent,hCross->GetBinContent(iPt+1));
    gCrossSyst->SetPointError(iPt,0.15,0.15,errsystlow[iPt]*hCross->GetBinContent(iPt+1),errsysthigh[iPt]*hCross->GetBinContent(iPt+1));
  }
  gCrossSyst->SetLineWidth(2);

  TH1F *hRatio = (TH1F*)hCross->Clone();
  hRatio->Divide(hCross,hPub,1.,1.,"");
  hRatio->SetDirectory(0);
  hRatio->SetTitle("");
  hRatio->GetXaxis()->SetTitleFont(42);
  hRatio->GetYaxis()->SetTitleFont(42);
  hRatio->GetXaxis()->SetLabelFont(42);
  hRatio->GetYaxis()->SetLabelFont(42);
  hRatio->GetXaxis()->SetTitleSize(0.05);
  hRatio->GetYaxis()->SetTitleSize(0.05);
  hRatio->GetYaxis()->SetTitle("ratio of #frac{d#sigma}{d#it{p}_{T}} w.r.t. the published result");
  hRatio->GetXaxis()->SetTitle(hPub->GetXaxis()->GetTitle());

  //systematic uncertainty on the ratio
  TGraphAsymmErrors* gSystRatio = new TGraphAsymmErrors(nPtBins);
  gSystRatio->SetLineColor(colors[0]);
  gSystRatio->SetFillStyle(20);
  gSystRatio->SetLineWidth(2);
  gSystRatio->SetName("gSystRatio");

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    Double_t ptcent = (PtLims[iPt]+PtLims[iPt+1])/2;

    //error on fprompt and track efficiency completely correlated
    Double_t systerrnotracknofpromptlow = TMath::Sqrt(gCutVarSyst->GetErrorYlow(iPt)*gCutVarSyst->GetErrorYlow(iPt)+
                                                      gRawSyst->GetErrorYlow(iPt)*gRawSyst->GetErrorYlow(iPt)+
                                                      gEffSyst->GetErrorYlow(iPt)*gEffSyst->GetErrorYlow(iPt)+
                                                      gPIDSyst->GetErrorYlow(iPt)*gPIDSyst->GetErrorYlow(iPt));
    Double_t systerrnotracknofprompthigh = TMath::Sqrt(gCutVarSyst->GetErrorYhigh(iPt)*gCutVarSyst->GetErrorYhigh(iPt)+
                                                       gRawSyst->GetErrorYhigh(iPt)*gRawSyst->GetErrorYhigh(iPt)+
                                                       gEffSyst->GetErrorYhigh(iPt)*gEffSyst->GetErrorYhigh(iPt)+
                                                       gPIDSyst->GetErrorYhigh(iPt)*gPIDSyst->GetErrorYhigh(iPt));
    Double_t xpub,ypub;
    gPub->GetPoint(iPt,xpub,ypub);

    Double_t puberrlow = TMath::Sqrt(gPub->GetErrorYlow(iPt)*gPub->GetErrorYlow(iPt)/(ypub*ypub)-0.09*0.09-gFracSystPub->GetErrorYlow(iPt+1)*gFracSystPub->GetErrorYlow(iPt+1));
    Double_t puberrhigh = TMath::Sqrt(gPub->GetErrorYhigh(iPt)*gPub->GetErrorYhigh(iPt)/(ypub*ypub)-0.09*0.09-gFracSystPub->GetErrorYhigh(iPt+1)*gFracSystPub->GetErrorYhigh(iPt+1));

    //ratio of dsdpt+sigma(fprompt) and dsdpt-sigma(fprompt) 
    Double_t ratioupedgesfprompt = (hCross->GetBinContent(iPt+1)*(1+gFracSyst->GetErrorYhigh(iPt+1)))/(hPub->GetBinContent(iPt+1)*(1+gFracSystPub->GetErrorYhigh(iPt+1)));
    Double_t ratiolowedgesfprompt = (hCross->GetBinContent(iPt+1)*(1-gFracSyst->GetErrorYlow(iPt+1)))/(hPub->GetBinContent(iPt+1)*(1-gFracSystPub->GetErrorYlow(iPt+1)));

    Double_t uperrfrac = ratioupedgesfprompt-hRatio->GetBinContent(iPt+1);
    Double_t lowerrfrac = hRatio->GetBinContent(iPt+1)-ratiolowedgesfprompt;
    Double_t uperrfrac2=0;
    Double_t lowerrfrac2=0;
    
    if(uperrfrac<0 && lowerrfrac>0) {
      lowerrfrac2=-uperrfrac;
      uperrfrac=0;
    }
    if(lowerrfrac<0 && uperrfrac>0) {
      uperrfrac2=-lowerrfrac;
      lowerrfrac=0;
    }
    if(lowerrfrac<0 && uperrfrac<0) {
      uperrfrac2=-lowerrfrac;
      lowerrfrac=-uperrfrac;
      uperrfrac=uperrfrac2;
      uperrfrac2=0;
    }

    if(lowerrfrac2>lowerrfrac) {
      lowerrfrac=lowerrfrac2;
    }
    if(uperrfrac2>uperrfrac) {
      uperrfrac=uperrfrac2;
    }
        
    //ratio of dsdpt+sigma(track) and dsdpt-sigma(track) 
    Double_t ratioupedgestrack = (hCross->GetBinContent(iPt+1)*(1+0.09))/(hPub->GetBinContent(iPt+1)*(1+0.09));
    Double_t ratiolowedgestrack = (hCross->GetBinContent(iPt+1)*(1-0.09))/(hPub->GetBinContent(iPt+1)*(1-0.09));
   
    gSystRatio->SetPoint(iPt,ptcent,hRatio->GetBinContent(iPt+1));
    gSystRatio->SetPointError(iPt,0.15,0.15,TMath::Sqrt(systerrnotracknofpromptlow*systerrnotracknofpromptlow+lowerrfrac*lowerrfrac+puberrlow*puberrlow),TMath::Sqrt(systerrnotracknofprompthigh*systerrnotracknofprompthigh+uperrfrac*uperrfrac+puberrhigh*puberrhigh));

    
  }
  
  TLegend* l = new TLegend(0.35,0.73,0.8,0.85);
  l->SetTextSize(0.045);
  TLegend* l2 = (TLegend*)l->Clone();
  l2->SetY1(0.25);
  l2->SetY2(0.5);
  l2->SetX1(0.55);
  l2->SetX2(0.85);
  l->AddEntry(hCross,"KFParticle","lpe");
  l->AddEntry(hPub,"Published Cross Section","lpe");
  l2->AddEntry(gFracKF,"KFParticle","lpe");
  l2->AddEntry(gFracPub,"Published","lpe");

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.045);
  latex->SetTextFont(132);

  TLine* line = new TLine(1,1,24,1);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->SetLineColor(colors[1]);
  
  TCanvas *c = new TCanvas("c","",900,1000);
  c->cd();
  TPad* pad1 = new TPad("pad1","",0,0.3,1.,1.);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->SetLogy();
  pad1->cd();
  hPub->GetYaxis()->SetRangeUser(0.2,50000);
  hPub->Draw();
  gPub->Draw("2");
  hCross->Draw("same");
  gCrossSyst->Draw("2");
  l->Draw("same");
  latex->DrawLatex(PtLims[2]-0.5,7.5,"Prompt D^{+}");
  latex->DrawLatex(PtLims[2]-0.5,3,"pPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->DrawLatex(PtLims[2]-0.5,1.2,"#pm 2.1% BR unc. not shown");
  latex->DrawLatex(PtLims[2]-0.5,0.5,"#pm 3.7% norm. unc. not shown");
  c->cd();
  TPad* pad2 = new TPad("pad2","",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();
  pad2->cd();
  TH1F* hRatioClone=(TH1F*)hRatio->Clone();
  hRatioClone->GetYaxis()->SetRangeUser(0.45,1.45);
  hRatioClone->GetYaxis()->SetTitle("#frac{KFParticle}{Published}   ");
  hRatioClone->GetYaxis()->SetNdivisions(505);
  hRatioClone->GetYaxis()->SetTitleOffset(0.7);
  hRatioClone->GetYaxis()->SetTitleSize(0.11);
  hRatioClone->GetXaxis()->SetTitleSize(0.11);
  hRatioClone->GetYaxis()->SetLabelSize(0.10);
  hRatioClone->GetXaxis()->SetLabelSize(0.10);
  hRatioClone->Draw();
  gSystRatio->Draw("2");
  line->Draw("same");
  
  TCanvas *cFrac = new TCanvas("cFrac","",800,600);
  gFracPub->GetYaxis()->SetRangeUser(0.5,1.);
  gFracPub->Draw("AP");
  gFracKF->Draw("P");
  l2->Draw("same");
  
  TCanvas *cRatio = new TCanvas("cRatio","",800,800);
  hRatio->GetYaxis()->SetRangeUser(0.45,1.5);
  hRatio->Draw();
  gSystRatio->Draw("2");

  c->SaveAs("CrossSection_KF.eps");
  cFrac->SaveAs("PromptFrac_KF.eps");
  cRatio->SaveAs("Ratio_KF_pub.eps");

  TFile outfile("DplusCrossSectionKF_pPb5TeV.root","RECREATE");
  gCutVarSyst->Write();
  gRawSyst->Write();
  gEffSyst->Write();
  gTrackSyst->Write();
  gPIDSyst->Write();
  gFracSyst->Write();
  gTotSyst->Write();

  hCross->Write();
  gCrossSyst->Write();
  outfile.Close();
  

  for(Int_t iPt=1; iPt<=nPtBins; iPt++) {

    Double_t fKF;
    Double_t fPub;
    Double_t pt;
    gFracKF->GetPoint(iPt,pt,fKF);
    gFracPub->GetPoint(iPt,pt,fPub);
    
    cout << "KFParticle: " << fKF << " + " << gFracSyst->GetErrorYhigh(iPt)<< " - " <<gFracSyst->GetErrorYlow(iPt) << "     Published: " << fPub << " + " << gFracSystPub->GetErrorYhigh(iPt)<< " - " <<gFracSystPub->GetErrorYlow(iPt)<<endl;
  }
  
}
