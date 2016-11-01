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
#include <TPad.h>
#include <TGaxis.h>

#endif

const Int_t colors[] = {kRed,kBlack,kBlue,kGreen,kOrange+7,kGreen+3,kMagenta+1,kYellow+2,kCyan+2,kOrange+3};

void CutVarSyst(Int_t var=1, Double_t BR=.0913);
void ExtimateUncertainty();

void CutVarSyst(Int_t var, Double_t BR) {

  Int_t nSets=10;
  if(var==3)
    nSets=9;
  if(var==2)
    nSets=7;
  if(var==4 || var==5)
    nSets=7;
  
  Int_t *setsno = new Int_t[nSets];
  TString *legendnames = new TString[nSets];
  TString *plotnames = new TString[11];
  setsno[0]=1; ///reference set
  legendnames[0]="cental set";
  if(var==1) { ///chi square variations
    setsno[1]=2;
    setsno[2]=3;
    setsno[3]=4;
    setsno[4]=5;
    setsno[5]=28;
    setsno[6]=6;
    setsno[7]=29;
    setsno[8]=7;
    setsno[9]=8;
    legendnames[1]="#chi^{2} < central val + 1.50";
    legendnames[2]="#chi^{2} < central val + 1.00";
    legendnames[3]="#chi^{2} < central val + 0.50";
    legendnames[4]="#chi^{2} < central val - 0.50";
    legendnames[5]="#chi^{2} < central val - 0.75";
    legendnames[6]="#chi^{2} < central val - 1.00";
    legendnames[7]="#chi^{2} < central val - 1.25";
    legendnames[8]="#chi^{2} < central val - 1.50";
    legendnames[9]="central set w/o #chi^{2} cut";
    plotnames[0]="KF_CutVarSyst_chiS.eps";
    plotnames[1]="KF_CutVarSyst_sigmaonly_chiS.eps";
    plotnames[2]="KF_CutVarSyst_ratioonly_chiS.eps";    
    plotnames[3]="KF_CutVarSyst_effprompt_chiS.eps";
    plotnames[4]="KF_CutVarSyst_effFD_chiS.eps";
    plotnames[5]="KF_CutVarSyst_rawyield_chiS.eps";
    plotnames[6]="KF_CutVarSyst_sigma_chiS.eps";
    plotnames[7]="KF_CutVarSyst_ratioeffprompt_chiS.eps";
    plotnames[8]="KF_CutVarSyst_ratioeffFD_chiS.eps";
    plotnames[9]="KF_CutVarSyst_ratioraw_chiS.eps";
    plotnames[10]="KF_CutVarSyst_promptfrac_chiS.eps";
  }
  else if(var==2) { ///PV chi square variations
    setsno[0]=36;
    setsno[1]=13;
    setsno[2]=12;
    setsno[3]=9;
    setsno[4]=10;
    setsno[5]=11;
    setsno[6]=8;
    legendnames[0]="PV #chi^{2} < 3.0 w/o KF #chi^{2} cut";
    legendnames[1]="PV #chi^{2} < 4.0 w/o KF #chi^{2} cut";
    legendnames[2]="PV #chi^{2} < 3.5 w/o KF #chi^{2} cut";
    legendnames[3]="PV #chi^{2} < 2.5 w/o KF #chi^{2} cut";
    legendnames[4]="PV #chi^{2} < 2.0 w/o KF #chi^{2} cut";
    legendnames[5]="PV #chi^{2} < 1.5 w/o KF #chi^{2} cut";      
    legendnames[6]="central set w/o #chi^{2} cut";      
    plotnames[0]="KF_CutVarSyst_PVchi.eps";
    plotnames[1]="KF_CutVarSyst_sigmaonly_PVchi.eps";
    plotnames[2]="KF_CutVarSyst_ratioonly_PVchi.eps";
    plotnames[3]="KF_CutVarSyst_effprompt_PVchi.eps";
    plotnames[4]="KF_CutVarSyst_effFD_PVchi.eps";
    plotnames[5]="KF_CutVarSyst_rawyield_PVchi.eps";
    plotnames[6]="KF_CutVarSyst_sigma_PVchi.eps";
    plotnames[7]="KF_CutVarSyst_ratioeffprompt_PVchi.eps";
    plotnames[8]="KF_CutVarSyst_ratioeffFD_PVchi.eps";
    plotnames[9]="KF_CutVarSyst_ratioraw_PVchi.eps";
    plotnames[10]="KF_CutVarSyst_promptfrac_PVchi.eps";
  }
  else if(var==3) { ///norm DecL XY variations
    setsno[1]=31;
    setsno[2]=19;
    setsno[3]=18;
    setsno[4]=17;
    setsno[5]=14;
    setsno[6]=15;
    setsno[7]=16;
    setsno[8]=31;
    legendnames[1]="nDecLXY > central val - 2.0";
    legendnames[2]="nDecLXY > central val - 1.5";      
    legendnames[3]="nDecLXY > central val - 1.0";      
    legendnames[4]="nDecLXY > central val - 0.5";      
    legendnames[5]="nDecLXY > central val + 0.5";
    legendnames[6]="nDecLXY > central val + 1.0";
    legendnames[7]="nDecLXY > central val + 1.5";
    legendnames[8]="nDecLXY > central val + 2.0";
    plotnames[0]="KF_CutVarSyst_nDecLXY.eps";
    plotnames[1]="KF_CutVarSyst_sigmaonly_nDecLXY.eps";
    plotnames[2]="KF_CutVarSyst_ratioonly_nDecLXY.eps";
    plotnames[3]="KF_CutVarSyst_effprompt_nDecLXY.eps";
    plotnames[4]="KF_CutVarSyst_effFD_nDecLXY.eps";
    plotnames[5]="KF_CutVarSyst_rawyield_nDecLXY.eps";
    plotnames[6]="KF_CutVarSyst_sigma_nDecLXY.eps";
  }
  else if(var==4) { ///sigvtx variations
    setsno[1]=33;
    setsno[2]=23;
    setsno[3]=22;
    setsno[4]=21;
    setsno[5]=20;
    setsno[6]=32;
    legendnames[1]="sigvtx < central val + 0.006";
    legendnames[2]="sigvtx < central val + 0.004";
    legendnames[3]="sigvtx < central val + 0.002";
    legendnames[4]="sigvtx < central val - 0.002";
    legendnames[5]="sigvtx < central val - 0.004";
    legendnames[6]="sigvtx < central val - 0.006";
    plotnames[0]="KF_CutVarSyst_sigvtx.eps";
    plotnames[1]="KF_CutVarSyst_sigmaonly_sigvtx.eps";
    plotnames[2]="KF_CutVarSyst_ratioonly_sigvtx.eps";
    plotnames[3]="KF_CutVarSyst_effprompt_sigvtx.eps";
    plotnames[4]="KF_CutVarSyst_effFD_sigvtx.eps";
    plotnames[5]="KF_CutVarSyst_rawyield_sigvtx.eps";
    plotnames[6]="KF_CutVarSyst_sigma_sigvtx.eps";
  }
  else if(var==5) { ///cosp variations
    setsno[1]=34;
    setsno[2]=25;
    setsno[3]=24;
    setsno[4]=26;
    setsno[5]=27;
    setsno[6]=35;
    legendnames[1]="cosp < central val - 0.015";
    legendnames[2]="cosp < central val - 0.010";
    legendnames[3]="cosp < central val - 0.005";
    legendnames[4]="cosp < central val + 0.005";
    legendnames[5]="cosp < central val + 0.010";
    legendnames[6]="cosp < central val + 0.013";
    plotnames[0]="KF_CutVarSyst_cosp.eps";
    plotnames[1]="KF_CutVarSyst_sigmaonly_cosp.eps";
    plotnames[2]="KF_CutVarSyst_ratioonly_cosp.eps";
    plotnames[3]="KF_CutVarSyst_effprompt_cosp.eps";
    plotnames[4]="KF_CutVarSyst_effFD_cosp.eps";
    plotnames[5]="KF_CutVarSyst_rawyield_cosp.eps";
    plotnames[6]="KF_CutVarSyst_sigma_cosp.eps";
  }
  else {
    cerr<< "no more variations"<< endl;
    return;
  }
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLineWidth(1);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLabelFont(42,"xy");
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleSize(0.05,"xzt");
  gStyle->SetTitleSize(0.045,"y");
  gStyle->SetLegendBorderSize(0);
  TGaxis::SetMaxDigits(2);
  
  TLegend* l = new TLegend(0.35,0.5,0.85,0.85);
  TLegend* lRatio = new TLegend(0.2,0.5,0.85,0.85);
  l->SetTextSize(0.045);
  lRatio->SetTextSize(0.045);
  TLine* line = new TLine(1,1,24,1);
  line->SetLineColor(colors[0]);
  line->SetLineWidth(2);
  lRatio->AddEntry(line,"central set","l");
  
  TH1F** hCross = new TH1F*[nSets];
  TH1F** hRatio = new TH1F*[nSets-1];
  TGraphAsymmErrors** gFrac = new TGraphAsymmErrors*[nSets];

  for(Int_t iSet=0; iSet<nSets; iSet++) {
    TFile infile(Form("HFPtSpectrum_combinedFD_cutset%d.root",setsno[iSet]),"UPDATE");
    hCross[iSet] = (TH1F*)infile.Get("histoSigmaCorr");
    hCross[iSet]->SetDirectory(0);
    hCross[iSet]->SetMarkerStyle(20+iSet);
    hCross[iSet]->SetMarkerSize(1.5);
    hCross[iSet]->SetMarkerColor(colors[iSet]);
    hCross[iSet]->SetLineColor(colors[iSet]);
    hCross[iSet]->SetLineWidth(2);
    hCross[iSet]->SetTitle("");
    hCross[iSet]->GetXaxis()->SetLabelFont(42);
    hCross[iSet]->GetYaxis()->SetLabelFont(42);
    hCross[iSet]->GetXaxis()->SetTitleFont(42);
    hCross[iSet]->GetYaxis()->SetTitleFont(42);
    hCross[iSet]->GetXaxis()->SetTitleOffset(1.2);
    hCross[iSet]->GetYaxis()->SetTitleOffset(1.7);
    hCross[iSet]->GetXaxis()->SetTitleSize(0.05);
    hCross[iSet]->GetYaxis()->SetTitleSize(0.045);
    hCross[iSet]->GetXaxis()->SetLabelSize(0.05);
    hCross[iSet]->GetYaxis()->SetLabelSize(0.05);
    hCross[iSet]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hCross[iSet]->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}} (#mub c/GeV)");
    hCross[iSet]->Scale(1./(1000000*BR));
    l->AddEntry(hCross[iSet],Form("%s",legendnames[iSet].Data()),"lpe");
 
    if(iSet>0) {
      hRatio[iSet-1] = (TH1F*)hCross[iSet]->Clone();
      hRatio[iSet-1]->Divide(hCross[iSet],hCross[0],1.,1.,"");
      hRatio[iSet-1]->SetDirectory(0);
      hRatio[iSet-1]->GetYaxis()->SetTitle("ratio of #frac{d#sigma}{d#it{p}_{T}} w.r.t. the central value");
      for(Int_t iPt=0; iPt<hRatio[iSet-1]->GetNbinsX(); iPt++) {
        hRatio[iSet-1]->SetBinError(iPt+1,1.e-10);
      }
      lRatio->AddEntry(hRatio[iSet-1],Form("%s",legendnames[iSet].Data()),"lp");
    }

    if(var==1 || var==2) {
      gFrac[iSet] = (TGraphAsymmErrors*)infile.Get("gFcCorrConservative");
      gFrac[iSet]->SetMarkerStyle(20+iSet);
      gFrac[iSet]->SetMarkerSize(1.5);
      gFrac[iSet]->SetMarkerColor(colors[iSet]);
      gFrac[iSet]->SetLineColor(colors[iSet]);
      gFrac[iSet]->SetLineWidth(2);
      gFrac[iSet]->SetTitle("");
      gFrac[iSet]->GetXaxis()->SetLabelFont(42);
      gFrac[iSet]->GetYaxis()->SetLabelFont(42);
      gFrac[iSet]->GetXaxis()->SetTitleFont(42);
      gFrac[iSet]->GetYaxis()->SetTitleFont(42);
      gFrac[iSet]->GetXaxis()->SetTitleOffset(1.2);
      gFrac[iSet]->GetYaxis()->SetTitleOffset(1.5);
      gFrac[iSet]->GetXaxis()->SetTitleSize(0.05);
      gFrac[iSet]->GetYaxis()->SetTitleSize(0.05);
      gFrac[iSet]->GetXaxis()->SetLabelSize(0.05);
      gFrac[iSet]->GetYaxis()->SetLabelSize(0.05);
      gFrac[iSet]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
      gFrac[iSet]->GetYaxis()->SetTitle("#it{f}_{prompt}");   
    }
    infile.Close();
  }
  
  TH1F** hEffPrompt = new TH1F*[nSets];
  TH1F** hEffFD = new TH1F*[nSets];
  TH1F** hRatioEffPrompt = new TH1F*[nSets-1];
  TH1F** hRatioEffFD = new TH1F*[nSets-1];
  
  TLegend *lEff=(TLegend*)l->Clone();
  lEff->SetY1(0.2);
  lEff->SetY2(0.5);
  lEff->SetX1(0.4);
  
  for(Int_t iSet=0; iSet<nSets; iSet++) {
    TFile efffile(Form("efficiency_Dplus_cutset%d.root",setsno[iSet]),"UPDATE");
    hEffPrompt[iSet] = (TH1F*)efffile.Get("hEffD");
    hEffFD[iSet] = (TH1F*)efffile.Get("hEffB");
    hEffPrompt[iSet]->SetDirectory(0);
    hEffPrompt[iSet]->SetMarkerStyle(20+iSet);
    hEffPrompt[iSet]->SetMarkerSize(1.5);
    hEffPrompt[iSet]->SetMarkerColor(colors[iSet]);
    hEffPrompt[iSet]->SetLineColor(colors[iSet]);
    hEffPrompt[iSet]->SetLineWidth(2);
    hEffPrompt[iSet]->SetTitle("");
    hEffPrompt[iSet]->GetXaxis()->SetLabelFont(42);
    hEffPrompt[iSet]->GetYaxis()->SetLabelFont(42);
    hEffPrompt[iSet]->GetXaxis()->SetTitleFont(42);
    hEffPrompt[iSet]->GetYaxis()->SetTitleFont(42);
    hEffPrompt[iSet]->GetXaxis()->SetTitleOffset(1.2);
    hEffPrompt[iSet]->GetYaxis()->SetTitleOffset(1.4);
    hEffPrompt[iSet]->GetXaxis()->SetLabelSize(0.05);
    hEffPrompt[iSet]->GetYaxis()->SetLabelSize(0.05);
    hEffPrompt[iSet]->GetXaxis()->SetTitleSize(0.05);
    hEffPrompt[iSet]->GetYaxis()->SetTitleSize(0.045);
    hEffPrompt[iSet]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hEffPrompt[iSet]->GetYaxis()->SetTitle("#epsilon_{prompt}");
    hEffFD[iSet]->SetDirectory(0);
    hEffFD[iSet]->SetMarkerStyle(20+iSet);
    hEffFD[iSet]->SetMarkerSize(1.5);
    hEffFD[iSet]->SetMarkerColor(colors[iSet]);
    hEffFD[iSet]->SetLineColor(colors[iSet]);
    hEffFD[iSet]->SetLineWidth(2);
    hEffFD[iSet]->SetTitle("");
    hEffFD[iSet]->GetXaxis()->SetLabelFont(42);
    hEffFD[iSet]->GetYaxis()->SetLabelFont(42);
    hEffFD[iSet]->GetXaxis()->SetTitleFont(42);
    hEffFD[iSet]->GetYaxis()->SetTitleFont(42);
    hEffFD[iSet]->GetXaxis()->SetTitleOffset(1.2);
    hEffFD[iSet]->GetYaxis()->SetTitleOffset(1.3);
    hEffFD[iSet]->GetXaxis()->SetLabelSize(0.05);
    hEffFD[iSet]->GetYaxis()->SetLabelSize(0.05);
    hEffFD[iSet]->GetXaxis()->SetTitleSize(0.05);
    hEffFD[iSet]->GetYaxis()->SetTitleSize(0.045);
    hEffFD[iSet]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hEffFD[iSet]->GetYaxis()->SetTitle("#epsilon_{FD}");

    efffile.Close();
  }

  TLegend *lRatioPromptEff = new TLegend(0.35,0.2,0.85,0.55); 
  TLegend *lRatioFDEff = new TLegend(0.35,0.5,0.85,0.85);
  lRatioPromptEff->SetTextSize(0.045);
  lRatioFDEff->SetTextSize(0.045);
 
  if(var==1 || var==2) {
    for(Int_t iSet=0; iSet<nSets-1; iSet++) {
      hRatioEffPrompt[iSet] = (TH1F*)hEffPrompt[iSet]->Clone();
      hRatioEffPrompt[iSet]->Divide(hEffPrompt[iSet],hEffPrompt[nSets-1],1.,1.,"");
      hRatioEffPrompt[iSet]->SetDirectory(0);
      hRatioEffPrompt[iSet]->GetYaxis()->SetTitle("ratio of #epsilon_{prompt}  w.r.t. the central set w/o #chi^{2} cut");
      hRatioEffFD[iSet] = (TH1F*)hEffFD[iSet]->Clone();
      hRatioEffFD[iSet]->Divide(hEffFD[iSet],hEffFD[nSets-1],1.,1.,"");
      hRatioEffFD[iSet]->SetDirectory(0);
      hRatioEffFD[iSet]->GetYaxis()->SetTitle("ratio of #epsilon_{FD}  w.r.t. the central set w/o #chi^{2} cut");
      for(Int_t iPt=0; iPt<hRatioEffPrompt[iSet]->GetNbinsX(); iPt++) {
        hRatioEffPrompt[iSet]->SetBinError(iPt+1,1.e-10);
        hRatioEffFD[iSet]->SetBinError(iPt+1,1.e-10);
      }
      lRatioPromptEff->AddEntry(hRatioEffPrompt[iSet],Form("%s",legendnames[iSet].Data()),"lp");
      lRatioFDEff->AddEntry(hRatioEffPrompt[iSet],Form("%s",legendnames[iSet].Data()),"lp");
    }
  }
  
  TH1F** hRaw = new TH1F*[nSets];
  TH1F** hSigma = new TH1F*[nSets];
  TH1F** hRatioRaw = new TH1F*[nSets-1];
 
  for(Int_t iSet=0; iSet<nSets; iSet++) {
    TFile rawfile(Form("rawyields_Dplus_cutset%d.root",setsno[iSet]),"UPDATE");
    hRaw[iSet] = (TH1F*)rawfile.Get("hSignal");
    hSigma[iSet] = (TH1F*)rawfile.Get("hSigma");
    hRaw[iSet]->SetDirectory(0);
    hRaw[iSet]->SetMarkerStyle(20+iSet);
    hRaw[iSet]->SetMarkerSize(1.5);
    hRaw[iSet]->SetMarkerColor(colors[iSet]);
    hRaw[iSet]->SetLineColor(colors[iSet]);
    hRaw[iSet]->SetLineWidth(2);
    hRaw[iSet]->SetTitle("");
    hRaw[iSet]->GetXaxis()->SetLabelFont(42);
    hRaw[iSet]->GetYaxis()->SetLabelFont(42);
    hRaw[iSet]->GetXaxis()->SetTitleFont(42);
    hRaw[iSet]->GetYaxis()->SetTitleFont(42);
    hRaw[iSet]->GetXaxis()->SetTitleOffset(1.2);
    hRaw[iSet]->GetYaxis()->SetTitleOffset(1.3);
    hRaw[iSet]->GetXaxis()->SetTitleSize(0.05);
    hRaw[iSet]->GetYaxis()->SetTitleSize(0.05);
    hRaw[iSet]->GetXaxis()->SetLabelSize(0.05);
    hRaw[iSet]->GetYaxis()->SetLabelSize(0.05);
    hRaw[iSet]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hRaw[iSet]->GetYaxis()->SetTitle("Raw Yield");
    hSigma[iSet]->SetDirectory(0);
    hSigma[iSet]->SetMarkerStyle(20+iSet);
    hSigma[iSet]->SetMarkerSize(1.5);
    hSigma[iSet]->SetMarkerColor(colors[iSet]);
    hSigma[iSet]->SetLineColor(colors[iSet]);
    hSigma[iSet]->SetLineWidth(2);
    hSigma[iSet]->SetTitle("");
    hSigma[iSet]->GetXaxis()->SetLabelFont(42);
    hSigma[iSet]->GetYaxis()->SetLabelFont(42);
    hSigma[iSet]->GetXaxis()->SetTitleFont(42);
    hSigma[iSet]->GetYaxis()->SetTitleFont(42);
    hSigma[iSet]->GetXaxis()->SetTitleOffset(1.2);
    hSigma[iSet]->GetYaxis()->SetTitleOffset(1.3);
    hSigma[iSet]->GetXaxis()->SetLabelSize(0.05);
    hSigma[iSet]->GetYaxis()->SetLabelSize(0.05);
    hSigma[iSet]->GetXaxis()->SetTitleSize(0.05);
    hSigma[iSet]->GetYaxis()->SetTitleSize(0.045);
    hSigma[iSet]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hSigma[iSet]->GetYaxis()->SetTitle("Sigma");
  }

  if(var==1 || var==2) {
    for(Int_t iSet=0; iSet<nSets-1; iSet++) {
      hRatioRaw[iSet]=(TH1F*)hRaw[iSet]->Clone();
      hRatioRaw[iSet]->Divide(hRaw[iSet],hRaw[nSets-1],1.,1.,"");
      hRatioRaw[iSet]->SetDirectory(0);
      hRatioRaw[iSet]->GetYaxis()->SetTitle("ratio of Y w.r.t. the central set w/o #chi^{2} cut");
      for(Int_t iPt=0; iPt<hRatioEffPrompt[iSet]->GetNbinsX(); iPt++) {
        hRatioRaw[iSet]->SetBinError(iPt+1,1.e-10);
      }
    }
  }
  
  TCanvas *cCross = new TCanvas("cCross","",800,800);
  cCross->SetLogy();
  hCross[0]->Draw();
  l->Draw("same");
  
  TCanvas *cRatio = new TCanvas("cRatio","",800,800);
  hRatio[0]->GetYaxis()->SetRangeUser(0.8,1.7);
  hRatio[0]->Draw();
  line->Draw("same");
  lRatio->Draw("same");
  
  TCanvas *c = new TCanvas("c","",1200,600);
  c->Divide(2,1);
  c->cd(1);
  hCross[0]->Draw();
  l->Draw("same");
  c->cd(2);
  hRatio[0]->Draw();
  line->Draw("same");
  lRatio->Draw("same");
  
  for(Int_t iSet=1; iSet<nSets; iSet++) {
    cCross->cd();
    hCross[iSet]->Draw("same");
    c->cd(1)->SetLogy();
    hCross[iSet]->Draw("same");

    if(iSet<nSets-1) {
      cRatio->cd();
      hRatio[iSet]->Draw("same");
      c->cd(2);
      hRatio[iSet]->Draw("same");
    }
  }

  c->SaveAs(plotnames[0].Data());
  cCross->SaveAs(plotnames[1].Data());
  cRatio->SaveAs(plotnames[2].Data());

  TCanvas* cEffPrompt = new TCanvas("cEffPrompt","",800,800);
  cEffPrompt->SetLogy();
  hEffPrompt[0]->GetYaxis()->SetRangeUser(0.001,0.5);
  hEffPrompt[0]->Draw();
  for(Int_t iSet=1; iSet<nSets; iSet++) {
    hEffPrompt[iSet]->Draw("same");
  }
  lEff->Draw("same");
  
  TCanvas* cEffFD = new TCanvas("cEffFD","",800,800);
  cEffFD->SetLogy();
  hEffFD[0]->GetYaxis()->SetRangeUser(0.001,0.5);
  hEffFD[0]->Draw();
  for(Int_t iSet=1; iSet<nSets; iSet++) {
    hEffFD[iSet]->Draw("same");
  }
  lEff->Draw("same");

  TCanvas* cRaw = new TCanvas("cRaw","",800,800);
  hRaw[0]->Draw();
  hRaw[0]->GetYaxis()->SetRangeUser(0.,2500.);
  for(Int_t iSet=1; iSet<nSets; iSet++) {
    hRaw[iSet]->Draw("same");
  }
  l->Draw("same");
  
  TCanvas* cSigma = new TCanvas("cSigma","",800,800);
  hSigma[0]->Draw();
  for(Int_t iSet=1; iSet<nSets; iSet++) {
    hSigma[iSet]->Draw("same");
  }
  lEff->Draw("same");

  cEffPrompt->SaveAs(plotnames[3].Data());
  cEffFD->SaveAs(plotnames[4].Data());
  cRaw->SaveAs(plotnames[5].Data());
  cSigma->SaveAs(plotnames[6].Data());

  if(var==1 || var==2) {
    TCanvas* cRatioEffPrompt = new TCanvas("cRatioEffPrompt","",800,800);
    hRatioEffPrompt[0]->GetYaxis()->SetRangeUser(0.5,1.02);
    hRatioEffPrompt[0]->GetYaxis()->SetTitleOffset(1.5);
    hRatioEffPrompt[0]->Draw();
    for(Int_t iSet=1; iSet<nSets-1; iSet++) {
      hRatioEffPrompt[iSet]->Draw("same");
    }
    lRatioPromptEff->Draw("same");
    
    TCanvas* cRatioEffFD = new TCanvas("cRatioEffFD","",800,800);
    if(var==1) {
      hRatioEffFD[0]->GetYaxis()->SetRangeUser(0.3,1.2);
      hRatioEffFD[0]->GetYaxis()->SetTitleOffset(1.5);
    }
    else
      hRatioEffFD[0]->GetYaxis()->SetRangeUser(0.7,1.25);
    hRatioEffFD[0]->Draw();
    for(Int_t iSet=1; iSet<nSets-1; iSet++) {
      hRatioEffFD[iSet]->Draw("same");
    }
    lRatioFDEff->Draw("same");
    
    cRatioEffPrompt->SaveAs(plotnames[7].Data());
    cRatioEffFD->SaveAs(plotnames[8].Data());

    TCanvas* cRatioRaw = new TCanvas("cRatioRaw","",800,800);
    hRatioRaw[0]->GetYaxis()->SetRangeUser(0.2,1.05);
    hRatioRaw[0]->Draw();
    for(Int_t iSet=1; iSet<nSets-1; iSet++) {
      hRatioRaw[iSet]->Draw("same");
    }
    lRatioPromptEff->Draw("same");
    
    cRatioRaw->SaveAs(plotnames[9].Data());

    TLegend* lfrac = (TLegend*)lRatioPromptEff->Clone();
    lfrac->AddEntry(gFrac[nSets-1],"central w/o #chi^{2} cut","lpe");
    lfrac->SetY1(0.25);
    lfrac->SetFillStyle(0);
    TCanvas* cFrac = new TCanvas("cFrac","",800,800);
    gFrac[nSets-1]->GetYaxis()->SetRangeUser(0.5,1.);
    gFrac[nSets-1]->Draw("AP");
    for(Int_t iSet=0; iSet<nSets-1; iSet++) {
      gFrac[iSet]->Draw("P");
    }
    lfrac->Draw("same");
    cFrac->SaveAs(plotnames[10].Data());    
  }
}

void ExtimateUncertainty() {

  gStyle->SetOptStat(1111);
  gStyle->SetStatFontSize(0.06);
  gStyle->SetOptFit(0);
  gStyle->SetLineWidth(1);
  gStyle->SetLabelFont(42,"xy");
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleSize(0.05,"xzt");
  gStyle->SetTitleSize(0.045,"y");
  gStyle->SetLegendBorderSize(0);
  TGaxis::SetMaxDigits(2);
  Int_t nSets=29;
  Int_t varsets[nSets] = {2,3,4,5,6,7,8,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35};
  TH1F** hRatio = new TH1F*[nSets];
  TH1F** hCross = new TH1F*[nSets];
  
  TFile reffile("HFPtSpectrum_combinedFD_cutset1.root","UPDATE");
  TH1F* hCrossRef = (TH1F*)reffile.Get("histoSigmaCorr");
  hCrossRef->SetDirectory(0);
  reffile.Close();

  Int_t nPtBins = hCrossRef->GetNbinsX();
  TH1F** hDisp = new TH1F*[nPtBins];

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    hDisp[iPt] = new TH1F(Form("hDisp_Pt%d",iPt),"",50,0.75,1.25);
    hDisp[iPt]->SetDirectory(0);
  }

  for(Int_t iSet=0; iSet<nSets; iSet++) {
    TFile infile(Form("HFPtSpectrum_combinedFD_cutset%d.root",varsets[iSet]),"UPDATE"); 
    hCross[iSet] = (TH1F*)infile.Get("histoSigmaCorr");
    hCross[iSet]->SetDirectory(0);
    hRatio[iSet]=(TH1F*)hCross[iSet]->Clone();
    hRatio[iSet]->SetDirectory(0);
    hRatio[iSet]->Divide(hCross[iSet],hCrossRef,1.,1.,"");
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      hDisp[iPt]->Fill(hRatio[iSet]->GetBinContent(iPt+1));
    }
    infile.Close();
  }
  TCanvas* c = new TCanvas("c","",1200,900);
  c->Divide(4,3);
  TCanvas** c2 = new TCanvas*[nPtBins];
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    c->cd(iPt+1);
    hDisp[iPt]->SetName("Cut variation");
    hDisp[iPt]->SetStats(1111);
    hDisp[iPt]->SetFillColor(kBlue);
    hDisp[iPt]->SetFillStyle(3004);
    hDisp[iPt]->GetYaxis()->SetTitle("Entries");
    hDisp[iPt]->GetXaxis()->SetTitle("ratio of #frac{d#sigma}{d#it{p}_{T}} w.r.t. the central value");
    hDisp[iPt]->GetYaxis()->SetTitleSize(0.05);
    hDisp[iPt]->GetXaxis()->SetTitleSize(0.05);
    hDisp[iPt]->GetYaxis()->SetLabelSize(0.05);
    hDisp[iPt]->GetXaxis()->SetLabelSize(0.05);
    hDisp[iPt]->GetYaxis()->SetTitleOffset(1.2);
    hDisp[iPt]->GetXaxis()->SetTitleOffset(1.5);
    Double_t centre = hRatio[0]->GetBinCenter(iPt+1);
    Double_t step = hRatio[0]->GetBinWidth(iPt+1)/2;
    hDisp[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",centre-step,centre+step));
    hDisp[iPt]->GetXaxis()->SetNdivisions(508);
    hDisp[iPt]->Draw();
    c2[iPt] = new TCanvas(Form("c2_%d",iPt),"",800,800);
    c2[iPt]->SetBottomMargin(0.18);
    hDisp[iPt]->Draw();
    c2[iPt]->SaveAs(Form("KF_CutVarSyst_Disp_Pt%d.eps",iPt));
  }

  for(Int_t iPt=0; iPt<nPtBins; iPt++)
    delete c2[iPt];
  delete[] c2;
    
  TFile outfile("KF_CutVarSyst_RMS.root","RECREATE");
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    hDisp[iPt]->Write();
  }
  c->Write();
  outfile.Close();
  
}
