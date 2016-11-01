#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <string>
#include <sstream>
#include <fstream>
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
#include <TDatabasePDG.h>
#include <TString.h>
#include <TLine.h>
#include <TList.h>
#include <TDirectoryFile.h>
#include <TLatex.h>

#include "AliHFMassFitter.h"
#include "AliHFCutVarFDsubAxis.h"
#include "AliHFCutVarFDsubCut.h"
#include "AliHFCutVarFDsubCutSet.h"
#include "AliHFCutVarFDsubEfficiency.h"
#include "AliHFCutVarFDsubMassFitter.h"
#include "AliHFCutVarFDsubMinimiser.h"
#include "AliHFCutVarFDsubAnalysisManager.h"
#include "AliHFCutVarFDsubAnalysisManagerDplus.h"

#endif

enum method {kMinimise,kIncentre};
Int_t promptcolors[] = {kBlue,kGreen+3};
Int_t FDcolors[] = {kRed,kMagenta};

void EfficiencySyst();
void ReadAxes(TString FileName, vector<string> &axesanmes, vector<int> &axesno);
void ReadSet(TString FileName, vector<string> &varnames, vector<double> &cutset);

void EfficiencySyst() {

  //________________________________________________________________________________________________________________
  //INPUT FILES
  //Input files and output directories
  TString datafile = "/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/Data/AnalysisResultspPbData.root";
  TString MCfile = "/home/fabrizio/tesi/Trains/ImpPar_and_CutVar/pPb/MC/AnalysisResultspPbMC.root";
  TString datalist = "coutputDplus_CutVarpPbData0100";
  TString MClist = "coutputDplus_CutVarpPbMC0100";
  
  TString datadir = "PWG3_D2H_InvMassDplus";
  TString MCdir = "PWG3_D2H_InvMassDplus";

  TString axesfile = "../../Axes/Axes";
  TString cutfile1 = "../../cutSets/set1_a";
  TString cutfile2 = "../../cutSets/set2_a";
  TString cutfile3 = "../../cutSets/set3_a";

  TString NotRewDir = "NotRewEff";
  TString RewDir = "RewEff";
  
  //Pt reweigth files
  TFile Dfile("PtD.root","UPDATE");
  TF1* funcWeightsD=(TF1*)Dfile.Get("fFuncWeight"); 
  Dfile.Close();
  TFile Bfile("PtB.root","UPDATE");
  TF1* funcWeightsB=(TF1*)Bfile.Get("fFuncWeight"); 
  Bfile.Close(); 

  //_________________________________________________________________________________________________________________
  //Import cuts from files
  const Int_t nSets=3;
  
  vector<string> axesnames;
  vector<int> axesno;
  ReadAxes(axesfile,axesnames,axesno);
  const Int_t nCutVars=axesnames.size(); //pt included

  UInt_t dataAxesNo[nCutVars];
  UInt_t MCgenAxesNo[nCutVars];
  UInt_t MCcutAxesNo[nCutVars];
  TString AxesNames[nCutVars];

  cout << "\n________________________________ AXES ___________________________________ \n\n Axis Name" <<setw(20) <<"DataAxis" << setw(20) << "MCGenAxis" << setw(20) << " MCRecoAxis\n"<< "_________________________________________________________________________\n" <<endl;
  
  for(Int_t iCutVar=0; iCutVar<nCutVars; iCutVar++) {
    AxesNames[iCutVar] = axesnames[iCutVar];
    dataAxesNo[iCutVar] = axesno[iCutVar];
    MCgenAxesNo[iCutVar] = axesno[iCutVar+nCutVars];
    MCcutAxesNo[iCutVar] = axesno[iCutVar+2*nCutVars];
    cout.width(5); cout <<AxesNames[iCutVar];
    cout.width(20); cout << dataAxesNo[iCutVar];
    cout.width(25); cout << MCgenAxesNo[iCutVar];
    cout.width(15); cout << MCcutAxesNo[iCutVar] << endl;
  }

  cout << "_____________________________________________________________________________\n" << endl;

  vector<string> varnames;
  vector<double> cutset1;
  vector<double> cutset2;
  vector<double> cutset3;
  ReadSet(cutfile1,varnames,cutset1);
  ReadSet(cutfile2,varnames,cutset2);
  ReadSet(cutfile3,varnames,cutset3);
  Int_t nPtBins= cutset1.size()/(2*nCutVars);
  const Int_t nptbins = nPtBins;
  
  Double_t*** cutlowset = new Double_t**[nSets]; //first: set, second: pt bin, third: cut variable
  Double_t*** cuthighset = new Double_t**[nSets];

  for(Int_t iSet=0; iSet<nSets; iSet++) {
    cutlowset[iSet] = new Double_t*[nPtBins];
    cuthighset[iSet] = new Double_t*[nPtBins];
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      cutlowset[iSet][iPt] = new Double_t[nCutVars];
      cuthighset[iSet][iPt] = new Double_t[nCutVars];
    }    
  }

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    for(Int_t iCutVar=0; iCutVar<nCutVars; iCutVar++) {
      cutlowset[0][iPt][iCutVar] = cutset1[(2*iCutVar)+(iPt*nCutVars*2)];
      cuthighset[0][iPt][iCutVar] = cutset1[(2*iCutVar+1)+(iPt*nCutVars*2)];      
      cutlowset[1][iPt][iCutVar] = cutset2[(2*iCutVar)+(iPt*nCutVars*2)];
      cuthighset[1][iPt][iCutVar] = cutset2[(2*iCutVar+1)+(iPt*nCutVars*2)];      
      cutlowset[2][iPt][iCutVar] = cutset3[(2*iCutVar)+(iPt*nCutVars*2)];
      cuthighset[2][iPt][iCutVar] = cutset3[(2*iCutVar+1)+(iPt*nCutVars*2)];      
    }
  }

  for(Int_t iSet=0; iSet<nSets; iSet++) {
    cout << "\n" <<Form("_______________________________________________ CUT SET %d ________________________________________________",iSet+1)<<"\n" << endl;
    for(Int_t iVarName=0; iVarName < nCutVars*2; iVarName++)
      cout << varnames[iVarName]<<"    ";//setw(11);
    cout << endl;
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      for(Int_t iCutVar=0; iCutVar<nCutVars; iCutVar++) {
        if(iCutVar>0) {
          cout.width(11); cout << cutlowset[iSet][iPt][iCutVar];
          cout.width(11); cout << cuthighset[iSet][iPt][iCutVar];
        }
        else {
          cout << cutlowset[iSet][iPt][iCutVar];
          cout.width(11); cout << cuthighset[iSet][iPt][iCutVar];
        }
      }
      cout << endl;
    }
    cout << "\n__________________________________________________________________________________________________________\n\n" << endl;
  }

 //_________________________________________________________________________________________________________________
  //means and sigmas for fit
  Double_t massD = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  
  Double_t **means = new Double_t*[nSets];
  Double_t **sigmasfromMC = new Double_t*[nSets];
  
  for(Int_t iSet=0; iSet<nSets; iSet++) {
    means[iSet] = new Double_t[nPtBins];
    sigmasfromMC[iSet] = new Double_t[nPtBins];    
  }

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    means[0][iPt] = -massD;
    means[1][iPt] = -massD;
    means[2][iPt] = -massD;
  }

  TFile *massfitsfile = TFile::Open("../../massfitsMCtruth/RawYields.root","READ");
  if(massfitsfile) {
    for(Int_t iSet = 0; iSet<nSets; iSet++) {
      TH1F* hSigma = (TH1F*)massfitsfile->Get(Form("hRawYieldSigmas_Set%d",iSet+1));
      for(Int_t iPt=0; iPt<nPtBins; iPt++) {
        sigmasfromMC[iSet][iPt] = hSigma->GetBinContent(iPt+1);
      }
    }
  }
  else {
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      sigmasfromMC[0][iPt] = -0.008;
      sigmasfromMC[1][iPt] = -0.008;
      sigmasfromMC[2][iPt] = -0.008;
    }  
  }

  //________________________________________________________________________________________________________________
  //Analysis
  TH1F*** hCorrYieldsPrompt = new TH1F**[2];
  TH1F*** hCorrYieldsFD = new TH1F**[2];
  TH1F** hRatioPrompt = new TH1F*[2];
  TH1F** hRatioFD = new TH1F*[2];
  
    for(Int_t iMethod=kMinimise; iMethod<=kIncentre; iMethod++) {
    hCorrYieldsPrompt[iMethod] = new TH1F*[2];
    hCorrYieldsFD[iMethod] = new TH1F*[2];
    
    for(Int_t iRew=0; iRew<=1; iRew++) {     

      AliHFCutVarFDsubAnalysisManagerDplus *AnalysisManagerDplus = new AliHFCutVarFDsubAnalysisManagerDplus();  
      AnalysisManagerDplus->SetPID(kTRUE,3);
      Int_t loadTH = AnalysisManagerDplus->GetTHnSparses(MCfile,datafile,MCdir,datadir,MClist,datalist,kFALSE);
      
      if(loadTH>0)
        return;
      
      AnalysisManagerDplus->GetAxes(dataAxesNo,MCgenAxesNo,MCcutAxesNo,AxesNames,nCutVars);
      AnalysisManagerDplus->GetCuts(cutlowset,cuthighset,means,sigmasfromMC,4,0,0,1.72,2.05,nSets,nPtBins,nCutVars);
      AnalysisManagerDplus->GetXaxisInformation();
      if(iRew==0) {
        AnalysisManagerDplus->GetEfficiencies(NotRewDir);
        AnalysisManagerDplus->DrawEfficiencies(NotRewDir);
      }
      else {
        AnalysisManagerDplus->GetEfficiencies(RewDir,kTRUE,funcWeightsD,funcWeightsB);
        AnalysisManagerDplus->DrawEfficiencies(RewDir);
      }
      AnalysisManagerDplus->GetRawYields(kFALSE);
      Bool_t min = AnalysisManagerDplus->Minimise(iMethod,10,kTRUE,0.,10000);

      hCorrYieldsPrompt[iMethod][iRew] = (TH1F*)AnalysisManagerDplus->GetYieldsPrompt();
      hCorrYieldsFD[iMethod][iRew] = (TH1F*)AnalysisManagerDplus->GetYieldsFD();
      
      hCorrYieldsPrompt[iMethod][iRew]->SetDirectory(0);
      hCorrYieldsFD[iMethod][iRew]->SetDirectory(0);
      
      hCorrYieldsPrompt[iMethod][iRew]->GetYaxis()->SetTitle("#frac{dN_{Prompt}}{d#it{p}_{T}} (c/GeV)");
      hCorrYieldsFD[iMethod][iRew]->GetYaxis()->SetTitle("#frac{dN_{FD}}{d#it{p}_{T}} (c/GeV)");
      hCorrYieldsPrompt[iMethod][iRew]->GetYaxis()->SetTitleSize(0.05);
      hCorrYieldsFD[iMethod][iRew]->GetYaxis()->SetTitleSize(0.05);
      hCorrYieldsPrompt[iMethod][iRew]->GetXaxis()->SetTitleSize(0.05);
      hCorrYieldsFD[iMethod][iRew]->GetXaxis()->SetTitleSize(0.05);
      hCorrYieldsPrompt[iMethod][iRew]->GetYaxis()->SetLabelSize(0.05);
      hCorrYieldsFD[iMethod][iRew]->GetYaxis()->SetLabelSize(0.05);
      hCorrYieldsPrompt[iMethod][iRew]->GetXaxis()->SetLabelSize(0.05);
      hCorrYieldsFD[iMethod][iRew]->GetXaxis()->SetLabelSize(0.05);
      hCorrYieldsPrompt[iMethod][iRew]->GetYaxis()->SetTitleOffset(1.1);
      hCorrYieldsFD[iMethod][iRew]->GetYaxis()->SetTitleOffset(1.1);
      hCorrYieldsPrompt[iMethod][iRew]->GetXaxis()->SetTitleOffset(1.1);
      hCorrYieldsFD[iMethod][iRew]->GetXaxis()->SetTitleOffset(1.1);
      
      hCorrYieldsPrompt[iMethod][iRew]->SetLineColor(promptcolors[iRew]);
      hCorrYieldsFD[iMethod][iRew]->SetLineColor(FDcolors[iRew]);
      hCorrYieldsPrompt[iMethod][iRew]->SetMarkerColor(promptcolors[iRew]);
      hCorrYieldsFD[iMethod][iRew]->SetMarkerColor(FDcolors[iRew]);
      hCorrYieldsPrompt[iMethod][iRew]->SetLineWidth(2);
      hCorrYieldsFD[iMethod][iRew]->SetLineWidth(2);
      hCorrYieldsPrompt[iMethod][iRew]->SetMarkerSize(1.5);
      hCorrYieldsFD[iMethod][iRew]->SetMarkerSize(1.5);
      hCorrYieldsPrompt[iMethod][iRew]->SetMarkerStyle(20);
      hCorrYieldsFD[iMethod][iRew]->SetMarkerStyle(21);
    }
    hRatioPrompt[iMethod]=(TH1F*)hCorrYieldsPrompt[iMethod][0]->Clone();
    hRatioPrompt[iMethod]->Divide(hCorrYieldsPrompt[iMethod][1],hCorrYieldsPrompt[iMethod][0],1.,1.);
    hRatioPrompt[iMethod]->SetDirectory(0);
    hRatioPrompt[iMethod]->GetYaxis()->SetTitle("Prompt corrected yields ratio");
    hRatioPrompt[iMethod]->GetXaxis()->SetTitle(hCorrYieldsPrompt[iMethod][0]->GetXaxis()->GetTitle());
    hRatioPrompt[iMethod]->GetYaxis()->SetTitleSize(0.05);
    hRatioPrompt[iMethod]->GetXaxis()->SetTitleSize(0.05);
    hRatioPrompt[iMethod]->GetYaxis()->SetLabelSize(0.05);
    hRatioPrompt[iMethod]->GetXaxis()->SetLabelSize(0.05);
    hRatioPrompt[iMethod]->GetYaxis()->SetTitleOffset(1.5);
    hRatioPrompt[iMethod]->GetXaxis()->SetTitleOffset(1.);
    hRatioPrompt[iMethod]->SetDirectory(0);
    hRatioFD[iMethod]=(TH1F*)hCorrYieldsFD[iMethod][0]->Clone();
    hRatioFD[iMethod]->Divide(hCorrYieldsFD[iMethod][1],hCorrYieldsFD[iMethod][0],1.,1.);
    hRatioFD[iMethod]->SetDirectory(0);
    hRatioFD[iMethod]->GetYaxis()->SetTitle("FD corrected yields ratio");
    hRatioFD[iMethod]->GetXaxis()->SetTitle(hCorrYieldsFD[iMethod][0]->GetXaxis()->GetTitle());
    hRatioFD[iMethod]->GetYaxis()->SetTitleSize(0.05);
    hRatioFD[iMethod]->GetXaxis()->SetTitleSize(0.05);
    hRatioFD[iMethod]->GetYaxis()->SetLabelSize(0.05);
    hRatioFD[iMethod]->GetXaxis()->SetLabelSize(0.05);
    hRatioFD[iMethod]->GetYaxis()->SetTitleOffset(1.5);
    hRatioFD[iMethod]->GetXaxis()->SetTitleOffset(1.);
  }

  //________________________________________________________________________________________________________
  //Plots
  
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetOptStat(0);

  TLine* line = new TLine(2,1,16,1);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->SetLineColor(kBlack);
  
  TLegend *lPrompt = new TLegend(0.45,0.65,0.8,0.89);
  lPrompt->SetTextSize(0.05);
  lPrompt->SetBorderSize(0);
  lPrompt->AddEntry(hCorrYieldsPrompt[0][0], "w/o #it{p}_{T} reweight","lpe");
  lPrompt->AddEntry(hCorrYieldsPrompt[0][1], "with #it{p}_{T} reweight","lpe");

  TLegend *lPromptRatio = new TLegend(0.4,0.7,0.85,0.89); 
  lPromptRatio->SetTextSize(0.045);
  lPromptRatio->SetBorderSize(0);
  lPromptRatio->AddEntry(hRatioPrompt[0], "#frac{with #it{p}_{T} reweight}{w/o #it{p}_{T} reweight}","l");

  TLegend *lFD = new TLegend(0.45,0.65,0.8,0.89); 
  lFD->SetTextSize(0.05);
  lFD->SetBorderSize(0);
  lFD->AddEntry(hCorrYieldsFD[0][0], "w/o #it{p}_{T} reweight","lpe");
  lFD->AddEntry(hCorrYieldsFD[0][1], "with #it{p}_{T} reweight","lpe");

  TLegend *lFDRatio = new TLegend(0.4,0.7,0.85,0.89); 
  lFDRatio->SetTextSize(0.05);
  lFDRatio->SetBorderSize(0);
  lFDRatio->AddEntry(hRatioFD[0], "#frac{with #it{p}_{T} reweight}{w/o #it{p}_{T} reweight}","l");

  for(Int_t iPt=0; iPt<hCorrYieldsPrompt[0][0]->GetNbinsX(); iPt++) {
    hRatioPrompt[0]->SetBinError(iPt+1,1.e-10);
    hRatioFD[0]->SetBinError(iPt+1,1.e-10);
    hRatioPrompt[1]->SetBinError(iPt+1,1.e-10);
    hRatioFD[1]->SetBinError(iPt+1,1.e-10);
  }
    
  TCanvas* cPromptMin = new TCanvas("cPromptMin","",1200,600);
  cPromptMin->Divide(2,1);
  cPromptMin->cd(1)->SetLogy();
  hCorrYieldsPrompt[0][0]->Draw();
  hCorrYieldsPrompt[0][1]->Draw("same");
  lPrompt->Draw("same");
  cPromptMin->cd(2);
  hRatioPrompt[0]->GetYaxis()->SetRangeUser(0.92,1.06);
  hRatioPrompt[0]->Draw("hist");
  lPromptRatio->Draw("same");
  cPromptMin->SaveAs("CorrYieldsPromptMin_syst_eff.eps");
  cPromptMin->SaveAs("CorrYieldsPromptMin_syst_eff.root");

  TCanvas* cFDMin = new TCanvas("cFDMin","",1200,600);
  cFDMin->Divide(2,1);
  cFDMin->cd(1)->SetLogy();
  hCorrYieldsFD[0][0]->Draw();
  hCorrYieldsFD[0][1]->Draw("same");
  lFD->Draw("same");
  cFDMin->cd(2);
  hRatioFD[0]->GetYaxis()->SetRangeUser(0.94,1.15);
  hRatioFD[0]->Draw("E");
  lFDRatio->Draw("same");
  cFDMin->SaveAs("CorrYieldsFDMin_syst_eff.eps");
  cFDMin->SaveAs("CorrYieldsFDMin_syst_eff.root");

  TLegend* lRatio = new TLegend(0.2,0.7,0.8,0.89);
  lRatio->SetTextSize(0.05);
  lRatio->SetBorderSize(0);
  lRatio->SetFillStyle(0);
  lRatio->AddEntry(hRatioPrompt[0],"Prompt","lpe");
  lRatio->AddEntry(hRatioFD[0],"Feed-down","lpe");
  line->Draw("same");
  
  TCanvas* cRatiosMin = new TCanvas("cRatiosMin","",800,800);
  hRatioPrompt[0]->SetTitle("Analytic minimisation");
  hRatioPrompt[0]->GetYaxis()->SetTitle("#frac{dN}{d#it{p}_{T}}(FONLL)/#frac{dN}{d#it{p}_{T}}(MC)");
  hRatioPrompt[0]->GetYaxis()->SetRangeUser(0.895,1.155);
  hRatioPrompt[0]->Draw("E");
  hRatioFD[0]->Draw("Esame");
  lRatio->Draw("same");
  cRatiosMin->SaveAs("CorrYieldsMin_syst_eff_onlyratio.eps");
  
  TCanvas* cPromptInc = new TCanvas("cPromptInc","",1200,600);
  cPromptInc->Divide(2,1);
  cPromptInc->cd(1)->SetLogy();
  hCorrYieldsPrompt[1][0]->Draw();
  hCorrYieldsPrompt[1][1]->Draw("same");
  lPrompt->Draw("same");
  cPromptInc->cd(2);
  hRatioPrompt[1]->GetYaxis()->SetRangeUser(0.92,1.06);
  hRatioPrompt[1]->Draw("hist");
  lPromptRatio->Draw("same");
  cPromptInc->SaveAs("CorrYieldsPromptInc_syst_eff.eps");
  cPromptInc->SaveAs("CorrYieldsPromptInc_syst_eff.root");

  TCanvas* cFDInc = new TCanvas("cFDInc","",1200,600);
  cFDInc->Divide(2,1);
  cFDInc->cd(1)->SetLogy();
  hCorrYieldsFD[1][0]->Draw();
  hCorrYieldsFD[1][1]->Draw("same");
  lFD->Draw("same");
  cFDInc->cd(2);
  hRatioFD[1]->GetYaxis()->SetRangeUser(0.94,1.15);
  hRatioFD[1]->Draw("hist");
  lFDRatio->Draw("same");
  cFDInc->SaveAs("CorrYieldsFDInc_syst_eff.eps");
  cFDInc->SaveAs("CorrYieldsFDInc_syst_eff.root");

  TCanvas* cRatiosInc = new TCanvas("cRatiosInc","",800,800);
  hRatioPrompt[1]->SetTitle("Incentre minimisation");
  hRatioPrompt[1]->GetYaxis()->SetTitle("#frac{dN}{d#it{p}_{T}}(FONLL)/#frac{dN}{d#it{p}_{T}}(MC)");
  hRatioPrompt[1]->GetYaxis()->SetRangeUser(0.895,1.205);
  hRatioPrompt[1]->Draw("E");
  hRatioFD[1]->Draw("Esame");
  lRatio->Draw("same");
  cRatiosMin->SaveAs("CorrYieldsMin_syst_eff_onlyratio.eps");
  
}

void ReadAxes(TString FileName, vector<string> &axesnames, vector<int> &axesno) {  
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
