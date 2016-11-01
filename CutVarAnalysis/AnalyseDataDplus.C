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
#include <TDirectory.h>

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

//macro for the Cut Variation Analyis for the D+
//author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino

//****************************************//
//                                        //
//    Main Function: AnalyseDataDplus     //
//                                        // 
//****************************************//

Int_t AnalyseDataDplus(Bool_t MCtruth=kFALSE, Bool_t isIncentre=kFALSE, Bool_t fixsigmafromMC=kTRUE, Bool_t fixmeanfromPDG=kFALSE, Int_t Rebin=4, Int_t funsig=0, Int_t funbkg=0, Bool_t PID=kTRUE);
void ReadAxes(TString FileName, vector<string> &axesanmes, vector<int> &axesno);
void ReadSet(TString FileName, vector<string> &varnames, vector<double> &cutset);

Int_t AnalyseDataDplus(Bool_t MCtruth,Bool_t isIncentre, Bool_t fixsigmafromMC, Bool_t fixmeanfromPDG, Int_t Rebin,Int_t funsig,Int_t funback,Bool_t PID) {

  TString AccFileName = "/home/fabrizio/ALICE/Files/Acceptance/Acceptance_Toy_DplusKpipi_yfidPtDep_etaDau09_ptDau100_FONLL5ptshape.root"; //acceptance file

  //________________________________________________________________________________________________________________
  //Input files and output directories
  
  TString datafile = "/home/fabrizio/ALICE/Files/Trains/LHC13/AnalysisResultspPbData.root";
  TString MCfile = "/home/fabrizio/ALICE/Files/Trains/LHC13/AnalysisResultspPbMC.root";
  TString datalist = "coutputDplus_CutVarpPbData0100";
  TString MClist = "coutputDplus_CutVarpPbMC0100";
  
  TString datadir = "PWG3_D2H_InvMassDplus";
  TString MCdir = "PWG3_D2H_InvMassDplus";

  TString axesfile = "AxesSymm";
  TString cutfile1 = "set1_a";
  TString cutfile2 = "set2_a";
  TString cutfile3 = "set3_a";
  
  TString distdir = ".";
  TString linesdir = ".";
  TString massfitsdir;
  if(MCtruth)
    massfitsdir = ".";
  else
    massfitsdir = ".";
  TString effdir = ".";

  TString crossdir = ".";
  TString corryielddir = ".";
  TString fracdir = ".";
  
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
  Bool_t isCutSymm[nCutVars];
  
  cout << "\n________________________________ AXES ___________________________________ \n\n Axis Name" <<setw(20) <<"DataAxis" << setw(20) << "MCGenAxis" << setw(20) << " MCRecoAxis\n"<< "_________________________________________________________________________\n" <<endl;
  
  for(Int_t iCutVar=0; iCutVar<nCutVars; iCutVar++) {
    AxesNames[iCutVar] = axesnames[iCutVar];
    dataAxesNo[iCutVar] = axesno[iCutVar];
    MCgenAxesNo[iCutVar] = axesno[iCutVar+nCutVars];
    MCcutAxesNo[iCutVar] = axesno[iCutVar+2*nCutVars];
    if(axesno[iCutVar+3*nCutVars]==0) isCutSymm[iCutVar]=kFALSE;
    else isCutSymm[iCutVar]=kTRUE;
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
  vector< vector<double> > cutsets;
  cutsets.push_back(cutset1);
  cutsets.push_back(cutset2);
  cutsets.push_back(cutset3);

  Int_t nPtBins= cutsets[0].size()/(2*nCutVars);
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
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      for(Int_t iCutVar=0; iCutVar<nCutVars; iCutVar++) {
        cutlowset[iSet][iPt][iCutVar] = cutsets[iSet][(2*iCutVar)+(iPt*nCutVars*2)];
        cuthighset[iSet][iPt][iCutVar] = cutsets[iSet][(2*iCutVar+1)+(iPt*nCutVars*2)];      
      }
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
  Double_t **sigmas = new Double_t*[nSets];
  
  for(Int_t iSet=0; iSet<nSets; iSet++) {
    means[iSet] = new Double_t[nPtBins];
    sigmas[iSet] = new Double_t[nPtBins];
    
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      if(fixmeanfromPDG) {
        means[iSet][iPt] = massD;
      }
      else {
        means[iSet][iPt] = -massD;
      }
      sigmas[iSet][iPt] = -0.008;
    }
  }
  
  TFile *massfitsfile = TFile::Open("massfitsMCtruth/RawYields.root","READ");
  if(massfitsfile) {
    for(Int_t iSet = 0; iSet<nSets; iSet++) {
      TH1F* hSigma = (TH1F*)massfitsfile->Get("hRawYieldSigmas_Set2"); //all the sigmas set/fixed to the MC values for the set 2
      if(fixsigmafromMC) {
        for(Int_t iPt=0; iPt<nPtBins; iPt++) {
          sigmas[iSet][iPt] = hSigma->GetBinContent(iPt+1);
        }
      }
      else {
        for(Int_t iPt=0; iPt<nPtBins; iPt++) {
          sigmas[iSet][iPt] = -hSigma->GetBinContent(iPt+1);
        } 
      }
    }
  }
  
  //________________________________________________________________________________________________________________
  //Analysis    
  
  AliHFCutVarFDsubAnalysisManagerDplus *AnalysisManagerDplus = new AliHFCutVarFDsubAnalysisManagerDplus();  
  AnalysisManagerDplus->SetPID(PID,3);
  Int_t loadTH = AnalysisManagerDplus->GetTHnSparses(MCfile,datafile,MCdir,datadir,MClist,datalist,MCtruth);

  if(loadTH>0)
    return 1;

  AnalysisManagerDplus->GetAxes(dataAxesNo,MCgenAxesNo,MCcutAxesNo,AxesNames,nCutVars,isCutSymm);
  AnalysisManagerDplus->GetCuts(cutlowset,cuthighset,means,sigmas,Rebin,funsig,funback,1.72,2.05,nSets,nPtBins,nCutVars);
    
  AnalysisManagerDplus->DrawDistributions(distdir);
  AnalysisManagerDplus->GetXaxisInformation();
  AnalysisManagerDplus->GetEfficiencies();
  AnalysisManagerDplus->DrawEfficiencies(effdir);
  AnalysisManagerDplus->GetRawYields(kTRUE,massfitsdir);

  Int_t method = 0;
  if(isIncentre)
    method = 1;
  Bool_t min = AnalysisManagerDplus->Minimise(method,10,kTRUE,0.,10000);
  AnalysisManagerDplus->DrawLines(linesdir,isIncentre);

  TString strmethod = "Min";
  if(isIncentre)
    strmethod = "Inc";
  
  //_____________________________________________________________________________________________________________________
  //Plots
  TCanvas *cCorrPrompt = new TCanvas("cCorrPrompt","",10,10,1920,1080);
  cCorrPrompt->Clear();
  TCanvas *cCorrFD = new TCanvas("cCorrFD","",10,10,1920,1080);
  cCorrFD->Clear();
  
  TH1F *hCorrPrompt = (TH1F*)AnalysisManagerDplus->GetYieldsPrompt();
  TH1F *hCorrFD = (TH1F*)AnalysisManagerDplus->GetYieldsFD();
  hCorrPrompt->SetStats(kFALSE);
  hCorrPrompt->SetDirectory(0);
  hCorrPrompt->GetYaxis()->SetTitle("#frac{dN_{prompt}}{d#it{p}_{T}} (c/GeV)");
  hCorrPrompt->SetStats(kFALSE);
  hCorrPrompt->SetLineColor(kBlack);
  hCorrPrompt->SetMarkerStyle(20);
  hCorrPrompt->SetMarkerSize(1);
  hCorrPrompt->SetMarkerColor(kBlack);
  hCorrFD->GetYaxis()->SetTitle("#frac{dN_{FD}}{d#it{p}_{T}} (c/GeV)");
  hCorrFD->SetStats(kFALSE);
  hCorrFD->SetDirectory(0);
  hCorrFD->SetLineColor(kBlack);
  hCorrFD->SetMarkerStyle(20);
  hCorrFD->SetMarkerSize(1);
  hCorrFD->SetMarkerColor(kBlack);

  TLegend *lPrompt = new TLegend(0.45,0.68,0.80,0.87);
  lPrompt->SetBorderSize(0);
  lPrompt->SetFillColor(kWhite);
  lPrompt->SetTextSizePixels(1000);
  lPrompt->SetTextFont(42);
  lPrompt->SetTextSize(0.05);
  
  TLegend *lFD = new TLegend(0.45,0.68,0.80,0.87);
  lFD->SetBorderSize(0);
  lFD->SetFillColor(kWhite);
  lFD->SetTextSizePixels(1000);
  lFD->SetTextFont(42);
  lFD->SetTextSize(0.05);
  
  lPrompt->AddEntry(hCorrPrompt,"Cut Variation Method - Prompt","lpe");    
  lFD->AddEntry(hCorrFD,"Cut Variation Method - FD","lpe");
    
  cCorrPrompt->Clear();
  cCorrPrompt->SetLogy();
  hCorrPrompt->Draw("E1");
  lPrompt->Draw("same");

  hCorrPrompt->SaveAs(Form("%s/CorrYieldsPrompt_%s.root",corryielddir.Data(),strmethod.Data()));

  cCorrFD->Clear();
  cCorrFD->SetLogy();
  hCorrFD->Draw("E1");
  lFD->Draw("same");

  hCorrFD->SaveAs(Form("%s/CorrYieldsFD_%s.root",corryielddir.Data(),strmethod.Data()));                      

  TFile corryieldfile(Form("%s/CorrYields_%s.root",corryielddir.Data(),strmethod.Data()),"RECREATE");
  hCorrPrompt->Write();
  hCorrFD->Write();
  corryieldfile.Close();
  
  if(MCtruth) {
    TFile file(MCfile.Data(),"UPDATE");
    TDirectoryFile* MCdir = (TDirectoryFile*)file.Get("PWG3_D2H_InvMassDplus");
    TList* MClist = (TList*)MCdir->Get("coutputDplus_CutVarpPbMC0100");
    THnSparseF* hMCGenAccPrompt = (THnSparseF*)MClist->FindObject("hMCAccPrompt");
    THnSparseF* hMCGenAccFD = (THnSparseF*)MClist->FindObject("hMCAccBFeed");
    TH1F* hPtMCGenAccPrompt = (TH1F*)hMCGenAccPrompt->Projection(0);
    TH1F* hPtMCGenAccFD = (TH1F*)hMCGenAccFD->Projection(0);
    hPtMCGenAccPrompt->Sumw2();
    hPtMCGenAccFD->Sumw2();
    
    TArrayD* ptarray = (TArrayD*)hCorrPrompt->GetXaxis()->GetXbins();
    Double_t* PtLims = (Double_t*)ptarray->GetArray();

    TH1F* hPtMCGenAccPromptReb = (TH1F*)hPtMCGenAccPrompt->Rebin(nPtBins,"hPtMCGenAccPromptReb",PtLims);
    TH1F* hPtMCGenAccFDReb = (TH1F*)hPtMCGenAccFD->Rebin(nPtBins,"hPtMCGenAccFDReb",PtLims);
    hPtMCGenAccPromptReb->Sumw2();
    hPtMCGenAccFDReb->Sumw2();

    hPtMCGenAccPromptReb->SetDirectory(0);
    hPtMCGenAccPromptReb->SetMarkerStyle(20);
    hPtMCGenAccPromptReb->SetMarkerColor(kBlue);
    hPtMCGenAccPromptReb->SetLineColor(kBlue);
    hPtMCGenAccFDReb->SetDirectory(0);
    hPtMCGenAccFDReb->SetMarkerStyle(20);
    hPtMCGenAccFDReb->SetMarkerColor(kBlue);
    hPtMCGenAccFDReb->SetLineColor(kBlue);

    cCorrPrompt->cd();
    hPtMCGenAccPromptReb->Draw("E1same");
    cCorrPrompt->SaveAs(Form("%s/CorrYieldsPrompt_%s.eps",corryielddir.Data(),strmethod.Data()));
    
    cCorrFD->cd();
    hPtMCGenAccFDReb->Draw("E1same");
    cCorrFD->SaveAs(Form("%s/CorrYieldsFD_%s.eps",corryielddir.Data(),strmethod.Data()));

    delete cCorrPrompt;
    delete cCorrFD;
  }
  else {
    cCorrPrompt->SaveAs(Form("%s/CorrYieldsPrompt_%s.eps",corryielddir.Data(),strmethod.Data()));
    cCorrFD->SaveAs(Form("%s/CorrYieldsFD_%s.eps",corryielddir.Data(),strmethod.Data()));
    delete cCorrPrompt;
    delete cCorrFD;

    TCanvas *cCrossSectionPrompt = new TCanvas("cCrossSectionPrompt","cCrossSectionPrompt",10,10,80,800);
    cCrossSectionPrompt->Clear();
    TCanvas *cCrossSectionFD = new TCanvas("cCrossSectionFD","cCrossSectionFD",10,10,800,800);
    cCrossSectionFD->Clear();
    
    TAxis* PtAxis = (TAxis*)hCorrPrompt->GetXaxis();
    TAxis* YAxis = (TAxis*)hCorrPrompt->GetYaxis();
    TArrayD* PtBinsArray = (TArrayD*)PtAxis->GetXbins();
    Double_t* PtLimsArray = (Double_t*)PtBinsArray->GetArray();
    Int_t nPtBins = PtBinsArray->GetSize()-1;
    
    TH1F *hCrossSecPrompt=(TH1F*)AnalysisManagerDplus->GetCrossSecPrompt(AccFileName,"hPtGenLimAcc","hPtGenAcc",0.0913,2.09);
    TH1F *hCrossSecFD=(TH1F*)AnalysisManagerDplus->GetCrossSecFD(AccFileName,"hPtGenLimAcc","hPtGenAcc",0.0913,2.09);
    
    hCrossSecPrompt->SetStats(kFALSE);
    hCrossSecPrompt->SetLineColor(kBlack);
    hCrossSecPrompt->SetMarkerStyle(20);
    hCrossSecPrompt->SetMarkerSize(1);
    hCrossSecPrompt->SetMarkerColor(kBlack);
    hCrossSecFD->SetStats(kFALSE);
    hCrossSecFD->SetLineColor(kBlack);
    hCrossSecFD->SetMarkerStyle(20);
    hCrossSecFD->SetMarkerSize(1);
    hCrossSecFD->SetMarkerColor(kBlack);
    
    cCrossSectionPrompt->Clear();
    cCrossSectionPrompt->SetLogy();
    hCrossSecPrompt->Draw("E1");
    lPrompt->Draw("same");
    
    hCrossSecPrompt->SaveAs(Form("%s/CrossPrompt_%s.root",crossdir.Data(),strmethod.Data()));
    
    cCrossSectionFD->Clear();
    cCrossSectionFD->SetLogy();
    hCrossSecFD->Draw("E1");
    lFD->Draw("same");
    
    hCrossSecFD->SaveAs(Form("%s/CrossFD_%s.root",crossdir.Data(),strmethod.Data()));                      

    TFile crossecfile(Form("%s/CrossSec_%s.root",crossdir.Data(),strmethod.Data()),"RECREATE");
    hCrossSecPrompt->Write();
    hCrossSecFD->Write();
    crossecfile.Close();

    delete hCrossSecPrompt;
    delete hCrossSecFD;
    delete lPrompt;
    delete lFD;
    delete cCrossSectionPrompt;
    delete cCrossSectionFD;
    
    //fraction - residuals - pulls
    TFile outfracfile(Form("%s/PromptFraction_Res_Pulls_%s.root",fracdir.Data(),strmethod.Data()),"RECREATE");
    TList* lFPrompt = (TList*)AnalysisManagerDplus->GetPromptFraction(); 
    TList* lFPromptRaw = (TList*)AnalysisManagerDplus->GetPromptFractionRaw(); 
    TList* lResiduals = (TList*)AnalysisManagerDplus->GetResiduals();
    TList* lPulls = (TList*)AnalysisManagerDplus->GetPulls();
    
    TH1F** hFPrompt = new TH1F*[nPtBins];
    TH1F** hFPromptRaw = new TH1F*[nPtBins];
    TH1F** hResiduals = new TH1F*[nPtBins];
    TH1F** hPulls = new TH1F*[nPtBins];
    
    for(Int_t iPt=0; iPt<nPtBins; ++iPt) {
      hFPrompt[iPt]=(TH1F*)lFPrompt->At(iPt);
      hFPrompt[iPt]->SetName(Form("hFPrompt_%d",iPt));
      hFPrompt[iPt]->Write();
      hFPromptRaw[iPt]=(TH1F*)lFPromptRaw->At(iPt);
      hFPromptRaw[iPt]->SetName(Form("hFPromptRaw_%d",iPt));
      hFPromptRaw[iPt]->Write();
      hResiduals[iPt]=(TH1F*)lResiduals->At(iPt);
      hResiduals[iPt]->SetName(Form("hResiduals_%d",iPt));
      hResiduals[iPt]->Write(); 
      hPulls[iPt]=(TH1F*)lPulls->At(iPt);
      hPulls[iPt]->SetName(Form("hPulls_%d",iPt));
      hPulls[iPt]->Write();
    }
    outfracfile.Close();

    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      delete hFPrompt[iPt];
      delete hFPromptRaw[iPt];
      delete hPulls[iPt];
      delete hResiduals[iPt];
    }
    delete[] hFPrompt;
    delete[] hFPromptRaw;
    delete[] hPulls;
    delete[] hResiduals;
  }

  for(Int_t iSet=0; iSet<nSets; iSet++) {
    delete means[iSet];
    delete sigmas[iSet];
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      delete cutlowset[iSet][iPt];
      delete cuthighset[iSet][iPt];
    }   
  }
  delete[] means;
  delete[] sigmas;
  delete[] cutlowset;
  delete[] cuthighset;
  
  delete AnalysisManagerDplus;
  
  return 0;
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
