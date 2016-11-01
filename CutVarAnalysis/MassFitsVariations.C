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
#include <TGaxis.h>

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

enum MassSigma{kFixed,kLow,kHigh};
enum MassMean{kFree,kFix};
enum Method{kMinimisation,kIncentre};

void MassFitsVariations(Int_t iMethod=0);
void ExtimateUncertainty(TString PromptOrFD="Prompt", TString Method="Min");
void ReadAxes(TString FileName, vector<string> &axesanmes, vector<int> &axesno);
void ReadSet(TString FileName, vector<string> &varnames, vector<double> &cutset);
Double_t GetMin(TH1F* h);
Double_t GetMax(TH1F* h);

void MassFitsVariations(Int_t iMethod) {

 //________________________________________________________________________________________________________________
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
  
  Double_t **meansfromPDG = new Double_t*[nSets];
  Double_t **sigmasfromMC = new Double_t*[nSets];
  
  for(Int_t iSet=0; iSet<nSets; iSet++) {
    meansfromPDG[iSet] = new Double_t[nPtBins];
    sigmasfromMC[iSet] = new Double_t[nPtBins];    
  }

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    meansfromPDG[0][iPt] = massD;
    meansfromPDG[1][iPt] = massD;
    meansfromPDG[2][iPt] = massD;
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
  
  TFile corrrefminfile("../../corryields/CorrYields_Min.root","UPDATE");
  TH1F* hCorrYieldPromptRefMin = (TH1F*)corrrefminfile.Get("hCorrYieldPrompt"); 
  TH1F* hCorrYieldFDRefMin = (TH1F*)corrrefminfile.Get("hCorrYieldFD"); 
  hCorrYieldPromptRefMin->SetDirectory(0);
  hCorrYieldFDRefMin->SetDirectory(0);
  corrrefminfile.Close();
  
  TH1F** hRawYieldRef = new TH1F*[nSets];
  TFile rawreffile("../../massfits/RawYields.root","UPDATE");
  for(Int_t iSet=0; iSet<nSets; iSet++) {
    hRawYieldRef[iSet] = (TH1F*)rawreffile.Get(Form("hRawYields_Set%d",iSet+1)); 
    hRawYieldRef[iSet]->SetDirectory(0);
  }
  rawreffile.Close();

  //________________________________________________________________________________________________________________
  //Analysis    
  const Int_t nReb = 5;
  Int_t rebin[nReb] = {2,4,5,6,8};
  const Int_t nSigma = 3;
  Int_t sigma[nSigma] = {kFixed,kLow,kHigh};
  Double_t **sigmas = new Double_t*[nSets];
  const Int_t nMean = 2;
  Int_t mean[nMean] = {kFree,kFix};
  Double_t **means = new Double_t*[nSets];
  const Int_t nBkgFunc = 2;
  Int_t funback[nBkgFunc] = {AliHFMassFitter::kExpo,AliHFMassFitter::kPol2};
  Int_t funsig = AliHFMassFitter::kGaus;
  const Int_t nRanges = 4;
  Double_t min[nRanges] = {1.72,1.73,1.74,1.75};
  Double_t max[nRanges] = {2.05,2.04,2.03,2.02};  

  Int_t nTrials = nReb*nSigma*nMean*nRanges*nRanges*nBkgFunc;
  
  TH1F*** hRawYields = new TH1F**[nSets];
  TH1F*** hRawYieldsVsTrial = new TH1F**[nSets];
  for(Int_t iSet=0; iSet<nSets; iSet++) {
    hRawYields[iSet] = new TH1F*[nPtBins];
    hRawYieldsVsTrial[iSet] = new TH1F*[nPtBins];
  }
  
  Int_t nBins = 40;
  Double_t MinsPrompt[nPtBins];
  Double_t MaxsPrompt[nPtBins];
  Double_t MinsFD[nPtBins];
  Double_t MaxsFD[nPtBins];
  
  Double_t MinsRaw[nSets][nPtBins];
  Double_t MaxsRaw[nSets][nPtBins];

  TString MethodName[2] = {"Min","Inc"};
  
  TH1F** hCorrYieldsPrompt = new TH1F*[nPtBins];
  TH1F** hCorrYieldsFD = new TH1F*[nPtBins];
  TH1F** hCorrYieldsPromptVsTrial = new TH1F*[nPtBins];
  TH1F** hCorrYieldsFDVsTrial = new TH1F*[nPtBins];

  for(Int_t iPt = 0; iPt<nPtBins; iPt++) {
    MinsPrompt[iPt] = hCorrYieldPromptRefMin->GetBinContent(iPt+1)-hCorrYieldPromptRefMin->GetBinContent(iPt+1)*0.5;
    MaxsPrompt[iPt] = hCorrYieldPromptRefMin->GetBinContent(iPt+1)+hCorrYieldPromptRefMin->GetBinContent(iPt+1)*0.5;
    MinsFD[iPt] = hCorrYieldFDRefMin->GetBinContent(iPt+1)-hCorrYieldFDRefMin->GetBinContent(iPt+1)*1.;
    MaxsFD[iPt] = hCorrYieldFDRefMin->GetBinContent(iPt+1)+hCorrYieldFDRefMin->GetBinContent(iPt+1)*1.;

    hCorrYieldsPrompt[iPt] = new TH1F(Form("hCorrYieldsPrompt%s_Pt%d",MethodName[iMethod].Data(),iPt),"",nBins,MinsPrompt[iPt],MaxsPrompt[iPt]);
    hCorrYieldsFD[iPt] = new TH1F(Form("hCorrYieldsFD%s_Pt%d",MethodName[iMethod].Data(),iPt),"",nBins,MinsFD[iPt],MaxsFD[iPt]);
    hCorrYieldsPromptVsTrial[iPt] = new TH1F(Form("hCorrYieldsPromptVsTrial%s_Pt%d",MethodName[iMethod].Data(),iPt),"",nTrials,0.5,nTrials+0.5);
    hCorrYieldsFDVsTrial[iPt] = new TH1F(Form("hCorrYieldsFDVsTrial%s_Pt%d",MethodName[iMethod].Data(),iPt),"",nTrials,0.5,nTrials+0.5); 
    hCorrYieldsPrompt[iPt]->GetXaxis()->SetTitle("N_{prompt}");
    hCorrYieldsPrompt[iPt]->GetYaxis()->SetTitle("Entries");
    hCorrYieldsFD[iPt]->GetXaxis()->SetTitle("N_{feed-down}");
    hCorrYieldsFD[iPt]->GetYaxis()->SetTitle("Entries");
    hCorrYieldsPrompt[iPt]->GetXaxis()->SetNdivisions(10,10,0);
    hCorrYieldsPrompt[iPt]->GetYaxis()->SetNdivisions(6,10,0);
    hCorrYieldsFD[iPt]->GetXaxis()->SetNdivisions(10,10,0);
    hCorrYieldsFD[iPt]->GetYaxis()->SetNdivisions(6,10,0);
    hCorrYieldsPromptVsTrial[iPt]->GetYaxis()->SetTitle("N_{prompt}");
    hCorrYieldsPromptVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hCorrYieldsFDVsTrial[iPt]->GetYaxis()->SetTitle("N_{feed-down}");
    hCorrYieldsFDVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hCorrYieldsPromptVsTrial[iPt]->GetXaxis()->SetNdivisions(10,10,0);
    hCorrYieldsPromptVsTrial[iPt]->GetYaxis()->SetNdivisions(20,10,0);
    hCorrYieldsFDVsTrial[iPt]->GetXaxis()->SetNdivisions(10,10,0);
    hCorrYieldsFDVsTrial[iPt]->GetYaxis()->SetNdivisions(20,10,0);

    for(Int_t iSet=0; iSet<nSets; iSet++) {
      MinsRaw[iSet][iPt] = hRawYieldRef[iSet]->GetBinContent(iPt+1)-hRawYieldRef[iSet]->GetBinContent(iPt+1)*0.5;
      MaxsRaw[iSet][iPt] = hRawYieldRef[iSet]->GetBinContent(iPt+1)+hRawYieldRef[iSet]->GetBinContent(iPt+1)*0.5;

      hRawYields[iSet][iPt] = new TH1F(Form("hRawYields%d_Pt%d",iSet+1,iPt),"",nBins,MinsRaw[iSet][iPt],MaxsRaw[iSet][iPt]);
      hRawYieldsVsTrial[iSet][iPt] = new TH1F(Form("hRawYieldsVsTrial%d_Pt%d",iSet+1,iPt),"",nTrials,0.5,nTrials+0.5);
      hRawYields[iSet][iPt]->SetDirectory(0);
      hRawYields[iSet][iPt]->SetDirectory(0);
      hRawYields[iSet][iPt]->SetTitle(Form("Set %d",iSet+1)); 
      hRawYields[iSet][iPt]->GetXaxis()->SetTitle("Raw Yield"); 
      hRawYields[iSet][iPt]->GetYaxis()->SetTitle("Entries"); 
      hRawYields[iSet][iPt]->GetXaxis()->SetNdivisions(10,10,0);
      hRawYields[iSet][iPt]->GetYaxis()->SetNdivisions(6,10,0);
      hRawYieldsVsTrial[iSet][iPt]->SetTitle(Form("Set %d",iSet+1)); 
      hRawYieldsVsTrial[iSet][iPt]->GetYaxis()->SetTitle("Raw Yield"); 
      hRawYieldsVsTrial[iSet][iPt]->GetXaxis()->SetTitle("Trial #"); 
      hRawYieldsVsTrial[iSet][iPt]->GetXaxis()->SetNdivisions(10,10,0);
      hRawYieldsVsTrial[iSet][iPt]->GetYaxis()->SetNdivisions(20,10,0);
    }
  }
  
  Int_t TrialCounter = 0;

  AliHFCutVarFDsubAnalysisManagerDplus *AnalysisManagerDplus = new AliHFCutVarFDsubAnalysisManagerDplus();  
  AnalysisManagerDplus->SetPID(kTRUE,3);
  Int_t loadTH = AnalysisManagerDplus->GetTHnSparses(MCfile,datafile,MCdir,datadir,MClist,datalist,kFALSE);
  
  if(loadTH>0)
    return;

  for(Int_t iSigma=0; iSigma<nSigma; iSigma++) {
    for(Int_t iMean=0; iMean<nMean; iMean++) {
      for(Int_t iReb=0; iReb<nReb; iReb++) {
        for(Int_t iLowRange=0; iLowRange<nRanges; iLowRange++) {
          for(Int_t iHighRange=0; iHighRange<nRanges; iHighRange++) {
            for(Int_t iFun=0; iFun<nBkgFunc; iFun++) {
              for(Int_t iSet=0; iSet<nSets; iSet++) {
                means[iSet] = new Double_t[nPtBins];
                
                for(Int_t iPt=0; iPt<nPtBins; iPt++) {
                  if(iMean==kFix) {
                    means[iSet][iPt] = meansfromPDG[iSet][iPt];
                  }
                  else {
                    means[iSet][iPt] = -meansfromPDG[iSet][iPt];
                  }
                }
              }
              
              for(Int_t iSet=0; iSet<nSets; iSet++) {
                sigmas[iSet] = new Double_t[nPtBins];
                
                for(Int_t iPt=0; iPt<nPtBins; iPt++) {
                  if(iSigma==kFixed) {
                    sigmas[iSet][iPt] = sigmasfromMC[iSet][iPt];
                  }
                  else if(iSigma==kLow) {
                    sigmas[iSet][iPt] = sigmasfromMC[iSet][iPt]-0.15*sigmasfromMC[iSet][iPt];
                  }
                  else {
                    sigmas[iSet][iPt] = sigmasfromMC[iSet][iPt]+0.15*sigmasfromMC[iSet][iPt];
                  }
                }
              }

              cout << "\n\n\\_______________________________________TRIAL NO "<< TrialCounter+1 << "/" << nTrials << " - " << MethodName[iMethod].Data() << " ________________________________________\\\n\n"<<endl;
              
              AnalysisManagerDplus->GetAxes(dataAxesNo,MCgenAxesNo,MCcutAxesNo,AxesNames,nCutVars);
              AnalysisManagerDplus->GetCuts(cutlowset,cuthighset,means,sigmas,rebin[iReb],funsig,funback[iFun],min[iLowRange],max[iHighRange],nSets,nPtBins,nCutVars);
              AnalysisManagerDplus->GetXaxisInformation();
              AnalysisManagerDplus->GetEfficiencies();
              AnalysisManagerDplus->GetRawYields(kFALSE);
              AnalysisManagerDplus->Minimise(iMethod,10,kTRUE,0.,10000);
              
              TH1F *hCorrPrompt = (TH1F*)AnalysisManagerDplus->GetYieldsPrompt();
              TH1F *hCorrFD = (TH1F*)AnalysisManagerDplus->GetYieldsFD();
              
              TFile rawyieldsfile("RawYields.root","UPDATE");
              TH1F* hPtRaw1 = (TH1F*)rawyieldsfile.Get("hRawYields_Set1");
              TH1F* hPtRaw2 = (TH1F*)rawyieldsfile.Get("hRawYields_Set2");
              TH1F* hPtRaw3 = (TH1F*)rawyieldsfile.Get("hRawYields_Set3");
              hPtRaw1->SetDirectory(0);
              hPtRaw2->SetDirectory(0);
              hPtRaw3->SetDirectory(0);
              rawyieldsfile.Close();
              
              for(Int_t iPt = 0; iPt<nPtBins; iPt++) {
                hCorrYieldsPrompt[iPt]->Fill(hCorrPrompt->GetBinContent(iPt+1));
                hCorrYieldsFD[iPt]->Fill(hCorrFD->GetBinContent(iPt+1));
                hCorrYieldsPromptVsTrial[iPt]->SetBinContent(TrialCounter+1,hCorrPrompt->GetBinContent(iPt+1));
                hCorrYieldsFDVsTrial[iPt]->SetBinContent(TrialCounter+1,hCorrFD->GetBinContent(iPt+1));      
                hCorrYieldsPromptVsTrial[iPt]->SetBinError(TrialCounter+1,hCorrPrompt->GetBinError(iPt+1));
                hCorrYieldsFDVsTrial[iPt]->SetBinError(TrialCounter+1,hCorrFD->GetBinError(iPt+1));      

                hRawYields[0][iPt]->Fill(hPtRaw1->GetBinContent(iPt+1));
                hRawYields[1][iPt]->Fill(hPtRaw2->GetBinContent(iPt+1));
                hRawYields[2][iPt]->Fill(hPtRaw3->GetBinContent(iPt+1));
                hRawYieldsVsTrial[0][iPt]->SetBinContent(TrialCounter+1,hPtRaw1->GetBinContent(iPt+1));
                hRawYieldsVsTrial[1][iPt]->SetBinContent(TrialCounter+1,hPtRaw2->GetBinContent(iPt+1));
                hRawYieldsVsTrial[2][iPt]->SetBinContent(TrialCounter+1,hPtRaw3->GetBinContent(iPt+1));
                hRawYieldsVsTrial[0][iPt]->SetBinError(TrialCounter+1,hPtRaw1->GetBinError(iPt+1));
                hRawYieldsVsTrial[1][iPt]->SetBinError(TrialCounter+1,hPtRaw2->GetBinError(iPt+1));
                hRawYieldsVsTrial[2][iPt]->SetBinError(TrialCounter+1,hPtRaw3->GetBinError(iPt+1));
              }
              TrialCounter++;
            }
          }
        }
      }
    }
  }
  
  //_____________________________________________________________________________________________________________________
  //Plots
  gStyle->SetPadLeftMargin(1.3);
  gStyle->SetPadBottomMargin(1.1);
  gStyle->SetTitleSize(0.045,"xy");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetTitleOffset(1.,"x");
    
  TCanvas *cPrompt = new TCanvas("cPrompt","Prompt",1200,900);
  cPrompt->Divide(nPtBins/2+1,2);
  TCanvas *cFD = new TCanvas("cFD","FD",1200,900);
  cFD->Divide(nPtBins/2+1,2);
  TCanvas *cPromptVsTrial = new TCanvas("cPromptVsTrial","Prompt",1200,900);
  cPromptVsTrial->Divide(nPtBins/2+1,2);
  TCanvas *cFDVsTrial = new TCanvas("cFDVsTrial","FD",1200,900);
  cFDVsTrial->Divide(nPtBins/2+1,2);
  TCanvas** cRaw = new TCanvas*[nPtBins];
  TCanvas** cRawVsTrial = new TCanvas*[nPtBins];
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    cPrompt->cd(iPt+1);
    hCorrYieldsPrompt[iPt]->Draw();
    cFD->cd(iPt+1);
    hCorrYieldsFD[iPt]->Draw();

    hCorrYieldsPromptVsTrial[iPt]->GetYaxis()->SetRangeUser(hCorrYieldsPromptVsTrial[iPt]->GetMinimum()*0.5,hCorrYieldsPromptVsTrial[iPt]->GetMaximum()*1.7);
    hCorrYieldsFDVsTrial[iPt]->GetYaxis()->SetRangeUser(hCorrYieldsFDVsTrial[iPt]->GetMinimum()*0.5,hCorrYieldsFDVsTrial[iPt]->GetMaximum()*1.7);

    cPromptVsTrial->cd(iPt+1);
    hCorrYieldsPromptVsTrial[iPt]->Draw("E");
    cFDVsTrial->cd(iPt+1);
    hCorrYieldsFDVsTrial[iPt]->Draw("E");
    
    cRaw[iPt] = new TCanvas(Form("cRaw_Pt%d",iPt),"Raw Yields",1200,900);
    cRaw[iPt]->Divide(3,1);
    cRawVsTrial[iPt] = new TCanvas(Form("cRawVsTrial_Pt%d",iPt),"Raw Yields",1200,900);
    cRawVsTrial[iPt]->Divide(3,1);
    
    for(Int_t iSet=0; iSet<nSets; iSet++) {
      hRawYieldsVsTrial[iSet][iPt]->GetYaxis()->SetRangeUser(hRawYieldsVsTrial[iSet][iPt]->GetMinimum()*0.5,hRawYieldsVsTrial[iSet][iPt]->GetMaximum()*1.5);
      cRaw[iPt]->cd(iSet+1);
      hRawYields[iSet][iPt]->Draw();
      cRawVsTrial[iPt]->cd(iSet+1);
      hRawYieldsVsTrial[iSet][iPt]->Draw();
    }
  }

  TFile outfile(Form("MassFitVariations_%s.root",MethodName[iMethod].Data()),"RECREATE");
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    hCorrYieldsPrompt[iPt]->Write();
    hCorrYieldsFD[iPt]->Write();
    hCorrYieldsPromptVsTrial[iPt]->Write();
    hCorrYieldsFDVsTrial[iPt]->Write();
    hRawYields[0][iPt]->Write();
    hRawYields[1][iPt]->Write();
    hRawYields[2][iPt]->Write();
    hRawYieldsVsTrial[0][iPt]->Write();
    hRawYieldsVsTrial[1][iPt]->Write();
    hRawYieldsVsTrial[2][iPt]->Write();
    cRaw[iPt]->Write();
    cRawVsTrial[iPt]->Write();
  }
  cPrompt->Write();
  cFD->Write();
  cPromptVsTrial->Write();
  cFDVsTrial->Write();
  outfile.Close();

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    delete hCorrYieldsPrompt[iPt];
    delete hCorrYieldsFD[iPt];
    delete hCorrYieldsPromptVsTrial[iPt];
    delete hCorrYieldsFDVsTrial[iPt];
  }
  
  delete[] hCorrYieldsPrompt;
  delete[] hCorrYieldsFD;
  delete[] hCorrYieldsPromptVsTrial;
  delete[] hCorrYieldsFDVsTrial;
  delete cPrompt;
  delete cFD;
  delete cPromptVsTrial;
  delete cFDVsTrial;
  
  for(Int_t iSet=0; iSet<nSets; iSet++) {
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      delete hRawYields[iSet][iPt];
      delete hRawYieldsVsTrial[iSet][iPt];
    }
    delete[] hRawYields[iSet];
    delete[] hRawYieldsVsTrial[iSet];
  }
  delete[] hRawYields;
  delete[] hRawYieldsVsTrial;

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    delete cRaw[iPt];
    delete cRawVsTrial[iPt];
  }
  delete[] cRaw;
  delete[] cRawVsTrial;
  
}

void ExtimateUncertainty(TString PromptOrFD, TString Method) {

  const Int_t nPtBins=7;
  const Int_t nPtLims=nPtBins+1;
  Double_t PtLims[nPtLims]={2,3,4,5,6,8,12,16};
  
  const Int_t nSets=3;
  
  gStyle->SetTitleSize(0.045,"xy");
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);

  TFile corrreffile(Form("../../corryields/CorrYields_%s.root",Method.Data()),"UPDATE");
  TH1F* hCorrYieldRef = (TH1F*)corrreffile.Get(Form("hCorrYield%s",PromptOrFD.Data())); 
  hCorrYieldRef->SetDirectory(0);
  corrreffile.Close();
  TLine** line = new TLine*[nPtBins];
  TLine** lineSigmaPlus = new TLine*[nPtBins];
  TLine** lineSigmaMinus = new TLine*[nPtBins];
  TLine** linevstrial = new TLine*[nPtBins];

  TH1F** hRawYieldRef = new TH1F*[nSets];
  TLine*** linerawvstrial = new TLine**[nSets];
  TFile rawreffile("../../massfits/RawYields.root","UPDATE");
  for(Int_t iSet=0; iSet<nSets; iSet++) {
    hRawYieldRef[iSet] = (TH1F*)rawreffile.Get(Form("hRawYields_Set%d",iSet+1)); 
    hRawYieldRef[iSet]->SetDirectory(0);
    linerawvstrial[iSet] = new TLine*[nPtBins];
  }
  rawreffile.Close();

  TH1F** hCorrYields = new TH1F*[nPtBins];
  TH1F** hCorrYieldsVsTrial = new TH1F*[nPtBins];
  TH1F*** hRawYields = new TH1F**[nSets];
  TH1F*** hRawYieldsVsTrial = new TH1F**[nSets];
  for(Int_t iSet=0; iSet<nSets; iSet++) {
    hRawYields[iSet] = new TH1F*[nPtBins];
    hRawYieldsVsTrial[iSet] = new TH1F*[nPtBins];
  }
  TCanvas** cCorrYields = new TCanvas*[nPtBins];
  TCanvas** cCorrYieldsVsTrial = new TCanvas*[nPtBins];
  TCanvas** cCorrYieldsDisp = new TCanvas*[nPtBins];
  TCanvas** cRawYields = new TCanvas*[nPtBins];
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.045);
  latex->SetTextFont(132);

  TLegend* l = new TLegend(0.18,0.8,0.48,0.89);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextSize(0.045);

  TString title;
  if(Method=="Min")
    title = "Analytic Minimisation";
  else
    title = "Incentre Minimisation";
  
  TFile infile(Form("MassFitVariations_%s.root",Method.Data()),"UPDATE");
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    hCorrYields[iPt]=(TH1F*)infile.Get(Form("hCorrYields%s%s_Pt%d",PromptOrFD.Data(),Method.Data(),iPt)); 
    hCorrYields[iPt]->SetDirectory(0);
    hCorrYields[iPt]->SetStats(kFALSE);
    hCorrYields[iPt]->GetXaxis()->SetTitleSize(0.05);
    hCorrYields[iPt]->GetYaxis()->SetTitleSize(0.05);
    hCorrYields[iPt]->GetXaxis()->SetLabelSize(0.05);
    hCorrYields[iPt]->GetYaxis()->SetLabelSize(0.05);
    hCorrYields[iPt]->GetYaxis()->SetTitleOffset(1.3);
    hCorrYields[iPt]->GetYaxis()->SetTitle("Entries");
    hCorrYields[iPt]->SetTitle(title.Data());
    hCorrYields[iPt]->SetFillStyle(3004);
    hCorrYields[iPt]->SetFillColor(kBlue);    
    hCorrYieldsVsTrial[iPt]=(TH1F*)infile.Get(Form("hCorrYields%sVsTrial%s_Pt%d",PromptOrFD.Data(),Method.Data(),iPt));
    hCorrYieldsVsTrial[iPt]->SetDirectory(0);
    hCorrYieldsVsTrial[iPt]->SetStats(kFALSE);
    hCorrYieldsVsTrial[iPt]->GetXaxis()->SetTitleSize(0.05);
    hCorrYieldsVsTrial[iPt]->GetYaxis()->SetTitleSize(0.05);
    hCorrYieldsVsTrial[iPt]->GetXaxis()->SetLabelSize(0.05);
    hCorrYieldsVsTrial[iPt]->GetYaxis()->SetLabelSize(0.05);
    hCorrYieldsVsTrial[iPt]->GetYaxis()->SetTitleOffset(1.3);
    hCorrYieldsVsTrial[iPt]->GetYaxis()->SetNdivisions(505);
    hCorrYieldsVsTrial[iPt]->GetXaxis()->SetNdivisions(505);
    hCorrYieldsVsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hCorrYields[iPt]->GetXaxis()->SetRangeUser(hCorrYieldsVsTrial[iPt]->GetMinimum()-hCorrYieldsVsTrial[iPt]->GetMinimum()*1.5,hCorrYieldsVsTrial[iPt]->GetMaximum()+hCorrYieldsVsTrial[iPt]->GetMaximum()*1.5);
    hCorrYields[iPt]->GetYaxis()->SetRangeUser(1,hCorrYields[iPt]->GetMaximum()*1.6);
    if(PromptOrFD=="FD") {
      hCorrYieldsVsTrial[iPt]->GetYaxis()->SetRangeUser(hCorrYieldRef->GetBinContent(iPt+1)*0.,hCorrYieldRef->GetBinContent(iPt+1)*2.5);
    } 
    TGaxis::SetMaxDigits(4);
 
    Double_t min = GetMin(hCorrYields[iPt]);
    Double_t max = GetMax(hCorrYields[iPt]);
    Double_t RMS = hCorrYields[iPt]->GetRMS();
    Double_t mean = hCorrYields[iPt]->GetMean();
    Double_t extimator = (max-min)/TMath::Sqrt(12);

    Double_t relRMS = RMS/mean*100;
    Double_t relext = extimator/mean*100;

    line[iPt] = new TLine(hCorrYieldRef->GetBinContent(iPt+1),1,hCorrYieldRef->GetBinContent(iPt+1),hCorrYields[iPt]->GetMaximum()/1.6);
    line[iPt]->SetLineColor(kRed);
    line[iPt]->SetLineWidth(2);

    lineSigmaPlus[iPt] = new TLine(320,hCorrYieldsVsTrial[iPt]->GetMinimum(),320,hCorrYieldsVsTrial[iPt]->GetMaximum());
    lineSigmaPlus[iPt]->SetLineColor(kBlack);
    lineSigmaPlus[iPt]->SetLineStyle(7);
    lineSigmaPlus[iPt]->SetLineWidth(2);

    lineSigmaMinus[iPt] = new TLine(640,hCorrYieldsVsTrial[iPt]->GetMinimum(),640,hCorrYieldsVsTrial[iPt]->GetMaximum());
    lineSigmaMinus[iPt]->SetLineColor(kBlack);
    lineSigmaMinus[iPt]->SetLineStyle(7);
    lineSigmaMinus[iPt]->SetLineWidth(2);

    linevstrial[iPt] = new TLine(0.5,hCorrYieldRef->GetBinContent(iPt+1),hCorrYieldsVsTrial[iPt]->GetNbinsX()-0.5,hCorrYieldRef->GetBinContent(iPt+1));
    linevstrial[iPt]->SetLineColor(kRed);
    linevstrial[iPt]->SetLineWidth(2);

    if(iPt==0)
      l->AddEntry(line[iPt],"Reference value","l");
    
    cCorrYields[iPt] = new TCanvas(Form("cCorrYieldsPt%d",iPt),"",10,10,1200,600);  
    cCorrYields[iPt]->Divide(2,1);
    cCorrYields[iPt]->cd(1);
    hCorrYieldsVsTrial[iPt]->Draw();
    linevstrial[iPt]->Draw("same");
    lineSigmaPlus[iPt]->Draw("same");
    lineSigmaMinus[iPt]->Draw("same");
    
    cCorrYields[iPt]->cd(2);
    hCorrYields[iPt]->Draw();
    line[iPt]->Draw("same");
    latex->DrawLatex(mean-0.37*mean,hCorrYields[iPt]->GetMaximum()*0.83,Form("mean = %0.f",mean));
    latex->DrawLatex(mean-0.37*mean,hCorrYields[iPt]->GetMaximum()*0.75,Form("RMS = %0.1f (%0.1f%)",RMS,relRMS));
    latex->DrawLatex(mean-0.37*mean,hCorrYields[iPt]->GetMaximum()*0.67,Form("(x_{max}-x_{min})/#sqrt{12} = %0.1f (%0.1f%)",extimator,relext));
    l->Draw("same");

    cCorrYieldsVsTrial[iPt] = new TCanvas(Form("cCorrYieldsVsTrialPt%d",iPt),"",10,10,800,800);  
    hCorrYieldsVsTrial[iPt]->Draw();
    lineSigmaPlus[iPt]->Draw("same");
    lineSigmaMinus[iPt]->Draw("same");
    linevstrial[iPt]->Draw("same");
    latex->DrawLatex(hCorrYieldsVsTrial[iPt]->GetNbinsX()*0.05,hCorrYieldsVsTrial[iPt]->GetMaximum()*0.92,"Mass width");
    latex->DrawLatex(hCorrYieldsVsTrial[iPt]->GetNbinsX()*0.05,hCorrYieldsVsTrial[iPt]->GetMaximum()*0.85,"fixed to MC");
    latex->DrawLatex(hCorrYieldsVsTrial[iPt]->GetNbinsX()*0.38,hCorrYieldsVsTrial[iPt]->GetMaximum()*0.92,"Mass width");
    latex->DrawLatex(hCorrYieldsVsTrial[iPt]->GetNbinsX()*0.38,hCorrYieldsVsTrial[iPt]->GetMaximum()*0.85,"fixed to MC");
    latex->DrawLatex(hCorrYieldsVsTrial[iPt]->GetNbinsX()*0.38,hCorrYieldsVsTrial[iPt]->GetMaximum()*0.78,"    -15%");
    latex->DrawLatex(hCorrYieldsVsTrial[iPt]->GetNbinsX()*0.71,hCorrYieldsVsTrial[iPt]->GetMaximum()*0.92,"Mass width");
    latex->DrawLatex(hCorrYieldsVsTrial[iPt]->GetNbinsX()*0.71,hCorrYieldsVsTrial[iPt]->GetMaximum()*0.85,"fixed to MC");
    latex->DrawLatex(hCorrYieldsVsTrial[iPt]->GetNbinsX()*0.71,hCorrYieldsVsTrial[iPt]->GetMaximum()*0.78,"    +15%");

    Double_t shift = 0.37;
    if(PromptOrFD=="FD")
      shift=0.8;
      
    hCorrYields[iPt]->GetYaxis()->SetNdivisions(505);
    hCorrYields[iPt]->GetXaxis()->SetNdivisions(505);
    cCorrYieldsDisp[iPt] = new TCanvas(Form("cCorrYieldsDispPt%d",iPt),"",10,10,800,800);  
    hCorrYields[iPt]->Draw();
    latex->DrawLatex(mean-shift*mean,hCorrYields[iPt]->GetMaximum()*0.83,Form("mean = %0.f",mean));
    latex->DrawLatex(mean-shift*mean,hCorrYields[iPt]->GetMaximum()*0.75,Form("RMS = %0.1f (%0.1f%)",RMS,relRMS));
    latex->DrawLatex(mean-shift*mean,hCorrYields[iPt]->GetMaximum()*0.67,Form("(x_{max}-x_{min})/#sqrt{12} = %0.1f (%0.1f%)",extimator,relext));
    line[iPt]->Draw("same");
    l->Draw("same");

    cCorrYields[iPt]->SaveAs(Form("CorrYieldsVar%s%s_Pt%d.eps",PromptOrFD.Data(),Method.Data(),iPt));
    cCorrYieldsVsTrial[iPt]->SaveAs(Form("CorrYieldsVsTrial%s%s_Pt%d.eps",PromptOrFD.Data(),Method.Data(),iPt));
    cCorrYieldsDisp[iPt]->SaveAs(Form("CorrYieldsDisp%s%s_Pt%d.eps",PromptOrFD.Data(),Method.Data(),iPt));

    Double_t meanraw[nSets];
    Double_t RMSraw[nSets];
    Double_t extimatorraw[nSets];
    Double_t relRMSraw[nSets];
    Double_t relextraw[nSets];

    cRawYields[iPt] = new TCanvas(Form("cRawYieldsPt%d",iPt),"",10,10,1200,400);  
    cRawYields[iPt]->Divide(3,1);
    
    for(Int_t iSet=0; iSet<nSets; iSet++) {

      hRawYields[iSet][iPt]=(TH1F*)infile.Get(Form("hRawYields%d_Pt%d",iSet+1,iPt)); 
      hRawYields[iSet][iPt]->SetDirectory(0);
      hRawYields[iSet][iPt]->SetStats(kFALSE);
      hRawYields[iSet][iPt]->GetYaxis()->SetTitle("Entries");
      hRawYieldsVsTrial[iSet][iPt]=(TH1F*)infile.Get(Form("hRawYieldsVsTrial%d_Pt%d",iSet+1,iPt));
      if(iSet==nSets-1)
        hRawYieldsVsTrial[iSet][iPt]->GetYaxis()->SetRangeUser(hRawYieldRef[iSet]->GetBinContent(iPt+1)*0.,hRawYieldRef[iSet]->GetBinContent(iPt+1)*2.5);
      hRawYieldsVsTrial[iSet][iPt]->SetDirectory(0);
      hRawYieldsVsTrial[iSet][iPt]->SetStats(kFALSE);
      hRawYieldsVsTrial[iSet][iPt]->SetTitle(Form("Set %d",iSet+1));
      hRawYieldsVsTrial[iSet][iPt]->GetXaxis()->SetTitleSize(0.04);
      hRawYieldsVsTrial[iSet][iPt]->GetYaxis()->SetTitleSize(0.04);
      hRawYieldsVsTrial[iSet][iPt]->GetYaxis()->SetTitleOffset(1.6);
      hRawYieldsVsTrial[iSet][iPt]->GetYaxis()->SetNdivisions(10,10,0);
      hRawYieldsVsTrial[iSet][iPt]->GetXaxis()->SetNdivisions(20,10,0);

      meanraw[iSet] = hRawYields[iSet][iPt]->GetMean();
      RMSraw[iSet] = hRawYields[iSet][iPt]->GetRMS();
      extimatorraw[iSet] = (GetMax(hRawYields[iSet][iPt])-GetMin(hRawYields[iSet][iPt]))/TMath::Sqrt(12);
      relRMSraw[iSet] = RMSraw[iSet]/meanraw[iSet]*100;
      relextraw[iSet] = extimatorraw[iSet]/meanraw[iSet]*100;
    
      linerawvstrial[iSet][iPt] = new TLine(0.5,hRawYieldRef[iSet]->GetBinContent(iPt+1),hRawYieldsVsTrial[iSet][iPt]->GetNbinsX()-0.5,hRawYieldRef[iSet]->GetBinContent(iPt+1));
      linerawvstrial[iSet][iPt]->SetLineColor(kRed);
      linerawvstrial[iSet][iPt]->SetLineWidth(2);
      
      cRawYields[iPt]->cd(iSet+1);
      hRawYieldsVsTrial[iSet][iPt]->Draw();
      latex->DrawLatex(hRawYieldsVsTrial[iSet][iPt]->GetNbinsX()*0.3,hRawYieldsVsTrial[iSet][iPt]->GetMaximum()*0.95,Form("mean = %0.f",meanraw[iSet]));
      latex->DrawLatex(hRawYieldsVsTrial[iSet][iPt]->GetNbinsX()*0.3,hRawYieldsVsTrial[iSet][iPt]->GetMaximum()*0.87,Form("RMS = %0.1f (%0.1f%)",RMSraw[iSet],relRMSraw[iSet]));
      latex->DrawLatex(hRawYieldsVsTrial[iSet][iPt]->GetNbinsX()*0.3,hRawYieldsVsTrial[iSet][iPt]->GetMaximum()*0.79,Form("(x_{max}-x_{min})/#sqrt{12} = %0.1f (%0.1f%)",extimatorraw[iSet],relextraw[iSet]));
      linerawvstrial[iSet][iPt]->Draw("same");
    }
    
    cRawYields[iPt]->SaveAs(Form("RawYieldsVar_Pt%d.eps",iPt));
  }
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    delete hCorrYields[iPt];
    delete hCorrYieldsVsTrial[iPt];
    delete cCorrYieldsVsTrial[iPt];
    delete cCorrYieldsDisp[iPt];
    delete cCorrYields[iPt];
    delete hRawYields[0][iPt];
    delete hRawYieldsVsTrial[0][iPt];
    delete hRawYields[1][iPt];
    delete hRawYieldsVsTrial[1][iPt];
    delete hRawYields[2][iPt];
    delete hRawYieldsVsTrial[2][iPt];
    delete cRawYields[iPt];
    delete line[iPt];

    for(Int_t iSet=0; iSet<nSets; iSet++) {
      delete linerawvstrial[iSet][iPt];
    }
  }
  delete[] cCorrYieldsVsTrial;
  delete[] cCorrYieldsDisp;
  delete[] cCorrYields;
  delete[] hCorrYields;
  delete[] hCorrYieldsVsTrial;
  delete[] hRawYields;
  delete[] hRawYieldsVsTrial;
  delete[] cRawYields;
  delete[] line;
  delete[] linerawvstrial;
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
