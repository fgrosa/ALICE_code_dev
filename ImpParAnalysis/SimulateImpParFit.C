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
#include <TTree.h>
#include <TNtuple.h>

#include "AliDplusCharmFractionIPfitter.h"

#endif

//_____________________________________________________________________________________________
//GLOBAL VARIABLES

//PtBins of the analysis
const Int_t nPtBins = 7;
const Int_t nPtLims = nPtBins+1;
const Double_t PtLims[nPtLims] = {2,3,4,5,6,8,11,16};
//IP range limits
Double_t d0limFD[nPtBins] = {200,300,300,300,300,300,400};
Double_t d0lim[nPtBins] = {1000,1000,1000,1000,1000,1000,1000};
Double_t d0limPrompt[nPtBins] = {300,300,300,300,300,300,300};

//PtBins for the FONLL fraction
const Int_t FONLLptbins[nPtBins] = {2,3,4,5,6,8,9};
TString FONLLfilename="HFPtSpectrum_combinedFD.root";


//input file names
const TString infileMCname="$HOME/ALICE_WORK/Files/Trains/Run1/LHC13/AnalysisResultspPbMC.root";
const TString dirMCname="PWG3_D2H_InvMassDplus";
const TString listMCname="coutputDplus_ImpParpPbMC0100";
const TString infileDataname="$HOME/ALICE_WORK/Files/Trains/Run2/LHC15/LHC15o/AnalysisResults_3050_central_imppar.root";//"$HOME/ALICE_WORK/Files/Trains/Run1/LHC13/AnalysisResultspPbData.root";
const TString dirDataname=dirMCname;
const TString listDataname="coutputDplus_3050_CentralCuts_kINT73050";//"coutputDplus_ImpParpPbData0100";

///generation parameters
Double_t genfraction=0.8;
Double_t S=10000;
Double_t B=10000;

enum {kUnbinned,kBinned,kBinnedBkgSub,kBinnedVarBin};

//_____________________________________________________________________________________________
//FUNCTION PROTOTYPES
Int_t SimulateImpParFit(Int_t Nfits=50,
                        Int_t method=kUnbinned,
                        Double_t nSigmas=2,
                        Bool_t isbkg=kTRUE,
                        Bool_t isSigmaFixed=kFALSE,
                        Int_t genversion=AliDplusCharmFractionIPfitter::kFromHisto,
                        Bool_t realfractions=kFALSE,
                        Bool_t realstatistics=kTRUE,
                        Bool_t applyd0cut=kTRUE,
                        Double_t d0cut=60);

Int_t LoadMCSparses(THnSparseF *&promptsparse, THnSparseF *&trueFDsparse, THnSparseF *&recoFDsparse, THnSparseF *&bkgsparse);
Int_t LoadDataSparsesAndTree(Int_t method, THnSparseF *&datasparse, TNtuple *&tree);
void SetStyle();

//_____________________________________________________________________________________________
//SIMULATION FUNCTION
Int_t SimulateImpParFit(Int_t Nfits,Int_t method,Double_t nSigmas,Bool_t isbkg,Bool_t isSigmaFixed,Int_t genversion,Bool_t realfractions,Bool_t realstatistics,Bool_t applyd0cut,Double_t d0cut) {
  
  //input files
  THnSparseF* hMassPtImpParPrompt=0x0;
  THnSparseF* hMassPtImpParRecoFD=0x0;
  THnSparseF* hMassPtImpParTrueFD=0x0;
  THnSparseF* hMassPtImpParBkg=0x0;
  Int_t loadMC=LoadMCSparses(hMassPtImpParPrompt,hMassPtImpParTrueFD,hMassPtImpParRecoFD,hMassPtImpParBkg);
  if(loadMC>0) {return 1;}
  
  THnSparseF* hMassPtImpParAll=0x0;
  TNtuple* dataTree=0x0;
  Int_t loadData=LoadDataSparsesAndTree(method,hMassPtImpParAll,dataTree);
  if(loadData>0) {return 2;}
  
  Double_t FONLLfractions[nPtBins];
  TGraphAsymmErrors* gfPromptFONLL = 0x0;
  if(realfractions) {
    TFile* FONLLfile=TFile::Open(FONLLfilename.Data(),"READ");
    if(FONLLfile) {
      gfPromptFONLL = (TGraphAsymmErrors*)FONLLfile->Get("gFcCorrConservative");
      FONLLfile->Close();
      Double_t pt;
      for(Int_t iBin=0; iBin<nPtBins; iBin++) {
        gfPromptFONLL->GetPoint(FONLLptbins[iBin],pt,FONLLfractions[iBin]);
      }
    }
    else {
      cout << "Warning: unable to load FONLL fraction! Set gen frac = "<< genfraction << endl;
    }
  }
  
  //analysis
  AliDplusCharmFractionIPfitter *ImpParFitter = new AliDplusCharmFractionIPfitter();
  ImpParFitter->SetDataSparse(hMassPtImpParAll);
  ImpParFitter->SetMCPromptSparse(hMassPtImpParPrompt);
  ImpParFitter->SetMCTrueFDSparse(hMassPtImpParTrueFD);
  ImpParFitter->SetMCRecoFDSparse(hMassPtImpParRecoFD);

  ImpParFitter->SetFDFunction(AliDplusCharmFractionIPfitter::kConvolution);    
  ImpParFitter->SetBkgFunction(AliDplusCharmFractionIPfitter::kDoubleGaussExpoSymm);   
  ImpParFitter->SetNSigmas(nSigmas);
  ImpParFitter->SetGenPromptFraction(genfraction);
  ImpParFitter->SetPID(kTRUE);
  
  Double_t initparprompt[5] = {0.9,0.,35,300,1};
  Double_t initparFD[5] = {0.5,0.,100,10,1};
  Double_t initparBkg[14] = {0.3,-30.,50.,50.,100.,0.8,0.4,30.,50.,50.,100.,0.8,0.5,1.};
  
  SetStyle();
  
  TH1F* hBiasVsPt = new TH1F("hBiasVsPt","",nPtBins,PtLims);
  hBiasVsPt->SetLineColor(kRed);
  hBiasVsPt->SetLineWidth(2);
  hBiasVsPt->SetMarkerStyle(20);
  hBiasVsPt->SetMarkerSize(1.5);
  hBiasVsPt->SetMarkerColor(kRed);
  hBiasVsPt->GetYaxis()->SetTitle("< #it{f}_{prompt}^{meas} - #it{f}_{prompt}^{true} >");
  hBiasVsPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  TH1F* hResVsPt = new TH1F("hResVsPt","",nPtBins,PtLims);
  hResVsPt->SetLineColor(kRed);
  hResVsPt->SetLineWidth(2);
  hResVsPt->SetMarkerStyle(20);
  hResVsPt->SetMarkerSize(1.5);
  hResVsPt->SetMarkerColor(kRed);
  hResVsPt->GetYaxis()->SetTitle("#sigma(#it{f}_{prompt}^{meas} - #it{f}_{prompt}^{true})");
  hResVsPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  TH1F* hPullsVsPt = new TH1F("hPullsVsPt","",nPtBins,PtLims);
  hPullsVsPt->SetLineColor(kRed);
  hPullsVsPt->SetLineWidth(2);
  hPullsVsPt->SetMarkerStyle(20);
  hPullsVsPt->SetMarkerSize(1.5);
  hPullsVsPt->SetMarkerColor(kRed);
  hPullsVsPt->GetYaxis()->SetTitle("#sigma((#it{f}_{prompt}^{meas} - #it{f}_{prompt}^{true})/#sigma_{#it{f}_{prompt}^{meas}})");
  hPullsVsPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  TH1F* hSigmaTrueVsPt = new TH1F("hSigmaTrueVsPt","",nPtBins,PtLims);
  hSigmaTrueVsPt->SetLineColor(kBlue);
  hSigmaTrueVsPt->SetLineWidth(2);
  hSigmaTrueVsPt->SetMarkerStyle(21);
  hSigmaTrueVsPt->SetMarkerSize(1.5);
  hSigmaTrueVsPt->SetMarkerColor(kBlue);
  hSigmaTrueVsPt->GetYaxis()->SetTitle("#sigma_{prompt} (#mum)");
  hSigmaTrueVsPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  TH1F* hSigmaVsPt = new TH1F("hSigmaVsPt","",nPtBins,PtLims);
  hSigmaVsPt->SetLineColor(kRed);
  hSigmaVsPt->SetLineWidth(2);
  hSigmaVsPt->SetMarkerStyle(20);
  hSigmaVsPt->SetMarkerSize(1.5);
  hSigmaVsPt->SetMarkerColor(kRed);
  hSigmaVsPt->GetYaxis()->SetTitle("<#sigma_{prompt}^{meas}> (#mum)");
  hSigmaVsPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  TH1F* hSigmaResVsPt = new TH1F("hSigmaResVsPt","",nPtBins,PtLims);
  hSigmaResVsPt->SetLineColor(kRed);
  hSigmaResVsPt->SetLineWidth(2);
  hSigmaResVsPt->SetMarkerStyle(20);
  hSigmaResVsPt->SetMarkerSize(1.5);
  hSigmaResVsPt->SetMarkerColor(kRed);
  hSigmaResVsPt->GetYaxis()->SetTitle("#sigma(#sigma_{prompt}^{meas} - #sigma_{prompt}^{true}) (#mum)");
  hSigmaResVsPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  TH1F* hSigmaPullsVsPt = new TH1F("hSigmaPullsVsPt","",nPtBins,PtLims);
  hSigmaPullsVsPt->SetLineColor(kRed);
  hSigmaPullsVsPt->SetLineWidth(2);
  hSigmaPullsVsPt->SetMarkerStyle(20);
  hSigmaPullsVsPt->SetMarkerSize(1.5);
  hSigmaPullsVsPt->SetMarkerColor(kRed);
  hSigmaPullsVsPt->GetYaxis()->SetTitle("#sigma(#sigma_{prompt}^{meas} - #sigma_{prompt}^{true})/#sigma_{#sigma_{prompt}^{meas}})");
  hSigmaPullsVsPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  
  TLine* line = new TLine(PtLims[0],0,PtLims[nPtBins],0);
  line->SetLineColor(kBlue);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  
  TH1F** hRes = new TH1F*[nPtBins];
  TH1F** hPulls = new TH1F*[nPtBins];
  TH1F** hSigma = new TH1F*[nPtBins];
  TH1F** hSigmaRes = new TH1F*[nPtBins];
  TH1F** hSigmaPulls = new TH1F*[nPtBins];
  Bool_t print = kTRUE;

  for(Int_t iBin=0; iBin<nPtBins; iBin++) {
    hRes[iBin] = new TH1F(Form("hRes%d",iBin),"",100,-0.5,0.5);
    hRes[iBin]->SetStats(kTRUE);
    hRes[iBin]->GetYaxis()->SetTitle("Entries");
    hRes[iBin]->GetXaxis()->SetTitle("#it{f}_{prompt}^{meas}-#it{f}_{prompt}^{true}");
    hRes[iBin]->SetFillColor(kBlue);
    hRes[iBin]->SetFillStyle(3004);
    hPulls[iBin] = new TH1F(Form("hPulls%d",iBin),"",100,-10,10);
    hPulls[iBin]->GetXaxis()->SetTitle("(#it{f}_{prompt}^{meas} - #it{f}_{prompt}^{true})/#sigma_{#it{f}_{prompt}^{meas}}");
    hPulls[iBin]->GetYaxis()->SetTitle("Entries");
    hPulls[iBin]->SetFillColor(kBlue);
    hPulls[iBin]->SetFillStyle(3004);
    hSigma[iBin] = new TH1F(Form("hSigma%d",iBin),"",100,10,60);
    hSigma[iBin]->GetXaxis()->SetTitle("#sigma_{prompt}^{meas} (#mum)");
    hSigma[iBin]->GetYaxis()->SetTitle("Entries");
    hSigma[iBin]->SetFillColor(kBlue);
    hSigma[iBin]->SetFillStyle(3004);
    hSigmaRes[iBin] = new TH1F(Form("hSigmaRes%d",iBin),"",100,-25,25);
    hSigmaRes[iBin]->GetXaxis()->SetTitle("#sigma_{prompt}^{meas} - #sigma_{prompt}^{true} (#mum)");
    hSigmaRes[iBin]->GetYaxis()->SetTitle("Entries");
    hSigmaRes[iBin]->SetFillColor(kBlue);
    hSigmaRes[iBin]->SetFillStyle(3004);
    hSigmaPulls[iBin] = new TH1F(Form("hSigmaPulls%d",iBin),"",100,-10,10);
    hSigmaPulls[iBin]->GetXaxis()->SetTitle("(#sigma_{prompt}^{meas} - #sigma_{prompt}^{true})/#sigma_{#sigma_{prompt}^{meas}}");
    hSigmaPulls[iBin]->GetYaxis()->SetTitle("Entries");
    hSigmaPulls[iBin]->SetFillColor(kBlue);
    hSigmaPulls[iBin]->SetFillStyle(3004);
    
    if(realfractions && gfPromptFONLL) {ImpParFitter->SetGenPromptFraction(FONLLfractions[iBin]);}
    
    ImpParFitter->SetRebinImpParHistos(2);
    ImpParFitter->SetPtLims(PtLims[iBin],PtLims[iBin+1]);
    ImpParFitter->SetNSigmas(nSigmas);
    ImpParFitter->SetNSigmaSBLimits(4,15);
    if(realstatistics)
      ImpParFitter->GetSignal(4,0,0,1.68,2.05,AliDplusCharmFractionIPfitter::kCentralValue);
    else {
      ImpParFitter->SetSignal(S);
      ImpParFitter->SetBkg(B);
    }
    ImpParFitter->SetFitOptions("RLEM0");
    ImpParFitter->SetRelativeLimitSigmaPrompt(0.2);
    ImpParFitter->SetLimitsPromptFrac(0.2,1.5);
    ImpParFitter->SetInitialParameters(initparprompt,initparFD,initparBkg);
    ImpParFitter->PrefitStep(-d0limPrompt[iBin],d0limPrompt[iBin],-d0limFD[iBin],d0limFD[iBin],-1000,1000);

    Double_t sigmatrue = ImpParFitter->GetPromptSigmaMC();
    Double_t sigmatrueerr = ImpParFitter->GetPromptSigmaMCErr();
    hSigmaTrueVsPt->SetBinContent(iBin+1,sigmatrue);
    hSigmaTrueVsPt->SetBinError(iBin+1,sigmatrueerr);
    
    ImpParFitter->FixSigmaPromptFromMC(isSigmaFixed);
  
    for(Int_t iGen=0; iGen<Nfits; iGen++) {
      if(realfractions) 
        cout << "\n\n ***** Generation number " << iGen << "   Bin number " << iBin << "   f = " << FONLLfractions[iBin] << " *****\n\n"<<endl;
      else
        cout << "\n\n ***** Generation number " << iGen << "   Bin number " << iBin << "   f = " << genfraction << " *****\n\n"<<endl;
      
      if(iGen>0)
        print = kFALSE;
      if(method==kBinned || method==kBinnedBkgSub || method==kBinnedVarBin) {
        ImpParFitter->GenerateEntries(isbkg,genversion,AliDplusCharmFractionIPfitter::kHisto);
        if(method==kBinnedBkgSub) {ImpParFitter->SetBkgSubtraction(kTRUE);}
        if(method==kBinnedVarBin) {ImpParFitter->SetVariableBinningHisto(kTRUE,5);}
        ImpParFitter->FitHisto(-d0lim[iBin],d0lim[iBin],print);
      }
      else if(method==kUnbinned) {
        ImpParFitter->GenerateEntries(isbkg,genversion,AliDplusCharmFractionIPfitter::kTree);        
        ImpParFitter->FitTree(-d0lim[iBin],d0lim[iBin],print);
      }
      else {
        cout << "Error: only binned or unbinned fit are supported!" << endl;
        return 3;
      }
   
      Double_t fraction;
      Double_t fractionerr;
      if(!applyd0cut) {
        fraction = ImpParFitter->GetPromptFraction();
        fractionerr = ImpParFitter->GetPromptFractionErr();
      }
      else {
        ImpParFitter->GetPromptFractionWithIPCut(d0cut, fraction, fractionerr);
      }
      Double_t sigmaprompt = ImpParFitter->GetPromptSigma();
      Double_t sigmaprompterr = ImpParFitter->GetPromptSigmaErr();
      
      Double_t truefraction;
      Double_t trueerr;//dummy variable -> always 0
      if(!applyd0cut) {
        if(realfractions) {
          truefraction = FONLLfractions[iBin];
        }
        else {
          truefraction = genfraction;
        }
      }
      else {
        ImpParFitter->GetPromptFractionWithIPCut(d0cut, truefraction, trueerr, kTRUE);
      }
      
      hRes[iBin]->Fill(fraction-truefraction);
      hPulls[iBin]->Fill((fraction-truefraction)/fractionerr);
      hSigma[iBin]->Fill(sigmaprompt);
      hSigmaRes[iBin]->Fill(sigmaprompt-sigmatrue);
      if(sigmaprompterr!=0)
        hSigmaPulls[iBin]->Fill((sigmaprompt-sigmatrue)/sigmaprompterr);
      else
        hSigmaPulls[iBin]->Fill(0);        
    }
    
    print=kTRUE;
    Double_t resmean = hRes[iBin]->GetMean();
    Double_t resmeanerr = hRes[iBin]->GetMeanError();
    Double_t RMS = hRes[iBin]->GetRMS();
    Double_t RMSerr = hRes[iBin]->GetRMSError();
    Double_t pull = hPulls[iBin]->GetRMS();
    Double_t pullerr = hPulls[iBin]->GetRMSError();
    Double_t meansigma = hSigma[iBin]->GetMean();
    Double_t meansigmaerr = hSigma[iBin]->GetMeanError();
    Double_t resmeansigma = hSigmaRes[iBin]->GetMean();
    Double_t resmeansigmaerr = hSigmaRes[iBin]->GetMeanError();
    Double_t RMSsigma = hSigmaRes[iBin]->GetRMS();
    Double_t RMSsigmaerr = hSigmaRes[iBin]->GetRMSError();
    Double_t pullsigma = hSigmaPulls[iBin]->GetRMS();
    Double_t pullsigmaerr = hSigmaPulls[iBin]->GetRMSError();
    hBiasVsPt->SetBinContent(iBin+1,resmean);
    hBiasVsPt->SetBinError(iBin+1,resmeanerr);
    hResVsPt->SetBinContent(iBin+1,RMS);
    hResVsPt->SetBinError(iBin+1,RMSerr);
    hPullsVsPt->SetBinContent(iBin+1,pull);
    hPullsVsPt->SetBinError(iBin+1,pullerr);
    hSigmaVsPt->SetBinContent(iBin+1,meansigma);
    hSigmaVsPt->SetBinError(iBin+1,meansigmaerr);
    hSigmaVsPt->SetBinContent(iBin+1,meansigma);
    hSigmaResVsPt->SetBinContent(iBin+1,RMSsigma);
    hSigmaResVsPt->SetBinError(iBin+1,RMSsigmaerr);
    hSigmaPullsVsPt->SetBinContent(iBin+1,pullsigma);
    hSigmaPullsVsPt->SetBinError(iBin+1,pullsigmaerr);
  }

  TLegend* sigmalegend = new TLegend(0.6,0.6,0.8,0.8);
  sigmalegend->SetTextSize(0.05);
  sigmalegend->SetFillStyle(0);
  sigmalegend->AddEntry(hSigmaTrueVsPt," #sigma_{prompt}^{MC prefit} ","lpe");
  sigmalegend->AddEntry(hSigmaVsPt,"< #sigma_{prompt}^{meas} >","lpe");  
  
  TCanvas* cBiasVsPt = new TCanvas("cBiasVsPt","cBiasVsPt",800,800);
  hBiasVsPt->Draw("E1");
  line->Draw("same");
  TCanvas* cResVsPt = new TCanvas("cResVsPt","cResVsPt",800,800);
  hResVsPt->Draw("E1");
  TCanvas* cPullsVsPt = new TCanvas("cPullsVsPt","cPullsVsPt",800,800);
  hPullsVsPt->Draw("E1"); 
  TCanvas *cSigmaVsPt = new TCanvas("cSigmaVsPt","cSigmaVsPt",800,800);
  hSigmaTrueVsPt->Draw();
  hSigmaTrueVsPt->GetYaxis()->SetRangeUser(hSigmaTrueVsPt->GetMinimum()*0.5,hSigmaTrueVsPt->GetMaximum()*1.5);
  hSigmaVsPt->Draw("same");
  sigmalegend->Draw("same");
  TCanvas *cSigmaResVsPt = new TCanvas("cSigmaResVsPt","cSigmaResVsPt",800,800);
  hSigmaResVsPt->Draw();
  TCanvas* cSigmaPullsVsPt = new TCanvas("cSigmaPullsVsPt","cSigmaPullsVsPt",800,800);
  hSigmaPullsVsPt->Draw("E1");

  //output files
  TString bkg = "bkg";
  if(!isbkg)
    bkg="nobkg";
  TString sigmapar = "sigmafree";
  if(isSigmaFixed)
    bkg="sigmafixed";
  
  TString methodname="unbinned";
  if(method==kBinned) {methodname="binned";}
  if(method==kBinnedBkgSub) {methodname="binnedbkgsub";}
  if(method==kBinnedVarBin) {methodname="binnedvarbin";}
  
  TString outfracname = Form("PromptFraction_TOYMC_Nbfrac_%s_%s_%s.root",bkg.Data(),methodname.Data(),sigmapar.Data());
  TString outsigmaname = Form("SigmaPrompt_TOYMC_Nbfrac_%s_%s_%s.root",bkg.Data(),methodname.Data(),sigmapar.Data());
  if(!realfractions) {
    outfracname.ReplaceAll("Nbfrac",Form("%0.2f",genfraction));
    outsigmaname.ReplaceAll("Nbfrac",Form("%0.2f",genfraction));
  }
  TFile outfracfile(outfracname,"RECREATE");
  hBiasVsPt->Write();
  hResVsPt->Write();
  hPullsVsPt->Write();
  for(Int_t iBin=0; iBin<nPtBins; iBin++) {
    hRes[iBin]->Write();
    hPulls[iBin]->Write();
  }
  outfracfile.Close();

  TFile outsigmafile(outsigmaname,"RECREATE");
  hSigmaVsPt->Write();
  hSigmaTrueVsPt->Write();
  hSigmaResVsPt->Write();
  hSigmaPullsVsPt->Write();
  for(Int_t iBin=0; iBin<nPtBins; iBin++) {
    hSigma[iBin]->Write();
    hSigmaRes[iBin]->Write();
    hSigmaPulls[iBin]->Write();
  }
  outsigmafile.Close();
  
  return 0;
}

//_____________________________________________________________________________________________
//LOAD MC SPARSES FUNCTION
Int_t LoadMCSparses(THnSparseF *&promptsparse, THnSparseF *&trueFDsparse, THnSparseF *&recoFDsparse, THnSparseF *&bkgsparse) {
  
  cout<<"Opening MC file " <<infileMCname<< "..." << endl;
  TFile* infileMC = TFile::Open(infileMCname.Data(),"READ");
  TDirectoryFile* dirMC=0x0;
  TList* listMC=0x0;
  
  if(infileMC) {dirMC=(TDirectoryFile*)infileMC->Get(dirMCname.Data()); cout << "MC file opened!" << endl;}
  else {cerr << "Error: File " << infileMCname << " not found. Exit." << endl; return 1;}
  if(dirMC) listMC=(TList*)dirMC->Get(listMCname.Data());
  else {cerr << "Error: Wrong TDirectoryFile name " << dirMCname << ". Exit." << endl; return 2;}
  if(listMC) {
    promptsparse=(THnSparseF*)listMC->FindObject("hMassPtImpParPrompt");
    recoFDsparse=(THnSparseF*)listMC->FindObject("hMassPtImpParBfeed");
    trueFDsparse=(THnSparseF*)listMC->FindObject("hMassPtImpParTrueBfeed");
    bkgsparse=(THnSparseF*)listMC->FindObject("hMassPtImpParBkg");
    cout << "Sparses got!" << endl;
  }
  else {cerr << "Error: Wrong TList name " << listMCname << ". Exit." << endl; return 3;}
  if(!promptsparse) {cerr << "Error: No MC sparse for prompt D+. Check if the name hMassPtImpParPrompt is right!" << endl; return 4;}
  if(!recoFDsparse) {cerr << "Error: No MC sparse for feed-down D+. Check if the name hMassPtImpParBfeed is right!" << endl; return 5;}
  if(!trueFDsparse) {cerr << "Error: No MC sparse for true feed-down D+. Check if the name hMassPtImpParTrueBfeed is right!" << endl; return 6;}
  if(!bkgsparse) {cout << "Warning: No MC sparse for the background. Check if the name hMassPtImpParBkg is right!" << endl;}
  
  infileMC->Close();
  cout<<"MC file closed."<< endl;
  
  return 0;
}

//_____________________________________________________________________________________________
//LOAD DATA SPARSE AND TREE FUNCTION
Int_t LoadDataSparsesAndTree(Int_t method, THnSparseF *&datasparse, TNtuple *&tree) {
  
  cout<<"Opening data file " <<infileMCname<< "..." << endl;
  TFile* infileData = TFile::Open(infileDataname.Data(),"READ");
  TDirectoryFile* dirData=0x0;
  TList* listData=0x0;
  THnSparseF* hMassPtImpParAll=0x0;
  TNtuple* dataTree=0x0;
  
  if(infileData) {dirData=(TDirectoryFile*)infileData->Get(dirDataname.Data()); cout << "Data file opened!" << endl;}
  else {cerr << "Error: File " << infileDataname << " not found. Exit." << endl; return 1;}
  if(dirData) listData=(TList*)dirData->Get(listDataname.Data());
  else {cerr << "Error: Wrong TDirectoryFile name " << dirDataname << ". Exit." << endl; return 2;}
  if(listData) {
    datasparse=(THnSparseF*)listData->FindObject("hMassPtImpParAll");
    cout << "Sparse got!" << endl;
    if(method==kUnbinned) {tree=(TNtuple*)dirData->Get("fNtupleDplus"); cout << "TNtuple got!" << endl;}
  }
  else {cerr << "Error: Wrong TList name " << listDataname << ". Exit." << endl; return 3;}
  if(!datasparse && method==kUnbinned) {cout << "Warning: No data sparse. Check if the name hMassPtImpParAll is right!" << endl;}
  if(!datasparse && method==kBinned) {cout << "Error: No data sparse. Check if the name hMassPtImpParAll is right!" << endl; return 4;}
  if(!tree && method==kUnbinned) {cout << "Error: No data tree. Check if the name fNtupleDplus is right!" << endl; return 5;}
  
  infileData->Close();
  cout<<"Data file closed."<< endl;
  
  return 0;
}

//_____________________________________________________________________________________________
//DRAW STYLE FUNCTION
void SetStyle() {
  
  cout << "Setting drawing style!" << endl;
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.05,"xyzt");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.16);
}

