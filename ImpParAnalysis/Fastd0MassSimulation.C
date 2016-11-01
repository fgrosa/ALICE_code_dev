#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TF2.h>
#include <THnSparse.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TBranch.h>
#include <TColor.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include <TLine.h>
#include <TDatabasePDG.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>

//#include "AliDplusCharmFractionIPfitter.h"

#endif

//______________________________________________________________________________________________________________________
//global variables
const Double_t massD = 1.869;
const Double_t masssignalpars[2] = {massD,0.010};
const Double_t massbkgpars[1] = {-2};
const Double_t impparpromptpars[4] = {0.,40.,0.92,43};
const Double_t impparFDpars[4] = {0.88,0.,70,8.};
const Double_t impparbkgpars[5] = {-60,0.66,40,115,60.};

enum {kBinned,kUnbinned};
//______________________________________________________________________________________________________________________

Int_t Fastd0MassSimulation(Int_t fitmeth=kBinned, Bool_t generatetree=kTRUE, Int_t nSgn=50000, Int_t nBkg=500000, Double_t promptfraction=0.95);
void GenerateTree(Int_t nSgn,Int_t nBkg,Double_t promptfraction,TF1* fMassSignal,TF1* fMassBkg,TF1* fImpParPrompt,TF1* fImpParRecoFD,TF1* fImpParBkg,TString outfilename="Massd0tree_FastSim.root", TString treename="Massd0tree");
TTree* OpenTreeFile(TString infilename="Massd0tree_FastSim.root", TString treename="Massd0tree");
void CloseFile(TString infilename="Massd0tree_FastSim.root");
Double_t Gauss(Double_t var, Double_t mean, Double_t sigma);
Double_t ExpoDouble(Double_t d0, Double_t mean, Double_t lambda);
Double_t GetImpParPromptProb(Double_t* d0, Double_t* pars);
Double_t GetImpParTrueFDProb(Double_t* d0, Double_t* pars);
Double_t GetImpParRecoFDProb(Double_t* d0, Double_t* pars);
Double_t GetImpParBkgProb(Double_t* d0, Double_t* pars);
Double_t ResolutionConvolution(Double_t d0, Double_t *parprompt, Double_t *parFD, Double_t d0min, Double_t d0max, Int_t nsteps);
Double_t GetImpParProb(Double_t* d0, Double_t* pars);
Double_t GetMassPeakProb(Double_t* mass, Double_t* pars);
Double_t GetMassBkgProb(Double_t* mass, Double_t* pars);
Double_t GetMassProb(Double_t* m, Double_t* pars);
Double_t GetMassImpParTotProb(Double_t* vars, Double_t* pars);

Int_t Fastd0MassSimulation(Int_t fitmeth, Bool_t generatetree, Int_t nSgn, Int_t nBkg, Double_t promptfraction)
{
  Double_t signalfrac = Double_t(nSgn)/(nSgn+nBkg); ///BEWARE: IN ALL THE MASS RANGE, NOT ONLY THE SIGNAL REGION!!!!
  
  TF1* fMassSignal = new TF1("fMassSignal",GetMassPeakProb,1.66,2.06,1);
  fMassSignal->SetParameter(0,1.);
  TF1* fMassBkg = new TF1("fMassBkg",GetMassBkgProb,1.66,2.06,1);
  fMassBkg->SetParameter(0,1.);
  fMassBkg->SetLineColor(kBlack);
  TF1* fMass = new TF1("fMass",GetMassProb,1.66,2.06,2);
  fMass->SetParameters(1.,signalfrac);
  fMass->SetLineColor(kBlue);
  
  TCanvas* cMassFunc = new TCanvas("cMassFunc","cMassFunc",1920,1080);
  cMassFunc->Divide(2,1);
  cMassFunc->cd(1);
  fMassSignal->Draw();
  fMassBkg->Draw("same");
  cMassFunc->cd(2);
  fMass->Draw();

  TF1* fImpParPrompt = new TF1("fImpParPrompt",GetImpParPromptProb,-1000,1000,1);
  fImpParPrompt->SetParameter(0.,1.);
  fImpParPrompt->SetLineColor(kGreen+3);
  TF1* fImpParTrueFD = new TF1("fImpParTrueFD",GetImpParTrueFDProb,-1000,1000,1);
  fImpParTrueFD->SetParameter(0,1.);
  fImpParTrueFD->SetLineColor(kBlue+3);
  TF1* fImpParRecoFD = new TF1("fImpParRecoFD",GetImpParRecoFDProb,-1000,1000,1);
  fImpParRecoFD->SetParameter(0,1.);
  fImpParRecoFD->SetLineColor(kBlue+3);
  TF1* fImpParBkg = new TF1("fImpParBkg",GetImpParBkgProb,-1000,1000,1);
  fImpParBkg->SetParameter(0,1.);
  fImpParBkg->SetLineColor(kMagenta+2);
  TF1* fImpPar = new TF1("fImpPar",GetImpParProb,-1000,1000,3);
  fImpPar->SetParameters(1.,signalfrac,promptfraction);
  
  TCanvas* cImpParFunc = new TCanvas("cImpParFunc","cImpParFunc",1920,1080);
  cImpParFunc->Divide(2,1);
  cImpParFunc->cd(1)->SetLogy();
  fImpParPrompt->GetYaxis()->SetRangeUser(1.e-8,3.e-2);
  fImpParPrompt->Draw();
  fImpParRecoFD->Draw("same");
  fImpParBkg->Draw("same");
  cImpParFunc->cd(2)->SetLogy();
  fImpPar->Draw();

  TF2* fMassImpPar = new TF2("fMassImpPar",GetMassImpParTotProb,1.66,2.06,-1000,1000,3);
  fMassImpPar->SetParameters(nSgn+nBkg,nSgn,promptfraction);

  TCanvas* cMassImpParFunc = new TCanvas("cMassImpParFunc","cMassImpParFunc",1920,1080);
  fMassImpPar->SetNpy(200);
  fMassImpPar->Draw("surf4");
  
  if(generatetree)
    GenerateTree(nSgn,nBkg,promptfraction,fMassSignal,fMassBkg,fImpParPrompt,fImpParRecoFD,fImpParBkg,"Massd0tree_FastSim.root");

  TTree* tree=0x0;
  tree=OpenTreeFile();
  if(!tree) {return 1;}
  
  TH2F* hMassImpPar = new TH2F("hMassImpPar","",50,1.66,2.06,200,-1000,1000);
  
  TCanvas* cMassImpPar = new TCanvas("cMassImpPar","",1920,1080);
  tree->Project("hMassImpPar","d0:mass");
  hMassImpPar->Draw("LEGO");
  Double_t binarea = hMassImpPar->GetXaxis()->GetBinWidth(10)*hMassImpPar->GetYaxis()->GetBinWidth(10);
  
  if(fitmeth==kBinned) {fMassImpPar->SetParameters((nSgn+nBkg)*binarea,nSgn*binarea,promptfraction); hMassImpPar->Fit("fMassImpPar","RLEM");}
  else if(fitmeth==kUnbinned) {tree->UnbinnedFit("fMassImpPar","d0:mass","d0<300 && d0>-300 && mass<1.95 && mass>1.75");}
  else {CloseFile(); return 2;}
  
  CloseFile();
  
  return 0;
 }

//______________________________________________________________________________
void GenerateTree(Int_t nSgn,Int_t nBkg,Double_t promptfraction,TF1* fMassSignal,TF1* fMassBkg,TF1* fImpParPrompt,TF1* fImpParRecoFD,TF1* fImpParBkg,TString outfilename,TString treename)
{
  Double_t mass;
  Double_t d0;
  TTree *tree = new TTree(treename.Data(),treename.Data());
  TBranch* massbranch = (TBranch*)tree->Branch("mass",&mass);
  TBranch* d0branch = (TBranch*)tree->Branch("d0",&d0);
  
  cout << "Prompt generation started." << endl;
  for(Int_t iPrompt=0; iPrompt<nSgn*promptfraction; iPrompt++) {
    mass = fMassSignal->GetRandom();
    d0 = fImpParPrompt->GetRandom();
    tree->Fill();
    if((iPrompt+1)%(Int_t)(nSgn*promptfraction/2)==0) {cout << "Generated " << iPrompt+1 << " prompt entries"<<endl;}
  }
  cout << "Prompt generation completed. Skip to FD." << endl;
  
  for(Int_t iFD=0; iFD<nSgn*(1-promptfraction); iFD++) {
    mass = fMassSignal->GetRandom();
    d0 = fImpParRecoFD->GetRandom();
    tree->Fill();
    if((iFD+1)%(Int_t)(nSgn*(1-promptfraction)/2)==0) {cout << "Generated " << iFD+1 << " FD entries"<<endl;}
  }

  cout << "FD generation completed. Skip to bkg." << endl;
  for(Int_t iBkg=0; iBkg<nBkg; iBkg++) {
    mass = fMassBkg->GetRandom();
    d0 = fImpParBkg->GetRandom();
    tree->Fill();
    if((iBkg+1)%(Int_t)(nBkg/2)==0) {cout << "Generated " << iBkg+1 << " bkg entries"<<endl;}
  }

  TFile outfile(outfilename.Data(),"RECREATE");
  tree->Write();
  outfile.Close();
}

//______________________________________________________________________________
TTree* OpenTreeFile(TString infilename, TString treename)
{
  TFile* infile=TFile::Open(infilename.Data(),"READ");
  TTree* tree = (TTree*)infile->Get(treename.Data());
  if(!tree) {infile->Close(); return 0x0;}
  return tree;
}

//______________________________________________________________________________
void CloseFile(TString infilename) {
  TFile* infile=TFile::Open(infilename.Data(),"READ");
  infile->Close();
}

//______________________________________________________________________________
Double_t Gauss(Double_t var,Double_t mean,Double_t sigma)
{
  return TMath::Gaus(var,mean,sigma,kTRUE); //kTRUE -> normalised at 1
}

//______________________________________________________________________________
Double_t ExpoDouble(Double_t d0, Double_t mean, Double_t lambda)
{
  return 1./(2.*lambda)*TMath::Exp(-TMath::Abs(d0-mean)/lambda);
}

//______________________________________________________________________________
Double_t GetImpParPromptProb(Double_t* d0, Double_t* pars)
{  
  Double_t d00 = d0[0];
  Double_t norm = pars[0];
  Double_t meanD = impparpromptpars[0];
  Double_t sigma = impparpromptpars[1];
  Double_t fractionG = impparpromptpars[2];
  Double_t lambdaD = impparpromptpars[3];
    
  return norm*((1.-fractionG)*ExpoDouble(d00,meanD,lambdaD)+fractionG*Gauss(d00,meanD,sigma));
}

//______________________________________________________________________________
Double_t GetImpParTrueFDProb(Double_t* d0, Double_t* pars)
{
  Double_t d00 = d0[0];
  Double_t norm = pars[0];
  Double_t fraction1 = impparFDpars[0];
  Double_t meanB = impparFDpars[1];
  Double_t lambda1 = impparFDpars[2];
  Double_t lambda2 = impparFDpars[3];
    
  return norm*(fraction1*ExpoDouble(d00,meanB,lambda1)+(1.-fraction1)*ExpoDouble(d00,meanB,lambda2));
}

//______________________________________________________________________________
Double_t GetImpParRecoFDProb(Double_t* d0, Double_t* pars)
{
  Double_t d00 = d0[0];
  Double_t norm = pars[0];
  Double_t parprompt[1] = {1.};
  Double_t parFD[1] = {1.};

  return norm*ResolutionConvolution(d00,parprompt,parFD,-1000,1000,500);
}

//______________________________________________________________________________
Double_t GetImpParBkgProb(Double_t* d0, Double_t* pars)
{
  Double_t d00 = d0[0];
  Double_t norm = pars[0];
  Double_t mean1 = impparbkgpars[0];
  Double_t fraction1 = impparbkgpars[1];
  Double_t sigma1 = impparbkgpars[2];
  Double_t lambda1 = impparbkgpars[3];
  Double_t mean2 =  impparbkgpars[4];

  return norm/2*(((1.-fraction1)*ExpoDouble(d00,mean1,lambda1)+fraction1*Gauss(d00,mean1,sigma1))+
                 ((1.-fraction1)*ExpoDouble(d00,mean2,lambda1)+fraction1*Gauss(d00,mean2,sigma1)));
}

//______________________________________________________________________________
Double_t ResolutionConvolution(Double_t d0, Double_t *parprompt, Double_t *parFD, Double_t d0min, Double_t d0max, Int_t nsteps)
{
  Double_t d0true[1];
  Double_t diffd0[1];

  Double_t sum=0.;
  Double_t dd0=(d0max-d0min)/(Double_t)nsteps;

  for(Int_t iStep=0;iStep<nsteps;iStep++){
    d0true[0]=d0min+dd0*iStep;
    diffd0[0]=d0-d0true[0];
    Double_t molt = GetImpParTrueFDProb(d0true,parFD)*GetImpParPromptProb(diffd0,parprompt);
    sum = sum+molt*dd0;
  }
    
  return sum;
}

//______________________________________________________________________________
Double_t GetImpParProb(Double_t* d0, Double_t* pars)
{
  Double_t d00[1] = {d0[0]};
  Double_t norm = pars[0];
  Double_t signal = pars[1];
  Double_t promptfraction = pars[2];
  Double_t parprompt[1] = {1.};
  Double_t parFD[1] = {1.};
  Double_t parbkg[1] = {1.};
  
  return norm*(signal/norm*(promptfraction*GetImpParPromptProb(d00,parprompt)+(1-promptfraction)*GetImpParRecoFDProb(d00,parFD))+(norm-signal)/norm*GetImpParBkgProb(d00,parbkg));

}

//______________________________________________________________________________
Double_t GetMassPeakProb(Double_t* mass, Double_t* pars)
{
  Double_t mmass = mass[0];
  Double_t norm = pars[0];
  Double_t mean = masssignalpars[0];
  Double_t sigma = masssignalpars[1];
  
  return norm*Gauss(mmass,mean,sigma);
}

//______________________________________________________________________________
Double_t GetMassBkgProb(Double_t* mass, Double_t* pars)
{
  Double_t mmass = mass[0];
  Double_t norm = pars[0];
  Double_t coef1 = massbkgpars[0];
    
  return norm*coef1/(TMath::Exp(coef1*2.06)-TMath::Exp(coef1*1.66))*TMath::Exp(coef1*mmass);
}

//______________________________________________________________________________
Double_t GetMassProb(Double_t* m, Double_t* pars)
{
  Double_t mass[1] = {m[0]};
  Double_t norm = pars[0];
  Double_t signal = pars[1];
  Double_t massmean = masssignalpars[0];
  Double_t masssigma = masssignalpars[1];
  Double_t bkgcoef = massbkgpars[0];
    
  Double_t parMpeak[1] = {1.};
  Double_t parMbkg[1]= {1.};

  return norm*((norm-signal)/norm*GetMassBkgProb(mass,parMbkg)+signal/norm*GetMassPeakProb(mass,parMpeak));
}

//______________________________________________________________________________
Double_t GetMassImpParTotProb(Double_t* vars, Double_t* pars)
{
  Double_t mass[1] = {vars[0]};
  Double_t d0[1] = {vars[1]};
  Double_t norm = pars[0];
  Double_t S = pars[1];
  Double_t promptfraction = pars[2];
  
  Double_t parprompt[1] = {1.};
  Double_t parFD[1] = {1.};
  Double_t parbkg[1] = {1.};
  Double_t parMpeak[1] = {1.};
  Double_t parMbkg[1]= {1.};
  
  return S*(promptfraction*GetImpParPromptProb(d0,parprompt)+(1-promptfraction)*GetImpParRecoFDProb(d0,parFD))*GetMassPeakProb(mass,parMpeak)+(norm-S)*GetImpParBkgProb(d0,parbkg)*GetMassBkgProb(mass,parMbkg);
  
}

