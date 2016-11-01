#ifndef ALIDPLUSCHARMFRACTIONIPFITTER_H
#define ALIDPLUSCHARMFRACTIONIPFITTER_H
/// \class AliDplusCharmFractionIPfitter
/// \class that performs the impact-parameter unbinned (binned) fit on a tree (histogram) for the determination of the fraction of prompt D+ 
/// \author Fabrizio Grosa <grosa@to.infn.it>, INFN Torino

#include <TObject.h>
#include <THnSparse.h>
#include <TH1F.h>
#include <TTree.h>
#include <TAxis.h>

class AliDplusCharmFractionIPfitter : public TObject {
  
public:
  AliDplusCharmFractionIPfitter();
  ~AliDplusCharmFractionIPfitter();

  enum SovT {kCentralValue,kLowerLimit,kLowerLimitOnlyStat,kLowerLimitOnlySyst,kUpperLimit,kUpperLimitOnlyStat,kUpperLimitOnlySyst};
  enum bkgtype {kSingleGaussExpo, kDoubleGaussExpo, kDoubleGaussExpoSymm};
  enum Fdtype {kConvolution, kGaussExpo};
  enum SBregion {kBoth, kLeft, kRight};
  enum TreeOrHisto {kTree, kHisto};
  enum Meson {kB, kD};
  enum GenFrom {kFromHisto, kFromFunc};
  
  void SetMCPromptSparse(THnSparseF* sparse) {fMCPromptSparse=(THnSparseF*)sparse->Clone();}
  void SetMCTrueFDSparse(THnSparseF* sparse) {fMCTrueFDSparse=(THnSparseF*)sparse->Clone();}
  void SetMCRecoFDSparse(THnSparseF* sparse) {fMCRecoFDSparse=(THnSparseF*)sparse->Clone();}
  void SetMCBkgSparse(THnSparseF* sparse) {fMCBkgSparse=(THnSparseF*)sparse->Clone();}
  void SetDataSparse(THnSparseF* sparse){fDataSparse=(THnSparseF*)sparse->Clone();}
  void SetDataTree(TTree* tree, TString massbranch="InvMass", TString ptbranch="Pt", TString impparbranch="d0") {fDataTree=(TTree*)tree->Clone(); fMassBranch=massbranch; fPtBranch=ptbranch; fImpParBranch=impparbranch;}

  void SetNumAxes(Int_t massaxis, Int_t ptaxis, Int_t impparaxis, Int_t PIDaxis) {fMassAxis=massaxis; fPtAxis=ptaxis; fImpParAxis=impparaxis; fPIDAxis = PIDaxis; SetPID(fPID);}
  void SetPtLims(Double_t ptmin, Double_t ptmax) {fPtMin=ptmin; fPtMax=ptmax;}
  void SetGenPromptFraction(Double_t fraction) {fGenPromptFraction=fraction;}
  void SetSgnAndBkg(Double_t sig, Double_t bkg) {fSig = sig; fBkg = bkg;}
  void SetNSigmas(Double_t nsigma) {fNSigma = nsigma;}
  void SetSignal(Int_t sig) {fSig = sig;}
  void SetBkg(Int_t bkg) {fBkg = bkg;}
  void SetFitOptions(TString fitoptions) {fFitOptions = fitoptions; if(!fFitOptions.Contains("0")) fFitOptions += "0";}
  void SetFitOptionsPrefit(TString fitoptions) {fFitOptionsPrefit = fitoptions;}
  void SetRebinImpParHistos(Int_t rebin) {fReb = rebin;}
  void SetNSigmaSBLimits(Int_t low, Int_t high) {fNSigmaSBLow=low;fNSigmaSBHigh=high;}
  void SetFDFunction(Int_t FDtype = kConvolution) {fFDType = FDtype;}
  void SetBkgFunction(Int_t bkgtype = kDoubleGaussExpoSymm) {fBkgType = bkgtype;}
  void SetSideBandsRegion(Int_t sbregion = kBoth) {fSBRegion = sbregion;}
  void SetPID(Bool_t isPIDon);
  void SetGaussOnlyForPrompt(Bool_t gaussonly) {fGaussOnlyForPrompt = gaussonly;}
  void SetBkgFromMC(Bool_t MCbkg) {fBkgFromMC = MCbkg;}
  void SetInitialParameters(Double_t* parprompt, Double_t* parFD, Double_t* parbkg) {
    fInitParamPrompt=parprompt;
    fInitParamFD=parFD;
    fInitParamBkg=parbkg;
  }
  void FixSigmaPromptFromMC(Bool_t isSigmaFixed=kFALSE) {fSigmaPromptFixed=isSigmaFixed;}
  void SetRelativeLimitSigmaPrompt(Double_t lim=0.2) {fSigmaPromptLim=lim;}
  void SetLimitsPromptFrac(Double_t min=0.2, Double_t max=1.5) {fPromptFracMin=min; fPromptFracMax=max;}

  void SetVariableBinningHisto(Bool_t varbin, Int_t ncountsmin) {fVariableBinning=varbin; fnCountsMin=ncountsmin;}
  void SetBkgSubtraction(Bool_t subbkg) {fSubBkg=subbkg;}
  
  void SetPtBWeightsHisto(TH1F *ptBweights) {fPtBWeightsHisto = (TH1F*)ptBweights->Clone(); fPtBWeight = kTRUE;}
  void SetPtDWeightsHisto(TH1F *ptDweights) {fPtDWeightsHisto = (TH1F*)ptDweights->Clone(); fPtDWeight = kTRUE;}
  void SetPtBWeights(Bool_t ptBweights) {fPtBWeight = ptBweights;}
  void SetPtDWeights(Bool_t ptDweights) {fPtDWeight = ptDweights;}
  
  void GetSignal(Int_t rebin=1, Int_t fsig=0, Int_t fbkg=0, Double_t min=1.68, Double_t max=2.05, Int_t sovert=kCentralValue);
  void CheckSideBandsImpParDist(Int_t nSigmaWidth);
  void PrefitOnPrompt(Double_t d0minPrompt=-1000, Double_t d0maxPrompt=1000);
  void PrefitOnFD(Double_t d0minFD=-1000,Double_t d0maxFD=1000);
  void PrefitOnBkg(Double_t d0minBkg=-1000,Double_t d0maxBkg=1000,Bool_t printSB=kTRUE);
  void PrefitStep(Double_t d0minPrompt=-1000,
                  Double_t d0maxPrompt=1000,
                  Double_t d0minFD=-1000,
                  Double_t d0maxFD=1000,
                  Double_t d0minBkg=-1000,
                  Double_t d0maxBkg=1000);
  void GenerateEntries(Bool_t background=kTRUE, Int_t gentype=kFromFunc, Int_t treeOrhisto=kTree);
  void FitTree(Double_t d0min=-1000, Double_t d0max=1000, Bool_t print=kTRUE);
  void FitHisto(Double_t d0min=-1000, Double_t d0max=1000, Bool_t print=kTRUE);
  
  void GetPromptFractionWithIPCut(Double_t d0cut, Double_t &promptfrac, Double_t &err, Bool_t genfrac=kFALSE);
  Double_t GetPromptFraction() {return fPromptFraction;}
  Double_t GetPromptFractionErr() {return fPromptFractionErr;}
  Double_t GetPromptFractionGauss() {return fPromptFractionGauss;}
  Double_t GetPromptFractionGaussErr() {return fPromptFractionGaussErr;}
  Double_t GetPromptMean() {return fPromptMean;}
  Double_t GetPromptMeanErr() {return fPromptMeanErr;}
  Double_t GetPromptSigma() {return fPromptSigma;}
  Double_t GetPromptSigmaErr() {return fPromptSigmaErr;}
  Double_t GetPromptSigmaMC() {return fPromptSigmaMC;}
  Double_t GetPromptSigmaMCErr() {return fPromptSigmaMCErr;}
  Double_t GetPromptLambda() {return fPromptLambda;}
  Double_t GetPromptLambdaErr() {return fPromptLambdaErr;}
  Double_t GetFDFraction1() {return fFDFraction1;}
  Double_t GetFDFraction1Err() {return fFDFraction1Err;}
  Double_t GetFDMean() {return fFDMean;}
  Double_t GetFDMeanErr() {return fFDMeanErr;}
  Double_t GetFDLambda1() {return fFDLambda1;}
  Double_t GetFDLambda1Err() {return fFDLambda1Err;}
  Double_t GetFDLambda2() {return fFDLambda2;}
  Double_t GetFDLambda2Err() {return fFDLambda2Err;}
  Double_t GetBkgFractionGauss1() {return fBkgFractionGauss1;}
  Double_t GetBkgFractionGauss1Err() {return fBkgFractionGauss1Err;}
  Double_t GetBkgMean1() {return fBkgMean1;}
  Double_t GetBkgMean1Err() {return fBkgMean1Err;}
  Double_t GetBkgSigma1() {return fBkgSigma1;}
  Double_t GetBkgSigma1Err() {return fBkgSigma1Err;}
  Double_t GetBkgLambda1() {return fBkgLambda1;}
  Double_t GetBkgLambda1Err() {return fBkgLambda1Err;}
  Double_t GetBkgFractionGauss2() {return fBkgFractionGauss2;}
  Double_t GetBkgFractionGauss2Err() {return fBkgFractionGauss2Err;}
  Double_t GetBkgMean2() {return fBkgMean2;}
  Double_t GetBkgMean2Err() {return fBkgMean2Err;}
  Double_t GetBkgSigma2() {return fBkgSigma2;}
  Double_t GetBkgSigma2Err() {return fBkgSigma2Err;}
  Double_t GetBkgLambda2() {return fBkgLambda2;}
  Double_t GetBkgLambda2Err() {return fBkgLambda2Err;}
  
  Double_t GetChi() {return fChiSquare;}
  Double_t GetRedChi() {return fChiSquare/fNDF;}
  Double_t GetDegOfFreedom() {return fNDF;}
  
private:
  
  ///fit functions
  Double_t Gauss(Double_t d0, Double_t mean, Double_t sigma);
  Double_t ExpoDouble(Double_t d0, Double_t mean, Double_t lambda);
  Double_t FunctionImpParPrompt(Double_t* x, Double_t* par);
  Double_t FunctionTrueImpParFD(Double_t* x, Double_t* par);
  Double_t FunctionRecoImpParFD(Double_t *x, Double_t *par);
  Double_t FunctionImpParBkg(Double_t *x, Double_t *par);
  Double_t Convolution(Double_t x, Double_t xmin, Double_t xmax, Int_t nsteps);
  Double_t FitFunction(Double_t* x, Double_t *par);
  
  ///derivative functions with respect to sigma (for the error propagation to fprompt within [-d0cut,d0cut])
  Double_t GaussSigmaDerivative(Double_t *d0,Double_t *pars);
  Double_t ConvolutionSigmaDerivative(Double_t x, Double_t xmin, Double_t xmax, Int_t nsteps);
  Double_t RecoFDSigmaDerivative(Double_t *d0, Double_t *pars);

  ///function to reweight impact parameter FD histos for the pT of the B mesons
  TH1F* GetPtreweightedHisto(THnSparseF* sparse, Int_t meson);
  
  ///function to draw the histogram with the fit functions
  void DrawResult(Bool_t isFromTree=kTRUE);

  ///function to rebin with variable bin width the impact-parameter histo
  void RebinVariableWidth();
  
  void ResetAxes(THnSparseF* sparse);
  TH1F* GetSidebandsDist(Bool_t printSB=kFALSE);
  void SetPtRange(THnSparseF* sparse);
  void SetMassRange(THnSparseF* sparse);
  
  ///data members
  THnSparseF* fMCPromptSparse; //!<! MC prompt sparse
  THnSparseF* fMCTrueFDSparse; //!<! MC FD sparse with true impact parameter
  THnSparseF* fMCRecoFDSparse; //!<! MC FD sparse with reco impact parameter 
  THnSparseF* fMCBkgSparse; //!<! MC Bkg sparse
  THnSparseF* fDataSparse; //!<! Data sparse
  TTree* fDataTree; //!<! Data tree (evenctually generated with MC toy)
  TH1F* fImpParHisto; //!<! Data histo with only impact parameter (eventually generated with MC toy)
  TString fMassBranch; ///mass branch name in data tree
  TString fPtBranch; ///pt branch name in data tree
  TString fImpParBranch; ///imp par branch name in data tree
  Double_t fPtMin; ///lower pT limit
  Double_t fPtMax; ///upper pT limit
  Double_t fMassMean; ///D+ mass obtained with massfit
  Double_t fMassSigma; ///D+ mass sigma obtained with mass fit
  Double_t fNSigma; ///number of sigma for the mass cut
  Double_t fGenPromptFraction; ///prompt fraction for tree generation
  TString fFitOptions; ///options for fit on data
  TString fFitOptionsPrefit; ///options for the prefits
  Int_t fReb; ///rebinning value for d0 histos
  Int_t fNSigmaSBLow;///lower value in numer of mass sigmas of SB region
  Int_t fNSigmaSBHigh;///upper value in numer of mass sigmas of SB region
  Int_t fMassAxis;///Mass axis in THnSparse
  Int_t fPtAxis;///pT axis in THnSparse
  Int_t fImpParAxis;///Imp Par axis in THnSparse
  Int_t fPIDAxis;///PID axis in THnSparses
  Int_t fPIDbin;///PID bin in the PID axis
  Int_t fPID;///flag to activate PID

  //initial parameters
  Double_t *fInitParamPrompt; ///prompt template initial parameters
  Double_t *fInitParamFD; ///FD template initial parameters
  Double_t *fInitParamBkg; ///bkg template initial parameters
 
  //parameters from mass fit
  Double_t fSig; ///signal obtained with mass fit (n sigma)
  Double_t fSigErr; ///signal obtained with mass fit (n sigma)
  Double_t fSigErrStat; ///statistic error on fSig
  Double_t fSigRelErrSyst; ///relative systematic error on fSig
  Double_t fBkg; ///background obtained with mass fit (n sigma)
  Double_t fBkgErrStat; ///statistic error on fBkg
  Double_t fIntegral; ///number of entries within the mass range (signal+bkg)
  
  //fit parameters
  Double_t fPromptFraction; ///prompt fraction
  Double_t fPromptFractionGauss;///
  Double_t fPromptMean;///
  Double_t fPromptSigma;///
  Double_t fPromptSigmaMC;///
  Double_t fPromptLambda;///
  Double_t fFDFraction1;///
  Double_t fFDMean;///
  Double_t fFDLambda1;///
  Double_t fFDLambda2;///
  Double_t fBkgFractionGauss1;///
  Double_t fBkgMean1;///
  Double_t fBkgSigma1;///
  Double_t fBkgLambda1;///
  Double_t fBkgFractionGauss2;///
  Double_t fBkgMean2;///
  Double_t fBkgSigma2;///
  Double_t fBkgLambda2;///
  Double_t fBkgFracFunc1;///

  //limits & constraints
  Bool_t fSigmaPromptFixed; ///fix sigma prompt from MC in fit on data
  Double_t fSigmaPromptLim; ///set relative limit (w.r.t. MC value) for the sigma prompt parameter
  Double_t fPromptFracMin; ///lower limit for the prompt fractioin
  Double_t fPromptFracMax; ///upper limit for the prompt fractioin
  
  //chi square and NDF
  Double_t fChiSquare;///
  Double_t fNDF;///
  
  //error on fit parameters
  Double_t fPromptFractionErr;///
  Double_t fPromptFractionGaussErr;///
  Double_t fPromptMeanErr;///
  Double_t fPromptSigmaErr;///
  Double_t fPromptSigmaMCErr;///
  Double_t fPromptLambdaErr;///
  Double_t fFDFraction1Err;///
  Double_t fFDMeanErr;///
  Double_t fFDLambda1Err;///
  Double_t fFDLambda2Err;///
  Double_t fBkgFractionGauss1Err;///
  Double_t fBkgMean1Err;///
  Double_t fBkgSigma1Err;///
  Double_t fBkgLambda1Err;///
  Double_t fBkgFractionGauss2Err;///
  Double_t fBkgMean2Err;///
  Double_t fBkgSigma2Err;///
  Double_t fBkgLambda2Err;//
  Double_t fBkgFracFunc1Err;///
  
  Double_t fCovFracSigmaPrompt; ///covariance between fprompt and sigmaprompt obtained from the fit
  
  Bool_t fMCTest;///flag or MC tests or real data fit
  Bool_t fGaussOnlyForPrompt;///flag to set the prefit on prompt distribution only with a gaussian
  Bool_t fBkgFromMC;///flag for taking the bkg parametrization from MC instead of SB
  
  Int_t fSBRegion;///only left, only right, or both
  Int_t fBkgType;///gauss+expo, double gauss+expo, double gauss+expo with same parameters
  Int_t fFDType;///convolution or gauss+expo

  //pTB reweight
  Bool_t fPtBWeight;///flag to activate the reweight for pT of the B mesons for FD
  Bool_t fPtDWeight;///flag to activate the reweight for pT of the D mesons for prompt
  TH1F* fPtBWeightsHisto;//!<! pTB weights histogram
  TH1F* fPtDWeightsHisto;//!<! pTD weights histogram

  //variable binning 
  Bool_t fVariableBinning;///flag to activate the variable binning histogram
  Int_t fnCountsMin;///flag to activate the minimum number of entries in each bin in case of variable binning

  //background subtraction
  Bool_t fSubBkg;///flag to activate the background subtraction in case of binned fit
  
  /// \cond CLASSDEF
  ClassDef(AliDplusCharmFractionIPfitter,1);
  /// \endcond
};
#endif //ALIDPLUSCHARMFRACTIONIPFITTER
