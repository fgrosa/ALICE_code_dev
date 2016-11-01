/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//*************************************************************************
/// \class Class AliAnalysisTaskDmesonKF
/// \brief AliAnalysisTaskSE for D0 and D+ mesons using the Kalman Filter
/// \author Fabrizio Grosa, INFN Turin, grosa@to.infn.it
/// \date: April 2016
//*************************************************************************

#ifndef HomogeneousField 
#define HomogeneousField
#endif

#define KFParticleStandalone

#include <Riostream.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TVector3.h>

#include "AliAnalysisTaskDmesonKF.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliCentrality.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliKFVertex.h"
#include "AliESDtrack.h"
#include "AliESDUtils.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisUtils.h"
#include "AliVertexerTracks.h"
#include "AliAODVertex.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonKF);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskDmesonKF::AliAnalysisTaskDmesonKF():
  AliAnalysisTaskSE()
  ,fOutput(0x0)
  ,fHistNEvents(0x0)
  ,fMCGenAccPrompt(0x0)
  ,fMCGenAccFD(0x0)
  ,fHistPVres(0x0)
  ,fHistPVpulls(0x0)
  ,fHistPVchiS(0x0)
  ,fHistPVprob(0x0)
  ,fHistPVHFres(0x0)
  ,fHistPVHFpulls(0x0)
  ,fHistPVHFchiS(0x0)
  ,fHistPVHFprob(0x0)
  ,fHistPVHFremDaures(0x0)
  ,fHistPVHFremDaupulls(0x0)
  ,fHistPVHFremDauchiS(0x0)
  ,fHistPVHFremDauprob(0x0)
  ,fHistPVESDres(0x0)
  ,fHistPVESDpulls(0x0)
  ,fHistPVESDchiS(0x0)
  ,fHistPVESDprob(0x0)
  ,fHistPVESDHFres(0x0)
  ,fHistPVESDHFpulls(0x0)
  ,fHistPVESDHFchiS(0x0)
  ,fHistPVESDHFprob(0x0)
  ,fHistPVESDHFremDaures(0x0)
  ,fHistPVESDHFremDaupulls(0x0)
  ,fHistPVESDHFremDauchiS(0x0)
  ,fHistPVESDHFremDauprob(0x0)
  ,fHistMres(0x0)
  ,fHistMpulls(0x0)
  ,fHistMresTopo(0x0)
  ,fHistMpullsTopo(0x0)
  ,fHistMresMass(0x0)
  ,fHistMpullsMass(0x0)
  ,fHistMchiS(0x0)
  ,fHistMprob(0x0)
  ,fHistMchiSTopo(0x0)
  ,fHistMprobTopo(0x0)
  ,fHistMchiSMass(0x0)
  ,fHistMprobMass(0x0)
  ,fHistAliVertMres(0x0)
  ,fHistAliVertMpulls(0x0)
  ,fHistDres(0x0)
  ,fHistDpulls(0x0)
  ,fHistDchiS(0x0)
  ,fHistDprob(0x0)
  ,fHistTrackres(0x0)
  ,fHistTrackpulls(0x0)
  ,fReadMC(kTRUE)
  ,fQAonKF(kTRUE)
  ,fPDGcode(421)
  ,fMassMean(0)
  ,fMassMin(0)
  ,fMassMax(0)
  ,fPtMin(0.2)
  ,fNDau(2)
  ,fESDtrackCuts(0x0)
  ,fPIDResponse(0x0)
  ,fMaxRapidCand(-999)
  ,fSearchUpToQuark(kTRUE)
  ,fTriggerMask(AliVEvent::kINT7)
  ,fUsePID(kTRUE)
  ,fUseStrongPID(kTRUE)
  ,fMaxPtForStrongPID(2.)
  ,fOptPileup(kRejectPileupEvent) 
  ,fMaxVtxZ(10)
  ,fMinVtxType(3) 
  ,fMinVtxContr(1)   
  ,fCutOnzVertexSPD(0)
  ,fMinNClustersTPCPID(0)
  ,fCutTOFmismatch(0.01)
  ,fLooseCuts(kFALSE)
  ,fPtLims(0)
  ,fMinDauPt(0) 
  ,fCospMin(0)
  ,fDecLMin(0)
  ,fNDecLXYMin(0)
  ,fSigVtxMax(0)
  ,fChiMax(0)
  ,fPVChiMax(0)
  ,fPosTracksArray(0x0)
  ,fNegTracksArray(0x0)
  ,fTransportToReco(kFALSE)
{
  /// Default constructor
  
  fMassMean = TDatabasePDG::Instance()->GetParticle(fPDGcode)->Mass();
  fMassMin = fMassMean-0.15;
  fMassMax = fMassMean+0.15;
  fTriggerClass[0] = "";
  fTriggerClass[1] = "";

  for(Int_t iSparse=0; iSparse<3; iSparse++) {
    fSparse[iSparse] = 0x0;
  }

  fPtLims.push_back(0);
  fPtLims.push_back(100000000000.);
  fCospMin.push_back(0.9);
  fDecLMin.push_back(0.);
  fNDecLXYMin.push_back(1.);
  fSigVtxMax.push_back(0.06);
}

//________________________________________________________________________
AliAnalysisTaskDmesonKF::AliAnalysisTaskDmesonKF(const char *name, Int_t meson, Bool_t readMC):
  AliAnalysisTaskSE(name)
  ,fOutput(0x0)
  ,fHistNEvents(0x0)
  ,fMCGenAccPrompt(0x0)
  ,fMCGenAccFD(0x0)
  ,fHistPVres(0x0)
  ,fHistPVpulls(0x0)
  ,fHistPVchiS(0x0)
  ,fHistPVprob(0x0)
  ,fHistPVHFres(0x0)
  ,fHistPVHFpulls(0x0)
  ,fHistPVHFchiS(0x0)
  ,fHistPVHFprob(0x0)
  ,fHistPVHFremDaures(0x0)
  ,fHistPVHFremDaupulls(0x0)
  ,fHistPVHFremDauchiS(0x0)
  ,fHistPVHFremDauprob(0x0)
  ,fHistPVESDres(0x0)
  ,fHistPVESDpulls(0x0)
  ,fHistPVESDchiS(0x0)
  ,fHistPVESDprob(0x0)
  ,fHistPVESDHFres(0x0)
  ,fHistPVESDHFpulls(0x0)
  ,fHistPVESDHFchiS(0x0)
  ,fHistPVESDHFprob(0x0)
  ,fHistPVESDHFremDaures(0x0)
  ,fHistPVESDHFremDaupulls(0x0)
  ,fHistPVESDHFremDauchiS(0x0)
  ,fHistPVESDHFremDauprob(0x0)
  ,fHistMres(0x0)
  ,fHistMpulls(0x0)
  ,fHistMresTopo(0x0)
  ,fHistMpullsTopo(0x0)
  ,fHistMresMass(0x0)
  ,fHistMpullsMass(0x0)
  ,fHistMchiS(0x0)
  ,fHistMprob(0x0)
  ,fHistMchiSTopo(0x0)
  ,fHistMprobTopo(0x0)
  ,fHistMchiSMass(0x0)
  ,fHistMprobMass(0x0)
  ,fHistAliVertMres(0x0)
  ,fHistAliVertMpulls(0x0)
  ,fHistDres(0x0)
  ,fHistDpulls(0x0)
  ,fHistDchiS(0x0)
  ,fHistDprob(0x0)
  ,fHistTrackres(0x0)
  ,fHistTrackpulls(0x0)
  ,fReadMC(readMC)
  ,fQAonKF(kTRUE)
  ,fPDGcode(421)
  ,fMassMean(0)
  ,fMassMin(0)
  ,fMassMax(0)
  ,fPtMin(0.2)
  ,fNDau(2)
  ,fESDtrackCuts(0x0)
  ,fPIDResponse(0x0)
  ,fMaxRapidCand(-999)
  ,fSearchUpToQuark(kTRUE)
  ,fTriggerMask(AliVEvent::kINT7)
  ,fUsePID(kTRUE)
  ,fUseStrongPID(kTRUE)
  ,fMaxPtForStrongPID(2.)
  ,fOptPileup(kRejectPileupEvent) 
  ,fMaxVtxZ(10)
  ,fMinVtxType(3) 
  ,fMinVtxContr(1)  
  ,fCutOnzVertexSPD(0)
  ,fMinNClustersTPCPID(0)
  ,fCutTOFmismatch(0.01)
  ,fLooseCuts(kFALSE)
  ,fPtLims(0)
  ,fMinDauPt(0) 
  ,fCospMin(0)
  ,fDecLMin(0)
  ,fNDecLXYMin(0)
  ,fSigVtxMax(0)
  ,fChiMax(0)
  ,fPVChiMax(0)
  ,fPosTracksArray(0x0)
  ,fNegTracksArray(0x0)
  ,fTransportToReco(kFALSE)
{
  /// Standard constructor

  if(meson==kDplusToKpipi) {
    fPDGcode=411;
    fNDau = 3;
  }
  
  fMassMean = TDatabasePDG::Instance()->GetParticle(fPDGcode)->Mass();
  fMassMin = fMassMean-0.15;
  fMassMax = fMassMean+0.15;
  fTriggerClass[0] = "CINT7";

  for(Int_t iSparse=0; iSparse<3; iSparse++) {
    fSparse[iSparse] = 0x0;
  }

  fPtLims.push_back(0);
  fPtLims.push_back(100000000000.);
  fCospMin.push_back(0.9);
  fDecLMin.push_back(0.);
  fNDecLXYMin.push_back(1.);
  fSigVtxMax.push_back(0.06);

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskDmesonKF::~AliAnalysisTaskDmesonKF()
{
  /// Destructor

  if(fPosTracksArray) delete fPosTracksArray;
  if(fNegTracksArray) delete fNegTracksArray;

  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    for(Int_t iSparse=0;iSparse<3;iSparse++){
      delete fSparse[iSparse];
    }
    delete fHistPVres;
    delete fHistPVpulls;
    delete fHistPVchiS;
    delete fHistPVprob;
    delete fHistPVHFres;
    delete fHistPVHFpulls;
    delete fHistPVHFchiS;
    delete fHistPVHFprob;
    delete fHistPVHFremDaures;
    delete fHistPVHFremDaupulls;
    delete fHistPVHFremDauchiS;
    delete fHistPVHFremDauprob;
    delete fHistPVESDres;
    delete fHistPVESDpulls;
    delete fHistPVESDchiS;
    delete fHistPVESDprob;
    delete fHistPVESDHFres;
    delete fHistPVESDHFpulls;
    delete fHistPVESDHFchiS;
    delete fHistPVESDHFprob;
    delete fHistPVESDHFremDaures;
    delete fHistPVESDHFremDaupulls;
    delete fHistPVESDHFremDauchiS;
    delete fHistPVESDHFremDauprob;
    delete fMCGenAccPrompt;
    delete fMCGenAccFD;
    delete fHistPVres;
    delete fHistPVpulls;
    delete fHistPVchiS;
    delete fHistPVprob;
    delete fHistMres;
    delete fHistMpulls;
    delete fHistMresTopo;
    delete fHistMpullsTopo;
    delete fHistMresMass;
    delete fHistMpullsMass;
    delete fHistMchiS;
    delete fHistMprob;
    delete fHistMchiSTopo;
    delete fHistMprobTopo;
    delete fHistMchiSMass;
    delete fHistMprobMass;
    delete fHistAliVertMres;
    delete fHistAliVertMpulls;
    delete fHistDres;
    delete fHistDpulls;
    delete fHistDchiS;
    delete fHistDprob;
    delete fHistTrackres;
    delete fHistTrackpulls;
    delete fOutput;   
  }
  if(fESDtrackCuts) delete fESDtrackCuts;
}  

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::UserCreateOutputObjects()
{
  ///tracks arrays
  fPosTracksArray = new TObjArray(10);
  fNegTracksArray = new TObjArray(10);

  ///Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("fHistNEvents", "number of events ",8,-0.5,7.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvents accepted");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"PV reco");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"no PV reco");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"PV outside range");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"PV with too few tracks");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"rejected due to pile-up");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"rejected due to PhysSel");

  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);  
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  CreateSparses();
  if(fReadMC) CreateMCGenAccHistos();
  if(fReadMC && fQAonKF) CreateQAhistos();

  PostData(1,fOutput);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::UserExec(Option_t */*option*/)
{
  ///==============ESD EVENT=============== 
  AliESDEvent *esdEvent = (AliESDEvent*) (InputEvent());
  
  if(!esdEvent) {
    printf("AliAnalysisTaskSDDRP::Exec(): bad ESD\n");
    return;
  } 

  if(!fESDtrackCuts){
    Int_t year=2011;
    if(esdEvent->GetRunNumber()<=139517) year=2010;
    if(year==2010) fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE);
    else fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
    fESDtrackCuts->SetMaxDCAToVertexXY(2.4);
    fESDtrackCuts->SetMaxDCAToVertexZ(3.2);
    fESDtrackCuts->SetDCAToVertex2D(kTRUE);
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* eventHandler = (AliInputEventHandler*)mgr->GetInputEventHandler();
  fPIDResponse=eventHandler->GetPIDResponse();  

  ///refit with multivertexer 
  AliESDUtils::RefitESDVertexTracks(esdEvent,6);

  Double_t Bz = esdEvent->GetMagneticField();
  KFParticle::SetField(Bz);
  const AliESDVertex *ESDprimary = (AliESDVertex*)esdEvent->GetPrimaryVertex(); 
  //if(!ESDprimary || TMath::Abs(Bz)<0.001) return; //?

  fHistNEvents->Fill(0);

  ///==============MC INFO===============  
  AliStack *stack=0x0;
  AliMCEvent* mcEvent=0x0;
  if(fReadMC){
   AliMCEventHandler* MCeventHandler = dynamic_cast<AliMCEventHandler*> (mgr->GetMCtruthEventHandler());
    if (!MCeventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }
    mcEvent = MCeventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
    stack = mcEvent->Stack();
    if (!stack) {
      Printf("ERROR: stack not available");
      return;
    }
  }
  
  ///==============EVENT SELECTION===============  
  Bool_t isAccepted=kTRUE;
  Int_t whyRejection=0;

  TString firedTriggerClasses=((AliESDEvent*)esdEvent)->GetFiredTriggerClasses();
  //don't do for MC and for PbPb 2010 data
  if(!fReadMC && (esdEvent->GetRunNumber()<136851 && esdEvent->GetRunNumber()>139517)) {
    if(firedTriggerClasses.Contains(fTriggerClass[0].Data()) && 
       (fTriggerClass[1].CompareTo("")==0 || firedTriggerClasses.Contains(fTriggerClass[1].Data()))) {
      whyRejection=5;
      isAccepted=kFALSE;
    }
  }

  ///physics selection
  UInt_t maskPhysSel = eventHandler->IsEventSelected();
  Bool_t isSelected = (maskPhysSel & fTriggerMask);
  if(!isSelected) {
    if(isAccepted) whyRejection=7;
    isAccepted=kFALSE;
  }

 ///vtx requirements
  if(!ESDprimary) {
    isAccepted=kFALSE;
  }
  else{
    TString title=ESDprimary->GetTitle();
    if(title.Contains("Z") && fMinVtxType>1){
      isAccepted=kFALSE;
    }
    else if(title.Contains("3D") && fMinVtxType>2){
      isAccepted=kFALSE;
    }
    if(ESDprimary->GetNContributors()<fMinVtxContr){
      if(isAccepted) whyRejection=2;
      isAccepted=kFALSE;
    }
    if(TMath::Abs(ESDprimary->GetZ())>fMaxVtxZ) {
      if(isAccepted) whyRejection=6;
      isAccepted=kFALSE;
    } 
  }
  if(fCutOnzVertexSPD>0){
    const AliVVertex *vSPD = ((AliESDEvent*)esdEvent)->GetPrimaryVertexSPD();
    if(!vSPD || (vSPD && vSPD->GetNContributors()<fMinVtxContr)){
      isAccepted=kFALSE;
    }
    else{
      if(fCutOnzVertexSPD==1 && TMath::Abs(vSPD->GetZ())>12.) {
	if(isAccepted) whyRejection=6;
	isAccepted=kFALSE;
      } 
      if(fCutOnzVertexSPD==2 && ESDprimary){
	if(TMath::Abs(vSPD->GetZ()-ESDprimary->GetZ())>0.5) {
	  if(isAccepted) whyRejection=6;
	  isAccepted=kFALSE;
	} 
      }
    }
  }
  
  ///pileup rejection
  if(!fReadMC) {
    if(fOptPileup==kRejectPileupEvent){
      if(esdEvent->IsPileupFromSPD(5,0.8,3.,2.,10.)){//hard coded (can be changed)
	if(isAccepted) whyRejection=1;
	isAccepted=kFALSE;
      }
    }
    else if(fOptPileup==kRejectMVPileupEvent){
      AliAnalysisUtils utils;
      Bool_t isPUMV = utils.IsPileUpMV(esdEvent);
      if(isPUMV) {
	if(isAccepted) whyRejection=1;
	isAccepted=kFALSE;
      }
    }
  }

  ///protection for events with empty trigger mask in p-Pb
  if(fReadMC) {
    if(esdEvent->GetTriggerMask()==0 && (esdEvent->GetRunNumber()>=195344 && esdEvent->GetRunNumber()<=195677))
      return; 
  }
  ///==============FILL MC HISTOS FOR EFFICIENCY===============  
  if(fReadMC) FillMCGenAccHistos(stack);  

  //post data already here
  PostData(1,fOutput);
  
  if(isAccepted || (!isAccepted && whyRejection==6)) fHistNEvents->Fill(2); 
  if(!isAccepted && (whyRejection==0 || whyRejection==2)) fHistNEvents->Fill(3); 
  if(!isAccepted && whyRejection==6) fHistNEvents->Fill(4); 
  if(!isAccepted &&  whyRejection==2) fHistNEvents->Fill(5); 
  if(!isAccepted &&  whyRejection==1) fHistNEvents->Fill(6); 
  if(!isAccepted &&  whyRejection==7) fHistNEvents->Fill(7); 

  if(!isAccepted) return;
  
  fHistNEvents->Fill(1);

  ///==============QA ANALYSIS ON KF===============  
  if(fReadMC && fQAonKF) 
    DoQAanalysis(esdEvent,stack);

  ///==============CONVERT PV===============  
  KFVertex pvKF;
  KFPVertex pvertex;
  ConvertVertex(ESDprimary, pvertex);
  pvKF = KFVertex(pvertex);

  ///==============D MESONS ANALYSIS===============  
  ///get D meson candidates
  GetTracksArrays(esdEvent);

  ///combine tracks to get the D candidates and fill sparses
  if(fPDGcode==421) {
    Combine2ProngsKF(pvKF,stack);
  }
  else {
    Combine3ProngsKF(pvKF,stack);
  }
  
  fPosTracksArray->Clear();
  fNegTracksArray->Clear();

  PostData(1,fOutput); 
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  if(fDebug > 1) printf("AnalysisTaskDmesonKF: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
  if(fHistNEvents){
    printf("Number of analyzed events = %d\n",(Int_t)fHistNEvents->GetBinContent(0));
  }else{
    printf("ERROR: fHistNEvents not available\n");
    return;
  }
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::CreateSparses()
{
  /// Sparses with physical variables
  ///dimensions for THnSparses
  Int_t nmassbins=(fMassMax-fMassMin)*1000/2;
  
  Int_t nptbins=80;
  Double_t ptmin=0.;
  Double_t ptmax=40.;
 
  Int_t nchibins=200;
  Double_t chimin=0.;
  Double_t chimax=50.;

  Int_t nd0XYbins=200;
  Double_t d0XYmin=-1000.;
  Double_t d0XYmax=1000.;

  Int_t nDLbins;
  Double_t DLmin=0.;
  Double_t DLmax;
  Int_t nDLXYbins;
  Double_t DLXYmin=0.;
  Double_t DLXYmax;
  if(fPDGcode==411) {
    nDLbins=50;
    nDLXYbins=50;
    DLmax=0.5;
    DLXYmax=0.5;
  }
  else {
    nDLbins=50;
    nDLXYbins=50;
    DLmax=0.5;
    DLXYmax=0.5;
  }

  Int_t nNDLXYbins=40;
  Double_t NDLXYmin=0.;
  Double_t NDLXYmax=20.;

  Int_t nCbins;
  Double_t Cmin;
  Double_t Cmax;
  Int_t nCXYbins;
  Double_t CXYmin;
  Double_t CXYmax;

  if(fLooseCuts) {
    nCbins=100;
    nCXYbins=100;
    Cmin=0.9;
    Cmax=1.;
    CXYmin=0.9;
    CXYmax=1.;
  }
  else {
    nCbins=200;
    nCXYbins=200;
    Cmin=-1.;
    Cmax=1.;
    CXYmin=-1.;
    CXYmax=1.;
  }
  Int_t nDauPtbins=40;
  Double_t DauPtmin=0.;
  Double_t DauPtmax=2.;

  Int_t nSigVtx;
  Double_t SigVtxmin;
  Double_t SigVtxmax;
  if(fLooseCuts) {
    nSigVtx=60;
    SigVtxmin=0.;
    SigVtxmax=0.06;
  }
  else {
    nSigVtx=100;
    SigVtxmin=0.;
    SigVtxmax=0.1;
  }
  Int_t nd0d0exp = 40;
  Double_t d0d0expmin = -10.;
  Double_t d0d0expmax = 10.;

  Int_t nPID;
  Double_t PIDmin=-0.5;
  Double_t PIDmax;

  if(fPDGcode==421) {
    nPID=4;
    PIDmax=3.5;
  }
  else {
    nPID=2;
    PIDmax=1.5;
  }

  TString axTit[fNVarsForSparses];
  if(fPDGcode==421)
    axTit[0]="M_{K#pi} (GeV/c^{2})";
  else
    axTit[0]="M_{K#pi#pi} (GeV/c^{2})";
  
  axTit[1]="#it{p}_{T} (GeV/c)";
  axTit[2]="ImpPar XY (#mum)";  	
  axTit[3]="PID";
  axTit[4]="decL (cm)";  
  axTit[5]="decL XY (cm)";  
  axTit[6]="norm. decL XY";  
  axTit[7]="Cos(#theta_{P})";  
  axTit[8]="Cos(#theta_{P}^{XY})"; 
  axTit[9]="min daughter #it{p}_{T} (GeV/c)"; 
  axTit[10]="sigma vtx (cm)"; 
  axTit[11]="Norm max d0-d0exp"; 
  axTit[12]="KF #chi^{2}";  
  axTit[13]="PV #chi^{2}";  

  Int_t nbins[fNVarsForSparses]={nmassbins,nptbins,nd0XYbins,nPID,nDLbins,nDLXYbins,nNDLXYbins,nCbins,nCXYbins,nDauPtbins,nSigVtx,nd0d0exp,nchibins,nchibins};
  Double_t min[fNVarsForSparses]={fMassMin,ptmin,d0XYmin,PIDmin,DLmin,DLXYmin,NDLXYmin,Cmin,CXYmin,DauPtmin,SigVtxmin,d0d0expmin,chimin,chimin};
  Double_t max[fNVarsForSparses]={fMassMax,ptmax,d0XYmax,PIDmax,DLmax,DLXYmax,NDLXYmax,Cmax,CXYmax,DauPtmax,SigVtxmax,d0d0expmax,chimax,chimax};
  
  fSparse[0]=new THnSparseF("fSparseAll","Mass vs. pt vs. #chi^{2} vs. all - All",fNVarsForSparses,nbins,min,max);
  fSparse[1]=new THnSparseF("fSparsePrompt","Mass vs. pt vs. #chi^{2} vs. all - Prompt",fNVarsForSparses,nbins,min,max);
  fSparse[2]=new THnSparseF("fSparseFD","Mass vs. pt vs. #chi^{2} vs. all - FD",fNVarsForSparses,nbins,min,max);  

  for(Int_t iSparse=0; iSparse<3; iSparse++){
    for(Int_t iAxis=0; iAxis<fNVarsForSparses; iAxis++) fSparse[iSparse]->GetAxis(iAxis)->SetTitle(axTit[iAxis].Data());
    fSparse[iSparse]->GetAxis(fNVarsForSparses-1)->SetBinLabel(1,"rejected");
    if(fPDGcode==421) {
      fSparse[iSparse]->GetAxis(fNVarsForSparses-1)->SetBinLabel(2,"D^{0}");
      fSparse[iSparse]->GetAxis(fNVarsForSparses-1)->SetBinLabel(3,"#overline{D}^{0}");
      fSparse[iSparse]->GetAxis(fNVarsForSparses-1)->SetBinLabel(4,"D^{0} or #overline{D}^{0}");
    }
    else 
      fSparse[iSparse]->GetAxis(fNVarsForSparses-1)->SetBinLabel(2,"D^{#pm}");
    fSparse[iSparse]->GetAxis(fNVarsForSparses-1)->SetNdivisions(1,kFALSE);
    
  fOutput->Add(fSparse[iSparse]);
  }  
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::CreateMCGenAccHistos()
{
  /// MC gen acc step histos

  const Int_t nVars = 2;

  Int_t nptbins=80;
  Double_t ptmin=0.;
  Double_t ptmax=40.;

  Int_t nybins=100;
  Double_t ymin=-1.;
  Double_t ymax=1.;

  Int_t nbins[nVars]={nptbins,nybins};
  Double_t min[nVars]={ptmin,ymin};
  Double_t max[nVars]={ptmax,ymax};
  TString axTit[nVars] = {"#it{p}_{T} (GeV/c)","y"};
 
  fMCGenAccPrompt=new THnSparseF("fMCGenAccPrompt","MCGenAcc - Prompt",nVars,nbins,min,max);
  fMCGenAccFD=new THnSparseF("fMCGenAccFD","MCGenAcc - Prompt",nVars,nbins,min,max);

  for(Int_t iAxis=0; iAxis<nVars; iAxis++) {
    fMCGenAccPrompt->GetAxis(iAxis)->SetTitle(axTit[iAxis].Data());
    fMCGenAccFD->GetAxis(iAxis)->SetTitle(axTit[iAxis].Data());
  }

  fOutput->Add(fMCGenAccPrompt);
  fOutput->Add(fMCGenAccFD);
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::CreateQAhistos()
{
  /// QA histograms
  Int_t nbins = 100;
  Double_t minres = -1000;
  Double_t maxres = 1000;
  Double_t minresPV = -600;
  Double_t maxresPV = 600;
  Double_t minrestopo = -400;
  Double_t maxrestopo = 400;
  Double_t minresE = -0.1;
  Double_t maxresE = 0.1;
  Double_t minresp = -0.1;
  Double_t maxresp = 0.1;
  Double_t minpull = -6.;
  Double_t maxpull = 6.;
  Double_t minchi = 0.;
  Double_t maxchi = 25.;
  Double_t minmult = 0.;
  Double_t maxmult = 100.;

  Int_t PVbins[4] = {nbins,nbins,nbins,maxmult};
  Double_t PVminres[4] = {minresPV,minresPV,minresPV,minmult-0.5};
  Double_t PVmaxres[4] = {maxresPV,maxresPV,maxresPV,maxmult-0.5};
  Double_t PVminpulls[4] = {minpull,minpull,minpull,minmult-0.5};
  Double_t PVmaxpulls[4] = {maxpull,maxpull,maxpull,maxmult-0.5};
  TString PVrestit[4] = {"X_{reco}-X_{MC} (#mum)","Y_{reco}-Y_{MC} (#mum)","Z_{reco}-Z_{MC} (#mum)","N_{tracks}"}; 
  TString PVpullstit[4] = {"(X_{reco}-X_{MC})/#sigma_{X_{reco}}","(Y_{reco}-Y_{MC})/#sigma_{Y_{reco}}","(Z_{reco}-Z_{MC})/#sigma_{Z_{reco}}","N_{tracks}"}; 

  fHistPVres = new THnSparseF("fHistPVres","PV - res",4,PVbins,PVminres,PVmaxres);
  fHistPVpulls = new THnSparseF("fHistPVpulls","PV - pulls",4,PVbins,PVminpulls,PVmaxpulls);
  fHistPVHFres = new THnSparseF("fHistPVHFres","PV - res",4,PVbins,PVminres,PVmaxres);
  fHistPVHFpulls = new THnSparseF("fHistPVHFpulls","PV - pulls",4,PVbins,PVminpulls,PVmaxpulls);
  fHistPVHFremDaures = new THnSparseF("fHistPVHFremDaures","PV dau removed & mother added- res",4,PVbins,PVminres,PVmaxres);
  fHistPVHFremDaupulls = new THnSparseF("fHistPVHFremDaupulls","PV dau removed & mother added - pulls",4,PVbins,PVminpulls,PVmaxpulls);
  for(Int_t iAxis=0; iAxis<4; iAxis++) {
    fHistPVres->GetAxis(iAxis)->SetTitle(PVrestit[iAxis].Data());
    fHistPVpulls->GetAxis(iAxis)->SetTitle(PVpullstit[iAxis].Data());
    fHistPVHFres->GetAxis(iAxis)->SetTitle(PVrestit[iAxis].Data());
    fHistPVHFpulls->GetAxis(iAxis)->SetTitle(PVpullstit[iAxis].Data());
    fHistPVHFremDaures->GetAxis(iAxis)->SetTitle(PVrestit[iAxis].Data());
    fHistPVHFremDaupulls->GetAxis(iAxis)->SetTitle(PVpullstit[iAxis].Data());
   }
  fHistPVchiS = new TH1F("fHistPVchiS","PV - #chi^{2}",nbins,minchi,maxchi);
  fHistPVchiS->GetXaxis()->SetTitle("#chi^{2}");
  fHistPVprob = new TH1F("fHistPVprob","PV - probability",nbins,0,1);
  fHistPVprob->GetXaxis()->SetTitle("prob");
  fHistPVHFchiS = new TH1F("fHistPVHFchiS","PV - #chi^{2}",nbins,minchi,maxchi);
  fHistPVHFchiS->GetXaxis()->SetTitle("#chi^{2}");
  fHistPVHFprob = new TH1F("fHistPVHFprob","PV - probability",nbins,0,1);
  fHistPVHFprob->GetXaxis()->SetTitle("prob");
  fHistPVHFremDauchiS = new TH1F("fHistPVHFremDauchiS","PV dau removed & mother added - #chi^{2}",nbins,minchi,maxchi);
  fHistPVHFremDauchiS->GetXaxis()->SetTitle("#chi^{2}");
  fHistPVHFremDauprob = new TH1F("fHistPVHFremDauprob","PV dau removed & mother added - probability",nbins,0,1);
  fHistPVHFremDauprob->GetXaxis()->SetTitle("prob");
  
  fHistPVESDres = new THnSparseF("fHistPVESDres","ESD PV - res",4,PVbins,PVminres,PVmaxres);
  fHistPVESDpulls = new THnSparseF("fHistPVESDpulls","ESD PV - pulls",4,PVbins,PVminpulls,PVmaxpulls);
  fHistPVESDHFres = new THnSparseF("fHistPVESDHFres","ESD PV - res",4,PVbins,PVminres,PVmaxres);
  fHistPVESDHFpulls = new THnSparseF("fHistPVESDHFpulls","ESD PV - pulls",4,PVbins,PVminpulls,PVmaxpulls);
  fHistPVESDHFremDaures = new THnSparseF("fHistPVESDHFremDaures","ESD PV dau removed & mother added- res",4,PVbins,PVminres,PVmaxres);
  fHistPVESDHFremDaupulls = new THnSparseF("fHistPVESDHFremDaupulls","ESD PV dau removed & mother added - pulls",4,PVbins,PVminpulls,PVmaxpulls);
  for(Int_t iAxis=0; iAxis<4; iAxis++) {
    fHistPVESDres->GetAxis(iAxis)->SetTitle(PVrestit[iAxis].Data());
    fHistPVESDpulls->GetAxis(iAxis)->SetTitle(PVpullstit[iAxis].Data());
    fHistPVESDHFres->GetAxis(iAxis)->SetTitle(PVrestit[iAxis].Data());
    fHistPVESDHFpulls->GetAxis(iAxis)->SetTitle(PVpullstit[iAxis].Data());
    fHistPVESDHFremDaures->GetAxis(iAxis)->SetTitle(PVrestit[iAxis].Data());
    fHistPVESDHFremDaupulls->GetAxis(iAxis)->SetTitle(PVpullstit[iAxis].Data());
   }
  fHistPVESDchiS = new TH1F("fHistPVESDchiS","ESD PV - #chi^{2}",nbins,minchi,maxchi);
  fHistPVESDchiS->GetXaxis()->SetTitle("#chi^{2}");
  fHistPVESDprob = new TH1F("fHistPVESDprob","ESD PV - probability",nbins,0,1);
  fHistPVESDprob->GetXaxis()->SetTitle("prob");
  fHistPVESDHFchiS = new TH1F("fHistPVESDHFchiS","ESD PV - #chi^{2}",nbins,minchi,maxchi);
  fHistPVESDHFchiS->GetXaxis()->SetTitle("#chi^{2}");
  fHistPVESDHFprob = new TH1F("fHistPVESDHFprob","ESD PV - probability",nbins,0,1);
  fHistPVESDHFprob->GetXaxis()->SetTitle("prob");
  fHistPVESDHFremDauchiS = new TH1F("fHistPVESDHFremDauchiS","ESD PV dau removed - #chi^{2}",nbins,minchi,maxchi);
  fHistPVESDHFremDauchiS->GetXaxis()->SetTitle("#chi^{2}");
  fHistPVESDHFremDauprob = new TH1F("fHistPVESDHFremDauprob","ESD PV dau removed - probability",nbins,0,1);
  fHistPVESDHFremDauprob->GetXaxis()->SetTitle("prob");
  
  fOutput->Add(fHistPVres);
  fOutput->Add(fHistPVpulls); 
  fOutput->Add(fHistPVchiS);
  fOutput->Add(fHistPVprob);
  fOutput->Add(fHistPVHFres);
  fOutput->Add(fHistPVHFpulls); 
  fOutput->Add(fHistPVHFchiS);
  fOutput->Add(fHistPVHFprob);
  fOutput->Add(fHistPVHFremDaures);
  fOutput->Add(fHistPVHFremDaupulls);
  fOutput->Add(fHistPVHFremDauchiS);
  fOutput->Add(fHistPVHFremDauprob);
  fOutput->Add(fHistPVESDres);
  fOutput->Add(fHistPVESDpulls); 
  fOutput->Add(fHistPVESDchiS);
  fOutput->Add(fHistPVESDprob);
  fOutput->Add(fHistPVESDHFres);
  fOutput->Add(fHistPVESDHFpulls); 
  fOutput->Add(fHistPVESDHFchiS);
  fOutput->Add(fHistPVESDHFprob);
  fOutput->Add(fHistPVESDHFremDaures);
  fOutput->Add(fHistPVESDHFremDaupulls);
  fOutput->Add(fHistPVESDHFremDauchiS);
  fOutput->Add(fHistPVESDHFremDauprob);
 
  Int_t nptbins = 80;
  Double_t minpt = 0.;
  Double_t maxpt = 40.;

  Double_t minsparse[fNVarsForResHistos] = {minres,minres,minres,minresp,minresp,minresp,minresE,minresE,minresp,minres/2,minres/2,minres/2,minres,minres,minres/2,minpt};
  Double_t maxsparse[fNVarsForResHistos] = {maxres,maxres,maxres,maxresp,maxresp,maxresp,maxresE,maxresE,maxresp,maxres/2,maxres/2,maxres/2,maxres,maxres,maxres/2,maxpt};
  Double_t minsparsetopo[fNVarsForResHistos] = {minres,minres,minres,minresp,minresp,minresp,minresE,minresE,minresp,minrestopo,minrestopo,minres/2,minres,minres,minrestopo,minpt};
  Double_t maxsparsetopo[fNVarsForResHistos] = {maxres,maxres,maxres,maxresp,maxresp,maxresp,maxresE,maxresE,maxresp,maxrestopo,maxrestopo,maxres/2,maxres,maxres,maxrestopo,maxpt};
  Double_t minsparsemass[fNVarsForResHistos] = {minres,minres,minres,minresp/2,minresp/2,minresp/2,minresE,minresE/2,minres,minres/2,minres/2,minres/2,minres,minres,minres/2,minpt};
  Double_t maxsparsemass[fNVarsForResHistos] = {maxres,maxres,maxres,maxresp/2,maxresp/2,maxresp/2,maxresE,maxresE/2,maxres,maxres/2,maxres/2,maxres/2,maxres,maxres,maxres/2,maxpt};
  Double_t minsparsepull[fNVarsForResHistos] = {minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpt};
  Double_t maxsparsepull[fNVarsForResHistos] = {maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpt};
  Int_t nbinsparse[fNVarsForResHistos] = {nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nptbins};
  
  TString axtit[fNVarsForResHistos] = {"decay vtx X_{reco}-X_{MC} (#mum)",
                                       "decay vtx Y_{reco}-Y_{MC} (#mum)",
                                       "decay vtx Z_{reco}-Z_{MC} (#mum)",
                                       "(p_{X}^{reco}-p_{X}^{MC})/p_{X}^{reco}",
                                       "(p_{Y}^{reco}-p_{Y}^{MC})/p_{Y}^{reco}",
                                       "(p_{Z}^{reco}-p_{Z}^{MC})/p_{Z}^{reco}",
                                       "E_{reco}-E_{MC} (GeV)",
                                       "M_{reco}-M_{MC} (GeV/c^{2})",
                                       "(#it{p}_{T}^{reco}-#it{p}_{T}^{MC})/#it{p}_{T}^{reco}",
                                       "prod vtx X_{reco}-X_{MC} (#mum)",
                                       "prod vtx Y_{reco}-Y_{MC} (#mum)",
                                       "prod vtx Z_{reco}-Z_{MC} (#mum)",
                                       "decL_{reco} - decL_{MC}) (#mum)",
                                       "decL_{XY}^{reco} - decL_{XY}^{MC} (#mum)",
                                       "d0_{XY}^{reco}-d0_{XY}^{MC} (#mum)",
                                       "#it{p}_{T} (GeV/c)"};
  
  TString axtitpull[fNVarsForResHistos] = {"decay vtx (X_{reco}-X_{MC})/#sigma_{X_{reco}}",
                                           "decay vtx (Y_{reco}-Y_{MC})/#sigma_{Y_{reco}}",
                                           "decay vtx (Z_{reco}-Z_{MC})/#sigma_{Z_{reco}}",
                                           "(p_{X}^{reco}-p_{X}^{MC})/#sigma_{p_{X}^{reco}}",
                                           "(p_{Y}^{reco}-p_{Y}^{MC})/#sigma_{p_{Y}^{reco}}",
                                           "(p_{Z}^{reco}-p_{Z}^{MC})/#sigma_{p_{Z}^{reco}}",
                                           "(E_{reco}-E_{reco})/#sigma_{E_{MC}}",
                                           "(M_{reco}-M_{reco})/#sigma_{M_{MC}}",
                                           "(#it{p}_{T}^{reco}-#it{p}_{T}^{MC})/#sigma_{#it{p}_{T}^{reco}}",
                                           "prod vtx (X_{reco}-X_{MC})/#sigma_{X_{reco}}",
                                           "prod vtx (Y_{reco}-Y_{MC})/#sigma_{Y_{reco}}",
                                           "prod vtx (Z_{reco}-Z_{MC})/#sigma_{Z_{reco}}",
                                           "(decL_{reco} - decL_{MC})/#sigma_{decL_{reco}}",
                                           "(decL_{XY}^{reco} - decL_{XY}^{MC})/#sigma_{decL_{XY}^{reco}}",
                                           "(d0_{XY}^{reco}-d0_{XY}^{MC})/#sigma_{d0_{XY}^{reco}}",
                                           "#it{p}_{T} (GeV/c)"};
  
  fHistMres = new THnSparseF("fHistMres","Mother - res",fNVarsForResHistos,nbinsparse,minsparse,maxsparse);
  fHistMpulls = new THnSparseF("fHistMpulls","Mother - pull",fNVarsForResHistos,nbinsparse,minsparsepull,maxsparsepull);
  fHistMresTopo = new THnSparseF("fHistMresTopo","Mother - res (topo const)",fNVarsForResHistos,nbinsparse,minsparsetopo,maxsparsetopo);
  fHistMpullsTopo = new THnSparseF("fHistMpullsTopo","Mother - pull (topo const)",fNVarsForResHistos,nbinsparse,minsparsepull,maxsparsepull);
  fHistMresMass = new THnSparseF("fHistMresMass","Mother - res (mass const)",fNVarsForResHistos,nbinsparse,minsparsemass,maxsparsemass);
  fHistMpullsMass = new THnSparseF("fHistMpullsMass","Mother - pull (mass const)",fNVarsForResHistos,nbinsparse,minsparsepull,maxsparsepull);
  for(Int_t iAxis=0; iAxis<fNVarsForResHistos; iAxis++) {
    fHistMres->GetAxis(iAxis)->SetTitle(axtit[iAxis].Data());
    fHistMpulls->GetAxis(iAxis)->SetTitle(axtitpull[iAxis].Data());
    fHistMresTopo->GetAxis(iAxis)->SetTitle(axtit[iAxis].Data());
    fHistMpullsTopo->GetAxis(iAxis)->SetTitle(axtitpull[iAxis].Data());
    fHistMresMass->GetAxis(iAxis)->SetTitle(axtit[iAxis].Data());
    fHistMpullsMass->GetAxis(iAxis)->SetTitle(axtitpull[iAxis].Data());
  }
  fHistMchiS = new TH2F("fHistMchiS","Mother - #chi^{2}",nbins,minchi,maxchi,nptbins,minpt,maxpt);
  fHistMchiSTopo = new TH2F("fHistMchiSTopo","Mother - #chi^{2} (topo const)",nbins,minchi,maxchi,nptbins,minpt,maxpt);
  fHistMchiSMass = new TH2F("fHistMchiSMass","Mother - #chi^{2} (mass const)",nbins,minchi,maxchi,nptbins,minpt,maxpt);
  fHistMchiS->GetXaxis()->SetTitle("#chi^{2}");
  fHistMchiS->GetYaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  fHistMchiSTopo->GetXaxis()->SetTitle("#chi^{2}");
  fHistMchiSTopo->GetYaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  fHistMchiSMass->GetXaxis()->SetTitle("#chi^{2}");
  fHistMchiSMass->GetYaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  fHistMprob = new TH2F("fHistMprob","Mother - probability",nbins,0,1,nptbins,minpt,maxpt);
  fHistMprobTopo = new TH2F("fHistMprobTopo","Mother - probability (topo const)",nbins,0,1,nptbins,minpt,maxpt);
  fHistMprobMass = new TH2F("fHistMprobMass","Mother - probability (mass const)",nbins,0,1,nptbins,minpt,maxpt);
  fHistMprob->GetXaxis()->SetTitle("prob");
  fHistMprob->GetYaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  fHistMprobTopo->GetXaxis()->SetTitle("prob");
  fHistMprobTopo->GetYaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  fHistMprobMass->GetXaxis()->SetTitle("prob");
  fHistMprobMass->GetYaxis()->SetTitle("#it{p}_{T} (GeV/c)");

  fOutput->Add(fHistMres);
  fOutput->Add(fHistMpulls);
  fOutput->Add(fHistMresTopo);
  fOutput->Add(fHistMpullsTopo);
  fOutput->Add(fHistMresMass);
  fOutput->Add(fHistMpullsMass);
  fOutput->Add(fHistMchiS);
  fOutput->Add(fHistMprob);
  fOutput->Add(fHistMchiSTopo);
  fOutput->Add(fHistMprobTopo);
  fOutput->Add(fHistMchiSMass);
  fOutput->Add(fHistMprobMass);

  Double_t minsparseAV[13] = {minres,minres,minres,minresp,minresp,minresp,minresE,minresE,minresp,minres,minres,minres/2,minpt};
  Double_t maxsparseAV[13] = {maxres,maxres,maxres,maxresp,maxresp,maxresp,maxresE,maxresE,maxresp,maxres,maxres,maxres/2,maxpt};
  Int_t nbinsparseAV[13] = {nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nptbins};
  Double_t minsparsepullAV[9] = {minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpt};
  Double_t maxsparsepullAV[9] = {maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpt};
  Int_t nbinsparsepullAV[9] = {nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nptbins};

  fHistAliVertMres = new THnSparseF("fHistAliVertMres","AliVert - res",13,nbinsparseAV,minsparseAV,maxsparseAV);
  fHistAliVertMpulls = new THnSparseF("fHistAliVertMpulls","AliVert - pull",9,nbinsparsepullAV,minsparsepullAV,maxsparsepullAV);
  for(Int_t iAxis=0; iAxis<9; iAxis++) {
    fHistAliVertMres->GetAxis(iAxis)->SetTitle(axtit[iAxis].Data());
  }
  fHistAliVertMres->GetAxis(9)->SetTitle(axtit[fNVarsForResHistos-4].Data());
  fHistAliVertMres->GetAxis(10)->SetTitle(axtit[fNVarsForResHistos-3].Data());
  fHistAliVertMres->GetAxis(11)->SetTitle(axtit[fNVarsForResHistos-2].Data());
  fHistAliVertMres->GetAxis(12)->SetTitle(axtit[fNVarsForResHistos-1].Data());

  for(Int_t iAxis=0; iAxis<6; iAxis++) {
    fHistAliVertMpulls->GetAxis(iAxis)->SetTitle(axtitpull[iAxis].Data());
  }
  fHistAliVertMpulls->GetAxis(6)->SetTitle(axtit[fNVarsForResHistos-4].Data());
  fHistAliVertMpulls->GetAxis(7)->SetTitle(axtitpull[fNVarsForResHistos-3].Data());
  fHistAliVertMpulls->GetAxis(8)->SetTitle(axtitpull[fNVarsForResHistos-1].Data());

  fOutput->Add(fHistAliVertMres);
  fOutput->Add(fHistAliVertMpulls);

  Double_t minsparsedau[fNVarsForResHistosDau] = {minres,minres,minres,minresp,minresp,minresp,minresE,minresE/2,minres,minpt};
  Double_t maxsparsedau[fNVarsForResHistosDau] = {maxres,maxres,maxres,maxresp,maxresp,maxresp,maxresE,maxrestopoE/2,maxres,maxpt};
  Double_t minsparsepulldau[fNVarsForResHistosDau] = {minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpull,minpt};
  Double_t maxsparsepulldau[fNVarsForResHistosDau] = {maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpull,maxpt};
  Int_t nbinsparsedau[fNVarsForResHistosDau] = {nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nbins,nptbins};

  fHistDres = new THnSparseF("fHistDres","Daughters - res",fNVarsForResHistosDau,nbinsparsedau,minsparsedau,maxsparsedau);
  fHistDpulls = new THnSparseF("fHistDpulls","Daughters - pull",fNVarsForResHistosDau,nbinsparsedau,minsparsepulldau,maxsparsepulldau);
  for(Int_t iAxis=0; iAxis<fNVarsForResHistosDau-1; iAxis++) {
    fHistDres->GetAxis(iAxis)->SetTitle(axtit[iAxis].Data());
    fHistDpulls->GetAxis(iAxis)->SetTitle(axtitpull[iAxis].Data());
  }
  fHistDres->GetAxis(fNVarsForResHistosDau-1)->SetTitle(axtit[fNVarsForResHistos-1].Data());
  fHistDpulls->GetAxis(fNVarsForResHistosDau-1)->SetTitle(axtitpull[fNVarsForResHistos-1].Data());
 
  fHistDchiS = new TH2F("fHistDchiS","Daughters - #chi^{2}",nbins,minchi,maxchi,nptbins,minpt,maxpt);
  fHistDchiS->GetXaxis()->SetTitle("#chi^{2}");
  fHistDchiS->GetYaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  fHistDprob = new TH2F("fHistDprob","Daughters - probability",nbins,0,1,nptbins,minpt,maxpt);
  fHistDprob->GetXaxis()->SetTitle("prob");
  fHistDprob->GetYaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  
  fOutput->Add(fHistDres);
  fOutput->Add(fHistDpulls);
  fOutput->Add(fHistDchiS);
  fOutput->Add(fHistDprob);

  Int_t ntrbins[6] = {nbins,nbins,nbins,nbins,nbins,nptbins};
  Double_t minresTr[6] = {minres,minres,minresp/2,minresp/2,minresp,minpt};
  Double_t maxresTr[6] = {maxres,maxres,maxresp/2,maxresp/2,maxresp,maxpt};
  Double_t minpullsTr[6] = {minpull,minpull,minpull,minpull,minpull,minpt};
  Double_t maxpullsTr[6] = {maxpull,maxpull,maxpull,maxpull,maxpull,maxpt};
  TString axtitresTr[6] = {"Y_{reco}-Y_{MC} (#mum)",
                           "Z_{reco}-Z_{MC} (#mum)",
                           "sin(#phi_{reco})-sin(#phi_{MC})",
                           "#frac{p_{Z}}{p_{T}}_{reco}-#frac{p_{Z}}{p_{T}}_{MC}",
                           "#frac{q}{p_{T}}_{reco}-#frac{q}{p_{T}}_{MC} (c/GeV)",
                           "#it{p}_{T} (GeV/c)"};
  TString axtitpullsTr[6] = {"(Y_{reco}-Y_{MC})/#sigma_{Y_{reco}}",
                             "(Z_{reco}-Z_{MC})/#sigma_{Z_{reco}}",
                             "(sin(#phi_{reco})-sin(#phi_{MC}))/#sigma_{sin(#phi_{reco})}",
                             "(#frac{p_{Z}}{p_{T}}_{reco}-#frac{p_{Z}}{p_{T}}_{MC})/#sigma_{#frac{p_{Z}}{p_{T}}_{reco}}",
                             "(#frac{q}{p_{T}}_{reco}-#frac{q}{p_{T}}_{MC})/#sigma_{#frac{q}{p_{T}}_{reco}})","#it{p}_{T} (GeV/c)"};

  fHistTrackres = new THnSparseF("fHistTrackres","ESD track params - res",6,ntrbins,minresTr,maxresTr);
  fHistTrackpulls = new THnSparseF("fHistTrackpulls","ESD track params - pulls",6,ntrbins,minpullsTr,maxpullsTr);
  for(Int_t iAxis=0; iAxis<6; iAxis++) {
    fHistTrackres->GetAxis(iAxis)->SetTitle(axtitresTr[iAxis].Data());
    fHistTrackpulls->GetAxis(iAxis)->SetTitle(axtitpullsTr[iAxis].Data());
  }

  fOutput->Add(fHistTrackres);
  fOutput->Add(fHistTrackpulls);
   
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::DoQAanalysis(AliESDEvent *esdEvent, AliStack *stack)
{
  ///Number of particles in the stack
  const Int_t npart = stack->GetNtrack();

  ///Get magnetic field for KF
  Double_t Bz = esdEvent->GetMagneticField();
  KFParticle::SetField(Bz);
  
  ///check PV quality
  Double_t resPV[4] = {0., 0., 0., 0.};
  Double_t pullPV[4] = {0., 0., 0., 0.};
  Double_t resESDPV[4] = {0., 0., 0., 0.};
  Double_t pullESDPV[4] = {0., 0., 0., 0.};
  Double_t chiSESDPV;
  Int_t ndfESDPV;
  Double_t chiSPV;
  Int_t ndfPV;

  //convert primary vertex
  AliESDVertex *ESDprimary = (AliESDVertex*)esdEvent->GetPrimaryVertex(); 
  if(!ESDprimary) return;
  RecoESDPVQuality(resESDPV,pullESDPV,chiSESDPV,ndfESDPV,ESDprimary,stack);

  KFVertex pvKF;
  KFPVertex pvertex;
  ConvertVertex(ESDprimary, pvertex);
  pvKF = KFVertex(pvertex);
  RecoPVQuality(resPV,pullPV,chiSPV,ndfPV,pvKF,stack);

  ///fill histos for primary vertex QA 
  fHistPVres->Fill(resPV);
  fHistPVpulls->Fill(pullPV);
  fHistPVchiS->Fill(chiSPV/ndfPV);
  fHistPVprob->Fill(TMath::Prob(chiSPV,ndfPV));
  fHistPVESDres->Fill(resESDPV);
  fHistPVESDpulls->Fill(pullESDPV);
  fHistPVESDchiS->Fill(chiSESDPV/ndfESDPV);
  fHistPVESDprob->Fill(TMath::Prob(chiSESDPV,ndfESDPV));
 
  Int_t nDmesons=0;
  Double_t resESDnoDau[4] = {0., 0., 0., 0.};
  Double_t pullESDnoDau[4] = {0., 0., 0., 0.};
  Double_t chiSESDnoDau;
  Int_t ndfESDnoDau;
  for(Int_t iPart=0; iPart<npart; iPart++) { 
    /// Get Particle from Stack
    TParticle * mcpart = (TParticle*)stack->Particle(iPart);
    if (!mcpart)  {
      printf("Stack loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iPart);
      continue;
    }
    
    Int_t pdg = mcpart->GetPdgCode();
    if(TMath::Abs(pdg)!=fPDGcode){continue;}
    
    Int_t daulabel[4];
    Double_t decaytype;
    
    ///checks if are the corrected daughters for the selected decay channels: D0->Kpi, D+->Kpipi
    if(fPDGcode==421) {
      decaytype=AliVertexingHFUtils::CheckD0Decay(stack,iPart,daulabel);
      if(decaytype!=1) {continue;}
    }
    else if(fPDGcode==411) {
      decaytype=AliVertexingHFUtils::CheckDplusDecay(stack,iPart,daulabel);
      if(decaytype<0) {continue;}
    }
    else {
      cerr << "ERROR: Only D+ and D0 implemented!" << endl; return;
    }
      
    /// Create Daughter Particles
    TParticle *dau1 = (TParticle*)stack->Particle(daulabel[0]);
    TParticle *dau2 = (TParticle*)stack->Particle(daulabel[1]);
    TParticle *dau3 = 0x0;
    
    if(fNDau==3) dau3 = (TParticle*)stack->Particle(daulabel[2]);
    
    if(!dau1){continue;}
    if(!dau2){continue;}
    if(fPDGcode==411 && !dau3){continue;}
        
    Int_t daupdg[3];
    daupdg[0] = dau1->GetPdgCode();
    daupdg[1] = dau2->GetPdgCode();
    if(fPDGcode==411)
      daupdg[2] = dau3->GetPdgCode();

    Int_t orig = AliVertexingHFUtils::CheckOrigin(stack,mcpart,fSearchUpToQuark); //4->prompt, 5->feed-down   
    
    if(orig!=4) {continue;} //takes only prompt D mesons 

    const Int_t ntracks = esdEvent->GetNumberOfTracks();
    /// flags to check for the daughters ESD-tracks and save the tracks in new pointers
    Bool_t track[3] = {kFALSE, kFALSE, kFALSE};
    AliESDtrack* dautrack[3];
    
    for(Int_t iTrack=0; iTrack<ntracks; iTrack++) {
      /// pointer to reconstructed to track
      AliESDtrack* esdtrack = (AliESDtrack*)esdEvent->GetTrack(iTrack);
      if(!esdtrack) {
        AliError(Form("ERROR: Could not retrieve esdtrack %d",iTrack));
        continue;
      }
      
      /// if this is not a primary track, skip to the next one
      if(!fESDtrackCuts->AcceptTrack(esdtrack)) {continue;}
      
      /// read the label
      Int_t lab=esdtrack->GetLabel();
      
      if(lab==-1.) {continue;}
      
      /// compare labels of particle daughters and label of tracks => daughter gets its track
      for(Int_t iDau=0; iDau<fNDau; iDau++) {
        if(lab == daulabel[iDau]) {
          track[iDau]=kTRUE;
          dautrack[iDau]=esdtrack;
        }
      }     
      /// if all the daughters have a track -> leave loop
      if (fPDGcode == 411 && track[0] == kTRUE && track[1] == kTRUE && track[2] == kTRUE) {break;}
      if (fPDGcode == 421 && track[0] == kTRUE && track[1] == kTRUE) {break;}
    } /// loop over tracks
    
    /// only go on if every daughter has a track
    if (fPDGcode==411 && (track[0] == kFALSE || track[1] == kFALSE || track[2] == kFALSE)) {continue;}
    if (fPDGcode==421 && (track[0] == kFALSE || track[1] == kFALSE)) {continue;}

    nDmesons++;
    
    ///mother QA  
    Double_t vecforMsparseres[fNVarsForResHistos];
    Double_t vecforMsparsepulls[fNVarsForResHistos];
    Double_t chiSM;
    Int_t ndfM;
    Double_t vecforMsparserestopo[fNVarsForResHistos];
    Double_t vecforMsparsepullstopo[fNVarsForResHistos];
    Double_t chiSMtopo;
    Int_t ndfMtopo;
    Double_t vecforMsparseresmass[fNVarsForResHistos];
    Double_t vecforMsparsepullsmass[fNVarsForResHistos];
    Double_t chiSMmass;
    Int_t ndfMmass;
    Double_t vecforMsparseresAV[11];
    Double_t vecforMsparsepullsAV[7];
    float mcdecvtx[3] = {dau1->Vx(),dau1->Vy(),dau1->Vz()};    

    float decvtx[3] = {0};
    float decvtxC[6] = {0};
    ///no constraints
    RecoKFMotherQuality(vecforMsparseres,vecforMsparsepulls,chiSM,ndfM,mcpart,dautrack,daupdg,pvKF,mcdecvtx,kFALSE,kFALSE,decvtx,decvtxC); 
    ///topological constraint
    RecoKFMotherQuality(vecforMsparserestopo,vecforMsparsepullstopo,chiSMtopo,ndfMtopo,mcpart,dautrack,daupdg,pvKF,mcdecvtx,kFALSE,kTRUE,decvtx,decvtxC); 
    ///mass constraint
    RecoKFMotherQuality(vecforMsparseresmass,vecforMsparsepullsmass,chiSMmass,ndfMmass,mcpart,dautrack,daupdg,pvKF,mcdecvtx,kTRUE,kFALSE,decvtx,decvtxC);
    ///AliVertexer
    RecoESDMotherQuality(vecforMsparseresAV,vecforMsparsepullsAV,mcpart,dautrack,daupdg,ESDprimary,mcdecvtx,esdEvent);
  
    Double_t probM = TMath::Prob(chiSM, ndfM);
    if(ndfM==0) ndfM=1;
    Double_t probMtopo = TMath::Prob(chiSMtopo, ndfMtopo);
    if(ndfMtopo==0) ndfMtopo=1;
    Double_t probMmass = TMath::Prob(chiSMmass, ndfMmass);
    if(ndfMmass==0) ndfMmass=1;
    
    Double_t PtM = vecforMsparseres[fNVarsForResHistos-1];
    
    ///fill histos for mother QA
    fHistMres->Fill(vecforMsparseres);
    fHistMpulls->Fill(vecforMsparsepulls);
    fHistMresTopo->Fill(vecforMsparserestopo);
    fHistMpullsTopo->Fill(vecforMsparsepullstopo);
    fHistMresMass->Fill(vecforMsparseresmass);
    fHistMpullsMass->Fill(vecforMsparsepullsmass);
    fHistMchiS->Fill(chiSM/ndfM,PtM);
    fHistMprob->Fill(probM,PtM);
    fHistMchiSTopo->Fill(chiSMtopo/ndfMtopo,PtM);
    fHistMprobTopo->Fill(probMtopo,PtM);
    fHistMchiSMass->Fill(chiSMmass/ndfMmass,PtM);
    fHistMprobMass->Fill(probMmass,PtM);
    fHistAliVertMres->Fill(vecforMsparseresAV);
    fHistAliVertMpulls->Fill(vecforMsparsepullsAV);

    ///Vertex manipulation
    KFParticle DMesonKF;
    KFPTrack pdaughterfirst;
    KFPTrack pdaughtersecond;
    ConvertAliExternalTrackParamToKFPTrack(dautrack[0],pdaughterfirst);
    ConvertAliExternalTrackParamToKFPTrack(dautrack[1],pdaughtersecond);
    KFParticle daughterfirst(pdaughterfirst,daupdg[0]);
    KFParticle daughtersecond(pdaughtersecond,daupdg[1]);
    KFParticle daughterthird;
    if(fNDau==2) {
      DMesonKF = KFParticle(daughterfirst,daughtersecond);
    }
    else {
      KFPTrack pdaughterthird;
      ConvertAliExternalTrackParamToKFPTrack(dautrack[2],pdaughterthird);
      daughterthird = KFParticle(pdaughterthird,daupdg[2]);
      DMesonKF = KFParticle(daughterfirst,daughtersecond,daughterthird);
    }

    if(nDmesons==1) {
      fHistPVHFres->Fill(resPV);
      fHistPVHFpulls->Fill(pullPV);
      fHistPVHFchiS->Fill(chiSPV/ndfPV);
      fHistPVHFprob->Fill(TMath::Prob(chiSPV,ndfPV));
      
      fHistPVESDHFres->Fill(resESDPV);
      fHistPVESDHFpulls->Fill(pullESDPV);
      fHistPVESDHFchiS->Fill(chiSESDPV/ndfESDPV);
      fHistPVESDHFprob->Fill(TMath::Prob(chiSESDPV,ndfESDPV));
    }

    ///Remove the daughters from the vertex     
    KFVertex pvcopy = pvKF;                                                                                               
    if(ESDprimary->UsesTrack(dautrack[0]->GetID())) pvcopy -= daughterfirst;
    if(ESDprimary->UsesTrack(dautrack[1]->GetID())) pvcopy -= daughtersecond;
    if(fNDau==3 && ESDprimary->UsesTrack(dautrack[2]->GetID())) pvcopy -= daughterthird;
    ///set primary vertex constraint
    DMesonKF.SetProductionVertex(pvcopy);
    /// Add mother to the primary                                                                                                               
    pvKF += DMesonKF;
    ///remove daughters from ESD

    //check PV quality after removing the D daughters and adding the mother (only for KF) to the PV
    AliESDVertex* ESDprimaryNoDau= (AliESDVertex*)RemoveDaughtersFromESDPrimary(ESDprimary, dautrack, esdEvent);
    RecoESDPVQuality(resESDnoDau,pullESDnoDau,chiSESDnoDau,ndfESDnoDau,ESDprimaryNoDau,stack);    
    RecoPVQuality(resPV,pullPV,chiSPV,ndfPV,pvKF,stack);
    
    fHistPVHFremDaures->Fill(resPV);
    fHistPVHFremDaupulls->Fill(pullPV);
    fHistPVHFremDauchiS->Fill(chiSPV/ndfPV);
    fHistPVHFremDauprob->Fill(TMath::Prob(chiSPV,ndfPV));

    fHistPVESDHFremDaures->Fill(resESDnoDau);
    fHistPVESDHFremDaupulls->Fill(pullESDnoDau);
    fHistPVESDHFremDauchiS->Fill(chiSESDnoDau/ndfESDnoDau);
    fHistPVESDHFremDauprob->Fill(TMath::Prob(chiSESDnoDau,ndfESDnoDau));

    ///daughters QA
    TParticle* mcTrack[3];
    mcTrack[0] = dau1;
    mcTrack[1] = dau2;
    if(fNDau==3)
      mcTrack[2] = dau3;
    
    for(int iDau=0; iDau<fNDau; iDau++) {
      Double_t vecforDsparseres[fNVarsForResHistosDau];
      Double_t vecforDsparsepulls[fNVarsForResHistosDau];
      
      Double_t chiSD;
      Int_t ndfD;        
      Double_t probD;
      
      RecoKFDaughterQuality(vecforDsparseres,vecforDsparsepulls,chiSD,ndfD,mcTrack[iDau],dautrack[iDau],decvtx,decvtxC);
      
      probD = TMath::Prob(chiSD, ndfD);
      if(ndfD == 0) ndfD = 1; 

      vecforDsparseres[fNVarsForResHistosDau-1] = PtM;
      vecforDsparsepulls[fNVarsForResHistosDau-1] = PtM;
      
      ///fill histos for daughters QA      
      fHistDres->Fill(vecforDsparseres);
      fHistDpulls->Fill(vecforDsparsepulls);
      fHistDchiS->Fill(chiSD/ndfD,PtM);
      fHistDprob->Fill(probD,PtM);
     
      ///ESD tracks QA
      ///fill histos for esd tracks QA
      Double_t vecforTrsparseres[6];
      Double_t vecforTrsparsepulls[6];
      RecoESDTracksQuality(vecforTrsparseres,vecforTrsparsepulls,mcTrack[iDau],dautrack[iDau],esdEvent,ESDprimary);
      vecforTrsparseres[5]=PtM;
      vecforTrsparsepulls[5]=PtM;

      fHistTrackres->Fill(vecforTrsparseres);
      fHistTrackpulls->Fill(vecforTrsparsepulls);
    }  
  } /// loop over the stack
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::ConvertAliExternalTrackParamToKFPTrack(const AliExternalTrackParam* track, KFPTrack &kfptrack)
{
  Double_t xyztrack[3], pxpypztrack[3], covtrack[21] = {0.};
  track->GetXYZ(xyztrack);
  track->PxPyPz(pxpypztrack);
  float Charge = track->Charge();

  kfptrack.SetParameters((Float_t) xyztrack[0],(Float_t) xyztrack[1],(Float_t) xyztrack[2],(Float_t) pxpypztrack[0],(Float_t) pxpypztrack[1],(Float_t) pxpypztrack[2]);
  kfptrack.SetCharge(Charge);
  
  Double_t pt=1./TMath::Abs(track->GetParameter()[4]) * TMath::Abs(Charge);
  Double_t cs=TMath::Cos(track->GetAlpha()), sn=TMath::Sin(track->GetAlpha());
  Double_t r=TMath::Sqrt((1.-track->GetParameter()[2])*(1.+track->GetParameter()[2]));

  Double_t m00=-sn, m10=cs;
  Double_t m23=-pt*(sn + track->GetParameter()[2]*cs/r), m43=-pt*pt*(r*cs - track->GetParameter()[2]*sn);
  Double_t m24= pt*(cs - track->GetParameter()[2]*sn/r), m44=-pt*pt*(r*sn + track->GetParameter()[2]*cs);
  Double_t m35=pt, m45=-pt*pt*track->GetParameter()[3];

  m43*=track->GetSign();
  m44*=track->GetSign();
  m45*=track->GetSign();

  const Double_t *cTr = track->GetCovariance();

  Float_t covtrackf[21] = {0.f};
  covtrackf[0] = cTr[0]*m00*m00;
  covtrackf[1] = cTr[0]*m00*m10; 
  covtrackf[2] = cTr[0]*m10*m10;
  covtrackf[3] = cTr[1]*m00; 
  covtrackf[4] = cTr[1]*m10; 
  covtrackf[5] = cTr[2];
  covtrackf[6] = m00*(cTr[3]*m23 + cTr[10]*m43); 
  covtrackf[7] = m10*(cTr[3]*m23 + cTr[10]*m43); 
  covtrackf[8] = cTr[4]*m23 + cTr[11]*m43; 
  covtrackf[9] = m23*(cTr[5]*m23 + cTr[12]*m43)  +  m43*(cTr[12]*m23 + cTr[14]*m43);
  covtrackf[10] = m00*(cTr[3]*m24 + cTr[10]*m44); 
  covtrackf[11] = m10*(cTr[3]*m24 + cTr[10]*m44); 
  covtrackf[12] = cTr[4]*m24 + cTr[11]*m44; 
  covtrackf[13] = m23*(cTr[5]*m24 + cTr[12]*m44)  +  m43*(cTr[12]*m24 + cTr[14]*m44);
  covtrackf[14] = m24*(cTr[5]*m24 + cTr[12]*m44)  +  m44*(cTr[12]*m24 + cTr[14]*m44);
  covtrackf[15] = m00*(cTr[6]*m35 + cTr[10]*m45); 
  covtrackf[16] = m10*(cTr[6]*m35 + cTr[10]*m45); 
  covtrackf[17] = cTr[7]*m35 + cTr[11]*m45; 
  covtrackf[18] = m23*(cTr[8]*m35 + cTr[12]*m45)  +  m43*(cTr[13]*m35 + cTr[14]*m45);
  covtrackf[19] = m24*(cTr[8]*m35 + cTr[12]*m45)  +  m44*(cTr[13]*m35 + cTr[14]*m45); 
  covtrackf[20] = m35*(cTr[9]*m35 + cTr[13]*m45)  +  m45*(cTr[13]*m35 + cTr[14]*m45);
  
  kfptrack.SetCovarianceMatrix(covtrackf);
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::ConvertVertex(const AliVVertex* alipvertex, KFPVertex &pvertex)
{
  // Convertx double to float
  Double_t covmatrvtx[6];
  alipvertex->GetCovarianceMatrix(covmatrvtx); 
  
  Float_t covmatrvtxf[6] = {(Float_t) covmatrvtx[0],(Float_t) covmatrvtx[1],(Float_t) covmatrvtx[2],(Float_t) covmatrvtx[3],(Float_t) covmatrvtx[4],(Float_t) covmatrvtx[5]}; 
  pvertex.SetCovarianceMatrix(covmatrvtxf); 
  
  Double_t vtxpos[3]; 
  alipvertex->GetXYZ(vtxpos);
  Float_t vtxposf[3] = {(Float_t) vtxpos[0],(Float_t) vtxpos[1],(Float_t) vtxpos[2]};
  pvertex.SetXYZ(vtxposf);
  
  pvertex.SetChi2((Float_t)alipvertex->GetChi2());
  pvertex.SetNDF(alipvertex->GetNDF());
  pvertex.SetNContributors(alipvertex->GetNContributors());
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::RecoKFMotherQuality(Double_t *resvec,Double_t *pullsvec,Double_t &chiS,Int_t &ndf,TParticle *mcpart,AliESDtrack **dau,Int_t *daupdg,KFVertex &pv, float *mcdecvtx, Bool_t massconstraint,Bool_t topoconstraint, float *decvtx, float *decvtxC) 
{
  KFParticle DMesonKF;   
  float Mass, ErrMass;
  float mcE = TMath::Sqrt(fMassMean*fMassMean + mcpart->Px()*mcpart->Px() + mcpart->Py()*mcpart->Py()+ mcpart->Pz()*mcpart->Pz());
  float mcParam[fNVarsForResHistos-1]={mcpart->Vx(),mcpart->Vy(),mcpart->Vz(),mcpart->Px(),mcpart->Py(),mcpart->Pz(),mcE,fMassMean,mcpart->Pt()};     
  float recParam[fNVarsForResHistos-1] = {0};
  float errParam[fNVarsForResHistos-1] = {0};
  Double_t PtM, ErrPtM;

  KFPTrack pdaughterfirst;
  KFPTrack pdaughtersecond;
  ConvertAliExternalTrackParamToKFPTrack(dau[0],pdaughterfirst);
  ConvertAliExternalTrackParamToKFPTrack(dau[1],pdaughtersecond);
  KFParticle daughterfirst(pdaughterfirst,daupdg[0]);
  KFParticle daughtersecond(pdaughtersecond,daupdg[1]);
  KFParticle daughterthird;
  if(fNDau==2) {
    DMesonKF = KFParticle(daughterfirst,daughtersecond);
  }
  else {
    KFPTrack pdaughterthird;
    ConvertAliExternalTrackParamToKFPTrack(dau[2],pdaughterthird);
    daughterthird = KFParticle(pdaughterthird,daupdg[2]);
    DMesonKF = KFParticle(daughterfirst,daughtersecond,daughterthird);
  }

  ///primary vertex and primary vertex covariance
  float pvtx[3] = {pv.X(),pv.Y(),pv.Z()};
  float pvtxC[6] = {pv.GetCovariance(0,0),
		    pv.GetCovariance(1,1),
		    pv.GetCovariance(2,2),
		    pv.GetCovariance(3,3),
		    pv.GetCovariance(4,4),
		    pv.GetCovariance(5,5)};

  ///decay vertex and decay vertex covariance
  decvtx[0] = DMesonKF.X();
  decvtx[1] = DMesonKF.Y();
  decvtx[2] = DMesonKF.Z();
  decvtxC[0] = DMesonKF.GetCovariance(0,0);
  decvtxC[1] = DMesonKF.GetCovariance(1,1);
  decvtxC[2] = DMesonKF.GetCovariance(2,2);
  decvtxC[3] = DMesonKF.GetCovariance(3,3);
  decvtxC[4] = DMesonKF.GetCovariance(4,4);
  decvtxC[5] = DMesonKF.GetCovariance(5,5);

  float mcprodvtx[3] = {mcpart->Vx(),mcpart->Vy(),mcpart->Vz()};

  if(topoconstraint) {
    //Constraint at the primary vertex (D meson at primary vertex)
    DMesonKF.SetProductionVertex(pv);
  }    
  else {
    //transport to primary vertex to check the resolution and pulls at the primary vertex
    if(massconstraint)
      DMesonKF.SetMassConstraint(fMassMean);
    if(!fTransportToReco) 
      DMesonKF.TransportToPoint(mcprodvtx); ///MC primary vertex
  }
  
  if(topoconstraint || !fTransportToReco) {
    for(int iPar=0; iPar<7; iPar++) {
      recParam[iPar] = DMesonKF.GetParameter(iPar);
      Double_t error= DMesonKF.GetCovariance(iPar,iPar);
      if(error < 0.) {error = 1.e20;}
      errParam[iPar] = TMath::Sqrt(error);
    }    
  }
  else { ///reco decay vertex -> error on reconstructed point taken into account with GetParametersAtPoint(pvtx,pvtxC,m,V)
    float m[8];
    float V[36];
    for(Int_t iPar=0; iPar<8; iPar++)
      m[iPar]=DMesonKF.GetParameter(iPar);
    for(Int_t iC=0; iC<36; iC++) 
      V[iC] = DMesonKF.GetCovariance(iC);
    DMesonKF.GetParametersAtPoint(pvtx,pvtxC,m,V);
    for(Int_t iPar=0; iPar<7; iPar++) {
      recParam[iPar]=m[iPar];
    }
    errParam[0] = TMath::Sqrt(V[0]);
    errParam[1] = TMath::Sqrt(V[2]);
    errParam[2] = TMath::Sqrt(V[5]);
    errParam[3] = TMath::Sqrt(V[9]);
    errParam[4] = TMath::Sqrt(V[14]);
    errParam[5] = TMath::Sqrt(V[20]);
    errParam[6] = TMath::Sqrt(V[27]);
  }
  
  for(Int_t iPar=3; iPar<7; iPar++){
    if(iPar<6) { 
      resvec[iPar] = (recParam[iPar]-mcParam[iPar])/recParam[iPar]; 
      pullsvec[iPar] = (recParam[iPar]-mcParam[iPar])/errParam[iPar]; 
    }
    else {
      resvec[iPar] = recParam[iPar]-mcParam[iPar]; 
      pullsvec[iPar] = resvec[iPar]/errParam[iPar]; 
    }
  }
  
  DMesonKF.GetMass(Mass,ErrMass);
  PtM = DMesonKF.GetPt();
  ErrPtM=DMesonKF.GetErrPt();

  resvec[7] = Mass - mcParam[7];
  if(TMath::Abs(ErrMass) > 1.e-20) pullsvec[7] = resvec[7]/ErrMass;    
  resvec[8] = (PtM - mcParam[8])/PtM;
  if(TMath::Abs(ErrPtM) > 1.e-20) pullsvec[8] = (PtM - mcParam[8])/ErrPtM;    
  //put the residuals and pulls at PV at the end of the array
  resvec[9] = (recParam[0] - mcParam[0])*10000; //microns
  if(TMath::Abs(errParam[0]) > 1.e-20) pullsvec[9] = resvec[9]/(errParam[0]*10000);    
  resvec[10] = (recParam[1] - mcParam[1])*10000; //microns
  if(TMath::Abs(errParam[1]) > 1.e-20) pullsvec[10] = resvec[10]/(errParam[1]*10000);    
  resvec[11] = (recParam[2] - mcParam[2])*10000; //microns
  if(TMath::Abs(errParam[2]) > 1.e-20) pullsvec[11] = resvec[11]/(errParam[2]*10000);    

  float IPxy;
  float IPxyerr;
  DMesonKF.GetDistanceFromVertexXY(pv,IPxy,IPxyerr);
  resvec[14] = IPxy;
  pullsvec[14] = IPxy/IPxyerr;

  //if topological constraint transport to decay vertex in order to take the residuals and pulls at decay vertex
  if(topoconstraint) {
    if(!fTransportToReco) {
      DMesonKF.TransportToPoint(mcdecvtx); ///MC decay vertex   
      for(Int_t iPar=0; iPar<3; iPar++) {
        decvtx[iPar] = DMesonKF.GetParameter(iPar);
        decvtxC[iPar] = DMesonKF.GetCovariance(iPar,iPar);
        if(decvtxC[iPar] < 0.) {decvtxC[iPar] = 1.e20;}
        errParam[iPar] = TMath::Sqrt(decvtxC[iPar]);
      }
    } 
    else { ///reco decay vertex -> error on reconstructed point taken into account with GetParametersAtPoint(decvtx,decvtxC,m,V)
      float m[8];
      float V[36];
      for(Int_t iPar=0; iPar<8; iPar++)
        m[iPar]=DMesonKF.GetParameter(iPar);
      for(Int_t iC=0; iC<36; iC++) 
        V[iC] = DMesonKF.GetCovariance(iC);
      DMesonKF.GetParametersAtPoint(decvtx,decvtxC,m,V);
      errParam[0] = TMath::Sqrt(V[0]);
      errParam[1] = TMath::Sqrt(V[2]);
      errParam[2] = TMath::Sqrt(V[5]);
    }
    DMesonKF.TransportToPoint(decvtx); ///reco decay vertex   
    for(Int_t iPar=0; iPar<3; iPar++) {
      decvtx[iPar]=DMesonKF.GetParameter(iPar);
    }
  }
  else {///if not topo constraint
    for(Int_t iPar=0; iPar<3; iPar++)
      errParam[iPar] = TMath::Sqrt(decvtxC[iPar]);
  }

  ///store decay vertex residuals and pulls
  for(Int_t iPar=0; iPar<3; iPar++) {
    resvec[iPar] = (decvtx[iPar]-mcdecvtx[iPar])*10000; //microns 
    pullsvec[iPar] = resvec[iPar]/(errParam[iPar]*10000); 
  }

  //decay length
  float decL;
  float errdecL;
  DMesonKF.GetDecayLength(decL,errdecL);
  Double_t mcdecL = TMath::Sqrt((mcdecvtx[0]-mcpart->Vx())*(mcdecvtx[0]-mcpart->Vx())+
                                (mcdecvtx[1]-mcpart->Vy())*(mcdecvtx[1]-mcpart->Vy())+
                                (mcdecvtx[2]-mcpart->Vz())*(mcdecvtx[2]-mcpart->Vz()));
  resvec[12] = (decL-mcdecL)*10000; //microns
  if(errdecL<0) {errdecL=1.e20;}
  pullsvec[12] = resvec[12]/(errdecL*10000);

  //decay length XY
  float decLXY;
  float errdecLXY;
  DMesonKF.GetDecayLengthXY(decLXY,errdecLXY);
  Double_t mcdecLXY = TMath::Sqrt((mcdecvtx[0]-mcpart->Vx())*(mcdecvtx[0]-mcpart->Vx())+
                                  (mcdecvtx[1]-mcpart->Vy())*(mcdecvtx[1]-mcpart->Vy()));

  resvec[13] = (decLXY-mcdecLXY)*10000; //microns
  if(errdecLXY<0) {errdecLXY=1.e20;}
  pullsvec[13] = resvec[13]/(errdecLXY*10000);

  resvec[fNVarsForResHistos-1] = PtM;
  pullsvec[fNVarsForResHistos-1] = PtM;

  chiS = DMesonKF.GetChi2();
  ndf = DMesonKF.GetNDF();
}  

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::RecoKFDaughterQuality(Double_t *resvec,Double_t *pullsvec,Double_t &chiS,Int_t &ndf,TParticle *mctrack,AliESDtrack *esdtrack, float *recodecvtx, float *recodecvtxC) 
{
  TParticlePDG* particlePDG = TDatabasePDG::Instance()->GetParticle(mctrack->GetPdgCode());
  float mcVX =  mctrack->Vx();
  float mcVY =  mctrack->Vy();
  float mcVZ =  mctrack->Vz();
  float mcPx = mctrack->Px();
  float mcPy = mctrack->Py();
  float mcPz = mctrack->Pz();
  float mcPt = mctrack->Pt();
      
  float decayVtx[3] = {mcVX, mcVY, mcVZ};
  float recParam[fNVarsForResHistosDau] = {0};
  float errParam[fNVarsForResHistosDau] = {0};

  Double_t massMC = (particlePDG) ? particlePDG->Mass() : -1.; 
  Double_t Emc = TMath::Sqrt(mcPx*mcPx+ mcPy*mcPy + mcPz*mcPz + massMC*massMC);
  float mcParam[fNVarsForResHistosDau] = {mcVX,mcVY,mcVZ,mcPx,mcPy,mcPz,Emc,massMC,mcPt};
    
  float Mass, ErrMass;
  float Pt, ErrPt;

  KFPTrack pdaughter;
  ConvertAliExternalTrackParamToKFPTrack(esdtrack,pdaughter);
  KFParticle DaughterKF(pdaughter,mctrack->GetPdgCode()); 
  if(fTransportToReco) {
    float m[8];
    float V[36];
    for(Int_t iPar=0; iPar<8; iPar++)
      m[iPar]=DaughterKF.GetParameter(iPar);
    for(Int_t iC=0; iC<36; iC++) 
      V[iC] = DaughterKF.GetCovariance(iC);
    DaughterKF.GetParametersAtPoint(recodecvtx,recodecvtxC,m,V);
    for(Int_t iPar=0; iPar<7; iPar++) {
      recParam[iPar]=m[iPar];
    }
    errParam[0] = TMath::Sqrt(V[0]);
    errParam[1] = TMath::Sqrt(V[2]);
    errParam[2] = TMath::Sqrt(V[5]);
    errParam[3] = TMath::Sqrt(V[9]);
    errParam[4] = TMath::Sqrt(V[14]);
    errParam[5] = TMath::Sqrt(V[20]);
    errParam[6] = TMath::Sqrt(V[27]);
  }
  else {
    DaughterKF.TransportToPoint(decayVtx); ///MC decay vertex
    for(int iPar=0; iPar < 7; iPar++ ) {
      Double_t error = DaughterKF.GetCovariance(iPar,iPar);
      if(error < 0.) { error = 1.e20;}
      recParam[iPar] = DaughterKF.GetParameter(iPar);
      errParam[iPar] = TMath::Sqrt(error);
    }
  }

  for(Int_t iPar=0; iPar<7; iPar++) {
    if(iPar<3) {
      resvec[iPar]  = (recParam[iPar] - mcParam[iPar])*10000;
      if(errParam[iPar] > 1.e-20) pullsvec[iPar] = resvec[iPar]/(errParam[iPar]*10000);
    }
    else if(iPar>=3 && iPar<6) {
      resvec[iPar]  = (recParam[iPar] - mcParam[iPar])/recParam[iPar];
      if(errParam[iPar] > 1.e-20) pullsvec[iPar] = (recParam[iPar] - mcParam[iPar])/errParam[iPar];
    }
    else {
      resvec[iPar]  = recParam[iPar] - mcParam[iPar];
      if(errParam[iPar] > 1.e-20) pullsvec[iPar] = resvec[iPar]/errParam[iPar];
    }
  }

  DaughterKF.GetMass(Mass,ErrMass);
  Pt = DaughterKF.GetPt();
  ErrPt = DaughterKF.GetErrPt();
  chiS = DaughterKF.GetChi2();
  ndf = DaughterKF.GetNDF(); 
  
  resvec[7] = Mass - mcParam[7];
  if(TMath::Abs(ErrMass) > 1.e-20) pullsvec[7] = resvec[7]/ErrMass;

  resvec[8] = (Pt - mcParam[8])/Pt;
  if(TMath::Abs(ErrMass) > 1.e-20) pullsvec[8] = (Pt - mcParam[8])/ErrPt;
}
  
//________________________________________________________________________
void AliAnalysisTaskDmesonKF::RecoPVQuality(Double_t *resvec,Double_t *pullsvec,Double_t &chiS,Int_t &ndf,KFVertex pv,AliStack *stack)
{
  Int_t npart = stack->GetNtrack();
  Double_t mcPV[3] = {0., 0., 0.};
  for(Int_t iPartMC=0; iPartMC<npart; iPartMC++) {
    TParticle *mcpart = (TParticle*)stack->Particle(iPartMC);
    if(mcpart->GetFirstMother()<0) {//it means that is a first mother particle
      mcPV[0] = mcpart->Vx();
      mcPV[1] = mcpart->Vy();
      mcPV[2] = mcpart->Vz();
      break;
    }
  }
  Int_t nTracks = pv.GetNContributors();
  resvec[0] = (pv.X()-mcPV[0])*10000; //microns
  resvec[1] = (pv.Y()-mcPV[1])*10000; //microns
  resvec[2] = (pv.Z()-mcPV[2])*10000; //microns
  pullsvec[0] = resvec[0]/(TMath::Sqrt(pv.GetCovariance(0,0))*10000);
  pullsvec[1] = resvec[1]/(TMath::Sqrt(pv.GetCovariance(1,1))*10000);
  pullsvec[2] = resvec[2]/(TMath::Sqrt(pv.GetCovariance(2,2))*10000);
  chiS = pv.Chi2();
  ndf = pv.NDF();
  
  resvec[3] = nTracks;
  pullsvec[3] = nTracks;
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::RecoESDMotherQuality(Double_t *resvec,Double_t *pullsvec,/* Double_t &chiS, Int_t &ndf,*/TParticle *mcpart,AliESDtrack** dau, Int_t *daupdg, AliESDVertex* pv, float *mcdecvtx, AliESDEvent* esdEvent) 
{
  ///secondary vertex reconstruction 
  Double_t Bz = esdEvent->GetMagneticField();
  AliVertexerTracks* vertexer = new AliVertexerTracks(Bz);
  vertexer->SetVtxStart(pv);

  TObjArray trkarray(fNDau);
  for(Int_t iDau=0; iDau<fNDau; iDau++) {
    trkarray.AddLast(dau[iDau]);
  }

  AliESDVertex* secvtx = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(&trkarray);  
  trkarray.Clear();
  resvec[0] = (secvtx->GetX()-mcdecvtx[0])*10000; //microns
  resvec[1] = (secvtx->GetY()-mcdecvtx[1])*10000; //microns
  resvec[2] = (secvtx->GetZ()-mcdecvtx[2])*10000; //microns
  pullsvec[0] = resvec[0]/(secvtx->GetXRes()*10000);
  pullsvec[1] = resvec[1]/(secvtx->GetYRes()*10000);
  pullsvec[2] = resvec[2]/(secvtx->GetZRes()*10000);
  
  ///momentum, invariant mass, energy 
  Double_t px = 0.;
  Double_t py = 0.;
  Double_t pz = 0.;
  Double_t massdau[3];
  Double_t energydau[3];
  Double_t energy = 0.;

  Double_t cov[21];
  Double_t d0z0[2];
  Double_t covd0z0[3];

  for(Int_t iDau=0; iDau<fNDau; iDau++) {
    dau[iDau]->PropagateToDCA(pv,Bz,100000,d0z0,covd0z0);
    px += dau[iDau]->Px();
    py += dau[iDau]->Py();
    pz += dau[iDau]->Pz(); 
    massdau[iDau] = TDatabasePDG::Instance()->GetParticle(daupdg[iDau])->Mass();
    energydau[iDau] = TMath::Sqrt(massdau[iDau]*massdau[iDau]+dau[iDau]->P()*dau[iDau]->P());
    energy += energydau[iDau];

    Double_t daucov[21];
    dau[iDau]->GetCovarianceXYZPxPyPz(daucov);
    for(Int_t iCov=0; iCov<21; iCov++) {
      cov[iCov] += daucov[iCov];
    }
  }
  Double_t pt = TMath::Sqrt(px*px+py*py);
  Double_t p = TMath::Sqrt(px*px+py*py+pz*pz);
  Double_t mass = TMath::Sqrt(energy*energy-p*p);

  resvec[3] = (px-mcpart->Px())/px;
  resvec[4] = (py-mcpart->Py())/py;
  resvec[5] = (pz-mcpart->Pz())/pz;
  pullsvec[3] = (px-mcpart->Px())/TMath::Sqrt(cov[9]);
  pullsvec[4] = (py-mcpart->Py())/TMath::Sqrt(cov[14]);
  pullsvec[5] = (pz-mcpart->Pz())/TMath::Sqrt(cov[20]);

  ///only resolution
  resvec[6] = energy-mcpart->Energy();
  resvec[7] = mass-fMassMean;
  resvec[8] = (pt-mcpart->Pt())/pt; 

  ///decay length 
  Double_t decL = TMath::Sqrt((secvtx->GetX()-pv->GetX())*(secvtx->GetX()-pv->GetX())+
                              (secvtx->GetY()-pv->GetY())*(secvtx->GetY()-pv->GetY())+
                              (secvtx->GetZ()-pv->GetZ())*(secvtx->GetZ()-pv->GetZ()));

  Double_t Pchi2perNDF = pv->GetChi2toNDF();
  Double_t Ppos[3];
  Double_t Pcov[6];
  pv->GetXYZ(Ppos);
  pv->GetCovMatrix(Pcov);

  Double_t Schi2perNDF = secvtx->GetChi2toNDF();
  Double_t Spos[3];
  Double_t Scov[6];
  secvtx->GetXYZ(Spos);
  secvtx->GetCovMatrix(Scov);

  AliAODVertex* AODprimary = new AliAODVertex(Ppos,Pcov,Pchi2perNDF);
  AliAODVertex* AODsecondary = new AliAODVertex(Spos,Scov,Schi2perNDF);

  Double_t errdecL = TMath::Sqrt(AODsecondary->Error2DistanceToVertex(AODprimary));
 
  Double_t mcdecL = TMath::Sqrt((mcdecvtx[0]-mcpart->Vx())*(mcdecvtx[0]-mcpart->Vx())+
                                (mcdecvtx[1]-mcpart->Vy())*(mcdecvtx[1]-mcpart->Vy())+
                                (mcdecvtx[2]-mcpart->Vz())*(mcdecvtx[2]-mcpart->Vz()));
  resvec[9] = (decL - mcdecL)*10000; //microns 
  pullsvec[6] = resvec[9]/(errdecL*10000);

  ///decay lengthXY 
  Double_t decLXY = TMath::Sqrt((secvtx->GetX()-pv->GetX())*(secvtx->GetX()-pv->GetX())+
                                (secvtx->GetY()-pv->GetY())*(secvtx->GetY()-pv->GetY()));
  Double_t mcdecLXY = TMath::Sqrt((mcdecvtx[0]-mcpart->Vx())*(mcdecvtx[0]-mcpart->Vx())+
                                  (mcdecvtx[1]-mcpart->Vy())*(mcdecvtx[1]-mcpart->Vy()));
  Double_t errdecLXY = TMath::Sqrt(AODsecondary->Error2DistanceXYToVertex(AODprimary));
  
  resvec[10] = (decLXY - mcdecLXY)*10000; //microns 
  pullsvec[7] = resvec[10]/(errdecLXY*10000);
				
  ///imp par XY
  Double_t k = -(secvtx->GetX()-pv->GetX())*px-(secvtx->GetY()-pv->GetY())*py;
  k /= pt*pt;
  Double_t dx = secvtx->GetX()-pv->GetX()+k*px;
  Double_t dy = secvtx->GetY()-pv->GetY()+k*py;
  Double_t absImpPar = TMath::Sqrt(dx*dx+dy*dy);
  TVector3 mom(px,py,pz);
  TVector3 fline(secvtx->GetX()-pv->GetX(),secvtx->GetY()-pv->GetY(),secvtx->GetZ()-pv->GetZ());
  TVector3 cross = mom.Cross(fline);

  if(cross.Z()>0.) resvec[11] = absImpPar*10000; //microns
  else resvec[11] = -absImpPar*10000; //microns  

  resvec[12] = pt;
  pullsvec[8] = pt;

  delete vertexer;
  delete AODprimary;
  delete AODsecondary;
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::RecoESDTracksQuality(Double_t *resvec, Double_t *pullsvec, TParticle* mctrack, AliESDtrack* track, AliESDEvent* esdEvent, AliESDVertex* pv)
{
  Double_t cosA = TMath::Cos(track->GetAlpha());
  Double_t sinA = TMath::Sin(track->GetAlpha());
  //Double_t mcX = mctrack->Vx()*cosA + mctrack->Vy()*sinA;
  Double_t mcY = -mctrack->Vx()*sinA + mctrack->Vy()*cosA;
  Double_t mcZ = mctrack->Vz();
  Double_t mcEx = mctrack->Px()*cosA + mctrack->Py()*sinA;
  Double_t mcEy = -mctrack->Px()*sinA + mctrack->Py()*cosA;
  Double_t mcEz = mctrack->Pz();
  Double_t mcEt = TMath::Sqrt(mcEx*mcEx + mcEy*mcEy);
  Double_t mcSinPhi = mcEy/mcEt;
  Double_t mcDzDs = mcEz/mcEt;
  Double_t mcQ=0;
  
  if(mctrack->GetPdgCode() < 9999999)
    mcQ = TDatabasePDG::Instance()->GetParticle(mctrack->GetPdgCode())->Charge()/3.0;
  Double_t mcQPt = mcQ/mctrack->Pt();
  
  Double_t d0z0[2];
  Double_t covd0z0[3];
  track->PropagateToDCA(pv,esdEvent->GetMagneticField(),100000,d0z0,covd0z0);

  Double_t MCPar[5] = {mcY, mcZ, mcSinPhi, mcDzDs, mcQPt};
  Double_t Err[5] = {track->GetSigmaY2(),track->GetSigmaZ2(),
                     track->GetSigmaSnp2(),track->GetSigmaTgl2(),track->GetSigma1Pt2()};
  
  for(int iPar=0; iPar<5; iPar++) {
    if(iPar<2) {
      resvec[iPar] = (track->GetParameter()[iPar]-MCPar[iPar])*10000; //microns
      if(Err[iPar] >= 0.)
        pullsvec[iPar] = resvec[iPar]/(TMath::Sqrt(Err[iPar])*10000);
    }
    else {
      resvec[iPar] = track->GetParameter()[iPar]-MCPar[iPar];
      if(Err[iPar] >= 0.)
        pullsvec[iPar] = resvec[iPar]/TMath::Sqrt(Err[iPar]);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::RecoESDPVQuality(Double_t *resvec, Double_t *pullsvec, Double_t &chiS, Int_t &ndf, AliESDVertex* pv, AliStack *stack) 
{
  Int_t npart = stack->GetNtrack();
  Double_t mcPV[3];
  for(Int_t iPartMC=0; iPartMC<npart; iPartMC++) {
    TParticle *mcpart = (TParticle*)stack->Particle(iPartMC);
    if(mcpart->GetFirstMother()<0) {//it means that is a first mother particle
      mcPV[0] = mcpart->Vx();
      mcPV[1] = mcpart->Vy();
      mcPV[2] = mcpart->Vz();
      break;
    }
  }
  Int_t nTracks = pv->GetNContributors();
  resvec[0] = (pv->GetX()-mcPV[0])*10000; //microns
  resvec[1] = (pv->GetY()-mcPV[1])*10000; //microns
  resvec[2] = (pv->GetZ()-mcPV[2])*10000; //microns
  Double_t covmatrix[6];
  pv->GetCovarianceMatrix(covmatrix);
  pullsvec[0] = resvec[0]/(TMath::Sqrt(covmatrix[0])*10000);
  pullsvec[1] = resvec[1]/(TMath::Sqrt(covmatrix[2])*10000);
  pullsvec[2] = resvec[2]/(TMath::Sqrt(covmatrix[5])*10000);
  chiS = pv->GetChi2();
  ndf = pv->GetNDF();
  
  resvec[3] = nTracks;
  pullsvec[3] = nTracks;
}

//________________________________________________________________________
AliESDVertex* AliAnalysisTaskDmesonKF::RemoveDaughtersFromESDPrimary(AliESDVertex *pv, AliESDtrack** dautracks, AliESDEvent* esdEvent) {

  AliVertexerTracks *vertexer = new AliVertexerTracks(esdEvent->GetMagneticField());
  vertexer->SetITSMode();
  vertexer->SetMinClusters(3);
  vertexer->SetConstraintOff();
  vertexer->SetVtxStart(pv);

  Int_t id[3];
  Int_t ntracks = fNDau;
  for(Int_t iDau=0; iDau<fNDau; iDau++) {
    if(dautracks[iDau]->GetID()) ntracks--; 
    else id[iDau] = dautracks[iDau]->GetID();
  }

  vertexer->SetSkipTracks(ntracks,id);
  AliESDVertex *ESDprimaryDauRem = vertexer->FindPrimaryVertex(esdEvent);
  delete vertexer;
  if(!ESDprimaryDauRem) return 0;
  if(ESDprimaryDauRem->GetNContributors()<=0) {
    delete ESDprimaryDauRem;
    ESDprimaryDauRem = NULL;
    return 0;
  }

  return ESDprimaryDauRem;
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::GetTracksArrays(AliESDEvent *esdEvent) 
{  
  const Int_t nTracks = esdEvent->GetNumberOfTracks();

  for(Int_t iTrack=0; iTrack<nTracks; iTrack++) {
    AliESDtrack *track = (AliESDtrack*)esdEvent->GetTrack(iTrack);
    if(!fESDtrackCuts->AcceptTrack(track)) {continue;}
    
    Double_t charge=track->GetSign();
    if(charge>0)
      fPosTracksArray->AddLast(track);
    else
      fNegTracksArray->AddLast(track);
  }
}

//________________________________________________________________________
Int_t AliAnalysisTaskDmesonKF::IdentifyParticleTPC(AliESDtrack* track, AliPID::EParticleType species)
{
  if(!CheckTPCPIDStatus(track)) return 0;

  Double_t nTPCsigmas= TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,species));
  Double_t pt = track->Pt();
  
  Bool_t TPCcompatible=kFALSE;
  Bool_t TPCidentified=kFALSE;

  //TPC
  if(nTPCsigmas>3) {
    TPCcompatible = kFALSE;
    TPCidentified = kFALSE;
  }
  else if(nTPCsigmas>2 && nTPCsigmas<3) {
    TPCcompatible=kTRUE;
  } 
  else if(nTPCsigmas>1 && nTPCsigmas<2) {
    if(pt<0.6) TPCidentified=kTRUE;
    else TPCcompatible=kTRUE;
  }
  else {
    if(pt<0.8) TPCidentified=kTRUE;
    else TPCcompatible=kTRUE;
  } 
  
  if(TPCidentified)
    return 1;
  else if(TPCcompatible)
    return 0;
  else
    return -1;
}

//________________________________________________________________________
Int_t AliAnalysisTaskDmesonKF::IdentifyParticleTOF(AliESDtrack* track, AliPID::EParticleType species)
{
  if(!CheckTOFPIDStatus(track)) return 0;

  Double_t nTOFsigmas= TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,species));
  Double_t pt = track->Pt();
  
  Bool_t TOFcompatible=kFALSE;
  Bool_t TOFidentified=kFALSE;
  
  //TOF
  if(nTOFsigmas>3) {
    TOFcompatible = kFALSE; 
    TOFidentified = kFALSE;
  }
  else {
    if(pt<1.5) TOFidentified = kTRUE;
    else TOFcompatible = kTRUE;
  }
  
  if(TOFidentified)
    return 1;
  else if(TOFcompatible)
    return 0;
  else
    return -1;
}

//________________________________________________________________________
Int_t AliAnalysisTaskDmesonKF::DplusPID(AliESDtrack** tracks, Bool_t strongPID) 
{
  ///0->rejected, 1->D+
  ///first 2 tracks always like sign, third different sign
  Int_t protonCombined[3] = {0,0,0};
  Int_t kaonCombined[3] = {0,0,0};
  Int_t pionCombined[3] = {0,0,0};

  for(Int_t iTrack=0; iTrack<fNDau; iTrack++) {
    protonCombined[iTrack] = IdentifyParticleTPC(tracks[iTrack],AliPID::kProton) + IdentifyParticleTOF(tracks[iTrack],AliPID::kProton);
    kaonCombined[iTrack] = IdentifyParticleTPC(tracks[iTrack],AliPID::kKaon) + IdentifyParticleTOF(tracks[iTrack],AliPID::kKaon);
    pionCombined[iTrack] = IdentifyParticleTPC(tracks[iTrack],AliPID::kPion) + IdentifyParticleTOF(tracks[iTrack],AliPID::kPion);
  }

  ///convention: if TPC and TOF disagree (rejected Vs identified) -> unknown (compatible)
 
  ///if at least one track is identified as a proton and not as a pion or kaon reject
  Int_t protID=0;
  for(Int_t iTrack=0; iTrack<fNDau; iTrack++) {
    if(protonCombined[iTrack]>=1 && pionCombined[iTrack]<=-1 && kaonCombined[iTrack]<=-1)
      protID++;
  }
  if(protID>0) return 0;
  
  ///if more than one track is identified as a kaon and rejected as a pion, from at least TOF or TPC, reject
  Int_t accKaonRejPion=0;
  for(Int_t iTrack=0; iTrack<fNDau; iTrack++) {
    if(kaonCombined[iTrack]==2 && pionCombined[iTrack]<=-1)
      accKaonRejPion++;
  }
  if(accKaonRejPion>1) return 0;
  
  ///if more than 2 tracks are rejected as kaons from TPC or TOF reject
  Int_t rejKaon=0;
  for(Int_t iTrack=0; iTrack<fNDau; iTrack++) {
    if(kaonCombined[iTrack]<=-1)
      rejKaon++;
  }
  if(rejKaon>2) return 0;
  
  ///if the same sign tracks are not compatible with the pion hypothesis reject
  if(pionCombined[0]<=-1 || pionCombined[1]<=-1) return 0;
  ///if the different sign track is not compatible with the kaon hypothesis reject
  if(kaonCombined[2]<=-1) return 0;
   
  ///if strong PID the pions must be identified
  if(strongPID && (pionCombined[0]<=0 || pionCombined[1]<=0)) return 0;
  ///if strong PID the kaon must be identified
  if(strongPID && kaonCombined[2]<=0) return 0;

  return 1;
  
}

//________________________________________________________________________
Int_t AliAnalysisTaskDmesonKF::D0PID(AliESDtrack** tracks, Bool_t strongPID) 
{
  ///0->rejected, 1->D0, 2->D0bar, 3->D0 or D0bar
  Int_t kaonCombined[2] = {0,0};
  Int_t pionCombined[2] = {0,0};
  
  for(Int_t iTrack=0; iTrack<fNDau; iTrack++) {
    kaonCombined[iTrack] = IdentifyParticleTPC(tracks[iTrack],AliPID::kKaon) + IdentifyParticleTOF(tracks[iTrack],AliPID::kKaon);
    pionCombined[iTrack] = IdentifyParticleTPC(tracks[iTrack],AliPID::kPion) + IdentifyParticleTOF(tracks[iTrack],AliPID::kPion);
  }

  Int_t isD0D0barPID[2] = {1,2};

  //convention: if TPC and TOF disagree (rejected Vs identified) -> unknown (compatible)

  for(Int_t iTrack=0; iTrack<fNDau; iTrack++) {
    if(kaonCombined[iTrack]<=-1 && pionCombined[iTrack]<=-1){// if not a K- and not a pi- both D0 and D0bar excluded
      isD0D0barPID[0]=0;
      isD0D0barPID[1]=0;
    }
    else if(kaonCombined[iTrack]==2 && pionCombined[iTrack]>=1){// if in conflict (both pi- and K-), if k for both TPC and TOF -> is K
      if(tracks[iTrack]->Charge()==-1)isD0D0barPID[1]=0;//if K- D0bar excluded
      else isD0D0barPID[0]=0;// if K+ D0 excluded
    }
    else if(kaonCombined[iTrack]==1 && pionCombined[iTrack]>=1){// if in conflict (both pi- and K-) and k- only for TPC or TOF -> reject
      isD0D0barPID[0]=0;
      isD0D0barPID[1]=0;
    }
    else if(kaonCombined[iTrack]>=1 || pionCombined[iTrack]<=-1){
      if(tracks[iTrack]->Charge()==-1)isD0D0barPID[1]=0;// not a D0bar if K- or if pi- excluded
      else isD0D0barPID[0]=0;//  not a D0 if K+ or if pi+ excluded
    }
    else if(kaonCombined[iTrack]<=-1 || pionCombined[iTrack]>=1){
      if(tracks[iTrack]->Charge()==-1)isD0D0barPID[0]=0;// not a D0 if pi- or if K- excluded
      else isD0D0barPID[1]=0;// not a D0bar if pi+ or if K+ excluded
    }

    //strong PID -> kaon must be identified
    if(strongPID && kaonCombined[iTrack]<=0) {
      if(tracks[iTrack]->Charge()==-1) isD0D0barPID[0]=0; /// if is not identified as a K- is not a D0
      else isD0D0barPID[1]=0; /// if is not identified as a K+ is not a D0bar
    }
  }
  
  return isD0D0barPID[0]+isD0D0barPID[1];
}

//________________________________________________________________________
Bool_t AliAnalysisTaskDmesonKF::CheckTPCPIDStatus(AliESDtrack* track) 
{
  AliPIDResponse::EDetPidStatus status = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,track);
  if(status != AliPIDResponse::kDetPidOk) return kFALSE;
  UInt_t nclsTPCPID = track->GetTPCsignalN();
  if(nclsTPCPID<fMinNClustersTPCPID) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskDmesonKF::CheckTOFPIDStatus(AliESDtrack* track) 
{
  AliPIDResponse::EDetPidStatus status = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,track);
  if(status != AliPIDResponse::kDetPidOk) return kFALSE;
  Float_t probMis = fPIDResponse->GetTOFMismatchProbability(track);
  if(probMis>fCutTOFmismatch) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::Combine2ProngsKF(KFVertex pv, AliStack* stack)
{  
  Int_t nPos=fPosTracksArray->GetEntriesFast();
  Int_t nNeg=fNegTracksArray->GetEntriesFast();

  //combine positive pions with negative kaons or vice versa
  for(Int_t iPos=0; iPos<nPos; iPos++) {
    for(Int_t iNeg=0; iNeg<nNeg; iNeg++) {
      AliESDtrack* postrk = (AliESDtrack*)fPosTracksArray->UncheckedAt(iPos); 
      AliESDtrack* negtrk = (AliESDtrack*)fNegTracksArray->UncheckedAt(iNeg);
      AliESDtrack* tracks[2] = {postrk,negtrk};
      Int_t label1 = postrk->GetLabel();
      Int_t label2 = negtrk->GetLabel();
      if(fReadMC && (label1<0 || label2<0)) continue;
      
      KFPTrack ppospart;
      KFPTrack pnegpart;
      KFParticle PosPart;
      KFParticle NegPart;
      KFParticle DMesonKF;
      KFParticle PosPart2;
      KFParticle NegPart2;
      KFParticle DMesonKF2;

      Int_t daupdg[2] = {0,0};
      Int_t PIDrespo = D0PID(tracks);
      
      if(PIDrespo==1) { ///D0 
        daupdg[0] = 211; ///pion pdg
        daupdg[1] = -321; ///kaon pdg
      }
      else if(PIDrespo==2) { ///D0bar
        daupdg[0] = 321; ///kaon pdg 
        daupdg[1] = -211; ///pion pdg
      }
      else if(PIDrespo==3) { ///if compatible with both D0 and D0bar or no PID
        ///start with D0
        daupdg[0] = 211; ///pion pdg
        daupdg[1] = -321; ///kaon pdg
      }
      else if(fUsePID && PIDrespo==0) { ///if rejected and use PID continue
        continue;
      }
      else { ///if rejected and don't use PID start from D0
        daupdg[0] = 211; ///pion pdg
        daupdg[1] = -321; ///kaon pdg
      }
      
      ConvertAliExternalTrackParamToKFPTrack(postrk,ppospart);
      ConvertAliExternalTrackParamToKFPTrack(negtrk,pnegpart);
      PosPart = KFParticle(ppospart,daupdg[0]); 
      NegPart = KFParticle(pnegpart,daupdg[1]);
      KFParticle* dau[2] = {&PosPart,&NegPart};

      DMesonKF = KFParticle(PosPart,NegPart); ///D meson candidate

      KFParticle* dau2[2];
      if(PIDrespo==3 || PIDrespo==0) {///if dont't use PID or is both D0 and D0bar build also D0bar
        PosPart2 = KFParticle(ppospart,-daupdg[1]); 
        NegPart2 = KFParticle(pnegpart,-daupdg[0]);
        dau2[0] = &PosPart2;
        dau2[1] = &NegPart2;
        DMesonKF2 = KFParticle(PosPart2,NegPart2); 
      }
      
      Double_t mass = DMesonKF.GetMass();
      Double_t pt = DMesonKF.GetPt();
      Double_t y = GetRapidity(DMesonKF,fPDGcode);
      Double_t mass2 = -1.;
      Double_t pt2 = -1.;
      Double_t y2 = -1.;

      if(PIDrespo==3 || PIDrespo==0) {
        mass2 = DMesonKF2.GetMass();
        pt2 = DMesonKF2.GetPt();
        y2 = GetRapidity(DMesonKF2,fPDGcode);
      }
      
      if(fUseStrongPID && (pt<fMaxPtForStrongPID || ((PIDrespo==3 || PIDrespo==0) && pt2<fMaxPtForStrongPID))) {
        //apply strong PID
        PIDrespo = D0PID(tracks,kTRUE);
        if(fUsePID && PIDrespo==0)
          continue;
      }
      
      Bool_t D0acc=kFALSE;      
      Bool_t D0baracc=kFALSE;
      
      if(PIDrespo!=3 && PIDrespo!=0){ ///if use PID and is only a D0 or D0bar 
        if(mass>fMassMax || mass<fMassMin || pt<fPtMin || !IsInFiducialAcceptance(pt, y)) continue;
      }
      else {///if don't use PID or is both D0 and D0bar
        if(mass<fMassMax && mass>fMassMin && pt>fPtMin && IsInFiducialAcceptance(pt, y)) D0acc=kTRUE;
        if(mass2<fMassMax && mass2>fMassMin && pt2>fPtMin && IsInFiducialAcceptance(pt2, y2)) D0baracc=kTRUE;
        if(!D0acc && !D0baracc) continue;
      }
      
      KFVertex pvcopy = pv;

      if((PIDrespo!=3 && PIDrespo!=0) || D0acc) {
        ///Remove the daughters and add the mother in the primary vertex
        pvcopy -= PosPart;
        pvcopy -= NegPart;	
        pvcopy += DMesonKF;
   
        Double_t daupt1=PosPart.GetPt();
        Double_t daupt2=NegPart.GetPt();
        Double_t mindaupt=-1;
        
        if(daupt1<daupt2)
          mindaupt = daupt1;
        else
          mindaupt = daupt2;

        Double_t sigvtx = GetSigmaVtx(DMesonKF,dau);
        Double_t d0d0exp = GetNormMaxd0d0exp(DMesonKF,dau,pv);

        if(PIDrespo==3 && D0acc && !D0baracc) PIDrespo=1;

        FillSparses(DMesonKF,pv,label1,label2,-1,mindaupt,sigvtx,d0d0exp,PIDrespo,stack);
      }
      if(D0baracc) {
        ///Remove the daughters and add the mother in the primary vertex
        pvcopy -= PosPart2;
        pvcopy -= NegPart2;	
        pvcopy += DMesonKF2;

        Double_t daupt1=PosPart2.GetPt();
        Double_t daupt2=NegPart2.GetPt();
        Double_t mindaupt=-1.;
        
        if(daupt1<daupt2)
          mindaupt = daupt1;
        else
          mindaupt = daupt2;

        Double_t sigvtx = GetSigmaVtx(DMesonKF2,dau2);
        Double_t d0d0exp = GetNormMaxd0d0exp(DMesonKF,dau,pv);

        if(PIDrespo==3 && !D0acc && D0baracc) PIDrespo=2;

        FillSparses(DMesonKF,pv,label1,label2,-1,mindaupt,sigvtx,d0d0exp,PIDrespo,stack);        
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::Combine3ProngsKF(KFVertex pv, AliStack* stack) 
{  
  Int_t nPos=fPosTracksArray->GetEntriesFast();
  Int_t nNeg=fNegTracksArray->GetEntriesFast();

  KFParticle DMesonKF;  

  Double_t mass=-1;
  Double_t pt=-1;
  Double_t y=-1;
  Int_t label1=-1;
  Int_t label2=-1;
  Int_t label3=-1;
  Double_t daupt1=0;
  Double_t daupt2=0;
  Double_t daupt3=0;
  Double_t mindaupt=-1;
  Double_t sigvtx=-1;
  Double_t d0d0exp=-1000;
  Int_t PIDrespo=0;

  //combine two positive pions with one negative kaon 
  for(Int_t iPos1=0; iPos1<nPos; iPos1++) {
    for(Int_t iPos2=iPos1; iPos2<nPos; iPos2++) { ///avoid double counting
      for(Int_t iNeg=0; iNeg<nNeg; iNeg++) {
        
        if(iPos1==iPos2) continue;

        AliESDtrack* postrk1 = (AliESDtrack*)fPosTracksArray->UncheckedAt(iPos1); 
        AliESDtrack* postrk2 = (AliESDtrack*)fPosTracksArray->UncheckedAt(iPos2);
        AliESDtrack* negtrk = (AliESDtrack*)fNegTracksArray->UncheckedAt(iNeg);	
        AliESDtrack* tracks[3] = {postrk1,postrk2,negtrk};
        label1=postrk1->GetLabel();
        label2=postrk2->GetLabel();
        label3=negtrk->GetLabel();
  	if(fReadMC && (label1<0 || label2<0 || label3<0)) continue;
        
        KFPTrack ppospart1;
        KFPTrack ppospart2;
        KFPTrack pnegpart;
        KFParticle PosPart1;
        KFParticle PosPart2;
        KFParticle NegPart;
        
        PIDrespo=DplusPID(tracks);  
        if(fUsePID && PIDrespo==0) continue;
        
        ConvertAliExternalTrackParamToKFPTrack(postrk1,ppospart1);
        ConvertAliExternalTrackParamToKFPTrack(postrk2,ppospart2);
        ConvertAliExternalTrackParamToKFPTrack(negtrk,pnegpart);
        PosPart1 = KFParticle(ppospart1,211); //pion pdg
        PosPart2 = KFParticle(ppospart2,211); //pion pdg
        NegPart = KFParticle(pnegpart,-321); //kaon pdg
        KFParticle* dau[3] = {&PosPart1,&PosPart2,&NegPart};
        
        DMesonKF = KFParticle(PosPart1,PosPart2,NegPart); //D+ meson candidate
    
        mass = DMesonKF.GetMass();
        pt = DMesonKF.GetPt();
        y = GetRapidity(DMesonKF,fPDGcode);

        if(mass>fMassMax || mass<fMassMin || pt<fPtMin || !IsInFiducialAcceptance(pt, y)) continue;
        if(fUseStrongPID && pt<fMaxPtForStrongPID) {
          //apply strong PID
          PIDrespo = DplusPID(tracks,kTRUE);
          if(fUsePID && PIDrespo==0)
            continue;
        }
        
        KFVertex pvcopy;
        // Remove the daughters from the vertex
        pvcopy -= PosPart1;
        pvcopy -= PosPart2;
        pvcopy -= NegPart;
        pvcopy += DMesonKF;
        
        daupt1=PosPart1.GetPt();
        daupt2=PosPart2.GetPt();
        daupt3=NegPart.GetPt();	
        if(daupt1<daupt2 && daupt1<daupt3)
          mindaupt = daupt1;
        else if(daupt2<daupt1 && daupt2<daupt3)
          mindaupt = daupt2;
        else
          mindaupt = daupt3;
        
        sigvtx = GetSigmaVtx(DMesonKF,dau);
        d0d0exp = GetNormMaxd0d0exp(DMesonKF,dau,pv);

        FillSparses(DMesonKF,pv,label1,label2,label3,mindaupt,sigvtx,d0d0exp,PIDrespo,stack);
      }
    }
  }
  //combine two negative pions with one positive kaon 
  for(Int_t iNeg1=0; iNeg1<nNeg; iNeg1++) {
    for(Int_t iNeg2=iNeg1; iNeg2<nNeg; iNeg2++) { ///avoid double counting
      for(Int_t iPos=0; iPos<nPos; iPos++) {
        
        if(iNeg1==iNeg2) continue;
		
        AliESDtrack* negtrk1 = (AliESDtrack*)fNegTracksArray->UncheckedAt(iNeg1); 
        AliESDtrack* negtrk2 = (AliESDtrack*)fNegTracksArray->UncheckedAt(iNeg2);
        AliESDtrack* postrk = (AliESDtrack*)fPosTracksArray->UncheckedAt(iPos);	
        AliESDtrack* tracks[3] = {negtrk1,negtrk2,postrk};	
        label1=negtrk1->GetLabel();
        label2=negtrk2->GetLabel();
        label3=postrk->GetLabel();
	if(fReadMC && (label1<0 || label2<0 || label3<0)) continue;
        
        KFPTrack pnegpart1;
        KFPTrack pnegpart2;
        KFPTrack ppospart;
        KFParticle NegPart1;
        KFParticle NegPart2;
        KFParticle PosPart;
        
        PIDrespo=DplusPID(tracks);  
        if(fUsePID && PIDrespo==0) continue;
        
        ConvertAliExternalTrackParamToKFPTrack(negtrk1,pnegpart1);
        ConvertAliExternalTrackParamToKFPTrack(negtrk2,pnegpart2);
        ConvertAliExternalTrackParamToKFPTrack(postrk,ppospart);
        NegPart1 = KFParticle(pnegpart1,-211); //pion pdg
        NegPart2 = KFParticle(pnegpart2,-211); //pion pdg
        PosPart = KFParticle(ppospart,321); //kaon pdg
        KFParticle* dau[3] = {&NegPart1,&NegPart2,&PosPart};
        
        DMesonKF = KFParticle(NegPart1,NegPart2,PosPart); //D- meson candidate
        
        mass = DMesonKF.GetMass();
        pt = DMesonKF.GetPt();
        y = GetRapidity(DMesonKF,fPDGcode);

        if(mass>fMassMax || mass<fMassMin || pt<fPtMin || !IsInFiducialAcceptance(pt, y)) continue;
        if(fUseStrongPID && pt<fMaxPtForStrongPID) {
          //apply strong PID
          PIDrespo = DplusPID(tracks,kTRUE);
          if(fUsePID && PIDrespo==0)
            continue;
        }
        
        KFVertex pvcopy = pv;
        pvcopy -= NegPart1;
        pvcopy -= NegPart2;
        pvcopy -= PosPart;
        pvcopy += DMesonKF;
        
        daupt1=NegPart1.GetPt();
        daupt2=NegPart2.GetPt();
        daupt3=PosPart.GetPt();
        if(daupt1<daupt2 && daupt1<daupt3)
          mindaupt = daupt1;
        else if(daupt2<daupt1 && daupt2<daupt3)
          mindaupt = daupt2;
        else
          mindaupt = daupt3;
        
        sigvtx = GetSigmaVtx(DMesonKF,dau);
        d0d0exp = GetNormMaxd0d0exp(DMesonKF,dau,pv);

        FillSparses(DMesonKF,pv,label1,label2,label3,mindaupt,sigvtx,d0d0exp,PIDrespo,stack);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::FillSparses(KFParticle DMesonKF,KFVertex pv,Int_t label1,Int_t label2,Int_t label3,Double_t mindaupt,Double_t sigvtx,Double_t d0d0exp,Int_t PID,AliStack* stack)
{
  ///fill THnsparse for all reconstructed candidates and only prompt or feed-down in case of fReadMC

  Double_t mass = DMesonKF.GetMass();
  Double_t Pt = DMesonKF.GetPt();
  Double_t d0XY = DMesonKF.GetDistanceFromVertexXY(pv)*10000;    
  Double_t cosp = TMath::Cos(GetPointingAngle(DMesonKF,pv,kFALSE));
  Double_t cospXY = TMath::Cos(GetPointingAngle(DMesonKF,pv,kTRUE));
    
  ///add mother to the primary vertex and evaluate the PV chi square
  KFVertex pvcopy = pv;
  pvcopy += DMesonKF;
  Double_t chi2PV = pvcopy.GetChi2();
  Int_t ndfPV = pvcopy.GetNDF();
  Double_t PVredchi=chi2PV/ndfPV;
  
  ///set production vertex (without mother)
  KFParticle DMesonKFcopy = DMesonKF;
  DMesonKFcopy.SetProductionVertex(pv);
  Double_t chi2 = DMesonKFcopy.GetChi2();
  Int_t ndf = DMesonKFcopy.GetNDF();
  Double_t redchi=chi2/ndf;
  float decL;
  float decLerr;
  float decLXY;
  float decLXYerr;
  DMesonKFcopy.GetDecayLength(decL,decLerr);
  DMesonKFcopy.GetDecayLengthXY(decLXY,decLXYerr);
  
  Double_t NdecLXY=decLXY/decLXYerr;

  Int_t iPt=0;
  while(Pt>fPtLims[iPt+1])
    iPt++;

  if(!fLooseCuts || (fLooseCuts && mindaupt>fMinDauPt[iPt] && cosp>fCospMin[iPt] && decL>fDecLMin[iPt] && NdecLXY>fNDecLXYMin[iPt] && sigvtx<fSigVtxMax[iPt] && redchi<fChiMax[iPt] && PVredchi<fPVChiMax[iPt])) {
  
    Double_t vecforsparse[fNVarsForSparses] = {mass,Pt,d0XY,PID,decL,decLXY,NdecLXY,cosp,cospXY,mindaupt,sigvtx,d0d0exp,redchi,PVredchi};
    ///not to go overflow if cosp is exactly 1
    if(TMath::Abs(vecforsparse[7])==1) {
      vecforsparse[8]*=0.99999;
    }
    if(TMath::Abs(vecforsparse[8])==1) {
    vecforsparse[9]*=0.99999;
    }
    ///not to go underflow if chi2/ndf is exactly 0
    if(vecforsparse[12]==0) {
      vecforsparse[12] = 0.000001;
    }
    if(vecforsparse[13]==0) {
      vecforsparse[13] = 0.000001;
    }
    fSparse[0]->Fill(vecforsparse);
  
    if(fReadMC && stack) {
      Int_t orig;
      if(fPDGcode==421)
        orig=MatchToMC(stack,label1,label2,-1);    
      else
        orig=MatchToMC(stack,label1,label2,label3); 
      if(orig==4) 
        fSparse[1]->Fill(vecforsparse);
      if(orig==5)
        fSparse[2]->Fill(vecforsparse);
    }
}
} 

//________________________________________________________________________
Int_t AliAnalysisTaskDmesonKF::MatchToMC(AliStack* stack, Int_t firstlabel,Int_t secondlabel,Int_t thirdlabel) 
{  

  Int_t labDau[3]={firstlabel,secondlabel,thirdlabel};
  Int_t daupdg[3]={211,321,-1};
  if(fPDGcode==411)
    daupdg[2] = 211;
  Int_t labMom[3]={0.,0.,0.};
  Int_t lab, labMother, pdgMother, pdgPart;
  TParticle* part = 0x0;
  TParticle* mother = 0x0;
  Double_t pxSumDgs=0., pySumDgs=0, pzSumDgs=0.;
  Bool_t pdgUsed[3]={kFALSE,kFALSE,kFALSE};

  ///loop on daughters label
  for(Int_t iDau=0; iDau<fNDau; iDau++) {
    labMom[iDau]=-1.;
    lab = TMath::Abs(labDau[iDau]);
    if(lab<0) {
      cerr << Form("daughter with negative label %d",lab) <<endl;
      return 0;
    }
    part = (TParticle*)stack->Particle(lab);
    if(!part) {
      cerr << "no MC particle" <<endl;
      return 0;
    }

    ///check the pdg of the daughters
    pdgPart=TMath::Abs(part->GetPdgCode());
    for(Int_t iDauPdg=0; iDauPdg<fNDau; iDauPdg++) {
      if(!pdgUsed[iDauPdg] && pdgPart==daupdg[iDauPdg]) {
	pdgUsed[iDauPdg]=kTRUE;
	break;
      }
    }

    //check mothers
    mother=part;
    while(mother->GetFirstMother()>=0) {
      labMother=mother->GetFirstMother();
      mother=(TParticle*)stack->Particle(labMother);
      if(!mother) {
	cerr << "no mother particle" << endl;
	break;    
      }
      pdgMother=TMath::Abs(mother->GetPdgCode());
      if(pdgMother==fPDGcode) {
	labMom[iDau]=labMother;
	pxSumDgs += part->Px();
	pySumDgs += part->Py();
	pzSumDgs += part->Pz();
	break;
      }
      else if(pdgMother>fPDGcode || pdgMother<10) {
	break;
      }
    }
    if(labMom[iDau]==-1) return 0;
  }///end loop on daughters

  ///check if the candidate is signal
  labMother=labMom[0];
  for(Int_t iDau=0; iDau<fNDau; iDau++) {
    if(labMom[iDau]==-1) return 0;
    if(labMom[iDau]!=labMother) return 0;
  }

  ///check that all daughter PDGs are matched
  for(Int_t iDau=0; iDau<fNDau; iDau++) {
    if(pdgUsed[iDau]==kFALSE) return 0;
  } 

  ///checks mom conservation
  mother = (TParticle*)stack->Particle(labMother);
  Double_t pxMother = mother->Px();
  Double_t pyMother = mother->Py();
  Double_t pzMother = mother->Pz();
  ///within 0.1%
  if((TMath::Abs(pxMother-pxSumDgs)/(TMath::Abs(pxMother)+1.e-13))>0.00001 && 
     (TMath::Abs(pyMother-pySumDgs)/(TMath::Abs(pyMother)+1.e-13))>0.00001 && 
     (TMath::Abs(pzMother-pzSumDgs)/(TMath::Abs(pzMother)+1.e-13))>0.00001)
    return 0;

  Int_t orig=AliVertexingHFUtils::CheckOrigin(stack,mother,fSearchUpToQuark); ///4->prompt 5->feed-down

  /*
  if(firstlabel<0 || secondlabel<0) return 0;
  
  TParticle *mcdau1 = (TParticle*)stack->Particle(firstlabel);
  TParticle *mcdau2 = (TParticle*)stack->Particle(secondlabel);
  TParticle *mcdau3 = 0x0;
  Int_t motherlab1 = mcdau1->GetFirstMother();
  Int_t motherlab2 = mcdau2->GetFirstMother();
  Int_t motherlab3 = -1.;  
 
  if(motherlab1<0) return 0;
  if(motherlab1!=motherlab2) return 0;
  
  if(fPDGcode==411) {
    if(thirdlabel<0) return 0;
    mcdau3 = (TParticle*)stack->Particle(thirdlabel);
    motherlab3 = mcdau3->GetFirstMother();
    if(motherlab3!=motherlab1 || motherlab3!=motherlab2) return 0;
  }

  TParticle *mother= (TParticle*)stack->Particle(motherlab1);
  Int_t pdgcode = mother->GetPdgCode();
  if(TMath::Abs(pdgcode)!=fPDGcode) return 0;  
  
  Int_t orig=AliVertexingHFUtils::CheckOrigin(stack,mother,fSearchUpToQuark); ///4->prompt 5->feed-down
  */

  return orig;
}  

//________________________________________________________________________
Bool_t AliAnalysisTaskDmesonKF::IsInFiducialAcceptance(Double_t pt, Double_t y) 
{
  if(fMaxRapidCand>-998.) {
    if(TMath::Abs(y) > fMaxRapidCand) return kFALSE;
    else return kTRUE;
  }

  if(pt>5.) {
    if(TMath::Abs(y)>0.8) return kFALSE;
  }
  else { 
    //apply smooth cut for pt<5 GeV/c
    Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5;
    Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;
    if(y<minFiducialY || y>maxFiducialY) return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskDmesonKF::FillMCGenAccHistos(AliStack* stack)
{  
  Int_t npart=stack->GetNtrack();
  Double_t mcPvZ;
  for(Int_t iPartMC=0; iPartMC<npart; iPartMC++) {
    TParticle *mcpart = (TParticle*)stack->Particle(iPartMC);
    if(mcpart->GetFirstMother()<0) {//it means that is a first mother particle
      mcPvZ = mcpart->Vz();
      break;
    } 
  }
  
  for(Int_t iPartMC=0; iPartMC<npart; iPartMC++) {
    TParticle *mcpart = (TParticle*)stack->Particle(iPartMC);
    if(TMath::Abs(mcpart->GetPdgCode())==fPDGcode) {
      Int_t orig = AliVertexingHFUtils::CheckOrigin(stack,mcpart,fSearchUpToQuark); //4->prompt, 5->feed-down   
      
      Int_t daulabel[4];
      Double_t decaytype=0;
      Bool_t isGoodDecay=kFALSE;

      ///checks if are the corrected daughters for the selected decay channels: D0->Kpi, D+->Kpipi
      if(fPDGcode==421) {
	decaytype=AliVertexingHFUtils::CheckD0Decay(stack,iPartMC,daulabel);
	if(decaytype==1) isGoodDecay=kTRUE;
      }
      else if(fPDGcode==411) {
	decaytype=AliVertexingHFUtils::CheckDplusDecay(stack,iPartMC,daulabel);
	if(decaytype>0) isGoodDecay=kTRUE;
      }
      else {
	cerr << "ERROR: Only D+ and D0 implemented!" << endl; return;
      }
      if(daulabel[0]<0) continue;
      
      if(isGoodDecay) {
	///check if the D meson is in fiducial acceptance
	Bool_t isFidAcc=IsInFiducialAcceptance(mcpart->Pt(),mcpart->Y()); 
	///check if the daughters are in acceptance
	Bool_t isInAcc=CheckAcceptance(stack,fNDau,daulabel);
	
	if(TMath::Abs(mcPvZ)<fMaxVtxZ && isFidAcc && isInAcc) {
	  Double_t sparsearray[2]={mcpart->Pt(),mcpart->Y()};
	  if(orig==4) 
	    fMCGenAccPrompt->Fill(sparsearray);
	  else if(orig==5)
	  fMCGenAccFD->Fill(sparsearray);
	  else 
	    continue;
	}
      }
    }
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskDmesonKF::CheckAcceptance(AliStack *stack, Int_t nDau, Int_t *labdau) 
{  
  for(Int_t iDau=0; iDau<nDau; iDau++) {
    TParticle* mcDau = (TParticle*)stack->Particle(labdau[iDau]);
    if(!mcDau) return kFALSE;
    Double_t eta= mcDau->Eta();
    Double_t pt= mcDau->Pt();
    if(TMath::Abs(eta)>0.9 || pt<0.1) return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskDmesonKF::GetRapidity(KFParticle part, Int_t pdg) 
{
  Double_t px = part.GetPx();
  Double_t py = part.GetPy();
  Double_t pz = part.GetPz();

  Double_t mass = fMassMean=TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  Double_t E = TMath::Sqrt(mass*mass+px*px+py*py+pz*pz);

  Double_t y = 0.5*TMath::Log((E+pz)/(E-pz+1.e-13));

  return y;
}

//________________________________________________________________________
Double_t AliAnalysisTaskDmesonKF::GetPointingAngle(KFParticle part, KFVertex pv, Bool_t isXY) 
{
  //momentum of the particle
  TVector3 momentum(part.GetPx(),part.GetPy(),part.GetPz());
  if(isXY)
    momentum.SetZ(0.);

  part.TransportToDecayVertex();
  //flight line
  TVector3 fline(part.X()-pv.X(),part.Y()-pv.Y(),part.Z()-pv.Z());
  if(isXY)
    fline.SetZ(0.);

  Double_t pointangle = momentum.Angle(fline);

  return pointangle;
}

//________________________________________________________________________
Double_t AliAnalysisTaskDmesonKF::GetSigmaVtx(KFParticle part, KFParticle *dau[])
{
  part.TransportToDecayVertex();  
  float decayvtx[3] = {part.X(),part.Y(),part.Z()};
  Double_t d[3];

  for(Int_t iDau=0; iDau<fNDau; iDau++) {
    dau[iDau]->TransportToPoint(decayvtx);
    d[iDau] = dau[iDau]->GetDistanceFromVertex(decayvtx);
  }
  
  Double_t sigvtx;
  if(fNDau==2)
    sigvtx = TMath::Sqrt(d[0]*d[0]+d[1]*d[1]);
  else
    sigvtx = TMath::Sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
  
  return sigvtx;
}

//________________________________________________________________________
Double_t AliAnalysisTaskDmesonKF::GetNormMaxd0d0exp(KFParticle part, KFParticle *dau[], KFVertex pv)
{
  float dlxy;
  float dlxyerr; 
  part.SetProductionVertex(pv);
  part.GetDecayLength(dlxy,dlxyerr);
  part.TransportToDecayVertex();
  float secvtx[3] = {part.X(),part.Y(),part.Z()};

  Double_t normdiff[3];

  for(Int_t iDau=0; iDau<fNDau; iDau++) {
    float d0meas;
    float d0measerr; 
    dau[iDau]->GetDistanceFromVertexXY(pv,d0meas,d0measerr);    
    dau[iDau]->TransportToPoint(secvtx);
    Double_t px = dau[iDau]->Px();
    Double_t py = dau[iDau]->Py();
    Double_t pt = dau[iDau]->GetPt();
    Double_t sinthetap = (px*part.Py()-py*part.Px())/(part.GetPt()*pt);
    normdiff[iDau] = (d0meas-dlxy*sinthetap)/TMath::Sqrt(d0measerr*d0measerr+dlxyerr*dlxyerr*sinthetap*sinthetap);
  }  

  if(fNDau==2) {
    if(TMath::Abs(normdiff[0])>TMath::Abs(normdiff[1])) return normdiff[0];
    else return normdiff[1];
  }
  else {
    if(TMath::Abs(normdiff[0])>TMath::Abs(normdiff[1]) && TMath::Abs(normdiff[0])>TMath::Abs(normdiff[2])) return normdiff[0];
    else if(TMath::Abs(normdiff[1])>TMath::Abs(normdiff[0]) && TMath::Abs(normdiff[1])>TMath::Abs(normdiff[2])) return normdiff[1];
    else return normdiff[2];
  }
}
