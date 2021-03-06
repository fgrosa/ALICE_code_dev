#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenerator.h"
#include "AliVertexingHFUtils.h"
#include "AliMultiplicity.h"
#include <TParticle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TChain.h>
#include "AliESDInputHandlerRP.h"
#include "AliAnalysisTaskCheckHFMCProd2.h"

/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
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
// Implementation of class AliAnalysisTaskCheckHFMCProd2
// AliAnalysisTask to check MC production at ESD+Kine level
// 
//
// Authors: F. Prino, prino@to.infn.it
//          
//*************************************************************************

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCheckHFMCProd2);
/// \endcond

//______________________________________________________________________________
AliAnalysisTaskCheckHFMCProd2::AliAnalysisTaskCheckHFMCProd2() : AliAnalysisTaskSE("HFMCChecks"), 
  fOutput(0),
  fHistoNEvents(0),
  fHistoPhysPrim(0),
  fHistoTracks(0),
  fHistoSelTracks(0),
  fHistoTracklets(0),
  fHistoTrackletsEta1(0),
  fHistoPtPhysPrim(0),
  fHistoEtaPhysPrim(0),
  fHistoSPD3DVtxX(0),
  fHistoSPD3DVtxY(0),
  fHistoSPD3DVtxZ(0),
  fHistoSPDZVtxX(0),
  fHistoSPDZVtxY(0),
  fHistoSPDZVtxZ(0),
  fHistoTRKVtxX(0),
  fHistoTRKVtxY(0),
  fHistoTRKVtxZ(0),
  fHistoNcharmed(0),
  fHistoNbVsNc(0),
  fHistOriginPrompt(0),
  fHistOriginFeeddown(0),
  fHistMotherID(0),
  fHistDSpecies(0),
  fHistBSpecies(0),
  fHistLcDecayChan(0),
  fHistNcollHFtype(0),
  fHistEtaPhiPtGenEle(0),
  fHistEtaPhiPtGenPi(0),
  fHistEtaPhiPtGenK(0),
  fHistEtaPhiPtGenPro(0),
  fHistEtaPhiPtRecEle(0),
  fHistEtaPhiPtRecPi(0),
  fHistEtaPhiPtRecK(0),
  fHistEtaPhiPtRecPro(0),
  fSearchUpToQuark(kFALSE),
  fSystem(0),
  fESDtrackCuts(0x0),
  fReadMC(kTRUE)
{
  //
  for(Int_t i=0; i<5; i++){
    fHistBYPtAllDecay[i]=0x0;
    fHistYPtAllDecay[i]=0x0;
    fHistYPtPromptAllDecay[i]=0x0;
    fHistYPtFeeddownAllDecay[i]=0x0;
    fHistYPtPrompt[i]=0x0;
    fHistYPtFeeddown[i]=0x0;
  }
  for(Int_t i=0; i<2; i++){
    fHistYPtD0byDecChannel[i]=0x0;
    fHistYPtDplusbyDecChannel[i]=0x0;
    fHistYPtDsbyDecChannel[i]=0x0;
  }
  fHistYPtDplusbyDecChannel[2]=0x0;

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskCheckHFMCProd2::~AliAnalysisTaskCheckHFMCProd2(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  delete fESDtrackCuts;

}
   
//___________________________________________________________________________
void AliAnalysisTaskCheckHFMCProd2::UserCreateOutputObjects() {
  /// create output histos

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistoNEvents = new TH1F("hNEvents", "Number of processed events",3,-0.5,2.5);
  fHistoNEvents->Sumw2();
  fHistoNEvents->SetMinimum(0);
  fOutput->Add(fHistoNEvents);

  Double_t maxMult=100.;
  if(fSystem==1) maxMult=10000.;
  if(fSystem==2) maxMult=500.;
  fHistoPhysPrim = new TH1F("hPhysPrim","",100,-0.5,maxMult-0.5);
  fHistoPhysPrim->Sumw2();
  fOutput->Add(fHistoPhysPrim);
  fHistoTracks = new TH1F("hTracks","",100,-0.5,maxMult*2-0.5);
  fHistoTracks->Sumw2();
  fOutput->Add(fHistoTracks);
  fHistoSelTracks = new TH1F("hSelTracks","",100,-0.5,maxMult-0.5);
  fHistoSelTracks->Sumw2();
  fOutput->Add(fHistoSelTracks);
  fHistoTracklets = new TH1F("hTracklets","",100,-0.5,maxMult-0.5);
  fHistoTracklets->Sumw2();
  fOutput->Add(fHistoTracklets);
  fHistoTrackletsEta1 = new TH1F("hTrackletsEta1","",100,-0.5,maxMult-0.5);
  fHistoTrackletsEta1->Sumw2();
  fOutput->Add(fHistoTrackletsEta1);
  fHistoPtPhysPrim = new TH1F("hPtPhysPrim","",100,0.,20.);
  fHistoPtPhysPrim->Sumw2();
  fOutput->Add(fHistoPtPhysPrim);
  fHistoEtaPhysPrim = new TH1F("hEtaPhysPrim","",100,-10.,10.);
  fHistoEtaPhysPrim->Sumw2();
  fOutput->Add(fHistoEtaPhysPrim);

  fHistoSPD3DVtxX = new TH1F("hSPD3DvX","",100,-1.,1.);
  fHistoSPD3DVtxX->Sumw2();
  fOutput->Add(fHistoSPD3DVtxX);
  fHistoSPD3DVtxY = new TH1F("hSPD3DvY","",100,-1.,1.);
  fHistoSPD3DVtxY->Sumw2();
  fOutput->Add(fHistoSPD3DVtxY);
  fHistoSPD3DVtxZ = new TH1F("hSPD3DvZ","",100,-15.,15.);
  fHistoSPD3DVtxZ->Sumw2();
  fOutput->Add(fHistoSPD3DVtxZ);

  fHistoSPDZVtxX = new TH1F("hSPDZvX","",100,-1.,1.);
  fHistoSPDZVtxX->Sumw2();
  fOutput->Add(fHistoSPDZVtxX);
  fHistoSPDZVtxY = new TH1F("hSPDZvY","",100,-1.,1.);
  fHistoSPDZVtxY->Sumw2();
  fOutput->Add(fHistoSPDZVtxY);
  fHistoSPDZVtxZ = new TH1F("hSPDZvZ","",100,-15.,15.);
  fHistoSPDZVtxZ->Sumw2();
  fOutput->Add(fHistoSPDZVtxZ);


  fHistoTRKVtxX = new TH1F("hTRKvX","",100,-1.,1.);
  fHistoTRKVtxX->Sumw2();
  fOutput->Add(fHistoTRKVtxX);
  fHistoTRKVtxY = new TH1F("hTRKvY","",100,-1.,1.);
  fHistoTRKVtxY->Sumw2();
  fOutput->Add(fHistoTRKVtxY);
  fHistoTRKVtxZ = new TH1F("hTRKvZ","",100,-15.,15.);
  fHistoTRKVtxZ->Sumw2();
  fOutput->Add(fHistoTRKVtxZ);

  Int_t nBinscb=11;
  if(fSystem==1) nBinscb=200;
  if(fSystem==2) nBinscb=21;
  Double_t maxncn=nBinscb-0.5;
  fHistoNcharmed = new TH2F("hncharmed","",100,-0.5,maxMult-0.5,nBinscb,-0.5,maxncn);
  fHistoNcharmed->Sumw2();
  fOutput->Add(fHistoNcharmed);
  fHistoNbVsNc = new TH2F("hnbvsnc","",nBinscb,-0.5,maxncn,nBinscb,-0.5,maxncn);
  fHistoNbVsNc->Sumw2();
  fOutput->Add(fHistoNbVsNc);

  fHistYPtPrompt[0] = new TH2F("hyptD0prompt","D0 - Prompt",80,0.,40.,20,-2.,2.);
  fHistYPtPrompt[1] = new TH2F("hyptDplusprompt","Dplus - Prompt",80,0.,40.,20,-2.,2.);
  fHistYPtPrompt[2] = new TH2F("hyptDstarprompt","Dstar - Prompt",80,0.,40.,20,-2.,2.);
  fHistYPtPrompt[3] = new TH2F("hyptDsprompt","Ds - Prompt",80,0.,40.,20,-2.,2.);
  fHistYPtPrompt[4] = new TH2F("hyptLcprompt","Lc - Prompt",80,0.,40.,20,-2.,2.);

  fHistBYPtAllDecay[0] = new TH2F("hyptB0AllDecay","B0 - All",80,0.,40.,40,-2.,2.);
  fHistBYPtAllDecay[1] = new TH2F("hyptBplusAllDecay","Bplus - All",80,0.,40.,40,-2.,2.);
  fHistBYPtAllDecay[2] = new TH2F("hyptBstarAllDecay","Bstar - All",80,0.,40.,40,-2.,2.);
  fHistBYPtAllDecay[3] = new TH2F("hyptBsAllDecay","Bs - All",80,0.,40.,40,-2.,2.);
  fHistBYPtAllDecay[4] = new TH2F("hyptLbAllDecay","LB - All",40,0.,40.,40,-2.,2.);

  fHistYPtAllDecay[0] = new TH2F("hyptD0AllDecay","D0 - All",80,0.,40.,40,-2.,2.);
  fHistYPtAllDecay[1] = new TH2F("hyptDplusAllDecay","Dplus - All",80,0.,40.,40,-2.,2.);
  fHistYPtAllDecay[2] = new TH2F("hyptDstarAllDecay","Dstar - All",80,0.,40.,40,-2.,2.);
  fHistYPtAllDecay[3] = new TH2F("hyptDsAllDecay","Ds - All",80,0.,40.,40,-2.,2.);
  fHistYPtAllDecay[4] = new TH2F("hyptLcAllDecay","Lc - All",80,0.,40.,40,-2.,2.);

  fHistYPtPromptAllDecay[0] = new TH2F("hyptD0promptAllDecay","D0 - Prompt",80,0.,40.,40,-2.,2.);
  fHistYPtPromptAllDecay[1] = new TH2F("hyptDpluspromptAllDecay","Dplus - Prompt",80,0.,40.,40,-2.,2.);
  fHistYPtPromptAllDecay[2] = new TH2F("hyptDstarpromptAllDecay","Dstar - Prompt",80,0.,40.,40,-2.,2.);
  fHistYPtPromptAllDecay[3] = new TH2F("hyptDspromptAllDecay","Ds - Prompt",80,0.,40.,40,-2.,2.);
  fHistYPtPromptAllDecay[4] = new TH2F("hyptLcpromptAllDecay","Lc - Prompt",80,0.,40.,40,-2.,2.);

  fHistYPtFeeddownAllDecay[0] = new TH2F("hyptD0feeddownAllDecay","D0 - FromB",80,0.,40.,40,-2.,2.);
  fHistYPtFeeddownAllDecay[1] = new TH2F("hyptDplusfeeddownAllDecay","Dplus - FromB",80,0.,40.,40,-2.,2.);
  fHistYPtFeeddownAllDecay[2] = new TH2F("hyptDstarfeeddownAllDecay","Dstar - FromB",80,0.,40.,40,-2.,2.);
  fHistYPtFeeddownAllDecay[3] = new TH2F("hyptDsfeeddownAllDecay","Ds - FromB",80,0.,40.,40,-2.,2.);
  fHistYPtFeeddownAllDecay[4] = new TH2F("hyptLcfeeddownAllDecay","Lc - FromB",80,0.,40.,40,-2.,2.);


  fHistYPtFeeddown[0] = new TH2F("hyptD0feeddown","D0 - Feeddown",80,0.,40.,20,-2.,2.);
  fHistYPtFeeddown[1] = new TH2F("hyptDplusfeeddown","Dplus - Feeddown",80,0.,40.,20,-2.,2.);
  fHistYPtFeeddown[2] = new TH2F("hyptDstarfeedown","Dstar - Feeddown",80,0.,40.,20,-2.,2.);
  fHistYPtFeeddown[3] = new TH2F("hyptDsfeedown","Ds - Feeddown",80,0.,40.,20,-2.,2.);
  fHistYPtFeeddown[4] = new TH2F("hyptLcfeedown","Lc - Feeddown",80,0.,40.,20,-2.,2.);

  for(Int_t ih=0; ih<5; ih++){
    fHistBYPtAllDecay[ih]->Sumw2();
    fHistBYPtAllDecay[ih]->SetMinimum(0);
    fOutput->Add(fHistBYPtAllDecay[ih]);
    fHistYPtAllDecay[ih]->Sumw2();
    fHistYPtAllDecay[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtAllDecay[ih]);
    fHistYPtPromptAllDecay[ih]->Sumw2();
    fHistYPtPromptAllDecay[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtPromptAllDecay[ih]);
    fHistYPtFeeddownAllDecay[ih]->Sumw2();
    fHistYPtFeeddownAllDecay[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtFeeddownAllDecay[ih]);
    fHistYPtPrompt[ih]->Sumw2();
    fHistYPtPrompt[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtPrompt[ih]);
    fHistYPtFeeddown[ih]->Sumw2();
    fHistYPtFeeddown[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtFeeddown[ih]);
  }

  fHistYPtD0byDecChannel[0] = new TH2F("hyptD02","D0 - 2prong",80,0.,40.,20,-2.,2.);
  fHistYPtD0byDecChannel[1] = new TH2F("hyptD04","D0 - 4prong",80,0.,40.,20,-2.,2.);
  fHistYPtDplusbyDecChannel[0] = new TH2F("hyptDplusnonreson","Dplus - non reson",80,0.,40.,20,-2.,2.);
  fHistYPtDplusbyDecChannel[1] = new TH2F("hyptDplusreson","Dplus - reson via K0*",80,0.,40.,20,-2.,2.);
  fHistYPtDplusbyDecChannel[2] = new TH2F("hyptDplusKKpi","Dplus -> KKpi",80,0.,40.,20,-2.,2.);
  fHistYPtDsbyDecChannel[0] = new TH2F("hyptDsphi","Ds - vis Phi",80,0.,40.,20,-2.,2.);
  fHistYPtDsbyDecChannel[1] = new TH2F("hyptDsk0st","Ds - via k0*",80,0.,40.,20,-2.,2.);

  for(Int_t ih=0; ih<2; ih++){

    fHistYPtD0byDecChannel[ih]->Sumw2();
    fHistYPtD0byDecChannel[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtD0byDecChannel[ih]);
    fHistYPtDplusbyDecChannel[ih]->Sumw2();
    fHistYPtDplusbyDecChannel[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtDplusbyDecChannel[ih]);
    fHistYPtDsbyDecChannel[ih]->Sumw2();
    fHistYPtDsbyDecChannel[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtDsbyDecChannel[ih]);
  }
  fHistYPtDplusbyDecChannel[2]->Sumw2();
  fHistYPtDplusbyDecChannel[2]->SetMinimum(0);
  fOutput->Add(fHistYPtDplusbyDecChannel[2]);

  fHistOriginPrompt=new TH1F("hOriginPrompt","",100,0.,0.5);
  fHistOriginPrompt->Sumw2();
  fHistOriginPrompt->SetMinimum(0);
  fOutput->Add(fHistOriginPrompt);
  fHistOriginFeeddown=new TH1F("hOriginFeeddown","",100,0.,0.5);
  fHistOriginFeeddown->Sumw2();
  fHistOriginFeeddown->SetMinimum(0);
  fOutput->Add(fHistOriginFeeddown);
  fHistMotherID=new TH1F("hMotherID","",1000,-1.5,998.5);
  fHistMotherID->SetMinimum(0);
  fOutput->Add(fHistMotherID);
  fHistDSpecies=new TH1F("hDSpecies","",10,-0.5,9.5);
  fHistDSpecies->GetXaxis()->SetBinLabel(1,"D0");
  fHistDSpecies->GetXaxis()->SetBinLabel(2,"D0bar");
  fHistDSpecies->GetXaxis()->SetBinLabel(3,"D+");
  fHistDSpecies->GetXaxis()->SetBinLabel(4,"D-");
  fHistDSpecies->GetXaxis()->SetBinLabel(5,"D*+");
  fHistDSpecies->GetXaxis()->SetBinLabel(6,"D*-");
  fHistDSpecies->GetXaxis()->SetBinLabel(7,"Ds+");
  fHistDSpecies->GetXaxis()->SetBinLabel(8,"Ds-");
  fHistDSpecies->GetXaxis()->SetBinLabel(9,"Lc+");
  fHistDSpecies->GetXaxis()->SetBinLabel(10,"Lc-");
  fHistDSpecies->SetMinimum(0);
  fOutput->Add(fHistDSpecies);
  fHistBSpecies=new TH1F("hBSpecies","",10,-0.5,9.5);
  fHistBSpecies->GetXaxis()->SetBinLabel(1,"B0");
  fHistBSpecies->GetXaxis()->SetBinLabel(2,"B0bar");
  fHistBSpecies->GetXaxis()->SetBinLabel(3,"B+");
  fHistBSpecies->GetXaxis()->SetBinLabel(4,"B-");
  fHistBSpecies->GetXaxis()->SetBinLabel(5,"B*+");
  fHistBSpecies->GetXaxis()->SetBinLabel(6,"B*-");
  fHistBSpecies->GetXaxis()->SetBinLabel(7,"Bs+");
  fHistBSpecies->GetXaxis()->SetBinLabel(8,"Bs-");
  fHistBSpecies->GetXaxis()->SetBinLabel(9,"Lb+");
  fHistBSpecies->GetXaxis()->SetBinLabel(10,"Lb-");
  fHistBSpecies->SetMinimum(0);
  fOutput->Add(fHistBSpecies);
  fHistLcDecayChan=new TH1F("hLcDecayChan","",9,-2.5,6.5);
  fHistLcDecayChan->GetXaxis()->SetBinLabel(1,"Violates p cons");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(2,"Other decay");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(3,"Error");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(4,"pK#pi non res");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(5,"pK#pi via K*0");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(6,"pK#pi via #Delta++");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(7,"pK#pi via #Lambda1520");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(8,"pK0s#rightarrowp#pi#pi");
  fHistLcDecayChan->GetXaxis()->SetBinLabel(9,"#pi#Lambda#rightarrowp#pi#pi");
  fHistLcDecayChan->SetMinimum(0);
  fOutput->Add(fHistLcDecayChan);

  fHistNcollHFtype=new TH2F("hNcollHFtype","",5,-1.5,3.5,30,-0.5,29.5);
  fOutput->Add(fHistNcollHFtype);

  Double_t binseta[11]={-1.0,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.0};
  const Int_t nBinsPhi=40;
  Double_t binsphi[nBinsPhi+1];
  for(Int_t ib=0; ib<=nBinsPhi; ib++) binsphi[ib]=ib*TMath::Pi()/20.;
  const Int_t nBinsPt=24;  
  Double_t binspt[nBinsPt+1]={0.,0.10,0.15,0.2,0.25,
			      0.3,0.4,0.5,0.6,0.7,
			      0.8,0.9,1.,1.25,1.5,
			      1.75,2.,2.5,3.,4.,
			      5.,7.5,10.,15.,20.};

   fHistEtaPhiPtGenEle=new TH3F("hEtaPhiPtGenEle","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtGenEle);
  fHistEtaPhiPtGenPi=new TH3F("hEtaPhiPtGenPi","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtGenPi);
  fHistEtaPhiPtGenK=new TH3F("hEtaPhiPtGenK","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtGenK);
  fHistEtaPhiPtGenPro=new TH3F("hEtaPhiPtGenPro","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtGenPro);


  fHistEtaPhiPtRecEle=new TH3F("hEtaPhiPtRecEle","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtRecEle);
  fHistEtaPhiPtRecPi=new TH3F("hEtaPhiPtRecPi","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtRecPi);
  fHistEtaPhiPtRecK=new TH3F("hEtaPhiPtRecK","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtRecK);
  fHistEtaPhiPtRecPro=new TH3F("hEtaPhiPtRecPro","",10,binseta,nBinsPhi,binsphi,nBinsPt,binspt);
  fOutput->Add(fHistEtaPhiPtRecPro);

 

  PostData(1,fOutput);

}
//______________________________________________________________________________
void AliAnalysisTaskCheckHFMCProd2::UserExec(Option_t *)
{
  //

  AliESDEvent *esd = (AliESDEvent*) (InputEvent());


  if(!esd) {
    printf("AliAnalysisTaskSDDRP::Exec(): bad ESD\n");
    return;
  } 

  fHistoNEvents->Fill(0);

  if(!fESDtrackCuts){
    Int_t year=2011;
    if(esd->GetRunNumber()<=139517) year=2010;
    if(year==2010) fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE);
    else fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
    fESDtrackCuts->SetMaxDCAToVertexXY(2.4);
    fESDtrackCuts->SetMaxDCAToVertexZ(3.2);
    fESDtrackCuts->SetDCAToVertex2D(kTRUE);
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					    AliESDtrackCuts::kAny);
  }

  Int_t nTracks=esd->GetNumberOfTracks();
  fHistoTracks->Fill(nTracks);
  Int_t nSelTracks=0;
  for(Int_t it=0; it<nTracks; it++){
    AliESDtrack* tr=esd->GetTrack(it);
    UInt_t status=tr->GetStatus();
    if(!(status&AliESDtrack::kITSrefit)) continue;
    if(!(status&AliESDtrack::kTPCin)) continue;
    nSelTracks++;
  }
  fHistoSelTracks->Fill(nSelTracks);

  const AliMultiplicity* mult=esd->GetMultiplicity();
  Int_t nTracklets=mult->GetNumberOfTracklets();
  Int_t nTracklets1=0;
  for(Int_t it=0; it<nTracklets; it++){
    Double_t eta=TMath::Abs(mult->GetEta(it));
    if(eta<1) nTracklets1++;
  }
  fHistoTracklets->Fill(nTracklets);
  fHistoTrackletsEta1->Fill(nTracklets1);
  
  const AliESDVertex *spdv=esd->GetVertex();
  if(spdv && spdv->IsFromVertexer3D()){
    fHistoSPD3DVtxX->Fill(spdv->GetX());
    fHistoSPD3DVtxY->Fill(spdv->GetY());
    fHistoSPD3DVtxZ->Fill(spdv->GetZ());
  }
  if(spdv && spdv->IsFromVertexerZ()){
    fHistoSPDZVtxX->Fill(spdv->GetX());
    fHistoSPDZVtxY->Fill(spdv->GetY());
    fHistoSPDZVtxZ->Fill(spdv->GetZ());
  }
  const AliESDVertex *trkv=esd->GetPrimaryVertex();
  if(trkv && trkv->GetNContributors()>1){
    fHistoTRKVtxX->Fill(trkv->GetX());
    fHistoTRKVtxY->Fill(trkv->GetY());
    fHistoTRKVtxZ->Fill(trkv->GetZ());
  }

  AliStack* stack=0;
  if(fReadMC){
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
    stack = mcEvent->Stack();
    if (!stack) {
      Printf("ERROR: stack not available");
      return;
    }
    const AliVVertex* mcVert=mcEvent->GetPrimaryVertex();
    if(!mcVert){
      Printf("ERROR: generated vertex not available");
      return;
    }
    if(TMath::Abs(mcVert->GetZ())>10) return;

    //    const AliHeader* h=(AliHeader*)mcEvent->GetHeader();
    //    cout<<h<<endl;
    TString genname=mcEvent->GenEventHeader()->ClassName();
    Int_t nColl=-1;
    Int_t typeHF=-1;
    TList* lgen=0x0;
    if(genname.Contains("CocktailEventHeader")){
      AliGenCocktailEventHeader *cockhead=(AliGenCocktailEventHeader*)mcEvent->GenEventHeader();
      lgen=cockhead->GetHeaders();
      for(Int_t ig=0; ig<lgen->GetEntries(); ig++){
	AliGenerator* gen=(AliGenerator*)lgen->At(ig);
	TString title=gen->GetName();
	if(title.Contains("bchadr")) typeHF=1;
	else if(title.Contains("chadr")) typeHF=0;
	else if(title.Contains("bele")) typeHF=3;
 	else if(title.Contains("cele")) typeHF=2;
      }
      nColl=lgen->GetEntries();
      fHistNcollHFtype->Fill(typeHF,nColl);
    }else{
      TString genTitle=mcEvent->GenEventHeader()->GetTitle();
      if(genTitle.Contains("bchadr")) typeHF=1;
      else if(genTitle.Contains("chadr")) typeHF=0;
      else if(genTitle.Contains("bele")) typeHF=3;
      else if(genTitle.Contains("cele")) typeHF=2;
      fHistNcollHFtype->Fill(typeHF,1.);
    }
    Int_t nParticles=stack->GetNtrack();
    Double_t dNchdy = 0.;
    Int_t nb = 0, nc=0;
    Int_t nCharmed=0;
    Int_t nPhysPrim=0;
    for (Int_t i=0;i<nParticles;i++){
      TParticle* part = (TParticle*)stack->Particle(i);
      Int_t absPdg=TMath::Abs(part->GetPdgCode());
      Int_t pdg=part->GetPdgCode();
      if(absPdg==4) nc++;
      if(absPdg==5) nb++;
      if(stack->IsPhysicalPrimary(i)){
	Double_t eta=part->Eta();
	fHistoEtaPhysPrim->Fill(eta);
	if(absPdg==11) fHistEtaPhiPtGenEle->Fill(eta,part->Phi(),part->Pt());
	else if(absPdg==211) fHistEtaPhiPtGenPi->Fill(eta,part->Phi(),part->Pt());
	else if(absPdg==321) fHistEtaPhiPtGenK->Fill(eta,part->Phi(),part->Pt());
	else if(absPdg==2212) fHistEtaPhiPtGenPro->Fill(eta,part->Phi(),part->Pt());
	
	if(TMath::Abs(eta)<0.5){
	  dNchdy+=0.6666;   // 2/3 for the ratio charged/all
	  nPhysPrim++;
	}
	if(TMath::Abs(eta)<0.9){
	  fHistoPtPhysPrim->Fill(part->Pt());
	}
      }
      Float_t rapid=-999.;
      if (part->Energy() != TMath::Abs(part->Pz())){
	rapid=0.5*TMath::Log((part->Energy()+part->Pz())/(part->Energy()-part->Pz()));
      }

       
      Int_t iPart=-1;
      Int_t iType=0;
      Int_t iSpecies=-1;
      Int_t dummy[4];
      if(absPdg==421){
	iSpecies=0;
	iType=AliVertexingHFUtils::CheckD0Decay(stack,i,dummy); 
	if(iType>0) iPart=0;	
      }
      else if(absPdg==411){
	iSpecies=1;
	iType=AliVertexingHFUtils::CheckDplusDecay(stack,i,dummy);
	if(iType<0){
	  Int_t iTypeKKpi=AliVertexingHFUtils::CheckDplusKKpiDecay(stack,i,dummy);
	  if(iTypeKKpi>0) iType=3;
	}
	if(iType>0) iPart=1;
      }
      else if(absPdg==413){
	iSpecies=2;
	iType=AliVertexingHFUtils::CheckDstarDecay(stack,i,dummy);
	if(iType>0) iPart=2;
      }
      else if(absPdg==431){
	iSpecies=3;
	iType=AliVertexingHFUtils::CheckDsDecay(stack,i,dummy);
	if(iType==1 || iType==2) iPart=3;
      }
      else if(absPdg==4122){
	iSpecies=4;
	iType=AliVertexingHFUtils::CheckLcpKpiDecay(stack,i,dummy);
	if(iType<0){
	  Int_t iTypeV0=AliVertexingHFUtils::CheckLcV0bachelorDecay(stack,i,dummy);
	  if(iTypeV0==1) iType=5;
	  if(iTypeV0==2) iType=6;
	}
	fHistLcDecayChan->Fill(iType);
	if(iType>=0) iPart=4;
      }
      if(iSpecies>=0) fHistYPtAllDecay[iSpecies]->Fill(part->Pt(),rapid);

      // check beauty mesons
      if(absPdg==511) fHistBYPtAllDecay[0]->Fill(part->Pt(),rapid);
      else if(absPdg==521) fHistBYPtAllDecay[1]->Fill(part->Pt(),rapid);
      else if(absPdg==513) fHistBYPtAllDecay[2]->Fill(part->Pt(),rapid);
      else if(absPdg==531) fHistBYPtAllDecay[3]->Fill(part->Pt(),rapid);
      else if(absPdg==5122) fHistBYPtAllDecay[4]->Fill(part->Pt(),rapid);

      if(pdg==511) fHistBSpecies->Fill(0);
      else if(pdg==-511) fHistBSpecies->Fill(1);
      else if(pdg==521) fHistBSpecies->Fill(2);
      else if(pdg==-521) fHistBSpecies->Fill(3);
      else if(pdg==513) fHistBSpecies->Fill(4);
      else if(pdg==-513) fHistBSpecies->Fill(5);
      else if(pdg==531) fHistBSpecies->Fill(6);
      else if(pdg==-531) fHistBSpecies->Fill(7);
      else if(pdg==5122) fHistBSpecies->Fill(8);
      else if(pdg==-5122) fHistBSpecies->Fill(9);

     if(iSpecies<0) continue; // not a charm meson

      if(pdg==421) fHistDSpecies->Fill(0);
      else if(pdg==-421) fHistDSpecies->Fill(1);
      else if(pdg==411) fHistDSpecies->Fill(2);
      else if(pdg==-411) fHistDSpecies->Fill(3);
      else if(pdg==413) fHistDSpecies->Fill(4);
      else if(pdg==-413) fHistDSpecies->Fill(5);
      else if(pdg==431) fHistDSpecies->Fill(6);
      else if(pdg==-431) fHistDSpecies->Fill(7);
      else if(pdg==4122) fHistDSpecies->Fill(8);
      else if(pdg==-4122) fHistDSpecies->Fill(9);

      Double_t distx=part->Vx()-mcVert->GetX();
      Double_t disty=part->Vy()-mcVert->GetY();
      Double_t distz=part->Vz()-mcVert->GetZ();
      Double_t distToVert=TMath::Sqrt(distx*distx+disty*disty+distz*distz);
      fHistMotherID->Fill(part->GetFirstMother());
      Int_t iFromB=AliVertexingHFUtils::CheckOrigin(stack,part,fSearchUpToQuark);
      if(iFromB==4){
	fHistYPtPromptAllDecay[iSpecies]->Fill(part->Pt(),rapid);
	fHistOriginPrompt->Fill(distToVert);
      }
      else if(iFromB==5){
	fHistYPtFeeddownAllDecay[iSpecies]->Fill(part->Pt(),rapid);
	fHistOriginFeeddown->Fill(distToVert);
      }

      if(iPart<0) continue;
      if(iType<0) continue;
      nCharmed++;
      if(iPart==0 && iType>0 && iType<=2){
	fHistYPtD0byDecChannel[iType-1]->Fill(part->Pt(),rapid);
      }else if(iPart==1 && iType>0 && iType<=3){
	fHistYPtDplusbyDecChannel[iType-1]->Fill(part->Pt(),rapid);
      }else if(iPart==3 &&  iType>0 && iType<=2){
	fHistYPtDsbyDecChannel[iType-1]->Fill(part->Pt(),rapid);
      }
      
      if(iFromB==4 && iPart>=0 && iPart<5) fHistYPtPrompt[iPart]->Fill(part->Pt(),rapid);
      else if(iFromB==5 && iPart>=0 && iPart<5) fHistYPtFeeddown[iPart]->Fill(part->Pt(),rapid);      
    }

    for(Int_t i=0; i<nTracks; i++){
      AliESDtrack* track=esd->GetTrack(i);
      if(fESDtrackCuts->AcceptTrack(track)){
	Int_t label=TMath::Abs(track->GetLabel());
	if(stack->IsPhysicalPrimary(label)){
	  TParticle* part = (TParticle*)stack->Particle(label);
	  Int_t absPdg=TMath::Abs(part->GetPdgCode());
	  Double_t eta=part->Eta();
	  if(absPdg==11) fHistEtaPhiPtRecEle->Fill(eta,part->Phi(),part->Pt());
	  else if(absPdg==211) fHistEtaPhiPtRecPi->Fill(eta,part->Phi(),part->Pt());
	  else if(absPdg==321) fHistEtaPhiPtRecK->Fill(eta,part->Phi(),part->Pt());
	  else if(absPdg==2212) fHistEtaPhiPtRecPro->Fill(eta,part->Phi(),part->Pt());      
	}
      }
    }
    fHistoNcharmed->Fill(dNchdy,nCharmed);
    fHistoNbVsNc->Fill(nc,nb);
    fHistoPhysPrim->Fill(nPhysPrim);
  }

  PostData(1,fOutput);
  
}
//______________________________________________________________________________
void AliAnalysisTaskCheckHFMCProd2::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  return;
}




//______________________________________________________________________________
Int_t AliAnalysisTaskCheckHFMCProd2::CheckLcDecay(Int_t labLc, AliStack* stack) const{
  if(labLc<0) return -1;
  TParticle* dp = (TParticle*)stack->Particle(labLc);
  Int_t pdgdp=dp->GetPdgCode();
  Int_t nDau=dp->GetNDaughters();

  if(nDau==3){
    Int_t nKaons=0;
    Int_t nPions=0;
    Int_t nProtons=0;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
      TParticle* dau=(TParticle*)stack->Particle(iDau);
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	if(pdgdp>0 && pdgdau>0) return -1;
	if(pdgdp<0 && pdgdau<0) return -1;
	nKaons++;
      }else if(TMath::Abs(pdgdau)==211){
	if(pdgdp<0 && pdgdau>0) return -1;
	if(pdgdp>0 && pdgdau<0) return -1;
	nPions++;
      }else if(TMath::Abs(pdgdau)==2212){
	if(pdgdp<0 && pdgdau>0) return -1;
	if(pdgdp>0 && pdgdau<0) return -1;
	nProtons++;
      }else{
	return -1;
      }
    }
    if(nPions!=1) return -1;
    if(nKaons!=1) return -1;
    if(nProtons!=1) return -1;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
    }
    return 0;
  }

  if(nDau==2){
    Int_t nKaons=0;
    Int_t nPions=0;
    Int_t nProtons=0;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
      TParticle* dau=(TParticle*)stack->Particle(iDau);
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==313 || TMath::Abs(pdgdau)==3124 || TMath::Abs(pdgdau)==2224
	 || TMath::Abs(pdgdau)==3122 || TMath::Abs(pdgdau)==311){
	Int_t nDauRes=dau->GetNDaughters();
	if(nDauRes!=2)  return -1;
	for(Int_t resDau=dau->GetFirstDaughter(); resDau<=dau->GetLastDaughter(); resDau++){
	  if(resDau<0) return -1;
	  TParticle* resdau=(TParticle*)stack->Particle(resDau);
	  Int_t pdgresdau=resdau->GetPdgCode();
	  if(TMath::Abs(pdgresdau)==321){
	    if(pdgdp>0 && pdgresdau>0) return -1;
	    if(pdgdp<0 && pdgresdau<0) return -1;
	    nKaons++;
	  }
	  if(TMath::Abs(pdgresdau)==211){
	    if(pdgdp<0 && pdgresdau>0) return -1;
	    if(pdgdp>0 && pdgresdau<0) return -1;
	    nPions++;
	  }
	  if(TMath::Abs(pdgresdau)==2212){
	    if(pdgdp<0 && pdgresdau>0) return -1;
	    if(pdgdp>0 && pdgresdau<0) return -1;
	    nProtons++;
	  }
	}
      }else if(TMath::Abs(pdgdau)==321){
	if(pdgdp>0 && pdgdau>0) return -1;
	if(pdgdp<0 && pdgdau<0) return -1;
	nKaons++;
      }else if(TMath::Abs(pdgdau)==211){
	if(pdgdp<0 && pdgdau>0) return -1;
	if(pdgdp>0 && pdgdau<0) return -1;
	nPions++;
      }else if(TMath::Abs(pdgdau)==2212){
	if(pdgdp<0 && pdgdau>0) return -1;
	if(pdgdp>0 && pdgdau<0) return -1;
	nProtons++;
      }else{
	return -1;
      }
    }
    if(nPions!=1) return -1;
    if(nKaons!=1) return -1;
    if(nProtons!=1) return -1;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
    }
    return 1;
  } 
  return -1;
}

