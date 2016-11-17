#include <Riostream.h>
#include <TFile.h>
//#include <AliRDHFCutsDplustoKpipi.h>
#include <TClonesArray.h>
#include <TParameter.h>

//macro to make a .root file which contains an AliRDHFCutsDplustoKpipi with loose set of cuts (for significance maximization) and TParameter with the tighest value of these cuts
//Needed for AliAnalysisTaskSEDplus, AliCFTaskVertexingHF3Prong, AliAnalysisTaskSESignificance

void makeTFile4CutsDplustoKpipi_010(Bool_t fUseMC=kFALSE);
void makeTFile4CutsDplustoKpipi_020(Bool_t fUseMC=kFALSE);
void makeTFile4CutsDplustoKpipi_2040(Bool_t fUseMC=kFALSE);
void makeTFile4CutsDplustoKpipi_3050(Bool_t fUseMC=kFALSE);
void makeTFile4CutsDplustoKpipi_6080(Bool_t fUseMC=kFALSE);

//__________________________________________________________________________________________
void makeTFile4CutsDplustoKpipi_010(Bool_t fUseMC){

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4); 
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                         AliESDtrackCuts::kAny);
  // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.6,1.e10);
  
  TString cent="";
  //centrality selection (Pb-Pb)
  Float_t minc=0.,maxc=10;
  const Int_t nptbins=15;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  
  ptbins[0]=2.;
  ptbins[1]=3.;
  ptbins[2]=4.;
  ptbins[3]=5.;
  ptbins[4]=6.;
  ptbins[5]=7.;
  ptbins[6]=8.;
  ptbins[7]=9.;
  ptbins[8]=10.;
  ptbins[9]=11.;
  ptbins[10]=12.;
  ptbins[11]=14.;
  ptbins[12]=16.;
  ptbins[13]=24.;
  ptbins[14]=36.;
  ptbins[15]=100.;
  
  const Int_t nvars=14;
  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
  
  
  Int_t ic=0;//minv
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.2;
  }
  
  
  ic=1;//ptK
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[1][ipt]=0.6;
  }
  
  ic=2;//ptPi
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.6;
  }
  
  ic=3;//d0K
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=4;//d0Pi
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=5;//dist12
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  
  ic=6;//sigvert
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.03;
  }
  
  ic=7;//declen
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.08;
  }
  
  ic=8;//pM
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.0;
  }
  
  //cosp
  ic=9;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.99;
  }
  anacutsval[ic][0]=0.992;
  anacutsval[ic][1]=0.992;
  anacutsval[ic][2]=0.992;
  anacutsval[ic][10]=0.98;
  anacutsval[ic][11]=0.98;
  anacutsval[ic][12]=0.97;
  anacutsval[ic][13]=0.97;
  anacutsval[ic][14]=0.97;
  
  ic=10;//sumd02
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.0;
  }
  
  ic=11;//dca
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=10000000000.;
  }
  
  ic=12;//ndlXY
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=6.;
  }
  anacutsval[ic][0]=10.;
  anacutsval[ic][1]=10.;
  anacutsval[ic][2]=9.;
  anacutsval[ic][3]=8.;
  anacutsval[ic][4]=8.;
  anacutsval[ic][5]=8.;
  anacutsval[ic][13]=5.;
  anacutsval[ic][14]=5.;
  
  ic=13;//cospXY
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.99;
  }
  anacutsval[ic][0]=0.992;
  anacutsval[ic][1]=0.992;
  anacutsval[ic][2]=0.992;
  anacutsval[ic][3]=0.992;
  anacutsval[ic][4]=0.992;
  anacutsval[ic][5]=0.992;
  anacutsval[ic][10]=0.98;
  anacutsval[ic][11]=0.98;
  anacutsval[ic][12]=0.97;
  anacutsval[ic][13]=0.97;
  anacutsval[ic][14]=0.97;
  
  Float_t *d0cutsval=new Float_t[nptbins];
  for(Int_t ipt=0;ipt<nptbins;ipt++){ //d0
    d0cutsval[ipt]=60;
  }
  d0cutsval[0]=80;
  d0cutsval[11]=40;
  d0cutsval[12]=40;
  d0cutsval[13]=40;
  d0cutsval[14]=40;
  
  AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
  analysiscuts->SetName("AnalysisCuts");
  analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
  analysiscuts->SetPtBins(nptbins+1,ptbins);
  analysiscuts->SetCuts(nvars,nptbins,anacutsval);
  analysiscuts->Setd0Cut(nptbins,d0cutsval);
  analysiscuts->AddTrackCuts(esdTrackCuts);
  analysiscuts->SetScaleNormDLxyBypOverPt(kFALSE);
  //   AliAODPidHF* pidHF=(AliAODPidHF*)analysiscuts->GetPidHF();
  //pidHF->SetOldPid(kFALSE);
  //analysiscuts->SetPidHF(pidHF);
  //analysiscuts->GetPidHF()->SetOldPid(kFALSE);
  cout<<"************** checking old PID (it should be FALSE by default - July 10)--> "<<analysiscuts->GetPidHF()->GetOldPid()<<endl;
  
  // analysiscuts->SetUsePID(kFALSE);
  analysiscuts->SetUsePID(kTRUE);
  // analysiscuts->GetPidHF()->SetMaxTrackMomForCombinedPID(4.);
  //cout<<"++++++++ Max Pt for PID --> "<<analysiscuts->GetPidHF()->GetMaxTrackMomForCombinedPID()<<endl;
  analysiscuts->SetUseImpParProdCorrCut(kFALSE);
  analysiscuts->SetOptPileup(kFALSE);
  //analysiscuts->SetUseAOD049(kTRUE);
  analysiscuts->SetMinCentrality(minc);
  analysiscuts->SetMaxCentrality(maxc);
  
  analysiscuts->SetRemoveTrackletOutliers(kTRUE);//added on June 28
  analysiscuts->SetCutOnzVertexSPD(0);//added on July19
  
  cent=Form("%.0f%.0f",minc,maxc);
  analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
  if(fUseMC) {
    analysiscuts->SetTriggerClass("");
    analysiscuts->ResetMaskAndEnableMBTrigger();
    analysiscuts->SetTriggerMask(AliVEvent::kMB);
  }
  else {
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    analysiscuts->SetTriggerMask(AliVEvent::kINT7);
  }
  
  // analysiscuts->EnableSemiCentralTrigger();
  analysiscuts->SetMinPtCandidate(2.);
  analysiscuts->SetMaxPtCandidate(100.);
  
  cout<<"This is the object I'm going to save:"<<nptbins<<endl;
  
  analysiscuts->PrintAll();
  analysiscuts->PrintTrigger();
  TString filename="DplustoKpipiCuts_010_central_d0cut_kINT7.root";
  if(fUseMC) filename="DplustoKpipiCuts_010_central_d0cut_MC.root";
  TFile* fout=new TFile(filename.Data(),"RECREATE");
  fout->cd();
  analysiscuts->Write();
  fout->Close();
  
}

//__________________________________________________________________________________________
void makeTFile4CutsDplustoKpipi_020(Bool_t fUseMC){
  
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4); 
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                         AliESDtrackCuts::kAny);
  // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.6,1.e10);
  
  TString cent="";
  //centrality selection (Pb-Pb)
  Float_t minc=0.,maxc=20;
  const Int_t nptbins=15;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  
  ptbins[0]=2.;
  ptbins[1]=3.;
  ptbins[2]=4.;
  ptbins[3]=5.;
  ptbins[4]=6.;
  ptbins[5]=7.;
  ptbins[6]=8.;
  ptbins[7]=9.;
  ptbins[8]=10.;
  ptbins[9]=11.;
  ptbins[10]=12.;
  ptbins[11]=14.;
  ptbins[12]=16.;
  ptbins[13]=24.;
  ptbins[14]=36.;
  ptbins[15]=100.;
  
  const Int_t nvars=14;
  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
  
  
  Int_t ic=0;//minv
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.2;
  }
  
  
  ic=1;//ptK
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[1][ipt]=0.6;
  }
  
  ic=2;//ptPi
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.6;
  }
  
  ic=3;//d0K
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=4;//d0Pi
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=5;//dist12
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  
  ic=6;//sigvert
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.03;
  }
  
  ic=7;//declen
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.08;
  }
  
  ic=8;//pM
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.0;
  }
  
  //cosp
  ic=9;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.99;
  }
  anacutsval[ic][0]=0.992;
  anacutsval[ic][1]=0.992;
  anacutsval[ic][2]=0.992;
  anacutsval[ic][10]=0.98;
  anacutsval[ic][11]=0.98;
  anacutsval[ic][12]=0.97;
  anacutsval[ic][13]=0.97;
  anacutsval[ic][14]=0.97;
  
  ic=10;//sumd02
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.0;
  }
  
  ic=11;//dca
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=10000000000.;
  }
  
  ic=12;//ndlXY
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=6.;
  }
  anacutsval[ic][0]=10.;
  anacutsval[ic][1]=10.;
  anacutsval[ic][2]=9.;
  anacutsval[ic][3]=8.;
  anacutsval[ic][4]=8.;
  anacutsval[ic][5]=8.;
  anacutsval[ic][13]=5.;
  anacutsval[ic][14]=5.;
  
  ic=13;//cospXY
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.99;
  }
  anacutsval[ic][0]=0.992;
  anacutsval[ic][1]=0.992;
  anacutsval[ic][2]=0.992;
  anacutsval[ic][3]=0.992;
  anacutsval[ic][4]=0.992;
  anacutsval[ic][5]=0.992;
  anacutsval[ic][10]=0.98;
  anacutsval[ic][11]=0.98;
  anacutsval[ic][12]=0.97;
  anacutsval[ic][13]=0.97;
  anacutsval[ic][14]=0.97;
  
  Float_t *d0cutsval=new Float_t[nptbins];
  for(Int_t ipt=0;ipt<nptbins;ipt++){ //d0
    d0cutsval[ipt]=60;
  }
  d0cutsval[0]=80;
  d0cutsval[11]=40;
  d0cutsval[12]=40;
  d0cutsval[13]=40;
  d0cutsval[14]=40;
  analysiscuts->SetCuts(nvars,nptbins,anacutsval);

  AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
  analysiscuts->SetName("AnalysisCuts");
  analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
  analysiscuts->SetPtBins(nptbins+1,ptbins);
  analysiscuts->SetCuts(nvars,nptbins,anacutsval);
  analysiscuts->Setd0Cut(nptbins,d0cutsval);
  analysiscuts->AddTrackCuts(esdTrackCuts);
  analysiscuts->SetScaleNormDLxyBypOverPt(kFALSE);
  //   AliAODPidHF* pidHF=(AliAODPidHF*)analysiscuts->GetPidHF();
  //pidHF->SetOldPid(kFALSE);
  //analysiscuts->SetPidHF(pidHF);
  //analysiscuts->GetPidHF()->SetOldPid(kFALSE);
  cout<<"************** checking old PID (it should be FALSE by default - July 10)--> "<<analysiscuts->GetPidHF()->GetOldPid()<<endl;
  
  // analysiscuts->SetUsePID(kFALSE);
  analysiscuts->SetUsePID(kTRUE);
  // analysiscuts->GetPidHF()->SetMaxTrackMomForCombinedPID(4.);
  //cout<<"++++++++ Max Pt for PID --> "<<analysiscuts->GetPidHF()->GetMaxTrackMomForCombinedPID()<<endl;
  analysiscuts->SetUseImpParProdCorrCut(kFALSE);
  analysiscuts->SetOptPileup(kFALSE);
  //analysiscuts->SetUseAOD049(kTRUE);
  analysiscuts->SetMinCentrality(minc);
  analysiscuts->SetMaxCentrality(maxc);
  
  analysiscuts->SetRemoveTrackletOutliers(kTRUE);//added on June 28
  analysiscuts->SetCutOnzVertexSPD(0);//added on July19
  
  cent=Form("%.0f%.0f",minc,maxc);
  analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
  if(fUseMC) {
    analysiscuts->SetTriggerClass("");
    analysiscuts->ResetMaskAndEnableMBTrigger();
    analysiscuts->SetTriggerMask(AliVEvent::kMB);
  }
  else {
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    analysiscuts->SetTriggerMask(AliVEvent::kINT7);
  }
  
  // analysiscuts->EnableSemiCentralTrigger();
  analysiscuts->SetMinPtCandidate(2.);
  analysiscuts->SetMaxPtCandidate(100.);
  
  cout<<"This is the object I'm going to save:"<<nptbins<<endl;
  
  analysiscuts->PrintAll();
  analysiscuts->PrintTrigger();
  TString filename="DplustoKpipiCuts_020_central_d0cut_kINT7.root";
  if(fUseMC) filename="DplustoKpipiCuts_020_central_d0cut_MC.root";
  TFile* fout=new TFile(filename.Data(),"RECREATE");
  fout->cd();
  analysiscuts->Write();
  fout->Close();
  
}

//__________________________________________________________________________________________
void makeTFile4CutsDplustoKpipi_2040(Bool_t fUseMC){
  
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4); 
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                         AliESDtrackCuts::kAny);
  // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.4,1.e10);
  
  TString cent="";
  //centrality selection (Pb-Pb)
  Float_t minc=20.,maxc=40;
  const Int_t nptbins=15;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  
  ptbins[0]=2.;
  ptbins[1]=3.;
  ptbins[2]=4.;
  ptbins[3]=5.;
  ptbins[4]=6.;
  ptbins[5]=7.;
  ptbins[6]=8.;
  ptbins[7]=9.;
  ptbins[8]=10.;
  ptbins[9]=11.;
  ptbins[10]=12.;
  ptbins[11]=14.;
  ptbins[12]=16.;
  ptbins[13]=24.;
  ptbins[14]=36.;
  ptbins[15]=100.;
  
  const Int_t nvars=14;
  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
  
  
  Int_t ic=0;//minv
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.2;
  }
  
  ic=1;//ptK
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[1][ipt]=0.4;
  }
  
  ic=2;//ptPi
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.4;
  }
  
  ic=3;//d0K
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=4;//d0Pi
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=5;//dist12
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  
  ic=6;//sigvert
  anacutsval[ic][0]=0.022;
  anacutsval[ic][1]=0.022;
  anacutsval[ic][2]=0.022;
  anacutsval[ic][3]=0.022;
  anacutsval[ic][4]=0.022;
  anacutsval[ic][5]=0.022;
  anacutsval[ic][6]=0.025;
  anacutsval[ic][7]=0.025;
  anacutsval[ic][8]=0.025;
  anacutsval[ic][9]=0.025;
  anacutsval[ic][10]=0.025;
  anacutsval[ic][11]=0.025;
  anacutsval[ic][12]=0.025;
  anacutsval[ic][13]=0.035;
  anacutsval[ic][14]=0.035;
  
  ic=7;//declen
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.10;
  }
  anacutsval[ic][12]=0.15;
  anacutsval[ic][13]=0.15;
  anacutsval[ic][14]=0.15;
  
  ic=8;//pM
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.0;
  }
  
  //cosp
  ic=9;
  anacutsval[ic][0]=0.997;
  anacutsval[ic][1]=0.997;
  anacutsval[ic][2]=0.996;
  anacutsval[ic][3]=0.996;
  anacutsval[ic][4]=0.996;
  anacutsval[ic][5]=0.996;
  anacutsval[ic][6]=0.992;
  anacutsval[ic][7]=0.992;
  anacutsval[ic][8]=0.992;
  anacutsval[ic][9]=0.992;
  anacutsval[ic][10]=0.992;
  anacutsval[ic][11]=0.992;
  anacutsval[ic][12]=0.980;
  anacutsval[ic][13]=0.970;
  anacutsval[ic][14]=0.950;
  
  ic=10;//sumd02
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.0;
  }
  
  ic=11;//dca
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=10000000000.;
  }
  
  ic=12;//ndlXY
  anacutsval[ic][0]=14.;
  anacutsval[ic][1]=14.;
  anacutsval[ic][2]=14.;
  anacutsval[ic][3]=14.;
  anacutsval[ic][4]=11.;
  anacutsval[ic][5]=11.;
  anacutsval[ic][6]=11.;
  anacutsval[ic][7]=11.;
  anacutsval[ic][8]=10.;
  anacutsval[ic][9]=10.;
  anacutsval[ic][10]=10.;
  anacutsval[ic][11]=10.;
  anacutsval[ic][12]=8.;
  anacutsval[ic][13]=6.;
  anacutsval[ic][14]=6.;
  
  ic=13;//cospXY
  anacutsval[ic][0]=0.997;
  anacutsval[ic][1]=0.997;
  anacutsval[ic][2]=0.996;
  anacutsval[ic][3]=0.996;
  anacutsval[ic][4]=0.996;
  anacutsval[ic][5]=0.996;
  anacutsval[ic][6]=0.992;
  anacutsval[ic][7]=0.992;
  anacutsval[ic][8]=0.992;
  anacutsval[ic][9]=0.992;
  anacutsval[ic][10]=0.992;
  anacutsval[ic][11]=0.992;
  anacutsval[ic][12]=0.980;
  anacutsval[ic][13]=0.970;
  anacutsval[ic][14]=0.950;
  
  Float_t *d0cutsval=new Float_t[nptbins];
  for(Int_t ipt=0;ipt<nptbins;ipt++){ //d0
    d0cutsval[ipt]=60;
  }
  d0cutsval[0]=80;
  d0cutsval[11]=40;
  d0cutsval[12]=40;
  d0cutsval[13]=40;
  d0cutsval[14]=40;
  
  AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
  analysiscuts->SetName("AnalysisCuts");
  analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
  analysiscuts->SetPtBins(nptbins+1,ptbins);
  analysiscuts->SetCuts(nvars,nptbins,anacutsval);
  analysiscuts->Setd0Cut(nptbins,d0cutsval);
  analysiscuts->AddTrackCuts(esdTrackCuts);
  analysiscuts->SetScaleNormDLxyBypOverPt(kFALSE);
  //   AliAODPidHF* pidHF=(AliAODPidHF*)analysiscuts->GetPidHF();
  //pidHF->SetOldPid(kFALSE);
  //analysiscuts->SetPidHF(pidHF);
  //analysiscuts->GetPidHF()->SetOldPid(kFALSE);
  cout<<"************** checking old PID (it should be FALSE by default - July 10)--> "<<analysiscuts->GetPidHF()->GetOldPid()<<endl;
  
  analysiscuts->SetUsePID(kTRUE);
  // analysiscuts->GetPidHF()->SetMaxTrackMomForCombinedPID(4.);
  //cout<<"++++++++ Max Pt for PID --> "<<analysiscuts->GetPidHF()->GetMaxTrackMomForCombinedPID()<<endl;
  analysiscuts->SetUseImpParProdCorrCut(kFALSE);
  analysiscuts->SetOptPileup(kFALSE);
  analysiscuts->SetMinCentrality(minc);
  analysiscuts->SetMaxCentrality(maxc);
  
  analysiscuts->SetRemoveTrackletOutliers(kTRUE);//added on June 28
  analysiscuts->SetCutOnzVertexSPD(0);//added on July19
  
  cent=Form("%.0f%.0f",minc,maxc);
  analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
  if(fUseMC) {
    analysiscuts->SetTriggerClass("");
    analysiscuts->ResetMaskAndEnableMBTrigger();
    analysiscuts->SetTriggerMask(AliVEvent::kMB);
  }
  else {
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    analysiscuts->SetTriggerMask(AliVEvent::kINT7);
  }
  
  // analysiscuts->EnableSemiCentralTrigger();
  analysiscuts->SetMinPtCandidate(2.);
  analysiscuts->SetMaxPtCandidate(100.);
  
  cout<<"This is the object I'm going to save:"<<nptbins<<endl;
  
  analysiscuts->PrintAll();
  analysiscuts->PrintTrigger();
  TString filename="DplustoKpipiCuts_2040_central_d0cut_kINT7.root";
  if(fUseMC) filename="DplustoKpipiCuts_2040_central_d0cut_MC.root";
  TFile* fout=new TFile(filename.Data(),"RECREATE");
  fout->cd();
  analysiscuts->Write();
  fout->Close();
  
}

//__________________________________________________________________________________________
void makeTFile4CutsDplustoKpipi_3050(Bool_t fUseMC){
  
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4); 
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                         AliESDtrackCuts::kAny);
  // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.4,1.e10);
  
  TString cent="";
  //centrality selection (Pb-Pb)
  Float_t minc=30.,maxc=50;
  const Int_t nptbins=15;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  
  ptbins[0]=2.;
  ptbins[1]=3.;
  ptbins[2]=4.;
  ptbins[3]=5.;
  ptbins[4]=6.;
  ptbins[5]=7.;
  ptbins[6]=8.;
  ptbins[7]=9.;
  ptbins[8]=10.;
  ptbins[9]=11.;
  ptbins[10]=12.;
  ptbins[11]=14.;
  ptbins[12]=16.;
  ptbins[13]=24.;
  ptbins[14]=36.;
  ptbins[15]=100.;
  
  const Int_t nvars=14;
  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
  
  
  Int_t ic=0;//minv
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.2;
  }
  
  ic=1;//ptK
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[1][ipt]=0.4;
  }
  
  ic=2;//ptPi
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.4;
  }
  
  ic=3;//d0K
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=4;//d0Pi
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=5;//dist12
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  
  ic=6;//sigvert
  anacutsval[ic][0]=0.020;
  anacutsval[ic][1]=0.022;
  anacutsval[ic][2]=0.022;
  anacutsval[ic][3]=0.022;
  anacutsval[ic][4]=0.022;
  anacutsval[ic][5]=0.022;
  anacutsval[ic][6]=0.025;
  anacutsval[ic][7]=0.025;
  anacutsval[ic][8]=0.025;
  anacutsval[ic][9]=0.025;
  anacutsval[ic][10]=0.025;
  anacutsval[ic][11]=0.025;
  anacutsval[ic][12]=0.025;
  anacutsval[ic][13]=0.035;
  anacutsval[ic][14]=0.035;
  
  ic=7;//declen
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.12;
  }
  anacutsval[ic][0]=0.08;
  anacutsval[ic][12]=0.16;
  anacutsval[ic][13]=0.16;
  anacutsval[ic][14]=0.16;
  
  ic=8;//pM
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.0;
  }
  
  //cosp
  ic=9;
  anacutsval[ic][0]=0.996;
  anacutsval[ic][1]=0.996;
  anacutsval[ic][2]=0.995;
  anacutsval[ic][3]=0.995;
  anacutsval[ic][4]=0.995;
  anacutsval[ic][5]=0.995;
  anacutsval[ic][6]=0.990;
  anacutsval[ic][7]=0.990;
  anacutsval[ic][8]=0.990;
  anacutsval[ic][9]=0.990;
  anacutsval[ic][10]=0.990;
  anacutsval[ic][11]=0.990;
  anacutsval[ic][12]=0.980;
  anacutsval[ic][13]=0.970;
  anacutsval[ic][14]=0.950;
  
  ic=10;//sumd02
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.0;
  }
  
  ic=11;//dca
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=10000000000.;
  }
  
  ic=12;//ndlXY
  anacutsval[ic][0]=12.;
  anacutsval[ic][1]=12.;
  anacutsval[ic][2]=12.;
  anacutsval[ic][3]=10.;
  anacutsval[ic][4]=10.;
  anacutsval[ic][5]=10.;
  anacutsval[ic][6]=10.;
  anacutsval[ic][7]=10.;
  anacutsval[ic][8]=9.;
  anacutsval[ic][9]=9.;
  anacutsval[ic][10]=9.;
  anacutsval[ic][11]=9.;
  anacutsval[ic][12]=8.;
  anacutsval[ic][13]=8.;
  anacutsval[ic][14]=6.;
  
  ic=13;//cospXY
  anacutsval[ic][0]=0.996;
  anacutsval[ic][1]=0.996;
  anacutsval[ic][2]=0.995;
  anacutsval[ic][3]=0.995;
  anacutsval[ic][4]=0.995;
  anacutsval[ic][5]=0.995;
  anacutsval[ic][6]=0.990;
  anacutsval[ic][7]=0.990;
  anacutsval[ic][8]=0.990;
  anacutsval[ic][9]=0.990;
  anacutsval[ic][10]=0.990;
  anacutsval[ic][11]=0.990;
  anacutsval[ic][12]=0.980;
  anacutsval[ic][13]=0.970;
  anacutsval[ic][14]=0.950;
  
  Float_t *d0cutsval=new Float_t[nptbins];
  for(Int_t ipt=0;ipt<nptbins;ipt++){ //d0
    d0cutsval[ipt]=60;
  }
  d0cutsval[0]=80;
  d0cutsval[11]=40;
  d0cutsval[12]=40;
  d0cutsval[13]=40;
  d0cutsval[14]=40;
  
  AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
  analysiscuts->SetName("AnalysisCuts");
  analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
  analysiscuts->SetPtBins(nptbins+1,ptbins);
  analysiscuts->SetCuts(nvars,nptbins,anacutsval);
  analysiscuts->Setd0Cut(nptbins,d0cutsval);
  analysiscuts->AddTrackCuts(esdTrackCuts);
  analysiscuts->SetScaleNormDLxyBypOverPt(kFALSE);
  //   AliAODPidHF* pidHF=(AliAODPidHF*)analysiscuts->GetPidHF();
  //pidHF->SetOldPid(kFALSE);
  //analysiscuts->SetPidHF(pidHF);
  //analysiscuts->GetPidHF()->SetOldPid(kFALSE);
  cout<<"************** checking old PID (it should be FALSE by default - July 10)--> "<<analysiscuts->GetPidHF()->GetOldPid()<<endl;
  
  // analysiscuts->SetUsePID(kFALSE);
  analysiscuts->SetUsePID(kTRUE);
  // analysiscuts->GetPidHF()->SetMaxTrackMomForCombinedPID(4.);
  //cout<<"++++++++ Max Pt for PID --> "<<analysiscuts->GetPidHF()->GetMaxTrackMomForCombinedPID()<<endl;
  analysiscuts->SetUseImpParProdCorrCut(kFALSE);
  analysiscuts->SetOptPileup(kFALSE);
  //analysiscuts->SetUseAOD049(kTRUE);
  analysiscuts->SetMinCentrality(minc);
  analysiscuts->SetMaxCentrality(maxc);
  
  analysiscuts->SetRemoveTrackletOutliers(kTRUE);//added on June 28
  analysiscuts->SetCutOnzVertexSPD(0);//added on July19
  
  cent=Form("%.0f%.0f",minc,maxc);
  analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
  if(fUseMC) {
    analysiscuts->SetTriggerClass("");
    analysiscuts->ResetMaskAndEnableMBTrigger();
    analysiscuts->SetTriggerMask(AliVEvent::kMB);
  }
  else {
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    analysiscuts->SetTriggerMask(AliVEvent::kINT7);
  }
  
  // analysiscuts->EnableSemiCentralTrigger();
  analysiscuts->SetMinPtCandidate(2.);
  analysiscuts->SetMaxPtCandidate(100.);
  
  cout<<"This is the object I'm going to save:"<<nptbins<<endl;
  
  analysiscuts->PrintAll();
  analysiscuts->PrintTrigger();
  TString filename="DplustoKpipiCuts_3050_central_d0cut_kINT7.root";
  if(fUseMC) filename="DplustoKpipiCuts_3050_central_d0cut_MC.root";
  TFile* fout=new TFile(filename.Data(),"RECREATE");
  fout->cd();
  analysiscuts->Write();
  fout->Close();
  
}

//__________________________________________________________________________________________
void makeTFile4CutsDplustoKpipi_6080(Bool_t fUseMC){
  
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4); 
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                         AliESDtrackCuts::kAny);
  // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.4,1.e10);
  
  TString cent="";
  //centrality selection (Pb-Pb)
  Float_t minc=60.,maxc=80;
  const Int_t nptbins=15;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  
  ptbins[0]=2.;
  ptbins[1]=3.;
  ptbins[2]=4.;
  ptbins[3]=5.;
  ptbins[4]=6.;
  ptbins[5]=7.;
  ptbins[6]=8.;
  ptbins[7]=9.;
  ptbins[8]=10.;
  ptbins[9]=11.;
  ptbins[10]=12.;
  ptbins[11]=14.;
  ptbins[12]=16.;
  ptbins[13]=24.;
  ptbins[14]=36.;
  ptbins[15]=100.;
  
  const Int_t nvars=14;
  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
  
  
  Int_t ic=0;//minv
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.2;
  }
  
  ic=1;//ptK
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[1][ipt]=0.4;
  }
  
  ic=2;//ptPi
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.4;
  }
  
  ic=3;//d0K
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=4;//d0Pi
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=5;//dist12
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  
  ic=6;//sigvert
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.032;
  }
  
  ic=7;//declen
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.05;
  }
  
  ic=8;//pM
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.0;
  }
  
  //cosp
  ic=9;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.98;
  }
  anacutsval[ic][0]=0.99;
  anacutsval[ic][1]=0.99;
  anacutsval[ic][2]=0.99;
  anacutsval[ic][3]=0.99;
  anacutsval[ic][4]=0.99;
  anacutsval[ic][5]=0.99;
  anacutsval[ic][13]=0.97;
  anacutsval[ic][14]=0.97;
  
  ic=10;//sumd02
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.0;
  }
  
  ic=11;//dca
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=10000000000.;
  }
  
  ic=12;//ndlXY
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=5.;
  }
  anacutsval[ic][0]=8.;
  anacutsval[ic][1]=8.;
  anacutsval[ic][2]=7.;
  anacutsval[ic][3]=7.;
  anacutsval[ic][4]=6.;
  anacutsval[ic][5]=6.;
  anacutsval[ic][6]=6.;
  anacutsval[ic][7]=6.;
  anacutsval[ic][8]=6.;
  anacutsval[ic][9]=6.;
  
  ic=13;//cospXY
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.98;
  }
  anacutsval[ic][0]=0.99;
  anacutsval[ic][1]=0.99;
  anacutsval[ic][2]=0.99;
  anacutsval[ic][3]=0.99;
  anacutsval[ic][4]=0.99;
  anacutsval[ic][5]=0.99;
  anacutsval[ic][13]=0.97;
  anacutsval[ic][14]=0.97;
  
  Float_t *d0cutsval=new Float_t[nptbins];
  for(Int_t ipt=0;ipt<nptbins;ipt++){ //d0
    d0cutsval[ipt]=60;
  }
  d0cutsval[0]=80;
  d0cutsval[11]=40;
  d0cutsval[12]=40;
  d0cutsval[13]=40;
  d0cutsval[14]=40;
  
  AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
  analysiscuts->SetName("AnalysisCuts");
  analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
  analysiscuts->SetPtBins(nptbins+1,ptbins);
  analysiscuts->SetCuts(nvars,nptbins,anacutsval);
  analysiscuts->Setd0Cut(nptbins,d0cutsval);
  analysiscuts->AddTrackCuts(esdTrackCuts);
  analysiscuts->SetScaleNormDLxyBypOverPt(kFALSE);
  //   AliAODPidHF* pidHF=(AliAODPidHF*)analysiscuts->GetPidHF();
  //pidHF->SetOldPid(kFALSE);
  //analysiscuts->SetPidHF(pidHF);
  //analysiscuts->GetPidHF()->SetOldPid(kFALSE);
  cout<<"************** checking old PID (it should be FALSE by default - July 10)--> "<<analysiscuts->GetPidHF()->GetOldPid()<<endl;
  
  // analysiscuts->SetUsePID(kFALSE);
  analysiscuts->SetUsePID(kTRUE);
  // analysiscuts->GetPidHF()->SetMaxTrackMomForCombinedPID(4.);
  //cout<<"++++++++ Max Pt for PID --> "<<analysiscuts->GetPidHF()->GetMaxTrackMomForCombinedPID()<<endl;
  analysiscuts->SetUseImpParProdCorrCut(kFALSE);
  analysiscuts->SetOptPileup(kFALSE);
  //analysiscuts->SetUseAOD049(kTRUE);
  analysiscuts->SetMinCentrality(minc);
  analysiscuts->SetMaxCentrality(maxc);
  
  analysiscuts->SetRemoveTrackletOutliers(kTRUE);//added on June 28
  analysiscuts->SetCutOnzVertexSPD(0);//added on July19
  
  cent=Form("%.0f%.0f",minc,maxc);
  analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
  if(fUseMC) {
    analysiscuts->SetTriggerClass("");
    analysiscuts->ResetMaskAndEnableMBTrigger();
    analysiscuts->SetTriggerMask(AliVEvent::kMB);
  }
  else {
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    analysiscuts->SetTriggerMask(AliVEvent::kINT7);
  }
  
  // analysiscuts->EnableSemiCentralTrigger();
  analysiscuts->SetMinPtCandidate(2.);
  analysiscuts->SetMaxPtCandidate(100.);
  
  cout<<"This is the object I'm going to save:"<<nptbins<<endl;
  
  analysiscuts->PrintAll();
  analysiscuts->PrintTrigger();
  TString filename="DplustoKpipiCuts_6080_central_d0cut_kINT7.root";
  if(fUseMC) filename="DplustoKpipiCuts_6080_central_d0cut_MC.root";
  TFile* fout=new TFile(filename.Data(),"RECREATE");
  fout->cd();
  analysiscuts->Write();
  fout->Close();
  
}
