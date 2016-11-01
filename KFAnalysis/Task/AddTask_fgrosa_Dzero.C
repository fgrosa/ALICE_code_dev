AliAnalysisTask *AddTask_fgrosa_Dzero() {
  
  // Test macro for the AliAnalysisTaskDmesonKF for D0 and D+ candidates using the Kalman Filter
                                        
  //  Fabrizio Grosa, grosa@to.infn.it    
  // Get the pointer to the existing analysis manager via the static access method.                                                     
  //==============================================================================
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskDzeroKF", "No analysis manager to connect to.");
  }  

  TString taskName=("AliAnalysisTaskDmesonKF.cxx+");

  gROOT->LoadMacro(taskName.Data());
  if (gProof){
    TString taskSO=gSystem->pwd();
    taskSO+="/";
    taskSO+=taskName(0,taskName.First('.'))+"_cxx.so";
    gProof->Exec(Form("gSystem->Load(\"%s\")",taskSO.Data()),kTRUE);
  }
 
  Int_t meson = AliAnalysisTaskDmesonKF::kDzeroToKpi;
  Bool_t topoloosecuts=kTRUE;
  Bool_t readMC=kFALSE;
  Bool_t QAanalysis=kFALSE;
  Bool_t ITSrefit=kTRUE;
  Bool_t TPCrefit=kTRUE;
  Bool_t accKinks=kFALSE;
  Int_t minclsITS=2;
  Int_t minclsTPC=70;
  Double_t maxchiTPC=4;
  AliESDtrackCuts::Detector ITSreq1=AliESDtrackCuts::kSPD;
  AliESDtrackCuts::ITSClusterRequirement ITSreq2=AliESDtrackCuts::kAny;
  Bool_t SigToVtx=kFALSE;
  Double_t minDCA=0.;
  Int_t minpt;
  Int_t pdg=421; 
  Double_t mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  Double_t maxrapid=-999;
  Bool_t usePID=kTRUE;
  Bool_t useStrongPID=kTRUE;
  Double_t ptmaxForPID=2.;
  ULong64_t trmask= AliVEvent::kINT7;
  TString trclass[2] = {"",""};
  Int_t optpileup = 1;

  //topological cuts
  vector<Double_t> ptlims; //pt lims for different topological cuts values
  ptlims.push_back(0.);
  ptlims.push_back(2.);
  ptlims.push_back(5.);
  ptlims.push_back(1.e+6);
  vector<Double_t> mindaupt;
  mindaupt.push_back(0.2);
  mindaupt.push_back(0.3);
  mindaupt.push_back(0.3);
  vector<Double_t> cospmin;
  cospmin.push_back(0.9);
  cospmin.push_back(0.9);
  cospmin.push_back(0.9);
  vector<Double_t> decLmin;
  decLmin.push_back(0.);
  decLmin.push_back(0.);
  decLmin.push_back(0.);
  vector<Double_t> NdecLmin;
  NdecLmin.push_back(3.);
  NdecLmin.push_back(2.);
  NdecLmin.push_back(1.5);
  vector<Double_t> sigvtxmax;
  sigvtxmax.push_back(0.06);
  sigvtxmax.push_back(0.06);
  sigvtxmax.push_back(0.06);
  vector<Double_t> chimax;
  chimax.push_back(0.06);
  chimax.push_back(0.06);
  chimax.push_back(0.06);
  vector<Double_t> PVchimax;
  PVchimax.push_back(0.06);
  PVchimax.push_back(0.06);
  PVchimax.push_back(0.06);

  //Analysis Task
  AliAnalysisTaskDmesonKF *task = new AliAnalysisTaskDmesonKF("AnalysisTaskDmesonKF",meson,readMC);
  task->SetMeson(meson);
  task->SetReadMC(readMC);
  task->SetDoQAanalysis(QAanalysis);
  task->SetMinPt(minpt);
  task->SetMassWindow(mass-0.2,mass+0.2);
  task->SetTriggerMask(trmask);
  task->SetTriggerClass(trclass);
  task->SetMaxFiducialRapidity(maxrapid);
  task->SetUsePID(usePID);
  task->SetUseStrongPID(useStrongPID,ptmaxForPID);
  task->SetMinVtxContr(1);
  task->SetLooseCuts(topoloosecuts,ptlims,mindaupt,cospmin,decLmin,NdecLmin,sigvtxmax,chimax,PVchimax);
  task->SetOptPileup(optpileup);

  AliESDtrackCuts *cuts = new AliESDtrackCuts();
  cuts->SetRequireTPCRefit(TPCrefit);
  cuts->SetRequireITSRefit(ITSrefit);
  cuts->SetMinNClustersITS(minclsITS);
  cuts->SetClusterRequirementITS(ITSreq1,ITSreq2);
  cuts->SetMinNClustersTPC(minclsTPC);
  cuts->SetMaxChi2PerClusterTPC(maxchiTPC);   
  cuts->SetRequireSigmaToVertex(SigToVtx);
  cuts->SetMinDCAToVertexXY(minDCA);
  cuts->SetPtRange(0.2,1.e10);
  cuts->SetEtaRange(-0.8,0.8);

  task->SetESDtrackCuts(cuts);

  mgr->AddTask(task);

  //==============================================================================                                                        
  // Create containers for input/output 
  
  TString outname = "coutputDzeroKF";
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  TString outputfile = "fgrosa_Dzero_KF.root";
  AliAnalysisDataContainer *coutputDmesonKF = mgr->CreateContainer(outname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutputDmesonKF);

  return task;
}
