class AliAnalysisGrid;
class AliAnalysisAlien;

void RunAnalysisAODVertexingHFCFComp()
{
  //
  // Test macro for AliAnalysisTaskSE's for heavy-flavour candidates
  // It has the structure of a Analysis Train:
  // - in this macro, change things related to running mode
  //   and input preparation 
  // - add your task using a AddTaskXXX macro 
  //
  // A.Dainese, andrea.dainese@lnl.infn.it
  // "grid" mode added by R.Bala, bala@to.infn.it
  //
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libOADB.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGHFbase.so");
  gSystem->Load("libPWGflowBase.so");
  gSystem->Load("libPWGflowTasks.so");
  gSystem->Load("libPWGHFvertexingHF.so");
  
  TString trainName = "D2H";
  TString analysisMode = "grid"; // "local", "grid", or "proof"
  TString inputMode    = "list"; // "list", "xml", or "dataset"
  Long64_t nentries=123567890,firstentry=0;
  Bool_t useParFiles=kFALSE;
  Bool_t useAlienPlugin=kTRUE;
  TString pluginmode="terminate";// "test" "full" "terminate(merge dei file)"
  Bool_t saveProofToAlien=kFALSE;
  TString proofOutdir = "";
  TString loadMacroPath="$ALICE_PHYSICS/../src/PWGHF/vertexingHF/macros/";
  //TString loadMacroPath="./"; // this is normally needed for CAF
  //

  if(analysisMode=="grid") {
    // Connect to AliEn
    TGrid::Connect("alien://");
  } else if(analysisMode=="proof") {
    // Connect to the PROOF cluster
    if(inputMode!="dataset") {printf("Input mode must be dataset, for proof analysis\n"); return;}
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    TProof::Open("alicecaf");
    //TProof::Reset("alicecaf");
    if(saveProofToAlien) {
      TGrid::Connect("alien://");
      if(gGrid) {
	TString homedir = gGrid->GetHomeDirectory();
	TString workdir = homedir + trainName;
	if(!gGrid->Cd(workdir)) {
	  gGrid->Cd(homedir);
	  if(gGrid->Mkdir(workdir)) {
	    gGrid->Cd(trainName);
	    ::Info("VertexingTrain::Connect()", "Directory %s created", gGrid->Pwd());
	  }
	}	   
	gGrid->Mkdir("proof_output");
	gGrid->Cd("proof_output");
	proofOutdir = Form("alien://%s", gGrid->Pwd());
      } 
    }
  }
  
  // AliRoot libraries
  if(analysisMode=="local" || analysisMode=="grid") {
    TString loadLibraries="LoadLibraries.C";
    loadLibraries.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(loadLibraries.Data());
    LoadLibraries(useParFiles);
  }
  else if (analysisMode=="proof") {
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libVMC.so");    
    gSystem->Load("libMinuit.so");    
    // Enable the needed packages
    //gProof->ClearPackages();
    TString parDir="/afs/cern.ch/user/d/dainesea/code/";
    TString parFile;
    if(!useParFiles) {
      gProof->UploadPackage("AF-v4-17");
      gProof->EnablePackage("AF-v4-17");
      // --- Enable the PWGHFvertexingHF Package
      parFile="PWGHFvertexingHF.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("PWGHFvertexingHF");
    } else {
      // --- Enable the STEERBase Package
      parFile="STEERBase.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("STEERBase");
      // --- Enable the ESD Package
      parFile="ESD.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("ESD");
      // --- Enable the AOD Package
      parFile="AOD.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("AOD");
      // --- Enable the ANALYSIS Package
      parFile="ANALYSIS.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("ANALYSIS");
      // --- Enable the ANALYSISalice Package
      parFile="ANALYSISalice.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("ANALYSISalice");
      // --- Enable the CORRFW Package
      parFile="CORRFW.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("CORRFW");
      // --- Enable the PWGHFbase Package
      parFile="PWGHFbase.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("PWGHFbase");
      // --- Enable the PWGHFvertexingHF Package
      parFile="PWGHFvertexingHF.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("PWGHFvertexingHF");
    }
    gProof->ShowEnabledPackages(); // show a list of enabled packages
  }
  
  // Create Alien plugin, if requested
  if(useAlienPlugin) {  
    if(analysisMode!="grid") {printf("Analysis mode must be grid, to use alien plugin\n"); return;}
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,useParFiles);  
    if(!alienHandler) return;
  }

  //-------------------------------------------------------------------
  // Prepare input
  TChain *chainAOD = 0;
  TString dataset; // for proof

  if(!useAlienPlugin) {
    TString makeAODInputChain="../MakeAODInputChain.C"; makeAODInputChain.Prepend(loadMacroPath.Data());
    if(inputMode=="list") {
      // Local files
      gROOT->LoadMacro(makeAODInputChain.Data());
      chainAOD = MakeAODInputChain();// with this it reads ./AliAOD.root and ./AliAOD.VertexingHF.root
      //chainAOD = MakeAODInputChain("alien:///alice/cern.ch/user/r/rbala/newtrain/out_lhc08x/180100/",1,1);
      printf("ENTRIES %d\n",chainAOD->GetEntries());
    } else if(inputMode=="xml") {
      // xml
      gROOT->LoadMacro(makeAODInputChain.Data());
      chainAOD = MakeAODInputChain("collection_aod.xml","collection_aodHF.xml");
    } else if(inputMode=="dataset") {
      // CAF dataset
      //gProof->ShowDataSets();
      dataset="/ITS/dainesea/AODVertexingHF_LHC08x_180100";
    }
  }

  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
  // Connect plug-in to the analysis manager
  if(useAlienPlugin) mgr->SetGridHandler(alienHandler);

  // Input
  AliAODInputHandler *inputHandler = new AliAODInputHandler("handler","handler for D2H");
  if(analysisMode=="proof" ) {
    inputHandler->AddFriend("./AliAOD.VertexingHF.root");
    //inputHandler->AddFriend("deltas/AliAOD.VertexingHF.root");
    if(saveProofToAlien) mgr->SetSpecialOutputLocation(proofOutdir);
  }
  mgr->SetInputEventHandler(inputHandler);

  //-------------------------------------------------------------------

  
  //-------------------------------------------------------------------
  // Analysis tasks (wagons of the train)   
  //
  // First add the task for the PID response setting
  
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskSE *setupTask = AddTaskPIDResponse(kTRUE,kTRUE,kTRUE,2); //kTRUE,kTRUE per MC, oppure kFALSE,kTRUE per dati

  TString taskName;

  taskName="AddTaskCFVertexingHF3Prong.C";
  taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliCFTaskVertexingHF *task1 = AddTaskCFVertexingHF3Prong("ImpParpPbMC_CF","alien:///alice/cern.ch/user/f/fgrosa/DplustoKpipiCutspPbMC_ImpPar.root",AliCFTaskVertexingHF::kCheetah,kFALSE,kFALSE);
  AliCFTaskVertexingHF *task2 = AddTaskCFVertexingHF3Prong("ImpParpPbMC_CF","alien:///alice/cern.ch/user/f/fgrosa/DplustoKpipiCutspPbMC_ImpPar.root",AliCFTaskVertexingHF::kCheetah,kTRUE,kTRUE);

  taskName="AddTaskDplus.C";
  taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliAnalysisTaskSEDplus *dplustask = AddTaskDplus(1,0,100,2,kTRUE,kTRUE,"_ImpParpPbMC","alien:///alice/cern.ch/user/f/fgrosa/DplustoKpipiCutspPbMC_ImpPar.root","AnalysisCuts",0);//ultimo Bool_t -> readMC (kFALSE per dati)
  //-------------------------------------------------------------------
  
  //
  // Run the analysis
  //    
  if(chainAOD) printf("CHAIN HAS %d ENTRIES\n",(Int_t)chainAOD->GetEntries());
  
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(analysisMode=="grid" && !useAlienPlugin) analysisMode="local";
  if(analysisMode!="proof") {
    mgr->StartAnalysis(analysisMode.Data(),chainAOD,nentries,firstentry);
  } else {
    // proof
    mgr->StartAnalysis(analysisMode.Data(),dataset.Data(),nentries,firstentry);
  }
  
  return;
}
//_____________________________________________________________________________
//
AliAnalysisGrid* CreateAlienHandler(TString pluginmode="test",Bool_t useParFiles=kFALSE)
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(pluginmode.Data());
   plugin->SetUser("fgrosa");
   // Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-34-30-alice-3");
   plugin->SetAliROOTVersion("v5-07-02-1");
   plugin->SetAliPhysicsVersion("vAN-20151022-1");
   plugin->SetNtestFiles(1);
   gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGHF/vertexingHF/AddGoodRuns.C");

   // Declare input data to be processed.
   //************************************************
   // Set data search pattern for DATA
   //************************************************  
   //Method 1: To create automatically xml through plugin  
   // plugin->SetGridDataDir("/alice/data/2010/LHC10d"); // specify LHC period
   //plugin->SetDataPattern("pass2/AOD018/*AliAOD.root"); // specify reco pass and AOD set

   //PbPb
   // plugin->SetGridDataDir("/alice/sim/2012/LHC12a17a"); // specify LHC period
   // plugin->SetDataPattern("AOD110/*AliAOD.root"); // specify reco pass and AOD set

   //pp
   // plugin->SetGridDataDir("/alice/sim/LHC10f7a"); // specify LHC period
   // plugin->SetDataPattern("AOD136a/*AliAOD.root"); // specify reco pass and AOD set

   //pPb
   plugin->SetGridDataDir("/alice/sim/2013/LHC13d3"); // specify LHC period
   plugin->SetDataPattern("*/AliAOD.root"); // specify reco pass and AOD set

   
   plugin->SetFriendChainName("AliAOD.VertexingHF.root");

   // OR plugin->SetFriendChainName("deltas/AliAOD.VertexingHF.root");
   // Adds only the good runs from the Monalisa Run Condition Table
   // More than one period can be added but the period name has to be removed from GridDataDir (to be tested)
   Int_t totruns=0;

   //  totruns += AddGoodRuns(plugin,"LHC11h_pass2","LHC12a17a"); // specify LHC period (PbPb)
 
   //   totruns += AddGoodRuns(plugin,"LHC10c","LHC10f7a"); // specify LHC period (pp)

   // plugin->AddRunNumber(129744);

   //totruns += AddGoodRuns(plugin,"LHC10c"); // specify LHC period
   //   totruns += AddGoodRuns(plugin,"LHC11h_pass2","LHC12a17a"); // specify LHC period
   //   plugin->AddRunNumber(170572);

    //(pPb)
    plugin->AddRunNumber(195389);
    plugin->AddRunNumber(195390);
    plugin->AddRunNumber(195391);
    plugin->AddRunNumber(195478);
    plugin->AddRunNumber(195479);
    plugin->AddRunNumber(195480);
    plugin->AddRunNumber(195481);
    plugin->AddRunNumber(195482);
    plugin->AddRunNumber(195596);
    plugin->AddRunNumber(195673);
    plugin->AddRunNumber(195675);
    plugin->AddRunNumber(195677);
    
    plugin->SetNrunsPerMaster(12);
    
    //plugin->SetGridDataDir("/alice/sim/LHC10d4"); // specify LHC period
    //plugin->SetDataPattern("AOD056/*AliAOD.root"); // specify reco pass and AOD set
    
    plugin->SetFriendChainName("AliAOD.VertexingHF.root");
        
   // Method 2: Declare existing data files (e.g xml collections)

   //plugin->AddDataFile("/alice/cern.ch/user/r/rbala/000168068_000170593.xml");
   //  plugin->SetDataPattern("*AliAOD.root");
   //  plugin->SetFriendChainName("./AliAOD.VertexingHF.root"); 

   //************************************************
   // Set data search pattern for MONTECARLO
   //************************************************
   /* 
   plugin->SetGridDataDir("/alice/sim/LHC10d3"); // specify MC sample
   plugin->SetDataPattern("AOD005/*AliAOD.root"); // specify AOD set
   plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
   // OR plugin->SetFriendChainName("deltas/AliAOD.VertexingHF.root");
   // Adds only the good runs from the Monalisa Run Condition Table 
   // More than one period can be added!
   Int_t totruns=0;
   totruns += AddGoodRuns(plugin,"LHC10b","LHC10d3"); // specify LHC period for anchor runs; and the name of the MC production
   //totruns += AddGoodRuns(plugin,"LHC10c","LHC10f7"); // specify LHC period for anchor runs;  and the name of the MC production
   //totruns += AddGoodRuns(plugin,"LHC10d","LHC10f7"); // specify LHC period for anchor runs;  and the name of the MC production
   plugin->SetNrunsPerMaster(totruns);
   */
   //
   // Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir("ImpParpPbMC");//Dplus_Effic_LHC10d4
   // plugin->SetGridWorkingDir("Dplus_Effic_LHC12a17a");

   // Name of executable

   plugin->SetExecutable("ImpParpPbMC.sh");
   //  plugin->SetExecutable("Dplus_Effic_LHC12a17a.sh");

   // Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
   // Declare the analysis source files names separated by blancs. To be compiled runtime
   // using ACLiC on the worker nodes.
   //plugin->SetAnalysisSource("AliAnalysisTaskSEDplus2.cxx");
   // Declare all libraries (other than the default ones for the framework. These will be
   // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   //gROOT->ProcessLine(".L AliAnalysisTaskSEDplus2.cxx+g");

   plugin->SetAdditionalLibs("libSTEERBase.so libESD.so libPWGflowBase.so libPWGflowTasks.so libPWGHFbase.so libPWGHFvertexingHF.so libOADB.so");
   
   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS/include -g");
   //   plugin->SetAdditionalLibs("libPWGflowBase.so libPWGflowTasks.so libPWGHFbase.so libPWGHFvertexingHF.so AliAnalysisTaskSEDplus2.h AliAnalysisTaskSEDplus2.cxx");

  // use par files
   if(useParFiles) {
     plugin->EnablePackage("STEERBase.par");
     plugin->EnablePackage("ESD.par");
     plugin->EnablePackage("AOD.par");
     plugin->EnablePackage("ANALYSIS.par");
     plugin->EnablePackage("OADB.par");
     plugin->EnablePackage("ANALYSISalice.par");
     plugin->EnablePackage("CORRFW.par");
     plugin->EnablePackage("PWGHFbase.par");
     plugin->EnablePackage("PWGHFvertexingHF.par");
   }

    plugin->SetDefaultOutputs(kTRUE);
   // merging via jdl
   plugin->SetMergeViaJDL(kFALSE);
   plugin->SetOneStageMerging(kFALSE);
   plugin->SetMaxMergeStages(2);

   // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("AnalysisHFCuts.C");
   // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   // Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("TaskHFCuts.jdl");
   
   return plugin;
}
