class AliAnalysisGrid;
class AliAnalysisAlien;


void RunHFMCCheck(TString runmode="terminate"){

  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  
  gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS  -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/CDB -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS/include -g");
//  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
//  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  // Load common libraries
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");   
  gSystem->Load("libCORRFW.so");
  
  gROOT->LoadMacro("AliAnalysisTaskCheckHFMCProd2.cxx++g");   
  gROOT->LoadMacro("AddHFMCCheck.C");   

  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

  AliAnalysisGrid *alienHandler = 0x0;
  TChain *chainESD= 0x0;
  // Create and configure the alien handler plugin
  if(runmode=="local"){
    chainESD= new TChain("esdTree");
    chainESD->Add("AliESDs.root");
  }else{
    alienHandler = CreateAlienHandler(runmode);  
    if (!alienHandler) return;	
    mgr->SetGridHandler(alienHandler);
  }

  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  // Apply the event selection
  AliAnalysisTaskCheckHFMCProd2 *task = AddHFMCCheck(2);

  mgr->InitAnalysis();
  mgr->PrintStatus();
  // Start analysis in grid.
  if(runmode=="local"){
    mgr->StartAnalysis("local",chainESD);
  }else{
    mgr->StartAnalysis("grid");
  }
};


AliAnalysisGrid* CreateAlienHandler(TString runmode) {
  // Check if user has a valid token, otherwise make one. This has limitations
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
  //if (!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(runmode.Data());
  plugin->SetNtestFiles(1);
  plugin->SetUser("fgrosa");
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-30-alice-12");
  plugin->SetAliROOTVersion("v5-08-01-1");
  plugin->SetAliPhysicsVersion("vAN-20160310-1");
  // Declare input data to be processed.
  // Method 1: Create automatically XML collections using alien 'find' command
  // Define production directory LFN
  //plugin->SetGridDataDir("/alice/sim/2012/LHC12a17a_fix");
  //plugin->SetGridDataDir("/alice/data/2011/LHC11a");
  plugin->SetGridDataDir("/alice/sim/2013/LHC13d3"); // specify LHC period
  // Set data search pattern
  plugin->SetDataPattern("*AliESDs.root");
  //  plugin->SetDataPattern("pass2_without_SDD/*AliESDs.root");
  // ...then add run numbers to be considered
  //  plugin->SetRunPrefix("000");
  
  //  plugin->AddRunNumber(167915);
  //  plugin->AddRunNumber(167920);

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

  
  // Method 2: Declare existing data files (raw collections, xml collections, root file)
  // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
  // XML collections added via this method can be combined with the first method if
  // the content is compatible (using or not tags)
  //plugin->AddDataFile("alice/cern.ch/user/e/ebiolcat/work/000104073.xml");
  //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir("HFMCcheck_pPb_ptB");
  plugin->SetExecutable("HFMCcheck_pPb_ptB.sh");
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  plugin->SetAnalysisSource("AliAnalysisTaskCheckHFMCProd2.cxx"); 
  gROOT->ProcessLine(".L AliAnalysisTaskCheckHFMCProd2.cxx+g");
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS  -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/CDB -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS/include -g");
                         
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalLibs("libPWGflowBase.so libPWGflowTasks.so libPWGHFbase.so libPWGHFvertexingHF.so AliAnalysisTaskCheckHFMCProd2.h AliAnalysisTaskCheckHFMCProd2.cxx");// "libGui.so libProof.so libMinuit.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so");
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  plugin->SetDefaultOutputs(kFALSE);
  plugin->SetOutputFiles("AnalysisResults.root");
  // Optionally define the files to be archived.
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
  plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  //  plugin->SetAnalysisMacro("analysisITSsaTracks.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
//  plugin->SetSplitMaxInputFileNumber(100);
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//  plugin->SetMaxInitFailed(5);
  // Optionally resubmit threshold.
//  plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
//  plugin->SetTTL(20000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("addtaskHFMCCheck.jdl");
  // Optionally modify job price (default 1)
//  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");
  return plugin;
}
