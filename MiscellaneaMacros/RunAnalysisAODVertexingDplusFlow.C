class AliAnalysisGrid;
class AliAnalysisAlien;

void RunAnalysisAODVertexingDplusFlow(Bool_t fUseMC = kFALSE, TString pluginmode="full")
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

    gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS/include -g");

    TString trainName = "D2H";
    TString analysisMode = "grid"; // "local", "grid", or "proof"
    TString inputMode    = "list"; // "list", "xml", or "dataset"
    Long64_t nentries=123567890,firstentry=0;
    Bool_t useParFiles=kFALSE;
    Bool_t useAlienPlugin=kTRUE;
    Bool_t saveProofToAlien=kFALSE;
    TString proofOutdir = "";
    TString loadMacroPath="$ALICE_PHYSICS/PWGHF/vertexingHF/macros/";
  
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
        TString loadLibraries="LoadLibraries.C"; loadLibraries.Prepend(loadMacroPath.Data());
        gROOT->LoadMacro(loadLibraries.Data());
        LoadLibraries(useParFiles);
    } else if (analysisMode=="proof") {
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
        AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,useParFiles,fUseMC);
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
    Bool_t flagForMC = kFALSE;
    if(fUseMC) flagForMC = kTRUE;
    AliAnalysisTaskSE *setupTask = AddTaskPIDResponse(flagForMC,flagForMC);//kTRUE su MC!!

    TString taskName;
    TString suffix="3050";
    if(fUseMC) suffix="3050_MC";
  
    TString cutfilename="alien:///alice/cern.ch/user/f/fgrosa/DplustoKpipiCuts_3050_central_topod0cut_kINT7.root";
    if(fUseMC) cutfilename="alien:///alice/cern.ch/user/f/fgrosa/DplustoKpipiCuts_3050_central_topod0cut_MC.root";
    TString suffix="_Central_Topomatic_d0_Cut";
    if(fUseMC) suffix="_Central_Topomatic_d0_Cut_MC";

    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/AddTaskFlowQnVectorCorrectionsToLegoTrainNewDetConfig.C");
    AliAnalysisDataContainer *corrTask = AddTaskFlowQnVectorCorrectionsToLegoTrainNewDetConfig("alien:///alice/cern.ch/user/f/fgrosa/QnConfig");

    gROOT->LoadMacro("/$ALICE_PHYSICS/PWGHF/vertexingHF/charmFlow/AddTaskHFv2.C");
    AliAnalysisTaskSEHFv2 *task1 = AddTaskHFv2(cutfilename.Data(),AliAnalysisTaskSEHFv2::kDplustoKpipi,"AnalysisCuts",kFALSE,suffix.Data(),0,30.,50.,kTRUE,AliAnalysisTaskSEHFv2::kEP,"QoverM");
    AliAnalysisTaskSEHFv2 *task2 = AddTaskHFv2(cutfilename.Data(),AliAnalysisTaskSEHFv2::kDplustoKpipi,"AnalysisCuts",kFALSE,suffix.Data(),0,30.,50.,kTRUE,AliAnalysisTaskSEHFv2::kSP,"QoverM");
    AliAnalysisTaskSEHFv2 *task3 = AddTaskHFv2(cutfilename.Data(),AliAnalysisTaskSEHFv2::kDplustoKpipi,"AnalysisCuts",kFALSE,suffix.Data(),0,30.,50.,kTRUE,AliAnalysisTaskSEHFv2::kEvShape,"QoverM",AliAnalysisTaskSEHFv2::kq2TPC);
  
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskCleanupVertexingHF.C");
    AliAnalysisTaskSECleanupVertexingHF *taskclean = AddTaskCleanupVertexingHF();

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
AliAnalysisGrid* CreateAlienHandler(TString pluginmode="test",Bool_t useParFiles=kFALSE, Bool_t fUseMC)
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
    plugin->SetAliPhysicsVersion("vAN-20161109-1");
    plugin->SetNtestFiles(1);
    plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS/include -g");
    // gROOT->LoadMacro("AddGoodRunsSPLIT.C");
    gROOT->LoadMacro("$ALIPHYSICS/PWGHF/vertexingHF/AddGoodRuns.C");

    ///PER MC
    if(fUseMC==kTRUE) {
        //to be setted!
        plugin->SetGridDataDir("/alice/sim/2016/LHC16i2b");
        plugin->SetDataPattern("/AOD/*AliAOD.root");
        plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
        Int_t totruns=0;
        plugin->AddRunNumber(245833);
        plugin->SetNrunsPerMaster(1);
    }
    else {
        plugin->SetGridDataDir("/alice/data/2015/LHC15o/");
        plugin->SetDataPattern("/pass1/AOD/*AliAOD.root");
        plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
        plugin->SetRunPrefix("000"); //data
        plugin->AddRunNumber(246087); //very big run pass1
        //plugin->AddRunNumber(246982); //small pass1
        //plugin->AddRunNumber(244917);//low intensity pass2_lowIR
        plugin->SetNrunsPerMaster(1);
    }
  
    TString suffix="3050";
    if(fUseMC) suffix="3050_MC";
  
    plugin->SetGridWorkingDir(Form("Dplusv2_PbPb_%s",suffix.Data()));
    // Name of executable
    plugin->SetExecutable(Form("Dplusv2_PbPb_%s.sh",suffix.Data()));
    // Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
    // Declare the analysis source files names separated by blancs. To be compiled runtime
    // using ACLiC on the worker nodes.
    // Declare all libraries (other than the default ones for the framework. These will be
    // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    plugin->SetAdditionalLibs("libPWGflowBase.so libPWGflowTasks.so libPWGHFbase.so libPWGHFvertexingHF.so");

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
    plugin->SetMaxMergeStages(3);
    plugin->SetSplitMaxInputFileNumber(5);
    // Optionally set time to live (default 30000 sec)
    // plugin->SetTTL(30000);

    // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    plugin->SetAnalysisMacro("AnalysisHF.C");
    // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
    plugin->SetMaxMergeFiles(20);
    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName("TaskHF.jdl");

    return plugin;
}
