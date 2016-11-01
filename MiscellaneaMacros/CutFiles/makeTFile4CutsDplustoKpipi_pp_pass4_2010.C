#include <Riostream.h>
#include <TFile.h>
//#include <AliRDHFCutsDplustoKpipi.h>
#include <TClonesArray.h>
#include <TParameter.h>

//macro to make a .root file which contains an AliRDHFCutsDplustoKpipi with loose set of cuts (for significance maximization) and TParameter with the tighest value of these cuts
//Needed for AliAnalysisTaskSEDplus, AliCFTaskVertexingHF3Prong, AliAnalysisTaskSESignificance

//Use:
//Set hard coded commented with //set this!!

//.L makeTFile4CutsDplustoKpipi.C
// makeInputAliAnalysisTaskSEDplus()
// makeInputAliAnalysisTaskSESignificanceMaximization()


void makeTFile4CutsDplustoKpipi_pp(){
    
    //  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWG3 -I$ALICE_ROOT/PWG3/vertexingHF -I$ALICE_ROOT/PWG3/vertexingH/macros -g");
    
    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    //default
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    //esdTrackCuts->SetMinNClustersITS(4); // default is 5
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                           AliESDtrackCuts::kAny);
    // default is kBoth, otherwise kAny
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3,1.e10);
    
    TString cent="";
    //centrality selection (Pb-Pb)
    //Float_t minc=0.,maxc=20;
    const Int_t nptbins=16;
    Float_t* ptbins;
    ptbins=new Float_t[nptbins+1];

    ptbins[0]=1.;
    ptbins[1]=2.;
    ptbins[2]=3.;
    ptbins[3]=4.;
    ptbins[4]=5.;
    ptbins[5]=6.;
    ptbins[6]=7.;
    ptbins[7]=8.;
    ptbins[8]=9.;
    ptbins[9]=10.;
    ptbins[10]=11.;
    ptbins[11]=12.;
    ptbins[12]=14.;
    ptbins[13]=16.;
    ptbins[14]=24.;
    ptbins[15]=36.;
       
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
        anacutsval[1][ipt]=0.3;
    }
    anacutsval[1][0]=0.4;
    anacutsval[1][1]=0.4;
    anacutsval[1][2]=0.4;
//    anacutsval[1][2]=0.9;
//    anacutsval[1][3]=0.6;
//    anacutsval[1][4]=0.6;
//    anacutsval[1][5]=0.4;
//    anacutsval[1][6]=0.4;
//    anacutsval[1][7]=0.4;
//    anacutsval[1][8]=0.4;
//    anacutsval[1][9]=0.4;
//    anacutsval[1][10]=0.4;
    
    ic=2;//ptPi
    for(Int_t ipt=0;ipt<nptbins;ipt++){
        anacutsval[ic][ipt]=0.3;
    }
    anacutsval[1][0]=0.4;
    anacutsval[1][1]=0.4;
    anacutsval[1][2]=0.4;
//    anacutsval[2][0]=0.6;
//    anacutsval[2][1]=0.73;
//    anacutsval[2][2]=0.8;
//    anacutsval[2][3]=0.6;
//    anacutsval[2][4]=0.6;
//    
//    anacutsval[2][5]=0.4;//8-9
//    anacutsval[2][6]=0.4;//9 - 10
//    anacutsval[2][7]=0.4;//10-11
//    anacutsval[2][8]=0.4;//11 - 12
//    anacutsval[2][9]=0.4;//12 - 14
//    anacutsval[2][10]=0.4;//14-16
    
    
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
        anacutsval[ic][ipt]=0.07;
    }
    anacutsval[7][3]=0.1;
    anacutsval[7][4]=0.1;
    anacutsval[7][5]=0.1;
    anacutsval[7][6]=0.1;
    anacutsval[7][7]=0.1;
    anacutsval[7][8]=0.12;
    anacutsval[7][9]=0.12;
    anacutsval[7][10]=0.12;
    anacutsval[7][11]=0.12;
    anacutsval[7][12]=0.12;
    anacutsval[7][13]=0.12;
    anacutsval[7][14]=0.12;
    anacutsval[7][15]=0.2;

    ic=8;//pM
    for(Int_t ipt=0;ipt<nptbins;ipt++){
        anacutsval[ic][ipt]=0.0;
    }
    
    
    //cosp
    ic=9;
    for(Int_t ipt=2;ipt<nptbins;ipt++){
        anacutsval[ic][ipt]=0.97;
    }
    anacutsval[9][14]=0.98;
    anacutsval[9][15]=0.98;
    
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
    anacutsval[12][0]=6.;
    anacutsval[12][15]=8.;
   
    ic=13;//cospXY
    for(Int_t ipt=0;ipt<nptbins;ipt++){
        anacutsval[ic][ipt]=0.;
    }
    anacutsval[13][0]=0.985;
    anacutsval[13][1]=0.985;
    anacutsval[13][2]=0.98;  
    
    //topomatic variable (max norm d0meas-d0exp)
    Float_t topocutvalues[nptbins];
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      topocutvalues[ipt]=2.5;
    } 

    AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
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
    //analysiscuts->SetMinCentrality(minc);
    //analysiscuts->SetMaxCentrality(maxc);
    
    analysiscuts->SetRemoveTrackletOutliers(kTRUE);//added on June 28
    analysiscuts->SetCutOnzVertexSPD(0);//added on July19
    
    // cent=Form("%.0f%.0f",minc,maxc);
    // analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    // analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    // analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    // analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    
    // analysiscuts->EnableSemiCentralTrigger();
    analysiscuts->SetMinPtCandidate(1.);
    analysiscuts->SetMaxPtCandidate(36.);
    analysiscuts->Setd0MeasMinusExpCut(nptbins,topocutvalues);

    cout<<"This is the object I'm going to save:"<<nptbins<<endl;
    
    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TFile* fout=new TFile("DplustoKpipiCuts_pp_pass4.root","recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();
    
}
