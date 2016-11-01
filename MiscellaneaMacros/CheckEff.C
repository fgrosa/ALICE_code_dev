void CheckEff(Bool_t usePID=kTRUE, Int_t PIDaxis=7) {

  TString fileName="AnalysisResults.root";

  TString baseDirFilName="PWG3_D2H_CFtaskDplustoKpipi_PromptImpParpPbMC_CF";
  TString baseContainerName="CFHFcontainer_DplustoKpipi_PromptImpParpPbMC_CF";
  TString baseDirFilNameB="PWG3_D2H_CFtaskDplustoKpipi_FromBImpParpPbMC_CF";
  TString baseContainerNameB="CFHFcontainer_DplustoKpipi_FromBImpParpPbMC_CF";
  
  TFile* file = new TFile(fileName.Data());

  TDirectoryFile* df=(TDirectoryFile*)file->Get(Form("%s",baseDirFilName.Data()));
  AliCFContainer* cont=(AliCFContainer*)df->Get(Form("%s",baseContainerName.Data()));
  
  TDirectoryFile* dfB=(TDirectoryFile*)file->Get(Form("%s",baseDirFilNameB.Data()));
  AliCFContainer* contB=(AliCFContainer*)dfB->Get(Form("%s",baseContainerNameB.Data()));

  Int_t MCGenAccStep = 2;
  Int_t RecoStep = 8;
  Int_t RecoPIDStep = 9;
  
  TH1F *hptMCAccPromptCF = (TH1F*)cont->Project(MCGenAccStep,0);
  TH1F *hptMCAccBFeedCF = (TH1F*)contB->Project(MCGenAccStep,0);
  TH1F *hptRecoCutsPromptCF = 0x0;
  TH1F *hptRecoCutsBFeedCF = 0x0;
  
  if(usePID) {
    hptRecoCutsPromptCF = (TH1F*)cont->Project(RecoPIDStep,0);
    hptRecoCutsBFeedCF = (TH1F*)contB->Project(RecoPIDStep,0);
  }
  else {
    hptRecoCutsPromptCF = (TH1F*)cont->Project(RecoStep,0);
    hptRecoCutsBFeedCF = (TH1F*)contB->Project(RecoStep,0);
  }
    
  hptMCAccPromptCF->SetDirectory(0);
  hptMCAccBFeedCF->SetDirectory(0);
  hptRecoCutsPromptCF->SetDirectory(0);
  hptRecoCutsBFeedCF->SetDirectory(0);
  
  TDirectoryFile *DplusDir = (TDirectoryFile*)file->Get("PWG3_D2H_InvMassDplus");
  TList *List = (TList*)DplusDir->Get("coutputDplus_ImpParpPbMC0100");
  THnSparse *hMCAccPrompt = (THnSparse*)List->FindObject("hMCAccPrompt");
  THnSparse *hMCAccBFeed = (THnSparse*)List->FindObject("hMCAccBFeed");
  THnSparse *hSparsePrompt= (THnSparse*)List->FindObject("hMassPtImpParPrompt");
  THnSparse *hSparseBFeed= (THnSparse*)List->FindObject("hMassPtImpParBfeed");

  for(Int_t iAxis=0; iAxis<hSparsePrompt->GetNdimensions(); iAxis++) {
    hSparsePrompt->GetAxis(iAxis)->SetRange(-1,-1);
    hSparseBFeed->GetAxis(iAxis)->SetRange(-1,-1);
  }
  /*
  if(usePID) {
    hSparsePrompt->GetAxis(PIDaxis)->SetRange(3,3);
    hSparseBFeed->GetAxis(PIDaxis)->SetRange(3,3);
  }
  */
  TH1F *hptMCAccPromptCuts = (TH1F*)hMCAccPrompt->Projection(0);
  TH1F *hptMCAccBFeedCuts = (TH1F*)hMCAccBFeed->Projection(0);

  TH1F *hptRecoCutsPromptCuts = (TH1F*)hSparsePrompt->Projection(1);
  TH1F *hptRecoCutsBFeedCuts = (TH1F*)hSparseBFeed->Projection(1);

  hptMCAccPromptCuts->SetDirectory(0);
  hptMCAccBFeedCuts->SetDirectory(0);
  hptRecoCutsPromptCuts->SetDirectory(0);
  hptRecoCutsBFeedCuts->SetDirectory(0);

  file->Close();

  //rebin for same binning

  const Int_t nPtBins = 11;
  const Int_t nPtLims = nPtBins+1;

  Double_t ptBins[nPtLims] = {1,2,3,4,5,6,7,8,12,16,24,36};

  //CF rebinned
  TH1F* hptMCAccPromptCFReb = (TH1F*)hptMCAccPromptCF->Rebin(nPtBins,"hptMCAccPromptCFReb",ptBins);
  TH1F *hptRecoCutsPromptCFReb = (TH1F*)hptRecoCutsPromptCF->Rebin(nPtBins,"hptRecoCutsPromptCFReb",ptBins);
  TH1F* hptMCAccBFeedCFReb = (TH1F*)hptMCAccBFeedCF->Rebin(nPtBins,"hptMCAccBFeedCFReb",ptBins);
  TH1F *hptRecoCutsBFeedCFReb = (TH1F*)hptRecoCutsBFeedCF->Rebin(nPtBins,"hptRecoCutsBFeedCFReb",ptBins);

  hptMCAccPromptCFReb->Sumw2();
  hptRecoCutsPromptCFReb->Sumw2();
  hptMCAccBFeedCFReb->Sumw2();
  hptRecoCutsBFeedCFReb->Sumw2();
  
  //CutsTask rebinned
  TH1F* hptMCAccPromptCutsReb = (TH1F*)hptMCAccPromptCuts->Rebin(nPtBins,"hptMCAccPromptCutsReb",ptBins);
  TH1F *hptRecoCutsPromptCutsReb = (TH1F*)hptRecoCutsPromptCuts->Rebin(nPtBins,"hptRecoCutsPromptCutsReb",ptBins);
  TH1F* hptMCAccBFeedCutsReb = (TH1F*)hptMCAccBFeedCuts->Rebin(nPtBins,"hptMCAccBFeedCutsReb",ptBins);
  TH1F *hptRecoCutsBFeedCutsReb = (TH1F*)hptRecoCutsBFeedCuts->Rebin(nPtBins,"hptRecoCutsBFeedCutsReb",ptBins);

  hptMCAccPromptCutsReb->Sumw2();
  hptRecoCutsPromptCutsReb->Sumw2();
  hptMCAccBFeedCutsReb->Sumw2();
  hptRecoCutsBFeedCutsReb->Sumw2();
  
  //Efficiency for Prompt and FeedDown with CF
  TH1F *hEffPromptCF =new TH1F("hEffPromptCF","kStepRecoCuts/kStepMCAcc Prompt",nPtBins,ptBins);
  hEffPromptCF->Divide(hptRecoCutsPromptCFReb,hptMCAccPromptCFReb,1,1,"B");

  TH1F *hEffBFeedCF =new TH1F("hEffBFeedCF","kStepRecoCuts/kStepMCAcc FeedDown",nPtBins,ptBins);
  hEffBFeedCF->Divide(hptRecoCutsBFeedCFReb,hptMCAccBFeedCFReb,1,1,"B");

  //Efficiency for Prompt and FeedDown with CutsTask
  TH1F *hEffPromptCuts =new TH1F("hEffPromptCuts","kStepRecoCuts/kStepMCAcc Prompt",nPtBins,ptBins);
  hEffPromptCuts->Divide(hptRecoCutsPromptCutsReb,hptMCAccPromptCutsReb,1,1,"B");
  hEffPromptCuts->SetLineColor(kRed);
  
  TH1F *hEffBFeedCuts =new TH1F("hEffBFeedCuts","kStepRecoCuts/kStepMCAcc FeedDown",nPtBins,ptBins);
  hEffBFeedCuts->Divide(hptRecoCutsBFeedCutsReb,hptMCAccBFeedCutsReb,1,1,"B");
  hEffBFeedCuts->SetLineColor(kRed);

  TCanvas *cAllReb = new TCanvas("cAllReb","",10,10,1200,900);
  cAllReb->Clear();
  cAllReb->Divide(4,2);

  cAllReb->cd(1);
  hptMCAccPromptCFReb->Draw();
 
  cAllReb->cd(2);
  hptRecoCutsPromptCFReb->Draw();

  cAllReb->cd(3);
  hptMCAccBFeedCFReb->Draw();

  cAllReb->cd(4);
  hptRecoCutsBFeedCFReb->Draw();
  
  cAllReb->cd(5);
  hptMCAccPromptCutsReb->Draw("E0");
  hptMCAccPromptCutsReb->SetLineColor(kRed);
  
  cAllReb->cd(6);
  hptRecoCutsPromptCutsReb->Draw("E0");
  hptRecoCutsPromptCutsReb->SetLineColor(kRed);
  
  cAllReb->cd(7);
  hptMCAccBFeedCutsReb->Draw("E0");
  hptMCAccBFeedCutsReb->SetLineColor(kRed);
  
  cAllReb->cd(8);
  hptRecoCutsBFeedCutsReb->Draw("E0");
  hptRecoCutsBFeedCutsReb->SetLineColor(kRed);

  TCanvas *cAll = new TCanvas("cAll","",10,10,1200,900);
  cAll->Clear();
  cAll->Divide(4,2);

  cAll->cd(1);
  hptMCAccPromptCF->Draw();
 
  cAll->cd(2);
  hptRecoCutsPromptCF->Draw();

  cAll->cd(3);
  hptMCAccBFeedCF->Draw();

  cAll->cd(4);
  hptRecoCutsBFeedCF->Draw();
  
  cAll->cd(5);
  hptMCAccPromptCuts->Draw("E0");
  hptMCAccPromptCuts->SetLineColor(kRed);
  
  cAll->cd(6);
  hptRecoCutsPromptCuts->Draw("E0");
  hptRecoCutsPromptCuts->SetLineColor(kRed);
  
  cAll->cd(7);
  hptMCAccBFeedCuts->Draw("E0");
  hptMCAccBFeedCuts->SetLineColor(kRed);
  
  cAll->cd(8);
  hptRecoCutsBFeedCuts->Draw("E0");
  hptRecoCutsBFeedCuts->SetLineColor(kRed);
  
  TCanvas *cComparison = new TCanvas("cComparison","",10,10,1200,900);
  cComparison->Clear();
  cComparison->Divide(2,1);

  TLegend *legend = new TLegend(0.18,0.68,0.47,0.87);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->SetTextSizePixels(1000);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(hEffPromptCF,"Correction Framework","l");
  legend->AddEntry(hEffPromptCuts,"D^{+} Task","l");
  
  cComparison->cd(1);
  hEffPromptCF->GetYaxis()->SetRangeUser(0,1);
  hEffPromptCF->Draw();
  hEffPromptCF->SetStats(0);
  hEffPromptCF->GetYaxis()->SetTitle("Efficiency");
  hEffPromptCF->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hEffPromptCuts->Draw("E0same");
  legend->Draw("same");
  
  cComparison->cd(2);
  hEffBFeedCF->GetYaxis()->SetRangeUser(0,1);
  hEffBFeedCF->Draw();
  hEffBFeedCF->SetStats(0);
  hEffBFeedCF->GetYaxis()->SetTitle("Efficiency");
  hEffBFeedCF->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hEffBFeedCuts->Draw("E0same");
  legend->Draw("same");

  cComparison->Print("CFComp.pdf");

  TH1F *hRatioPrompt = new TH1F("hRatioPrompt","",nPtBins,ptBins);
  hRatioPrompt->Divide(hEffPromptCF,hEffPromptCuts,1.,1.,"B");
  TH1F *hRatioFD = new TH1F("hRatioFD","",nPtBins,ptBins);
  hRatioFD->Divide(hEffBFeedCF,hEffBFeedCuts,1.,1.,"B");
  
  TCanvas *cRatio = new TCanvas("cRatio","",10,10,1200,900);
  cRatio->Clear();
  cRatio->Divide(2,1);

  cRatio->cd(1);
  hRatioPrompt->Draw();
  cRatio->cd(2);
  hRatioFD->Draw();
 
  TFile *OutFile = new TFile("ratios.root","recreate");
  hRatioPrompt->Write();
  hRatioFD->Write();
  OutFile->Close();
}

