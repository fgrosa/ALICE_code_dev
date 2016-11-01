  void Eff() {

  TFile* file = new TFile("LHC13d3plus.root","update");
  TDirectoryFile *df = (TDirectoryFile*)file->Get("PWG3_D2H_CFtaskDplustoKpipi_PromptMultCNtrk_C");
  AliCFContainer* cont=(AliCFContainer*)df->Get("CFHFcontainer_DplustoKpipi_PromptMultCNtrk_C");
  TDirectoryFile *dfB = (TDirectoryFile*)file->Get("PWG3_D2H_CFtaskDplustoKpipi_FromBMultCNtrk_C"); 
  AliCFContainer* contB=(AliCFContainer*)df->Get("CFHFcontainer_DplustoKpipi_FromBMultCNtrk_C");
  
  TH1F *hpt0 = (TH1F*)cont->Project(0,0);
  TH1F *hpt1 = (TH1F*)cont->Project(1,0);
  TH1F *hpt2 = (TH1F*)cont->Project(2,0);
  TH1F *hpt3 = (TH1F*)cont->Project(3,0);
  TH1F *hpt4 = (TH1F*)cont->Project(4,0);
  TH1F *hpt5 = (TH1F*)cont->Project(5,0);
  TH1F *hpt6 = (TH1F*)cont->Project(6,0);
  TH1F *hpt7 = (TH1F*)cont->Project(7,0);
  TH1F *hpt8 = (TH1F*)cont->Project(8,0);
  TH1F *hpt9 = (TH1F*)cont->Project(9,0);
  
  TH1F *hy0 =(TH1F*)cont->Project(0,1); 
  TH1F *hy1 = (TH1F*)cont->Project(1,1);
  TH1F *hy2 = (TH1F*)cont->Project(2,1);
  TH1F *hy3 = (TH1F*)cont->Project(3,1);
  TH1F *hy4 = (TH1F*)cont->Project(4,1);
  TH1F *hy5 = (TH1F*)cont->Project(5,1);
  TH1F *hy6 = (TH1F*)cont->Project(6,1);
  TH1F *hy7 = (TH1F*)cont->Project(7,1);
  TH1F *hy8 = (TH1F*)cont->Project(8,1);
  TH1F *hy9 = (TH1F*)cont->Project(9,1);

  //Efficiency calculation
  //Efficiency vs pT
  Int_t Nptbins =hpt0->GetNbinsX();
  Float_t *ptbins = new Float_t[Nptbins+1];
  
  for(Int_t i = 0; i < Nptbins; i++) {
    ptbins[i] = hpt0->GetBinLowEdge(i+1);
  }

  ptbins[Nptbins]=hpt0->GetBinLowEdge(Nptbins)+hpt0->GetBinWidth(Nptbins);
  
  TH1F *hptEff28 = new TH1F("hptEff28","Efficiency vs pt (RecoCuts/MCAcc)",Nptbins,ptbins);
  TH1F *hptEff29 = new TH1F("hptEff29","Efficiency vs pt (RecoPID/MCAcc)",Nptbins,ptbins);

  for(Int_t i = 0; i< Nptbins; i++) {
    hptEff28->SetBinContent(i+1,(Double_t)hpt8->GetBinContent(i+1)/hpt2->GetBinContent(i+1));
    hptEff29->SetBinContent(i+1,(Double_t)hpt9->GetBinContent(i+1)/hpt2->GetBinContent(i+1));
  }
  
  //Efficiency vs y
  Int_t Nybins =hy0->GetNbinsX();
  Float_t *ybins = new Float_t[Nybins+1];
  
  for(Int_t i = 0; i < Nybins; i++) {
    ybins[i] = hy0->GetBinLowEdge(i+1);
  }
  
  ybins[Nybins]=hy0->GetBinLowEdge(Nybins)+hy0->GetBinWidth(Nybins);

  TH1F *hyEff28 = new TH1F("hyEff28","Efficiency vs y (RecoCuts/MCAcc)",Nybins,ybins);
  TH1F *hyEff29 = new TH1F("hyEff29","Efficiency vs y (RecoPID/MCAcc)",Nybins,ybins);
  
  for(Int_t i = 0; i< Nybins; i++) {
    if(hy2->GetBinContent(i+1)!=0) {
	  hyEff28->SetBinContent(i+1,(Double_t)hy8->GetBinContent(i+1)/hy2->GetBinContent(i+1));
	  hyEff29->SetBinContent(i+1,(Double_t)hy9->GetBinContent(i+1)/hy2->GetBinContent(i+1));
	}
  }
  
  //Acceptance calculation
  //Acceptance vs pt
  TH1F *hptAcc = new TH1F("hptAcc","Acceptance vs pt (MCAcc/MCLimAcc)",Nptbins,ptbins);
  
  for(Int_t i = 0; i< Nybins; i++) {
    hptAcc->SetBinContent(i+1,(Double_t)hpt2->GetBinContent(i+1)/hpt0->GetBinContent(i+1));
  }
  
  //Acceptance vs y
  TH1F *hyAcc = new TH1F("hyAcc","Acceptance vs y (MCAcc/MCLimAcc)",Nybins,ybins);
  
  for(Int_t i = 0; i< Nybins; i++) {
    if(hy0->GetBinContent(i+1)!=0){
	  hyAcc->SetBinContent(i+1,(Double_t)hy2->GetBinContent(i+1)/hy0->GetBinContent(i+1));
	}
  }
  
  TLegend *legend1 = new TLegend(0.58,0.48,0.87,0.87);
  legend1->SetBorderSize(0);
  legend1->SetFillColor(kWhite);
  legend1->SetTextSizePixels(1000);
  legend1->SetTextFont(42);
  legend1->SetTextSize(0.035);
  legend1->AddEntry(hpt0,"kStepMCLimAcc","lpe");
  legend1->AddEntry(hpt1,"kStepMC","lpe");
  legend1->AddEntry(hpt2,"kStepMCAcc","lpe");
  legend1->AddEntry(hpt3,"kStepRecoVertex","lpe");
  legend1->AddEntry(hpt4,"kStepRecoRefit","lpe");
  legend1->AddEntry(hpt5,"kStepReco","lpe");
  legend1->AddEntry(hpt6,"kStepRecoAcc","lpe");
  legend1->AddEntry(hpt7,"kStepRecoITSCluster","lpe");
  legend1->AddEntry(hpt8,"kStepRecoCuts","lpe");
  legend1->AddEntry(hpt9,"kStepRecoPID","lpe");

  TLegend *legend2 = new TLegend(0.35,0.15,0.69,0.54);
  legend2->SetBorderSize(0);
  legend2->SetFillColor(kWhite);
  legend2->SetTextSizePixels(1000);
  legend2->SetTextFont(42);
  legend2->SetTextSize(0.035);
  legend2->AddEntry(hy0,"kStepMCLimAcc","lpe");
  legend2->AddEntry(hy1,"kStepMC","lpe");
  legend2->AddEntry(hy2,"kStepMCAcc","lpe");
  legend2->AddEntry(hy3,"kStepRecoVertex","lpe");
  legend2->AddEntry(hy4,"kStepRecoRefit","lpe");
  legend2->AddEntry(hy5,"kStepReco","lpe");
  legend2->AddEntry(hy6,"kStepRecoAcc","lpe");
  legend2->AddEntry(hy7,"kStepRecoITSCluster","lpe");
  legend2->AddEntry(hy8,"kStepRecoCuts","lpe");
  legend2->AddEntry(hy9,"kStepRecoPID","lpe");
  
  TCanvas *cPt = new TCanvas("cPt","cPt",10,10,1200,900);
  cPt->Clear();
  cPt->SetLogy();

  hpt0->GetYaxis()->SetRangeUser(50,600000);
  hpt0->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hpt0->GetYaxis()->SetTitle("dN/dp_{T}");
  hpt0->SetTitle("");
  hpt0->SetStats(0);
  hpt1->SetStats(0);
  hpt2->SetStats(0);
  hpt3->SetStats(0);
  hpt4->SetStats(0);
  hpt5->SetStats(0);
  hpt6->SetStats(0);
  hpt7->SetStats(0);
  hpt8->SetStats(0);
  hpt9->SetStats(0);
  
  hpt0->SetLineColor(kBlue);
  hpt0->SetMarkerColor(kBlue);
  hpt1->SetLineColor(kRed);
  hpt1->SetMarkerColor(kRed);
  hpt2->SetLineColor(kGreen+3);
  hpt2->SetMarkerColor(kGreen+3);
  hpt3->SetLineColor(kBlack);
  hpt3->SetMarkerColor(kBlack);
  hpt4->SetLineColor(kYellow);
  hpt4->SetMarkerColor(kYellow);
  hpt5->SetLineColor(kMagenta);
  hpt5->SetMarkerColor(kMagenta);
  hpt6->SetLineColor(kAzure);
  hpt6->SetMarkerColor(kAzure);
  hpt7->SetLineColor(kGreen);
  hpt7->SetMarkerColor(kGreen);
  hpt8->SetLineColor(kOrange+10);
  hpt8->SetMarkerColor(kOrange+10);
  hpt9->SetLineColor(kOrange+3);
  hpt9->SetMarkerColor(kOrange+3);

  hpt0->SetMarkerStyle(20);
  hpt1->SetMarkerStyle(20);
  hpt2->SetMarkerStyle(20);
  hpt3->SetMarkerStyle(20);
  hpt4->SetMarkerStyle(20);
  hpt5->SetMarkerStyle(20);
  hpt6->SetMarkerStyle(20);
  hpt7->SetMarkerStyle(20);
  hpt8->SetMarkerStyle(20);
  hpt9->SetMarkerStyle(20);

  hpt0->SetMarkerSize(0.8);
  hpt1->SetMarkerSize(0.8);
  hpt2->SetMarkerSize(0.8);
  hpt3->SetMarkerSize(0.8);
  hpt4->SetMarkerSize(0.8);
  hpt5->SetMarkerSize(0.8);
  hpt6->SetMarkerSize(0.8);
  hpt7->SetMarkerSize(0.8);
  hpt9->SetMarkerSize(0.8);
  hpt8->SetMarkerSize(0.8);
  
  hpt0->Draw();
  hpt1->Draw("same");
  hpt2->Draw("same");
  hpt3->Draw("same");
  hpt4->Draw("same");
  hpt5->Draw("same");
  hpt6->Draw("same");
  hpt7->Draw("same");
  hpt8->Draw("same");
  hpt9->Draw("same");
  legend1->Draw("same");
  
  TCanvas *cPt2 = new TCanvas("cPt2","cPt2",10,10,1200,900);
  cPt2->Clear();
  cPt2->Divide(5,2);
  
  cPt2->cd(1)->SetLogy();
  hpt0->Draw();
  cPt2->cd(2)->SetLogy();
  hpt1->Draw();
  cPt2->cd(3)->SetLogy();
  hpt2->Draw();
  cPt2->cd(4)->SetLogy();
  hpt3->Draw();
  cPt2->cd(5)->SetLogy();
  hpt4->Draw();
  cPt2->cd(6)->SetLogy();
  hpt5->Draw();
  cPt2->cd(7)->SetLogy();
  hpt6->Draw();
  cPt2->cd(8)->SetLogy();
  hpt7->Draw();
  cPt2->cd(9)->SetLogy();
  hpt8->Draw();
  cPt2->cd(10)->SetLogy();
  hpt9->Draw();
  
  TCanvas *cy = new TCanvas("cy","cy",10,10,1200,900);
  cy->Clear();
  cy->SetLogy();

  hy0->GetYaxis()->SetRangeUser(0.1,200000);
  hy0->GetXaxis()->SetTitle("y");
  hy0->GetYaxis()->SetTitle("dN/dy");
  hy0->SetTitle("");
  hy0->SetStats(0);
  hy1->SetStats(0);
  hy2->SetStats(0);
  hy3->SetStats(0);
  hy4->SetStats(0);
  hy5->SetStats(0);
  hy6->SetStats(0);
  hy7->SetStats(0);
  hy8->SetStats(0);
  hy9->SetStats(0);

  hy0->SetLineColor(kBlue);
  hy0->SetMarkerColor(kBlue);
  hy1->SetLineColor(kRed);
  hy1->SetMarkerColor(kRed);
  hy2->SetLineColor(kGreen+3);
  hy2->SetMarkerColor(kGreen+3);
  hy3->SetLineColor(kBlack);
  hy3->SetMarkerColor(kBlack);
  hy4->SetLineColor(kYellow);
  hy4->SetMarkerColor(kYellow);
  hy5->SetLineColor(kMagenta);
  hy5->SetMarkerColor(kMagenta);
  hy6->SetLineColor(kAzure);
  hy6->SetMarkerColor(kAzure);
  hy7->SetLineColor(kGreen);
  hy7->SetMarkerColor(kGreen);
  hy8->SetLineColor(kOrange+10);
  hy8->SetMarkerColor(kOrange+10);
  hy9->SetLineColor(kOrange+3);
  hy9->SetMarkerColor(kOrange+3);
  
  hy0->SetMarkerStyle(20);
  hy1->SetMarkerStyle(20);
  hy2->SetMarkerStyle(20);
  hy3->SetMarkerStyle(20);
  hy4->SetMarkerStyle(20);
  hy5->SetMarkerStyle(20);
  hy6->SetMarkerStyle(20);
  hy7->SetMarkerStyle(20);
  hy8->SetMarkerStyle(20);
  hy9->SetMarkerStyle(20);
  
  hy0->SetMarkerSize(0.8);
  hy1->SetMarkerSize(0.8);
  hy2->SetMarkerSize(0.8);
  hy3->SetMarkerSize(0.8);
  hy4->SetMarkerSize(0.8);
  hy5->SetMarkerSize(0.8);
  hy6->SetMarkerSize(0.8);
  hy7->SetMarkerSize(0.8);
  hy8->SetMarkerSize(0.8);
  hy9->SetMarkerSize(0.8);
  
  hy0->Draw();
  hy1->Draw("same");
  hy2->Draw("same");
  hy3->Draw("same");
  hy4->Draw("same");
  hy5->Draw("same");
  hy6->Draw("same");
  hy7->Draw("same");
  hy8->Draw("same");
  hy9->Draw("same");
  legend2->Draw("same");
  
  TCanvas *cy2 = new TCanvas("cy2","cy2",10,10,1200,900);
  cy2->Clear();
  cy2->Divide(5,2);
  
  cy2->cd(1)->SetLogy();
  hy0->Draw();
  cy2->cd(2)->SetLogy();
  hy1->Draw();
  cy2->cd(3)->SetLogy();
  hy2->Draw();
  cy2->cd(4)->SetLogy();
  hy3->Draw();
  cy2->cd(5)->SetLogy();
  hy4->Draw();
  cy2->cd(6)->SetLogy();
  hy5->Draw();
  cy2->cd(7)->SetLogy();
  hy6->Draw();
  cy2->cd(8)->SetLogy();
  hy7->Draw();
  cy2->cd(9)->SetLogy();
  hy8->Draw();
  cy2->cd(10)->SetLogy();
  hy9->Draw();

  TCanvas *cptEff = new TCanvas("cptEff","cptEff",10,10,1200,900);
  cptEff->Clear();
  cptEff->Divide(2,1);

  cptEff->cd(1)->SetLogy();
  hptEff28->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hptEff28->GetYaxis()->SetTitle("Efficiency");  
  hptEff28->Draw();
  
  cptEff->cd(2)->SetLogy();
  hptEff29->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hptEff29->GetYaxis()->SetTitle("Efficiency");
  hptEff29->Draw();

  TCanvas *cyEff = new TCanvas("cyEff","cyEff",10,10,1200,900);
  cyEff->Clear();
  cyEff->Divide(2,1);

  cyEff->cd(1)->SetLogy();
  hyEff28->GetXaxis()->SetTitle("y");
  hyEff28->GetYaxis()->SetTitle("Efficiency");
  hyEff28->Draw();
  
  cyEff->cd(2)->SetLogy();
  hyEff29->GetXaxis()->SetTitle("y");
  hyEff29->GetYaxis()->SetTitle("Efficiency");
  hyEff29->Draw();

  TCanvas *cptAcc = new TCanvas("cptAcc","cptAcc",10,10,1200,900);
  cptAcc->Clear();

  hptAcc->Draw();

  TCanvas *cyAcc = new TCanvas("cyAcc","cyAcc",10,10,1200,900);
  cyAcc->Clear();

  hyAcc->Draw();

  cPt->Print("PTsteps.pdf");
  cy->Print("Ysteps.pdf");
  cptEff->Print("PTeff.pdf");
  cyEff->Print("Yeff.pdf");
  
}
