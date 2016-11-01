void PlotEffTimesAcc(TString cutset="cutset1") {

  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetLegendBorderSize(0.);
  
  TFile infile(Form("efficiency_Dplus_%s.root",cutset.Data()),"UPDATE");
  TH1F* hPrompt = (TH1F*)infile.Get("hEffD");
  TH1F* hFD = (TH1F*)infile.Get("hEffB");
  TH1F* hAcc = (TH1F*)infile.Get("hAccToy");
  hPrompt->SetDirectory(0);
  hFD->SetDirectory(0);
  hAcc->SetDirectory(0);
  infile.Close();
  hPrompt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hPrompt->GetYaxis()->SetTitle("Efficiency, Acceptance");
  hPrompt->GetXaxis()->SetTitleFont(42);
  hPrompt->GetYaxis()->SetTitleFont(42);
  hPrompt->GetXaxis()->SetLabelFont(42);
  hPrompt->GetYaxis()->SetLabelFont(42); 
  hPrompt->GetXaxis()->SetTitleSize(0.045);
  hPrompt->GetYaxis()->SetTitleSize(0.045);
  hPrompt->GetXaxis()->SetTitleOffset(1.2);
  hPrompt->GetYaxis()->SetTitleOffset(1.4);
  hPrompt->SetMarkerStyle(20);
  hPrompt->SetLineWidth(1);
  hFD->SetMarkerStyle(20);
  hFD->SetLineWidth(1);  
  hAcc->SetMarkerStyle(20);
  hAcc->SetLineWidth(1);
  hAcc->SetMarkerColor(kGreen+3);
  hAcc->SetMarkerStyle(kGreen+3);
  hAcc->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hAcc->GetYaxis()->SetTitle("Acceptance");
  hAcc->GetXaxis()->SetTitleSize(0.045);
  hAcc->GetYaxis()->SetTitleSize(0.045);
  hAcc->GetXaxis()->SetTitleOffset(1.2);
  hAcc->GetYaxis()->SetTitleOffset(1.4);
  
  TLegend* l = new TLegend(0.4,0.3,0.85,0.5);
  l->SetTextSize(0.04);
  l->AddEntry(hPrompt,"Prompt efficiency","lpe");
  l->AddEntry(hFD,"Feed-down efficiency","lpe");
  l->AddEntry(hAcc,"Acceptance","lpe");
  TCanvas* c = new TCanvas("c","",800,800);
  c->SetLogy();
  hPrompt->GetYaxis()->SetRangeUser(0.002,1.8);
  hPrompt->Draw();
  hFD->Draw("same");
  hAcc->Draw("same");
  l->Draw("same");

  TCanvas* cAcc = new TCanvas("cAcc","",800,800);
  cAcc->SetLogy();
  hAcc->Draw();

  cAcc->SaveAs(Form("Acc_Dplus_%s.eps",cutset.Data()));
  delete cAcc;
  c->SaveAs(Form("AccEff_Dplus_%s.eps",cutset.Data()));
  delete c;
}
