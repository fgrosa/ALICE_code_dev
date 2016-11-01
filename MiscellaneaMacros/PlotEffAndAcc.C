void PlotEffTimesAcc() {

  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetLegendBorderSize(0.);
  
  TFile infile("efficiencies.root","UPDATE");
  TH1F* hPrompt = (TH1F*)infile.Get("hEff_C");
  TH1F* hFD = (TH1F*)infile.Get("hEff_B");
  TH1F* hAcc = (TH1F*)infile.Get("hAccToy");
  hPrompt->SetDirectory(0);
  hFD->SetDirectory(0);
  hAcc->SetDirectory(0);
  infile.Close();
  hPrompt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hPrompt->GetYaxis()->SetTitle("Efficiency");
  hPrompt->GetXaxis()->SetTitleFont(42);
  hPrompt->GetYaxis()->SetTitleFont(42);
  hPrompt->GetXaxis()->SetLabelFont(42);
  hPrompt->GetYaxis()->SetLabelFont(42); 
  hPrompt->GetXaxis()->SetTitleSize(0.05);
  hPrompt->GetYaxis()->SetTitleSize(0.05);
  hPrompt->GetXaxis()->SetLabelSize(0.05);
  hPrompt->GetYaxis()->SetLabelSize(0.05);
  hPrompt->GetXaxis()->SetTitleOffset(1.2);
  hPrompt->GetYaxis()->SetTitleOffset(1.4);
  hPrompt->SetMarkerStyle(20);
  hPrompt->SetLineWidth(2);
  hPrompt->SetMarkerSize(1.5);
  hFD->SetMarkerStyle(21);
  hFD->SetMarkerSize(1.5);
  hFD->SetLineWidth(2);  
  hAcc->SetMarkerStyle(20);
  hAcc->SetLineWidth(2);
  hAcc->SetMarkerSize(1.5);
  hAcc->SetMarkerColor(kGreen+3);
  hAcc->SetMarkerStyle(kGreen+3);
  hAcc->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hAcc->GetYaxis()->SetTitle("Acceptance");
  hAcc->GetXaxis()->SetTitleSize(0.05);
  hAcc->GetYaxis()->SetTitleSize(0.05);
  hAcc->GetXaxis()->SetLabelSize(0.05);
  hAcc->GetYaxis()->SetLabelSize(0.05);
  hAcc->GetXaxis()->SetTitleOffset(1.2);
  hAcc->GetYaxis()->SetTitleOffset(1.4);

  for(Int_t iPt=0; iPt<hAcc->GetNbinsX(); iPt++) {
    Double_t pt = hAcc->GetBinCenter(iPt+1);
    Double_t fidacc = 0.8;
    if(pt<=5)
      fidacc = -0.2/15*pt*pt+1.9/15*pt+0.5;
    cout << pt << "  " << fidacc <<endl;
    hAcc->SetBinContent(iPt+1,hAcc->GetBinContent(iPt+1)/(2*fidacc));
  }
  
  TLegend* l = new TLegend(0.35,0.3,0.85,0.5);
  l->SetTextSize(0.045);
  l->AddEntry(hPrompt,"Prompt efficiency","lpe");
  l->AddEntry(hFD,"Feed-down efficiency","lpe");
  TCanvas* c = new TCanvas("c","",800,800);
  c->SetLogy();
  hPrompt->GetYaxis()->SetRangeUser(0.002,1.);
  hPrompt->Draw();
  hFD->Draw("same");
  l->Draw("same");

  TCanvas* cAcc = new TCanvas("cAcc","",800,800);
//  cAcc->SetLogy();
  hAcc->GetYaxis()->SetRangeUser(0.1,1.);
  hAcc->Draw();

  cAcc->SaveAs("Acc_Dplus.eps");
  delete cAcc;
  c->SaveAs("Eff_Dplus.eps");
  delete c;
}
