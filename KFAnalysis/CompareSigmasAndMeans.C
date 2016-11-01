void CompareSigmasAndMeans() {

  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetTitleSize(0.05,"xyt");
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetTitleOffset(1.6,"y");
  gStyle->SetTitleOffset(1.3,"x");
  gStyle->SetLegendBorderSize(0);
  TGaxis::SetMaxDigits(2);
  
  TFile infileMC("rawyieldsMC_Dplus_cutset1.root","UPDATE");
  TH1F* hSigmaMC = (TH1F*)infileMC.Get("hSigma");
  TH1F* hMeanMC = (TH1F*)infileMC.Get("hMean");
  hSigmaMC->SetDirectory(0);
  hMeanMC->SetDirectory(0);
  infileMC.Close();
  
  TFile infile("rawyields_Dplus_cutset1.root","UPDATE");
  TH1F* hSigma = (TH1F*)infile.Get("hSigma");
  TH1F* hMean = (TH1F*)infile.Get("hMean");
  hSigma->SetDirectory(0);
  hMean->SetDirectory(0);
  infile.Close();

  Double_t massD = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  TLine* line = new TLine(1,massD,24,massD);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->SetLineColor(kGreen+3);
  
  TLegend *l = new TLegend(0.2,0.65,0.8,0.85);
  l->SetTextSize(0.05);
  l->AddEntry(hSigma,"fit on data","lpe");
  l->AddEntry(hSigmaMC,"fit on MC","lpe"); 
  TLegend *l2 = new TLegend(0.2,0.65,0.8,0.89);
  l2->SetTextSize(0.05);
  l2->AddEntry(hMean,"fit on data","lpe");
  l2->AddEntry(hMeanMC,"fit on MC","lpe");  
  l2->AddEntry(line,"PDG value","l");
  
  TCanvas* cSigma = new TCanvas("cSigma","",800,800);
  hSigmaMC->GetYaxis()->SetRangeUser(0.005,0.025);
  hSigmaMC->GetYaxis()->SetTitle("Sigma (GeV/c^{2})");
  hSigmaMC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hSigmaMC->GetYaxis()->SetTitleOffset(1.5);
  hSigmaMC->GetXaxis()->SetTitleOffset(1.2);
  hSigmaMC->GetYaxis()->SetLabelSize(0.05);
  hSigmaMC->GetXaxis()->SetLabelSize(0.05);
  hSigmaMC->GetYaxis()->SetTitleSize(0.05);
  hSigmaMC->GetXaxis()->SetTitleSize(0.05);
  hSigmaMC->SetLineColor(kRed);
  hSigmaMC->SetMarkerColor(kRed);
  hSigmaMC->SetMarkerStyle(21);
  hSigmaMC->SetMarkerSize(1.5);
  hSigmaMC->SetLineWidth(2);
  hSigma->SetLineColor(kBlue);
  hSigma->SetMarkerColor(kBlue);
  hSigma->SetMarkerStyle(20);
  hSigma->SetMarkerSize(1.5);
  hSigma->SetLineWidth(2);
  hSigmaMC->Draw();
  hSigma->Draw("same");
  l->Draw("same");

  TCanvas* cMean = new TCanvas("cMean","",800,800);
  hMeanMC->GetYaxis()->SetRangeUser(1.865,1.880);
  hMeanMC->GetYaxis()->SetTitle("Mean (GeV/c^{2})");
  hMeanMC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hMeanMC->GetYaxis()->SetTitleOffset(1.8);
  hMeanMC->GetXaxis()->SetTitleOffset(1.2);
  hMeanMC->SetMarkerSize(1.5);
  hMeanMC->SetLineWidth(2);
  hMeanMC->SetLineColor(kRed);
  hMeanMC->SetMarkerColor(kRed);
  hMeanMC->SetMarkerStyle(21);
  hMeanMC->GetYaxis()->SetLabelSize(0.05);
  hMeanMC->GetXaxis()->SetLabelSize(0.05);
  hMeanMC->GetYaxis()->SetTitleSize(0.05);
  hMeanMC->GetXaxis()->SetTitleSize(0.05);
  hMean->SetMarkerSize(1.5);
  hMean->SetLineWidth(2);
  hMean->SetLineColor(kBlue);
  hMean->SetMarkerColor(kBlue);
  hMean->SetMarkerStyle(20);
  hMeanMC->Draw();
  hMean->Draw("same");
  line->Draw("same");
  l2->Draw("same");
  
  cSigma->SaveAs("SigmaComparison.eps");
  delete cSigma;
  cMean->SaveAs("MeanComparison.eps");
  delete cMean;
}
