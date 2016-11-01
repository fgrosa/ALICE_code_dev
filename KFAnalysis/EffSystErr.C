void EffSystErr() {

  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHistLineWidth(2);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.05,"xyz");
  
  ///for prompt -> efficiency variation
  TFile NoWeightFile("efficiency_Dplus_cutset1_noPtReweight.root","UPDATE");
  TH1F* hEffPrompt = (TH1F*)NoWeightFile.Get("hEffD");
  hEffPrompt->SetDirectory(0);;
  TH1F* hEffFD = (TH1F*)NoWeightFile.Get("hEffB");
  hEffFD->SetDirectory(0);
  NoWeightFile.Close();
  
  TFile WeightFile("efficiency_Dplus_cutset1_PtReweight.root","UPDATE");
  TH1F* hEffPromptWeight = (TH1F*)WeightFile.Get("hEffD");
  hEffPromptWeight->SetDirectory(0);;
  TH1F* hEffFDWeight = (TH1F*)WeightFile.Get("hEffB");
  hEffFDWeight->SetDirectory(0);
  WeightFile.Close();

  TH1F* hRatioPrompt = (TH1F*)hEffPrompt->Clone();
  hRatioPrompt->Divide(hEffPrompt,hEffPromptWeight,1.,1.,"");
  TH1F* hRatioFD = (TH1F*)hEffFD->Clone();
  hRatioFD->Divide(hEffFD,hEffFDWeight,1.,1.,"");

  TGaxis::SetMaxDigits(2);

  hRatioPrompt->SetLineWidth(2);
  hRatioFD->SetLineWidth(2);
  hRatioPrompt->GetYaxis()->SetTitleSize(0.05);
  hRatioPrompt->GetXaxis()->SetTitleSize(0.05);
  hRatioPrompt->GetYaxis()->SetLabelSize(0.05);
  hRatioPrompt->GetXaxis()->SetLabelSize(0.05);
  hRatioPrompt->GetYaxis()->SetTitleOffset(1.6);
  hRatioPrompt->GetYaxis()->SetNdivisions(505);
  hRatioFD->GetYaxis()->SetTitleSize(0.05);
  hRatioFD->GetXaxis()->SetTitleSize(0.05);
  hRatioFD->GetYaxis()->SetLabelSize(0.05);
  hRatioFD->GetXaxis()->SetLabelSize(0.05);
  hRatioFD->GetYaxis()->SetTitleOffset(1.6);
  hRatioFD->GetYaxis()->SetNdivisions(505);
  hRatioFD->SetMarkerStyle(21);
  hRatioFD->SetMarkerSize(1.5);
  
  for(Int_t iPt=0; iPt<hRatioPrompt->GetNbinsX(); iPt++) {
    hRatioPrompt->SetBinError(iPt+1,1.e-10);
    hRatioFD->SetBinError(iPt+1,1.e-10);
  }
  
  ///for FD -> fprompt variation
  TFile NoWeightFpromptFile("HFPtSpectrum_combinedFD_cutset1_noPtReweight.root","UPDATE");
  TGraphAsymmErrors* gFprompt = (TGraphAsymmErrors*)NoWeightFpromptFile.Get("gFcCorrConservative");
  NoWeightFpromptFile.Close();
  
  TFile WeightFpromptFile("HFPtSpectrum_combinedFD_cutset1_PtReweight.root","UPDATE");
  TGraphAsymmErrors* gFpromptWeight = (TGraphAsymmErrors*)WeightFpromptFile.Get("gFcCorrConservative");
  WeightFpromptFile.Close();

  TH1F* hFpromptRatio=(TH1F*)hRatioFD->Clone();
  hFpromptRatio->GetYaxis()->SetRangeUser(0.99,1.01);
  hFpromptRatio->GetYaxis()->SetTitle("#it{f}_{prompt}(FONLL)/#it{f}_{prompt}(Pythia)");
  hFpromptRatio->GetYaxis()->SetNdivisions(503);
  hFpromptRatio->GetYaxis()->SetTitleOffset(1.4);
 
  for(Int_t iPt=0; iPt<hEffPrompt->GetNbinsX(); iPt++) {
    Double_t pt;
    Double_t frac;
    Double_t fracweight;
    gFprompt->GetPoint(iPt+1,pt,frac);
    gFpromptWeight->GetPoint(iPt+1,pt,fracweight);
    hFpromptRatio->SetBinContent(iPt+1,fracweight/frac);
    hFpromptRatio->SetBinError(iPt+1,1.e-10);
  }
  
  TCanvas* cPrompt = new TCanvas("cPrompt","",800,800);
  cPrompt->SetLogy();
  hEffPrompt->Draw();
  hEffPromptWeight->Draw("same");
  
  TCanvas* cFD = new TCanvas("cFD","",800,800);
  cFD->SetLogy();
  hEffFD->Draw();
  hEffFDWeight->Draw("same");
  
  TCanvas* cRatioPrompt = new TCanvas("cRatioPrompt","",800,800);
  hRatioPrompt->GetYaxis()->SetTitle("#epsilon_{prompt}(FONLL)/#epsilon_{prompt}(Pythia)");
  hRatioPrompt->GetYaxis()->SetRangeUser(0.98,1.01);
  hRatioPrompt->SetMarkerSize(1.5);
  hRatioPrompt->Draw("E");

  TCanvas* cRatioFD = new TCanvas("cRatioFD","",800,800);
  hRatioFD->GetYaxis()->SetTitle("#epsilon_{feed-down}(FONLL)/#epsilon_{feed-down}(Pythia)");
  hRatioFD->GetYaxis()->SetRangeUser(0.97,1.01);
  hRatioFD->SetMarkerStyle(21);
  hRatioFD->SetMarkerSize(1.5);
  hRatioFD->Draw("E");

  TCanvas* cRatioFprompt = new TCanvas("cRatioFprompt","",800,800);
  hFpromptRatio->Draw("E");
  
  TH1F* hRatioPromptCopy = (TH1F*)hRatioPrompt->Clone();
  hRatioPromptCopy->GetYaxis()->SetTitle("Efficiency ratio FONLL/MC");
  hRatioPromptCopy->SetMarkerSize(1.5);
  TH1F* hRatioFDCopy = (TH1F*)hRatioFD->Clone();
  hRatioFDCopy->SetMarkerStyle(21);
  hRatioFDCopy->SetMarkerSize(1.5);
  
  TLegend* l = new TLegend(0.45,0.7,0.89,0.85);
  l->SetTextSize(0.05);
  l->AddEntry(hRatioPromptCopy,"Prompt","lpe");
  l->AddEntry(hRatioFDCopy,"FD","lpe");
  
  TCanvas *cRatio = new TCanvas("cRatio","",800,800);
  hRatioPromptCopy->Draw("E");
  hRatioFDCopy->Draw("Esame");
  l->Draw("same");
  
  cRatioPrompt->SaveAs("RatioEffPrompt_KF.eps");
  cRatioFD->SaveAs("RatioEffFD_KF.eps");
  cRatio->SaveAs("KF_ratio_eff.eps");
  cRatioFprompt->SaveAs("KF_ratio_fprompt.eps");
}
