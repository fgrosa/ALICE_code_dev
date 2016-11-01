void RawYieldsComparison(TString filename1="RawYields.root", TString filename2="RawYields_d0cut.root", TString name1="w/o d_{0} cut", TString name2="|d_{0}| < 100 #mum") {
  
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.3);
  
  TFile file1(filename1,"READ");
  TH1F* hRawYields1 = (TH1F*)file1.Get("hRawYields");
  TH1F* hRawYieldsSigma1 = (TH1F*)file1.Get("hRawYieldsSigma");
  TH1F* hRawYieldsMean1 = (TH1F*)file1.Get("hRawYieldsMean");
  TH1F* hRawYieldsSignificance1 = (TH1F*)file1.Get("hRawYieldsSignificance");
  TH1F* hRawYieldsSoverB1 = (TH1F*)file1.Get("hRawYieldsSoverB");
  TH1F* hRawYieldsChiSquare1 = (TH1F*)file1.Get("hRawYieldsChiSquare");
  hRawYields1->SetDirectory(0);
  hRawYieldsSigma1->SetDirectory(0);
  hRawYieldsMean1->SetDirectory(0);
  hRawYieldsSignificance1->SetDirectory(0);
  hRawYieldsSoverB1->SetDirectory(0);
  hRawYieldsChiSquare1->SetDirectory(0);
  file1.Close();

  TFile file2(filename2,"READ");
  TH1F* hRawYields2 = (TH1F*)file2.Get("hRawYields");
  TH1F* hRawYieldsSigma2 = (TH1F*)file2.Get("hRawYieldsSigma");
  TH1F* hRawYieldsMean2 = (TH1F*)file2.Get("hRawYieldsMean");
  TH1F* hRawYieldsSignificance2 = (TH1F*)file2.Get("hRawYieldsSignificance");
  TH1F* hRawYieldsSoverB2 = (TH1F*)file2.Get("hRawYieldsSoverB");
  TH1F* hRawYieldsChiSquare2 = (TH1F*)file2.Get("hRawYieldsChiSquare");
  hRawYields2->SetDirectory(0);
  hRawYieldsSigma2->SetDirectory(0);
  hRawYieldsMean2->SetDirectory(0);
  hRawYieldsSignificance2->SetDirectory(0);
  hRawYieldsSoverB2->SetDirectory(0);
  hRawYieldsChiSquare2->SetDirectory(0);
  file2.Close();
  
  TLegend* l = new TLegend(0.5,0.6,0.89,0.85);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.06);
  l->AddEntry(hRawYields1,name1,"lpe");
  l->AddEntry(hRawYields2,name2,"lpe");
  TLegend* l2 = new TLegend(0.5,0.2,0.89,0.45);
  l2->SetFillColor(0);
  l2->SetBorderSize(0);
  l2->SetTextSize(0.06);
  l2->AddEntry(hRawYields1,name1,"lpe");
  l2->AddEntry(hRawYields2,name2,"lpe");
  TLegend* l3 = new TLegend(0.18,0.64,0.57,0.89);
  l3->SetFillColor(0);
  l3->SetFillStyle(0);
  l3->SetBorderSize(0);
  l3->SetTextSize(0.06);
  l3->AddEntry(hRawYields1,name1,"lpe");
  l3->AddEntry(hRawYields2,name2,"lpe");
  TLegend* l4 = new TLegend(0.5,0.6,0.89,0.85);
  l4->SetFillColor(0);
  l4->SetFillStyle(0);
  l4->SetBorderSize(0);
  l4->SetTextSize(0.05);
  l4->AddEntry(hRawYields1,name1,"lpe");
  l4->AddEntry(hRawYields2,name2,"lpe");
  
  TCanvas* cRawYields = new TCanvas("cRawYields","",10,10,600,800);
  cRawYields->cd();
  TPad *padRawYields1 = new TPad("padRawYields1","padRawYields1",0,0.35,1,1);
  padRawYields1->SetBottomMargin(0);
  padRawYields1->Draw();
  padRawYields1->cd();
  hRawYields1->GetYaxis()->SetTitleSize(0.07);
  hRawYields1->GetYaxis()->SetLabelSize(0.06);
  hRawYields1->GetYaxis()->SetTitleOffset(1.);
  hRawYields1->GetYaxis()->SetRangeUser(0.1,hRawYields1->GetMaximum()*1.1);
  hRawYields1->SetLineColor(kBlue);
  hRawYields1->SetMarkerColor(kBlue);
  hRawYields2->SetLineColor(kRed);
  hRawYields2->SetMarkerColor(kRed);
  hRawYields2->SetMarkerStyle(21);
  hRawYields1->Draw();
  hRawYields2->Draw("same");
  l->Draw("same");
  cRawYields->cd();
  TPad *padRawYields2 = new TPad("padRawYields2","padRawYields2",0.,0.,1,0.35);
  padRawYields2->SetTopMargin(0);
  padRawYields2->Draw();
  padRawYields2->cd();
  TH1F* hRawYieldsRatio = (TH1F*)hRawYields2->Clone();
  hRawYieldsRatio->SetDirectory(0);
  hRawYieldsRatio->Divide(hRawYields2,hRawYields1,1.,1.,"");
  for(Int_t iPt=0; iPt<hRawYieldsRatio->GetNbinsX(); iPt++) {
    hRawYieldsRatio->SetBinError(iPt+1,1.e-10);
  }
  hRawYieldsRatio->GetYaxis()->SetTitle("Ratio");
  hRawYieldsRatio->GetYaxis()->SetTitleSize(0.11);
  hRawYieldsRatio->GetXaxis()->SetTitleSize(0.11);
  hRawYieldsRatio->GetYaxis()->SetLabelSize(0.10);
  hRawYieldsRatio->GetXaxis()->SetLabelSize(0.10);
  hRawYieldsRatio->GetXaxis()->SetTitleOffset(1.2);
  hRawYieldsRatio->GetYaxis()->SetTitleOffset(0.6);
  hRawYieldsRatio->GetYaxis()->SetDecimals(2);
  hRawYieldsRatio->Draw();

  TCanvas* cRawYieldsSigma = new TCanvas("cRawYieldsSigma","",10,10,600,800);
  cRawYieldsSigma->cd();
  TPad *padRawYieldsSigma1 = new TPad("padRawYieldsSigma1","padRawYieldsSigma1",0,0.35,1,1);
  padRawYieldsSigma1->SetBottomMargin(0);
  padRawYieldsSigma1->Draw();
  padRawYieldsSigma1->cd();
  hRawYieldsSigma1->GetYaxis()->SetTitleSize(0.07);
  hRawYieldsSigma1->GetYaxis()->SetLabelSize(0.06);
  hRawYieldsSigma1->GetYaxis()->SetTitleOffset(1.2);
  hRawYieldsSigma1->SetLineColor(kBlue);
  hRawYieldsSigma1->SetMarkerColor(kBlue);
  hRawYieldsSigma2->SetLineColor(kRed);
  hRawYieldsSigma2->SetMarkerColor(kRed);
  hRawYieldsSigma2->SetMarkerStyle(21);
  hRawYieldsSigma1->Draw();
  hRawYieldsSigma2->Draw("same");
  l3->Draw("same");
  cRawYieldsSigma->cd();
  TPad *padRawYieldsSigma2 = new TPad("padRawYieldsSigma2","padRawYieldsSigma2",0.,0.,1,0.35);
  padRawYieldsSigma2->SetTopMargin(0);
  padRawYieldsSigma2->Draw();
  padRawYieldsSigma2->cd();
  TH1F* hRawYieldsSigmaRatio = (TH1F*)hRawYieldsSigma2->Clone();
  hRawYieldsSigmaRatio->SetDirectory(0);
  hRawYieldsSigmaRatio->Divide(hRawYieldsSigma2,hRawYieldsSigma1,1.,1.,"");
  for(Int_t iPt=0; iPt<hRawYieldsSigmaRatio->GetNbinsX(); iPt++) {
    hRawYieldsSigmaRatio->SetBinError(iPt+1,1.e-10);
  }
  hRawYieldsSigmaRatio->GetYaxis()->SetTitle("Ratio");
  hRawYieldsSigmaRatio->GetYaxis()->SetTitleSize(0.11);
  hRawYieldsSigmaRatio->GetXaxis()->SetTitleSize(0.11);
  hRawYieldsSigmaRatio->GetYaxis()->SetLabelSize(0.10);
  hRawYieldsSigmaRatio->GetXaxis()->SetLabelSize(0.10);
  hRawYieldsSigmaRatio->GetXaxis()->SetTitleOffset(1.2);
  hRawYieldsSigmaRatio->GetYaxis()->SetTitleOffset(0.8);
  hRawYieldsSigmaRatio->GetYaxis()->SetDecimals(2);
  hRawYieldsSigmaRatio->GetYaxis()->SetRangeUser(hRawYieldsSigmaRatio->GetMinimum()*0.9,hRawYieldsSigmaRatio->GetMaximum()*1.1);
  hRawYieldsSigmaRatio->Draw();
  
  TCanvas* cRawYieldsMean = new TCanvas("cRawYieldsMean","",10,10,600,800);
  cRawYieldsMean->cd();
  TPad *padRawYieldsMean1 = new TPad("padRawYieldsMean1","padRawYieldsMean1",0,0.35,1,1);
  padRawYieldsMean1->SetBottomMargin(0);
  padRawYieldsMean1->Draw();
  padRawYieldsMean1->cd();
  hRawYieldsMean1->GetYaxis()->SetRangeUser(hRawYieldsMean1->GetMinimum()*0.992,hRawYieldsMean1->GetMinimum()*1.005);
  hRawYieldsMean1->GetYaxis()->SetTitleSize(0.07);
  hRawYieldsMean1->GetYaxis()->SetLabelSize(0.06);
  hRawYieldsMean1->GetYaxis()->SetTitleOffset(1.2);
  hRawYieldsMean1->SetLineColor(kBlue);
  hRawYieldsMean1->SetMarkerColor(kBlue);
  hRawYieldsMean2->SetLineColor(kRed);
  hRawYieldsMean2->SetMarkerColor(kRed);
  hRawYieldsMean2->SetMarkerStyle(21);
  hRawYieldsMean1->Draw();
  hRawYieldsMean2->Draw("same");
  l2->Draw("same");
  cRawYieldsMean->cd();
  TPad *padRawYieldsMean2 = new TPad("padRawYieldsMean2","padRawYieldsMean2",0.,0.,1,0.35);
  padRawYieldsMean2->SetTopMargin(0);
  padRawYieldsMean2->Draw();
  padRawYieldsMean2->cd();
  TH1F* hRawYieldsMeanRatio = (TH1F*)hRawYieldsMean2->Clone();
  hRawYieldsMeanRatio->SetDirectory(0);
  hRawYieldsMeanRatio->Divide(hRawYieldsMean2,hRawYieldsMean1,1.,1.,"");
  for(Int_t iPt=0; iPt<hRawYieldsMeanRatio->GetNbinsX(); iPt++) {
    hRawYieldsMeanRatio->SetBinError(iPt+1,1.e-10);
  }
  hRawYieldsMeanRatio->GetYaxis()->SetTitle("Ratio");
  hRawYieldsMeanRatio->GetYaxis()->SetTitleSize(0.11);
  hRawYieldsMeanRatio->GetXaxis()->SetTitleSize(0.11);
  hRawYieldsMeanRatio->GetYaxis()->SetLabelSize(0.10);
  hRawYieldsMeanRatio->GetXaxis()->SetLabelSize(0.10);
  hRawYieldsMeanRatio->GetXaxis()->SetTitleOffset(1.2);
  hRawYieldsMeanRatio->GetYaxis()->SetTitleOffset(0.8);
  hRawYieldsMeanRatio->GetYaxis()->SetDecimals(2);
  hRawYieldsMeanRatio->GetYaxis()->SetRangeUser(hRawYieldsMeanRatio->GetMinimum()*0.999,hRawYieldsMeanRatio->GetMaximum()*1.001);
  hRawYieldsMeanRatio->Draw();
  
  TCanvas* cRawYieldsSignificance = new TCanvas("cRawYieldsSignificance","",10,10,600,800);
  cRawYieldsSignificance->cd();
  TPad *padRawYieldsSignificance1 = new TPad("padRawYieldsSignificance1","padRawYieldsSignificance1",0,0.35,1,1);
  padRawYieldsSignificance1->SetBottomMargin(0);
  padRawYieldsSignificance1->Draw();
  padRawYieldsSignificance1->cd();
  hRawYieldsSignificance1->GetYaxis()->SetRangeUser(hRawYieldsSignificance1->GetMinimum()*0.75,hRawYieldsSignificance1->GetMaximum()*1.3);
  hRawYieldsSignificance1->GetYaxis()->SetTitleSize(0.07);
  hRawYieldsSignificance1->GetYaxis()->SetLabelSize(0.06);
  hRawYieldsSignificance1->GetYaxis()->SetTitleOffset(1.);
  hRawYieldsSignificance1->SetLineColor(kBlue);
  hRawYieldsSignificance1->SetMarkerColor(kBlue);
  hRawYieldsSignificance2->SetLineColor(kRed);
  hRawYieldsSignificance2->SetMarkerColor(kRed);
  hRawYieldsSignificance2->SetMarkerStyle(21);
  hRawYieldsSignificance1->Draw();
  hRawYieldsSignificance2->Draw("same");
  l->Draw("same");
  cRawYieldsSignificance->cd();
  TPad *padRawYieldsSignificance2 = new TPad("padRawYieldsSignificance2","padRawYieldsSignificance2",0.,0.,1,0.35);
  padRawYieldsSignificance2->SetTopMargin(0);
  padRawYieldsSignificance2->Draw();
  padRawYieldsSignificance2->cd();
  TH1F* hRawYieldsSignificanceRatio = (TH1F*)hRawYieldsSignificance2->Clone();
  hRawYieldsSignificanceRatio->SetDirectory(0);
  hRawYieldsSignificanceRatio->Divide(hRawYieldsSignificance2,hRawYieldsSignificance1,1.,1.,"");
  for(Int_t iPt=0; iPt<hRawYieldsSignificanceRatio->GetNbinsX(); iPt++) {
    hRawYieldsSignificanceRatio->SetBinError(iPt+1,1.e-10);
  }
  hRawYieldsSignificanceRatio->GetYaxis()->SetTitle("Ratio");
  hRawYieldsSignificanceRatio->GetYaxis()->SetTitleSize(0.11);
  hRawYieldsSignificanceRatio->GetXaxis()->SetTitleSize(0.11);
  hRawYieldsSignificanceRatio->GetYaxis()->SetLabelSize(0.10);
  hRawYieldsSignificanceRatio->GetXaxis()->SetLabelSize(0.10);
  hRawYieldsSignificanceRatio->GetXaxis()->SetTitleOffset(1.2);
  hRawYieldsSignificanceRatio->GetYaxis()->SetTitleOffset(0.6);
  hRawYieldsSignificanceRatio->GetYaxis()->SetDecimals(2);
  hRawYieldsSignificanceRatio->GetYaxis()->SetRangeUser(hRawYieldsSignificanceRatio->GetMinimum()*0.8,hRawYieldsSignificanceRatio->GetMaximum()*1.2);
  hRawYieldsSignificanceRatio->Draw();

  TCanvas* cRawYieldsSoverB = new TCanvas("cRawYieldsSoverB","",10,10,600,800);
  cRawYieldsSoverB->cd();
  TPad *padRawYieldsSoverB1 = new TPad("padRawYieldsSoverB1","padRawYieldsSoverB1",0,0.35,1,1);
  padRawYieldsSoverB1->SetBottomMargin(0);
  padRawYieldsSoverB1->Draw();
  padRawYieldsSoverB1->cd();
  hRawYieldsSoverB1->GetYaxis()->SetRangeUser(hRawYieldsSoverB1->GetMinimum()*0.75,hRawYieldsSoverB1->GetMaximum()*1.8);
  hRawYieldsSoverB1->GetYaxis()->SetTitleSize(0.07);
  hRawYieldsSoverB1->GetYaxis()->SetLabelSize(0.06);
  hRawYieldsSoverB1->GetYaxis()->SetTitleOffset(1.);
  hRawYieldsSoverB1->SetLineColor(kBlue);
  hRawYieldsSoverB1->SetMarkerColor(kBlue);
  hRawYieldsSoverB2->SetLineColor(kRed);
  hRawYieldsSoverB2->SetMarkerColor(kRed);
  hRawYieldsSoverB2->SetMarkerStyle(21);
  hRawYieldsSoverB1->Draw();
  hRawYieldsSoverB2->Draw("same");
  l->Draw("same");
  cRawYieldsSoverB->cd();
  TPad *padRawYieldsSoverB2 = new TPad("padRawYieldsSoverB2","padRawYieldsSoverB2",0.,0.,1,0.35);
  padRawYieldsSoverB2->SetTopMargin(0);
  padRawYieldsSoverB2->Draw();
  padRawYieldsSoverB2->cd();
  TH1F* hRawYieldsSoverBRatio = (TH1F*)hRawYieldsSoverB2->Clone();
  hRawYieldsSoverBRatio->SetDirectory(0);
  hRawYieldsSoverBRatio->Divide(hRawYieldsSoverB2,hRawYieldsSoverB1,1.,1.,"");
  for(Int_t iPt=0; iPt<hRawYieldsSoverBRatio->GetNbinsX(); iPt++) {
    hRawYieldsSoverBRatio->SetBinError(iPt+1,1.e-10);
  }
  hRawYieldsSoverBRatio->GetYaxis()->SetTitle("Ratio");
  hRawYieldsSoverBRatio->GetYaxis()->SetTitleSize(0.11);
  hRawYieldsSoverBRatio->GetXaxis()->SetTitleSize(0.11);
  hRawYieldsSoverBRatio->GetYaxis()->SetLabelSize(0.10);
  hRawYieldsSoverBRatio->GetXaxis()->SetLabelSize(0.10);
  hRawYieldsSoverBRatio->GetXaxis()->SetTitleOffset(1.2);
  hRawYieldsSoverBRatio->GetYaxis()->SetTitleOffset(0.6);
  hRawYieldsSoverBRatio->GetYaxis()->SetDecimals(2);
  hRawYieldsSoverBRatio->GetYaxis()->SetRangeUser(hRawYieldsSoverBRatio->GetMinimum()*0.8,hRawYieldsSoverBRatio->GetMaximum()*1.2);
  hRawYieldsSoverBRatio->Draw();
  
  TCanvas* cChi = new TCanvas("cChi","",10,10,800,800);
  cChi->SetBottomMargin(0.12);
  hRawYieldsChiSquare1->GetYaxis()->SetTitleOffset(1.4);
  hRawYieldsChiSquare1->SetLineColor(kBlue);
  hRawYieldsChiSquare1->SetMarkerColor(kBlue);
  hRawYieldsChiSquare2->SetLineColor(kRed);
  hRawYieldsChiSquare2->SetMarkerColor(kRed);
  hRawYieldsChiSquare2->SetMarkerStyle(21);
  hRawYieldsChiSquare1->Draw();
  hRawYieldsChiSquare2->Draw("same");
  l4->Draw("same");

  cRawYields->SaveAs("Comparison_RawYields.eps");
  cRawYieldsSigma->SaveAs("Comparison_Sigma.eps");
  cRawYieldsMean->SaveAs("Comparison_Mean.eps");
  cRawYieldsSignificance->SaveAs("Comparison_Significance.eps");
  cRawYieldsSoverB->SaveAs("Comparison_SoverB.eps");
  cChi->SaveAs("Comparison_ChiSquare.eps");
  
}

void EfficiencyComparison(TString filename1="Efficiency.root", TString filename2="Efficiency_d0cut.root", TString name1="w/o d_{0} cut", TString name2="|d_{0}| < 100 #mum") {
  
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.3);
  
  TFile file1(filename1,"READ");
  TH1F* hEffPrompt1 = (TH1F*)file1.Get("hEffPrompt");
  TH1F* hEffFD1 = (TH1F*)file1.Get("hEffFD");
  hEffPrompt1->SetDirectory(0);
  hEffFD1->SetDirectory(0);
  file1.Close();
  
  TFile file2(filename2,"READ");
  TH1F* hEffPrompt2 = (TH1F*)file2.Get("hEffPrompt");
  TH1F* hEffFD2 = (TH1F*)file2.Get("hEffFD");
  hEffPrompt2->SetDirectory(0);
  hEffFD2->SetDirectory(0);
  file2.Close();
  
  TLegend* l = new TLegend(0.55,0.25,0.85,0.5);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.06);
  l->AddEntry(hEffPrompt1,name1,"lpe");
  l->AddEntry(hEffPrompt2,name2,"lpe");
  
  TCanvas* cEffPrompt = new TCanvas("cEffPrompt","",10,10,600,800);
  cEffPrompt->cd();
  TPad *padPrompt1 = new TPad("padPrompt1","padPrompt1",0,0.35,1,1);
  padPrompt1->SetBottomMargin(0);
  padPrompt1->Draw();
  padPrompt1->SetLogy();
  padPrompt1->cd();
  cEffPrompt->SetLogy();
  hEffPrompt1->GetYaxis()->SetTitle("#epsilon_{prompt}");
  hEffPrompt1->GetYaxis()->SetTitleSize(0.08);
  hEffPrompt1->GetYaxis()->SetLabelSize(0.06);
  hEffPrompt1->GetYaxis()->SetTitleOffset(0.8);
  hEffPrompt1->SetLineColor(kBlue);
  hEffPrompt1->SetMarkerColor(kBlue);
  hEffPrompt2->SetLineColor(kRed);
  hEffPrompt2->SetMarkerColor(kRed);
  hEffPrompt2->SetMarkerStyle(21);
  hEffPrompt1->Draw();
  hEffPrompt2->Draw("same");
  l->Draw("same");
  cEffPrompt->cd();
  TPad *padPrompt2 = new TPad("padPrompt2","padPrompt2",0.,0.,1,0.35);
  padPrompt2->SetTopMargin(0);
  padPrompt2->Draw();
  padPrompt2->cd();
  TH1F* hEffPromptRatio = (TH1F*)hEffPrompt2->Clone();
  hEffPromptRatio->SetDirectory(0);
  hEffPromptRatio->Divide(hEffPrompt2,hEffPrompt1,1.,1.,"");
  for(Int_t iPt=0; iPt<hEffPromptRatio->GetNbinsX(); iPt++) {
    hEffPromptRatio->SetBinError(iPt+1,1.e-10);
  }
  hEffPromptRatio->GetYaxis()->SetTitle("Ratio");
  hEffPromptRatio->GetYaxis()->SetTitleSize(0.10);
  hEffPromptRatio->GetXaxis()->SetTitleSize(0.10);
  hEffPromptRatio->GetYaxis()->SetLabelSize(0.10);
  hEffPromptRatio->GetXaxis()->SetLabelSize(0.10);
  hEffPromptRatio->GetXaxis()->SetTitleOffset(1.2);
  hEffPromptRatio->GetYaxis()->SetTitleOffset(0.6);
  hEffPromptRatio->GetYaxis()->SetDecimals(2);
  hEffPromptRatio->Draw();

  
  TCanvas* cEffFD = new TCanvas("cEffFD","",10,10,600,800);
  cEffFD->cd();
  TPad *padFD1 = new TPad("padFD1","padFD1",0,0.35,1,1);
  padFD1->SetBottomMargin(0);
  padFD1->Draw();
  padFD1->SetLogy();
  padFD1->cd();
  cEffFD->SetLogy();
  hEffFD1->GetYaxis()->SetTitle("#epsilon_{feed-down}");
  hEffFD1->GetYaxis()->SetTitleSize(0.08);
  hEffFD1->GetYaxis()->SetLabelSize(0.06);
  hEffFD1->GetYaxis()->SetTitleOffset(0.8);
  hEffFD1->SetLineColor(kBlue);
  hEffFD1->SetMarkerColor(kBlue);
  hEffFD2->SetLineColor(kRed);
  hEffFD2->SetMarkerColor(kRed);
  hEffFD2->SetMarkerStyle(21);
  hEffFD1->Draw();
  hEffFD2->Draw("same");
  l->Draw("same");
  cEffFD->cd();
  TPad *padFD2 = new TPad("padFD2","padFD2",0.,0.,1,0.35);
  padFD2->SetTopMargin(0);
  padFD2->Draw();
  padFD2->cd();
  TH1F* hEffFDRatio = (TH1F*)hEffFD2->Clone();
  hEffFDRatio->SetDirectory(0);
  hEffFDRatio->Divide(hEffFD2,hEffFD1,1.,1.,"");
  for(Int_t iPt=0; iPt<hEffFDRatio->GetNbinsX(); iPt++) {
    hEffFDRatio->SetBinError(iPt+1,1.e-10);
  }
  hEffFDRatio->GetYaxis()->SetTitle("Ratio");
  hEffFDRatio->GetYaxis()->SetTitleSize(0.10);
  hEffFDRatio->GetXaxis()->SetTitleSize(0.10);
  hEffFDRatio->GetYaxis()->SetLabelSize(0.10);
  hEffFDRatio->GetXaxis()->SetLabelSize(0.10);
  hEffFDRatio->GetXaxis()->SetTitleOffset(1.2);
  hEffFDRatio->GetYaxis()->SetTitleOffset(0.6);
  hEffFDRatio->GetYaxis()->SetDecimals(2);
  hEffFDRatio->Draw();
  
  cEffPrompt->SaveAs("Comparison_EffPrompt.eps");
  cEffFD->SaveAs("Comparison_EffFD.eps");
  
}

void FPromptComparison(TString filename1="fprompt_unbinned_sigmafree.root",
                       TString filename2="fprompt_unbinned_sigmafree_d0cut.root",
                       TString name1="w/o d_{0} cut", TString name2="|d_{0}| < 100 #mum") {
  
  
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.3);
  
  TFile file1(filename1,"READ");
  TH1F* hFPrompt1 = (TH1F*)file1.Get("hFrac");
  hFPrompt1->SetDirectory(0);
  file1.Close();
  
  TFile file2(filename2,"READ");
  TH1F* hFPrompt2 = (TH1F*)file2.Get("hFrac");
  hFPrompt2->SetDirectory(0);
  file2.Close();
  
  
  TLegend* l = new TLegend(0.59,0.05,0.89,0.3);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.06);
  l->AddEntry(hFPrompt1,name1,"lpe");
  l->AddEntry(hFPrompt2,name2,"lpe");
  
  TCanvas* cFPrompt = new TCanvas("cFPrompt","",10,10,600,800);
  cFPrompt->cd();
  TPad *padPrompt1 = new TPad("padPrompt1","padPrompt1",0,0.35,1,1);
  padPrompt1->SetBottomMargin(0);
  padPrompt1->Draw();
  padPrompt1->cd();
  hFPrompt1->GetYaxis()->SetRangeUser(0.51,1.2);
  hFPrompt1->GetYaxis()->SetTitleSize(0.08);
  hFPrompt1->GetYaxis()->SetLabelSize(0.06);
  hFPrompt1->GetYaxis()->SetTitleOffset(0.8);
  hFPrompt1->SetLineColor(kBlue);
  hFPrompt1->SetMarkerColor(kBlue);
  hFPrompt2->SetLineColor(kRed);
  hFPrompt2->SetMarkerColor(kRed);
  hFPrompt2->SetMarkerStyle(21);
  hFPrompt1->Draw();
  hFPrompt2->Draw("same");
  l->Draw("same");
  cFPrompt->cd();
  TPad *padPrompt2 = new TPad("padPrompt2","padPrompt2",0.,0.,1,0.35);
  padPrompt2->SetTopMargin(0);
  padPrompt2->Draw();
  padPrompt2->cd();
  TH1F* hFPromptRatio = (TH1F*)hFPrompt2->Clone();
  hFPromptRatio->SetDirectory(0);
  hFPromptRatio->Divide(hFPrompt2,hFPrompt1,1.,1.,"");
  for(Int_t iPt=0; iPt<hFPromptRatio->GetNbinsX(); iPt++) {
    hFPromptRatio->SetBinError(iPt+1,1.e-10);
  }
  hFPromptRatio->GetYaxis()->SetTitle("Ratio");
  hFPromptRatio->GetYaxis()->SetTitleSize(0.10);
  hFPromptRatio->GetXaxis()->SetTitleSize(0.10);
  hFPromptRatio->GetYaxis()->SetLabelSize(0.10);
  hFPromptRatio->GetXaxis()->SetLabelSize(0.10);
  hFPromptRatio->GetXaxis()->SetTitleOffset(1.2);
  hFPromptRatio->GetYaxis()->SetTitleOffset(0.6);
  hFPromptRatio->GetYaxis()->SetDecimals(2);
  hFPromptRatio->Draw();
  
  cFPrompt->SaveAs("Comparison_Fprompt.eps");
  
}

void CrossSectionComparison(TString filename1="DplusCrossSection.root",
                            TString filename2="DplusCrossSection_d0cut.root",
                            TString name1="w/o d_{0} cut", TString name2="|d_{0}| < 100 #mum")) {
  
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.3);
  
  TFile file1(filename1,"READ");
  TH1F* hCrossPrompt1 = (TH1F*)file1.Get("hPromptCrossSection");
  hCrossPrompt1->SetDirectory(0);
  file1.Close();
  
  TFile file2(filename2,"READ");
  TH1F* hCrossPrompt2 = (TH1F*)file2.Get("hPromptCrossSection");
  hCrossPrompt2->SetDirectory(0);
  file2.Close();
  
  TLegend* l = new TLegend(0.55,0.6,0.85,0.85);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.06);
  l->AddEntry(hCrossPrompt1,name1,"lpe");
  l->AddEntry(hCrossPrompt2,name2,"lpe");
  
  TCanvas* cCrossPrompt = new TCanvas("cCrossPrompt","",10,10,600,800);
  cCrossPrompt->cd();
  TPad *padPrompt1 = new TPad("padPrompt1","padPrompt1",0,0.35,1,1);
  padPrompt1->SetBottomMargin(0);
  padPrompt1->Draw();
  padPrompt1->cd();
  padPrompt1->SetLogy();
  hCrossPrompt1->GetYaxis()->SetRangeUser(hCrossPrompt1->GetMinimum()*0.7,hCrossPrompt1->GetMaximum());
  hCrossPrompt1->GetYaxis()->SetTitleSize(0.08);
  hCrossPrompt1->GetYaxis()->SetLabelSize(0.06);
  hCrossPrompt1->GetYaxis()->SetTitleOffset(0.8);
  hCrossPrompt1->SetLineColor(kBlue);
  hCrossPrompt1->SetMarkerColor(kBlue);
  hCrossPrompt2->SetLineColor(kRed);
  hCrossPrompt2->SetMarkerColor(kRed);
  hCrossPrompt2->SetMarkerStyle(21);
  hCrossPrompt1->Draw();
  hCrossPrompt2->Draw("same");
  l->Draw("same");
  cCrossPrompt->cd();
  TPad *padPrompt2 = new TPad("padPrompt2","padPrompt2",0.,0.,1,0.35);
  padPrompt2->SetTopMargin(0);
  padPrompt2->Draw();
  padPrompt2->cd();
  TH1F* hCrossPromptRatio = (TH1F*)hCrossPrompt2->Clone();
  hCrossPromptRatio->SetDirectory(0);
  hCrossPromptRatio->Divide(hCrossPrompt2,hCrossPrompt1,1.,1.,"");
  for(Int_t iPt=0; iPt<hCrossPromptRatio->GetNbinsX(); iPt++) {
    hCrossPromptRatio->SetBinError(iPt+1,1.e-10);
  }
  hCrossPromptRatio->GetYaxis()->SetTitle("Ratio");
  hCrossPromptRatio->GetYaxis()->SetTitleSize(0.10);
  hCrossPromptRatio->GetXaxis()->SetTitleSize(0.10);
  hCrossPromptRatio->GetYaxis()->SetLabelSize(0.10);
  hCrossPromptRatio->GetXaxis()->SetLabelSize(0.10);
  hCrossPromptRatio->GetXaxis()->SetTitleOffset(1.2);
  hCrossPromptRatio->GetYaxis()->SetTitleOffset(0.6);
  hCrossPromptRatio->GetYaxis()->SetDecimals(2);
  hCrossPromptRatio->Draw();
  
  cCrossPrompt->SaveAs("Comparison_PromptCrossSection.eps");
  
}

