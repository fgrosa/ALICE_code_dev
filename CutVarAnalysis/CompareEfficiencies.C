void CompareEfficiencies(Int_t iSetPrompt=1, Int_t iSetFD=3) {

  gStyle->SetOptStat(0);
  
  TString titlePrompt="max. prompt set";
  TString titleFD="max. prompt set";
  if(iSetPrompt==2)
    titlePrompt = "mixed set";
  if(iSetPrompt==3)
    titlePrompt = "max. feed-down set";
  if(iSetFD==2)
    titleFD = "mixed set";
  if(iSetFD==3)
    titleFD = "max. feed-down set";

  TFile rewfile("RewEff/Efficiencies.root","UPDATE");
  TH1F* hEffPromptRew = (TH1F*)rewfile.Get(Form("hEffPrompt_Set%d",iSetPrompt));
  TH1F* hEffFDRew = (TH1F*)rewfile.Get(Form("hEffFD_Set%d",iSetFD));
  hEffPromptRew->SetDirectory(0);
  hEffFDRew->SetDirectory(0);
  hEffPromptRew->SetLineColor(kGreen+3);
  hEffFDRew->SetLineColor(kMagenta);
  hEffPromptRew->SetMarkerColor(kGreen+3);
  hEffFDRew->SetMarkerColor(kMagenta);
  hEffPromptRew->SetLineWidth(2);
  hEffFDRew->SetLineWidth(2);
  hEffPromptRew->SetLineWidth(2);
  hEffFDRew->SetLineWidth(2);
  hEffPromptRew->SetMarkerSize(1.5);
  hEffFDRew->SetMarkerSize(1.5);
  hEffPromptRew->SetMarkerStyle(20);
  hEffFDRew->SetMarkerStyle(21);
  rewfile.Close();
  
  TFile notrewfile("NotRewEff/Efficiencies.root","UPDATE");
  TH1F* hEffPromptNotRew = (TH1F*)notrewfile.Get(Form("hEffPrompt_Set%d",iSetPrompt));
  TH1F* hEffFDNotRew = (TH1F*)notrewfile.Get(Form("hEffFD_Set%d",iSetFD));
  hEffPromptNotRew->SetDirectory(0);
  hEffPromptNotRew->SetLineWidth(2);
  hEffPromptNotRew->SetLineColor(kBlue);
  hEffFDNotRew->SetDirectory(0);
  hEffFDNotRew->SetLineWidth(2);
  hEffFDNotRew->SetLineColor(kRed);
  hEffPromptNotRew->SetMarkerSize(1.5);
  hEffFDNotRew->SetMarkerSize(1.5);
  hEffPromptNotRew->SetMarkerStyle(20);
  hEffFDNotRew->SetMarkerStyle(21);
  hEffPromptNotRew->SetMarkerColor(kBlue);
  hEffFDNotRew->SetMarkerColor(kRed);
  notrewfile.Close();
  
  TLine* line = new TLine(2,1,16,1);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->SetLineColor(kBlack);
  
  TH1F* hRatioPrompt = (TH1F*)hEffPromptNotRew->Clone();
  hRatioPrompt->Divide(hEffPromptRew,hEffPromptNotRew,1.,1.);

  TH1F* hRatioFD = (TH1F*)hEffFDNotRew->Clone();
  hRatioFD->Divide(hEffFDRew,hEffFDNotRew,1.,1.);

  hEffPromptRew->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hEffFDRew->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hEffPromptNotRew->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hEffFDNotRew->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
 
  hEffFDNotRew->GetXaxis()->SetTitle(hEffPromptNotRew->GetXaxis()->GetTitle());
  hRatioPrompt->GetXaxis()->SetTitle(hEffPromptNotRew->GetXaxis()->GetTitle());
  hRatioFD->GetXaxis()->SetTitle(hEffPromptNotRew->GetXaxis()->GetTitle());
  hRatioPrompt->SetTitle(titlePrompt.Data());
  hRatioFD->SetTitle(titleFD.Data());

  hRatioPrompt->GetXaxis()->SetTitleSize(0.05);
  hRatioPrompt->GetYaxis()->SetTitleSize(0.05);
  hRatioFD->GetXaxis()->SetTitleSize(0.05);
  hRatioFD->GetYaxis()->SetTitleSize(0.05);
  hRatioPrompt->GetXaxis()->SetLabelSize(0.05);
  hRatioPrompt->GetYaxis()->SetLabelSize(0.05);
  hRatioFD->GetXaxis()->SetLabelSize(0.05);
  hRatioFD->GetYaxis()->SetLabelSize(0.05);

  for(Int_t iPt=0; iPt<hRatioFD->GetNbinsX(); iPt++) {
    hRatioFD->SetBinError(iPt+1,1.e-10);
    hRatioPrompt->SetBinError(iPt+1,1.e-10);
  }

  hEffPromptNotRew->GetYaxis()->SetTitle("#epsilon_{Prompt}");
  hEffFDNotRew->GetYaxis()->SetTitle("#epsilon_{FD}");
  hRatioPrompt->GetYaxis()->SetTitle("Ratio of #epsilon_{Prompt}");
  hRatioFD->GetYaxis()->SetTitle("Ratio of #epsilon_{FD}");
  hRatioPrompt->GetYaxis()->SetTitleOffset(1.5);
  hRatioFD->GetYaxis()->SetTitleOffset(1.5);

  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.13);

  TLegend *lPrompt = new TLegend(0.55,0.15,0.75,0.35);
  if(iSetPrompt==3) {
    lPrompt->SetY1(0.7);
    lPrompt->SetY2(0.89);
    lPrompt->SetX1(0.2);
    lPrompt->SetX2(0.45);
  }
  lPrompt->SetTextSize(0.045);
  lPrompt->SetBorderSize(0);
  lPrompt->AddEntry(hEffPromptNotRew, "w/o #it{p}_{T} reweight","lpe");
  lPrompt->AddEntry(hEffPromptRew, "with #it{p}_{T} reweight","lpe");

  TLegend *lPromptRatio = new TLegend(0.4,0.7,0.75,0.89); 
  lPromptRatio->SetTextSize(0.045);
  lPromptRatio->SetBorderSize(0);
  lPromptRatio->AddEntry(hRatioPrompt, "#frac{FONLL}{MC}","lpe");

  TLegend *lFD = new TLegend(0.55,0.15,0.75,0.35); 
  if(iSetFD==3) {
    lFD->SetY1(0.2);
    lFD->SetY2(0.4);
  }
  lFD->SetTextSize(0.05);
  lFD->SetBorderSize(0);
  lFD->AddEntry(hEffFDNotRew, "w/o #it{p}_{T} reweight","lpe");
  lFD->AddEntry(hEffFDRew, "with #it{p}_{T} reweight","lpe");
                                                   
  TLegend *lFDRatio = new TLegend(0.4,0.7,0.75,0.89); 
  lFDRatio->SetTextSize(0.045);
  lFDRatio->SetBorderSize(0);
  lFDRatio->AddEntry(hRatioFD, "#frac{FONLL}{MC}","lpe");

  TLegend* lComb = new TLegend(0.2,0.7,0.6,0.89);
  lComb->SetTextSize(0.045);
  lComb->SetBorderSize(0);
  lComb->AddEntry(hRatioPrompt, Form("Prompt - %s",titlePrompt.Data()),"lpe");
  lComb->AddEntry(hRatioFD, Form("Feed-down - %s",titleFD.Data()),"lpe");
  
  TCanvas* cPrompt = new TCanvas("cPrompt","",800,800);
  cPrompt->SetLogy();
  hEffPromptNotRew->Draw();
  hEffPromptRew->Draw("same");
  lPrompt->Draw("same");
  
  TCanvas* cPromptRatio = new TCanvas("cPromptRatio","",800,800);
  hRatioPrompt->GetYaxis()->SetRangeUser(0.96,1.08);
  if(iSetPrompt==2)
    hRatioPrompt->GetYaxis()->SetRangeUser(0.94,1.06);  
  if(iSetPrompt==3)
    hRatioPrompt->GetYaxis()->SetRangeUser(0.94,1.04);
  hRatioPrompt->Draw("E");
  lPromptRatio->Draw("same");
  line->Draw("same");

  TCanvas* cFD = new TCanvas("cFD","",800,800);
  cFD->SetLogy();
  hEffFDNotRew->Draw();
  hEffFDRew->Draw("same");
  lFD->Draw("same");

  TCanvas* cFDRatio = new TCanvas("cFDRatio","",800,800);
  hRatioFD->GetYaxis()->SetRangeUser(0.97,1.05);
  hRatioFD->Draw("E");
  lFDRatio->Draw("same");
  line->Draw("same");
  
  cPrompt->SaveAs(Form("EffPromptCompSet%d.eps",iSetPrompt));
  cPromptRatio->SaveAs(Form("EffPromptRatioSet%d.eps",iSetPrompt));
  cFD->SaveAs(Form("EffFDCompSet%d.eps",iSetFD));
  cFDRatio->SaveAs(Form("EffFDRatioSet%d.eps",iSetFD));

  TCanvas* cComb = new TCanvas("cComb","",800,800);
  hRatioFD->SetTitle("");
  hRatioFD->GetYaxis()->SetRangeUser(0.915,1.125);
  hRatioFD->GetYaxis()->SetTitle("Efficiency(FONLL)/Efficiency(MC)");
  hRatioFD->Draw("E");
  hRatioPrompt->Draw("Esame");
  lComb->Draw("same");
  line->Draw("same");

  cComb->SaveAs(Form("EffRatioCombined_Set_%d-%d.eps",iSetPrompt,iSetFD));
}
