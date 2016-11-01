void PlotBkgFunc(Double_t min=-750, Double_t max=750, Double_t pTmin = 2, Double_t pTmax=3) {

  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  
  TFile symmfile(Form("ImpParBkg_%0.f-%0.f.root",pTmin,pTmax),"UPDATE");
  TCanvas *cSymm=(TCanvas*)symmfile.Get("cBkg");
  TH1F* hImpParBkg=(TH1F*)cSymm->GetPrimitive("hImpParBkg");
  hImpParBkg->SetDirectory(0);
  TF1* fSymm=(TF1*)cSymm->GetPrimitive("ImpParBkgFunc");
  fSymm->SetLineColor(kRed);
  fSymm->SetLineWidth(3);
  symmfile.Close();

  TFile asymmfile(Form("ImpParBkgAsymm_%0.f-%0.f.root",pTmin,pTmax),"UPDATE");
  TCanvas *cAsymm=(TCanvas*)asymmfile.Get("cBkg");
  TF1* fAsymm=(TF1*)cAsymm->GetPrimitive("ImpParBkgFunc");
  fAsymm->SetLineColor(kBlue);
  fAsymm->SetLineWidth(3);
  asymmfile.Close();

  hImpParBkg->GetXaxis()->SetTitle("Imp Par XY (#mum)");
  hImpParBkg->GetYaxis()->SetTitle(Form("Entries/(%0.f #mum)",hImpParBkg->GetBinWidth(10)));
  hImpParBkg->GetXaxis()->SetTitleSize(0.05);
  hImpParBkg->GetYaxis()->SetTitleSize(0.05);
  hImpParBkg->GetXaxis()->SetLabelSize(0.05);
  hImpParBkg->GetYaxis()->SetLabelSize(0.05);
  hImpParBkg->GetXaxis()->SetTitleOffset(1.2);
  hImpParBkg->GetYaxis()->SetTitleOffset(1.2);
  hImpParBkg->GetXaxis()->SetNdivisions(508,kTRUE);

  TLegend* l = new TLegend(0.15,0.65,0.4,0.8);
  l->AddEntry(fSymm,"symmetric","l");
  l->AddEntry(fAsymm,"asymmetric","l");
  l->SetTextSize(0.045);
  l->SetFillStyle(0);

  TLatex* lat = new TLatex();
  lat->SetTextColor(kBlack);
  lat->SetTextSize(0.045);
  lat->SetTextFont(132);  
  TLatex* latS = new TLatex();
  latS->SetTextColor(kRed);
  latS->SetTextSize(0.04);
  latS->SetTextFont(132);
  TLatex* latA = new TLatex();
  latA->SetTextColor(kBlue);
  latA->SetTextSize(0.04);
  latA->SetTextFont(132);

  Double_t chiS = fSymm->GetChisquare();
  Int_t ndfS = fSymm->GetNDF();
  Double_t chiA = fAsymm->GetChisquare();
  Int_t ndfA = fAsymm->GetNDF();
  
  TCanvas *c = new TCanvas("c","",800,800);
  c->SetLogy();
  hImpParBkg->SetStats(0);
  hImpParBkg->GetXaxis()->SetRangeUser(min,max);
  hImpParBkg->Draw("E1");
  fSymm->Draw("same");
  fAsymm->Draw("same");
  l->Draw("same");
  lat->DrawLatex(-650,700,Form("%0.f < #it{p}_{T} < %0.f GeV/c",pTmin,pTmax));
  latS->DrawLatex(190,420,Form("#chi^{2}/ndf = %0.1f/%d",chiS,ndfS));
  latA->DrawLatex(190,170,Form("#chi^{2}/ndf = %0.1f/%d",chiA,ndfA));

  c->SaveAs("BkgFitFuncComp.eps");
}
