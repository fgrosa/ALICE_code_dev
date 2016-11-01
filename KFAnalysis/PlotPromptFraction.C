void PlotPromptFraction(TString cutset="cutset1") {

  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetLegendBorderSize(0.);
  
  TFile fcfile(Form("HFPtSpectrum_fc_%s.root",cutset.Data()),"UPDATE");
  TGraphAsymmErrors* gfc = (TGraphAsymmErrors*)fcfile.Get("gFcConservative");
  fcfile.Close();
  TFile Nbfile(Form("HFPtSpectrum_Nb_%s.root",cutset.Data()),"UPDATE");
  TGraphAsymmErrors* gNb = (TGraphAsymmErrors*)Nbfile.Get("gFcConservative");
  Nbfile.Close();
  TFile Totfile(Form("HFPtSpectrum_combinedFD_%s.root",cutset.Data()),"UPDATE");
  TGraphAsymmErrors* gTot = (TGraphAsymmErrors*)Totfile.Get("gFcCorrConservative");
  Totfile.Close();
  
  TLegend* l = new TLegend(0.6,0.25,0.8,0.5);
  l->SetTextSize(0.04);
  l->AddEntry(gfc,"#it{f_{c}} method","lpe");
  l->AddEntry(gNb,"#it{N_{b}} method","lpe");
  l->AddEntry(gTot,"combined #it{f_{c}}+#it{N_{b}}","f");
 
  TCanvas* c = new TCanvas("c","",1200,900);
  gfc->GetYaxis()->SetTitle("#it{f}_{prompt}");
  gfc->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  gfc->GetYaxis()->SetTitleSize(0.05);
  gfc->GetXaxis()->SetTitleSize(0.05);
  gfc->GetYaxis()->SetTitleFont(42);
  gfc->GetXaxis()->SetTitleFont(42);
  gfc->GetYaxis()->SetLabelFont(42);
  gfc->GetXaxis()->SetLabelFont(42);
  gfc->GetYaxis()->SetRangeUser(0.2,1.2);
  gfc->SetTitle();
  gfc->SetLineColor(kRed);
  gNb->SetLineColor(kBlue);
  gfc->SetMarkerColor(kRed);
  gNb->SetMarkerColor(kBlue);
  gfc->SetMarkerStyle(20);
  gNb->SetMarkerStyle(20);
  gTot->SetFillStyle(20);
//  gTot->SetFillColor(kGreen+3);
  gTot->SetLineColor(kGreen+3);
  gTot->SetLineWidth(1);

  gfc->Draw("AP");
  gNb->Draw("P");  
  gTot->Draw("2");
  l->Draw("same");
    
  c->SaveAs(Form("Fprompt_%s.eps",cutset.Data()));
  delete c;
}
