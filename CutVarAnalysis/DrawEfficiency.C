void DrawEfficiency() {
    
  gStyle->SetOptStat(1);
  gStyle->SetTitleOffset(1.1,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetTitleSize(0.06,"t");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.14);
  
  //____________________________________________________________________________
  Int_t nPtBins = 7;
  const Int_t nPtLims = nPtBins+1;
  Double_t PtLims[nPtLims] = {2,3,4,5,6,8,12,16};
  //____________________________________________________________________________
    
  TH1F** hEffPrompt = new TH1F*[3];
  TH1F** hEffFD = new TH1F*[3];
  for(Int_t iSet=0; iSet<3; ++iSet) {
    hEffPrompt[iSet] = new TH1F(Form("hEffPromptSet_%d",iSet+1),"",nPtBins,PtLims);
    hEffFD[iSet] = new TH1F(Form("hEffFDSet_%d",iSet+1),"",nPtBins,PtLims);
  }
      
  for(Int_t iBin=0; iBin<nPtBins; ++iBin) {
    TFile infile(Form("eff_%d.root",iBin),"READ");
    TCanvas* c=(TCanvas*)infile.Get("cEff");
    TH1F* hPrompt=(TH1F*)c->GetPrimitive(Form("hEffPrompt_%d.0-%d.0",(Int_t)PtLims[iBin],(Int_t)PtLims[iBin+1]));
    TH1F* hFD=(TH1F*)c->GetPrimitive(Form("hEffFD_%d.0-%d.0",(Int_t)PtLims[iBin],(Int_t)PtLims[iBin+1]));
    for(Int_t iSet=0; iSet<3; ++iSet) {
      hEffPrompt[iSet]->SetBinContent(iBin+1,hPrompt->GetBinContent(iSet+1));
      hEffFD[iSet]->SetBinContent(iBin+1,hFD->GetBinContent(iSet+1));
      hEffPrompt[iSet]->SetBinError(iBin+1,hPrompt->GetBinError(iSet+1));
      hEffFD[iSet]->SetBinError(iBin+1,hFD->GetBinError(iSet+1));
    }
    infile.Close();
  }

  TLegend* l = new TLegend(0.2,0.7,0.59,0.89);
  l->SetTextSize(0.06);
  l->AddEntry(hEffPrompt[0],"Prompt","lpe");
  l->AddEntry(hEffFD[0],"Feed-down","lpe");  
  TLegend* l2 = new TLegend(0.4,0.25,0.79,0.44);
  l2->SetTextSize(0.06);
  l2->AddEntry(hEffPrompt[0],"Prompt","lpe");
  l2->AddEntry(hEffFD[0],"Feed-down","lpe");  
  TLegend *l3 = new TLegend(0.2,0.7,0.59,0.89);
  l3->SetTextSize(0.06);
  TLegend *l4 = new TLegend(0.4,0.25,0.79,0.44);
  l4->SetTextSize(0.06);
  
  TCanvas **cEff = new TCanvas*[3];
  for(Int_t iSet=0; iSet<3; ++iSet) {
    cout << "\n************************* Eff prompt SET " << iSet << " ************************* " << endl;
    for(Int_t iBin=0; iBin<nPtBins; iBin++) {
      cout << hEffPrompt[iSet]->GetBinContent(iBin+1)<< " +/- " << hEffPrompt[iSet]->GetBinError(iBin+1) <<" err rel " << hEffPrompt[iSet]->GetBinError(iBin+1)/hEffPrompt[iSet]->GetBinContent(iBin+1)<<endl;
    }
  }
  
  for(Int_t iSet=0; iSet<3; ++iSet) {
    cout << "\n************************* Eff FD SET " << iSet << " ************************* " << endl;
    for(Int_t iBin=0; iBin<nPtBins; iBin++) {
      cout << hEffFD[iSet]->GetBinContent(iBin+1)<< " +/- " << hEffFD[iSet]->GetBinError(iBin+1)<<" err rel " << hEffFD[iSet]->GetBinError(iBin+1)/hEffFD[iSet]->GetBinContent(iBin+1)<<endl;
    }
  }

  TCanvas* cEff3Set = new TCanvas("cEff3Set","",1200,600);
  cEff3Set->Divide(3,1);

  TString setname[3] = {"Prompt-enhanced","Mixed","Feed-down-enhanced"};
  
  for(Int_t iSet=0; iSet<3; ++iSet) {
    cout << "\n************************* Ratios  SET " << iSet << " ************************* " << endl;
    for(Int_t iBin=0; iBin<nPtBins; iBin++) {
      cout << hEffPrompt[iSet]->GetBinContent(iBin+1)/hEffFD[iSet]->GetBinContent(iBin+1)<< "  " << hEffFD[iSet]->GetBinContent(iBin+1)/hEffPrompt[iSet]->GetBinContent(iBin+1)<<endl;
    }
    
    cEff[iSet] = new TCanvas(Form("cEff_Set%d",iSet+1),"",1200,900);
    cEff[iSet]->cd();
    cEff[iSet]->Clear();
    cEff[iSet]->SetLogy();
    Double_t ymax=(hEffPrompt[iSet]->GetMaximum()>hEffFD[iSet]->GetMaximum()) ? hEffPrompt[iSet]->GetMaximum() : hEffFD[iSet]->GetMaximum();
    Double_t ymin=(hEffPrompt[iSet]->GetMinimum()<hEffFD[iSet]->GetMinimum()) ? hEffPrompt[iSet]->GetMinimum() : hEffFD[iSet]->GetMinimum();
    if(1.5*ymax>1) ymax=1/1.5;
    hEffPrompt[iSet]->SetStats(0);
    hEffPrompt[iSet]->GetYaxis()->SetRangeUser(5.e-4,1.);
    hEffPrompt[iSet]->SetTitle(Form("%s",setname[iSet].Data()));
    hEffPrompt[iSet]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hEffPrompt[iSet]->GetYaxis()->SetTitle("Efficiency");
    hEffPrompt[iSet]->GetXaxis()->SetTitleSize(0.06);
    hEffPrompt[iSet]->GetYaxis()->SetTitleSize(0.06);
    hEffPrompt[iSet]->GetXaxis()->SetLabelSize(0.05);
    hEffPrompt[iSet]->GetYaxis()->SetLabelSize(0.05);
    hEffPrompt[iSet]->SetMarkerColor(kBlue);
    hEffPrompt[iSet]->SetMarkerStyle(20);
    hEffPrompt[iSet]->SetMarkerSize(1.5);
    hEffPrompt[iSet]->SetLineWidth(2);
    hEffFD[iSet]->SetStats(0);
    hEffFD[iSet]->SetLineColor(kRed);
    hEffFD[iSet]->SetMarkerColor(kRed);
    hEffFD[iSet]->SetMarkerStyle(21);
    hEffFD[iSet]->SetMarkerSize(1.5);
    hEffFD[iSet]->SetLineWidth(2);
    hEffPrompt[iSet]->Draw("E1");
    hEffFD[iSet]->Draw("E1same");
    if(iSet==2)
      l->Draw("same");
    else
      l2->Draw("same");
    
    cEff3Set->cd(iSet+1)->SetLogy();
    TH1F* hEffPromptCopy=(TH1F*)hEffPrompt[iSet]->Clone();
    hEffPromptCopy->SetMarkerSize(1.);
    TH1F* hEffFDCopy=(TH1F*)hEffFD[iSet]->Clone();
    hEffFDCopy->SetMarkerSize(1.);
    if(iSet==0) {
      l3->AddEntry(hEffPromptCopy,"Prompt","lpe");
      l3->AddEntry(hEffFDCopy,"Feed-down","lpe");
    }
    if(iSet==0) {
      l4->AddEntry(hEffPromptCopy,"Prompt","lpe");
      l4->AddEntry(hEffFDCopy,"Feed-down","lpe");
    }
    hEffPromptCopy->Draw();
    hEffFDCopy->Draw("same");    
    if(iSet==2)
      l3->Draw("same");
    else
      l4->Draw("same");
    
    cEff[iSet]->SaveAs(Form("Eff_Set%d.eps",iSet+1));
    cEff[iSet]->SaveAs(Form("Eff_Set%d.root",iSet+1));

  }
    cEff3Set->SaveAs("Eff3Sets.eps");

}
  
