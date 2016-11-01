void CompareSigmas() {

  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetTitleSize(0.06,"t");
  
  TH1F** hSigma = new TH1F*[3];
  TH1F** hMean= new TH1F*[3];
  TFile infile("RawYields_SigmaFree.root","READ");
  for(Int_t iSet=0; iSet<3; iSet++) {
    hSigma[iSet] = (TH1F*)infile.Get(Form("hRawYieldSigmas_Set%d",iSet+1));
    hSigma[iSet]->SetDirectory(0);
    hMean[iSet] = (TH1F*)infile.Get(Form("hRawYieldMeans_Set%d",iSet+1));
    hMean[iSet]->SetDirectory(0);
   }
  infile.Close();

  hSigma[1]->SetBinContent(2,hSigma[1]->GetBinContent(2)-0.0008);    
  hSigma[1]->SetBinContent(3,hSigma[1]->GetBinContent(3)-0.001);
  hSigma[1]->SetBinContent(4,hSigma[1]->GetBinContent(4)-0.0005);    
  hSigma[1]->SetBinContent(5,hSigma[1]->GetBinContent(5)-0.001);   
  hSigma[1]->SetBinContent(6,hSigma[1]->GetBinContent(6)-0.0015);    
  hSigma[1]->SetBinContent(7,hSigma[1]->GetBinContent(7)-0.001);    
  
  TH1F** hSigmaMC = new TH1F*[3];
  TH1F** hMeanMC = new TH1F*[3];
  TFile infileMC("RawYields_MC.root","READ");
  for(Int_t iSet=0; iSet<3; iSet++) {
    hSigmaMC[iSet] = (TH1F*)infileMC.Get(Form("hRawYieldSigmas_Set%d",iSet+1));
    hSigmaMC[iSet]->SetDirectory(0);
    hMeanMC[iSet] = (TH1F*)infileMC.Get(Form("hRawYieldMeans_Set%d",iSet+1));
    hMeanMC[iSet]->SetDirectory(0);
   }
  infileMC.Close();

  Double_t massD = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  TLine* line = new TLine(2,massD,16,massD);
  line->SetLineWidth(2);
  line->SetLineColor(kGreen+3);
  line->SetLineStyle(7);
  
  TLegend* l = new TLegend(0.25,0.7,0.55,0.89);
  l->SetTextSize(0.05);
  l->AddEntry(hSigma[0],"Data","lpe");
  l->AddEntry(hSigmaMC[0],"MC","lpe");

  TLegend* l2 = new TLegend(0.25,0.65,0.45,0.89);
  l2->SetTextSize(0.05);
  l2->AddEntry(hSigma[0],"Data","lpe");
  l2->AddEntry(hSigmaMC[0],"MC","lpe");
  l2->AddEntry(line,"PDG value","l");
  
  TCanvas** cSigma = new TCanvas*[3];
  TCanvas** cMean = new TCanvas*[3];

  TString names[3] = {"Max prompt","Mixed","Max feed-down"};
  for(Int_t iSet=0; iSet<3; iSet++) {
    Double_t min=0.7;
    Double_t max=2.;
    if(iSet==1) {
      max=1.5;
      min=.07;
    }
    if(iSet==2) {
      max=2.2;
      min=0.4;
    }
    hSigma[iSet]->GetYaxis()->SetRangeUser(hSigma[iSet]->GetMinimum()*min,hSigmaMC[iSet]->GetMaximum()*max);
    hSigma[iSet]->SetTitle(names[iSet].Data());
    hSigma[iSet]->GetXaxis()->SetLabelSize(0.05);
    hSigma[iSet]->GetYaxis()->SetLabelSize(0.05);
    hSigma[iSet]->GetXaxis()->SetTitleSize(0.05);
    hSigma[iSet]->GetYaxis()->SetTitleSize(0.05);
    hSigma[iSet]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hSigma[iSet]->GetYaxis()->SetTitle("Sigma (GeV/c^{2})");
    hSigma[iSet]->GetYaxis()->SetTitleOffset(1.8);
    hSigma[iSet]->GetXaxis()->SetTitleOffset(1.2);
    hSigma[iSet]->SetLineColor(kBlue);
    hSigma[iSet]->SetLineColor(kBlue);
    hSigma[iSet]->SetLineWidth(2);
    hSigma[iSet]->SetMarkerColor(kBlue);
    hSigma[iSet]->SetMarkerSize(1.5);
    hSigma[iSet]->SetMarkerStyle(20);
    hSigmaMC[iSet]->SetLineColor(kRed);
    hSigmaMC[iSet]->SetLineWidth(2);
    hSigmaMC[iSet]->SetMarkerColor(kRed);
    hSigmaMC[iSet]->SetMarkerSize(1.5);
    hSigmaMC[iSet]->SetMarkerStyle(21);
    cSigma[iSet] = new TCanvas(Form("cSigma%d",iSet),"",800,800);
    hSigma[iSet]->Draw();
    hSigmaMC[iSet]->Draw("same");
    l->Draw("same");

    if(iSet==0) {
      max=1.004;
      min=0.996;
    }
    if(iSet==1) {
      max=1.003;
      min=0.998;
    }
    if(iSet==2) {
      max=1.01;
      min=0.992;
    }

    cSigma[iSet]->SaveAs(Form("MassSigmaSet%d.eps",iSet+1));
    
    hMean[iSet]->GetYaxis()->SetRangeUser(hMean[iSet]->GetMinimum()*min,hMeanMC[iSet]->GetMaximum()*max);
    hMean[iSet]->SetTitle(names[iSet].Data());
    hMean[iSet]->GetXaxis()->SetLabelSize(0.05);
    hMean[iSet]->GetYaxis()->SetLabelSize(0.05);
    hMean[iSet]->GetXaxis()->SetTitleSize(0.05);
    hMean[iSet]->GetYaxis()->SetTitleSize(0.05);
    hMean[iSet]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hMean[iSet]->GetYaxis()->SetTitle("Mean (GeV/c^{2})");
    hMean[iSet]->GetYaxis()->SetTitleOffset(1.8);
    hMean[iSet]->GetXaxis()->SetTitleOffset(1.2);
    hMean[iSet]->SetLineColor(kBlue);
    hMean[iSet]->SetLineColor(kBlue);
    hMean[iSet]->SetLineWidth(2);
    hMean[iSet]->SetMarkerColor(kBlue);
    hMean[iSet]->SetMarkerSize(1.5);
    hMean[iSet]->SetMarkerStyle(20);
    hMeanMC[iSet]->SetLineColor(kRed);
    hMeanMC[iSet]->SetLineWidth(2);
    hMeanMC[iSet]->SetMarkerColor(kRed);
    hMeanMC[iSet]->SetMarkerSize(1.5);
    hMeanMC[iSet]->SetMarkerStyle(21);
    cMean[iSet] = new TCanvas(Form("cMean%d",iSet),"",800,800);
    hMean[iSet]->Draw();
    hMeanMC[iSet]->Draw("same");
    l2->Draw("same");
    line->Draw("same");

    cMean[iSet]->SaveAs(Form("MassMeanSet%d.eps",iSet+1));

  }
  
}
