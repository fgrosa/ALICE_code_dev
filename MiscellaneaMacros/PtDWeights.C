Double_t PtWeightsFromFONLL5overLHC13d3(Double_t *x, Double_t *pars);
void CreatePtDWeightsHisto(Int_t nPtBins=80, Double_t PtMin=0, Double_t PtMax=40);

void CreatePtDWeightsHisto(Int_t nPtBins, Double_t PtMin, Double_t PtMax) {

  Int_t nPars=9;
  TF1* fFuncWeight = new TF1("fFuncWeight",PtWeightsFromFONLL5overLHC13d3,PtMin,PtMax,nPars);
  fFuncWeight->SetParameters(2.94999e+00,3.47032e+00,2.81278e+00,2.5,1.93370e-02,3.86865e+00,-1.54113e-01,8.86944e-02,2.56267e-02);

  TH1F* hWeights = new TH1F("hWeights","",nPtBins,PtMin,PtMax);
  Double_t PtBinWidth = hWeights->GetBinWidth(1);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    cout << iPt*PtBinWidth+PtBinWidth/2 << endl;
    hWeights->SetBinContent(iPt+1,fFuncWeight->Eval(iPt*PtBinWidth+PtBinWidth/2));
  }
  
  TCanvas *cPtDWeights = new TCanvas("cPtDWeights","",1200,900);
  cPtDWeights->Clear();
  hWeights->SetStats(0);
  hWeights->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hWeights->GetYaxis()->SetTitle("#frac{FONLL}{GenLimAcc}");
  hWeights->GetXaxis()->SetTitleOffset(1.3);
  hWeights->GetYaxis()->SetTitleOffset(1.3);
  hWeights->Draw();
  fFuncWeight->Draw("same");

  cPtDWeights->SaveAs("PtDWeights.eps");
  
  TFile outfile("PtD.root","RECREATE");
  hWeights->Write();
  fFuncWeight->Write();
  outfile.Close();
  
}

//_________________________________________________________________________
Double_t PtWeightsFromFONLL5overLHC13d3(Double_t *x, Double_t *pars){
  // weight function from the ratio of the LHC13d3 MC
  // and FONLL calculations for pp data
  Double_t pt = x[0];
  Double_t weight = (pars[0]*pt)/TMath::Power(pars[2],(1+TMath::Power(pars[3],pt/pars[1])))+pars[4]*TMath::Exp(pars[5]+pars[6]*pt)+pars[7]*TMath::Exp(pars[8]*pt);
  
  return weight;
}
