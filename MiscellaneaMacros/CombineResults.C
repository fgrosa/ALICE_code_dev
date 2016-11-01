const Int_t colors[] = {kBlack,kBlue,kGreen+3,kRed,kOrange+7};

void CombineCrossSection() {

  //________________________________________________________________________________________________
  //Input files
  
  TFile CutVarMinFile("CutVar/DplusCrossSection_CutVar_Min.root","READ");
  TFile CutVarIncFile("CutVar/DplusCrossSection_CutVar_Inc.root","READ");
  TFile ImpParFitFile("ImpPar/DplusCrossSection_ImpParFit.root","READ");
  TFile KFFile("KF/DplusCrossSectionKF_pPb5TeV.root","READ");
  TFile PublishedFile("published/DplusCrossSec_method2_fd2_br1.root","READ");
  
  //_________________________________________________________________________________________________
  //Cross sections
  Double_t BR = 0.0913;

  CutVarIncFile.cd();
  TH1F *hCrossCutVarInc=(TH1F*)CutVarIncFile.Get("hCrossSecPrompt");
  TGraphAsymmErrors *gCrossCutVarInc=(TGraphAsymmErrors*)CutVarIncFile.Get("gCrossSystPromptInc");
  hCrossCutVarInc->SetDirectory(0);
  CutVarIncFile.Close();
  
  CutVarMinFile.cd();
  TH1F *hCrossCutVarMin=(TH1F*)CutVarMinFile.Get("hCrossSecPrompt");
  TGraphAsymmErrors *gCrossCutVarMin=(TGraphAsymmErrors*)CutVarMinFile.Get("gCrossSystPromptMin");
  hCrossCutVarMin->SetDirectory(0);
  CutVarMinFile.Close();

  ImpParFitFile.cd();
  TH1F *hCrossImpParFit=(TH1F*)ImpParFitFile.Get("hPromptCrossSec");
  TGraphAsymmErrors *gCrossImpParFit=(TGraphAsymmErrors*)ImpParFitFile.Get("gPromptCrossSec");
  hCrossImpParFit->SetDirectory(0);
  ImpParFitFile.Close();

  KFFile.cd();
  TH1F *hCrossKF=(TH1F*)KFFile.Get("histoSigmaCorr");
  TGraphAsymmErrors *gCrossKF=(TGraphAsymmErrors*)KFFile.Get("gCrossSyst");
  hCrossKF->SetDirectory(0);
  KFFile.Close();
  
  PublishedFile.cd();
  TH1F* hCrossPub=(TH1F*)PublishedFile.Get("hAAC");
  TGraphAsymmErrors *gCrossPub=(TGraphAsymmErrors*)PublishedFile.Get("gaaCsystTot");
  hCrossPub->SetDirectory(0);
  PublishedFile.Close();

  hCrossImpParFit->SetMarkerColor(colors[1]);
  hCrossImpParFit->SetLineColor(colors[1]);
  hCrossImpParFit->SetLineWidth(2);  
  hCrossImpParFit->SetMarkerStyle(22);
  gCrossImpParFit->SetFillStyle(20);
  gCrossImpParFit->SetLineColor(colors[1]);
  gCrossImpParFit->SetLineWidth(2);
  hCrossCutVarMin->SetMarkerStyle(21);
  hCrossCutVarMin->SetMarkerColor(colors[2]);
  hCrossCutVarMin->SetLineColor(colors[2]);
  hCrossCutVarMin->SetLineWidth(2);  
  hCrossCutVarInc->SetMarkerStyle(23);
  hCrossCutVarInc->SetMarkerColor(colors[4]);
  hCrossCutVarInc->SetLineColor(colors[4]);
  hCrossCutVarInc->SetLineWidth(2);
  gCrossCutVarMin->SetFillStyle(20);
  gCrossCutVarMin->SetLineColor(colors[2]);
  gCrossCutVarMin->SetLineWidth(2);  
  gCrossCutVarInc->SetFillStyle(20);
  gCrossCutVarInc->SetLineColor(colors[4]);
  gCrossCutVarInc->SetLineWidth(2);  
  hCrossKF->SetMarkerStyle(33);
  hCrossKF->SetMarkerColor(colors[3]);
  hCrossKF->SetLineColor(colors[3]);
  hCrossKF->SetLineWidth(2);  
  gCrossKF->SetFillStyle(20);
  gCrossKF->SetLineColor(colors[3]);
  gCrossKF->SetLineWidth(2);
  hCrossPub->SetMarkerStyle(20);
  hCrossPub->SetMarkerColor(colors[0]);
  hCrossPub->SetLineColor(colors[0]);  
  hCrossPub->SetLineWidth(2);  
  gCrossPub->SetFillStyle(20);
  gCrossPub->SetLineColor(colors[0]);
  gCrossPub->SetLineWidth(2);
  hCrossPub->SetTitle("");  
  hCrossPub->GetXaxis()->SetTitleSize(0.05);
  hCrossPub->GetYaxis()->SetTitleSize(0.05);
  hCrossPub->GetXaxis()->SetLabelSize(0.05);
  hCrossPub->GetYaxis()->SetLabelSize(0.05);
  hCrossPub->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");  
  hCrossPub->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}} (#mub GeV/c)");  
  hCrossPub->GetYaxis()->SetTitleOffset(1.5);  
  hCrossPub->GetXaxis()->SetTitleOffset(1.);  

  //_________________________________________________________________________________________________
  //Ratios with respect to published (no error)
  const Int_t nPtBins = hCrossImpParFit->GetNbinsX();
  const Int_t nPtLims = nPtBins+1;
  TAxis* PtAxis = (TAxis*)hCrossImpParFit->GetXaxis();
  TArrayD* ptarray = (TArrayD*)PtAxis->GetXbins();
  Double_t* PtLims = (Double_t*)ptarray->GetArray();

  const Int_t nPtBinsPub = hCrossPub->GetNbinsX();
  const Int_t nPtLimsPub = nPtBins+1;
  TAxis* PtAxisPub = (TAxis*)hCrossPub->GetXaxis();
  TArrayD* ptarrayPub = (TArrayD*)PtAxisPub->GetXbins();
  Double_t* PtLimsPub = (Double_t*)ptarrayPub->GetArray();
  
  Bool_t IsBinningEqual = kTRUE;
  
  if(nPtBins!=nPtBinsPub) 
    IsBinningEqual = kFALSE;
  
  Int_t Ptcounter=0;
  while(Ptcounter<nPtBins && IsBinningEqual==kTRUE) {
    if(PtLims[Ptcounter]!=PtLimsPub[Ptcounter])
      IsBinningEqual = kFALSE;
    Ptcounter++;
  }
    
  TH1F* hRatioCutVarMin = (TH1F*)hCrossCutVarMin->Clone();
  TH1F* hRatioCutVarInc = (TH1F*)hCrossCutVarInc->Clone();
  TH1F* hRatioImpParFit = (TH1F*)hCrossImpParFit->Clone();
  TH1F* hRatioKF = (TH1F*)hCrossKF->Clone();
  hRatioKF->Divide(hCrossKF,hCrossPub,1.,1.);
  
  if(IsBinningEqual) {
    hRatioCutVarMin->Divide(hCrossCutVarMin,hCrossPub,1.,1.);
    hRatioCutVarInc->Divide(hCrossCutVarInc,hCrossPub,1.,1.);
    hRatioImpParFit->Divide(hCrossImpParFit,hCrossPub,1.,1.);
  }
  else {
    TH1F* hPubCrossSection2 = (TH1F*)hCrossPub->Clone();
    for(Int_t iPt=0; iPt<nPtBinsPub; iPt++) {
      hPubCrossSection2->SetBinContent(iPt+1,hPubCrossSection2->GetBinContent(iPt+1)*(PtLimsPub[iPt+1]-PtLimsPub[iPt]));
      hPubCrossSection2->SetBinError(iPt+1,hPubCrossSection2->GetBinError(iPt+1)*(PtLimsPub[iPt+1]-PtLimsPub[iPt]));
    }
    TH1F* hPubCrossSectionReb = (TH1F*)hPubCrossSection2->Rebin(nPtBins,"hPubCrossSectionReb",PtLims);
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      hPubCrossSectionReb->SetBinContent(iPt+1,hPubCrossSectionReb->GetBinContent(iPt+1)/(PtLims[iPt+1]-PtLims[iPt]));
      hPubCrossSectionReb->SetBinError(iPt+1,hPubCrossSectionReb->GetBinError(iPt+1)/(PtLims[iPt+1]-PtLims[iPt]));
    }
    hRatioCutVarMin->Divide(hCrossCutVarMin,hPubCrossSectionReb,1.,1.);
    hRatioCutVarInc->Divide(hCrossCutVarInc,hPubCrossSectionReb,1.,1.);
    hRatioImpParFit->Divide(hCrossImpParFit,hPubCrossSectionReb,1.,1.);
  }
  for(Int_t iPt=0; iPt<nPtBinsPub; iPt++) {
    hRatioCutVarMin->SetBinError(iPt+1,1.e-10);
    hRatioCutVarInc->SetBinError(iPt+1,1.e-10);
    hRatioImpParFit->SetBinError(iPt+1,1.e-10);
  }

  hRatioCutVarMin->SetDirectory(0);
  hRatioCutVarInc->SetDirectory(0);
  hRatioImpParFit->SetDirectory(0);
  hRatioKF->SetDirectory(0);

  //_________________________________________________________________________________________________
  //plots

  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetOptStat(0);

  hRatioCutVarMin->GetYaxis()->SetTitle("ratio w.r.t. the published result");
  hRatioCutVarMin->GetYaxis()->SetRangeUser(0.65,1.6);
  hRatioCutVarMin->GetXaxis()->SetTitleOffset(1.);  
  hRatioCutVarMin->GetYaxis()->SetTitleOffset(1.5);  

  TLegend* l = new TLegend(0.35,0.65,0.87,0.87);
  l->SetTextSize(0.04);
  l->SetFillStyle(0);
  l->AddEntry(hCrossPub,"Published Cross Section","lpe");
  l->AddEntry(hCrossImpParFit,"Impact-parameter fit","lpe");
  l->AddEntry(hCrossCutVarMin,"Cut-variation (analytic)","lpe");
  l->AddEntry(hCrossKF,"KFParticle","lpe");

  TLegend* lRatio = new TLegend(0.2,0.6,0.67,0.87);
  lRatio->SetTextSize(0.04);
  lRatio->SetFillStyle(0);
  lRatio->AddEntry(hRatioImpParFit,"#frac{Impact parameter}{Published}","lpe");
  lRatio->AddEntry(hRatioCutVarMin,"#frac{Minimisation}{Published}","lpe");
  lRatio->AddEntry(hRatioCutVarInc,"#frac{Incentre}{Published}","lpe");  

  TLine* line = new TLine(PtLims[0],1,PtLims[nPtBins],1);
  line->SetLineColor(kBlack);
  line->SetLineStyle(6);

  TLatex latex;
  latex.SetTextFont(132);
  latex.SetTextSize(0.04);

  TCanvas* cCross = new TCanvas("cCross","",800,800);
  cCross->SetLogy();
  hCrossPub->GetYaxis()->SetRangeUser(0.2,50000);
  hCrossPub->Draw("E");
  gCrossPub->Draw("2");
  hCrossImpParFit->Draw("Esame");
  gCrossImpParFit->Draw("2");
  hCrossCutVarMin->Draw("Esame");
  gCrossCutVarMin->Draw("2");
  hCrossKF->Draw("Esame");
  gCrossKF->Draw("2");
  l->Draw("same"); 
  latex->DrawLatex(PtLims[1]-0.5,7.5,"Prompt D^{+}");
  latex->DrawLatex(PtLims[1]-0.5,3,"pPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->DrawLatex(PtLims[1]-0.5,1.2,"#pm 2.1% BR unc. not shown");
  latex->DrawLatex(PtLims[1]-0.5,0.5,"#pm 3.7% norm. unc. not shown");
  
  TCanvas* cRatio = new TCanvas("cRatio","",800,600);
  hRatioCutVarMin->Draw("E");
  hRatioImpParFit->Draw("Esame");
  hRatioCutVarInc->Draw("Esame");
  lRatio->Draw("same");
  line->Draw("same");
  
  cCross->SaveAs("CrossSecComp.eps");
  cRatio->SaveAs("CrossSecRatios.eps");
}

