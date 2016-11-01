const Int_t nPtBins = 7;
const Int_t nPtLims = nPtBins+1;
Double_t PtLims[nPtLims] = {2,3,4,5,6,8,12,16};

TCanvas* PlotFinalFit(Double_t xLatex=80, Double_t yLatex=100, Double_t ymin=0.1, Double_t ymax=300, Double_t pTmin = 5, Double_t pTmax=6, TString filename="FitUnbinned_5-6");
void PlotAllFits();

void PlotAllFits() {

  TCanvas* cAll = new TCanvas("cAll","",800,1000);
  cAll->Divide(2,nPtBins/2+1);
  
  Double_t xLatex=70;
  Double_t yLatex=300;
  Double_t ymin=0.1;
  Double_t ymax=400; 

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    if(PtLims[iPt]==8) {
      xLatex=55;
      yLatex=160;
      ymin=0.01;
      ymax=150;
    }
    else if(PtLims[iPt]==12) {
      xLatex=30;
      yLatex=160;
      ymin=0.01;
      ymax=150;
    }
    else {
      xLatex=70;
      yLatex=550;
      ymin=0.1;
      ymax=400;
    }
    
    TCanvas* c=(TCanvas*)PlotFinalFit(xLatex,yLatex,ymin,ymax,PtLims[iPt],PtLims[iPt+1],Form("FitUnbinned_%0.f-%0.f_bkg",PtLims[iPt],PtLims[iPt+1]));
    cAll->cd(iPt+1);
    c->DrawClonePad();
  }

  cAll->SaveAs("FinalFits.eps");
  delete cAll;
}

TCanvas* PlotFinalFit(Double_t xLatex, Double_t yLatex, Double_t ymin, Double_t ymax, Double_t pTmin, Double_t pTmax, TString filename) {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);

  TFile infile(Form("%s.root",filename.Data()),"UPDATE");
  TH1F* hImpPar = (TH1F*)infile.Get("hImpPar");
  TF1* fImpParTot = (TF1*)infile.Get("ImpParTotFunc");
  TF1* fImpParPrompt = (TF1*)infile.Get("ImpParPromptFunc");
  TF1* fImpParRecoFD = (TF1*)infile.Get("ImpParRecoFDFunc");
  TF1* fImpParBkg = (TF1*)infile.Get("ImpParBkgFunc");
  hImpPar->SetDirectory(0);
  infile.Close();
  fImpParPrompt->SetLineStyle(7);
  fImpParRecoFD->SetLineStyle(10);
  fImpParBkg->SetLineStyle(9);
  fImpParTot->SetLineWidth(4);
  fImpParPrompt->SetLineWidth(4);
  fImpParRecoFD->SetLineWidth(4);
  fImpParBkg->SetLineWidth(4);
  fImpParPrompt->SetLineStyle(7);
  fImpParRecoFD->SetLineStyle(10);
  fImpParBkg->SetLineStyle(9);
  fImpParTot->SetNpx(200);
  fImpParPrompt->SetNpx(200);
  fImpParRecoFD->SetNpx(200);
  fImpParBkg->SetNpx(200);
  
  TFile fracfile("/home/fabrizio/tesi/ImpParAnalysis/analysis/pPb_LHC13/result/fprompt.root","UPDATE");
  TH1F* hFPrompt = (TH1F*)fracfile.Get("hPromptFraction");
  hFPrompt->SetDirectory(0);

  Int_t iPt = 0;
  for(Int_t iBin=0; iBin<hFPrompt->GetNbinsX(); iBin++) {
    if(hFPrompt->GetBinLowEdge(iBin+1)==pTmin)
      iPt=iBin;
  }
  cout << iPt << endl;
  Double_t fraction = hFPrompt->GetBinContent(iPt+1);
  Double_t fracerr = hFPrompt->GetBinError(iPt+1);
  
  TLatex lat;
  lat.SetTextFont(132);
  lat.SetTextSize(0.05);

  TLegend* l = new TLegend(0.14,0.66,0.45,0.89);
  l->SetTextFont(132);
  l->SetTextSize(0.05);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetFillColor(kWhite);
  l->AddEntry(hImpPar, "Projected histogram", "lpe");
  l->AddEntry(fImpParTot, "Total fit function", "l");
  l->AddEntry(fImpParPrompt, "Prompt", "l");
  l->AddEntry(fImpParRecoFD, "Feed-down", "l");
  l->AddEntry(fImpParBkg, "Background", "l");

  TCanvas* c = new TCanvas("c","",1000,800);
  c->SetLogy();
  hImpPar->GetYaxis()->SetRangeUser(ymin,ymax);
  hImpPar->GetXaxis()->SetRangeUser(-400,400);
  hImpPar->GetYaxis()->SetTitleSize(0.05);
  hImpPar->GetXaxis()->SetTitleSize(0.05);
  hImpPar->GetYaxis()->SetLabelSize(0.05);
  hImpPar->GetXaxis()->SetLabelSize(0.05);
  hImpPar->Draw("E1");
  hImpPar->GetXaxis()->SetTitle("Imp Par XY (#mum)");
  hImpPar->GetYaxis()->SetTitle(Form("Entries/(%0.f #mum)",hImpPar->GetBinWidth(10)));
  hImpPar->GetYaxis()->SetRangeUser(hImpPar->GetMinimum()+0.01,hImpPar->GetMaximum()*5);
  fImpParTot->Draw("same");
  fImpParPrompt->Draw("same");
  fImpParRecoFD->Draw("same");
  fImpParBkg->Draw("same");
  lat.DrawLatex(xLatex, yLatex,Form("%0.f < #it{p}_{T} < %0.f GeV/c", pTmin, pTmax));
  lat.DrawLatex(xLatex, yLatex-0.6*yLatex,Form("#it{f}_{prompt} = %.2f #pm %.2f", fraction, fracerr));
  l->Draw("same");
  
  c->SaveAs(Form("%s_plot.eps",filename.Data()));

  return c;
}
