const Int_t nPtBins = 7;
const Int_t nPtLims = nPtBins+1;
Double_t PtLims[nPtLims] = {2,3,4,5,6,8,12,16};
Double_t xmin[nPtBins] = {-21100,-10000,-2700,-1200,-2000,-2000,-400};
Double_t xmax[nPtBins] = {35000,12000,6000,3000,3000,2000,500};
Double_t ymin[nPtBins] = {10100,10100,5100,4100,4100,1100,10};
Double_t ymax[nPtBins] = {70000,50000,25000,12000,10000,6000,1600};

void PlotLines();
void PlotLinesErr();
void PlotLinesErrDisp();
TCanvas* PlotLinesAndProjections(Int_t iPt=2, Double_t PtMin=4, Double_t PtMax=5, Double_t xmin=-2000, Double_t xmax=6000, Double_t ymin=5000, Double_t ymax=25000, Int_t type=0);

void PlotLines() {
  TCanvas* cAll = new TCanvas("cAll","",800,1000);
  cAll->Divide(2,nPtBins/2+1);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    TCanvas* c=(TCanvas*)PlotLinesAndProjections(iPt,PtLims[iPt],PtLims[iPt+1],xmin[iPt],xmax[iPt],ymin[iPt],ymax[iPt],0);
    cAll->cd(iPt+1);
    c->DrawClonePad();
  }

  cAll->SaveAs("Lines.eps");
  delete cAll;
}

void PlotLinesErr() {
  TCanvas* cAll = new TCanvas("cAll","",800,1000);
  cAll->Divide(2,nPtBins/2+1);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    TCanvas* c=(TCanvas*)PlotLinesAndProjections(iPt,PtLims[iPt],PtLims[iPt+1],xmin[iPt],xmax[iPt],ymin[iPt],ymax[iPt],1);
    cAll->cd(iPt+1);
    c->DrawClonePad();
  }

  cAll->SaveAs("LinesErr.eps");
  delete cAll;
}

void PlotLinesErrDisp() {
  TCanvas* cAll = new TCanvas("cAll","",800,1000);
  cAll->Divide(2,nPtBins/2+1);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    TCanvas* c=(TCanvas*)PlotLinesAndProjections(iPt,PtLims[iPt],PtLims[iPt+1],xmin[iPt],xmax[iPt],ymin[iPt],ymax[iPt],2);
    cAll->cd(iPt+1);
    c->DrawClonePad();
  }

  cAll->SaveAs("LinesErrDisp.eps");
  delete cAll;
}

TCanvas* PlotLinesAndProjections(Int_t iPt, Double_t PtMin, Double_t PtMax, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Int_t type) {
  
  gStyle->SetOptStat(1);
  gStyle->SetTitleOffset(1.,"y");
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetStatFontSize(0.07);
  TGaxis::SetMaxDigits(3);
  
  TFile infile(Form("Lines_%d.root",iPt),"READ");
  TCanvas* c=(TCanvas*)infile.Get("cLines");
  TH2F* hIncDist=(TH2F*)c->GetPrimitive("hIncDistError");
  TH1F* hProjX=(TH1F*)hIncDist->ProjectionX();
  TH1F* hProjY=(TH1F*)hIncDist->ProjectionY();
  TF1* lSet1=(TF1*)c->GetPrimitive(Form("lineSet1_PtBin%d",iPt));
  TF1* lSet2=(TF1*)c->GetPrimitive(Form("lineSet2_PtBin%d",iPt));
  TF1* lSet3=(TF1*)c->GetPrimitive(Form("lineSet3_PtBin%d",iPt));
  TF1* lSet1ErrLow=(TF1*)c->GetPrimitive(Form("lineErrLowSet1_PtBin%d",iPt));
  TF1* lSet2ErrLow=(TF1*)c->GetPrimitive(Form("lineErrLowSet2_PtBin%d",iPt));
  TF1* lSet3ErrLow=(TF1*)c->GetPrimitive(Form("lineErrLowSet3_PtBin%d",iPt));
  TF1* lSet1ErrHigh=(TF1*)c->GetPrimitive(Form("lineErrHighSet1_PtBin%d",iPt));
  TF1* lSet2ErrHigh=(TF1*)c->GetPrimitive(Form("lineErrHighSet2_PtBin%d",iPt));
  TF1* lSet3ErrHigh=(TF1*)c->GetPrimitive(Form("lineErrHighSet3_PtBin%d",iPt));
  hIncDist->SetDirectory(0);
  hProjX->SetDirectory(0);
  hProjY->SetDirectory(0);
  hProjX->GetYaxis()->SetTitle("Entries");
  hProjY->GetYaxis()->SetTitle("Entries");
  infile.Close();

  TH2F* hDummy = new TH2F("hDummy","",10,xmin,xmax,10,ymin,ymax);

  TLegend* l = new TLegend(0.68,0.65,0.89,0.89);
  l->AddEntry(lSet1," Set 1","l");
  l->AddEntry(lSet2," Set 2","l");
  l->AddEntry(lSet3," Set 3","l");
  l->SetTextSize(0.05);
  l->SetTextFont(42);

  hProjX->GetXaxis()->SetRangeUser(xmin,xmax);
  hProjY->GetXaxis()->SetRangeUser(ymin,ymax);
  Double_t nprompt = hProjY->GetMean();
  Double_t nprompterr = hProjY->GetRMS();
  Double_t nFD = hProjX->GetMean();
  Double_t nFDerr = hProjX->GetRMS();

  Double_t x1=0.32;
  if(iPt<=2)
    x1=0.35;

  cout << nprompt << "  " << nprompterr << endl;
  cout << nFD << "  " << nFDerr << endl;  
  
  TPaveText* info = new TPaveText(x1,0.2,0.4,0.35,"NDC");
  info->SetTextFont(132);
  info->SetFillStyle(0);
  info->SetFillColor(0);
  info->SetBorderSize(0);
  info->SetTextSize(0.05);
  info->Clear();
  info->AddText(Form("N_{prompt} = %0.f #pm %0.f",nprompt[iPt],nprompterr[iPt]));
  info->AddText(Form("N_{feed-down} = %0.f #pm %0.f",nFD[iPt],nFDerr[iPt]));

  TPaveText* info2 = new TPaveText(x1,0.2,0.4,0.35,"NDC");
  info2->SetTextFont(132);
  info2->SetFillStyle(0);
  info2->SetFillColor(0);
  info2->SetBorderSize(0);
  info2->SetTextSize(0.05);
  info2->Clear();
  info2->AddText(Form("N_{prompt} = %0.f",nprompt[iPt]));
  info2->AddText(Form("N_{feed-down} = %0.f",nFD[iPt]));
  
  TCanvas* cLines = new TCanvas("cLines","",1000,800);
  TCanvas* cLinesErr = new TCanvas("cLinesErr","",1000,800);
  TCanvas* cLinesDisp = new TCanvas("cLinesDisp","",1000,800);
  TCanvas* cProjX = new TCanvas("cProjX","",1000,800);
  TCanvas* cProjY = new TCanvas("cProjY","",1000,800);

  hIncDist->GetXaxis()->SetRangeUser(xmin,xmax);
  hIncDist->GetYaxis()->SetRangeUser(ymin,ymax);
  hIncDist->GetXaxis()->SetTitleSize(0.06);
  hIncDist->GetYaxis()->SetTitleSize(0.06);
  hIncDist->GetXaxis()->SetTitleOffset(1.1);
  hIncDist->GetYaxis()->SetTitleOffset(1.1);
  hIncDist->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",PtMin,PtMax));
  hIncDist->GetYaxis()->SetTitle("N_{prompt}");
  hIncDist->GetXaxis()->SetTitle("N_{feed-down}");
  hIncDist->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",PtMin,PtMax));
  hIncDist->SetTitleSize(0.035);
  hDummy->GetXaxis()->SetRangeUser(xmin,xmax);
  hDummy->GetYaxis()->SetRangeUser(ymin,ymax);
  hDummy->GetXaxis()->SetTitleOffset(1.);
  hDummy->GetYaxis()->SetTitleOffset(1.4);
  hDummy->GetXaxis()->SetTitle(hIncDist->GetXaxis()->GetTitle());
  hDummy->GetYaxis()->SetTitle(hIncDist->GetYaxis()->GetTitle());
  hDummy->GetXaxis()->SetTitleSize(0.06);
  hDummy->GetYaxis()->SetTitleSize(0.06);
  hDummy->GetXaxis()->SetLabelSize(0.05);
  hDummy->GetYaxis()->SetLabelSize(0.05);
  hDummy->SetTitle(hIncDist->GetTitle());
  hDummy->SetTitleSize(0.035);
  hDummy->SetStats(0);
  lSet1->SetLineWidth(2);
  lSet2->SetLineWidth(2);
  lSet3->SetLineWidth(2);
  lSet1->SetNpx(500);
  lSet2->SetNpx(500);
  lSet3->SetNpx(500);
  lSet1ErrLow->SetLineWidth(2);
  lSet2ErrLow->SetLineWidth(2);
  lSet3ErrLow->SetLineWidth(2);
  lSet1ErrLow->SetLineStyle(2);
  lSet2ErrLow->SetLineStyle(2);
  lSet3ErrLow->SetLineStyle(2);
  lSet1ErrLow->SetNpx(200);
  lSet2ErrLow->SetNpx(200);
  lSet3ErrLow->SetNpx(200);
  lSet1ErrHigh->SetLineWidth(2);
  lSet2ErrHigh->SetLineWidth(2);
  lSet3ErrHigh->SetLineWidth(2);
  lSet1ErrHigh->SetLineStyle(2);
  lSet2ErrHigh->SetLineStyle(2);
  lSet3ErrHigh->SetLineStyle(2);
  lSet1ErrHigh->SetNpx(200);
  lSet2ErrHigh->SetNpx(200);
  lSet3ErrHigh->SetNpx(200);
  cLines->cd();
  hDummy->GetXaxis()->SetTitleSize(0.06);
  hDummy->Draw();
  lSet1->Draw("same");
  lSet2->Draw("same");
  lSet3->Draw("same");
  l->Draw("same");
  info2->Draw("same");
  cLinesErr->cd();
  hDummy->GetXaxis()->SetTitleSize(0.06);
  hDummy->Draw();
  lSet1->Draw("same");
  lSet2->Draw("same");
  lSet3->Draw("same");  
  lSet1ErrLow->Draw("same");
  lSet2ErrLow->Draw("same");
  lSet3ErrLow->Draw("same");
  lSet1ErrHigh->Draw("same");
  lSet2ErrHigh->Draw("same");
  lSet3ErrHigh->Draw("same");
  l->Draw("same");
  info->Draw("same");
  cLinesDisp->cd();
  hIncDist->GetXaxis()->SetTitleSize(0.06);
  hIncDist->GetYaxis()->SetTitleSize(0.06);
  hIncDist->GetXaxis()->SetTitleOffset(1.);
  hIncDist->GetYaxis()->SetTitleOffset(1.4);
  hIncDist->GetXaxis()->SetLabelSize(0.05);
  hIncDist->GetYaxis()->SetLabelSize(0.05);
  hIncDist->GetYaxis()->SetRangeUser(ymin,ymax);
  hIncDist->Draw();  
  lSet1->Draw("same");
  lSet2->Draw("same");
  lSet3->Draw("same");
  lSet1ErrLow->Draw("same");
  lSet2ErrLow->Draw("same");
  lSet3ErrLow->Draw("same");
  lSet1ErrHigh->Draw("same");
  lSet2ErrHigh->Draw("same");
  lSet3ErrHigh->Draw("same");
  l->Draw("same");
  info->Draw("same");
  cProjX->cd();
  hProjX->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",PtMin,PtMax));
  hProjX->GetYaxis()->SetRangeUser(1,hProjX->GetMaximum()*1.2);
  hProjX->GetYaxis()->SetTitleOffset(1.4);
  hProjX->GetXaxis()->SetTitleSize(0.06);
  hProjX->GetYaxis()->SetTitleSize(0.06);
  hProjX->GetXaxis()->SetLabelSize(0.05);
  hProjX->GetYaxis()->SetLabelSize(0.05);
  hProjX->SetName("N_{feed-down} projection");
  hProjX->SetFillStyle(3004);
  hProjX->SetFillColor(kBlue);
  hProjX->Draw();
  cProjY->cd();
  hProjY->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",PtMin,PtMax));
  hProjY->GetYaxis()->SetRangeUser(1,hProjY->GetMaximum()*1.2);
  hProjY->GetXaxis()->SetTitleSize(0.06);
  hProjY->GetYaxis()->SetTitleOffset(1.4);
  hProjY->GetYaxis()->SetTitleSize(0.06);
  hProjY->GetXaxis()->SetLabelSize(0.05);
  hProjY->GetYaxis()->SetLabelSize(0.05);
  hProjY->SetName("N_{prompt} projection");
  hProjY->SetFillStyle(3004);
  hProjY->SetFillColor(kBlue);
  hProjY->Draw();
  infile.Close();

  cLines->SaveAs(Form("Lines_%0.f-%0.f.eps",PtMin,PtMax));
  cLinesErr->SaveAs(Form("LinesErr_%0.f-%0.f.eps",PtMin,PtMax));
  cLinesDisp->SaveAs(Form("LinesDisp_%0.f-%0.f.eps",PtMin,PtMax));
  cProjX->SaveAs(Form("Nfeeddown_disp_%0.f-%0.f.eps",PtMin,PtMax));
  cProjY->SaveAs(Form("Nprompt_disp_%0.f-%0.f.eps",PtMin,PtMax));

  if(type==0)
    return cLines;
  else if(type==1)
    return cLinesErr;
  else
    return cLinesDisp;
}
