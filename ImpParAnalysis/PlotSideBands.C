const Int_t nPtBins = 7;
const Int_t nPtLims = nPtBins+1;
Double_t PtLims[nPtLims] = {2,3,4,5,6,8,12,16};

void PlotAllSideBands();
void PlotSideBands(Int_t nSigmaMin=4, Int_t nSigmaMax=15, Double_t pTmin = 4, Double_t pTmax=5, TString massfilename="Mass_4-5", TString masscanvasname = "cMass", TString impparbkgfilename="Sidebands_Pt_4-5");
void PlotDifferentRanges(Double_t pTmin = 4, Double_t pTmax=5, TString massfilename="Mass_4-5", TString masscanvasname = "cMass");

void PlotAllSideBands() {

  Double_t minSB = 4;
  Double_t maxSB = 15;
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    if(PtLims[iPt]>=4 && PtLims[iPt]<8) {
      minSB = 4;
      maxSB = 17;
    }
    else if(PtLims[iPt]>=8 && PtLims[iPt]<=16) {
      minSB = 4;
      maxSB = 12;
    }    
    else {
      minSB = 4;
      maxSB = 20;
    }
    
    PlotSideBands(minSB,maxSB,PtLims[iPt],PtLims[iPt+1],Form("Mass_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"cMass",Form("ImpParBkg_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]));
  }
  
}

void PlotSideBands(Int_t nSigmaMin, Int_t nSigmaMax, Double_t pTmin, Double_t pTmax, TString massfilename, TString masscanvasname, TString impparbkgfilename) {

  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetLegendBorderSize(0);
  
  TFile inmassfile(Form("%s.root",massfilename.Data()),"UPDATE");
  TCanvas* c=(TCanvas*)inmassfile.Get(masscanvasname.Data());
  TH1F* h = (TH1F*)c->GetPrimitive("fhistoInvMass");
  TF1* fs = (TF1*)c->GetPrimitive("funcmass");
  TF1* fb = (TF1*)c->GetPrimitive("funcbkgFullRange");
  TPaveText* info = new TPaveText(0.2,0.7,0.45,0.85,"NDC");
  TPaveText* info2 = new TPaveText(0.1,0.6,0.5,0.9,"NDC");
  TPaveText* info3 = new TPaveText(0.58,0.8,0.88,0.85,"NDC");
  TPaveText* info4 = new TPaveText(0.58,0.7,0.88,0.8,"NDC");
  TPaveText* info5 = new TPaveText(0.59,0.65,0.89,0.7,"NDC");
  h->SetDirectory(0);
  inmassfile.Close();

  TFile inimpparbkgfile(Form("%s.root",impparbkgfilename.Data()),"UPDATE");
  TH1F* hImpParTot = (TH1F*)inimpparbkgfile.Get("hSumDist");
  TH1F* hImpParLeft = (TH1F*)inimpparbkgfile.Get("hLeftDist");
  TH1F* hImpParRight = (TH1F*)inimpparbkgfile.Get("hRightDist");
  hImpParTot->SetDirectory(0);
  hImpParLeft->SetDirectory(0);
  hImpParRight->SetDirectory(0);
  hImpParTot->Rebin(2);
  hImpParLeft->Rebin(2);
  hImpParRight->Rebin(2);
  inimpparbkgfile.Close();
  
  Int_t nSigma=3;
  Double_t mean = fs->GetParameter(3);
  Double_t sigma = fs->GetParameter(4);
  Double_t errmean = fs->GetParError(3);
  Double_t errsigma = fs->GetParError(4);
  Double_t ints=fs->Integral(mean-nSigma*sigma,mean+nSigma*sigma)/h->GetBinWidth(4);
  Double_t intb=fb->Integral(mean-nSigma*sigma,mean+nSigma*sigma)/h->GetBinWidth(2);
  Double_t signal = ints-intb;
  Double_t signalerr = fs->GetParError(fs->GetNpar()-3)/fs->GetParameter(fs->GetNpar()-3)*signal;
  Double_t bkg = intb;
  Double_t bkgerr = fb->GetParError(0)/fb->GetParameter(0)*bkg;
  Double_t significance = signal/TMath::Sqrt(signal+bkg);
  Double_t significanceerr = significance*TMath::Sqrt((signalerr*signalerr+bkgerr*bkgerr)/(4.*(signal+bkg)*(signal+bkg))+(bkg/(signal+bkg))*(signalerr*signalerr)/signal/signal);
  Double_t signaloverbkg = signal/bkg;

  Double_t binwidth = h->GetBinWidth(10)*1000;//in MeV

  Double_t y2 = h->GetBinContent(5)*1.1;
  Double_t y1 = h->GetMaximum()*0.1;
  Double_t x1 = mean+nSigmaMin*sigma;
  Double_t x2 = mean+nSigmaMax*sigma;
  
  TLine *l1 = new TLine(x1,y1,x1,y2);
  l1->SetLineColor(kGreen+3);
  l1->SetLineWidth(2);
  TLine *l2 = new TLine(x2,y1,x2,y2);
  l2->SetLineColor(kGreen+3);
  l2->SetLineWidth(2);
  TLine *l3 = new TLine(x1,y1,x2,y1);
  l3->SetLineColor(kGreen+3);
  l3->SetLineWidth(2);
  TLine *l4 = new TLine(x1,y2,x2,y2);
  l4->SetLineColor(kGreen+3);
  l4->SetLineWidth(2);
  TBox *box1 = new TBox(x1,y1,x2,y2);
  box1->SetFillColorAlpha(kGreen+3,0.2);
  
  Double_t x3 = mean-nSigmaMin*sigma;
  Double_t x4 = mean-nSigmaMax*sigma;
  
  TLine *l5 = new TLine(x3,y1,x3,y2);
  l5->SetLineColor(kOrange+7);
  l5->SetLineWidth(2);
  TLine *l6 = new TLine(x4,y1,x4,y2);
  l6->SetLineColor(kOrange+7);
  l6->SetLineWidth(2);
  TLine *l7 = new TLine(x3,y1,x4,y1);
  l7->SetLineColor(kOrange+7);
  l7->SetLineWidth(2);
  TLine *l8 = new TLine(x3,y2,x4,y2);
  l8->SetLineColor(kOrange+7);
  l8->SetLineWidth(2);
  TBox *box2 = new TBox(x3,y1,x4,y2);
  box2->SetFillColorAlpha(kOrange+7,0.2);
  
  info->Clear();
  info->SetTextSize(0.05);
  info->SetBorderSize(0);
  info->SetFillStyle(0);
  info->SetTextColor(kBlue);
  info->SetTextFont(132);
  info->AddText(Form("#mu = %.4f #pm %.4f",mean,errmean));
  info->AddText(Form("#sigma = %.4f #pm %.4f",sigma,errsigma));
  info2->Clear();    
  info2->SetTextSize(0.04);
  info2->SetBorderSize(0);
  info2->SetTextFont(132);
  info2->SetFillStyle(0);
  info2->AddText(Form("Significance (%d#sigma) = %.1f #pm %.1f",nSigma,significance,significanceerr));
  info2->AddText(Form("S (%d#sigma) = %.0f #pm %.0f",nSigma,signal,signalerr));
  info2->AddText(Form("B (%d#sigma) = %.0f #pm %.0f",nSigma,bkg,bkgerr));
  info2->AddText(Form("S/B (%d#sigma) = %.4f",nSigma,signaloverbkg));
  info3->SetTextSize(0.05);
  info3->SetBorderSize(0);
  info3->SetTextFont(132);
  info3->SetFillStyle(0);
  info3->SetTextColor(kBlack);
  info3->AddText("side-band regions");
  info4->SetTextSize(0.05);
  info4->SetBorderSize(0);
  info4->SetTextFont(132);
  info4->SetFillStyle(0);
  info4->SetTextColor(kOrange+7);
  info4->AddText("-15#sigma < M-M_{peak} < -4#sigma");
  info5->SetTextSize(0.05);
  info5->SetBorderSize(0);
  info5->SetTextFont(132);
  info5->SetFillStyle(0);
  info5->SetTextColor(kGreen+3);
  info5->AddText("4#sigma < M-M_{peak} < 15#sigma");
  
  TCanvas* cMass = new TCanvas(masscanvasname.Data(),"",1000,800);
  h->SetLineColor(kBlack);
  h->Draw("E");
  h->SetStats(0);
  h->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c", pTmin, pTmax));
  h->SetTitleSize(0.04);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitle(Form("Entries/(%0.f Mev/c^{2})",binwidth));
  h->GetYaxis()->SetRangeUser(h->GetMaximum()*0.1,h->GetMaximum()*1.1);
  fs->Draw("same");
  l1->Draw("same");
  l2->Draw("same");
  l3->Draw("same");
  l4->Draw("same");
  l5->Draw("same");
  l6->Draw("same");
  l7->Draw("same");
  l8->Draw("same");
  info->Draw("same");
  info3->Draw("same");
  info4->Draw("same");
  info5->Draw("same");
  box1->Draw("same");
  box2->Draw("same");

  cMass->SaveAs(Form("%s.pdf",massfilename.Data()));

  TLatex lat;
  lat.SetTextFont(132);
  lat.SetTextSize(0.055);
  Double_t min=-1000;
  Double_t max=1000;

  TLegend* l = new TLegend(0.62,0.69,0.85,0.89);
  l->SetTextSize(0.045);
  l->AddEntry(hImpParTot,"Weighted sum","le");
  l->AddEntry(hImpParLeft,"Left region","le");
  l->AddEntry(hImpParRight,"Right region","le");
 
  TCanvas* cBkg = new TCanvas("cBkg","",1000,800);
  cBkg->SetLogy();
  hImpParTot->SetStats(0);
  hImpParTot->SetLineColor(kBlack);
  hImpParTot->SetTitle("");
  hImpParTot->GetXaxis()->SetTitle("Imp Par XY (#mum)");
  hImpParTot->GetYaxis()->SetTitle(Form("Entries/(%0.f #mum)",hImpParTot->GetBinWidth(10)));
  hImpParTot->GetXaxis()->SetTitleSize(0.05);
  hImpParTot->GetYaxis()->SetTitleSize(0.05);
  hImpParTot->GetXaxis()->SetLabelSize(0.05);
  hImpParTot->GetYaxis()->SetLabelSize(0.05);
  hImpParTot->GetYaxis()->SetTitleOffset(1.2);
  hImpParTot->GetXaxis()->SetTitleOffset(1.2);
  hImpParTot->SetNdivisions(508);
  hImpParTot->SetStats(0);
//  hImpParTot->GetFunction("ImpParBkgFunc")->SetBit(TF1::kNotDraw);
  hImpParTot->Draw("E");
  hImpParLeft->SetLineColor(kOrange+7);
  hImpParLeft->Draw("Esame");
  hImpParRight->SetLineColor(kGreen+3);
  hImpParRight->Draw("Esame");
  lat.DrawLatex(min+TMath::Abs(min)*0.1, hImpParTot->GetMaximum()*0.7,Form("%0.f < #it{p}_{T} < %0.f GeV/c", pTmin, pTmax));
  l->Draw("same");
  
  cBkg->SaveAs(Form("%s_NoFit.eps",impparbkgfilename.Data()));
}

void PlotDifferentRanges(Double_t pTmin, Double_t pTmax, TString massfilename, TString masscanvasname) {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);
  
  TFile inmassfile(Form("%s.root",massfilename.Data()),"UPDATE");
  TCanvas* c=(TCanvas*)inmassfile.Get(masscanvasname.Data());
  TH1F* h = (TH1F*)c->GetPrimitive("fhistoInvMass");
  TF1* fs = (TF1*)c->GetPrimitive("funcmass");
  TF1* fb = (TF1*)c->GetPrimitive("funcbkgFullRange");
  TPaveText* info = new TPaveText(0.6,0.7,0.85,0.85,"NDC");
  TPaveText* info2 = new TPaveText(0.6,0.45,0.85,0.6,"NDC");
  TPaveText** info3 = new TPaveText*[5];
  h->SetDirectory(0);
  inmassfile.Close();
    
  Double_t mean = fs->GetParameter(3);
  Double_t sigma = fs->GetParameter(4);
  Double_t errmean = fs->GetParError(3);
  Double_t errsigma = fs->GetParError(4);

  Double_t binwidth = h->GetBinWidth(10)*1000;//in MeV

  Double_t nSigmaMin[5] = {3,4,4,5,6};
  Double_t nSigmaMax[5] = {14,15,15,16,17};
  TLine **l1 = new TLine*[5];
  TLine **l2 = new TLine*[5];
  TLine **l3 = new TLine*[5];
  TLine **l4 = new TLine*[5];
  TBox **box = new TBox*[5];
  TLine **l5 = new TLine*[5];
  TLine **l6 = new TLine*[5];
  TLine **l7 = new TLine*[5];
  TLine **l8 = new TLine*[5];
  TBox **box2 = new TBox*[5];
  
  TPaveText** info3 = new TPaveText*[5];

  const Int_t colors[] = {kMagenta,kOrange+7,kGreen+3,kOrange+3,kAzure+2};
  Double_t xpav1[5] = {0.58,0.18,0.18,0.58,0.58};
  Double_t xpav2[5] = {0.88,0.48,0.48,0.88,0.88};
  Double_t ypav1[5] = {0.8,0.75,0.68,0.73,0.66};
  Double_t ypav2[5] = {0.85,0.8,0.73,0.78,0.71};

  for(Int_t iSigma=0; iSigma<5; iSigma++) {
    Double_t step=0.03*iSigma;
    if(iSigma>=2)
      step=0.03*(iSigma-1);
    Double_t x1 = mean-nSigmaMax[iSigma]*sigma;
    Double_t x2 = mean-nSigmaMin[iSigma]*sigma;    
    Double_t y2 = h->GetMaximum()*(0.5+step);
    Double_t y1 = h->GetMaximum()*(0.2.+step);
    Double_t x3 = mean+nSigmaMin[iSigma]*sigma;
    Double_t x4 = mean+nSigmaMax[iSigma]*sigma;
    Double_t y3 = h->GetMaximum()*(0.52-step);
    Double_t y4 = h->GetMaximum()*(0.22-step);
    box[iSigma] = new TBox(x1,y1,x2,y2);
    box[iSigma]->SetFillColorAlpha(colors[iSigma],0.1);
    box2[iSigma] = new TBox(x3,y3,x4,y4);
    box2[iSigma]->SetFillColorAlpha(colors[iSigma],0.1);
    l1[iSigma] = new TLine(x1,y1,x1,y2);
    l1[iSigma]->SetLineColor(colors[iSigma]);
    l1[iSigma]->SetLineWidth(2);
    l2[iSigma] = new TLine(x2,y1,x2,y2);
    l2[iSigma]->SetLineColor(colors[iSigma]);
    l2[iSigma]->SetLineWidth(2);
    l3[iSigma] = new TLine(x1,y1,x2,y1);
    l3[iSigma]->SetLineColor(colors[iSigma]);
    l3[iSigma]->SetLineWidth(2);
    l4[iSigma] = new TLine(x1,y2,x2,y2);
    l4[iSigma]->SetLineColor(colors[iSigma]);
    l4[iSigma]->SetLineWidth(2);
    l5[iSigma] = new TLine(x3,y3,x3,y4);
    l5[iSigma]->SetLineColor(colors[iSigma]);
    l5[iSigma]->SetLineWidth(2);
    l6[iSigma] = new TLine(x4,y3,x4,y4);
    l6[iSigma]->SetLineColor(colors[iSigma]);
    l6[iSigma]->SetLineWidth(2);
    l7[iSigma] = new TLine(x3,y3,x4,y3);
    l7[iSigma]->SetLineColor(colors[iSigma]);
    l7[iSigma]->SetLineWidth(2);
    l8[iSigma] = new TLine(x3,y4,x4,y4);
    l8[iSigma]->SetLineColor(colors[iSigma]);
    l8[iSigma]->SetLineWidth(2);
    info3[iSigma] = new TPaveText(xpav1[iSigma],ypav1[iSigma],xpav2[iSigma],ypav2[iSigma],"NDC");
    info3[iSigma]->SetTextSize(0.04);
    info3[iSigma]->SetBorderSize(0);
    info3[iSigma]->SetTextFont(132);
    info3[iSigma]->SetFillStyle(0);
    info3[iSigma]->SetTextColor(colors[iSigma]);
    if(iSigma==1) {
      info3[iSigma]->AddText(Form("-%0.f#sigma < M-M_{peak} < -%0.f#sigma",nSigmaMin[iSigma],nSigmaMax[iSigma]));
    }
    else if(iSigma==2) {
      info3[iSigma]->AddText(Form("%0.f#sigma < M-M_{peak} < %0.f#sigma",nSigmaMin[iSigma],nSigmaMax[iSigma]));
    }
    else{
      info3[iSigma]->AddText(Form("%0.f#sigma < |M-M_{peak}| < %0.f#sigma",nSigmaMin[iSigma],nSigmaMax[iSigma]));
    }
  }
  
  info->Clear();
  info->SetTextSize(0.05);
  info->SetBorderSize(0);
  info->SetFillStyle(0);
  info->SetTextColor(kBlue);
  info->SetTextFont(132);
  info->AddText(Form("#mu = %.4f #pm %.4f",mean,errmean));
  info->AddText(Form("#sigma = %.4f #pm %.4f",sigma,errsigma));
  gStyle->SetPadBottomMargin(0.12);
  TCanvas* cMass = new TCanvas(masscanvasname.Data(),"",800,800);
  h->SetLineColor(kBlack);
  h->Draw("E");
  h->SetStats(0);
  h->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c", pTmin, pTmax));
  h->SetTitleSize(0.04);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.45);
  h->GetXaxis()->SetTitleOffset(1.);
  h->GetYaxis()->SetTitle(Form("Entries/(%0.f Mev/c^{2})",binwidth));
  h->GetYaxis()->SetRangeUser(h->GetMaximum()*0.1,h->GetMaximum()*1.1);
  fs->Draw("same");
  for(Int_t iSigma=0; iSigma<5; iSigma++) {
    if(iSigma!=2) {
      l1[4-iSigma]->Draw("same");
      l2[4-iSigma]->Draw("same");
      l3[4-iSigma]->Draw("same");
      l4[4-iSigma]->Draw("same");
      box[4-iSigma]->Draw("same");
    }
    if(iSigma!=1) {
      l5[iSigma]->Draw("same");
      l6[iSigma]->Draw("same");
      l7[iSigma]->Draw("same");
      l8[iSigma]->Draw("same");
      box2[iSigma]->Draw("same");
    }
    info3[iSigma]->Draw("same");
  }

  cMass->SaveAs(Form("%s_SBranges.pdf",massfilename.Data()));
}
