const Int_t nPtBins = 7;
const Int_t nPtLims = nPtBins+1;
Double_t PtLims[nPtLims] = {2,3,4,5,6,8,12,16};

void PlotFitMassRange(Double_t nSigma=2, Double_t pTmin=4, Double_t pTmax=5, TString massfilename="Mass_4-5", TString IPfilename="FitUnbinned_4-5_bkg", TString masscanvasname="cMass");
void PlotDifferentRanges(Double_t pTmin=4, Double_t pTmax=5, TString massfilename="Mass_4-5", TString masscanvasname="cMass");
TCanvas* PlotFitMassWORange(Double_t nSigma=3,Double_t pTmin=4,Double_t pTmax=5,TString massfilename="Mass_4-5", TString masscanvasname="cMass");
void Plot3FitMassWORange();
void PlotAllMassFits();

void PlotAllMassFits() {
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    PlotFitMassRange(2,PtLims[iPt],PtLims[iPt+1],Form("Mass_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),Form("FitUnbinned_%0.f-%0.f_bkg",PtLims[iPt],PtLims[iPt+1]),"cMass");
  }
}

void Plot3FitMassWORange() {

  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetTitleOffset(1.6,"y");
  gStyle->SetTitleSize(0.08,"t");
  
  TCanvas* c3Mass = new TCanvas("c3Mass","",1600,800);
  c3Mass->Divide(3,1);

  Double_t PtMin[3]={2,4,12};
  Double_t PtMax[3]={3,5,16};
  
  for(Int_t iFit=0; iFit<3; iFit++) {
    TCanvas* c=(TCanvas*)PlotFitMassWORange(3,PtMin[iFit],PtMax[iFit],Form("Mass_%0.f-%0.f",PtMin[iFit],PtMax[iFit]),"cMass");

    c3Mass->cd(iFit+1)->SetLeftMargin(0.18);
    c3Mass->cd(iFit+1)->SetRightMargin(-0.25);

    c->DrawClonePad();
  }
  
  c3Mass->SaveAs("MassFit3Bin.eps");
  
}
void PlotFitMassRange(Double_t nSigma, Double_t pTmin, Double_t pTmax, TString massfilename, TString IPfilename, TString masscanvasname) {

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
  TPaveText* info3 = new TPaveText(0.15,0.55,0.5,0.7,"NDC");
  h->SetDirectory(0);
  inmassfile.Close();

  TFile inIPfile(Form("%s.root",IPfilename.Data()));
  TH1F* hImpPar = (TH1F*)inIPfile.Get("hImpPar");
  hImpPar->SetDirectory(0);
  hImpPar->Rebin(2);
  inIPfile.Close();
    
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

  Double_t y2 = h->GetMaximum()*1.03;
  Double_t y1 = h->GetMaximum()*0.1.;
  Double_t x1 = mean+nSigma*sigma;
  Double_t x2 = mean-nSigma*sigma;
  
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
  TBox* box = new TBox(x1,y1,x2,y2);
  box->SetFillColorAlpha(kGreen+3,0.2);
  
  info->Clear();
  info->SetTextSize(0.05);
  info->SetBorderSize(0);
  info->SetFillStyle(0);
  info->SetTextColor(kBlue);
  info->SetTextFont(132);
  info->AddText(Form("#mu = %.4f #pm %.4f",mean,errmean));
  info->AddText(Form("#sigma = %.4f #pm %.4f",sigma,errsigma));
  info2->Clear();    
  info2->SetTextSize(0.05);
  info2->SetBorderSize(0);
  info2->SetTextFont(132);
  info2->SetFillStyle(0);
  info2->AddText(Form("S (%0.f#sigma) = %.0f #pm %.0f",nSigma,signal,signalerr));
  info2->AddText(Form("B (%0.f#sigma) = %.0f #pm %.0f",nSigma,bkg,bkgerr));
  info3->SetTextSize(0.05);
  info3->SetBorderSize(0);
  info3->SetTextFont(132);
  info3->SetFillStyle(0);
  info3->SetTextColor(kGreen+3);
  info3->AddText("mass region");
  info3->AddText("|M-M_{peak}| < 2#sigma");
  
  gStyle->SetPadBottomMargin(0.12);
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
  info->Draw("same");
  info2->Draw("same");
  info3->Draw("same");
  box->Draw("same");
  
  cMass->SaveAs(Form("%s_fit.pdf",massfilename.Data()));

  TLatex lat;
  lat.SetTextFont(132);
  lat.SetTextSize(0.055);
  Double_t min=-1000;
  Double_t max=1000;
  
  TCanvas* cImpPar = new TCanvas("cImpPar","",1000,800);
  cImpPar->SetLogy();
  hImpPar->SetLineColor(kGreen+3);
  hImpPar->GetXaxis()->SetTitleSize(0.05);
  hImpPar->GetYaxis()->SetTitleSize(0.05);
  hImpPar->GetXaxis()->SetLabelSize(0.05);
  hImpPar->GetYaxis()->SetLabelSize(0.05);
  hImpPar->GetYaxis()->SetTitleOffset(1.2);
  hImpPar->GetXaxis()->SetTitleOffset(1.2);
  hImpPar->GetYaxis()->SetTitle(Form("Entries/(%0.f #mum)",hImpPar->GetBinWidth(10)));
  hImpPar->GetYaxis()->SetRangeUser(0.2,hImpPar->GetMaximum()*2.);
  hImpPar->GetXaxis()->SetNdivisions(508,kTRUE);
  hImpPar->Draw("E");
  lat.DrawLatex(min+TMath::Abs(min)*0.1, hImpPar->GetMaximum()*0.4,Form("%0.f < #it{p}_{T} < %0.f GeV/c", pTmin, pTmax));

  cImpPar->SaveAs(Form("ImpParData_%0.f-%0.f.eps",pTmin,pTmax));
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
  h->SetDirectory(0);
  inmassfile.Close();
    
  Double_t mean = fs->GetParameter(3);
  Double_t sigma = fs->GetParameter(4);
  Double_t errmean = fs->GetParError(3);
  Double_t errsigma = fs->GetParError(4);

  Double_t binwidth = h->GetBinWidth(10)*1000;//in MeV

  Double_t y2 = h->GetMaximum()*1.03;
  Double_t y1 = h->GetMaximum()*0.1.;

  Double_t nSigma[4] = {1,1.5,2.5,3};
  TLine **l1 = new TLine*[4];
  TLine **l2 = new TLine*[4];
  TLine **l3 = new TLine*[4];
  TLine **l4 = new TLine*[4];
  TBox **box = new TBox*[4];
  
  TPaveText** info3 = new TPaveText*[4];

  const Int_t colors[] = {kGreen+3,kMagenta,kOrange+7,kAzure+2};
  Double_t xpav1[4] = {0.2,0.2,0.65,0.65};
  Double_t xpav2[4] = {0.4,0.4,0.85,0.85};
  Double_t ypav1[4] = {0.7,0.55,0.7,0.55};
  Double_t ypav2[4] = {0.8,0.65,0.8,0.65};
  
  for(Int_t iSigma=0; iSigma<4; iSigma++) {
    Double_t x1 = mean+nSigma[iSigma]*sigma;
    Double_t x2 = mean-nSigma[iSigma]*sigma;
    Double_t y2 = h->GetMaximum()*(1.03-0.01*iSigma);
    Double_t y1 = h->GetMaximum()*(0.1.+0.01*iSigma);
    box[iSigma] = new TBox(x1,y1,x2,y2);
    box[iSigma]->SetFillColorAlpha(colors[iSigma],0.15);
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
    info3[iSigma] = new TPaveText(xpav1[iSigma],ypav1[iSigma],xpav2[iSigma],ypav2[iSigma],"NDC");
    info3[iSigma]->SetTextSize(0.04);
    info3[iSigma]->SetBorderSize(0);
    info3[iSigma]->SetTextFont(132);
    info3[iSigma]->SetFillStyle(0);
    info3[iSigma]->SetTextColor(colors[iSigma]);
    info3[iSigma]->AddText(Form("|M-M_{peak}| < %0.1f#sigma",nSigma[iSigma]));
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
  h->GetXaxis()->SetRangeUser(1.75,1.99);
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
  for(Int_t iSigma=0; iSigma<4; iSigma++) {
    l1[3-iSigma]->Draw("same");
    l2[3-iSigma]->Draw("same");
    l3[3-iSigma]->Draw("same");
    l4[3-iSigma]->Draw("same");
    info3[iSigma]->Draw("same");
    box[3-iSigma]->Draw("same");
  }

  cMass->SaveAs(Form("%s_ranges.pdf",massfilename.Data()));
}

TCanvas* PlotFitMassWORange(Double_t nSigma, Double_t pTmin, Double_t pTmax, TString massfilename, TString masscanvasname) {
  
  TFile inmassfile(Form("%s.root",massfilename.Data()),"UPDATE");
  TCanvas* c=(TCanvas*)inmassfile.Get(masscanvasname.Data());
  TH1F* h = (TH1F*)c->GetPrimitive("fhistoInvMass");
  TF1* fs = (TF1*)c->GetPrimitive("funcmass");
  TF1* fb = (TF1*)c->GetPrimitive("funcbkgFullRange");
  TPaveText* info = new TPaveText(0.35,0.79,0.75,0.89,"NDC");
  Double_t y1=0.15;
  Double_t y2=0.32;
  Double_t x1=0.33;
  Double_t x2=0.63;
  if(pTmin<=3) {
    x1=0.38;
    x2=0.68;
    y1=0.22;
    y2=0.39;
  }
  if(pTmin<=6 && pTmin>3) {
    x1=0.38;
    x2=0.68;
    y1=0.18;
    y2=0.35;
  }  
  TPaveText* info2 = new TPaveText(x1,y1,x2,y2,"NDC");
  h->SetDirectory(0);
  inmassfile.Close();
  
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

  Double_t y2 = h->GetMaximum()*1.03;
  Double_t y1 = 0.;
  Double_t x1 = mean+nSigma*sigma;
  Double_t x2 = mean-nSigma*sigma;
  
  info->Clear();
  info->SetTextSize(0.055);
  info->SetBorderSize(0);
  info->SetFillStyle(0);
  info->SetTextColor(kBlue);
  info->SetTextFont(132);
  if(pTmin>=12) {
    info->AddText(Form("#mu = %.3f #pm %.3f",mean,errmean));
    info->AddText(Form("#sigma = %.3f #pm %.3f",sigma,errsigma));
  } 
  else {
    info->AddText(Form("#mu = %.4f #pm %.4f",mean,errmean));
    info->AddText(Form("#sigma = %.4f #pm %.4f",sigma,errsigma));
  }
  info2->Clear();    
  info2->SetTextSize(0.055);
  info2->SetBorderSize(0);
  info2->SetTextFont(132);
  info2->SetFillStyle(0);
  info2->AddText(Form("Signif. (%0.f#sigma) = %.1f #pm %.1f",nSigma,significance,significanceerr));
  info2->AddText(Form("S (%0.f#sigma) = %.0f #pm %.0f",nSigma,signal,signalerr));
  info2->AddText(Form("B (%0.f#sigma) = %.0f #pm %.0f",nSigma,bkg,bkgerr));
  info2->AddText(Form("S/B (%0.f#sigma) = %.4f",nSigma,signaloverbkg));
  
  TCanvas* cMass = new TCanvas(masscanvasname.Data(),"",1300,900);
  h->SetLineColor(kBlack);
  h->Draw("E");
  h->SetStats(0);
  h->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c", pTmin, pTmax));
  h->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitle(Form("Entries/(%0.f Mev/c^{2})",binwidth));
  h->GetYaxis()->SetRangeUser(h->GetMaximum()*0.,h->GetMaximum()*1.25);
  fs->Draw("same");
  info->Draw("same");
  info2->Draw("same");

  cMass->SaveAs(Form("%s_fit_norange.eps",massfilename.Data()));

  return cMass;
}
