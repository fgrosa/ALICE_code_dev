const Int_t nPtBins = 10;
const Int_t nPtLims = nPtBins+1;
Double_t PtLims[nPtLims] = {1,2,3,4,5,6,7,8,12,16,24};

void PlotAllMassFits(TString cutset="cutset1", Double_t height=1000, Double_t width=800, Bool_t isOneFig=kTRUE, Int_t nrows=2, Int_t ncols=5);
TCanvas* PlotMassFit(Double_t nSigma=3, Double_t pTmin = 1, Double_t pTmax=2, Double_t ymin=300, Double_t ymax=1000, TString massfilename="rawyields_Dplus_cutset1", TString masscanvasname = "cRaw_Pt0",Double_t x1Info=0.4,Double_t x2Info=0.6,Double_t y1Info=0.15,Double_t y2Info=0.3, Double_t infoheight=0.1);
void Plot3MassFits(TString cutset="cutset1",Int_t iPt1=0,Int_t iPt2=4, Int_t iPt3=9);
void CompareWOchicut(Double_t nSigma=3, TString cutset="cutset1", Double_t pTmin=5, Double_t pTmax=6, TString masscanvasname="cMass_Pt4", Double_t rebin=5);


void Plot3MassFits(TString cutset,Int_t iPt1,Int_t iPt2, Int_t iPt3) {

  Double_t ymin[nPtBins] = {0,600,400,0,0,0,0,0,0,0};
  Double_t ymax[nPtBins] = {280,2200,2000,1400,1300,360,300,375,120,50};
  Double_t x1[nPtBins] = {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25};
  Double_t x2[nPtBins] = {0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7};
  Double_t y1[nPtBins] = {0.78,0.78,0.78,0.78,0.78,0.78,0.78,0.78,0.78,0.78};
  Double_t y2[nPtBins] = {0.88,0.88,0.88,0.88,0.88,0.88,0.88,0.88,0.88,0.88};

  TCanvas *c3Mass = new TCanvas("c3Mass","",1600,800);
  c3Mass->Divide(3,1);
  Int_t ptbins[3] = {iPt1,iPt2,iPt3};
  
  for(Int_t iPt=0; iPt<3; iPt++) {
    TCanvas* c1=(TCanvas*)PlotMassFit(3,PtLims[ptbins[iPt]],PtLims[ptbins[iPt]+1],ymin[ptbins[iPt]],ymax[ptbins[iPt]],Form("rawyields_Dplus_%s",cutset.Data()),Form("cRaw_Pt%d",ptbins[iPt]),x1[ptbins[iPt]],x2[ptbins[iPt]],y1[ptbins[iPt]],y2[ptbins[iPt]]);
    c3Mass->cd(iPt+1);
    c1->DrawClonePad();
  }

  c3Mass->SaveAs(Form("MassFits_KF_Pt_%d_%d_%d.eps",iPt1,iPt2,iPt3));
}

void PlotAllMassFits(TString cutset, Double_t height, Double_t width, Bool_t isOneFig, Int_t nrows, Int_t ncols) {
  TCanvas* cAll = 0x0;
  TCanvas** cMass = new TCanvas*[nPtBins];

  if(isOneFig) {
    cAll = new TCanvas("cAll","",width,height);
    cAll->Divide(nrows,ncols);
  }
  else {
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      cMass[iPt] = new TCanvas(Form("cMass_%d",iPt),"",width,height);
    }
  }
    
  Double_t ymin[nPtBins] = {0,600,400,0,0,0,0,0,0,0};
  Double_t ymax[nPtBins] = {160,2200,2000,1400,700,360,300,375,120,30};
  Double_t x1[nPtBins] = {0.3,0.3,0.3,0.3,0.3,0.3,0.25,0.3,0.25,0.25};
  Double_t x2[nPtBins] = {0.6,0.6,0.6,0.6,0.6,0.6,0.45,0.6,0.45,0.45};
  Double_t y1[nPtBins] = {0.15,0.15,0.15,0.15,0.15,0.15,0.7,0.14,0.5,0.5};
  Double_t y2[nPtBins] = {0.3,0.3,0.3,0.3,0.3,0.3,0.85,0.29,0.65,0.65};
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    TCanvas* c1=(TCanvas*)PlotMassFit(3,PtLims[iPt],PtLims[iPt+1],ymin[iPt],ymax[iPt],Form("rawyields_Dplus_%s",cutset.Data()),Form("cRaw_Pt%d",iPt),x1[iPt],x2[iPt],y1[iPt],y2[iPt],0.1);
    if(isOneFig)
      cAll->cd(iPt+1);
    else
      cMass[iPt]->cd();
    c1->DrawClonePad();

    if(!isOneFig) {
      cMass[iPt]->SaveAs(Form("MassFits_%s_%d.eps",cutset.Data(),iPt));
      delete cMass[iPt];
    }
  }
  
  if(isOneFig) {
    cAll->SaveAs(Form("MassFits_%s.eps",cutset.Data()));
    delete cAll;
  }
  else
    delete[] cMass;
}

TCanvas* PlotMassFit(Double_t nSigma, Double_t pTmin, Double_t pTmax, Double_t ymin, Double_t ymax, TString massfilename, TString masscanvasname,Double_t x1Info,Double_t x2Info,Double_t y1Info,Double_t y2Info, Double_t infoheight) {

  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.19);
  gStyle->SetPadRightMargin(-0.15);
  gStyle->SetTitleSize(0.06,"xy");
  gStyle->SetTitleSize(0.07,"t");
  
  TFile inmassfile(Form("%s.root",massfilename.Data()),"UPDATE");
  TCanvas* c=(TCanvas*)inmassfile.Get(masscanvasname.Data());
  TList* l=(TList*)c->GetListOfPrimitives();
  
  TH1F* h = (TH1F*)c->GetPrimitive("fhistoInvMass");
  TF1* fs = (TF1*)c->GetPrimitive("funcmass");
  TF1* fb = (TF1*)c->GetPrimitive("funcbkgFullRange");
  h->SetDirectory(0);
  inmassfile.Close();

  Double_t y1=0.68;

  if(pTmin==7)
    y1=0.5;
  
  Double_t y2=y1+infoheight;

  TPaveText *info1 = new TPaveText(0.25,y1,0.7,y2,"NDC");
  TPaveText *info2 = new TPaveText(x1Info,y1Info,x2Info,y2Info,"NDC"); 
  TPaveText *info3 = new TPaveText(0.45,0.55,0.8,0.55+infoheight,"NDC"); 
    
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
  
  info3->Clear();
  info3->SetTextSize(0.06);
  info3->SetBorderSize(0);
  info3->SetFillStyle(0);
  info3->SetTextColor(kBlue);
  info3->SetTextFont(132);
  if(pTmin<2 || pTmin>12) {
    info3->AddText(Form("#mu = %.3f #pm %.3f",mean,errmean));
    info3->AddText(Form("#sigma = %.3f #pm %.3f",sigma,errsigma));
  }
  else {
    info3->AddText(Form("#mu = %.4f #pm %.4f",mean,errmean));
    info3->AddText(Form("#sigma = %.4f #pm %.4f",sigma,errsigma));
  }
  
  info1->Clear();
  info1->SetTextSize(0.06);
  info1->SetBorderSize(0);
  info1->SetTextFont(132);
  info1->SetFillStyle(0);
  info1->AddText(Form("S (%0.f#sigma) = %.0f #pm %.0f",nSigma,signal,signalerr));
  info1->AddText(Form("B (%0.f#sigma) = %.0f #pm %.0f",nSigma,bkg,bkgerr));
  
  info2->Clear();
  info2->SetTextSize(0.06);
  info2->SetBorderSize(0);
  info2->SetTextFont(132);
  info2->SetFillStyle(0);
  info2->AddText(Form("Signif. (%0.f#sigma) = %.1f #pm %.1f",nSigma,significance,significanceerr));
  info2->AddText(Form("S/B (%0.f#sigma) = %.4f",nSigma,signaloverbkg));

  Double_t binwidth = h->GetBinWidth(10)*1000;//in MeV

  TCanvas* cMass = new TCanvas("cMass","",1200,900);
  
  cMass->Clear();
  h->SetLineColor(kBlack);
  h->SetStats(0);
  h->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",pTmin,pTmax));
  h->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(1.);
  h->GetYaxis()->SetTitleOffset(1.6);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitle(Form("Entries/(%0.f Mev/c^{2})",binwidth));
  h->GetYaxis()->SetRangeUser(ymin,ymax);
  h->Draw("E");  
  fs->Draw("same");
  info1->Draw("same");
  info2->Draw("same");
  info3->Draw("same");
  
  return cMass;

}

void CompareWOchicut(Double_t nSigma, TString cutset, Double_t pTmin, Double_t pTmax, TString masscanvasname, Double_t rebin) {

  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleSize(0.06,"xyt");
  gStyle->SetLabelSize(0.05,"xy");

  TString masscanvasnameclone=masscanvasname;
  
  TFile inmassfileMC(Form("rawyieldsMC_Dplus_%s.root",cutset.Data()),"UPDATE");
  TH1F* hSigmaMC = (TH1F*)inmassfileMC.Get("hSigma");
  hSigmaMC->SetDirectory(0);
  inmassfileMC.Close();
  Int_t iPt = hSigmaMC->GetXaxis()->FindBin(pTmin*1.0001);
  Double_t sigmaMC = hSigmaMC->GetBinContent(iPt+1);
  
  TFile inmassfile(Form("rawyields_Dplus_%s.root",cutset.Data()),"UPDATE");
  TCanvas* c=(TCanvas*)inmassfile.Get(masscanvasname.Data());
  TList* l=(TList*)c->GetListOfPrimitives();

  TString histoname = masscanvasname.ReplaceAll("cMass","hMass");
  TString histonamechi = masscanvasnameclone.ReplaceAll("cMass","hMassChi2Cut");
  
  TH1F* hMass = (TH1F*)l->FindObject(histoname.Data());
  TH1F* hMassChi = (TH1F*)l->FindObject(histonamechi.Data());
  hMass->Rebin(rebin);
  hMassChi->Rebin(rebin);
  
  hMass->GetYaxis()->SetTitle(Form("Entries/(%.0f MeV)",hMass->GetBinWidth(20)*1000));
  hMass->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/c",pTmin,pTmax));
  hMass->GetYaxis()->SetRangeUser(hMassChi->GetMinimum()*0., hMass->GetBinContent(100/rebin)*1.2);
  hMass->GetYaxis()->SetTitleOffset(1.4);
  hMass->GetXaxis()->SetTitleOffset(1.);
  hMass->GetYaxis()->SetTitleSize(0.05);
  hMass->GetXaxis()->SetTitleSize(0.05);
  hMass->GetYaxis()->SetLabelSize(0.05);
  hMass->GetXaxis()->SetLabelSize(0.05);
  hMass->SetLineColor(kRed);
  hMassChi->SetLineColor(kBlue);

  hMass->SetMarkerStyle(20);
  hMassChi->SetMarkerStyle(20);

  AliHFMassFitter* fitter = new AliHFMassFitter(hMass,1.67,2.05,1,0,0);
  fitter->SetFixGaussianSigma(sigmaMC);
  fitter->SetInitialGaussianMean(1.869);
  fitter->SetUseLikelihoodFit();
  fitter->MassFitter(kFALSE);
  AliHFMassFitter* fitterchi = new AliHFMassFitter(hMassChi,1.67,2.05,1,0,0);
  fitterchi->SetFixGaussianSigma(sigmaMC);
  fitterchi->SetInitialGaussianMean(1.869);
  fitterchi->MassFitter(kFALSE);
  fitterchi->SetUseLikelihoodFit();
  TF1* func = (TF1*)fitter->GetMassFunc();
  func->SetLineColor(kRed);
  TF1* funcchi = (TF1*)fitterchi->GetMassFunc();

  TPaveText* info = new TPaveText(0.24,0.55,0.44,0.8,"NDC");
  info->SetFillStyle(0);
  info->SetBorderSize(0);
  info->SetTextFont(132);
  info->SetTextSize(0.045);
  info->SetTextColor(kRed);
  info->AddText("W/o #chi^{2} cut");
  Double_t signal, signalerr;
  Double_t bkg, bkgerr;
  fitter->Signal(3,signal,signalerr);
  fitter->Background(3,bkg,bkgerr);
  info->AddText(Form("S(3#sigma) = %0.f #pm %0.f",signal,signalerr));
  info->AddText(Form("B(3#sigma) = %0.f #pm %0.f",bkg,bkgerr));
  TPaveText* infochi = new TPaveText(0.58,0.55,0.89,0.8,"NDC");
  infochi->SetFillStyle(0);
  infochi->SetBorderSize(0);
  infochi->SetTextFont(132);
  infochi->SetTextSize(0.045);
  infochi->SetTextColor(kBlue);
  Double_t signalchi, signalerrchi;
  Double_t bkgchi, bkgerrchi;
  fitterchi->Signal(3,signalchi,signalerrchi);
  fitterchi->Background(3,bkgchi,bkgerrchi);
  infochi->AddText("With #chi^{2} cut");
  infochi->AddText(Form("S(3#sigma) = %0.f #pm %0.f",signalchi,signalerrchi));
  infochi->AddText(Form("B(3#sigma) = %0.f #pm %0.f",bkgchi,bkgerrchi));
  
  TCanvas* cComp = new TCanvas("cComp","",800,800);
  hMass->Draw("E");
  hMassChi->Draw("Esame");
  func->Draw("same");
  funcchi->Draw("same");
  info->Draw("same");
  infochi->Draw("same");

  cComp->SaveAs(Form("KFMassComp_Pt%d.eps",iPt));
}
