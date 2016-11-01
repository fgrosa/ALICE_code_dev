const Int_t nPtBins = 7;
const Int_t nPtLims = nPtBins+1;
Double_t PtLims[nPtLims] = {2,3,4,5,6,8,12,16};

void PlotAllMassFits(Bool_t isOneFig=kTRUE);
TCanvas* PlotFitMass3Set(Double_t nSigma=3,
                         Double_t pTmin = 4, Double_t pTmax=5,
                         TString massfilename="Fitter_2",
                         TString masscanvasname = "cFitter",
                         Int_t num=1,
                         Int_t titletype=0, //0->set 1->pt bin
                         Double_t min1=500,Double_t min2=0, Double_t min3=0,
                         Double_t max1=1800,Double_t max2=1000, Double_t max3=100,
                         Double_t x1Info1=0.4, Double_t x2Info1=0.6, Double_t y1Info1=0.78, Double_t y2Info1=0.88,
                         Double_t x1Info2=0.4, Double_t x2Info2=0.6, Double_t y1Info2=0.15, Double_t y2Info2=0.3,
                         Bool_t DrawInfo3=kFALSE,
                         Double_t x1Info3=0, Double_t x2Info3=0, Double_t y1Info3=0, Double_t y2Info3=0,
                         Bool_t Draw3Sets=kTRUE,
                         Double_t xLatex1=1.77, Double_t xLatex2=1.77, Double_t xLatex3=1.77,
                         Double_t yLatex1=800, Double_t yLatex2=790, Double_t yLatex3=75);

void PlotAllMassFits(Bool_t isOneFig) {
  TCanvas* cAllSet1=0x0;
  TCanvas* cAllSet2=0x0;
  TCanvas* cAllSet3=0x0;
  TCanvas** cSet1 = new TCanvas*[nPtBins];
  TCanvas** cSet2 = new TCanvas*[nPtBins];
  TCanvas** cSet3 = new TCanvas*[nPtBins];
  
  if(isOneFig) {
    cAllSet1 = new TCanvas("cAllSet1","",800,1000);
    cAllSet1->Divide(2,nPtBins/2+1);
    cAllSet2 = new TCanvas("cAllSet2","",800,1000);
    cAllSet2->Divide(2,nPtBins/2+1);
    cAllSet3 = new TCanvas("cAllSet3","",800,1000);
    cAllSet3->Divide(2,nPtBins/2+1);
  }
  else {
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      cSet1[iPt] = new TCanvas(Form("cSet1_Pt%d",iPt),"",1000,800);
      cSet2[iPt] = new TCanvas(Form("cSet2_Pt%d",iPt),"",1000,800);
      cSet3[iPt] = new TCanvas(Form("cSet3_Pt%d",iPt),"",1000,800);
    }
  }
  
  Double_t min1[nPtBins] = {0,0,700,300,300,20,0};
  Double_t min2[nPtBins] = {0,0,0,0,0,0,0};
  Double_t min3[nPtBins] = {0,40,0,0,0,0,0};
  Double_t max1[nPtBins] = {600,700,1800,1000,900,200,70};
  Double_t max2[nPtBins] = {800,1200,1000,670,700,450,160};
  Double_t max3[nPtBins] = {140,200,100,90,50,45,22};

  Double_t x1Info1[nPtBins] = {0.35,0.35,0.35,0.35,0.35,0.35,0.25};
  Double_t y1Info1[nPtBins] = {0.73,0.73,0.73,0.73,0.73,0.73,0.73};
  Double_t x2Info1[nPtBins] = {0.45,0.45,0.45,0.45,0.45,0.45,0.45};
  Double_t y2Info1[nPtBins] = {0.88,0.88,0.88,0.88,0.88,0.88,0.88};

  Double_t x1Info2[nPtBins] = {0.4,0.3,0.4,0.25,0.4,0.4,0.4};
  Double_t y1Info2[nPtBins] = {0.15,0.15,0.15,0.15,0.15,0.15,0.15};
  Double_t x2Info2[nPtBins] = {0.6,0.6,0.6,0.6,0.6,0.6,0.6};
  Double_t y2Info2[nPtBins] = {0.30,0.30,0.30,0.30,0.30,0.30,0.30};

  Double_t x1Info3[nPtBins] = {0.56,0.56,0.56,0.56,0.56,0.56,0.56};
  Double_t y1Info3[nPtBins] = {0.6,0.53,0.6,0.53,0.53,0.6,0.73};
  Double_t x2Info3[nPtBins] = {0.89,0.89,0.89,0.89,0.89,0.89,0.89};
  Double_t y2Info3[nPtBins] = {0.75,0.68,0.75,0.68,0.68,0.75,0.88};

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    TCanvas* c1=(TCanvas*)PlotFitMass3Set(3,PtLims[iPt],PtLims[iPt+1],Form("Fitter_%d",iPt),"cFitter",1,1,min1[iPt],min2[iPt],min3[iPt],max1[iPt],max2[iPt],max3[iPt],x1Info1[iPt],x2Info1[iPt],y1Info1[iPt],y2Info1[iPt],x1Info2[iPt],x2Info2[iPt],y1Info2[iPt],y2Info2[iPt],kTRUE,x1Info3[iPt],x2Info3[iPt],y1Info3[iPt],y2Info3[iPt]);
    if(isOneFig)
      cAllSet1->cd(iPt+1);
    else
      cSet1[iPt]->cd();
    c1->DrawClonePad();
    TCanvas* c2=(TCanvas*)PlotFitMass3Set(3,PtLims[iPt],PtLims[iPt+1],Form("Fitter_%d",iPt),"cFitter",2,1,min1[iPt],min2[iPt],min3[iPt],max1[iPt],max2[iPt],max3[iPt],x1Info1[iPt],x2Info1[iPt],y1Info1[iPt],y2Info1[iPt],x1Info2[iPt],x2Info2[iPt],y1Info2[iPt],y2Info2[iPt],kTRUE,x1Info3[iPt],x2Info3[iPt],y1Info3[iPt],y2Info3[iPt]);
    if(isOneFig)
      cAllSet2->cd(iPt+1);
    else
      cSet2[iPt]->cd();      
    c2->DrawClonePad();
    TCanvas* c3=(TCanvas*)PlotFitMass3Set(3,PtLims[iPt],PtLims[iPt+1],Form("Fitter_%d",iPt),"cFitter",3,1,min1[iPt],min2[iPt],min3[iPt],max1[iPt],max2[iPt],max3[iPt],x1Info1[iPt],x2Info1[iPt],y1Info1[iPt],y2Info1[iPt],x1Info2[iPt],x2Info2[iPt],y1Info2[iPt],y2Info2[iPt],kTRUE,x1Info3[iPt],x2Info3[iPt],y1Info3[iPt],y2Info3[iPt]);
    if(isOneFig)
      cAllSet3->cd(iPt+1);
    else
      cSet3[iPt]->cd();
    c3->DrawClonePad();

    if(!isOneFig) {
      cSet1[iPt]->SaveAs(Form("MassFitSet1_Pt%d.eps",iPt));
      cSet2[iPt]->SaveAs(Form("MassFitSet2_Pt%d.eps",iPt));
      cSet3[iPt]->SaveAs(Form("MassFitSet3_Pt%d.eps",iPt));

      delete cSet1[iPt];
      delete cSet2[iPt];
      delete cSet3[iPt];
    }
  }
  
  if(isOneFig) {
    cAllSet1->SaveAs("MassFitsSet1.eps");
    cAllSet2->SaveAs("MassFitsSet2.eps");
    cAllSet3->SaveAs("MassFitsSet3.eps");
    delete cAllSet1;
    delete cAllSet2;
    delete cAllSet3;
  }
  else {
    delete[] cSet1;
    delete[] cSet2;
    delete[] cSet3;
  }
}

TCanvas* PlotFitMass3Set(Double_t nSigma,
                         Double_t pTmin, Double_t pTmax,
                         TString massfilename,
                         TString masscanvasname,
                         Int_t num,
                         Int_t titletype,
                         Double_t min1,Double_t min2, Double_t min3,
                         Double_t max1,Double_t max2, Double_t max3,
                         Double_t x1Info1, Double_t x2Info1, Double_t y1Info1, Double_t y2Info1,
                         Double_t x1Info2, Double_t x2Info2, Double_t y1Info2, Double_t y2Info2,
                         Bool_t DrawInfo3,
                         Double_t x1Info3, Double_t x2Info3, Double_t y1Info3, Double_t y2Info3,
                         Bool_t Draw3Sets,
                         Double_t xLatex1, Double_t xLatex2, Double_t xLatex3,
                         Double_t yLatex1, Double_t yLatex2, Double_t yLatex3) {

  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetTitleSize(0.06,"t");

  TFile inmassfile(Form("%s.root",massfilename.Data()),"UPDATE");
  TCanvas* c=(TCanvas*)inmassfile.Get(masscanvasname.Data());
  TList* l=(TList*)c->GetListOfPrimitives();

  TH1F** h = new TH1F*[3];
  TF1** fs = new TF1*[3];
  TF1** fb = new TF1*[3];
  TPaveText** info1 = new TPaveText*[3];
  TPaveText** info2 = new TPaveText*[3];
  TPaveText** info3 = new TPaveText*[3];
  TLatex** latex = new TLatex*[3];
  
  for(Int_t iSet=0; iSet<3; iSet++) {
    h[iSet] = (TH1F*)c->GetPad(iSet+1)->GetPrimitive("fhistoInvMass");
    fs[iSet] = (TF1*)c->GetPad(iSet+1)->GetPrimitive("funcmass");
    fb[iSet] = (TF1*)c->GetPad(iSet+1)->GetPrimitive("funcbkgFullRange");
    h[iSet]->SetDirectory(0);
  }
  inmassfile.Close();

  for(Int_t iSet=0; iSet<3; iSet++) {
    latex[iSet] = new TLatex();
    latex[iSet]->SetTextFont(132);
    latex[iSet]->SetTextSize(0.06);
    
    info1[iSet] = new TPaveText(x1Info1,y1Info1,x2Info1,y2Info1,"NDC");
    info2[iSet] = new TPaveText(x1Info2,y1Info2,x2Info2,y2Info2,"NDC");
    info3[iSet] = new TPaveText(x1Info3,y1Info3,x2Info3,y2Info3,"NDC");
    
    Double_t mean = fs[iSet]->GetParameter(3);
    Double_t sigma = fs[iSet]->GetParameter(4);
    Double_t errmean = fs[iSet]->GetParError(3);
    Double_t errsigma = fs[iSet]->GetParError(4);
    Double_t ints=fs[iSet]->Integral(mean-nSigma*sigma,mean+nSigma*sigma)/h[iSet]->GetBinWidth(4);
    Double_t intb=fb[iSet]->Integral(mean-nSigma*sigma,mean+nSigma*sigma)/h[iSet]->GetBinWidth(2);
    Double_t signal = ints-intb;
    Double_t signalerr = fs[iSet]->GetParError(fs[iSet]->GetNpar()-3)/fs[iSet]->GetParameter(fs[iSet]->GetNpar()-3)*signal;
    Double_t bkg = intb;
    Double_t bkgerr = fb[iSet]->GetParError(0)/fb[iSet]->GetParameter(0)*bkg;
    Double_t significance = signal/TMath::Sqrt(signal+bkg);
    Double_t significanceerr = significance*TMath::Sqrt((signalerr*signalerr+bkgerr*bkgerr)/(4.*(signal+bkg)*(signal+bkg))+(bkg/(signal+bkg))*(signalerr*signalerr)/signal/signal);
    Double_t signaloverbkg = signal/bkg;

    info3[iSet]->Clear();
    info3[iSet]->SetTextSize(0.06);
    info3[iSet]->SetBorderSize(0);
    info3[iSet]->SetFillStyle(0);
    info3[iSet]->SetTextColor(kBlue);
    info3[iSet]->SetTextFont(132);
    info3[iSet]->AddText(Form("#mu = %.3f #pm %.3f",mean,errmean));
    info3[iSet]->AddText(Form("#sigma = %.3f #pm %.3f",sigma,errsigma));

    info1[iSet]->Clear();
    info1[iSet]->SetTextSize(0.06);
    info1[iSet]->SetBorderSize(0);
    info1[iSet]->SetTextFont(132);
    info1[iSet]->SetFillStyle(0);
    info1[iSet]->AddText(Form("S (%0.f#sigma) = %.0f #pm %.0f",nSigma,signal,signalerr));
    info1[iSet]->AddText(Form("B (%0.f#sigma) = %.0f #pm %.0f",nSigma,bkg,bkgerr));

    info2[iSet]->Clear();
    info2[iSet]->SetTextSize(0.06);
    info2[iSet]->SetBorderSize(0);
    info2[iSet]->SetTextFont(132);
    info2[iSet]->SetFillStyle(0);
    info2[iSet]->AddText(Form("Signif. (%0.f#sigma) = %.1f #pm %.1f",nSigma,significance,significanceerr));
    info2[iSet]->AddText(Form("S/B (%0.f#sigma) = %.4f",nSigma,signaloverbkg));
  }

  Double_t binwidth = h[0]->GetBinWidth(10)*1000;//in MeV

  TCanvas** cMass = new TCanvas*[3];
  for(Int_t iSet=0; iSet<3; iSet++) {
    cMass[iSet] = new TCanvas(Form("cMass%d",iSet+1),"",1200,900);
  }

  TCanvas* cMass3Set = new TCanvas("cMass3Set","",1200,600);
  cMass3Set->Divide(3,1);

  TString titname[3];
  if(titletype==0) {
    titname[0] = "Prompt-enhanced";
    titname[1] = "Mixed";
    titname[2] = "Feed-down-enhanced";
  }
  else {
    titname[0] = Form("%0.f < #it{p}_{T} < %0.f GeV/c",pTmin,pTmax);
    titname[1] = Form("%0.f < #it{p}_{T} < %0.f GeV/c",pTmin,pTmax);
    titname[2] = Form("%0.f < #it{p}_{T} < %0.f GeV/c",pTmin,pTmax);
  }
  
  for(Int_t iSet=0; iSet<3; iSet++) {
    h[iSet]->SetLineColor(kBlack);
    h[iSet]->SetStats(0);
    h[iSet]->SetTitle(titname[iSet].Data());
    h[iSet]->GetXaxis()->SetTitleSize(0.06);
    h[iSet]->GetYaxis()->SetTitleSize(0.06);
    h[iSet]->GetXaxis()->SetTitleOffset(1.);
    h[iSet]->GetYaxis()->SetTitleOffset(1.5);
    h[iSet]->GetXaxis()->SetLabelSize(0.05);
    h[iSet]->GetYaxis()->SetLabelSize(0.05);
    h[iSet]->GetYaxis()->SetTitle(Form("Entries/(%0.f Mev/c^{2})",binwidth));

    if(iSet==0)
      h[iSet]->GetYaxis()->SetRangeUser(min1,max1);
    else if(iSet==1)
      h[iSet]->GetYaxis()->SetRangeUser(min2,max2);
    else
      h[iSet]->GetYaxis()->SetRangeUser(min3,max3);

    cMass[iSet]->cd();
    cMass[iSet]->Clear();
    h[iSet]->Draw("E");
    fs[iSet]->Draw("same");
    info1[iSet]->Draw("same");
    info2[iSet]->Draw("same");
    if(DrawInfo3)
      info3[iSet]->Draw("same");

    cMass[iSet]->Update();

    if(Draw3Sets) {
      cMass3Set->cd(iSet+1);
      h[iSet]->Draw("E");
      fs[iSet]->Draw("same");
      info1[iSet]->Draw("same");
      info2[iSet]->Draw("same");
      if(DrawInfo3)
        info3[iSet]->Draw("same"); 
      
      if(titletype==0) {
        if(iSet==0) 
          latex[iSet]->DrawLatex(xLatex1,yLatex1,Form("%0.f < #it{p}_{T} < %0.f GeV/c",pTmin,pTmax));
        else if(iSet==1) 
          latex[iSet]->DrawLatex(xLatex2,yLatex2,Form("%0.f < #it{p}_{T} < %0.f GeV/c",pTmin,pTmax)); 
        else
          latex[iSet]->DrawLatex(xLatex3,yLatex3,Form("%0.f < #it{p}_{T} < %0.f GeV/c",pTmin,pTmax));
      }
    }
  }
 
  cMass3Set->SaveAs(Form("Mass3Set_%0.f-%0.f.eps",pTmin,pTmax));
    
  return cMass[num-1];
}
