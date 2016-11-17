#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TString.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TASImage.h>
#include <TPad.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>

#endif

const Int_t colors[] = {kRed+1,kBlue+1};

Int_t PlotMassFitsInOutOfPlane(TString filename="InvMassDeltaPhi_fs_Topod0Cut_VZERO_EP.root", TString masscanvasname="cinvmassdeltaphifs", Double_t nSigma=3) {

  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetTitleSize(0.07,"t");
  
  TFile* inmassfile=TFile::Open(filename,"READ");
  TCanvas* c1=0x0;
  TList* l1=0x0;
  if(inmassfile) {
    c1=(TCanvas*)inmassfile->Get(masscanvasname.Data());
    if(c1) {l1=(TList*)c1->GetListOfPrimitives();}
    else {return 2;}
    inmassfile->Close();
  }
  else {return 1;}
  
  const Int_t nPtBins = l1->GetSize()/2;
  cout << nPtBins << endl;
  
  TH1F** hMassInPlane = new TH1F*[nPtBins];
  TF1** fsInPlane = new TF1*[nPtBins];
  TF1** fbInPlane = new TF1*[nPtBins];
  TPaveText** infoInPlane1 = new TPaveText*[nPtBins];
  TPaveText** infoInPlane2 = new TPaveText*[nPtBins];
  TH1F** hMassOutOfPlane = new TH1F*[nPtBins];
  TF1** fsOutOfPlane = new TF1*[nPtBins];
  TF1** fbOutOfPlane = new TF1*[nPtBins];
  TPaveText** infoOutOfPlane1 = new TPaveText*[nPtBins];
  TPaveText** infoOutOfPlane2 = new TPaveText*[nPtBins];
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    hMassInPlane[iPt] = (TH1F*)c1->GetPad(iPt+1)->GetPrimitive("fhistoInvMass");
    fsInPlane[iPt] = (TF1*)c1->GetPad(iPt+1)->GetPrimitive("funcmass");
    fbInPlane[iPt] = (TF1*)c1->GetPad(iPt+1)->GetPrimitive("funcbkgRecalc");
    hMassInPlane[iPt]->SetDirectory(0);
    hMassOutOfPlane[iPt] = (TH1F*)c1->GetPad(nPtBins+iPt+1)->GetPrimitive("fhistoInvMass");
    fsOutOfPlane[iPt] = (TF1*)c1->GetPad(nPtBins+iPt+1)->GetPrimitive("funcmass");
    fbOutOfPlane[iPt] = (TF1*)c1->GetPad(nPtBins+iPt+1)->GetPrimitive("funcbkgRecalc");
    hMassOutOfPlane[iPt]->SetDirectory(0);
  }
  inmassfile->Close();
  
  TCanvas* cMassFits = new TCanvas("cMassFits","",10,10,1920,1080);
  cMassFits->Divide(4,3);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    infoInPlane1[iPt] = new TPaveText(0.25,0.58,0.5,0.70,"NDC");
    infoInPlane2[iPt] = new TPaveText(0.25,0.72,0.5,0.84,"NDC");
    
    Double_t meanInPlane = fsInPlane[iPt]->GetParameter(3);
    Double_t sigmaInPlane = fsInPlane[iPt]->GetParameter(4);
    Double_t errmeanInPlane = fsInPlane[iPt]->GetParError(3);
    Double_t errsigmaInPlane = fsInPlane[iPt]->GetParError(4);
    Double_t intsInPlane=fsInPlane[iPt]->Integral(meanInPlane-nSigma*sigmaInPlane,meanInPlane+nSigma*sigmaInPlane)/hMassInPlane[iPt]->GetBinWidth(4);
    Double_t intbInPlane=fbInPlane[iPt]->Integral(meanInPlane-nSigma*sigmaInPlane,meanInPlane+nSigma*sigmaInPlane)/hMassInPlane[iPt]->GetBinWidth(2);
    Double_t signalInPlane = intsInPlane-intbInPlane;
    Double_t signalerrInPlane = fsInPlane[iPt]->GetParError(fsInPlane[iPt]->GetNpar()-3)/fsInPlane[iPt]->GetParameter(fsInPlane[iPt]->GetNpar()-3)*signalInPlane;
    Double_t bkgInPlane = intbInPlane;
    Double_t bkgerrInPlane = fbInPlane[iPt]->GetParError(0)/fbInPlane[iPt]->GetParameter(0)*bkgInPlane;
    Double_t significanceInPlane = signalInPlane/TMath::Sqrt(signalInPlane+bkgInPlane);
    Double_t significanceerrInPlane = significanceInPlane*TMath::Sqrt((signalerrInPlane*signalerrInPlane+bkgerrInPlane*bkgerrInPlane)/(4.*(signalInPlane+bkgInPlane)*(signalInPlane+bkgInPlane))+(bkgInPlane/(signalInPlane+bkgInPlane))*(signalerrInPlane*signalerrInPlane)/signalInPlane/signalInPlane);
    Double_t signaloverbkgInPlane = signalInPlane/bkgInPlane;
    
    infoInPlane1[iPt]->Clear();
    infoInPlane1[iPt]->SetTextSize(0.045);
    infoInPlane1[iPt]->SetBorderSize(0);
    infoInPlane1[iPt]->SetTextFont(132);
    infoInPlane1[iPt]->SetFillStyle(0);
    infoInPlane1[iPt]->AddText(Form("S (%0.f#sigma) = %.0f #pm %.0f",nSigma,signalInPlane,signalerrInPlane));
    infoInPlane1[iPt]->AddText(Form("B (%0.f#sigma) = %.0f #pm %.0f",nSigma,bkgInPlane,bkgerrInPlane));
    
    infoInPlane2[iPt]->Clear();
    infoInPlane2[iPt]->SetTextSize(0.045);
    infoInPlane2[iPt]->SetBorderSize(0);
    infoInPlane2[iPt]->SetTextFont(132);
    infoInPlane2[iPt]->SetFillStyle(0);
    infoInPlane2[iPt]->AddText(Form("Signif. (%0.f#sigma) = %.1f #pm %.1f",nSigma,significanceInPlane,significanceerrInPlane));
    infoInPlane2[iPt]->AddText(Form("S/B (%0.f#sigma) = %.4f",nSigma,signaloverbkgInPlane));

    infoOutOfPlane1[iPt] = new TPaveText(0.65,0.58,0.89,0.70,"NDC");
    infoOutOfPlane2[iPt] = new TPaveText(0.65,0.72,0.89,0.84,"NDC");
    
    Double_t meanOutOfPlane = fsOutOfPlane[iPt]->GetParameter(3);
    Double_t sigmaOutOfPlane = fsOutOfPlane[iPt]->GetParameter(4);
    Double_t errmeanOutOfPlane = fsOutOfPlane[iPt]->GetParError(3);
    Double_t errsigmaOutOfPlane = fsOutOfPlane[iPt]->GetParError(4);
    Double_t intsOutOfPlane=fsOutOfPlane[iPt]->Integral(meanOutOfPlane-nSigma*sigmaOutOfPlane,meanOutOfPlane+nSigma*sigmaOutOfPlane)/hMassOutOfPlane[iPt]->GetBinWidth(4);
    Double_t intbOutOfPlane=fbOutOfPlane[iPt]->Integral(meanOutOfPlane-nSigma*sigmaOutOfPlane,meanOutOfPlane+nSigma*sigmaOutOfPlane)/hMassOutOfPlane[iPt]->GetBinWidth(2);
    Double_t signalOutOfPlane = intsOutOfPlane-intbOutOfPlane;
    Double_t signalerrOutOfPlane = fsOutOfPlane[iPt]->GetParError(fsOutOfPlane[iPt]->GetNpar()-3)/fsOutOfPlane[iPt]->GetParameter(fsOutOfPlane[iPt]->GetNpar()-3)*signalOutOfPlane;
    Double_t bkgOutOfPlane = intbOutOfPlane;
    Double_t bkgerrOutOfPlane = fbOutOfPlane[iPt]->GetParError(0)/fbOutOfPlane[iPt]->GetParameter(0)*bkgOutOfPlane;
    Double_t significanceOutOfPlane = signalOutOfPlane/TMath::Sqrt(signalOutOfPlane+bkgOutOfPlane);
    Double_t significanceerrOutOfPlane = significanceOutOfPlane*TMath::Sqrt((signalerrOutOfPlane*signalerrOutOfPlane+bkgerrOutOfPlane*bkgerrOutOfPlane)/(4.*(signalOutOfPlane+bkgOutOfPlane)*(signalOutOfPlane+bkgOutOfPlane))+(bkgOutOfPlane/(signalOutOfPlane+bkgOutOfPlane))*(signalerrOutOfPlane*signalerrOutOfPlane)/signalOutOfPlane/signalOutOfPlane);
    Double_t signaloverbkgOutOfPlane = signalOutOfPlane/bkgOutOfPlane;
    
    infoOutOfPlane1[iPt]->Clear();
    infoOutOfPlane1[iPt]->SetTextSize(0.05);
    infoOutOfPlane1[iPt]->SetBorderSize(0);
    infoOutOfPlane1[iPt]->SetTextFont(132);
    infoOutOfPlane1[iPt]->SetFillStyle(0);
    infoOutOfPlane1[iPt]->AddText(Form("S (%0.f#sigma) = %.0f #pm %.0f",nSigma,signalOutOfPlane,signalerrOutOfPlane));
    infoOutOfPlane1[iPt]->AddText(Form("B (%0.f#sigma) = %.0f #pm %.0f",nSigma,bkgOutOfPlane,bkgerrOutOfPlane));
    
    infoOutOfPlane2[iPt]->Clear();
    infoOutOfPlane2[iPt]->SetTextSize(0.05);
    infoOutOfPlane2[iPt]->SetBorderSize(0);
    infoOutOfPlane2[iPt]->SetTextFont(132);
    infoOutOfPlane2[iPt]->SetFillStyle(0);
    infoOutOfPlane2[iPt]->AddText(Form("Signif. (%0.f#sigma) = %.1f #pm %.1f",nSigma,significanceOutOfPlane,significanceerrOutOfPlane));
    infoOutOfPlane2[iPt]->AddText(Form("S/B (%0.f#sigma) = %.4f",nSigma,signaloverbkgOutOfPlane));
  }

  Double_t binwidth = hMassInPlane[0]->GetBinWidth(10)*1000;//in MeV
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    hMassInPlane[iPt]->GetFunction("funcbkgFullRange")->SetBit(TF1::kNotDraw);
    hMassOutOfPlane[iPt]->GetFunction("funcbkgFullRange")->SetBit(TF1::kNotDraw);
    hMassInPlane[iPt]->GetFunction("funcbkgRecalc")->SetBit(TF1::kNotDraw);
    hMassOutOfPlane[iPt]->GetFunction("funcbkgRecalc")->SetBit(TF1::kNotDraw);
    hMassInPlane[iPt]->SetMarkerColor(colors[1]);
    hMassInPlane[iPt]->SetMarkerSize(0.8);
    hMassInPlane[iPt]->SetLineColor(colors[1]);
    fsInPlane[iPt]->SetLineColor(colors[1]);
    fbInPlane[iPt]->SetLineColor(colors[1]);
    infoInPlane1[iPt]->SetTextColor(colors[1]);
    infoInPlane2[iPt]->SetTextColor(colors[1]);
    hMassInPlane[iPt]->SetStats(0);
    hMassInPlane[iPt]->GetXaxis()->SetTitleSize(0.06);
    hMassInPlane[iPt]->GetYaxis()->SetTitleSize(0.06);
    hMassInPlane[iPt]->GetXaxis()->SetTitleOffset(1.);
    hMassInPlane[iPt]->GetYaxis()->SetTitleOffset(1.6);
    hMassInPlane[iPt]->GetXaxis()->SetLabelSize(0.055);
    hMassInPlane[iPt]->GetYaxis()->SetLabelSize(0.055);
    binwidth = hMassOutOfPlane[iPt]->GetBinWidth(15)*1000;
    hMassInPlane[iPt]->GetYaxis()->SetTitle(Form("Entries/(%0.f Mev/c^{2})",binwidth));
    hMassInPlane[iPt]->GetXaxis()->SetTitle("M_{K#pi#pi} (GeV/c^{2})");
    TString title=hMassInPlane[iPt]->GetTitle();
    if(title.Contains(", #phi0")) {title.ReplaceAll(", #phi0","");}
    hMassInPlane[iPt]->SetTitle(title.Data());
    
    hMassOutOfPlane[iPt]->SetMarkerColor(colors[0]);
    hMassOutOfPlane[iPt]->SetMarkerSize(0.8);
    hMassOutOfPlane[iPt]->SetLineColor(colors[0]);
    fsOutOfPlane[iPt]->SetLineColor(colors[0]);
    fbOutOfPlane[iPt]->SetLineColor(colors[0]);
    infoOutOfPlane1[iPt]->SetTextColor(colors[0]);
    infoOutOfPlane2[iPt]->SetTextColor(colors[0]);
    
    Double_t min = hMassInPlane[iPt]->GetMinimum()*0.8;
    if(hMassInPlane[iPt]->GetMinimum()>hMassOutOfPlane[iPt]->GetMinimum())
      min = hMassOutOfPlane[iPt]->GetMinimum()*0.8;
    Double_t max = hMassInPlane[iPt]->GetMaximum()*2.;
    if(hMassInPlane[iPt]->GetMaximum()<hMassOutOfPlane[iPt]->GetMaximum())
      max = hMassOutOfPlane[iPt]->GetMaximum()*2.;
  
    hMassInPlane[iPt]->GetYaxis()->SetRangeUser(min,max);

    cMassFits->cd(iPt+1);
    hMassInPlane[iPt]->Draw("E");
    hMassOutOfPlane[iPt]->Draw("Esame");
    fsInPlane[iPt]->Draw("same");
    fsOutOfPlane[iPt]->Draw("same");
    //fbInPlane[iPt]->Draw("same");
    //fbOutOfPlane[iPt]->Draw("same");
    infoInPlane1[iPt]->Draw("same");
    infoInPlane2[iPt]->Draw("same");
    infoOutOfPlane1[iPt]->Draw("same");
    infoOutOfPlane2[iPt]->Draw("same");
  }
  
  cMassFits->SaveAs("MassFitsInOutOfPlane.pdf");
  return 0;
}
