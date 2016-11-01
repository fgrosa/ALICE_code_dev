#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <vector>
#include <TInterpreter.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TVirtualFitter.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <THnSparse.h>
#include <TString.h>
#include <TLine.h>
#include <TList.h>
#include <TF1.h>

#endif

void PlotPtBDistAndEff();
TH1F* GetFONLL(TString filename="FONLL_Bmesons_5TeV_80bin");
Double_t PtBWeightsFromFONLL5overLHC13d3(Double_t *x, Double_t *pars);

void PlotPtBDistAndEff() {

  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetPadBottomMargin(0.4);
  gStyle->SetPadLeftMargin(0.18);
  
  const Int_t nPtBins = 10;
  const Int_t nPtLims = nPtBins+1;
  Double_t PtLims[nPtLims] = {1,2,3,4,5,6,7,8,12,16,24};
  
  //load Reco and MCGenAcc histos
  TFile MCfile("AnalysisResultspPbMC.root","UPDATE");
  TDirectoryFile* dir=(TDirectoryFile*)MCfile.Get("PWG3_D2H_InvMassDplus");
  TList* list=(TList*)dir->Get("coutputDplus_ImpParpPbMC0100");
  THnSparseF* hMassPtImpParRecoFD=(THnSparseF*)list->FindObject("hMassPtImpParBfeed");
  hMassPtImpParRecoFD->GetAxis(7)->SetRange(3,3);
  THnSparseF* hPtYGenAccFD=(THnSparseF*)list->FindObject("hMCAccBFeed");
  TH1F* hPtGenAcc=(TH1F*)hPtYGenAccFD->Projection(0);
  TH1F* hPtReco=(TH1F*)hMassPtImpParRecoFD->Projection(1);
  hPtGenAcc->SetDirectory(0);
  hPtReco->SetDirectory(0);
  hPtGenAcc->Sumw2();
  hPtReco->Sumw2(); 
  
  //load ESDGenLimAcc histos for B0, Bplus and Bs
  TFile ESDfile("AnalysisResultspPbMCESD.root","UPDATE");
  TDirectory *ESDdir = (TDirectory*)ESDfile.Get("HFMCCheck");
  TList *ESDlist = (TList*)ESDdir->Get("clistHFMCCheck");
  TH2F *hPtYB0 = (TH2F*)ESDlist->FindObject("hyptB0AllDecay");
  TH2F *hPtYBplus = (TH2F*)ESDlist->FindObject("hyptBplusAllDecay");
  TH2F *hPtYBs = (TH2F*)ESDlist->FindObject("hyptBsAllDecay");
  hPtYB0->SetDirectory(0);
  hPtYBplus->SetDirectory(0);
  hPtYBs->SetDirectory(0);

  TH1F *hPtB0 = (TH1F*)hPtYB0->ProjectionX();
  TH1F *hPtBplus = (TH1F*)hPtYBplus->ProjectionX();
  TH1F *hPtBs = (TH1F*)hPtYBs->ProjectionX();
  hPtB0->SetDirectory(0);
  hPtBplus->SetDirectory(0);
  hPtBs->SetDirectory(0);
  hPtB0->Sumw2();
  hPtBplus->Sumw2();
  hPtBs->Sumw2();
  
  const Int_t nbins = hPtB0->GetNbinsX();
  Double_t ptmin = hPtB0->GetBinLowEdge(1);
  Double_t ptmax = hPtB0->GetBinLowEdge(hPtB0->GetNbinsX())+hPtB0->GetBinWidth(hPtB0->GetNbinsX());
  Double_t step = (ptmax-ptmin)/nbins;
  const Int_t nlims = nbins+1;
  Double_t ptlims[nlims];
  
  for(Int_t iLim=0; iLim<nlims; iLim++) {
    ptlims[iLim] = iLim*step; 
  }
  
  TH1F *hPtBGenLimAcc= new TH1F("hPtBGenLimAcc","",nbins,ptmin,ptmax);

  for(Int_t iBin=0; iBin<hPtB0->GetNbinsX(); iBin++)
    hPtBGenLimAcc->SetBinContent(iBin+1,hPtB0->GetBinContent(iBin+1)+hPtBplus->GetBinContent(iBin+1)+hPtBs->GetBinContent(iBin+1));
  hPtBGenLimAcc->SetDirectory(0);
  hPtBGenLimAcc->Sumw2();
  hPtBGenLimAcc->Scale(1./hPtBGenLimAcc->Integral());
  
  //load FONLL histo
  TH1F* hFONLL=GetFONLL();
  hFONLL->SetDirectory(0);

  //calculate weights and fit histo
  TH1F* hWeights = new TH1F("hWeights","",nbins,ptmin,ptmax);
  hWeights->Divide(hFONLL,hPtBGenLimAcc,1.,1.,"");
  hWeights->SetDirectory(0);
  TH1F *hWeightsCopy=(TH1F*)hWeights->Clone();
  hWeightsCopy->SetDirectory(0);
  TF1* fFuncWeight = new TF1("fFuncWeight",PtBWeightsFromFONLL5overLHC13d3,0,40,8);
  fFuncWeight->SetParameters(0.983543,2.72857,1.49631,-82.4067,51.359,5.05355,0.00136537,2.13765e-05);
  
  //reweight MCGenAcc and Reco histograms
  TH1F* hPtRecoWeight = new TH1F("hPtRecoWeight","",nbins,ptmin,ptmax);
  hPtRecoWeight->Multiply(hPtReco,hWeights,1.,1.,"");
  hPtRecoWeight->SetDirectory(0);
  hPtRecoWeight->Sumw2();
  TH1F* hPtGenAccWeight=(TH1F*)hPtGenAcc->Rebin(nbins,"hPtGenAccWeight",ptlims);
  hPtGenAccWeight->Multiply(hPtGenAccWeight,hWeights,1.,1.,"");
  hPtGenAccWeight->SetDirectory(0);
  hPtGenAccWeight->Sumw2();
  
  //calculate efficiencies
  TH1F* hPtEff =new TH1F("hEff","",nPtBins,PtLims);
  TH1F* hPtGenAccReb=(TH1F*)hPtGenAcc->Rebin(nPtBins,"hPtGenAccReb",PtLims);
  TH1F* hPtRecoReb=(TH1F*)hPtReco->Rebin(nPtBins,"hPtRecoReb",PtLims);
  hPtGenAccReb->SetDirectory(0);
  hPtRecoReb->SetDirectory(0);
  hPtGenAccReb->Sumw2();
  hPtRecoReb->Sumw2();
  hPtEff->Divide(hPtRecoReb,hPtGenAccReb,1.,1.,"B");
  hPtEff->SetDirectory(0);
  hPtEff->Sumw2();

  TH1F* hPtEffWeight = new TH1F("hPtEffWeight","",nPtBins,PtLims);
  TH1F* hPtGenAccWeightReb=(TH1F*)hPtGenAccWeight->Rebin(nPtBins,"hPtGenAccWeightReb",PtLims);
  TH1F* hPtRecoWeightReb=(TH1F*)hPtRecoWeight->Rebin(nPtBins,"hPtRecoWeightReb",PtLims);
  hPtGenAccWeightReb->SetDirectory(0);
  hPtRecoWeightReb->SetDirectory(0);
  hPtGenAccWeightReb->Sumw2();
  hPtRecoWeightReb->Sumw2();
  hPtEffWeight->Divide(hPtRecoWeightReb,hPtGenAccWeightReb,1.,1.,"B");
  hPtEffWeight->SetDirectory(0);
  hPtEffWeight->Sumw2();
  
  //ratio between efficiency and efficiency reweighted
  TH1F* hRatio = new TH1F("hRatio","",nPtBins,PtLims);
  hRatio->Divide(hPtEffWeight,hPtEff,1.,1.,"");
  hRatio->SetDirectory(0);

  TLegend* l = new TLegend(0.65,0.65,0.85,0.85);
  l->SetTextSize(0.06);
  l->AddEntry(hPtBGenLimAcc,"Pythia","lpe");
  l->AddEntry(hFONLL,"FONLL","lpe");

  TLegend* l2 = new TLegend(0.65,0.45,0.85,0.65);
  l2->SetTextSize(0.045);
  l2->AddEntry(hPtEff,"Not reweighted","lpe");
  l2->AddEntry(hPtEffWeight,"Reweighted","lpe");

  TLatex* lat = new TLatex();
  lat->SetTextSize(0.07);
  lat->SetTextFont(132);
  
  TCanvas* cSpec = new TCanvas("cSpec","",800,900);
  cSpec->Clear();
  //cSpec->SetLogy();
  TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->SetLogy();
  pad1->cd();
  hPtBGenLimAcc->GetXaxis()->SetTitleSize(0.07);
  hPtBGenLimAcc->GetYaxis()->SetTitleSize(0.07);
  hPtBGenLimAcc->GetXaxis()->SetLabelSize(0.06);
  hPtBGenLimAcc->GetYaxis()->SetLabelSize(0.06);
  hPtBGenLimAcc->GetYaxis()->SetTitleOffset(1.2);
  hPtBGenLimAcc->GetXaxis()->SetTitle("#it{p}_{T}^{B} (GeV/c)");
  hPtBGenLimAcc->GetYaxis()->SetTitle("Normalised entires");
  hPtBGenLimAcc->SetLineColor(kRed);
  hWeightsCopy->GetXaxis()->SetTitleSize(0.14);
  hWeightsCopy->GetYaxis()->SetTitleSize(0.13);
  hWeightsCopy->GetYaxis()->SetTitleOffset(0.6);
  hWeightsCopy->GetXaxis()->SetLabelSize(0.12);
  hWeightsCopy->GetYaxis()->SetLabelSize(0.12);
  hPtBGenLimAcc->Draw("E1");
  hFONLL->Draw("E1same");
  l->Draw("same");
  lat->DrawLatex(4,0.00005,"B-mesons decaying to D^{+}");
  cSpec->cd();
  TPad *pad2 = new TPad("pad2","pad2",0.,0.,1,0.35);
  pad2->SetTopMargin(0);
  pad2->Draw();
  pad2->cd();
  hWeightsCopy->GetXaxis()->SetTitle(hPtBGenLimAcc->GetXaxis()->GetTitle());
  hWeightsCopy->GetYaxis()->SetTitle("#frac{FONLL}{Pythia}");
  hWeightsCopy->Draw("E1");

  TCanvas* cWeights = new TCanvas("cWeights","",1200,900);
  cWeights->Clear();
  hWeights->SetStats(0);
  hWeights->Draw("E1");
  TVirtualFitter::SetDefaultFitter("Minuit2");
  hWeights->Fit("fFuncWeight","RE","",0,40);
  hWeights->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hWeights->GetYaxis()->SetTitle("#frac{FONLL}{GenLimAcc}");
  fFuncWeight->Draw("same");
  cWeights->SaveAs("PtBweights.eps");

  TCanvas* cEff = new TCanvas("cEff","",600,900);
  cEff->Clear();
  cEff->Divide(1,2);
  cEff->cd(1);
  hPtEff->SetLineColor(kRed);
  hPtEff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hPtEff->GetYaxis()->SetTitle("Efficiency");
  hPtEff->Draw("E1");
  hPtEffWeight->SetLineColor(kBlue);
  hPtEffWeight->Draw("E1same");
  l2->Draw("same");

  cEff->cd(2);
  hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hRatio->GetYaxis()->SetTitle("#frac{reweighted}{not reweighted}");
  hRatio->Draw("E1");

  cSpec->SaveAs("Bmesons_PtSpetrum.eps");
  cSpec->SaveAs("Bmesons_PtSpetrum.root");
  
  TFile outfile("PtB.root","RECREATE");
  hFONLL->Write();
  hPtBGenLimAcc->Write();
  hWeights->Write();
  fFuncWeight->Write();
  hPtEff->Write();
  hPtEffWeight->Write();
  hRatio->Write();
  outfile.Close();
}

TH1F* GetFONLL(TString filename) {

  ifstream infile(filename.Data());
  if(!infile){
    cerr<<"file "<<filename.Data() <<" does not exist "<<endl;
    return 0x0;
  }

  Double_t pt=0;
  Double_t dsigmadpt=0;
  Double_t min=0;
  Double_t max=0;
  Double_t min_sc=0;
  Double_t max_sc=0;
  Double_t min_mass=0;
  Double_t max_mass=0;
  vector<Double_t> Pt;
  vector<Double_t> Dsigmadpt;
  vector<Double_t> Min;
  vector<Double_t> Max;
  vector<Double_t> Min_sc;
  vector<Double_t> Max_sc;
  vector<Double_t> Min_mass;
  vector<Double_t> Max_mass;

  Int_t nbins=0;
  
  while(infile>>pt>>dsigmadpt>>min>>max>>min_sc>>max_sc>>min_mass>>max_mass) {
    Pt.push_back(pt);
    Dsigmadpt.push_back(dsigmadpt);
    nbins++;
  }

  Double_t step=(Pt[1]-Pt[0])/2; 
  Double_t ptmin=Pt[0]-step;
  Double_t ptmax=Pt[nbins-1]+step;

  TH1F* hFONLL = new TH1F("hFONLL","",nbins,ptmin,ptmax);

  for(Int_t iBin=0; iBin<nbins; iBin++)
    hFONLL->SetBinContent(iBin+1,Dsigmadpt[iBin]);

  hFONLL->Sumw2();
  hFONLL->Scale(1./hFONLL->Integral());

  infile.close();
  return hFONLL;
}

//_________________________________________________________________________
Double_t PtBWeightsFromFONLL5overLHC13d3(Double_t *x, Double_t *pars){
  // weight function from the ratio of the LHC13d3 MC
  // and FONLL calculations for pp data rescaled for A
  Double_t pt = x[0];
  Double_t weight = pars[0]*TMath::Gaus(pt,pars[1],pars[2],kTRUE)+pars[3]+pars[4]*TMath::Log(pars[5]+pars[6]*pt+pars[7]*pt*pt);
  
  return weight;
}
