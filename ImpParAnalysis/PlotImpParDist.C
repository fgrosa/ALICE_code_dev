const Int_t nPtBins = 7;
const Int_t nPtLims = nPtBins+1;
Double_t PtLims[nPtLims] = {2,3,4,5,6,8,12,16};

void PlotPrompt();
void PlotFD();
void PlotBkg();
TCanvas* PlotImpParDist(Double_t min=-250, Double_t max=250, Double_t pTmin = 4, Double_t pTmax=5, TString filename="ImpParPrompt_4-5", TString canvasname = "cPrompt", TString histoname="hMassPtImpParPrompt_proj_2", TString funcname = "ImpParPromptFunc");

void PlotPrompt() {
  TCanvas* cAll = new TCanvas("cAll","",800,1000);
  cAll->Divide(2,nPtBins/2+1);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    TCanvas* c=(TCanvas*)PlotImpParDist(-300,300,PtLims[iPt],PtLims[iPt+1],Form("ImpParPrompt_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"cPrompt","hMassPtImpParPrompt_proj_2","ImpParPromptFunc");
    cAll->cd(iPt+1);
    c->DrawClonePad();
  }

  cAll->SaveAs("PromptFits.eps");
  delete cAll;
}

void PlotFD() {
  TCanvas* cAllTrue = new TCanvas("cAllTrue","",800,1000);
  cAllTrue->Divide(2,nPtBins/2+1);
  TCanvas* cAllReco = new TCanvas("cAllReco","",800,1000);
  cAllReco->Divide(2,nPtBins/2+1);

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    TCanvas* cTrue=(TCanvas*)PlotImpParDist(-500,500,PtLims[iPt],PtLims[iPt+1],Form("ImpParTrueFD_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"cTrueFD","hMassPtImpParTrueBfeed_proj_2","ImpParTrueFDFunc");
    cAllTrue->cd(iPt+1);
    cTrue->DrawClonePad();
    TCanvas* cReco=(TCanvas*)PlotImpParDist(-500,500,PtLims[iPt],PtLims[iPt+1],Form("ImpParRecoFD_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"cRecoFD","hMassPtImpParBfeed_proj_2","ImpParRecoFDFunc");
    cAllReco->cd(iPt+1);
    cReco->DrawClonePad();
  }

  cAllTrue->SaveAs("TrueFDFits.eps");
  cAllReco->SaveAs("RecoFDFits.eps");
  delete cAllTrue;
  delete cAllReco;

}

void PlotBkg() {
  TCanvas* cAll = new TCanvas("cAll","",800,1000);
  cAll->Divide(2,nPtBins/2+1);
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    TCanvas* c=(TCanvas*)PlotImpParDist(-600,600,PtLims[iPt],PtLims[iPt+1],Form("ImpParBkg_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"cBkg","hImpParBkg","ImpParBkgFunc");
    cAll->cd(iPt+1);
    c->DrawClonePad();   
  }

  cAll->SaveAs("BkgFits.eps");
  delete cAll;
}

TCanvas* PlotImpParDist(Double_t min,Double_t max,Double_t pTmin,Double_t pTmax,TString filename,TString canvasname,TString histoname,TString funcname) {

  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);  
  
  cout << Form("%s.root",filename.Data()) << endl;
  
  TFile infile(Form("%s.root",filename.Data()),"UPDATE");
  TCanvas* cInput = (TCanvas*)infile.Get(canvasname.Data());
  TH1F* hImpPar = (TH1F*)cInput->GetPrimitive(histoname.Data());
  TF1* fImpPar = (TF1*)cInput->GetPrimitive(funcname.Data());

  Double_t ChiS =  fImpPar->GetChisquare();
  if(ChiS==0) {
    hImpPar->SetStats(0);
  }
    
  TLatex lat;
  lat.SetTextFont(132);
  lat.SetTextSize(0.06);
  
  TCanvas* c = new TCanvas("c","",1200,1000);
  c->SetLogy();
  hImpPar->GetXaxis()->SetRangeUser(min,max);
  hImpPar->GetYaxis()->SetRangeUser(0.1,hImpPar->GetMaximum()*7);
  hImpPar->Draw("E1");
  hImpPar->GetXaxis()->SetTitle("Imp Par XY (#mum)");
  hImpPar->GetYaxis()->SetTitle(Form("Entries/(%0.f #mum)",hImpPar->GetBinWidth(10)));
  hImpPar->GetXaxis()->SetTitleSize(0.05);
  hImpPar->GetYaxis()->SetTitleSize(0.05);
  hImpPar->GetXaxis()->SetLabelSize(0.05);
  hImpPar->GetYaxis()->SetLabelSize(0.05);
  hImpPar->GetXaxis()->SetTitleOffset(1.2);
  hImpPar->GetYaxis()->SetTitleOffset(1.2);
  hImpPar->GetXaxis()->SetNdivisions(508,kTRUE);
  fImpPar->Draw("same");
  lat.DrawLatex(min+TMath::Abs(min)*0.1, hImpPar->GetMaximum()*0.3,Form("%0.f < #it{p}_{T} < %0.f GeV/c", pTmin, pTmax));

  c->SaveAs(Form("%s.eps",filename.Data()));

  return c;
}
