void PlotKFDist(Int_t iPt=0, TString filename="KFchi_Dplus.root", TString canvasname="KFchi_pT", Int_t iProj=12) {

  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.14);
  
  TFile infile(filename.Data(),"READ");
  TCanvas* cin=(TCanvas*)infile.Get(Form("%s%d",canvasname.Data(),iPt));
  TH1F* hPrompt=(TH1F*)cin->GetPrimitive(Form("fSparsePrompt_proj_%d",iProj));
  TH1F* hFD=(TH1F*)cin->GetPrimitive(Form("fSparseFD_proj_%d",iProj));
  TH1F* hBkg=(TH1F*)cin->GetPrimitive(Form("fSparseAll_proj_%d",iProj));
  hPrompt->GetXaxis()->SetTitle("#chi^{2}/ndf");
  hPrompt->GetXaxis()->SetTitleSize(0.05);
  hPrompt->GetYaxis()->SetTitleSize(0.05);
  hPrompt->GetYaxis()->SetTitleOffset(1.4);
  hPrompt->GetXaxis()->SetLabelSize(0.05);
  hPrompt->GetYaxis()->SetLabelSize(0.05);
  hPrompt->SetLineWidth(2);
  hFD->SetLineWidth(2);
  hFD->SetLineStyle(9);
  hBkg->SetLineWidth(2);
  hBkg->SetLineStyle(2); 
  Double_t y1=0.75;
  Double_t y2=0.89;
    
  TLegend* l = new TLegend(0.35,y1,0.8,y2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextSize(0.05);
  l->AddEntry(hPrompt,"Prompt","l");
  l->AddEntry(hFD,"Feed-down","l");
  l->AddEntry(hBkg,"Background","l");

  TPaveText* latPrompt = new TPaveText(0.5,0.65,0.7,0.72,"NDC");
  latPrompt->SetTextColor(kBlue);
  latPrompt->SetFillStyle(0);
  latPrompt->SetFillColor(0);
  latPrompt->SetBorderSize(0);
  latPrompt->SetTextFont(132);
  latPrompt->SetTextSize(0.05);
  
  TPaveText* latFD = 0x0;
  TPaveText* latBkg = 0x0;
  
  if(iPt<6) {
    latBkg = new TPaveText(0.5,0.50,0.7,0.57,"NDC");
    latFD = new TPaveText(0.54,0.58,0.7,0.65,"NDC");
  }
  else {
    latBkg = new TPaveText(0.48,0.50,0.7,0.57,"NDC");
    latFD = new TPaveText(0.51,0.58,0.7,0.65,"NDC");
  }
  
  latFD->SetTextColor(kRed);
  latFD->SetFillStyle(0);
  latFD->SetFillColor(0);
  latFD->SetBorderSize(0);
  latFD->SetTextFont(132);
  latFD->SetTextSize(0.05);

  latBkg->SetTextColor(kGreen+3);
  latBkg->SetFillStyle(0);
  latBkg->SetFillColor(0);
  latBkg->SetBorderSize(0);
  latBkg->SetTextFont(132);
  latBkg->SetTextSize(0.05);
  
  latPrompt->AddText(Form("<#chi^{2}/ndf>_{prompt} = %0.2f #pm %0.2f",hPrompt->GetMean(),hPrompt->GetMeanError()));
  if(iPt<6) {
    latBkg->AddText(Form("<#chi^{2}/ndf>_{bkg}= %0.3f #pm %0.3f",hBkg->GetMean(),hBkg->GetMeanError()));
    latFD->AddText(Form("<#chi^{2}/ndf>_{feed-down}= %0.2f #pm %0.2f",hFD->GetMean(),hFD->GetMeanError()));
  }
  else {
    latBkg->AddText(Form("<#chi^{2}/ndf>_{bkg}= %0.2f #pm %0.2f",hBkg->GetMean(),hBkg->GetMeanError()));
    latFD->AddText(Form("<#chi^{2}/ndf>_{feed-down}= %0.1f #pm %0.1f",hFD->GetMean(),hFD->GetMeanError()));
  }
  
  TCanvas *cout = new TCanvas("cout","",800,800);
  cout->SetLogy();
  hPrompt->Draw();
  hFD->Draw("same");
  hBkg->Draw("same");
  l->Draw("same");
  if(iProj==12 || iProj==13) {
    latPrompt->Draw("same");
    latFD->Draw("same");
    latBkg->Draw("same");
  }
    
  TString outname = filename.ReplaceAll(".root",Form("_pT%d.eps",iPt));
  cout->SaveAs(outname.Data());
}
