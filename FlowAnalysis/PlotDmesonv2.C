#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TString.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TASImage.h>
#include <TPad.h>
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>

#endif

//_________________________________________________________________________________________________
//global variables
const TString infilename = "$HOME/ALICE_WORK/AnalysisPbPb2015/v2/v2Output_30_50_anis_Topod0Cut_VZERO_EP.root";
const TString graphstatname = "gav2fs";
const Double_t yieldsystunc[] = {0.03,0.03,0.02,0.03,0.02,0.03,0.03,0.02,0.02,0.04,0.03};

void PlotDmesonv2() {
  
  TGraphAsymmErrors* gstat = 0x0;
  TGraphAsymmErrors* gsyst = 0x0;

  TFile* infile = TFile::Open(infilename.Data(),"READ");
  if(infile) {
    gstat = (TGraphAsymmErrors*)infile->Get(graphstatname);
    infile->Close();
    gsyst = (TGraphAsymmErrors*)gstat->Clone();
    for(Int_t iPt=0; iPt<gstat->GetN(); iPt++) {
      gsyst->SetPointEYhigh(iPt,yieldsystunc[iPt]);
      gsyst->SetPointEYlow(iPt,yieldsystunc[iPt]);
      gsyst->SetPointEXhigh(iPt,0.25);
      gsyst->SetPointEXlow(iPt,0.25);
    }
  }
  
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.14);
  TCanvas* cv2 = new TCanvas("cv2","",800,800);
  TLine* line = new TLine(0.,0.,26.3,0.);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  gstat->SetMarkerSize(1);
  gstat->GetYaxis()->SetTitleOffset(1.4);
  gstat->GetYaxis()->SetRangeUser(-0.1,0.35);
  gsyst->SetFillStyle(20);
  gsyst->SetLineWidth(2);
  gstat->Draw("AP");
  gsyst->Draw("2");
  line->Draw("same");
  
  cv2->SaveAs("v2_result_3050.pdf");
}
