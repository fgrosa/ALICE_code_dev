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
const Int_t nFiles=4;
const TString filenames[nFiles] = {"$HOME/ALICE_WORK/Files/PublishedResults/PbPb_2.76TeV_2011/AverDv2-HalfTPCEP.root","$HOME/ALICE_WORK/AnalysisPbPb2015/v2/v2Output_30_50_anis__3050_step2_QoverM_VZEROEPVZERO_EP_ptbinning3.root","$HOME/ALICE_WORK/AnalysisPbPb2015/v2/v2Output_30_50_anis_Topod0Cut_VZERO_EP.root", "$HOME/ALICE_WORK/AnalysisPbPb2015/v2/v2Output_30_50_anis_VZERO_EP_Ds_Step2.root"};
const TString legendnames[nFiles] = {"D-meson average, #sqrt{s_{NN}} = 2.76 TeV","D^{0}, #sqrt{s_{NN}} = 5.02 TeV ","D^{+}, #sqrt{s_{NN}} = 5.02 TeV ","D_{s}, #sqrt{s_{NN}} = 5.02 TeV "};
const TString graphstatnames[nFiles] = {"gv2aveep3050","gav2fs","gav2fs","gav2fs"};
const TString graphsystnames[nFiles] = {"","","",""};
const Int_t colors[] = {kBlack,kRed,kBlue,kGreen+2};
const Int_t markers[] = {kFullCircle,kFullDiamond,kFullSquare,kFullTriangleUp};

TString outfilename = "Dmeson_FlowComparison_2.76-5TeV.pdf";

const Int_t startfile=0;
const Int_t stopfile=2;

//_________________________________________________________________________________________________
Int_t FlowComparison() {
  
  TGraphAsymmErrors** gstat = new TGraphAsymmErrors*[nFiles];
  TGraphAsymmErrors** gsyst = new TGraphAsymmErrors*[nFiles];
  
  for(Int_t iFile=startfile; iFile<nFiles; iFile++) {
    TFile* infile = TFile::Open(filenames[iFile].Data(),"READ");
    if(infile && graphstatnames[iFile]!="") gstat[iFile] = (TGraphAsymmErrors*)infile->Get(graphstatnames[iFile]);
    else gstat[iFile] = 0x0;
    if(infile && graphsystnames[iFile]!="") gsyst[iFile] = (TGraphAsymmErrors*)infile->Get(graphsystnames[iFile]);
    else gsyst[iFile] = 0x0;
    if(infile) infile->Close();
  }
  
  TLegend *leg = new TLegend(0.15,0.7,0.45,0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.04);
  for(Int_t iFile=startfile; iFile<stopfile; iFile++) {
    if(gstat[iFile]) leg->AddEntry(gstat[iFile],legendnames[iFile],"lpe");
  }
  
  TLine* zeroline = new TLine(0.,0.,50.,0.);
  zeroline->SetLineWidth(2);
  zeroline->SetLineStyle(7);
  
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetPadBottomMargin(0.14);
  
  TH1F* hDummy = new TH1F("hDummy","",25,0.,50.);
  
  TCanvas* cv2 = new TCanvas("cv2","",1000,1000);
  cv2->Clear();
  hDummy->GetYaxis()->SetRangeUser(-0.4,1.);
  hDummy->SetStats(kFALSE);
  hDummy->SetLineColor(kWhite);
  if(gstat[1]) hDummy->GetYaxis()->SetTitle(gstat[2]->GetYaxis()->GetTitle());
  if(gstat[1]) hDummy->GetXaxis()->SetTitle(gstat[2]->GetXaxis()->GetTitle());
  hDummy->Draw();
  TString drawopt="P";
  for(Int_t iFile=startfile; iFile<stopfile; iFile++) {
    if(gstat[iFile]) {
      gstat[iFile]->SetTitle("");
      gstat[iFile]->GetYaxis()->SetTitle("v_{2}");
      gstat[iFile]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
      gstat[iFile]->SetLineColor(colors[iFile]);
      gstat[iFile]->SetMarkerColor(colors[iFile]);
      gstat[iFile]->SetMarkerStyle(markers[iFile]);
      gstat[iFile]->SetLineWidth(2);
      gstat[iFile]->SetMarkerSize(1.5);
      gstat[iFile]->Draw(drawopt);
      drawopt = "P";
    }
    if(gsyst[iFile]) {
      gsyst[iFile]->SetLineColor(colors[iFile]);
      gsyst[iFile]->SetMarkerColor(colors[iFile]);
      gsyst[iFile]->Draw(drawopt);
      drawopt = "P";
    }
  }
  leg->Draw("same");
  zeroline->Draw("same");

  if(outfilename.Contains("root")) {
    TFile outfile(outfilename.Data(),"RECREATE");
    for(Int_t iFile=startfile; iFile<stopfile; iFile++) {
      if(gstat[iFile]) gstat[iFile]->Write();
      if(gsyst[iFile]) gsyst[iFile]->Write();
    }
    cv2->Write();
    outfile.Close();
    outfilename.ReplaceAll("root","pdf");
    cv2->SaveAs(outfilename.Data());
  }
  else {
    cv2->SaveAs(outfilename.Data());
  }
  
  return 0;
}
