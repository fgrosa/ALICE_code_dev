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
const Int_t nFiles=2;
const TString filenames[nFiles] = {"$HOME/ALICE_WORK/AnalysisPbPb2015/v2/v2Output_30_50_anis_TPC.root"
                                  ,"$HOME/ALICE_WORK/PublishedResults/PbPb_2.76TeV_2011/AverDv2-HalfTPCEP.root"};
const TString legendnames[nFiles] = {"D^{+}, #sqrt{s_{NN}} = 5.02 TeV ","D-meson average, #sqrt{s_{NN}} = 2.76 TeV"};
const TString graphstatnames[nFiles] = {"gav2fs","gv2aveep3050"};
const TString graphsystnames[nFiles] = {"",""};
const Int_t colors[] = {kRed,kBlack};
const Int_t markers[] = {kFullSquare,kFullCircle};

TString outfilename = "DmesonFlow_Comparison_2.76-5TeV.pdf";

//_________________________________________________________________________________________________
Int_t FlowComparison() {
  
  TGraphAsymmErrors** gstat = new TGraphAsymmErrors*[nFiles];
  TGraphAsymmErrors** gsyst = new TGraphAsymmErrors*[nFiles];
  
  for(Int_t iFile=0; iFile<nFiles; iFile++) {
    TFile* infile = TFile::Open(filenames[iFile].Data(),"READ");
    if(graphstatnames[iFile]!="") gstat[iFile] = (TGraphAsymmErrors*)infile->Get(graphstatnames[iFile]);
    else gstat[iFile] = 0x0;
    if(graphsystnames[iFile]!="") gsyst[iFile] = (TGraphAsymmErrors*)infile->Get(graphsystnames[iFile]);
    else gsyst[iFile] = 0x0;
    infile->Close();
  }
  
  TLegend *leg = new TLegend(0.35,0.7,0.89,0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.045);
  for(Int_t iFile=0; iFile<nFiles; iFile++) {
    if(gstat[iFile]) leg->AddEntry(gstat[iFile],legendnames[iFile],"lpe");
  }
  
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetPadBottomMargin(0.14);
  
  TCanvas* cv2 = new TCanvas("cv2","",1920,1080);
  cv2->Clear();
  TString drawopt="AP";
  for(Int_t iFile=0; iFile<nFiles; iFile++) {
    if(gstat[iFile]) {
      gstat[iFile]->SetTitle("");
      gstat[iFile]->GetYaxis()->SetRangeUser(-0.1,0.6);
      gstat[iFile]->GetYaxis()->SetTitle("v_{2}");
      gstat[iFile]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
      gstat[iFile]->SetLineColor(colors[iFile]);
      gstat[iFile]->SetMarkerColor(colors[iFile]);
      gstat[iFile]->SetMarkerStyle(markers[iFile]);
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
  
  if(outfilename.Contains("root")) {
    TFile outfile(outfilename.Data(),"RECREATE");
    for(Int_t iFile=0; iFile<nFiles; iFile++) {
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
