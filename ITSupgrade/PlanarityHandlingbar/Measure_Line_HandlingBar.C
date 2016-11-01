#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <TInterpreter.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <THnSparse.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TPad.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TLine.h>
#include <TList.h>
#include <TDirectory.h>
#include <TBox.h>

#endif

Int_t Measure_Line_HandlingBar(TString FileName="misure_handlingbar_3_10_2016.txt",
                               TString OutFileNamewoext="misure_handlingbar_3_10_2016_corr");
void ReadFile(TString FileName, vector<double> &y, vector<double> &z);
void CreateTxtFile(TString FileName, vector<double> &y, vector<double> &z);

Int_t Measure_Line_HandlingBar(TString FileName, TString OutFileNamewoext) {

  vector<double> y;
  vector<double> z;

  ReadFile(FileName,y,z);
  if(y.size()<=0 || z.size()<=0) {
    cerr << "Il numero di punti è <=0" << endl;
    return 1;
  }
  if(y.size()!=z.size()) {
    cerr << "Il numero di coordinate y è diversa dal numero di coordinate z" << endl;
    return 2;
  }

  TGraph *g = new TGraph(y.size());
  for(UInt_t iEntry=0; iEntry<y.size(); iEntry++) {
    g->SetPoint(iEntry,y[iEntry],z[iEntry]);
  }

  TF1* fline = new TF1("fline","pol1",-1000,1000);

  TCanvas* c = new TCanvas("c","",1200,900);
  g->SetMarkerStyle(20);
  g->GetXaxis()->SetTitle("y (mm)");
  g->GetYaxis()->SetTitle("z (mm)");
  g->SetTitle("");
  g->Draw("AP");
  g->Fit("fline");

  Double_t constant = fline->GetParameter(0);
  Double_t slope = fline->GetParameter(1);

  vector<double> zcorr;
  TGraph* gcorr = new TGraph(z.size());
  for(UInt_t iEntry=0; iEntry<y.size(); iEntry++) {
    zcorr.push_back(z[iEntry]-(constant+slope*y[iEntry]));
    gcorr->SetPoint(iEntry,y[iEntry],zcorr[iEntry]);
  }


  TBox* box1 = new TBox(y[0]-1,-0.200,y[0]+1,0.200);
  box1->SetLineColor(kRed);
  box1->SetLineWidth(1);
  box1->SetFillStyle(0);
  TBox* box2 = new TBox(y[6],-0.200,y[9],0.200);
  box2->SetLineColor(kRed);
  box2->SetLineWidth(1);
  box2->SetFillStyle(0);
  TBox* box3 = new TBox(y[16],-0.200,y[19],0.200);
  box3->SetLineColor(kRed);
  box3->SetLineWidth(1);
  box3->SetFillStyle(0);
  TBox* box4 = new TBox(y[31],-0.200,y[35],0.200);
  box4->SetLineColor(kRed);
  box4->SetLineWidth(1);
  box4->SetFillStyle(0);
  TBox* box5 = new TBox(y[51],-0.200,y[54],0.200);
  box5->SetLineColor(kRed);
  box5->SetLineWidth(1);
  box5->SetFillStyle(0);
  TBox* box6 = new TBox(y[67],-0.200,y[71],0.200);
  box6->SetLineColor(kRed);
  box6->SetLineWidth(1);
  box6->SetFillStyle(0);
  TBox* box7 = new TBox(y[78],-0.200,y[82],0.200);
  box7->SetLineColor(kRed);
  box7->SetLineWidth(1);
  box7->SetFillStyle(0);
  TBox* box8 = new TBox(y[94],-0.200,y[97],0.200);
  box8->SetLineColor(kRed);
  box8->SetLineWidth(1);
  box8->SetFillStyle(0);
  TBox* box9 = new TBox(y[114],-0.200,y[118],0.200);
  box9->SetLineColor(kRed);
  box9->SetLineWidth(1);
  box9->SetFillStyle(0);
  TBox* box10 = new TBox(y[128],-0.200,y[133],0.200);
  box10->SetLineColor(kRed);
  box10->SetLineWidth(1);
  box10->SetFillStyle(0);
  TBox* box11 = new TBox(y[139],-0.200,y[143],0.200);
  box11->SetLineColor(kRed);
  box11->SetLineWidth(1);
  box11->SetFillStyle(0);
  TBox* box12 = new TBox(y[149]-1,-0.200,y[149]+1,0.200);
  box12->SetLineColor(kRed);
  box12->SetLineWidth(1);
  box12->SetFillStyle(0);

  TCanvas *ccorr = new TCanvas("ccorr","",1200,900);
  gcorr->SetTitle("");
  gcorr->GetXaxis()->SetTitle("y (mm)");
  gcorr->GetYaxis()->SetTitle("z (mm)");
  gcorr->GetYaxis()->SetRangeUser(-0.250,0.250);
  gcorr->SetMarkerStyle(20);
  g->SetTitle("");
  gcorr->Draw("AP");
  box1->Draw("same");
  box2->Draw("same");
  box3->Draw("same");
  box4->Draw("same");
  box5->Draw("same");
  box6->Draw("same");
  box7->Draw("same");
  box8->Draw("same");
  box9->Draw("same");
  box10->Draw("same");
  box11->Draw("same");
  box12->Draw("same");

  ccorr->SaveAs(Form("%s.pdf",OutFileNamewoext.Data()));
  TString outfiletxt = Form("%s.txt",OutFileNamewoext.Data());
  CreateTxtFile(outfiletxt,y,zcorr);

  return 0;
}

void ReadFile(TString FileName, vector<double> &y, vector<double> &z) {
  ifstream inSet(FileName.Data());
  if(!inSet){
    cerr<<"Il file "<<FileName.Data() <<" non esiste "<<endl;
    return;
  }

  Double_t ycoord;
  Double_t zcoord;
  while(inSet>>ycoord>>zcoord) {
    y.push_back(ycoord);
    z.push_back(zcoord);
    cout << ycoord << " " << zcoord << endl;
  }

  inSet.close();
}

void CreateTxtFile(TString FileName, vector<double> &y, vector<double> &z) {
  ofstream outSet;
  outSet.open(FileName.Data());

  for(UInt_t iEntry=0; iEntry<z.size(); iEntry++) {
      outSet << y[iEntry] << " " << z[iEntry] << endl;
  }

  outSet.close();
}
