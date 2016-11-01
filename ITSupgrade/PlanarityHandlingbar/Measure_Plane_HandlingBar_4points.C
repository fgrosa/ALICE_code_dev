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
#include <TF2.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <THnSparse.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TPad.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TLine.h>
#include <TList.h>
#include <TDirectory.h>
#include <TBox.h>
#include <TVirtualFitter.h>

#endif

Int_t Measure_Plane_HandlingBar(TString FileName="misure_piano_barra_piedi_7_10_2016.txt",
                               TString OutFileNamewoext="misure_piano_barra_piedi_7_10_2016_corr",
                                Double_t x1=80, Double_t x2=120, Double_t x3=74, Double_t x4=125, Double_t toll=0.02);
void ReadFile(TString FileName, vector<double> &x, vector<double> &y, vector<double> &z);
void CreateTxtFile(TString FileName, vector<double> &x, vector<double> &y, vector<double> &z);

Int_t Measure_Plane_HandlingBar(TString FileName, TString OutFileNamewoext, Double_t x1, Double_t x2,  Double_t x3, Double_t x4, Double_t toll) {

  vector<double> x;
  vector<double> y;
  vector<double> z;

  ReadFile(FileName,x,y,z);
  if(x.size()<=0 || y.size()<=0 || z.size()<=0) {
    cerr << "The number of points is <=0. Exit" << endl;
    return 1;
  }
  if(y.size()!=z.size() || x.size()!=z.size()) {
    cerr << "There is a discrepacy between the number of entries for the different coordinates. Exit." << endl;
    return 2;
  }

  TGraph2D *g = new TGraph2D(y.size());
  for(UInt_t iEntry=0; iEntry<y.size(); iEntry++) {
    g->SetPoint(iEntry,x[iEntry],y[iEntry],z[iEntry]);
  }

  TF2* fplane = new TF2("fplane","[0]+[1]*x+[2]*y",85,115,-1000,1000);

  TCanvas* c = new TCanvas("c","",1200,900);
  g->SetMarkerStyle(20);
  fplane->GetXaxis()->SetTitle("x (mm)");
  fplane->GetYaxis()->SetTitle("y (mm)");
  fplane->GetZaxis()->SetTitle("z (mm)");
  g->SetTitle("");
  g->Fit("fplane");
  fplane->SetTitle("");
  fplane->Draw("surf");
  g->Draw("P0 same");
  
  
  Double_t para = fplane->GetParameter(0);
  Double_t parb = fplane->GetParameter(1);
  Double_t parc = fplane->GetParameter(2);

  Int_t npoints1=0;
  Int_t npoints2=0;
  Int_t npoints3=0;
  Int_t npoints4=0;
  for(UInt_t iEntry=0; iEntry<y.size(); iEntry++) {
    if(x[iEntry]>83.27 && x[iEntry]<87.8)
      npoints1++;
    else if(x[iEntry]>123.0 && x[iEntry]<127.2 && (iEntry<10 || iEntry>16))
      npoints2++;
    else if(x[iEntry]>82.0 && x[iEntry]<83.10)
      npoints3++;
    else if(x[iEntry]>126.9 && x[iEntry]<127.1 && (iEntry>10 && iEntry<16))
      npoints4++;
    else
      cerr << "x coordinate out of range. Change tollerance." << endl;
  }
  
  vector<double> zcorr;
  TGraph2D* gcorr = new TGraph2D(z.size());
  TGraph* gcorr_1 = new TGraph(npoints1);
  TGraph* gcorr_2 = new TGraph(npoints2);
  TGraph* gcorr_3 = new TGraph(npoints3);
  TGraph* gcorr_4 = new TGraph(npoints4);
  for(UInt_t iEntry=0; iEntry<y.size(); iEntry++) {
    zcorr.push_back(z[iEntry]-(para+parb*x[iEntry]+parc*y[iEntry]));
    
    gcorr->SetPoint(iEntry,x[iEntry],y[iEntry],zcorr[iEntry]);
    
    if(x[iEntry]>83.27 && x[iEntry]<87.8)
      gcorr_1->SetPoint(iEntry,y[iEntry],zcorr[iEntry]);
    else if(x[iEntry]>123.0 && x[iEntry]<127.2 && (iEntry<10 || iEntry>16))
      gcorr_2->SetPoint(iEntry,y[iEntry],zcorr[iEntry]);
    else if(x[iEntry]>82.0 && x[iEntry]<83.10)
      gcorr_3->SetPoint(iEntry,y[iEntry],zcorr[iEntry]);
    else if(x[iEntry]>126.9 && x[iEntry]<127.1 && (iEntry>10 && iEntry<16))
      gcorr_4->SetPoint(iEntry,y[iEntry],zcorr[iEntry]);
    else
      cerr << "x coordinate out of range. Change tollerance." << endl;
  }

  TBox* box1 = new TBox(-745.07,-0.200,-744.07,0.200);
  box1->SetLineColor(kRed);
  box1->SetLineWidth(1);
  box1->SetFillStyle(0);
  TBox* box2 = new TBox(-685.07,-0.200,-655.07,0.200);
  box2->SetLineColor(kRed);
  box2->SetLineWidth(1);
  box2->SetFillStyle(0);
  TBox* box3 = new TBox(-585.07,-0.200,-555.07,0.200);
  box3->SetLineColor(kRed);
  box3->SetLineWidth(1);
  box3->SetFillStyle(0);
  TBox* box4 = new TBox(-435.07,-0.200,-395.07,0.200);
  box4->SetLineColor(kRed);
  box4->SetLineWidth(1);
  box4->SetFillStyle(0);
  TBox* box5 = new TBox(-235.07,-0.200,-205.07,0.200);
  box5->SetLineColor(kRed);
  box5->SetLineWidth(1);
  box5->SetFillStyle(0);
  TBox* box6 = new TBox(-75.07,-0.200,-35.07,0.200);
  box6->SetLineColor(kRed);
  box6->SetLineWidth(1);
  box6->SetFillStyle(0);
  TBox* box7 = new TBox(34.93,-0.200,74.93,0.200);
  box7->SetLineColor(kRed);
  box7->SetLineWidth(1);
  box7->SetFillStyle(0);
  TBox* box8 = new TBox(194.93,-0.200,224.93,0.200);
  box8->SetLineColor(kRed);
  box8->SetLineWidth(1);
  box8->SetFillStyle(0);
  TBox* box9 = new TBox(394.93,-0.200,434.93,0.200);
  box9->SetLineColor(kRed);
  box9->SetLineWidth(1);
  box9->SetFillStyle(0);
  TBox* box10 = new TBox(534.93,-0.200,584.93,0.200);
  box10->SetLineColor(kRed);
  box10->SetLineWidth(1);
  box10->SetFillStyle(0);
  TBox* box11 = new TBox(644.93,-0.200,684.93,0.200);
  box11->SetLineColor(kRed);
  box11->SetLineWidth(1);
  box11->SetFillStyle(0);
  TBox* box12 = new TBox(744.93,-0.200,745.93,0.200);
  box12->SetLineColor(kRed);
  box12->SetLineWidth(1);
  box12->SetFillStyle(0);

  TCanvas *ccorr2D = new TCanvas("ccorr2D","",1200,900);
  gcorr->SetTitle("");
  gcorr->GetXaxis()->SetTitle("x (mm)");
  gcorr->GetYaxis()->SetTitle("y (mm)");
  gcorr->GetZaxis()->SetTitle("y (mm)");
  gcorr->GetZaxis()->SetRangeUser(-0.250,0.250);
  gcorr->SetMarkerStyle(20);
  gcorr->SetTitle("");
  gcorr->Draw("AP");
  
  TLegend* l = new TLegend(0.15,0.75,0.45,0.89);
  l->SetBorderSize(0);
  l->SetFillColor(kWhite);
  l->SetFillStyle(0);
  l->SetTextSize(0.045);
  l->AddEntry(gcorr_1,Form("x = %0.f #mum",x1),"p");
  l->AddEntry(gcorr_2,Form("x = %0.f #mum",x2),"p");
  l->AddEntry(gcorr_3,Form("x = %0.f #mum",x3),"p");
  l->AddEntry(gcorr_4,Form("x = %0.f #mum",x4),"p");
  
  TCanvas *ccorr = new TCanvas("ccorr","",1200,900);
  gcorr_1->SetTitle("");
  gcorr_1->GetXaxis()->SetTitle("y (mm)");
  gcorr_1->GetYaxis()->SetTitle("z (mm)");
  gcorr_1->GetYaxis()->SetRangeUser(-0.300,0.350);
  gcorr_1->SetMarkerStyle(20);
  gcorr_2->SetMarkerStyle(20);
  gcorr_3->SetMarkerStyle(20);
  gcorr_4->SetMarkerStyle(20);
  gcorr_1->SetMarkerColor(kBlack);
  gcorr_2->SetMarkerColor(kBlue);
  gcorr_3->SetMarkerColor(kGreen+2);
  gcorr_4->SetMarkerColor(kOrange+7);
  gcorr_1->SetTitle("");
  gcorr_1->Draw("AP");
  gcorr_2->Draw("P");
  gcorr_3->Draw("P");
  gcorr_4->Draw("P");
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
  l->Draw("same");

  ccorr2D->SaveAs(Form("%s2D.pdf",OutFileNamewoext.Data()));
  ccorr->SaveAs(Form("%s.pdf",OutFileNamewoext.Data()));
  TString outfiletxt = Form("%s.txt",OutFileNamewoext.Data());
  CreateTxtFile(outfiletxt,x,y,zcorr);
  
  return 0;
}

void ReadFile(TString FileName, vector<double> &x, vector<double> &y, vector<double> &z) {
  ifstream inSet(FileName.Data());
  if(!inSet){
    cerr<<"File "<<FileName.Data() <<" does not exists. "<<endl;
    return;
  }

  Double_t xcoord;
  Double_t ycoord;
  Double_t zcoord;
  while(inSet>>xcoord>>ycoord>>zcoord) {
    x.push_back(xcoord);
    y.push_back(ycoord);
    z.push_back(zcoord);
  }

  inSet.close();
}

void CreateTxtFile(TString FileName, vector<double> &x, vector<double> &y, vector<double> &z) {
  ofstream outSet;
  outSet.open(FileName.Data());

  for(UInt_t iEntry=0; iEntry<z.size(); iEntry++) {
      outSet << x[iEntry] << " " << y[iEntry] << " " << z[iEntry] << endl;
  }

  outSet.close();
}
