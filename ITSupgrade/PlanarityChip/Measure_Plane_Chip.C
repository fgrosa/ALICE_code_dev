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
#include <TPaveText.h>

#endif

//****************************************************************************//
//                                                                            //
//    Main Functions: Measure_Plane_Chip(), GetRelativePlanesInclination()    //
//                                                                            //
//****************************************************************************//

//_____________________________________________________________________________________________
//GLOBAL VARIBALES
const Double_t ymin=-10;
const Double_t ymax=20;
const Double_t xmin=-15;
const Double_t xmax=5;

//_____________________________________________________________________________________________
//FUNCTION PROTOTYPES
Int_t Measure_Plane_Chip(TString FileName="misure_handlingbar_plane_6_10_2016.txt",
                                TString OutFileNamewoext="misure_handlingbar_plane_6_10_2016_corr");
Int_t GetRelativePlanesInclination(TString FileName1="misure_handlingbar_plane_6_10_2016_corr.root", TString FileName2="misure_handlingbar_plane_7_10_2016_corr.root", TString OutFileNamewoext="comp_handlingbar_plane_7_10_2016_6_10_2016");

Int_t ReadFile(TString FileName, vector<double> &x, vector<double> &y, vector<double> &z);
void CreateTxtFile(TString FileName, vector<double> &x, vector<double> &y, vector<double> &z);
void GetMeanSigmaAndPlanarity(vector<double> z, Double_t &mean, Double_t &sigma, Double_t &planarity);
void GetDirCosine(TF2* fplane, Double_t dircos[3]);

//_____________________________________________________________________________________________
//FUNCTION FOR PLANARITY MEASUREMENT
Int_t Measure_Plane_Chip(TString FileName, TString OutFileNamewoext) {

  vector<double> x;
  vector<double> y;
  vector<double> z;

  Int_t read=ReadFile(FileName,x,y,z);
  if(read<0) {
    cerr << "Impossibile to find the input file. Exit" << endl;
    return 1;
  }
  if(x.size()<=0 || y.size()<=0 || z.size()<=0) {
    cerr << "The number of points is <=0. Exit" << endl;
    return 2;
  }
  if(y.size()!=z.size() || x.size()!=z.size()) {
    cerr << "There is a discrepacy between the number of entries for the different coordinates. Exit." << endl;
    return 3;
  }

  UInt_t nPoints=z.size();
  
  TGraph2D *g = new TGraph2D(nPoints);
  for(UInt_t iEntry=0; iEntry<nPoints; iEntry++) {
    g->SetPoint(iEntry,x[iEntry],y[iEntry],z[iEntry]);
  }

  TF2* fplane = new TF2("fplane","[0]+[1]*x+[2]*y",xmin,xmax,ymin,ymax);

  TCanvas* c = new TCanvas("c","",1200,900);
  g->SetMarkerStyle(20);
  fplane->GetXaxis()->SetTitle("x (mm)");
  fplane->GetYaxis()->SetTitle("y (mm)");
  fplane->GetZaxis()->SetTitle("z (mm)");
  g->SetTitle("");
  g->Fit("fplane");
  fplane->SetTitle("");
  fplane->GetZaxis()->SetRange(-16,-15);
  fplane->Draw("surf");
  g->Draw("P0 same");
  
  Double_t para = fplane->GetParameter(0);
  Double_t parb = fplane->GetParameter(1);
  Double_t parc = fplane->GetParameter(2);

  TGraph2D* gcorr = new TGraph2D(nPoints);
  gcorr->SetName("gcorr");
  TGraph* gcorr_dx = new TGraph(nPoints/3);
  gcorr->SetName("gcorr_dx");
  TGraph* gcorr_cent = new TGraph(nPoints/3);
  gcorr->SetName("gcorr_cent");
  TGraph* gcorr_sx = new TGraph(nPoints/3);
  gcorr->SetName("gcorr_sx");
  
  vector<double> zcorr;
  UInt_t poscounter=0;
  UInt_t iPoint=0;
  for(UInt_t iEntry=0; iEntry<nPoints; iEntry++) {
    if(poscounter>2) poscounter=0;
    zcorr.push_back(z[iEntry]-(para+parb*x[iEntry]+parc*y[iEntry]));

    gcorr->SetPoint(iEntry,x[iEntry],y[iEntry],zcorr[iEntry]);

    if(poscounter==0)
      gcorr_dx->SetPoint(iPoint,y[iEntry],zcorr[iEntry]);
    else if(poscounter==1)
      gcorr_cent->SetPoint(iPoint,y[iEntry],zcorr[iEntry]);
    else {
      gcorr_sx->SetPoint(iPoint,y[iEntry],zcorr[iEntry]);
      iPoint++;
    }
    poscounter++;
  }

  Double_t mean=0;
  Double_t sigma=0;
  Double_t planarity=0;
  GetMeanSigmaAndPlanarity(zcorr,mean,sigma,planarity);
  
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
  l->AddEntry(gcorr_dx,Form("x = %0.1f mm",x[0]),"p");
  l->AddEntry(gcorr_cent,Form("x = %0.1f mm",x[1]),"p");
  l->AddEntry(gcorr_sx,Form("x = %0.1f mm",x[2]),"p");
  
  TPaveText* info = new TPaveText(0.6,0.75,0.89,0.89,"NDC");
  info->SetTextSize(0.045);
  info->SetTextFont(42);
  info->SetFillColor(0);
  info->SetFillStyle(0);
  info->SetBorderSize(0);
  info->AddText(Form("mean = %0.3f mm",mean));
  info->AddText(Form("RMS = %0.3f mm",sigma));
  info->AddText(Form("planarity = %0.3f mm",planarity));
  
  TCanvas *ccorr = new TCanvas("ccorr","",1200,900);
  gcorr_dx->SetTitle("");
  gcorr_dx->GetXaxis()->SetTitle("y (mm)");
  gcorr_dx->GetYaxis()->SetTitle("z (mm)");
  gcorr_dx->GetYaxis()->SetRangeUser(-0.300,0.350);
  gcorr_dx->SetMarkerStyle(20);
  gcorr_cent->SetMarkerStyle(20);
  gcorr_sx->SetMarkerStyle(20);
  gcorr_dx->SetMarkerColor(kBlack);
  gcorr_cent->SetMarkerColor(kBlue);
  gcorr_sx->SetMarkerColor(kGreen+2);
  gcorr_dx->SetTitle("");
  gcorr_dx->Draw("AP");
  gcorr_cent->Draw("P");
  gcorr_sx->Draw("P");
  l->Draw("same");
  info->Draw("same");

  ccorr2D->SaveAs(Form("%s.pdf",OutFileNamewoext.Data()));
  ccorr->SaveAs(Form("%s.pdf",OutFileNamewoext.Data()));
  TString outfiletxt = Form("%s.txt",OutFileNamewoext.Data());
  CreateTxtFile(outfiletxt,x,y,zcorr);
  
  TFile outfile(Form("%s.root",OutFileNamewoext.Data()),"RECREATE");
  gcorr->Write();
  gcorr_sx->Write();
  gcorr_cent->Write();
  gcorr_dx->Write();
  fplane->Write();
  outfile.Close();
  
  return 0;
}

Int_t GetRelativePlanesInclination(TString FileName1, TString FileName2, TString OutFileNamewoext) {
  
  TFile* infile1 = TFile::Open(FileName1.Data(),"READ");
  TF2* fplane1=0x0;
  if(infile1) {fplane1=(TF2*)infile1->Get("fplane");}
  else {cerr << "File " << FileName1 << " does not exists. Exit" << endl; return 1;}
  infile1->Close();
  
  TFile* infile2 = TFile::Open(FileName2.Data(),"READ");
  TF2* fplane2=0x0;
  if(infile2) {fplane2=(TF2*)infile2->Get("fplane");fplane2->SetLineColor(kBlue);}
  else {cerr << "File " << FileName1 << " does not exists. Exit" << endl; return 2;}
  infile2->Close();
  
  TCanvas* cPlanes = new TCanvas("cPlanes","",1200,900);
  fplane1->Draw("surf");
  fplane2->Draw("surf same");
  
  Double_t* pars1 = (Double_t*)fplane1->GetParameters();
  Double_t* pars2 = (Double_t*)fplane2->GetParameters();
  
  Double_t xref = (xmax+xmin)/2;
  Double_t yref = (ymax+ymin)/2;
  
  TF1* flinezx1 = new TF1("flinezx1","pol1",0,xmax-xmin);
  flinezx1->SetParameters(0,pars1[1]);
  flinezx1->GetYaxis()->SetTitle("z(mm)");
  flinezx1->GetXaxis()->SetTitle("x(mm)");
  TF1* flinezx2 = new TF1("flinezx2","pol1",0,xmax-xmin);
  flinezx2->SetParameters(0,pars2[1]);
  flinezx2->GetYaxis()->SetTitle("z(mm)");
  flinezx2->GetXaxis()->SetTitle("x(mm)");
  flinezx2->SetLineColor(kBlue);
  
  Double_t toll=0.1;
  Double_t minzx=flinezx1->Eval(xmax-xmin);
  if(flinezx2->Eval(xmax-xmin)<flinezx1->Eval(xmax-xmin)) minzx=flinezx2->Eval(xmax-xmin);
  Double_t maxzx=flinezx1->Eval(xmax-xmin);
  if(flinezx2->Eval(xmax-xmin)>flinezx1->Eval(xmax-xmin)) maxzx=flinezx2->Eval(xmax-xmin);
  
  if(minzx>0) minzx=-0.;
  if(maxzx<0) maxzx=0.;
  
  Double_t anglezx1 = TMath::ATan(flinezx1->GetParameter(1));
  Double_t anglezx2 = TMath::ATan(flinezx2->GetParameter(1));
  Double_t deltaanglezx = anglezx1-anglezx2;
  if(deltaanglezx<0) deltaanglezx = TMath::Abs(deltaanglezx);
  
  TF1* flinezy1 = new TF1("flinezy1","pol1",0,ymax-ymin);
  flinezy1->SetParameters(0,pars1[2]);
  flinezy1->GetYaxis()->SetTitle("z(mm)");
  flinezy1->GetXaxis()->SetTitle("y(mm)");
  TF1* flinezy2 = new TF1("flinezy2","pol1",0,ymax-ymin);
  flinezy2->SetParameters(0,pars2[2]);
  flinezy2->GetYaxis()->SetTitle("z(mm)");
  flinezy2->GetXaxis()->SetTitle("y(mm)");
  flinezy2->SetLineColor(kBlue);
  
  Double_t minzy=flinezy1->Eval(ymax-ymin);
  if(flinezy2->Eval(ymax-ymin)<flinezy1->Eval(ymax-ymin)) minzy=flinezy2->Eval(ymax-ymin);
  Double_t maxzy=flinezy1->Eval(ymax-ymin);
  if(flinezy2->Eval(ymax-ymin)>flinezy1->Eval(ymax-ymin)) maxzy=flinezy2->Eval(ymax-ymin);
  
  if(minzy>0) minzy=0.;
  if(maxzy<0) maxzy=0.;
  
  Double_t anglezy1 = TMath::ATan(flinezy1->GetParameter(1));
  Double_t anglezy2 = TMath::ATan(flinezy2->GetParameter(1));
  Double_t deltaanglezy = anglezy1-anglezy2;
  if(deltaanglezy<0) deltaanglezy = TMath::Abs(deltaanglezy);
  
  TCanvas* cZX = new TCanvas("cZX","",1200,900);
  flinezx1->SetTitle(Form("y = %0.2f mm",yref));
  flinezx1->GetYaxis()->SetTitleOffset(1.3);
  flinezx1->GetYaxis()->SetRangeUser(minzx-toll,maxzx+toll);
  TPaveText* infozx_1 = new TPaveText(0.38,0.78,0.89,0.85,"NDC");
  infozx_1->SetTextSize(0.045);
  infozx_1->SetTextFont(42);
  infozx_1->SetTextColor(kRed);
  infozx_1->SetFillColor(0);
  infozx_1->SetFillStyle(0);
  infozx_1->SetBorderSize(0);
  infozx_1->AddText(Form("#alpha_{plane 1} = %0.4f rad ( = %0.2f deg)",anglezx1,anglezx1*180/TMath::Pi()));
  TPaveText* infozx_2 = new TPaveText(0.38,0.7,0.89,0.78,"NDC");
  infozx_2->SetTextSize(0.045);
  infozx_2->SetTextFont(42);
  infozx_2->SetTextColor(kBlue);
  infozx_2->SetFillColor(0);
  infozx_2->SetFillStyle(0);
  infozx_2->SetBorderSize(0);
  infozx_2->AddText(Form("#alpha_{plane 2} = %0.4f rad ( = %0.2f deg)",anglezx2,anglezx2*180/TMath::Pi()));
  TPaveText* infozx_3 = new TPaveText(0.38,0.62,0.89,0.7,"NDC");
  infozx_3->SetTextSize(0.045);
  infozx_3->SetTextFont(42);
  infozx_3->SetFillColor(0);
  infozx_3->SetFillStyle(0);
  infozx_3->SetBorderSize(0);
  infozx_3->AddText(Form("#Delta#alpha = |#alpha_{plane 2}-#alpha_{plane 1}| = %0.4f rad ( = %0.2f deg)",deltaanglezx,deltaanglezx*180/TMath::Pi()));
  flinezx1->Draw();
  flinezx2->Draw("same");
  infozx_1->Draw("same");
  infozx_2->Draw("same");
  infozx_3->Draw("same");
  
  TCanvas* cZY = new TCanvas("cZY","",1200,900);
  flinezy1->SetTitle(Form("x = %0.2f mm",xref));
  flinezy1->GetYaxis()->SetTitleOffset(1.2);
  flinezy1->GetYaxis()->SetRangeUser(minzy-toll,maxzy+toll);
  TPaveText* infozy_1 = new TPaveText(0.38,0.78,0.89,0.85,"NDC");
  infozy_1->SetTextSize(0.045);
  infozy_1->SetTextFont(42);
  infozy_1->SetTextColor(kRed);
  infozy_1->SetFillColor(0);
  infozy_1->SetFillStyle(0);
  infozy_1->SetBorderSize(0);
  infozy_1->AddText(Form("#alpha_{plane 1} = %0.4f rad ( = %0.2f deg)",anglezy1,anglezy1*180/TMath::Pi()));
  TPaveText* infozy_2 = new TPaveText(0.38,0.7,0.89,0.78,"NDC");
  infozy_2->SetTextSize(0.045);
  infozy_2->SetTextFont(42);
  infozy_2->SetTextColor(kBlue);
  infozy_2->SetFillColor(0);
  infozy_2->SetFillStyle(0);
  infozy_2->SetBorderSize(0);
  infozy_2->AddText(Form("#alpha_{plane 2} = %0.4f rad ( = %0.2f deg)",anglezy2,anglezy2*180/TMath::Pi()));
  TPaveText* infozy_3 = new TPaveText(0.38,0.62,0.89,0.7,"NDC");
  infozy_3->SetTextSize(0.045);
  infozy_3->SetTextFont(42);
  infozy_3->SetFillColor(0);
  infozy_3->SetFillStyle(0);
  infozy_3->SetBorderSize(0);
  infozy_3->AddText(Form("#Delta#alpha = |#alpha_{plane 2}-#alpha_{plane 1}| = %0.4f rad ( = %0.2f deg)",deltaanglezy,deltaanglezy*180/TMath::Pi()));
  flinezy1->Draw();
  flinezy2->Draw("same");
  infozy_1->Draw("same");
  infozy_2->Draw("same");
  infozy_3->Draw("same");
  
  cZX->SaveAs(Form("%s_ZX.pdf",OutFileNamewoext.Data()));
  cZY->SaveAs(Form("%s_ZY.pdf",OutFileNamewoext.Data()));
  
  return 0;
}

Int_t ReadFile(TString FileName, vector<double> &x, vector<double> &y, vector<double> &z) {
  ifstream inSet(FileName.Data());
  if(!inSet){
    cerr<<"File "<<FileName.Data() <<" does not exists. "<<endl;
    return -1;
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
  
  return 0;
}

void CreateTxtFile(TString FileName, vector<double> &x, vector<double> &y, vector<double> &z) {
  ofstream outSet;
  outSet.open(FileName.Data());

  for(UInt_t iEntry=0; iEntry<z.size(); iEntry++) {
      outSet << x[iEntry] << " " << y[iEntry] << " " << z[iEntry] << endl;
  }

  outSet.close();
}

void GetMeanSigmaAndPlanarity(vector<double> z, Double_t &mean, Double_t &sigma, Double_t &planarity) {
  
  mean=0;
  sigma=0;
  Double_t max=z[0];
  Double_t min=z[0];
  
  for(UInt_t iEntry=0; iEntry<z.size(); iEntry++) {
    mean += z[iEntry];
    if(z[iEntry]>max)
      max=z[iEntry];
    if(z[iEntry]<min)
      min=z[iEntry];
  }
  mean /= z.size();
  
  for(UInt_t iEntry=0; iEntry<z.size(); iEntry++) {
    sigma += (z[iEntry]-mean)*(z[iEntry]-mean);
  }
  sigma /= (z.size()-1);
  sigma = TMath::Sqrt(sigma);
  
  planarity = max-min;
  
}

void GetDirCosine(TF2* fplane, Double_t dircos[3]) {
  
  Double_t a = fplane->GetParameter(0);
  Double_t b = fplane->GetParameter(1);
  Double_t c = fplane->GetParameter(2);
  
  Double_t norm = TMath::Sqrt(a*a+b*b+c*c);
  dircos[0] = a/norm;
  dircos[1] = b/norm;
  dircos[2] = c/norm;
  
}

