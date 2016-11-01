/////////////////////////////////////////////////////////////////
// Macro to Merge Stage by Stage                               //
// Need mergeFiles.dat with link to files.                     //
// in alien:                                                   //
// find SourceDIR  *nameFILE  > dataset.dat                    //
/////////////////////////////////////////////////////////////////

void MergeSbyS(Int_t start=0, Int_t pack=15, Int_t offset=0, char* fl ="AnalysisResultsFiles.dat" ){
  TGrid::Connect("alien:",0,0,"t");
  Int_t iMerge=0;
  Int_t imf=0;
  FILE* lFiles=fopen(fl,"r");
  Char_t fName[200];
  TFileMerger* mg= new TFileMerger();
  while(!feof(lFiles)){
    iMerge++;
    fscanf(lFiles,"%s\n",fName);
    if (fName[0]=='#') continue;
    printf("\n %d File name %s \n", iMerge, fName);
    if (iMerge > start){ 
      mg->AddFile(Form("alien://%s",fName));
     // mg->AddFile(Form("%s",fName));
      if( (iMerge-start)%pack ==0){
        mg->OutputFile(Form("AnalysisResults_%03d.root",offset+(iMerge-start-1)/pack));
        printf("Merging up to: %s\n",fName);
        mg->Merge();
        mg->Reset();
        gROOT->Reset();
      }
    }

  }
  if(mg->GetMergeList()->GetEntries()){
    cout<<"Merging last "<<mg->GetMergeList()->GetEntries()<<" files \n";
    mg->OutputFile(Form("AnalysisResults_%03d.root",offset+(iMerge-start)/pack));
    printf("Merging up to: %s\n",fName);
    mg->Merge();
    mg->Reset();
  }
  delete mg;
  return;
}
