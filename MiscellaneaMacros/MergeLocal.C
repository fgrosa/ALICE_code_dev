/////////////////////////////////////////////////////////////////
// Macro to Merge Local Files                                  //
/////////////////////////////////////////////////////////////////

void MergeLocalFirstSetp(Int_t nFiles=70, Int_t startingfile=0, Int_t part=0){

  TFileMerger* mg= new TFileMerger();  
  for(Int_t iFile=0; iFile<nFiles; iFile++)
      mg->AddFile(Form("AnalysisResults_%03d.root",startingfile+iFile));
  mg->OutputFile(Form("AnalysisResultsPbPb2015_020_%d.root",part));
  mg->Merge();
  return;
}
void MergeLocalSecondSetp(Int_t nFiles=7){
  TFileMerger* mg= new TFileMerger();  
  for(Int_t iFile=0; iFile<nFiles; iFile++)
      mg->AddFile(Form("AnalysisResultsPbPb2015_020_%d.root",iFile));
  mg->OutputFile("AnalysisResultsPbPb2015_020.root");
  mg->Merge();
  return;
}
