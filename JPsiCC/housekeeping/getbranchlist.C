/*
Housekeeping script for getting a list of branches from a ROOT file.
Example usage:
    root -l -b -q 'getbranchlist.C("<fileName>", "<treeName>")' > branchlist.txt
*/

void getbranchlist(const char* fileName, const char* treeName)
{
    TFile* f = TFile::Open(fileName, "READ");
    TTree* tree = (TTree*) f->Get(treeName);
    TObjArray* branchlist = tree->GetListOfBranches();
    branchlist->Print();
}