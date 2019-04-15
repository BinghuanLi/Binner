//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr 14 14:38:37 2019 by ROOT version 6.16/00
// from TTree mvaTree/mvaTree
// found on file: ttHnobb_SigRegion.root
//////////////////////////////////////////////////////////

#ifndef mvaTree_h
#define mvaTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class mvaTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Long64_t        nEvent;
   Float_t         Bin2l;
   Float_t         mvaOutput_2lss_ttV;
   Float_t         mvaOutput_2lss_ttbar;
   Float_t         n_presel_jet;
   Float_t         EventWeight;
   Float_t         passTrigCut;
   Float_t         passMassllCut;
   Float_t         passTauNCut;
   Float_t         passZvetoCut;
   Float_t         passMetLDCut;
   Float_t         passTightChargeCut;
   Float_t         passLepTightNCut;
   Float_t         passGenMatchCut;

   // List of branches
   TBranch        *b_nEvent;   //!
   TBranch        *b_Bin2l;   //!
   TBranch        *b_mvaOutput_2lss_ttV;   //!
   TBranch        *b_mvaOutput_2lss_ttbar;   //!
   TBranch        *b_n_presel_jet;   //!
   TBranch        *b_EventWeight;   //!
   TBranch        *b_passTrigCut;   //!
   TBranch        *b_passMassllCut;   //!
   TBranch        *b_passTauNCut;   //!
   TBranch        *b_passZvetoCut;   //!
   TBranch        *b_passMetLDCut;   //!
   TBranch        *b_passTightChargeCut;   //!
   TBranch        *b_passLepTightNCut;   //!
   TBranch        *b_passGenMatchCut;   //!

   mvaTree(TTree *tree=0);
   virtual ~mvaTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mvaTree_cxx
mvaTree::mvaTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ttHnobb_SigRegion.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ttHnobb_SigRegion.root");
      }
      f->GetObject("mvaTree",tree);

   }
   Init(tree);
}

mvaTree::~mvaTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mvaTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mvaTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void mvaTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nEvent", &nEvent, &b_nEvent);
   fChain->SetBranchAddress("Bin2l", &Bin2l, &b_Bin2l);
   fChain->SetBranchAddress("mvaOutput_2lss_ttV", &mvaOutput_2lss_ttV, &b_mvaOutput_2lss_ttV);
   fChain->SetBranchAddress("mvaOutput_2lss_ttbar", &mvaOutput_2lss_ttbar, &b_mvaOutput_2lss_ttbar);
   fChain->SetBranchAddress("n_presel_jet", &n_presel_jet, &b_n_presel_jet);
   fChain->SetBranchAddress("EventWeight", &EventWeight, &b_EventWeight);
   fChain->SetBranchAddress("passTrigCut", &passTrigCut, &b_passTrigCut);
   fChain->SetBranchAddress("passMassllCut", &passMassllCut, &b_passMassllCut);
   fChain->SetBranchAddress("passTauNCut", &passTauNCut, &b_passTauNCut);
   fChain->SetBranchAddress("passZvetoCut", &passZvetoCut, &b_passZvetoCut);
   fChain->SetBranchAddress("passMetLDCut", &passMetLDCut, &b_passMetLDCut);
   fChain->SetBranchAddress("passTightChargeCut", &passTightChargeCut, &b_passTightChargeCut);
   fChain->SetBranchAddress("passLepTightNCut", &passLepTightNCut, &b_passLepTightNCut);
   fChain->SetBranchAddress("passGenMatchCut", &passGenMatchCut, &b_passGenMatchCut);
   Notify();
}

Bool_t mvaTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mvaTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mvaTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mvaTree_cxx
