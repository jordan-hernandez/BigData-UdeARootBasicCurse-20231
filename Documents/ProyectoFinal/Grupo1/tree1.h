//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun  3 17:56:45 2023 by ROOT version 6.28/00
// from TTree tree1/tree1
// found on file: ToyMC.root
//////////////////////////////////////////////////////////

#ifndef tree1_h
#define tree1_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class tree1 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        Mupull;
   Double_t        mu;
   Double_t        muE;
   Double_t        g1;
   Double_t        g1E;
   Double_t        g2;
   Double_t        g2E;
   Int_t           status;
   Int_t           covQual;
   Int_t           badNll;

   // List of branches
   TBranch        *b_Mupull;   //!
   TBranch        *b_mu;   //!
   TBranch        *b_muE;   //!
   TBranch        *b_g1;   //!
   TBranch        *b_g1E;   //!
   TBranch        *b_g2;   //!
   TBranch        *b_g2E;   //!
   TBranch        *b_status;   //!
   TBranch        *b_covQual;   //!
   TBranch        *b_badNll;   //!

   tree1(TTree *tree=0);
   virtual ~tree1();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tree1_cxx
tree1::tree1(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ToyMC.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ToyMC.root");
      }
      f->GetObject("tree1",tree);

   }
   Init(tree);
}

tree1::~tree1()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tree1::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tree1::LoadTree(Long64_t entry)
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

void tree1::Init(TTree *tree)
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

   fChain->SetBranchAddress("Mupull", &Mupull, &b_Mupull);
   fChain->SetBranchAddress("mu", &mu, &b_mu);
   fChain->SetBranchAddress("muE", &muE, &b_muE);
   fChain->SetBranchAddress("g1", &g1, &b_g1);
   fChain->SetBranchAddress("g1E", &g1E, &b_g1E);
   fChain->SetBranchAddress("g2", &g2, &b_g2);
   fChain->SetBranchAddress("g2E", &g2E, &b_g2E);
   fChain->SetBranchAddress("status", &status, &b_status);
   fChain->SetBranchAddress("covQual", &covQual, &b_covQual);
   fChain->SetBranchAddress("badNll", &badNll, &b_badNll);
   Notify();
}

Bool_t tree1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tree1::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tree1::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tree1_cxx
