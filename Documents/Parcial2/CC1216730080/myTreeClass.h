//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat May 27 11:55:13 2023 by ROOT version 6.28/00
// from TTree tree/Bc mass data
// found on file: Bc_data.root
//////////////////////////////////////////////////////////

#ifndef myTreeClass_h
#define myTreeClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class myTreeClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        Bcmass;
   Double_t        Tau;
   Double_t        TauError;

   // List of branches
   TBranch        *b_Bcmass;   //!
   TBranch        *b_Tau;   //!
   TBranch        *b_TauError;   //!

   myTreeClass(TTree *tree=0);
   virtual ~myTreeClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef myTreeClass_cxx
myTreeClass::myTreeClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Bc_data.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Bc_data.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

myTreeClass::~myTreeClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myTreeClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myTreeClass::LoadTree(Long64_t entry)
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

void myTreeClass::Init(TTree *tree)
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

   fChain->SetBranchAddress("Bcmass", &Bcmass, &b_Bcmass);
   fChain->SetBranchAddress("Tau", &Tau, &b_Tau);
   fChain->SetBranchAddress("TauError", &TauError, &b_TauError);
   Notify();
}

Bool_t myTreeClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myTreeClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myTreeClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef myTreeClass_cxx
