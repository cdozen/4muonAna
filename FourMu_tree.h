//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 20 21:16:21 2019 by ROOT version 6.16/00
// from TTree FourMu_tree/A ROOT tree
// found on file: fourMuMass_tree.root
//////////////////////////////////////////////////////////

#ifndef FourMu_tree_h
#define FourMu_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"

class FourMu_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Lumi;
   Int_t           Event;
   Double_t        Mass4mu;
   Double_t        Pt4mu;
   Double_t        Eta4mu;
   Double_t        Phi4mu;
   Double_t        Mass_upsilon;
   Double_t        Pt_upsilon;
   Double_t        Eta_upsilon;
   Double_t        Y_upsilon;
   Double_t        Phi_upsilon;
   TLorentzVector  *fourMuFit;
   TLorentzVector  *mu1;
   TLorentzVector  *mu2;
   TLorentzVector  *mu3;
   TLorentzVector  *mu4;
   TLorentzVector  *mu12;
   TLorentzVector  *mu34;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Lumi;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_Mass4mu;   //!
   TBranch        *b_Pt4mu;   //!
   TBranch        *b_Eta4mu;   //!
   TBranch        *b_Phi4mu;   //!
   TBranch        *b_Mass_upsilon;   //!
   TBranch        *b_Pt_upsilon;   //!
   TBranch        *b_Eta_upsilon;   //!
   TBranch        *b_Y_upsilon;   //!
   TBranch        *b_Phi_upsilon;   //!
   TBranch        *b_fourMuFit;   //!
   TBranch        *b_mu1;   //!
   TBranch        *b_mu2;   //!
   TBranch        *b_mu3;   //!
   TBranch        *b_mu4;   //!
   TBranch        *b_mu12;   //!
   TBranch        *b_mu34;   //!

   FourMu_tree(TTree *tree=0);
   virtual ~FourMu_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef FourMu_tree_cxx
FourMu_tree::FourMu_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("fourMuMass_tree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("fourMuMass_tree.root");
      }
      f->GetObject("FourMu_tree",tree);

   }
   Init(tree);
}

FourMu_tree::~FourMu_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FourMu_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FourMu_tree::LoadTree(Long64_t entry)
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

void FourMu_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   fourMuFit = 0;
   mu1 = 0;
   mu2 = 0;
   mu3 = 0;
   mu4 = 0;
   mu12 = 0;
   mu34 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Lumi", &Lumi, &b_Lumi);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("Mass4mu", &Mass4mu, &b_Mass4mu);
   fChain->SetBranchAddress("Pt4mu", &Pt4mu, &b_Pt4mu);
   fChain->SetBranchAddress("Eta4mu", &Eta4mu, &b_Eta4mu);
   fChain->SetBranchAddress("Phi4mu", &Phi4mu, &b_Phi4mu);
   fChain->SetBranchAddress("Mass_upsilon", &Mass_upsilon, &b_Mass_upsilon);
   fChain->SetBranchAddress("Pt_upsilon", &Pt_upsilon, &b_Pt_upsilon);
   fChain->SetBranchAddress("Eta_upsilon", &Eta_upsilon, &b_Eta_upsilon);
   fChain->SetBranchAddress("Y_upsilon", &Y_upsilon, &b_Y_upsilon);
   fChain->SetBranchAddress("Phi_upsilon", &Phi_upsilon, &b_Phi_upsilon);
   fChain->SetBranchAddress("fourMuFit", &fourMuFit, &b_fourMuFit);
   fChain->SetBranchAddress("mu1", &mu1, &b_mu1);
   fChain->SetBranchAddress("mu2", &mu2, &b_mu2);
   fChain->SetBranchAddress("mu3", &mu3, &b_mu3);
   fChain->SetBranchAddress("mu4", &mu4, &b_mu4);
   fChain->SetBranchAddress("mu12", &mu12, &b_mu12);
   fChain->SetBranchAddress("mu34", &mu34, &b_mu34);
   Notify();
}

Bool_t FourMu_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void FourMu_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FourMu_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef FourMu_tree_cxx
