//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct  9 14:51:52 2018 by ROOT version 6.12/07
// from TTree oniaTree/Tree of MuMuGamma
// found on file: root://cmsxrootd.fnal.gov//store/user/zhenhu/MuOnia/Rootuple_2018.root
//////////////////////////////////////////////////////////

#ifndef readTree2018_h
#define readTree2018_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"
#include "vector"
#include "vector"

class readTree2018 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           event;
   Int_t           irank;
   Int_t           trigger;
   Int_t           numPrimaryVertices;
   Float_t         pv_x;
   Float_t         pv_y;
   Float_t         pv_z;
   TLorentzVector  *mu1_p4;
   TLorentzVector  *mu2_p4;
   Int_t           mu1Charge;
   Int_t           mu2Charge;
   Float_t         mu1_d0;
   Float_t         mu1_d0err;
   Float_t         mu2_d0;
   Float_t         mu2_d0err;
   Float_t         mu1_dz;
   Float_t         mu1_dzerr;
   Float_t         mu2_dz;
   Float_t         mu2_dzerr;
   Float_t         mu1_vz;
   Float_t         mu2_vz;
   TLorentzVector  *dimuon_p4;
   Float_t         mumufit_Mass;
   Float_t         mumufit_MassErr;
   Float_t         mumufit_VtxCL;
   Float_t         mumufit_VtxCL2;
   Float_t         mumufit_DecayVtxX;
   Float_t         mumufit_DecayVtxY;
   Float_t         mumufit_DecayVtxZ;
   Float_t         mumufit_DecayVtxXE;
   Float_t         mumufit_DecayVtxYE;
   Float_t         mumufit_DecayVtxZE;
   TLorentzVector  *mumufit_p4;
   vector<float>   *fourMuFit_Mass_allComb_mix;
   Float_t         fourMuFit_Mass_mix;
   Float_t         fourMuFit_MassErr_mix;
   Float_t         fourMuFit_VtxX_mix;
   Float_t         fourMuFit_VtxY_mix;
   Float_t         fourMuFit_VtxZ_mix;
   Float_t         fourMuFit_VtxProb_mix;
   Float_t         fourMuFit_Chi2_mix;
   Int_t           fourMuFit_ndof_mix;
   Int_t           fourMuFit_3plus1_mix;
   TLorentzVector  *fourMuFit_p4_mix;
   TLorentzVector  *fourMuFit_mu1p4_mix;
   TLorentzVector  *fourMuFit_mu2p4_mix;
   TLorentzVector  *fourMuFit_mu3p4_mix;
   TLorentzVector  *fourMuFit_mu4p4_mix;
   Int_t           mu3Charge_mix;
   Int_t           mu4Charge_mix;
   TLorentzVector  *mu3_p4_mix;
   TLorentzVector  *mu4_p4_mix;
   Float_t         mu3_d0_mix;
   Float_t         mu3_d0err_mix;
   Float_t         mu4_d0_mix;
   Float_t         mu4_d0err_mix;
   Float_t         mu3_dz_mix;
   Float_t         mu3_dzerr_mix;
   Float_t         mu4_dz_mix;
   Float_t         mu4_dzerr_mix;
   Float_t         fourMuFit_Mass_mix3evts;
   Float_t         fourMuFit_VtxProb_mix3evts;
   TLorentzVector  *fourMuFit_p4_mix3evts;
   vector<float>   *fourMuFit_Mass_allComb;
   vector<float>   *v_mumufit_Mass;
   vector<float>   *mu1_trg_dR;
   vector<float>   *mu2_trg_dR;
   vector<float>   *fourMuFit_Mass;
   vector<float>   *fourMuFit_MassErr;
   vector<float>   *fourMuFit_Pt;
   vector<float>   *fourMuFit_Eta;
   vector<float>   *fourMuFit_Phi;
   vector<float>   *fourMuFit_VtxX;
   vector<float>   *fourMuFit_VtxY;
   vector<float>   *fourMuFit_VtxZ;
   vector<float>   *fourMuFit_VtxProb;
   vector<float>   *fourMuFit_Chi2;
   vector<int>     *fourMuFit_ndof;
   vector<float>   *fourMuFit_mu1Pt;
   vector<float>   *fourMuFit_mu1Eta;
   vector<float>   *fourMuFit_mu1Phi;
   vector<float>   *fourMuFit_mu1E;
   vector<float>   *fourMuFit_mu2Pt;
   vector<float>   *fourMuFit_mu2Eta;
   vector<float>   *fourMuFit_mu2Phi;
   vector<float>   *fourMuFit_mu2E;
   vector<float>   *fourMuFit_mu3Pt;
   vector<float>   *fourMuFit_mu3Eta;
   vector<float>   *fourMuFit_mu3Phi;
   vector<float>   *fourMuFit_mu3E;
   vector<float>   *fourMuFit_mu4Pt;
   vector<float>   *fourMuFit_mu4Eta;
   vector<float>   *fourMuFit_mu4Phi;
   vector<float>   *fourMuFit_mu4E;
   vector<float>   *fourMuFit_mu3_trg_dR;
   vector<float>   *fourMuFit_mu4_trg_dR;
   vector<float>   *mu3_Pt;
   vector<float>   *mu3_Eta;
   vector<float>   *mu3_Phi;
   vector<float>   *mu3_E;
   vector<float>   *mu4_Pt;
   vector<float>   *mu4_Eta;
   vector<float>   *mu4_Phi;
   vector<float>   *mu4_E;
   vector<int>     *mu3Charge;
   vector<int>     *mu4Charge;
   vector<float>   *mu3_d0;
   vector<float>   *mu3_d0err;
   vector<float>   *mu4_d0;
   vector<float>   *mu4_d0err;
   vector<float>   *mu3_dz;
   vector<float>   *mu3_dzerr;
   vector<float>   *mu4_dz;
   vector<float>   *mu4_dzerr;
   vector<int>     *mu1_Tight;
   vector<int>     *mu2_Tight;
   vector<int>     *mu3_Tight;
   vector<int>     *mu4_Tight;
   vector<int>     *mu1_Medium;
   vector<int>     *mu2_Medium;
   vector<int>     *mu3_Medium;
   vector<int>     *mu4_Medium;
   vector<int>     *mu1_pdgID;
   vector<int>     *mu2_pdgID;
   vector<int>     *mu3_pdgID;
   vector<int>     *mu4_pdgID;
   TLorentzVector  *mu1_p4_bestYMass;
   TLorentzVector  *mu2_p4_bestYMass;
   Int_t           mu1Charge_bestYMass;
   Int_t           mu2Charge_bestYMass;
   Float_t         mu1_d0_bestYMass;
   Float_t         mu1_d0err_bestYMass;
   Float_t         mu2_d0_bestYMass;
   Float_t         mu2_d0err_bestYMass;
   Float_t         mu1_dz_bestYMass;
   Float_t         mu1_dzerr_bestYMass;
   Float_t         mu2_dz_bestYMass;
   Float_t         mu2_dzerr_bestYMass;
   TLorentzVector  *dimuon_p4_bestYMass;
   Float_t         mumufit_Mass_bestYMass;
   Float_t         mumufit_MassErr_bestYMass;
   Float_t         mumufit_VtxCL_bestYMass;
   Float_t         mumufit_VtxCL2_bestYMass;
   Float_t         mumufit_DecayVtxX_bestYMass;
   Float_t         mumufit_DecayVtxY_bestYMass;
   Float_t         mumufit_DecayVtxZ_bestYMass;
   Float_t         mumufit_DecayVtxXE_bestYMass;
   Float_t         mumufit_DecayVtxYE_bestYMass;
   Float_t         mumufit_DecayVtxZE_bestYMass;
   TLorentzVector  *mumufit_p4_bestYMass;
   Int_t           bestVertex_and_bestYMass;
   vector<float>   *fourMuFit_Mass_allComb_mix_bestYMass;
   Float_t         fourMuFit_Mass_mix_bestYMass;
   Float_t         fourMuFit_MassErr_mix_bestYMass;
   Float_t         fourMuFit_VtxX_mix_bestYMass;
   Float_t         fourMuFit_VtxY_mix_bestYMass;
   Float_t         fourMuFit_VtxZ_mix_bestYMass;
   Float_t         fourMuFit_VtxProb_mix_bestYMass;
   Float_t         fourMuFit_Chi2_mix_bestYMass;
   Int_t           fourMuFit_ndof_mix_bestYMass;
   Int_t           fourMuFit_3plus1_mix_bestYMass;
   TLorentzVector  *fourMuFit_p4_mix_bestYMass;
   TLorentzVector  *fourMuFit_mu1p4_mix_bestYMass;
   TLorentzVector  *fourMuFit_mu2p4_mix_bestYMass;
   TLorentzVector  *fourMuFit_mu3p4_mix_bestYMass;
   TLorentzVector  *fourMuFit_mu4p4_mix_bestYMass;
   Int_t           mu3Charge_mix_bestYMass;
   Int_t           mu4Charge_mix_bestYMass;
   TLorentzVector  *mu3_p4_mix_bestYMass;
   TLorentzVector  *mu4_p4_mix_bestYMass;
   Float_t         mu3_d0_mix_bestYMass;
   Float_t         mu3_d0err_mix_bestYMass;
   Float_t         mu4_d0_mix_bestYMass;
   Float_t         mu4_d0err_mix_bestYMass;
   Float_t         mu3_dz_mix_bestYMass;
   Float_t         mu3_dzerr_mix_bestYMass;
   Float_t         mu4_dz_mix_bestYMass;
   Float_t         mu4_dzerr_mix_bestYMass;
   vector<float>   *fourMuFit_Mass_allComb_bestYMass;
   Float_t         fourMuFit_Mass_bestYMass;
   Float_t         fourMuFit_MassErr_bestYMass;
   Float_t         fourMuFit_VtxX_bestYMass;
   Float_t         fourMuFit_VtxY_bestYMass;
   Float_t         fourMuFit_VtxZ_bestYMass;
   Float_t         fourMuFit_VtxProb_bestYMass;
   Float_t         fourMuFit_Chi2_bestYMass;
   Int_t           fourMuFit_ndof_bestYMass;
   TLorentzVector  *fourMuFit_p4_bestYMass;
   TLorentzVector  *fourMuFit_mu1p4_bestYMass;
   TLorentzVector  *fourMuFit_mu2p4_bestYMass;
   TLorentzVector  *fourMuFit_mu3p4_bestYMass;
   TLorentzVector  *fourMuFit_mu4p4_bestYMass;
   Int_t           mu3Charge_bestYMass;
   Int_t           mu4Charge_bestYMass;
   TLorentzVector  *mu3_p4_bestYMass;
   TLorentzVector  *mu4_p4_bestYMass;
   Float_t         mu3_d0_bestYMass;
   Float_t         mu3_d0err_bestYMass;
   Float_t         mu4_d0_bestYMass;
   Float_t         mu4_d0err_bestYMass;
   Float_t         mu3_dz_bestYMass;
   Float_t         mu3_dzerr_bestYMass;
   Float_t         mu4_dz_bestYMass;
   Float_t         mu4_dzerr_bestYMass;
   Int_t           mu1_Tight_bestYMass;
   Int_t           mu2_Tight_bestYMass;
   Int_t           mu3_Tight_bestYMass;
   Int_t           mu4_Tight_bestYMass;
   Int_t           mu3_pdgID_bestYMass;
   Int_t           mu4_pdgID_bestYMass;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_irank;   //!
   TBranch        *b_trigger;   //!
   TBranch        *b_numPrimaryVertices;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_mu1_p4;   //!
   TBranch        *b_mu2_p4;   //!
   TBranch        *b_mu1Charge;   //!
   TBranch        *b_mu2Charge;   //!
   TBranch        *b_mu1_d0;   //!
   TBranch        *b_mu1_d0err;   //!
   TBranch        *b_mu2_d0;   //!
   TBranch        *b_mu2_d0err;   //!
   TBranch        *b_mu1_dz;   //!
   TBranch        *b_mu1_dzerr;   //!
   TBranch        *b_mu2_dz;   //!
   TBranch        *b_mu2_dzerr;   //!
   TBranch        *b_mu1_vz;   //!
   TBranch        *b_mu2_vz;   //!
   TBranch        *b_dimuon_p4;   //!
   TBranch        *b_mumufit_Mass;   //!
   TBranch        *b_mumufit_MassErr;   //!
   TBranch        *b_mumufit_VtxCL;   //!
   TBranch        *b_mumufit_VtxCL2;   //!
   TBranch        *b_mumufit_DecayVtxX;   //!
   TBranch        *b_mumufit_DecayVtxY;   //!
   TBranch        *b_mumufit_DecayVtxZ;   //!
   TBranch        *b_mumufit_DecayVtxXE;   //!
   TBranch        *b_mumufit_DecayVtxYE;   //!
   TBranch        *b_mumufit_DecayVtxZE;   //!
   TBranch        *b_mumufit_p4;   //!
   TBranch        *b_fourMuFit_Mass_allComb_mix;   //!
   TBranch        *b_fourMuFit_Mass_mix;   //!
   TBranch        *b_fourMuFit_MassErr_mix;   //!
   TBranch        *b_fourMuFit_VtxX_mix;   //!
   TBranch        *b_fourMuFit_VtxY_mix;   //!
   TBranch        *b_fourMuFit_VtxZ_mix;   //!
   TBranch        *b_fourMuFit_VtxProb_mix;   //!
   TBranch        *b_fourMuFit_Chi2_mix;   //!
   TBranch        *b_fourMuFit_ndof_mix;   //!
   TBranch        *b_fourMuFit_3plus1_mix;   //!
   TBranch        *b_fourMuFit_p4_mix;   //!
   TBranch        *b_fourMuFit_mu1p4_mix;   //!
   TBranch        *b_fourMuFit_mu2p4_mix;   //!
   TBranch        *b_fourMuFit_mu3p4_mix;   //!
   TBranch        *b_fourMuFit_mu4p4_mix;   //!
   TBranch        *b_mu3Charge_mix;   //!
   TBranch        *b_mu4Charge_mix;   //!
   TBranch        *b_mu3_p4_mix;   //!
   TBranch        *b_mu4_p4_mix;   //!
   TBranch        *b_mu3_d0_mix;   //!
   TBranch        *b_mu3_d0err_mix;   //!
   TBranch        *b_mu4_d0_mix;   //!
   TBranch        *b_mu4_d0err_mix;   //!
   TBranch        *b_mu3_dz_mix;   //!
   TBranch        *b_mu3_dzerr_mix;   //!
   TBranch        *b_mu4_dz_mix;   //!
   TBranch        *b_mu4_dzerr_mix;   //!
   TBranch        *b_fourMuFit_Mass_mix3evts;   //!
   TBranch        *b_fourMuFit_VtxProb_mix3evts;   //!
   TBranch        *b_fourMuFit_p4_mix3evts;   //!
   TBranch        *b_fourMuFit_Mass_allComb;   //!
   TBranch        *b_v_mumufit_Mass;   //!
   TBranch        *b_mu1_trg_dR;
   TBranch        *b_mu2_trg_dR;
   TBranch        *b_fourMuFit_Mass;   //!
   TBranch        *b_fourMuFit_MassErr;   //!
   TBranch        *b_fourMuFit_Pt;   //!
   TBranch        *b_fourMuFit_Eta;   //!
   TBranch        *b_fourMuFit_Phi;   //!
   TBranch        *b_fourMuFit_VtxX;   //!
   TBranch        *b_fourMuFit_VtxY;   //!
   TBranch        *b_fourMuFit_VtxZ;   //!
   TBranch        *b_fourMuFit_VtxProb;   //!
   TBranch        *b_fourMuFit_Chi2;   //!
   TBranch        *b_fourMuFit_ndof;   //!
   TBranch        *b_fourMuFit_mu1Pt;   //!
   TBranch        *b_fourMuFit_mu1Eta;   //!
   TBranch        *b_fourMuFit_mu1Phi;   //!
   TBranch        *b_fourMuFit_mu1E;   //!
   TBranch        *b_fourMuFit_mu2Pt;   //!
   TBranch        *b_fourMuFit_mu2Eta;   //!
   TBranch        *b_fourMuFit_mu2Phi;   //!
   TBranch        *b_fourMuFit_mu2E;   //!
   TBranch        *b_fourMuFit_mu3Pt;   //!
   TBranch        *b_fourMuFit_mu3Eta;   //!
   TBranch        *b_fourMuFit_mu3Phi;   //!
   TBranch        *b_fourMuFit_mu3E;   //!
   TBranch        *b_fourMuFit_mu4Pt;   //!
   TBranch        *b_fourMuFit_mu4Eta;   //!
   TBranch        *b_fourMuFit_mu4Phi;   //!
   TBranch        *b_fourMuFit_mu4E;   //!
   TBranch        *b_fourMuFit_mu3_trg_dR;
   TBranch        *b_fourMuFit_mu4_trg_dR;
   TBranch        *b_mu3_Pt;   //!
   TBranch        *b_mu3_Eta;   //!
   TBranch        *b_mu3_Phi;   //!
   TBranch        *b_mu3_E;   //!
   TBranch        *b_mu4_Pt;   //!
   TBranch        *b_mu4_Eta;   //!
   TBranch        *b_mu4_Phi;   //!
   TBranch        *b_mu4_E;   //!
   TBranch        *b_mu3Charge;   //!
   TBranch        *b_mu4Charge;   //!
   TBranch        *b_mu3_d0;   //!
   TBranch        *b_mu3_d0err;   //!
   TBranch        *b_mu4_d0;   //!
   TBranch        *b_mu4_d0err;   //!
   TBranch        *b_mu3_dz;   //!
   TBranch        *b_mu3_dzerr;   //!
   TBranch        *b_mu4_dz;   //!
   TBranch        *b_mu4_dzerr;   //!
   TBranch        *b_mu1_Tight;   //!
   TBranch        *b_mu2_Tight;   //!
   TBranch        *b_mu3_Tight;   //!
   TBranch        *b_mu4_Tight;   //!
   TBranch        *b_mu1_Medium;   //!
   TBranch        *b_mu2_Medium;   //!
   TBranch        *b_mu3_Medium;   //!
   TBranch        *b_mu4_Medium;   //!
   TBranch        *b_mu1_pdgID;   //!
   TBranch        *b_mu2_pdgID;   //!
   TBranch        *b_mu3_pdgID;   //!
   TBranch        *b_mu4_pdgID;   //!
   TBranch        *b_mu1_p4_bestYMass;   //!
   TBranch        *b_mu2_p4_bestYMass;   //!
   TBranch        *b_mu1Charge_bestYMass;   //!
   TBranch        *b_mu2Charge_bestYMass;   //!
   TBranch        *b_mu1_d0_bestYMass;   //!
   TBranch        *b_mu1_d0err_bestYMass;   //!
   TBranch        *b_mu2_d0_bestYMass;   //!
   TBranch        *b_mu2_d0err_bestYMass;   //!
   TBranch        *b_mu1_dz_bestYMass;   //!
   TBranch        *b_mu1_dzerr_bestYMass;   //!
   TBranch        *b_mu2_dz_bestYMass;   //!
   TBranch        *b_mu2_dzerr_bestYMass;   //!
   TBranch        *b_dimuon_p4_bestYMass;   //!
   TBranch        *b_mumufit_Mass_bestYMass;   //!
   TBranch        *b_mumufit_MassErr_bestYMass;   //!
   TBranch        *b_mumufit_VtxCL_bestYMass;   //!
   TBranch        *b_mumufit_VtxCL2_bestYMass;   //!
   TBranch        *b_mumufit_DecayVtxX_bestYMass;   //!
   TBranch        *b_mumufit_DecayVtxY_bestYMass;   //!
   TBranch        *b_mumufit_DecayVtxZ_bestYMass;   //!
   TBranch        *b_mumufit_DecayVtxXE_bestYMass;   //!
   TBranch        *b_mumufit_DecayVtxYE_bestYMass;   //!
   TBranch        *b_mumufit_DecayVtxZE_bestYMass;   //!
   TBranch        *b_mumufit_p4_bestYMass;   //!
   TBranch        *b_bestVertex_and_bestYMass;   //!
   TBranch        *b_fourMuFit_Mass_allComb_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_Mass_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_MassErr_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_VtxX_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_VtxY_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_VtxZ_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_VtxProb_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_Chi2_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_ndof_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_3plus1_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_p4_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_mu1p4_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_mu2p4_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_mu3p4_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_mu4p4_mix_bestYMass;   //!
   TBranch        *b_mu3Charge_mix_bestYMass;   //!
   TBranch        *b_mu4Charge_mix_bestYMass;   //!
   TBranch        *b_mu3_p4_mix_bestYMass;   //!
   TBranch        *b_mu4_p4_mix_bestYMass;   //!
   TBranch        *b_mu3_d0_mix_bestYMass;   //!
   TBranch        *b_mu3_d0err_mix_bestYMass;   //!
   TBranch        *b_mu4_d0_mix_bestYMass;   //!
   TBranch        *b_mu4_d0err_mix_bestYMass;   //!
   TBranch        *b_mu3_dz_mix_bestYMass;   //!
   TBranch        *b_mu3_dzerr_mix_bestYMass;   //!
   TBranch        *b_mu4_dz_mix_bestYMass;   //!
   TBranch        *b_mu4_dzerr_mix_bestYMass;   //!
   TBranch        *b_fourMuFit_Mass_allComb_bestYMass;   //!
   TBranch        *b_fourMuFit_Mass_bestYMass;   //!
   TBranch        *b_fourMuFit_MassErr_bestYMass;   //!
   TBranch        *b_fourMuFit_VtxX_bestYMass;   //!
   TBranch        *b_fourMuFit_VtxY_bestYMass;   //!
   TBranch        *b_fourMuFit_VtxZ_bestYMass;   //!
   TBranch        *b_fourMuFit_VtxProb_bestYMass;   //!
   TBranch        *b_fourMuFit_Chi2_bestYMass;   //!
   TBranch        *b_fourMuFit_ndof_bestYMass;   //!
   TBranch        *b_fourMuFit_p4_bestYMass;   //!
   TBranch        *b_fourMuFit_mu1p4_bestYMass;   //!
   TBranch        *b_fourMuFit_mu2p4_bestYMass;   //!
   TBranch        *b_fourMuFit_mu3p4_bestYMass;   //!
   TBranch        *b_fourMuFit_mu4p4_bestYMass;   //!
   TBranch        *b_mu3Charge_bestYMass;   //!
   TBranch        *b_mu4Charge_bestYMass;   //!
   TBranch        *b_mu3_p4_bestYMass;   //!
   TBranch        *b_mu4_p4_bestYMass;   //!
   TBranch        *b_mu3_d0_bestYMass;   //!
   TBranch        *b_mu3_d0err_bestYMass;   //!
   TBranch        *b_mu4_d0_bestYMass;   //!
   TBranch        *b_mu4_d0err_bestYMass;   //!
   TBranch        *b_mu3_dz_bestYMass;   //!
   TBranch        *b_mu3_dzerr_bestYMass;   //!
   TBranch        *b_mu4_dz_bestYMass;   //!
   TBranch        *b_mu4_dzerr_bestYMass;   //!
   TBranch        *b_mu1_Tight_bestYMass;   //!
   TBranch        *b_mu2_Tight_bestYMass;   //!
   TBranch        *b_mu3_Tight_bestYMass;   //!
   TBranch        *b_mu4_Tight_bestYMass;   //!
   TBranch        *b_mu3_pdgID_bestYMass;   //!
   TBranch        *b_mu4_pdgID_bestYMass;   //!

   readTree2018(TTree *tree=0);
   virtual ~readTree2018();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef readTree2018_cxx
readTree2018::readTree2018(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cmsxrootd.fnal.gov//store/user/zhenhu/MuOnia/Rootuple_2018.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://cmsxrootd.fnal.gov//store/user/zhenhu/MuOnia/Rootuple_2018.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://cmsxrootd.fnal.gov//store/user/zhenhu/MuOnia/Rootuple_2018.root:/rootuple");
      dir->GetObject("oniaTree",tree);

   }
   Init(tree);
}

readTree2018::~readTree2018()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t readTree2018::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t readTree2018::LoadTree(Long64_t entry)
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

void readTree2018::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   mu1_p4 = 0;
   mu2_p4 = 0;
   dimuon_p4 = 0;
   mumufit_p4 = 0;
   fourMuFit_Mass_allComb_mix = 0;
   fourMuFit_p4_mix = 0;
   fourMuFit_mu1p4_mix = 0;
   fourMuFit_mu2p4_mix = 0;
   fourMuFit_mu3p4_mix = 0;
   fourMuFit_mu4p4_mix = 0;
   mu3_p4_mix = 0;
   mu4_p4_mix = 0;
   fourMuFit_p4_mix3evts = 0;
   fourMuFit_Mass_allComb = 0;
   v_mumufit_Mass = 0;
   mu1_trg_dR = 0;
   mu2_trg_dR = 0;
   fourMuFit_Mass = 0;
   fourMuFit_MassErr = 0;
   fourMuFit_Pt = 0;
   fourMuFit_Eta = 0;
   fourMuFit_Phi = 0;
   fourMuFit_VtxX = 0;
   fourMuFit_VtxY = 0;
   fourMuFit_VtxZ = 0;
   fourMuFit_VtxProb = 0;
   fourMuFit_Chi2 = 0;
   fourMuFit_ndof = 0;
   fourMuFit_mu1Pt = 0;
   fourMuFit_mu1Eta = 0;
   fourMuFit_mu1Phi = 0;
   fourMuFit_mu1E = 0;
   fourMuFit_mu2Pt = 0;
   fourMuFit_mu2Eta = 0;
   fourMuFit_mu2Phi = 0;
   fourMuFit_mu2E = 0;
   fourMuFit_mu3Pt = 0;
   fourMuFit_mu3Eta = 0;
   fourMuFit_mu3Phi = 0;
   fourMuFit_mu3E = 0;
   fourMuFit_mu4Pt = 0;
   fourMuFit_mu4Eta = 0;
   fourMuFit_mu4Phi = 0;
   fourMuFit_mu4E = 0;
   fourMuFit_mu3_trg_dR = 0;
   fourMuFit_mu4_trg_dR = 0;
   mu3_Pt = 0;
   mu3_Eta = 0;
   mu3_Phi = 0;
   mu3_E = 0;
   mu4_Pt = 0;
   mu4_Eta = 0;
   mu4_Phi = 0;
   mu4_E = 0;
   mu3Charge = 0;
   mu4Charge = 0;
   mu3_d0 = 0;
   mu3_d0err = 0;
   mu4_d0 = 0;
   mu4_d0err = 0;
   mu3_dz = 0;
   mu3_dzerr = 0;
   mu4_dz = 0;
   mu4_dzerr = 0;
   mu1_Tight = 0;
   mu2_Tight = 0;
   mu3_Tight = 0;
   mu4_Tight = 0;
   mu1_Medium = 0;
   mu2_Medium = 0;
   mu3_Medium = 0;
   mu4_Medium = 0;
   mu1_pdgID = 0;
   mu2_pdgID = 0;
   mu3_pdgID = 0;
   mu4_pdgID = 0;
   mu1_p4_bestYMass = 0;
   mu2_p4_bestYMass = 0;
   dimuon_p4_bestYMass = 0;
   mumufit_p4_bestYMass = 0;
   fourMuFit_Mass_allComb_mix_bestYMass = 0;
   fourMuFit_p4_mix_bestYMass = 0;
   fourMuFit_mu1p4_mix_bestYMass = 0;
   fourMuFit_mu2p4_mix_bestYMass = 0;
   fourMuFit_mu3p4_mix_bestYMass = 0;
   fourMuFit_mu4p4_mix_bestYMass = 0;
   mu3_p4_mix_bestYMass = 0;
   mu4_p4_mix_bestYMass = 0;
   fourMuFit_Mass_allComb_bestYMass = 0;
   fourMuFit_p4_bestYMass = 0;
   fourMuFit_mu1p4_bestYMass = 0;
   fourMuFit_mu2p4_bestYMass = 0;
   fourMuFit_mu3p4_bestYMass = 0;
   fourMuFit_mu4p4_bestYMass = 0;
   mu3_p4_bestYMass = 0;
   mu4_p4_bestYMass = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("irank", &irank, &b_irank);
   fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
   fChain->SetBranchAddress("numPrimaryVertices", &numPrimaryVertices, &b_numPrimaryVertices);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("mu1_p4", &mu1_p4, &b_mu1_p4);
   fChain->SetBranchAddress("mu2_p4", &mu2_p4, &b_mu2_p4);
   fChain->SetBranchAddress("mu1Charge", &mu1Charge, &b_mu1Charge);
   fChain->SetBranchAddress("mu2Charge", &mu2Charge, &b_mu2Charge);
   fChain->SetBranchAddress("mu1_d0", &mu1_d0, &b_mu1_d0);
   fChain->SetBranchAddress("mu1_d0err", &mu1_d0err, &b_mu1_d0err);
   fChain->SetBranchAddress("mu2_d0", &mu2_d0, &b_mu2_d0);
   fChain->SetBranchAddress("mu2_d0err", &mu2_d0err, &b_mu2_d0err);
   fChain->SetBranchAddress("mu1_dz", &mu1_dz, &b_mu1_dz);
   fChain->SetBranchAddress("mu1_dzerr", &mu1_dzerr, &b_mu1_dzerr);
   fChain->SetBranchAddress("mu2_dz", &mu2_dz, &b_mu2_dz);
   fChain->SetBranchAddress("mu2_dzerr", &mu2_dzerr, &b_mu2_dzerr);
   fChain->SetBranchAddress("mu1_vz", &mu1_vz, &b_mu1_vz);
   fChain->SetBranchAddress("mu2_vz", &mu2_vz, &b_mu2_vz);
   fChain->SetBranchAddress("dimuon_p4", &dimuon_p4, &b_dimuon_p4);
   fChain->SetBranchAddress("mumufit_Mass", &mumufit_Mass, &b_mumufit_Mass);
   fChain->SetBranchAddress("mumufit_MassErr", &mumufit_MassErr, &b_mumufit_MassErr);
   fChain->SetBranchAddress("mumufit_VtxCL", &mumufit_VtxCL, &b_mumufit_VtxCL);
   fChain->SetBranchAddress("mumufit_VtxCL2", &mumufit_VtxCL2, &b_mumufit_VtxCL2);
   fChain->SetBranchAddress("mumufit_DecayVtxX", &mumufit_DecayVtxX, &b_mumufit_DecayVtxX);
   fChain->SetBranchAddress("mumufit_DecayVtxY", &mumufit_DecayVtxY, &b_mumufit_DecayVtxY);
   fChain->SetBranchAddress("mumufit_DecayVtxZ", &mumufit_DecayVtxZ, &b_mumufit_DecayVtxZ);
   fChain->SetBranchAddress("mumufit_DecayVtxXE", &mumufit_DecayVtxXE, &b_mumufit_DecayVtxXE);
   fChain->SetBranchAddress("mumufit_DecayVtxYE", &mumufit_DecayVtxYE, &b_mumufit_DecayVtxYE);
   fChain->SetBranchAddress("mumufit_DecayVtxZE", &mumufit_DecayVtxZE, &b_mumufit_DecayVtxZE);
   fChain->SetBranchAddress("mumufit_p4", &mumufit_p4, &b_mumufit_p4);
   fChain->SetBranchAddress("fourMuFit_Mass_allComb_mix", &fourMuFit_Mass_allComb_mix, &b_fourMuFit_Mass_allComb_mix);
   fChain->SetBranchAddress("fourMuFit_Mass_mix", &fourMuFit_Mass_mix, &b_fourMuFit_Mass_mix);
   fChain->SetBranchAddress("fourMuFit_MassErr_mix", &fourMuFit_MassErr_mix, &b_fourMuFit_MassErr_mix);
   fChain->SetBranchAddress("fourMuFit_VtxX_mix", &fourMuFit_VtxX_mix, &b_fourMuFit_VtxX_mix);
   fChain->SetBranchAddress("fourMuFit_VtxY_mix", &fourMuFit_VtxY_mix, &b_fourMuFit_VtxY_mix);
   fChain->SetBranchAddress("fourMuFit_VtxZ_mix", &fourMuFit_VtxZ_mix, &b_fourMuFit_VtxZ_mix);
   fChain->SetBranchAddress("fourMuFit_VtxProb_mix", &fourMuFit_VtxProb_mix, &b_fourMuFit_VtxProb_mix);
   fChain->SetBranchAddress("fourMuFit_Chi2_mix", &fourMuFit_Chi2_mix, &b_fourMuFit_Chi2_mix);
   fChain->SetBranchAddress("fourMuFit_ndof_mix", &fourMuFit_ndof_mix, &b_fourMuFit_ndof_mix);
   fChain->SetBranchAddress("fourMuFit_3plus1_mix", &fourMuFit_3plus1_mix, &b_fourMuFit_3plus1_mix);
   fChain->SetBranchAddress("fourMuFit_p4_mix", &fourMuFit_p4_mix, &b_fourMuFit_p4_mix);
   fChain->SetBranchAddress("fourMuFit_mu1p4_mix", &fourMuFit_mu1p4_mix, &b_fourMuFit_mu1p4_mix);
   fChain->SetBranchAddress("fourMuFit_mu2p4_mix", &fourMuFit_mu2p4_mix, &b_fourMuFit_mu2p4_mix);
   fChain->SetBranchAddress("fourMuFit_mu3p4_mix", &fourMuFit_mu3p4_mix, &b_fourMuFit_mu3p4_mix);
   fChain->SetBranchAddress("fourMuFit_mu4p4_mix", &fourMuFit_mu4p4_mix, &b_fourMuFit_mu4p4_mix);
   fChain->SetBranchAddress("mu3Charge_mix", &mu3Charge_mix, &b_mu3Charge_mix);
   fChain->SetBranchAddress("mu4Charge_mix", &mu4Charge_mix, &b_mu4Charge_mix);
   fChain->SetBranchAddress("mu3_p4_mix", &mu3_p4_mix, &b_mu3_p4_mix);
   fChain->SetBranchAddress("mu4_p4_mix", &mu4_p4_mix, &b_mu4_p4_mix);
   fChain->SetBranchAddress("mu3_d0_mix", &mu3_d0_mix, &b_mu3_d0_mix);
   fChain->SetBranchAddress("mu3_d0err_mix", &mu3_d0err_mix, &b_mu3_d0err_mix);
   fChain->SetBranchAddress("mu4_d0_mix", &mu4_d0_mix, &b_mu4_d0_mix);
   fChain->SetBranchAddress("mu4_d0err_mix", &mu4_d0err_mix, &b_mu4_d0err_mix);
   fChain->SetBranchAddress("mu3_dz_mix", &mu3_dz_mix, &b_mu3_dz_mix);
   fChain->SetBranchAddress("mu3_dzerr_mix", &mu3_dzerr_mix, &b_mu3_dzerr_mix);
   fChain->SetBranchAddress("mu4_dz_mix", &mu4_dz_mix, &b_mu4_dz_mix);
   fChain->SetBranchAddress("mu4_dzerr_mix", &mu4_dzerr_mix, &b_mu4_dzerr_mix);
   fChain->SetBranchAddress("fourMuFit_Mass_mix3evts", &fourMuFit_Mass_mix3evts, &b_fourMuFit_Mass_mix3evts);
   fChain->SetBranchAddress("fourMuFit_VtxProb_mix3evts", &fourMuFit_VtxProb_mix3evts, &b_fourMuFit_VtxProb_mix3evts);
   fChain->SetBranchAddress("fourMuFit_p4_mix3evts", &fourMuFit_p4_mix3evts, &b_fourMuFit_p4_mix3evts);
   fChain->SetBranchAddress("fourMuFit_Mass_allComb", &fourMuFit_Mass_allComb, &b_fourMuFit_Mass_allComb);
   fChain->SetBranchAddress("v_mumufit_Mass", &v_mumufit_Mass, &b_v_mumufit_Mass);
   fChain->SetBranchAddress("mu1_trg_dR",&mu1_trg_dR,&b_mu1_trg_dR);
   fChain->SetBranchAddress("mu2_trg_dR",&mu2_trg_dR,&b_mu2_trg_dR);
   fChain->SetBranchAddress("fourMuFit_Mass", &fourMuFit_Mass, &b_fourMuFit_Mass);
   fChain->SetBranchAddress("fourMuFit_MassErr", &fourMuFit_MassErr, &b_fourMuFit_MassErr);
   fChain->SetBranchAddress("fourMuFit_Pt", &fourMuFit_Pt, &b_fourMuFit_Pt);
   fChain->SetBranchAddress("fourMuFit_Eta", &fourMuFit_Eta, &b_fourMuFit_Eta);
   fChain->SetBranchAddress("fourMuFit_Phi", &fourMuFit_Phi, &b_fourMuFit_Phi);
   fChain->SetBranchAddress("fourMuFit_VtxX", &fourMuFit_VtxX, &b_fourMuFit_VtxX);
   fChain->SetBranchAddress("fourMuFit_VtxY", &fourMuFit_VtxY, &b_fourMuFit_VtxY);
   fChain->SetBranchAddress("fourMuFit_VtxZ", &fourMuFit_VtxZ, &b_fourMuFit_VtxZ);
   fChain->SetBranchAddress("fourMuFit_VtxProb", &fourMuFit_VtxProb, &b_fourMuFit_VtxProb);
   fChain->SetBranchAddress("fourMuFit_Chi2", &fourMuFit_Chi2, &b_fourMuFit_Chi2);
   fChain->SetBranchAddress("fourMuFit_ndof", &fourMuFit_ndof, &b_fourMuFit_ndof);
   fChain->SetBranchAddress("fourMuFit_mu1Pt", &fourMuFit_mu1Pt, &b_fourMuFit_mu1Pt);
   fChain->SetBranchAddress("fourMuFit_mu1Eta", &fourMuFit_mu1Eta, &b_fourMuFit_mu1Eta);
   fChain->SetBranchAddress("fourMuFit_mu1Phi", &fourMuFit_mu1Phi, &b_fourMuFit_mu1Phi);
   fChain->SetBranchAddress("fourMuFit_mu1E", &fourMuFit_mu1E, &b_fourMuFit_mu1E);
   fChain->SetBranchAddress("fourMuFit_mu2Pt", &fourMuFit_mu2Pt, &b_fourMuFit_mu2Pt);
   fChain->SetBranchAddress("fourMuFit_mu2Eta", &fourMuFit_mu2Eta, &b_fourMuFit_mu2Eta);
   fChain->SetBranchAddress("fourMuFit_mu2Phi", &fourMuFit_mu2Phi, &b_fourMuFit_mu2Phi);
   fChain->SetBranchAddress("fourMuFit_mu2E", &fourMuFit_mu2E, &b_fourMuFit_mu2E);
   fChain->SetBranchAddress("fourMuFit_mu3Pt", &fourMuFit_mu3Pt, &b_fourMuFit_mu3Pt);
   fChain->SetBranchAddress("fourMuFit_mu3Eta", &fourMuFit_mu3Eta, &b_fourMuFit_mu3Eta);
   fChain->SetBranchAddress("fourMuFit_mu3Phi", &fourMuFit_mu3Phi, &b_fourMuFit_mu3Phi);
   fChain->SetBranchAddress("fourMuFit_mu3E", &fourMuFit_mu3E, &b_fourMuFit_mu3E);
   fChain->SetBranchAddress("fourMuFit_mu4Pt", &fourMuFit_mu4Pt, &b_fourMuFit_mu4Pt);
   fChain->SetBranchAddress("fourMuFit_mu4Eta", &fourMuFit_mu4Eta, &b_fourMuFit_mu4Eta);
   fChain->SetBranchAddress("fourMuFit_mu4Phi", &fourMuFit_mu4Phi, &b_fourMuFit_mu4Phi);
   fChain->SetBranchAddress("fourMuFit_mu4E", &fourMuFit_mu4E, &b_fourMuFit_mu4E);
   fChain->SetBranchAddress("fourMuFit_mu3_trg_dR",&fourMuFit_mu3_trg_dR,&b_fourMuFit_mu3_trg_dR);
   fChain->SetBranchAddress("fourMuFit_mu4_trg_dR",&fourMuFit_mu4_trg_dR,&b_fourMuFit_mu4_trg_dR);
   fChain->SetBranchAddress("mu3_Pt", &mu3_Pt, &b_mu3_Pt);
   fChain->SetBranchAddress("mu3_Eta", &mu3_Eta, &b_mu3_Eta);
   fChain->SetBranchAddress("mu3_Phi", &mu3_Phi, &b_mu3_Phi);
   fChain->SetBranchAddress("mu3_E", &mu3_E, &b_mu3_E);
   fChain->SetBranchAddress("mu4_Pt", &mu4_Pt, &b_mu4_Pt);
   fChain->SetBranchAddress("mu4_Eta", &mu4_Eta, &b_mu4_Eta);
   fChain->SetBranchAddress("mu4_Phi", &mu4_Phi, &b_mu4_Phi);
   fChain->SetBranchAddress("mu4_E", &mu4_E, &b_mu4_E);
   fChain->SetBranchAddress("mu3Charge", &mu3Charge, &b_mu3Charge);
   fChain->SetBranchAddress("mu4Charge", &mu4Charge, &b_mu4Charge);
   fChain->SetBranchAddress("mu3_d0", &mu3_d0, &b_mu3_d0);
   fChain->SetBranchAddress("mu3_d0err", &mu3_d0err, &b_mu3_d0err);
   fChain->SetBranchAddress("mu4_d0", &mu4_d0, &b_mu4_d0);
   fChain->SetBranchAddress("mu4_d0err", &mu4_d0err, &b_mu4_d0err);
   fChain->SetBranchAddress("mu3_dz", &mu3_dz, &b_mu3_dz);
   fChain->SetBranchAddress("mu3_dzerr", &mu3_dzerr, &b_mu3_dzerr);
   fChain->SetBranchAddress("mu4_dz", &mu4_dz, &b_mu4_dz);
   fChain->SetBranchAddress("mu4_dzerr", &mu4_dzerr, &b_mu4_dzerr);
   fChain->SetBranchAddress("mu1_Tight", &mu1_Tight, &b_mu1_Tight);
   fChain->SetBranchAddress("mu2_Tight", &mu2_Tight, &b_mu2_Tight);
   fChain->SetBranchAddress("mu3_Tight", &mu3_Tight, &b_mu3_Tight);
   fChain->SetBranchAddress("mu4_Tight", &mu4_Tight, &b_mu4_Tight);
   fChain->SetBranchAddress("mu1_Medium", &mu1_Medium, &b_mu1_Medium);
   fChain->SetBranchAddress("mu2_Medium", &mu2_Medium, &b_mu2_Medium);
   fChain->SetBranchAddress("mu3_Medium", &mu3_Medium, &b_mu3_Medium);
   fChain->SetBranchAddress("mu4_Medium", &mu4_Medium, &b_mu4_Medium);
   fChain->SetBranchAddress("mu1_pdgID", &mu1_pdgID, &b_mu1_pdgID);
   fChain->SetBranchAddress("mu2_pdgID", &mu2_pdgID, &b_mu2_pdgID);
   fChain->SetBranchAddress("mu3_pdgID", &mu3_pdgID, &b_mu3_pdgID);
   fChain->SetBranchAddress("mu4_pdgID", &mu4_pdgID, &b_mu4_pdgID);
   fChain->SetBranchAddress("mu1_p4_bestYMass", &mu1_p4_bestYMass, &b_mu1_p4_bestYMass);
   fChain->SetBranchAddress("mu2_p4_bestYMass", &mu2_p4_bestYMass, &b_mu2_p4_bestYMass);
   fChain->SetBranchAddress("mu1Charge_bestYMass", &mu1Charge_bestYMass, &b_mu1Charge_bestYMass);
   fChain->SetBranchAddress("mu2Charge_bestYMass", &mu2Charge_bestYMass, &b_mu2Charge_bestYMass);
   fChain->SetBranchAddress("mu1_d0_bestYMass", &mu1_d0_bestYMass, &b_mu1_d0_bestYMass);
   fChain->SetBranchAddress("mu1_d0err_bestYMass", &mu1_d0err_bestYMass, &b_mu1_d0err_bestYMass);
   fChain->SetBranchAddress("mu2_d0_bestYMass", &mu2_d0_bestYMass, &b_mu2_d0_bestYMass);
   fChain->SetBranchAddress("mu2_d0err_bestYMass", &mu2_d0err_bestYMass, &b_mu2_d0err_bestYMass);
   fChain->SetBranchAddress("mu1_dz_bestYMass", &mu1_dz_bestYMass, &b_mu1_dz_bestYMass);
   fChain->SetBranchAddress("mu1_dzerr_bestYMass", &mu1_dzerr_bestYMass, &b_mu1_dzerr_bestYMass);
   fChain->SetBranchAddress("mu2_dz_bestYMass", &mu2_dz_bestYMass, &b_mu2_dz_bestYMass);
   fChain->SetBranchAddress("mu2_dzerr_bestYMass", &mu2_dzerr_bestYMass, &b_mu2_dzerr_bestYMass);
   fChain->SetBranchAddress("dimuon_p4_bestYMass", &dimuon_p4_bestYMass, &b_dimuon_p4_bestYMass);
   fChain->SetBranchAddress("mumufit_Mass_bestYMass", &mumufit_Mass_bestYMass, &b_mumufit_Mass_bestYMass);
   fChain->SetBranchAddress("mumufit_MassErr_bestYMass", &mumufit_MassErr_bestYMass, &b_mumufit_MassErr_bestYMass);
   fChain->SetBranchAddress("mumufit_VtxCL_bestYMass", &mumufit_VtxCL_bestYMass, &b_mumufit_VtxCL_bestYMass);
   fChain->SetBranchAddress("mumufit_VtxCL2_bestYMass", &mumufit_VtxCL2_bestYMass, &b_mumufit_VtxCL2_bestYMass);
   fChain->SetBranchAddress("mumufit_DecayVtxX_bestYMass", &mumufit_DecayVtxX_bestYMass, &b_mumufit_DecayVtxX_bestYMass);
   fChain->SetBranchAddress("mumufit_DecayVtxY_bestYMass", &mumufit_DecayVtxY_bestYMass, &b_mumufit_DecayVtxY_bestYMass);
   fChain->SetBranchAddress("mumufit_DecayVtxZ_bestYMass", &mumufit_DecayVtxZ_bestYMass, &b_mumufit_DecayVtxZ_bestYMass);
   fChain->SetBranchAddress("mumufit_DecayVtxXE_bestYMass", &mumufit_DecayVtxXE_bestYMass, &b_mumufit_DecayVtxXE_bestYMass);
   fChain->SetBranchAddress("mumufit_DecayVtxYE_bestYMass", &mumufit_DecayVtxYE_bestYMass, &b_mumufit_DecayVtxYE_bestYMass);
   fChain->SetBranchAddress("mumufit_DecayVtxZE_bestYMass", &mumufit_DecayVtxZE_bestYMass, &b_mumufit_DecayVtxZE_bestYMass);
   fChain->SetBranchAddress("mumufit_p4_bestYMass", &mumufit_p4_bestYMass, &b_mumufit_p4_bestYMass);
   fChain->SetBranchAddress("bestVertex_and_bestYMass", &bestVertex_and_bestYMass, &b_bestVertex_and_bestYMass);
   fChain->SetBranchAddress("fourMuFit_Mass_allComb_mix_bestYMass", &fourMuFit_Mass_allComb_mix_bestYMass, &b_fourMuFit_Mass_allComb_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_Mass_mix_bestYMass", &fourMuFit_Mass_mix_bestYMass, &b_fourMuFit_Mass_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_MassErr_mix_bestYMass", &fourMuFit_MassErr_mix_bestYMass, &b_fourMuFit_MassErr_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_VtxX_mix_bestYMass", &fourMuFit_VtxX_mix_bestYMass, &b_fourMuFit_VtxX_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_VtxY_mix_bestYMass", &fourMuFit_VtxY_mix_bestYMass, &b_fourMuFit_VtxY_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_VtxZ_mix_bestYMass", &fourMuFit_VtxZ_mix_bestYMass, &b_fourMuFit_VtxZ_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_VtxProb_mix_bestYMass", &fourMuFit_VtxProb_mix_bestYMass, &b_fourMuFit_VtxProb_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_Chi2_mix_bestYMass", &fourMuFit_Chi2_mix_bestYMass, &b_fourMuFit_Chi2_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_ndof_mix_bestYMass", &fourMuFit_ndof_mix_bestYMass, &b_fourMuFit_ndof_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_3plus1_mix_bestYMass", &fourMuFit_3plus1_mix_bestYMass, &b_fourMuFit_3plus1_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_p4_mix_bestYMass", &fourMuFit_p4_mix_bestYMass, &b_fourMuFit_p4_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_mu1p4_mix_bestYMass", &fourMuFit_mu1p4_mix_bestYMass, &b_fourMuFit_mu1p4_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_mu2p4_mix_bestYMass", &fourMuFit_mu2p4_mix_bestYMass, &b_fourMuFit_mu2p4_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_mu3p4_mix_bestYMass", &fourMuFit_mu3p4_mix_bestYMass, &b_fourMuFit_mu3p4_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_mu4p4_mix_bestYMass", &fourMuFit_mu4p4_mix_bestYMass, &b_fourMuFit_mu4p4_mix_bestYMass);
   fChain->SetBranchAddress("mu3Charge_mix_bestYMass", &mu3Charge_mix_bestYMass, &b_mu3Charge_mix_bestYMass);
   fChain->SetBranchAddress("mu4Charge_mix_bestYMass", &mu4Charge_mix_bestYMass, &b_mu4Charge_mix_bestYMass);
   fChain->SetBranchAddress("mu3_p4_mix_bestYMass", &mu3_p4_mix_bestYMass, &b_mu3_p4_mix_bestYMass);
   fChain->SetBranchAddress("mu4_p4_mix_bestYMass", &mu4_p4_mix_bestYMass, &b_mu4_p4_mix_bestYMass);
   fChain->SetBranchAddress("mu3_d0_mix_bestYMass", &mu3_d0_mix_bestYMass, &b_mu3_d0_mix_bestYMass);
   fChain->SetBranchAddress("mu3_d0err_mix_bestYMass", &mu3_d0err_mix_bestYMass, &b_mu3_d0err_mix_bestYMass);
   fChain->SetBranchAddress("mu4_d0_mix_bestYMass", &mu4_d0_mix_bestYMass, &b_mu4_d0_mix_bestYMass);
   fChain->SetBranchAddress("mu4_d0err_mix_bestYMass", &mu4_d0err_mix_bestYMass, &b_mu4_d0err_mix_bestYMass);
   fChain->SetBranchAddress("mu3_dz_mix_bestYMass", &mu3_dz_mix_bestYMass, &b_mu3_dz_mix_bestYMass);
   fChain->SetBranchAddress("mu3_dzerr_mix_bestYMass", &mu3_dzerr_mix_bestYMass, &b_mu3_dzerr_mix_bestYMass);
   fChain->SetBranchAddress("mu4_dz_mix_bestYMass", &mu4_dz_mix_bestYMass, &b_mu4_dz_mix_bestYMass);
   fChain->SetBranchAddress("mu4_dzerr_mix_bestYMass", &mu4_dzerr_mix_bestYMass, &b_mu4_dzerr_mix_bestYMass);
   fChain->SetBranchAddress("fourMuFit_Mass_allComb_bestYMass", &fourMuFit_Mass_allComb_bestYMass, &b_fourMuFit_Mass_allComb_bestYMass);
   fChain->SetBranchAddress("fourMuFit_Mass_bestYMass", &fourMuFit_Mass_bestYMass, &b_fourMuFit_Mass_bestYMass);
   fChain->SetBranchAddress("fourMuFit_MassErr_bestYMass", &fourMuFit_MassErr_bestYMass, &b_fourMuFit_MassErr_bestYMass);
   fChain->SetBranchAddress("fourMuFit_VtxX_bestYMass", &fourMuFit_VtxX_bestYMass, &b_fourMuFit_VtxX_bestYMass);
   fChain->SetBranchAddress("fourMuFit_VtxY_bestYMass", &fourMuFit_VtxY_bestYMass, &b_fourMuFit_VtxY_bestYMass);
   fChain->SetBranchAddress("fourMuFit_VtxZ_bestYMass", &fourMuFit_VtxZ_bestYMass, &b_fourMuFit_VtxZ_bestYMass);
   fChain->SetBranchAddress("fourMuFit_VtxProb_bestYMass", &fourMuFit_VtxProb_bestYMass, &b_fourMuFit_VtxProb_bestYMass);
   fChain->SetBranchAddress("fourMuFit_Chi2_bestYMass", &fourMuFit_Chi2_bestYMass, &b_fourMuFit_Chi2_bestYMass);
   fChain->SetBranchAddress("fourMuFit_ndof_bestYMass", &fourMuFit_ndof_bestYMass, &b_fourMuFit_ndof_bestYMass);
   fChain->SetBranchAddress("fourMuFit_p4_bestYMass", &fourMuFit_p4_bestYMass, &b_fourMuFit_p4_bestYMass);
   fChain->SetBranchAddress("fourMuFit_mu1p4_bestYMass", &fourMuFit_mu1p4_bestYMass, &b_fourMuFit_mu1p4_bestYMass);
   fChain->SetBranchAddress("fourMuFit_mu2p4_bestYMass", &fourMuFit_mu2p4_bestYMass, &b_fourMuFit_mu2p4_bestYMass);
   fChain->SetBranchAddress("fourMuFit_mu3p4_bestYMass", &fourMuFit_mu3p4_bestYMass, &b_fourMuFit_mu3p4_bestYMass);
   fChain->SetBranchAddress("fourMuFit_mu4p4_bestYMass", &fourMuFit_mu4p4_bestYMass, &b_fourMuFit_mu4p4_bestYMass);
   fChain->SetBranchAddress("mu3Charge_bestYMass", &mu3Charge_bestYMass, &b_mu3Charge_bestYMass);
   fChain->SetBranchAddress("mu4Charge_bestYMass", &mu4Charge_bestYMass, &b_mu4Charge_bestYMass);
   fChain->SetBranchAddress("mu3_p4_bestYMass", &mu3_p4_bestYMass, &b_mu3_p4_bestYMass);
   fChain->SetBranchAddress("mu4_p4_bestYMass", &mu4_p4_bestYMass, &b_mu4_p4_bestYMass);
   fChain->SetBranchAddress("mu3_d0_bestYMass", &mu3_d0_bestYMass, &b_mu3_d0_bestYMass);
   fChain->SetBranchAddress("mu3_d0err_bestYMass", &mu3_d0err_bestYMass, &b_mu3_d0err_bestYMass);
   fChain->SetBranchAddress("mu4_d0_bestYMass", &mu4_d0_bestYMass, &b_mu4_d0_bestYMass);
   fChain->SetBranchAddress("mu4_d0err_bestYMass", &mu4_d0err_bestYMass, &b_mu4_d0err_bestYMass);
   fChain->SetBranchAddress("mu3_dz_bestYMass", &mu3_dz_bestYMass, &b_mu3_dz_bestYMass);
   fChain->SetBranchAddress("mu3_dzerr_bestYMass", &mu3_dzerr_bestYMass, &b_mu3_dzerr_bestYMass);
   fChain->SetBranchAddress("mu4_dz_bestYMass", &mu4_dz_bestYMass, &b_mu4_dz_bestYMass);
   fChain->SetBranchAddress("mu4_dzerr_bestYMass", &mu4_dzerr_bestYMass, &b_mu4_dzerr_bestYMass);
   fChain->SetBranchAddress("mu1_Tight_bestYMass", &mu1_Tight_bestYMass, &b_mu1_Tight_bestYMass);
   fChain->SetBranchAddress("mu2_Tight_bestYMass", &mu2_Tight_bestYMass, &b_mu2_Tight_bestYMass);
   fChain->SetBranchAddress("mu3_Tight_bestYMass", &mu3_Tight_bestYMass, &b_mu3_Tight_bestYMass);
   fChain->SetBranchAddress("mu4_Tight_bestYMass", &mu4_Tight_bestYMass, &b_mu4_Tight_bestYMass);
   fChain->SetBranchAddress("mu3_pdgID_bestYMass", &mu3_pdgID_bestYMass, &b_mu3_pdgID_bestYMass);
   fChain->SetBranchAddress("mu4_pdgID_bestYMass", &mu4_pdgID_bestYMass, &b_mu4_pdgID_bestYMass);
   Notify();
}

Bool_t readTree2018::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void readTree2018::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t readTree2018::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef readTree2018_cxx
