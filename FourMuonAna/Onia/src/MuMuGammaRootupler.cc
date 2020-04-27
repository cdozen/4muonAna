// -*- C++ -*-
//
// Package:    MuMuGammaRootupler
// Class:      MuMuGammaRootupler
// 
// Description: Dump mumugamma decays
//
// Author:  Zhen Hu
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/Common/interface/View.h>
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace muon;
using namespace trigger;

std::vector<std::vector<pat::Muon>> muons_previousEvent;
std::vector<pat::Muon> muons_previousEvent_bestYMass;
const ParticleMass muonMass(0.10565837);
float muonSigma = muonMass*1E-6;
//
// class declaration
//

class MuMuGammaRootupler:public edm::EDAnalyzer {
	public:
		explicit MuMuGammaRootupler(const edm::ParameterSet &);
		~MuMuGammaRootupler();

		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

	private:
		UInt_t getTriggerBits(const edm::Event &);
                bool TriggerMatch(bool TriggerPassed, pat::CompositeCandidate dimuonCand);
                bool TriggerMatch_restMuons(TLorentzVector mu3p4, TLorentzVector mu4p4);
                bool findTrigger(edm::Handle<edm::TriggerResults> &hltR,
                                 std::vector < std::string > & triggersToCheck,
                                 std::vector < std::string > & triggerNameFound); 
                void analyzeTrigger(edm::Handle<edm::TriggerResults> &hltR,
                                       edm::Handle<trigger::TriggerEvent> &hltE,
                                         const std::string& triggerName);
                bool triggerDecision(edm::Handle<edm::TriggerResults> &hltR, int iTrigger);
		bool   isAncestor(const reco::Candidate *, const reco::Candidate *);
		const  reco::Candidate* GetAncestor(const reco::Candidate *);
		int   tightMuon(edm::View<pat::Muon>::const_iterator rmu, reco::Vertex vertex);
		int   mediumMuon(edm::View<pat::Muon>::const_iterator rmu);
		int   looseMuon(edm::View<pat::Muon>::const_iterator rmu);
		void fillUpsilonVector(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs);
		void fillUpsilonBestVertex(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs);
		void fillUpsilonBestMass(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs);
		void fourMuonFit(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);
		int fourMuonMixFit(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, std::vector<pat::Muon> muons_previous, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);
		int fourMuonMixFit(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, std::vector<pat::Muon> muons_previous1, std::vector<pat::Muon> muons_previous2, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);      
		void fourMuonFit_bestYMass(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);
		int  fourMuonMixFit_bestYMass(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, std::vector<pat::Muon> muons_previous, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);

		virtual void beginJob();
		virtual void analyze(const edm::Event &, const edm::EventSetup &);
		virtual void endJob(const edm::Event &);

		virtual void beginRun(edm::Run &, edm::EventSetup const&);
		virtual void endRun(edm::Run const &, edm::EventSetup const &);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
		virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);

		// ----------member data ---------------------------
		std::string file_name;
                std::string triggersPassed;
                vector<string> mu1_filtersMatched; // for mu1, all filters it is matched to
                vector<string> mu2_filtersMatched; // for mu2, all filters it is matched to
		edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
		edm::EDGetTokenT<pat::CompositeCandidateCollection> conversion_Label;
		edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
		edm::EDGetTokenT<reco::BeamSpot> bs_Label;
		edm::EDGetTokenT<edm::View<pat::Muon>> muon_Label;
		edm::EDGetTokenT<edm::TriggerResults> triggerResultsTok_;
                edm::EDGetTokenT<trigger::TriggerEvent>triggerEventTok_;
//                edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
		int  pdgid_;
                vector<string> triggersToApply;
		std::vector<double> OniaMassCuts_;
		bool isMC_;
		bool OnlyBest_;
		bool OnlyGen_;
		double upsilon_mass_;
		uint32_t triggerCuts_;
		bool best4muonCand_;

		TH1F* myFourMuM_fit;
		TH1F* myFourMuVtxP_fit;

		TH1F* myDimuonMass_all;
 		TH1F* myY1SFitMass_all;
                int nEventsAnalyzed;
		UInt_t run;
		UInt_t lumi;
		UInt_t event;
		Int_t  irank;
		UInt_t trigger;
		Int_t numPrimaryVertices;
		Float_t pv_x;
		Float_t pv_y;
		Float_t pv_z;

		std::vector<Int_t> v_mu1Charge;
		std::vector<Int_t> v_mu2Charge;
		std::vector<Float_t> v_mu1_d0;
		std::vector<Float_t> v_mu2_d0;
		std::vector<Float_t> v_mu1_d0err;
		std::vector<Float_t> v_mu2_d0err;
		std::vector<Float_t> v_mu1_dz;
		std::vector<Float_t> v_mu2_dz;
		std::vector<Float_t> v_mu1_dzerr;
		std::vector<Float_t> v_mu2_dzerr;
		std::vector<Float_t> v_mu1_vz;
		std::vector<Float_t> v_mu2_vz;
		std::vector<Float_t> v_mumufit_Mass;
		std::vector<Float_t> v_mumufit_MassErr;
		std::vector<Float_t> v_mumufit_VtxCL;
		std::vector<Float_t> v_mumufit_VtxCL2;
		std::vector<Float_t> v_mumufit_DecayVtxX;
		std::vector<Float_t> v_mumufit_DecayVtxY;
		std::vector<Float_t> v_mumufit_DecayVtxZ;
		std::vector<Float_t> v_mumufit_DecayVtxXE;
		std::vector<Float_t> v_mumufit_DecayVtxYE;
		std::vector<Float_t> v_mumufit_DecayVtxZE;

		TLorentzVector dimuon_p4;
		TLorentzVector mu1_p4;
		TLorentzVector mu2_p4;
		Int_t mu1Charge;
		Int_t mu2Charge;
		Float_t mu1_d0;
		Float_t mu2_d0;
		Float_t mu1_d0err;
		Float_t mu2_d0err;
		Float_t mu1_dz;
		Float_t mu2_dz;
		Float_t mu1_dzerr;
		Float_t mu2_dzerr;
		Float_t mu1_vz;
		Float_t mu2_vz;
		Float_t mumufit_Mass;
		Float_t mumufit_MassErr;
		Float_t mumufit_VtxCL;
		Float_t mumufit_VtxCL2;
		Float_t mumufit_DecayVtxX;
		Float_t mumufit_DecayVtxY;
		Float_t mumufit_DecayVtxZ;
		Float_t mumufit_DecayVtxXE;
		Float_t mumufit_DecayVtxYE;
		Float_t mumufit_DecayVtxZE;
		TLorentzVector mumufit_p4;

		std::vector<Float_t> fourMuFit_Mass_allComb_mix;
		Float_t fourMuFit_Mass_mix;
		Float_t fourMuFit_MassErr_mix;
		Float_t fourMuFit_VtxX_mix;
		Float_t fourMuFit_VtxY_mix;
		Float_t fourMuFit_VtxZ_mix;
		Float_t fourMuFit_VtxProb_mix;
		Float_t fourMuFit_Chi2_mix;
		TLorentzVector fourMuFit_p4_mix;
		Int_t fourMuFit_ndof_mix;
		Int_t fourMuFit_3plus1_mix;
		Int_t mu3Charge_mix;
		Int_t mu4Charge_mix;
		Float_t mu3_d0_mix;
		Float_t mu4_d0_mix;
		Float_t mu3_d0err_mix;
		Float_t mu4_d0err_mix;
		Float_t mu3_dz_mix;
		Float_t mu4_dz_mix;
		Float_t mu3_dzerr_mix;
		Float_t mu4_dzerr_mix;
		TLorentzVector fourMuFit_mu1p4_mix;
		TLorentzVector fourMuFit_mu2p4_mix;
		TLorentzVector fourMuFit_mu3p4_mix;
		TLorentzVector fourMuFit_mu4p4_mix;
		TLorentzVector mu3_p4_mix;
		TLorentzVector mu4_p4_mix;
		Float_t fourMuFit_Mass_mix3evts;
		Float_t fourMuFit_VtxProb_mix3evts;
		TLorentzVector fourMuFit_p4_mix3evts;

		std::vector<Float_t> genbkg_mu1_Pt;
		std::vector<Float_t> genbkg_mu1_Eta;
		std::vector<Float_t> genbkg_mu1_Phi;
		std::vector<Float_t> genbkg_mu1_Mass;
		std::vector<Float_t> genbkg_mu2_Pt;
		std::vector<Float_t> genbkg_mu2_Eta;
		std::vector<Float_t> genbkg_mu2_Phi;
		std::vector<Float_t> genbkg_mu2_Mass;
		std::vector<Float_t> genbkg_mu3_Pt;
		std::vector<Float_t> genbkg_mu3_Eta;
		std::vector<Float_t> genbkg_mu3_Phi;
		std::vector<Float_t> genbkg_mu3_Mass;
		std::vector<Float_t> genbkg_mu4_Pt;
		std::vector<Float_t> genbkg_mu4_Eta;
		std::vector<Float_t> genbkg_mu4_Phi;
		std::vector<Float_t> genbkg_mu4_Mass;

		std::vector<Float_t> fourMuFit_Mass_allComb;
		std::vector<Float_t> fourMuFit_Mass;
		std::vector<Float_t> fourMuFit_MassErr;
		std::vector<Float_t> fourMuFit_Pt;
		std::vector<Float_t> fourMuFit_Eta;
		std::vector<Float_t> fourMuFit_Phi;
		std::vector<Float_t> fourMuFit_VtxX;
		std::vector<Float_t> fourMuFit_VtxY;
		std::vector<Float_t> fourMuFit_VtxZ;
		std::vector<Float_t> fourMuFit_VtxProb;
		std::vector<Float_t> fourMuFit_Chi2;
		std::vector<Int_t> fourMuFit_ndof;
		std::vector<Float_t> fourMuFit_mu1Pt;
		std::vector<Float_t> fourMuFit_mu1Eta;
		std::vector<Float_t> fourMuFit_mu1Phi;
		std::vector<Float_t> fourMuFit_mu1E;
		std::vector<Float_t> fourMuFit_mu2Pt;
		std::vector<Float_t> fourMuFit_mu2Eta;
		std::vector<Float_t> fourMuFit_mu2Phi;
		std::vector<Float_t> fourMuFit_mu2E;
		std::vector<Float_t> fourMuFit_mu3Pt;
		std::vector<Float_t> fourMuFit_mu3Eta;
		std::vector<Float_t> fourMuFit_mu3Phi;
		std::vector<Float_t> fourMuFit_mu3E;
		std::vector<Float_t> fourMuFit_mu4Pt;
		std::vector<Float_t> fourMuFit_mu4Eta;
		std::vector<Float_t> fourMuFit_mu4Phi;
		std::vector<Float_t> fourMuFit_mu4E;
		std::vector<Float_t> mu3_Pt;
		std::vector<Float_t> mu3_Eta;
		std::vector<Float_t> mu3_Phi;
		std::vector<Float_t> mu3_E;
		std::vector<Float_t> mu4_Pt;
		std::vector<Float_t> mu4_Eta;
		std::vector<Float_t> mu4_Phi;
		std::vector<Float_t> mu4_E;
		std::vector<Int_t> mu3Charge;
		std::vector<Int_t> mu4Charge;
		std::vector<Float_t> mu3_d0;
		std::vector<Float_t> mu4_d0;
		std::vector<Float_t> mu3_d0err;
		std::vector<Float_t> mu4_d0err;
		std::vector<Float_t> mu3_dz;
		std::vector<Float_t> mu4_dz;
		std::vector<Float_t> mu3_dzerr;
		std::vector<Float_t> mu4_dzerr;
		std::vector<Int_t> mu1_Tight;
		std::vector<Int_t> mu2_Tight;
		std::vector<Int_t> mu3_Tight;
		std::vector<Int_t> mu4_Tight;
		std::vector<Int_t> mu1_Medium;
		std::vector<Int_t> mu2_Medium;
		std::vector<Int_t> mu3_Medium;
		std::vector<Int_t> mu4_Medium;
		std::vector<Int_t> mu1_Loose;
		std::vector<Int_t> mu2_Loose;
		std::vector<Int_t> mu3_Loose;
		std::vector<Int_t> mu4_Loose;
		std::vector<Int_t> mu1_pdgID;
		std::vector<Int_t> mu2_pdgID;
		std::vector<Int_t> mu3_pdgID;
		std::vector<Int_t> mu4_pdgID;

		/*
			std::vector<Float_t> fourMuFit_Mass_allComb;
			Float_t fourMuFit_Mass;
			Float_t fourMuFit_MassErr;
			Float_t fourMuFit_VtxX;
			Float_t fourMuFit_VtxY;
			Float_t fourMuFit_VtxZ;
			Float_t fourMuFit_VtxProb;
			Float_t fourMuFit_Chi2;
			Int_t fourMuFit_ndof;
			TLorentzVector fourMuFit_p4;
			Int_t mu3Charge;
			Int_t mu4Charge;
			Float_t mu3_d0;
			Float_t mu4_d0;
			Float_t mu3_d0err;
			Float_t mu4_d0err;
			Float_t mu3_dz;
			Float_t mu4_dz;
			Float_t mu3_dzerr;
			Float_t mu4_dzerr;
			TLorentzVector fourMuFit_mu1p4;
			TLorentzVector fourMuFit_mu2p4;
			TLorentzVector fourMuFit_mu3p4;
			TLorentzVector fourMuFit_mu4p4;
			TLorentzVector mu3_p4;
			TLorentzVector mu4_p4;
			Int_t mu1_Tight;
			Int_t mu2_Tight;
			Int_t mu3_Tight;
			Int_t mu4_Tight;
			Int_t mu3_pdgID;
			Int_t mu4_pdgID;
			std::vector<pat::CompositeCandidate> upsilonMuons;
			std::vector<pat::Muon> theRestMuons;
			*/

		TLorentzVector dimuon_p4_bestYMass;
		TLorentzVector mu1_p4_bestYMass;
		TLorentzVector mu2_p4_bestYMass;
		Int_t mu1Charge_bestYMass;
		Int_t mu2Charge_bestYMass;
		Float_t mu1_d0_bestYMass;
		Float_t mu2_d0_bestYMass;
		Float_t mu1_d0err_bestYMass;
		Float_t mu2_d0err_bestYMass;
		Float_t mu1_dz_bestYMass;
		Float_t mu2_dz_bestYMass;
		Float_t mu1_dzerr_bestYMass;
		Float_t mu2_dzerr_bestYMass;
		Float_t mumufit_Mass_bestYMass;
		Float_t mumufit_MassErr_bestYMass;
		Float_t mumufit_VtxCL_bestYMass;
		Float_t mumufit_VtxCL2_bestYMass;
		Float_t mumufit_DecayVtxX_bestYMass;
		Float_t mumufit_DecayVtxY_bestYMass;
		Float_t mumufit_DecayVtxZ_bestYMass;
		Float_t mumufit_DecayVtxXE_bestYMass;
		Float_t mumufit_DecayVtxYE_bestYMass;
		Float_t mumufit_DecayVtxZE_bestYMass;
		TLorentzVector mumufit_p4_bestYMass;
		Int_t bestVertex_and_bestYMass;

		std::vector<Float_t> fourMuFit_Mass_allComb_mix_bestYMass;
		Float_t fourMuFit_Mass_mix_bestYMass;
		Float_t fourMuFit_MassErr_mix_bestYMass;
		Float_t fourMuFit_VtxX_mix_bestYMass;
		Float_t fourMuFit_VtxY_mix_bestYMass;
		Float_t fourMuFit_VtxZ_mix_bestYMass;
		Float_t fourMuFit_VtxProb_mix_bestYMass;
		Float_t fourMuFit_Chi2_mix_bestYMass;
		Int_t fourMuFit_ndof_mix_bestYMass;
		Int_t fourMuFit_3plus1_mix_bestYMass;
		TLorentzVector fourMuFit_p4_mix_bestYMass;
		Int_t mu3Charge_mix_bestYMass;
		Int_t mu4Charge_mix_bestYMass;
		Float_t mu3_d0_mix_bestYMass;
		Float_t mu4_d0_mix_bestYMass;
		Float_t mu3_d0err_mix_bestYMass;
		Float_t mu4_d0err_mix_bestYMass;
		Float_t mu3_dz_mix_bestYMass;
		Float_t mu4_dz_mix_bestYMass;
		Float_t mu3_dzerr_mix_bestYMass;
		Float_t mu4_dzerr_mix_bestYMass;
		TLorentzVector fourMuFit_mu1p4_mix_bestYMass;
		TLorentzVector fourMuFit_mu2p4_mix_bestYMass;
		TLorentzVector fourMuFit_mu3p4_mix_bestYMass;
		TLorentzVector fourMuFit_mu4p4_mix_bestYMass;
		TLorentzVector mu3_p4_mix_bestYMass;
		TLorentzVector mu4_p4_mix_bestYMass;

		std::vector<Float_t> fourMuFit_Mass_allComb_bestYMass;
		Float_t fourMuFit_Mass_bestYMass;
		Float_t fourMuFit_MassErr_bestYMass;
		Float_t fourMuFit_VtxX_bestYMass;
		Float_t fourMuFit_VtxY_bestYMass;
		Float_t fourMuFit_VtxZ_bestYMass;
		Float_t fourMuFit_VtxProb_bestYMass;
		Float_t fourMuFit_Chi2_bestYMass;
		Int_t fourMuFit_ndof_bestYMass;
		TLorentzVector fourMuFit_p4_bestYMass;
		Int_t mu3Charge_bestYMass;
		Int_t mu4Charge_bestYMass;
		Float_t mu3_d0_bestYMass;
		Float_t mu4_d0_bestYMass;
		Float_t mu3_d0err_bestYMass;
		Float_t mu4_d0err_bestYMass;
		Float_t mu3_dz_bestYMass;
		Float_t mu4_dz_bestYMass;
		Float_t mu3_dzerr_bestYMass;
		Float_t mu4_dzerr_bestYMass;
		TLorentzVector fourMuFit_mu1p4_bestYMass;
		TLorentzVector fourMuFit_mu2p4_bestYMass;
		TLorentzVector fourMuFit_mu3p4_bestYMass;
		TLorentzVector fourMuFit_mu4p4_bestYMass;
		TLorentzVector mu3_p4_bestYMass;
		TLorentzVector mu4_p4_bestYMass;
		Int_t mu1_Tight_bestYMass;
		Int_t mu2_Tight_bestYMass;
		Int_t mu3_Tight_bestYMass;
		Int_t mu4_Tight_bestYMass;
		Int_t mu3_pdgID_bestYMass;
		Int_t mu4_pdgID_bestYMass;
                std::string rootFileName;
                edm::InputTag triggerEventTag_;   
                std::string hltName_;
                std::string   triggerName_;
                std::string hlTriggerSummaryAOD_; 
                edm::TriggerNames triggerNames;
                bool verbose;
                float trg_Match_dR_cut, trg_Match_dP_cut; 
                std::vector<std::string> triggerList;
                bool checkTrigger;
                int runNumber;
  		std::vector < reco::MuonCollection::const_iterator > allL1TrigMuons;
  		std::vector < reco::MuonCollection::const_iterator > allL2TrigMuons;
  		std::vector < reco::MuonCollection::const_iterator > allL3TrigMuons;
   		std::vector<TLorentzVector> allTrigMuons;
                std::vector<TLorentzVector> allRestTrigMuons;
  		std::vector < GlobalVector > allMuL1TriggerVectors;
  		std::vector < GlobalVector > allMuL2TriggerVectors;
  		std::vector < GlobalVector > allMuL3TriggerVectors_lowEff;
  		std::vector < GlobalVector > allMuL3TriggerVectors_highEff;
  		std::vector < GlobalVector > allMuHLTTriggerVectors;
                int lastTriggerModule;
                HLTConfigProvider hltConfig_;
		TTree *onia_tree;
		TTree *gen_tree;

		Int_t mother_pdgId;
		Int_t dimuon_pdgId;
		TLorentzVector gen_dimuon_p4;
		TLorentzVector gen_mu1_p4;
		TLorentzVector gen_mu2_p4;
		edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
		edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

//
// constructors and destructor
//
MuMuGammaRootupler::MuMuGammaRootupler(const edm::ParameterSet & iConfig):
	dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuons"))),
	conversion_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("conversions"))),
	primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
	bs_Label(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBeamSpot"))),
	muon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter< edm::InputTag>("muons"))),
	triggerResultsTok_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
        triggerEventTok_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerSummaryAOD"))),
	pdgid_(iConfig.getParameter<uint32_t>("onia_pdgid")),
	OniaMassCuts_(iConfig.getParameter<std::vector<double>>("onia_mass_cuts")),
	isMC_(iConfig.getParameter<bool>("isMC")),
	OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
	OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
	upsilon_mass_(iConfig.getParameter<double>("upsilon_mass")),
	triggerCuts_(iConfig.getParameter<uint32_t>("triggerCuts")),
	best4muonCand_(iConfig.getParameter<bool>("best4muonCand")),
        verbose(iConfig.getUntrackedParameter<bool>("VERBOSE",false)),
        trg_Match_dR_cut(iConfig.getUntrackedParameter<double>("TRG_Match_DR",0.2)),
        trg_Match_dP_cut(iConfig.getUntrackedParameter<double>("TRG_Match_DP",0.3)),
        triggerList(iConfig.getUntrackedParameter<std::vector<std::string>>("triggerList"))
{
	edm::Service < TFileService > fs;
	onia_tree = fs->make < TTree > ("oniaTree", "Tree of MuMuGamma");
	gen_tree = fs->make < TTree > ("genTree", "Tree of genCand");
        nEventsAnalyzed = 0;
        runNumber = -99;
        checkTrigger =true;
        hltName_ = "HLT";
        triggerName_ = "@"; // "@" means: analyze all triggers in config
	if (!OnlyGen_) {
		onia_tree->Branch("run",     &run,     "run/I");
		onia_tree->Branch("lumi",     &lumi,     "lumi/I");
		onia_tree->Branch("event",   &event,   "event/I");
		onia_tree->Branch("irank",   &irank,   "irank/I");
		onia_tree->Branch("trigger", &trigger, "trigger/I");
		onia_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
		onia_tree->Branch("pv_x",     &pv_x,     "pv_x/F");
		onia_tree->Branch("pv_y",     &pv_y,     "pv_y/F");
		onia_tree->Branch("pv_z",     &pv_z,     "pv_z/F");

		onia_tree->Branch("v_mu1Charge",   &v_mu1Charge);
		onia_tree->Branch("v_mu2Charge",   &v_mu2Charge);
		onia_tree->Branch("v_mu1_d0",   &v_mu1_d0);
		onia_tree->Branch("v_mu1_d0err",   &v_mu1_d0err);
		onia_tree->Branch("v_mu2_d0",   &v_mu2_d0);
		onia_tree->Branch("v_mu2_d0err",   &v_mu2_d0err);
		onia_tree->Branch("v_mu1_dz",   &v_mu1_dz);
		onia_tree->Branch("v_mu1_dzerr",   &v_mu1_dzerr);
		onia_tree->Branch("v_mu2_dz",   &v_mu2_dz);
		onia_tree->Branch("v_mu2_dzerr",   &v_mu2_dzerr);
		onia_tree->Branch("v_mu1_vz",   &v_mu1_vz);
		onia_tree->Branch("v_mu2_vz",   &v_mu2_vz);
		onia_tree->Branch("v_mumufit_Mass",&v_mumufit_Mass);
		onia_tree->Branch("v_mumufit_MassErr",&v_mumufit_MassErr);
		onia_tree->Branch("v_mumufit_VtxCL",&v_mumufit_VtxCL);
		onia_tree->Branch("v_mumufit_VtxCL2",&v_mumufit_VtxCL2);
		onia_tree->Branch("v_mumufit_DecayVtxX",&v_mumufit_DecayVtxX);
		onia_tree->Branch("v_mumufit_DecayVtxY",&v_mumufit_DecayVtxY);
		onia_tree->Branch("v_mumufit_DecayVtxZ",&v_mumufit_DecayVtxZ);
		onia_tree->Branch("v_mumufit_DecayVtxXE",&v_mumufit_DecayVtxXE);
		onia_tree->Branch("v_mumufit_DecayVtxYE",&v_mumufit_DecayVtxYE);
		onia_tree->Branch("v_mumufit_DecayVtxZE",&v_mumufit_DecayVtxZE);
                onia_tree->Branch("mu1_filtersMatched", & mu1_filtersMatched);
                onia_tree->Branch("mu2_filtersMatched", & mu2_filtersMatched);
		onia_tree->Branch("mu1_p4",  "TLorentzVector", &mu1_p4);
		onia_tree->Branch("mu2_p4",  "TLorentzVector", &mu2_p4);
		onia_tree->Branch("mu1Charge",   &mu1Charge,    "mu1Charge/I");
		onia_tree->Branch("mu2Charge",   &mu2Charge,    "mu2Charge/I");
		onia_tree->Branch("mu1_d0",   &mu1_d0,    "mu1_d0/F");
		onia_tree->Branch("mu1_d0err",   &mu1_d0err,    "mu1_d0err/F");
		onia_tree->Branch("mu2_d0",   &mu2_d0,    "mu2_d0/F");
		onia_tree->Branch("mu2_d0err",   &mu2_d0err,    "mu2_d0err/F");
		onia_tree->Branch("mu1_dz",   &mu1_dz,    "mu1_dz/F");
		onia_tree->Branch("mu1_dzerr",   &mu1_dzerr,    "mu1_dzerr/F");
		onia_tree->Branch("mu2_dz",   &mu2_dz,    "mu2_dz/F");
		onia_tree->Branch("mu2_dzerr",   &mu2_dzerr,    "mu2_dzerr/F");
		onia_tree->Branch("mu1_vz",   &mu1_vz,    "mu1_vz/F");
		onia_tree->Branch("mu2_vz",   &mu2_vz,    "mu2_vz/F");
		onia_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
		onia_tree->Branch("mumufit_Mass",&mumufit_Mass,"mumufit_Mass/F");
		onia_tree->Branch("mumufit_MassErr",&mumufit_MassErr,"mumufit_MassErr/F");
		onia_tree->Branch("mumufit_VtxCL",&mumufit_VtxCL,"mumufit_VtxCL/F");
		onia_tree->Branch("mumufit_VtxCL2",&mumufit_VtxCL2,"mumufit_VtxCL2/F");
		onia_tree->Branch("mumufit_DecayVtxX",&mumufit_DecayVtxX,"mumufit_DecayVtxX/F");
		onia_tree->Branch("mumufit_DecayVtxY",&mumufit_DecayVtxY,"mumufit_DecayVtxY/F");
		onia_tree->Branch("mumufit_DecayVtxZ",&mumufit_DecayVtxZ,"mumufit_DecayVtxZ/F");
		onia_tree->Branch("mumufit_DecayVtxXE",&mumufit_DecayVtxXE,"mumufit_DecayVtxXE/F");
		onia_tree->Branch("mumufit_DecayVtxYE",&mumufit_DecayVtxYE,"mumufit_DecayVtxYE/F");
		onia_tree->Branch("mumufit_DecayVtxZE",&mumufit_DecayVtxZE,"mumufit_DecayVtxZE/F");
		onia_tree->Branch("mumufit_p4",  "TLorentzVector", &mumufit_p4);

		onia_tree->Branch("fourMuFit_Mass_allComb_mix",&fourMuFit_Mass_allComb_mix);
		onia_tree->Branch("fourMuFit_Mass_mix",&fourMuFit_Mass_mix,"fourMuFit_Mass_mix/F");
		onia_tree->Branch("fourMuFit_MassErr_mix",&fourMuFit_MassErr_mix,"fourMuFit_MassErr_mix/F");
		onia_tree->Branch("fourMuFit_VtxX_mix",&fourMuFit_VtxX_mix,"fourMuFit_VtxX_mix/F");
		onia_tree->Branch("fourMuFit_VtxY_mix",&fourMuFit_VtxY_mix,"fourMuFit_VtxY_mix/F");
		onia_tree->Branch("fourMuFit_VtxZ_mix",&fourMuFit_VtxZ_mix,"fourMuFit_VtxZ_mix/F");
		onia_tree->Branch("fourMuFit_VtxProb_mix",&fourMuFit_VtxProb_mix,"fourMuFit_VtxProb_mix/F");
		onia_tree->Branch("fourMuFit_Chi2_mix",&fourMuFit_Chi2_mix,"fourMuFit_Chi2_mix/F");
		onia_tree->Branch("fourMuFit_ndof_mix",&fourMuFit_ndof_mix,"fourMuFit_ndof_mix/I");
		onia_tree->Branch("fourMuFit_3plus1_mix",&fourMuFit_3plus1_mix,"fourMuFit_3plus1_mix/I");
		onia_tree->Branch("fourMuFit_p4_mix",  "TLorentzVector", &fourMuFit_p4_mix);
		onia_tree->Branch("fourMuFit_mu1p4_mix",  "TLorentzVector", &fourMuFit_mu1p4_mix);
		onia_tree->Branch("fourMuFit_mu2p4_mix",  "TLorentzVector", &fourMuFit_mu2p4_mix);
		onia_tree->Branch("fourMuFit_mu3p4_mix",  "TLorentzVector", &fourMuFit_mu3p4_mix);
		onia_tree->Branch("fourMuFit_mu4p4_mix",  "TLorentzVector", &fourMuFit_mu4p4_mix);
		onia_tree->Branch("mu3Charge_mix",   &mu3Charge_mix,    "mu3Charge_mix/I");
		onia_tree->Branch("mu4Charge_mix",   &mu4Charge_mix,    "mu4Charge_mix/I");
		onia_tree->Branch("mu3_p4_mix",  "TLorentzVector", &mu3_p4_mix);
		onia_tree->Branch("mu4_p4_mix",  "TLorentzVector", &mu4_p4_mix);
		onia_tree->Branch("mu3_d0_mix",   &mu3_d0_mix,    "mu3_d0_mix/F");
		onia_tree->Branch("mu3_d0err_mix",   &mu3_d0err_mix,    "mu3_d0err_mix/F");
		onia_tree->Branch("mu4_d0_mix",   &mu4_d0_mix,    "mu4_d0_mix/F");
		onia_tree->Branch("mu4_d0err_mix",   &mu4_d0err_mix,    "mu4_d0err_mix/F");
		onia_tree->Branch("mu3_dz_mix",   &mu3_dz_mix,    "mu3_dz_mix/F");
		onia_tree->Branch("mu3_dzerr_mix",   &mu3_dzerr_mix,    "mu3_dzerr_mix/F");
		onia_tree->Branch("mu4_dz_mix",   &mu4_dz_mix,    "mu4_dz_mix/F");
		onia_tree->Branch("mu4_dzerr_mix",   &mu4_dzerr_mix,    "mu4_dzerr_mix/F");
		onia_tree->Branch("fourMuFit_Mass_mix3evts",&fourMuFit_Mass_mix3evts,"fourMuFit_Mass_mix3evts/F");
		onia_tree->Branch("fourMuFit_VtxProb_mix3evts",&fourMuFit_VtxProb_mix3evts,"fourMuFit_VtxProb_mix3evts/F");
		onia_tree->Branch("fourMuFit_p4_mix3evts",  "TLorentzVector", &fourMuFit_p4_mix3evts);

		gen_tree->Branch("genbkg_mu1_Pt",&genbkg_mu1_Pt);
		gen_tree->Branch("genbkg_mu1_Eta",&genbkg_mu1_Eta);
		gen_tree->Branch("genbkg_mu1_Phi",&genbkg_mu1_Phi);
		gen_tree->Branch("genbkg_mu1_Mass",&genbkg_mu1_Mass);
		gen_tree->Branch("genbkg_mu2_Pt",&genbkg_mu2_Pt);
		gen_tree->Branch("genbkg_mu2_Eta",&genbkg_mu2_Eta);
		gen_tree->Branch("genbkg_mu2_Phi",&genbkg_mu2_Phi);
		gen_tree->Branch("genbkg_mu2_Mass",&genbkg_mu2_Mass);
		gen_tree->Branch("genbkg_mu3_Pt",&genbkg_mu3_Pt);
		gen_tree->Branch("genbkg_mu3_Eta",&genbkg_mu3_Eta);
		gen_tree->Branch("genbkg_mu3_Phi",&genbkg_mu3_Phi);
		gen_tree->Branch("genbkg_mu3_Mass",&genbkg_mu3_Mass);
		gen_tree->Branch("genbkg_mu4_Pt",&genbkg_mu4_Pt);
		gen_tree->Branch("genbkg_mu4_Eta",&genbkg_mu4_Eta);
		gen_tree->Branch("genbkg_mu4_Phi",&genbkg_mu4_Phi);
		gen_tree->Branch("genbkg_mu4_Mass",&genbkg_mu4_Mass);

		onia_tree->Branch("fourMuFit_Mass_allComb",&fourMuFit_Mass_allComb);
		onia_tree->Branch("fourMuFit_Mass",&fourMuFit_Mass);
		onia_tree->Branch("fourMuFit_MassErr",&fourMuFit_MassErr);
		onia_tree->Branch("fourMuFit_Pt",&fourMuFit_Pt);
		onia_tree->Branch("fourMuFit_Eta",&fourMuFit_Eta);
		onia_tree->Branch("fourMuFit_Phi",&fourMuFit_Phi);
		onia_tree->Branch("fourMuFit_VtxX",&fourMuFit_VtxX);
		onia_tree->Branch("fourMuFit_VtxY",&fourMuFit_VtxY);
		onia_tree->Branch("fourMuFit_VtxZ",&fourMuFit_VtxZ);
		onia_tree->Branch("fourMuFit_VtxProb",&fourMuFit_VtxProb);
		onia_tree->Branch("fourMuFit_Chi2",&fourMuFit_Chi2);
		onia_tree->Branch("fourMuFit_ndof",&fourMuFit_ndof);
		onia_tree->Branch("fourMuFit_mu1Pt",&fourMuFit_mu1Pt);
		onia_tree->Branch("fourMuFit_mu1Eta",&fourMuFit_mu1Eta);
		onia_tree->Branch("fourMuFit_mu1Phi",&fourMuFit_mu1Phi);
		onia_tree->Branch("fourMuFit_mu1E",&fourMuFit_mu1E);
		onia_tree->Branch("fourMuFit_mu2Pt",&fourMuFit_mu2Pt);
		onia_tree->Branch("fourMuFit_mu2Eta",&fourMuFit_mu2Eta);
		onia_tree->Branch("fourMuFit_mu2Phi",&fourMuFit_mu2Phi);
		onia_tree->Branch("fourMuFit_mu2E",&fourMuFit_mu2E);
		onia_tree->Branch("fourMuFit_mu3Pt",&fourMuFit_mu3Pt);
		onia_tree->Branch("fourMuFit_mu3Eta",&fourMuFit_mu3Eta);
		onia_tree->Branch("fourMuFit_mu3Phi",&fourMuFit_mu3Phi);
		onia_tree->Branch("fourMuFit_mu3E",&fourMuFit_mu3E);
		onia_tree->Branch("fourMuFit_mu4Pt",&fourMuFit_mu4Pt);
		onia_tree->Branch("fourMuFit_mu4Eta",&fourMuFit_mu4Eta);
		onia_tree->Branch("fourMuFit_mu4Phi",&fourMuFit_mu4Phi);
		onia_tree->Branch("fourMuFit_mu4E",&fourMuFit_mu4E);
		onia_tree->Branch("mu3_Pt",   &mu3_Pt);
		onia_tree->Branch("mu3_Eta",   &mu3_Eta);
		onia_tree->Branch("mu3_Phi",   &mu3_Phi);
		onia_tree->Branch("mu3_E",   &mu3_E);
		onia_tree->Branch("mu4_Pt",   &mu4_Pt);
		onia_tree->Branch("mu4_Eta",   &mu4_Eta);
		onia_tree->Branch("mu4_Phi",   &mu4_Phi);
		onia_tree->Branch("mu4_E",   &mu4_E);
		onia_tree->Branch("mu3Charge",   &mu3Charge);
		onia_tree->Branch("mu4Charge",   &mu4Charge);
		onia_tree->Branch("mu3_d0",   &mu3_d0);
		onia_tree->Branch("mu3_d0err",   &mu3_d0err);
		onia_tree->Branch("mu4_d0",   &mu4_d0);
		onia_tree->Branch("mu4_d0err",   &mu4_d0err);
		onia_tree->Branch("mu3_dz",   &mu3_dz);
		onia_tree->Branch("mu3_dzerr",   &mu3_dzerr);
		onia_tree->Branch("mu4_dz",   &mu4_dz);
		onia_tree->Branch("mu4_dzerr",   &mu4_dzerr);
		onia_tree->Branch("mu1_Tight",   &mu1_Tight);
		onia_tree->Branch("mu2_Tight",   &mu2_Tight);
		onia_tree->Branch("mu3_Tight",   &mu3_Tight);
		onia_tree->Branch("mu4_Tight",   &mu4_Tight);
		onia_tree->Branch("mu1_Medium",   &mu1_Medium);
		onia_tree->Branch("mu2_Medium",   &mu2_Medium);
		onia_tree->Branch("mu3_Medium",   &mu3_Medium);
		onia_tree->Branch("mu4_Medium",   &mu4_Medium);
		onia_tree->Branch("mu1_Loose",   &mu1_Loose);
		onia_tree->Branch("mu2_Loose",   &mu2_Loose);
		onia_tree->Branch("mu3_Loose",   &mu3_Loose);
		onia_tree->Branch("mu4_Loose",   &mu4_Loose);
		onia_tree->Branch("mu1_pdgID",   &mu1_pdgID);
		onia_tree->Branch("mu2_pdgID",   &mu2_pdgID);
		onia_tree->Branch("mu3_pdgID",   &mu3_pdgID);
		onia_tree->Branch("mu4_pdgID",   &mu4_pdgID);

		/*
			onia_tree->Branch("fourMuFit_Mass_allComb",&fourMuFit_Mass_allComb);
			onia_tree->Branch("fourMuFit_Mass",&fourMuFit_Mass,"fourMuFit_Mass/F");
			onia_tree->Branch("fourMuFit_MassErr",&fourMuFit_MassErr,"fourMuFit_MassErr/F");
			onia_tree->Branch("fourMuFit_VtxX",&fourMuFit_VtxX,"fourMuFit_VtxX/F");
			onia_tree->Branch("fourMuFit_VtxY",&fourMuFit_VtxY,"fourMuFit_VtxY/F");
			onia_tree->Branch("fourMuFit_VtxZ",&fourMuFit_VtxZ,"fourMuFit_VtxZ/F");
			onia_tree->Branch("fourMuFit_VtxProb",&fourMuFit_VtxProb,"fourMuFit_VtxProb/F");
			onia_tree->Branch("fourMuFit_Chi2",&fourMuFit_Chi2,"fourMuFit_Chi2/F");
			onia_tree->Branch("fourMuFit_ndof",&fourMuFit_ndof,"fourMuFit_ndof/I");
			onia_tree->Branch("fourMuFit_p4",  "TLorentzVector", &fourMuFit_p4);
			onia_tree->Branch("fourMuFit_mu1p4",  "TLorentzVector", &fourMuFit_mu1p4);
			onia_tree->Branch("fourMuFit_mu2p4",  "TLorentzVector", &fourMuFit_mu2p4);
			onia_tree->Branch("fourMuFit_mu3p4",  "TLorentzVector", &fourMuFit_mu3p4);
			onia_tree->Branch("fourMuFit_mu4p4",  "TLorentzVector", &fourMuFit_mu4p4);
			onia_tree->Branch("mu3Charge",   &mu3Charge,    "mu3Charge/I");
			onia_tree->Branch("mu4Charge",   &mu4Charge,    "mu4Charge/I");
			onia_tree->Branch("mu3_p4",  "TLorentzVector", &mu3_p4);
			onia_tree->Branch("mu4_p4",  "TLorentzVector", &mu4_p4);
			onia_tree->Branch("mu3_d0",   &mu3_d0,    "mu3_d0/F");
			onia_tree->Branch("mu3_d0err",   &mu3_d0err,    "mu3_d0err/F");
			onia_tree->Branch("mu4_d0",   &mu4_d0,    "mu4_d0/F");
			onia_tree->Branch("mu4_d0err",   &mu4_d0err,    "mu4_d0err/F");
			onia_tree->Branch("mu3_dz",   &mu3_dz,    "mu3_dz/F");
			onia_tree->Branch("mu3_dzerr",   &mu3_dzerr,    "mu3_dzerr/F");
			onia_tree->Branch("mu4_dz",   &mu4_dz,    "mu4_dz/F");
			onia_tree->Branch("mu4_dzerr",   &mu4_dzerr,    "mu4_dzerr/F");
			onia_tree->Branch("mu1_Tight",   &mu1_Tight,    "mu1_Tight/I");
			onia_tree->Branch("mu2_Tight",   &mu2_Tight,    "mu2_Tight/I");
			onia_tree->Branch("mu3_Tight",   &mu3_Tight,    "mu3_Tight/I");
			onia_tree->Branch("mu4_Tight",   &mu4_Tight,    "mu4_Tight/I");
			onia_tree->Branch("mu3_pdgID",   &mu3_pdgID,    "mu3_pdgID/I");
			onia_tree->Branch("mu4_pdgID",   &mu4_pdgID,    "mu4_pdgID/I");
			onia_tree->Branch("upsilonMuons", &upsilonMuons);	
			onia_tree->Branch("theRestMuons", &theRestMuons); 
			*/
		onia_tree->Branch("mu1_p4_bestYMass",  "TLorentzVector", &mu1_p4_bestYMass);
		onia_tree->Branch("mu2_p4_bestYMass",  "TLorentzVector", &mu2_p4_bestYMass);
		onia_tree->Branch("mu1Charge_bestYMass",   &mu1Charge_bestYMass,    "mu1Charge_bestYMass/I");
		onia_tree->Branch("mu2Charge_bestYMass",   &mu2Charge_bestYMass,    "mu2Charge_bestYMass/I");
		onia_tree->Branch("mu1_d0_bestYMass",   &mu1_d0_bestYMass,    "mu1_d0_bestYMass/F");
		onia_tree->Branch("mu1_d0err_bestYMass",   &mu1_d0err_bestYMass,    "mu1_d0err_bestYMass/F");
		onia_tree->Branch("mu2_d0_bestYMass",   &mu2_d0_bestYMass,    "mu2_d0_bestYMass/F");
		onia_tree->Branch("mu2_d0err_bestYMass",   &mu2_d0err_bestYMass,    "mu2_d0err_bestYMass/F");
		onia_tree->Branch("mu1_dz_bestYMass",   &mu1_dz_bestYMass,    "mu1_dz_bestYMass/F");
		onia_tree->Branch("mu1_dzerr_bestYMass",   &mu1_dzerr_bestYMass,    "mu1_dzerr_bestYMass/F");
		onia_tree->Branch("mu2_dz_bestYMass",   &mu2_dz_bestYMass,    "mu2_dz_bestYMass/F");
		onia_tree->Branch("mu2_dzerr_bestYMass",   &mu2_dzerr_bestYMass,    "mu2_dzerr_bestYMass/F");
		onia_tree->Branch("dimuon_p4_bestYMass", "TLorentzVector", &dimuon_p4_bestYMass);
		onia_tree->Branch("mumufit_Mass_bestYMass",&mumufit_Mass_bestYMass,"mumufit_Mass_bestYMass/F");
		onia_tree->Branch("mumufit_MassErr_bestYMass",&mumufit_MassErr_bestYMass,"mumufit_MassErr_bestYMass/F");
		onia_tree->Branch("mumufit_VtxCL_bestYMass",&mumufit_VtxCL_bestYMass,"mumufit_VtxCL_bestYMass/F");
		onia_tree->Branch("mumufit_VtxCL2_bestYMass",&mumufit_VtxCL2_bestYMass,"mumufit_VtxCL2_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxX_bestYMass",&mumufit_DecayVtxX_bestYMass,"mumufit_DecayVtxX_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxY_bestYMass",&mumufit_DecayVtxY_bestYMass,"mumufit_DecayVtxY_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxZ_bestYMass",&mumufit_DecayVtxZ_bestYMass,"mumufit_DecayVtxZ_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxXE_bestYMass",&mumufit_DecayVtxXE_bestYMass,"mumufit_DecayVtxXE_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxYE_bestYMass",&mumufit_DecayVtxYE_bestYMass,"mumufit_DecayVtxYE_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxZE_bestYMass",&mumufit_DecayVtxZE_bestYMass,"mumufit_DecayVtxZE_bestYMass/F");
		onia_tree->Branch("mumufit_p4_bestYMass",  "TLorentzVector", &mumufit_p4_bestYMass);
		onia_tree->Branch("bestVertex_and_bestYMass", &bestVertex_and_bestYMass,"bestVertex_and_bestYMass/I");

		onia_tree->Branch("fourMuFit_Mass_allComb_mix_bestYMass",&fourMuFit_Mass_allComb_mix_bestYMass);
		onia_tree->Branch("fourMuFit_Mass_mix_bestYMass",&fourMuFit_Mass_mix_bestYMass,"fourMuFit_Mass_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_MassErr_mix_bestYMass",&fourMuFit_MassErr_mix_bestYMass,"fourMuFit_MassErr_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxX_mix_bestYMass",&fourMuFit_VtxX_mix_bestYMass,"fourMuFit_VtxX_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxY_mix_bestYMass",&fourMuFit_VtxY_mix_bestYMass,"fourMuFit_VtxY_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxZ_mix_bestYMass",&fourMuFit_VtxZ_mix_bestYMass,"fourMuFit_VtxZ_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxProb_mix_bestYMass",&fourMuFit_VtxProb_mix_bestYMass,"fourMuFit_VtxProb_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_Chi2_mix_bestYMass",&fourMuFit_Chi2_mix_bestYMass,"fourMuFit_Chi2_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_ndof_mix_bestYMass",&fourMuFit_ndof_mix_bestYMass,"fourMuFit_ndof_mix_bestYMass/I");
		onia_tree->Branch("fourMuFit_3plus1_mix_bestYMass",&fourMuFit_3plus1_mix_bestYMass,"fourMuFit_3plus1_mix_bestYMass/I");
		onia_tree->Branch("fourMuFit_p4_mix_bestYMass",  "TLorentzVector", &fourMuFit_p4_mix_bestYMass);
		onia_tree->Branch("fourMuFit_mu1p4_mix_bestYMass",  "TLorentzVector", &fourMuFit_mu1p4_mix_bestYMass);
		onia_tree->Branch("fourMuFit_mu2p4_mix_bestYMass",  "TLorentzVector", &fourMuFit_mu2p4_mix_bestYMass);
		onia_tree->Branch("fourMuFit_mu3p4_mix_bestYMass",  "TLorentzVector", &fourMuFit_mu3p4_mix_bestYMass);
		onia_tree->Branch("fourMuFit_mu4p4_mix_bestYMass",  "TLorentzVector", &fourMuFit_mu4p4_mix_bestYMass);
		onia_tree->Branch("mu3Charge_mix_bestYMass",   &mu3Charge_mix_bestYMass,    "mu3Charge_mix_bestYMass/I");
		onia_tree->Branch("mu4Charge_mix_bestYMass",   &mu4Charge_mix_bestYMass,    "mu4Charge_mix_bestYMass/I");
		onia_tree->Branch("mu3_p4_mix_bestYMass",  "TLorentzVector", &mu3_p4_mix_bestYMass);
		onia_tree->Branch("mu4_p4_mix_bestYMass",  "TLorentzVector", &mu4_p4_mix_bestYMass);
		onia_tree->Branch("mu3_d0_mix_bestYMass",   &mu3_d0_mix_bestYMass,    "mu3_d0_mix_bestYMass/F");
		onia_tree->Branch("mu3_d0err_mix_bestYMass",   &mu3_d0err_mix_bestYMass,    "mu3_d0err_mix_bestYMass/F");
		onia_tree->Branch("mu4_d0_mix_bestYMass",   &mu4_d0_mix_bestYMass,    "mu4_d0_mix_bestYMass/F");
		onia_tree->Branch("mu4_d0err_mix_bestYMass",   &mu4_d0err_mix_bestYMass,    "mu4_d0err_mix_bestYMass/F");
		onia_tree->Branch("mu3_dz_mix_bestYMass",   &mu3_dz_mix_bestYMass,    "mu3_dz_mix_bestYMass/F");
		onia_tree->Branch("mu3_dzerr_mix_bestYMass",   &mu3_dzerr_mix_bestYMass,    "mu3_dzerr_mix_bestYMass/F");
		onia_tree->Branch("mu4_dz_mix_bestYMass",   &mu4_dz_mix_bestYMass,    "mu4_dz_mix_bestYMass/F");
		onia_tree->Branch("mu4_dzerr_mix_bestYMass",   &mu4_dzerr_mix_bestYMass,    "mu4_dzerr_mix_bestYMass/F");

		onia_tree->Branch("fourMuFit_Mass_allComb_bestYMass",&fourMuFit_Mass_allComb_bestYMass);
		onia_tree->Branch("fourMuFit_Mass_bestYMass",&fourMuFit_Mass_bestYMass,"fourMuFit_Mass_bestYMass/F");
		onia_tree->Branch("fourMuFit_MassErr_bestYMass",&fourMuFit_MassErr_bestYMass,"fourMuFit_MassErr_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxX_bestYMass",&fourMuFit_VtxX_bestYMass,"fourMuFit_VtxX_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxY_bestYMass",&fourMuFit_VtxY_bestYMass,"fourMuFit_VtxY_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxZ_bestYMass",&fourMuFit_VtxZ_bestYMass,"fourMuFit_VtxZ_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxProb_bestYMass",&fourMuFit_VtxProb_bestYMass,"fourMuFit_VtxProb_bestYMass/F");
		onia_tree->Branch("fourMuFit_Chi2_bestYMass",&fourMuFit_Chi2_bestYMass,"fourMuFit_Chi2_bestYMass/F");
		onia_tree->Branch("fourMuFit_ndof_bestYMass",&fourMuFit_ndof_bestYMass,"fourMuFit_ndof_bestYMass/I");
		onia_tree->Branch("fourMuFit_p4_bestYMass",  "TLorentzVector", &fourMuFit_p4_bestYMass);
		onia_tree->Branch("fourMuFit_mu1p4_bestYMass",  "TLorentzVector", &fourMuFit_mu1p4_bestYMass);
		onia_tree->Branch("fourMuFit_mu2p4_bestYMass",  "TLorentzVector", &fourMuFit_mu2p4_bestYMass);
		onia_tree->Branch("fourMuFit_mu3p4_bestYMass",  "TLorentzVector", &fourMuFit_mu3p4_bestYMass);
		onia_tree->Branch("fourMuFit_mu4p4_bestYMass",  "TLorentzVector", &fourMuFit_mu4p4_bestYMass);
		onia_tree->Branch("mu3Charge_bestYMass",   &mu3Charge_bestYMass,    "mu3Charge_bestYMass/I");
		onia_tree->Branch("mu4Charge_bestYMass",   &mu4Charge_bestYMass,    "mu4Charge_bestYMass/I");
		onia_tree->Branch("mu3_p4_bestYMass",  "TLorentzVector", &mu3_p4_bestYMass);
		onia_tree->Branch("mu4_p4_bestYMass",  "TLorentzVector", &mu4_p4_bestYMass);
		onia_tree->Branch("mu3_d0_bestYMass",   &mu3_d0_bestYMass,    "mu3_d0_bestYMass/F");
		onia_tree->Branch("mu3_d0err_bestYMass",   &mu3_d0err_bestYMass,    "mu3_d0err_bestYMass/F");
		onia_tree->Branch("mu4_d0_bestYMass",   &mu4_d0_bestYMass,    "mu4_d0_bestYMass/F");
		onia_tree->Branch("mu4_d0err_bestYMass",   &mu4_d0err_bestYMass,    "mu4_d0err_bestYMass/F");
		onia_tree->Branch("mu3_dz_bestYMass",   &mu3_dz_bestYMass,    "mu3_dz_bestYMass/F");
		onia_tree->Branch("mu3_dzerr_bestYMass",   &mu3_dzerr_bestYMass,    "mu3_dzerr_bestYMass/F");
		onia_tree->Branch("mu4_dz_bestYMass",   &mu4_dz_bestYMass,    "mu4_dz_bestYMass/F");
		onia_tree->Branch("mu4_dzerr_bestYMass",   &mu4_dzerr_bestYMass,    "mu4_dzerr_bestYMass/F");
		onia_tree->Branch("mu1_Tight_bestYMass",   &mu1_Tight_bestYMass,    "mu1_Tight_bestYMass/I");
		onia_tree->Branch("mu2_Tight_bestYMass",   &mu2_Tight_bestYMass,    "mu2_Tight_bestYMass/I");
		onia_tree->Branch("mu3_Tight_bestYMass",   &mu3_Tight_bestYMass,    "mu3_Tight_bestYMass/I");
		onia_tree->Branch("mu4_Tight_bestYMass",   &mu4_Tight_bestYMass,    "mu4_Tight_bestYMass/I");
		onia_tree->Branch("mu3_pdgID_bestYMass",   &mu3_pdgID_bestYMass,    "mu3_pdgID_bestYMass/I");
		onia_tree->Branch("mu4_pdgID_bestYMass",   &mu4_pdgID_bestYMass,    "mu4_pdgID_bestYMass/I");

	}

	if (isMC_ || OnlyGen_) {
		onia_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/I");
		onia_tree->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/I");
		onia_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
		onia_tree->Branch("gen_mu1_p4",  "TLorentzVector",  &gen_mu1_p4);
		onia_tree->Branch("gen_mu2_p4",  "TLorentzVector",  &gen_mu2_p4);
	}
	genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
	packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
	//genCands_= consumes<reco::GenParticleCollection>((edm::InputTag)"genMuons");

	TFileDirectory hists = fs->mkdir( "hists" );
	myFourMuM_fit = hists.make<TH1F>("myFourMuM_fit","myFourMuM_fit",50000,0,500.0);
	myFourMuM_fit->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}#mu^{+}#mu^{-}} [GeV/c^{2}]");
	myFourMuVtxP_fit = hists.make<TH1F>("myFourMuVtxP_fit","myFourMuVtxP_fit",1000,0,1);
	myFourMuVtxP_fit->GetXaxis()->SetTitle("#mu^{+}#mu^{-}#mu^{+}#mu^{-} vertex fit #chi^{2} probability");
	myDimuonMass_all = hists.make<TH1F>("myDimuonMass_all","myDimuonMass_all",2000,0,20.0);
	myDimuonMass_all->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
	myY1SFitMass_all = hists.make<TH1F>("myY1SFitMass_all","myY1SFitMass_all",2000,0,20.0);
	myY1SFitMass_all->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
}

MuMuGammaRootupler::~MuMuGammaRootupler() {

}

//
// member functions
//

const reco::Candidate* MuMuGammaRootupler::GetAncestor(const reco::Candidate* p) {
	if (p->numberOfMothers()) {
		if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
		else return p->mother(0);
	}
	//std::cout << "GetAncestor: Inconsistet ancestor, particle does not have a mother " << std::endl;
	return p;
}

//Check recursively if any ancestor of particle is the given one
bool MuMuGammaRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
	if (ancestor == particle ) return true;
	for (size_t i=0; i< particle->numberOfMothers(); i++) {
		if (isAncestor(ancestor, particle->mother(i))) return true;
	}
	return false;
}

/* Grab Trigger information. Save it in variable trigger, trigger is an int between 0 and 127, in binary it is:
	(pass 2)(pass 1)(pass 0)
	ex. 7 = pass 0, 1 and 2
	ex. 6 = pass 2, 3
	ex. 1 = pass 0
	*/

bool MuMuGammaRootupler::triggerDecision(edm::Handle<edm::TriggerResults> &hltR, int iTrigger){
  bool triggerPassed = false;
  if(hltR->wasrun(iTrigger) &&
     hltR->accept(iTrigger) &&
     !hltR->error(iTrigger) ){
    triggerPassed = true;
  }
  return triggerPassed;
}


void MuMuGammaRootupler::analyzeTrigger(edm::Handle<edm::TriggerResults> &hltR,
                                       edm::Handle<trigger::TriggerEvent> &hltE,
                                       const std::string& triggerName) {
  using namespace trigger;

  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  if(verbose){
    std::cout<<" n = "<<n<<" triggerIndex = "<<triggerIndex<<" size = "<<hltConfig_.size()<<std::endl;
    std::cout<<" Analyze triggerName : "<<triggerName<<std::endl;
  }
  if (triggerIndex>=n) {
    if(verbose){
      cout << "FourmuonAnalyzer4::analyzeTrigger: path "
           << triggerName << " - not found!" << endl;
    }
    return;
  }
  const unsigned int moduleIndex(hltR->index(triggerIndex));
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));
  lastTriggerModule = moduleIndex;
  if(verbose){
    cout << "FourmuonAnalyzer::analyzeTrigger: path "
         << triggerName << " [" << triggerIndex << "]" << endl;
         
    std::cout<<"  n = "<< n<<" triggerIndex = "<<triggerIndex<<" m = "<<m<<std::endl;
    std::cout<<" moduleLabels = "<<moduleLabels.size()<<" moduleIndex = "<<moduleIndex<<std::endl;
    // Results from TriggerResults product
        cout << " Trigger path status:"
         << " WasRun=" << hltR->wasrun(triggerIndex)
         << " Accept=" << hltR->accept(triggerIndex)
         << " Error =" << hltR->error(triggerIndex)
         << endl;
    cout << " Last active module - label/type: "
         << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
         << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
         << endl;
  }
  assert (moduleIndex<m);
  // Results from TriggerEvent product - Attention: must look only for
  // modules actually run in this path for this event!
    std::vector < GlobalVector > passMomenta;
  for (unsigned int j=0; j<=moduleIndex; ++j) {
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
    // check whether the module is packed up in TriggerEvent product
        const unsigned int filterIndex(hltE->filterIndex(InputTag(moduleLabel,"",hltName_)));
    if(verbose){
      std::cout<<" j = "<<j<<" modLabel/moduleType = "<<moduleLabel<<"/"<<moduleType<<" filterIndex = "<<filterIndex<<" sizeF = "<<hltE->sizeFilters()<<std::endl;
    }
    if (filterIndex<hltE->sizeFilters()) { 
          if(verbose){
        cout << " 'L3' (or 'L1', 'L2') filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
      }
      const Vids& VIDS (hltE->filterIds(filterIndex));
      const Keys& KEYS(hltE->filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI==nK);
      const size_type n(max(nI,nK));
      if(verbose){
        cout << "   " << n  << " accepted 'L3' (or 'L1', 'L2') objects found: " << endl;
      }
      const TriggerObjectCollection& TOC(hltE->getObjects());
      for (size_type i=0; i!=n; ++i) {
        if(0==i){
          passMomenta.clear();
        }
        const TriggerObject& TO(TOC[KEYS[i]]);
        GlobalVector momentumT0(TO.px(),TO.py(),TO.pz());
        if(verbose){
          std::cout<<" i = "<<i<<" moduleLabel/moduleType : "<<moduleLabel<<"/"<<moduleType<<std::endl;
        }
        if(13==TO.id() || -13==TO.id() || 0==TO.id()){//TO.id() --> L1 Mu is always 0 (?)
          if(verbose){
            std::cout<<" current moduleType = "<<moduleType<<std::endl;
          }  
         if("HLTL1TSeed" == moduleType && "hltL1sTripleMu0orTripleMu500" == moduleLabel){ // HLT_Dimuon0_Jpsi_Muon_v5
            passMomenta.push_back(momentumT0);
            if(verbose){
              std::cout<<" L1 object found"<<std::endl;
            }
          }
          else if("HLTMuonL2FromL1TPreFilter"==moduleType){//HLT_Dimuon0_Jpsi_Muon_v5
            passMomenta.push_back(momentumT0);
            if(verbose){
              std::cout<<" L2 object found"<<std::endl;
            }
          }
          else if("HLTMuonL3PreFilter"==moduleType || "HLTMuonIsoFilter"==moduleType){
            passMomenta.push_back(momentumT0);
            if(verbose){
              std::cout<<" L3 (highEff) object found"<<std::endl;
            }
          }
          else if("HLTDiMuonGlbTrkFilter"==moduleType){
            passMomenta.push_back(momentumT0);
            if(verbose){
              std::cout<<" L3 (lowEff) object found"<<std::endl;
            }
          }
          else if("HLT2MuonMuonDZ"==moduleType || "HLTMuonIsoFilter"==moduleType || ("HLTMuonL3PreFilter"==moduleType)){
            passMomenta.push_back(momentumT0);
            if(verbose){
              std::cout<<" HLT object found"<<std::endl;
            }
          }
                   else if("HLTMuonDimuonL3Filter"==moduleType){//HLT_Dimuon0_Jpsi_Muon_v5
            passMomenta.push_back(momentumT0);
            if(verbose){
              std::cout<<" HLT L3 filter object found"<<std::endl;
            }
          }
        }
        if(verbose){
          cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
               << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
               << endl;
        }
      }                    
    }
    if("HLTL1TSeed" == moduleType && "hltL1sTripleMu0orTripleMu500" == moduleLabel){ // HLT_Dimuon0_Jpsi_Muon_v5
      for(unsigned int i=0;i<passMomenta.size();++i){
        allMuL1TriggerVectors.push_back(passMomenta[i]);
      }
      if(verbose){
        std::cout<<" L1 obj FOUND; current size = "<< allMuL1TriggerVectors.size()<<std::endl;
      }
    }
    if("HLTMuonL2FromL1TPreFilter"==moduleType){//HLT_Dimuon0_Jpsi_Muon_v5
      for(unsigned int i=0;i<passMomenta.size();++i){
        allMuL2TriggerVectors.push_back(passMomenta[i]);
      }
          if(verbose){
        std::cout<<" L2 obj FOUND; current size = "<< allMuL2TriggerVectors.size()<<std::endl;
      }
    }
       if("HLTMuonL3PreFilter" == moduleType){
      for(unsigned int i=0;i<passMomenta.size();++i){
        allMuL3TriggerVectors_highEff.push_back(passMomenta[i]);
      }
          if(verbose){
        std::cout<<" L3 (highEff) obj FOUND ; current size = " <<allMuL3TriggerVectors_highEff.size()<<std::endl;
      }
    }
    if("HLTDiMuonGlbTrkFilter"==moduleType){
      for(unsigned int i=0;i<passMomenta.size();++i){
        allMuL3TriggerVectors_lowEff.push_back(passMomenta[i]);
      }
          if(verbose){
        std::cout<<" L3 (lowEff) obj FOUND ; current size = " <<allMuL3TriggerVectors_lowEff.size()<<std::endl;
      }
    }
    if("HLTMuonDimuonL3Filter" == moduleType){
      for(unsigned int i=0;i<passMomenta.size();++i){
        allMuL3TriggerVectors_highEff.push_back(passMomenta[i]);
      }
           if(verbose){
        std::cout<<" HLT L3 filter object FOUND; current size = "<<allMuL3TriggerVectors_highEff.size()<<std::endl;
      }
    }
    if(("HLTMuonL3PreFilter"==moduleType) ){
      for(unsigned int i=0;i<passMomenta.size();++i){
        allMuHLTTriggerVectors.push_back(passMomenta[i]);
      }
            if(verbose){
        std::cout<<" HLT obj FOUND ; current size = " <<allMuHLTTriggerVectors.size()<<std::endl;
      }
    }
    passMomenta.clear();
  }
  return;
}

bool MuMuGammaRootupler::findTrigger(edm::Handle<edm::TriggerResults> &hltR,
                          std::vector < std::string > & triggersToCheck,
                          std::vector < std::string > & triggerNameFound)
{
  triggerNameFound.clear();
  if(verbose){
    std::cout<<" findTrigger()... "<<std::endl;
    if(1==nEventsAnalyzed){
      std::cout<<"   Request: "<<std::endl;
      for(unsigned int iT=0;iT<triggersToCheck.size();++iT){
        std::cout<<"       name ["<<iT<<"] = "<<triggersToCheck[iT]<<std::endl;
      }
    }
  }
  bool triggerFound = false;
  std::vector<std::string>  hlNames=triggerNames.triggerNames();
  for (uint iT=0; iT<hlNames.size(); ++iT) {
        if(verbose && 1==nEventsAnalyzed){
      std::cout<<" iT = "<<iT<<" hlNames[iT] = "<<hlNames[iT]<<
        " : wasrun = "<<hltR->wasrun(iT)<<" accept = "<<
        hltR->accept(iT)<<" !error = "<<
        !hltR->error(iT)<<std::endl;
    }
   for(uint imyT = 0;imyT<triggersToCheck.size();++imyT){
      if(string::npos!=hlNames[iT].find(triggersToCheck[imyT]))
         {
        if(verbose && 1==nEventsAnalyzed){
          std::cout<<" Trigger "<<hlNames[iT]<<" found to be compatible with the requested. "<<std::endl;
        }
        triggerNameFound.push_back(hlNames[iT]);
        if(triggerDecision(hltR, iT)){
          triggerFound = true;
        }
      }
    }
  }
  return triggerFound;
}

UInt_t MuMuGammaRootupler::getTriggerBits(const edm::Event& iEvent) {
	UInt_t itrigger = 0;
        // trigger collection
	edm::Handle<edm::TriggerResults> triggerResults_handle;
	iEvent.getByToken(triggerResultsTok_, triggerResults_handle);

	//   if ( triggerResults_handle.isValid() ) { 
	//     std::string testTriggerName;
	//     const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
	//     for(unsigned int trigIndex = 0; trigIndex < TheTriggerNames.size(); trigIndex++){
	//     testTriggerName = TheTriggerNames.triggerName(trigIndex);
	//     std::cout<<testTriggerName.c_str()<<std::endl;
	//     }
	//   }
	if ( triggerResults_handle.isValid() ) {
                // trigger collection 
		const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
		std::vector <unsigned int> bits_0, bits_1, bits_2, bits_3, bits_4, bits_5, bits_6, bits_7, bits_8, bits_9;
		for ( int version = 1; version<20; version ++ ) {
			std::stringstream ss0,ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9;
			ss0<<"HLT_Dimuon25_Jpsi_v"<<version;
			//ss0<<"HLT_Dimuon16_Jpsi_v"<<version;
			bits_0.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss0.str()).label().c_str()));

			ss1<<"HLT_Dimuon0_Upsilon_Muon_v"<<version; //2016
			//ss1<<"HLT_Dimuon13_PsiPrime_v"<<version;
			bits_1.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss1.str()).label().c_str()));

			ss2<<"HLT_Dimuon13_Upsilon_v"<<version; //2016
			bits_2.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss2.str()).label().c_str()));

			ss3<<"HLT_DiMuon0_Jpsi_Muon_v"<<version;
			//ss3<<"HLT_Dimuon10_Jpsi_Barrel_v"<<version;
			bits_3.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss3.str()).label().c_str()));

			ss4<<"HLT_TripleMu5_v"<<version; 
			//ss4<<"HLT_Dimuon8_PsiPrime_Barrel_v"<<version;
			bits_4.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss4.str()).label().c_str()));

			ss5<<"HLT_Dimuon12_Upsilon_y1p4_v"<<version;
			//ss5<<"HLT_Dimuon8_Upsilon_Barrel_v"<<version; //2016
			bits_5.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss5.str()).label().c_str()));

			ss6<<"HLT_Dimuon20_Jpsi_v"<<version;
			bits_6.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss6.str()).label().c_str()));

			ss7<<"HLT_Trimuon2_Upsilon5_Muon_v"<<version; //2017B
			//ss7<<"HLT_Dimuon14_Phi_Barrel_Seagulls_v"<<version;
			//ss7<<"HLT_Dimuon0_Phi_Barrel_v"<<version;//2016
			bits_7.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss7.str()).label().c_str()));

			ss8<<"HLT_Trimuon5_3p5_2_Upsilon_Muon_v"<<version; //2017C/D/E/F and 2018
			//ss8<<"HLT_DoubleMu4_3_Bs_v"<<version;
			bits_8.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss8.str()).label().c_str()));

			ss9<<"HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v"<<version;//2018
			//ss9<<"HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v"<<version;
			bits_9.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss9.str()).label().c_str()));
		}
		for (unsigned int i=0; i<bits_0.size(); i++) {
			unsigned int bit = bits_0[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 1;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_1.size(); i++) {
			unsigned int bit = bits_1[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 2;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_2.size(); i++) {
			unsigned int bit = bits_2[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 4;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_3.size(); i++) {
			unsigned int bit = bits_3[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 8;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_4.size(); i++) {
			unsigned int bit = bits_4[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 16;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_5.size(); i++) {
			unsigned int bit = bits_5[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 32;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_6.size(); i++) {
			unsigned int bit = bits_6[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 64;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_7.size(); i++) {
			unsigned int bit = bits_7[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 128;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_8.size(); i++) {
			unsigned int bit = bits_8[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) { 
					itrigger += 256;
					break;
				}   
			}   
		}   
		for (unsigned int i=0; i<bits_9.size(); i++) {
			unsigned int bit = bits_9[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) { 
					itrigger += 512;
					break;
				}   
			}   
		} 
            } 
	return itrigger;
}

bool MuMuGammaRootupler::TriggerMatch(bool TriggerPassed, pat::CompositeCandidate dimuonCand) {     
       allTrigMuons.clear(); 
       if (verbose) cout<< "Trigger mathcing for dimuon candidate"<<endl;
       double reco1_eta = dimuonCand.daughter("muon1")->eta();
       double reco1_phi = dimuonCand.daughter("muon1")->phi();
       double reco1_pt = dimuonCand.daughter("muon1")->pt();
       double reco1_mass = dimuonCand.daughter("muon1")->mass();
       double reco2_eta = dimuonCand.daughter("muon2")->eta();
       double reco2_phi = dimuonCand.daughter("muon2")->phi();
       double reco2_pt = dimuonCand.daughter("muon2")->pt();
       double reco2_mass = dimuonCand.daughter("muon2")->mass();
       if (verbose) cout<<"This Dimuon candidate"<<" mu1pt:" <<reco1_pt<<" mu1eta:"<<reco1_eta<<" mu1phi:"<<reco1_phi<<endl;
       if (verbose) cout<<"This Dimuon candidate"<<" mu2pt:" <<reco2_pt<<" mu1eta:"<<reco2_eta<<" mu1phi:"<<reco2_phi<<endl;
       float dR1 = -9999.;
       float dR2 = -9999.;
       float dR1_minimum = 99;
       float dR2_minimum = 99;
       float dPt1 = 999;
       float dPt2 = 999;
       TLorentzVector TempMomenta1;
       TLorentzVector TempMomenta2;
       if (verbose) cout<<"allMuHLTTriggerVectors.size():"<<allMuHLTTriggerVectors.size()<<endl;
       for(uint iTrig =0;iTrig<allMuHLTTriggerVectors.size();++iTrig){
           double hlt_pt = allMuHLTTriggerVectors[iTrig].perp();          
           double hlt_eta = allMuHLTTriggerVectors[iTrig].eta();
           double hlt_phi = allMuHLTTriggerVectors[iTrig].phi();
           double dR1 =  deltaR(reco1_eta,reco1_phi,hlt_eta,hlt_phi);
           double dR2 =  deltaR(reco2_eta,reco2_phi,hlt_eta,hlt_phi);
           if (verbose) cout<<"iTrig:"<<iTrig<<" iTrigpT:"<<hlt_pt<<" hlt_phi:"<<hlt_phi<<endl;
           if (verbose) cout<<"dR1:" <<dR1<<" dPt1:"<<dPt1<<endl;
           if (verbose) cout<<"dR2:" <<dR2<<" dPt2:"<<dPt2<<endl; 
           if (dR1<dR2)
           {
            if (dR1<dR1_minimum)       
            {
            dR1_minimum = dR1;   
            dPt1 = std::abs(reco1_pt - hlt_pt); 
            if(dR1 < trg_Match_dR_cut && dPt1 < trg_Match_dP_cut)      
              { 
               if (verbose) cout<<"Matching L3 sucessfull muPt1 = " <<reco1_pt<<" trigPt = "<<hlt_pt<<" dR = "<<dR1<<endl;
               TempMomenta1.SetPtEtaPhiM(reco1_pt,reco1_eta,reco1_phi,reco1_mass);
                } //Found matching to muon 1
               else cout<<"Matching failed -> iTrig = "<< hlt_pt<<" eta = "<<hlt_eta<<" phi = "<< hlt_phi<<" dR ="<<dR1<<endl;
               } // checking best matching object with muon 1 
            } // matching to muon1 
           if (dR2<dR1)
            {
             if (dR2<dR2_minimum)
              {
             dR2_minimum = dR2;
             dPt2 = std::abs(reco2_pt - hlt_pt);   
             if(dR2 < trg_Match_dR_cut && dPt2 < 0.3)
               {
               if (verbose) cout<<"Matching L3 sucessfull muPt2 = " <<reco2_pt<<" trigPt = "<<hlt_pt<<" dR = "<<dR2<<endl;
               TempMomenta2.SetPtEtaPhiM(reco2_pt,reco2_eta,reco2_phi,reco2_mass);
                }
             else cout<<"Matching failed -> iTrig = "<< hlt_pt<<" eta = "<<hlt_eta<<" phi = "<< hlt_phi<<" dR ="<<dR2<<endl;
                 } //checking best matching object  with muon 2
             } // matching to muon2
         } // loop over all HLT objects
        dR1 = dR1_minimum;
        dR2 = dR2_minimum;
        if( dR1 < trg_Match_dR_cut && dPt1 < trg_Match_dP_cut ){
          allTrigMuons.push_back(TempMomenta1);
           } // filling vector of matched muon 1
        if(dR2 < trg_Match_dR_cut && dPt2 < trg_Match_dP_cut ){
          allTrigMuons.push_back(TempMomenta2);
           } // filling vector of matched muon 2      
        if (verbose)
        {
        cout<<" AllTrigMu = "<<allTrigMuons.size()<<endl;
        cout<<" good fourMu  run/lumi/event : "<<run<<"/"<<lumi<<"/"<<event<<std::endl;
        if(TriggerPassed){
           cout<<"   Trigger passed "<<endl;
        if(allTrigMuons.size()>=2){
          cout<<"Matching good "<<endl;
             }
             }
          }
      if (allTrigMuons.size()>=2) return true;
      else return false;
     }
bool MuMuGammaRootupler::TriggerMatch_restMuons(TLorentzVector mu3p4, TLorentzVector mu4p4) {
       allRestTrigMuons.clear();
       if (verbose) cout<< "Trigger matching for Rest of muons candidates"<<endl;
       double reco1_eta = mu3p4.Eta();
       double reco1_phi = mu3p4.Phi();
       double reco1_pt = mu3p4.Pt();
       double reco1_mass = mu3p4.M();
       double reco2_eta = mu4p4.Eta();
       double reco2_phi = mu4p4.Phi();
       double reco2_pt = mu4p4.Pt();
       double reco2_mass = mu4p4.M();
       if (verbose) cout<<"First muon candidate"<<" mu1pt:" <<reco1_pt<<" mu1eta:"<<reco1_eta<<" mu1phi:"<<reco1_phi<<endl;
       if (verbose) cout<<"Second muon candidate"<<" mu2pt:" <<reco2_pt<<" mu2eta:"<<reco2_eta<<" mu2phi:"<<reco2_phi<<endl;
       float dR1 = -9999.;
       float dR2 = -9999.;
       float dR1_minimum = 99;
       float dR2_minimum = 99;
       float dPt1 = 999;
       float dPt2 = 999;
       TLorentzVector TempMomenta1;
       TLorentzVector TempMomenta2;
       if (verbose) cout<<"allMuHLTTriggerVectors.size():"<<allMuHLTTriggerVectors.size()<<endl;
       for(uint iTrig =0;iTrig<allMuHLTTriggerVectors.size();++iTrig){
           double hlt_pt = allMuHLTTriggerVectors[iTrig].perp();
           double hlt_eta = allMuHLTTriggerVectors[iTrig].eta();
           double hlt_phi = allMuHLTTriggerVectors[iTrig].phi();
           double dR1 =  deltaR(reco1_eta,reco1_phi,hlt_eta,hlt_phi);
           double dR2 =  deltaR(reco2_eta,reco2_phi,hlt_eta,hlt_phi);
           if (verbose) cout<<"iTrig:"<<iTrig<<" iTrigpT:"<<hlt_pt<<" hlt_phi:"<<hlt_phi<<endl;
           if (verbose) cout<<"dR1:" <<dR1<<" dPt1:"<<dPt1<<endl;
           if (verbose) cout<<"dR2:" <<dR2<<" dPt2:"<<dPt2<<endl;
           if (dR1<dR2)
           {
            if (dR1<dR1_minimum)
            {
            dR1_minimum = dR1;
            dPt1 = std::abs(reco1_pt - hlt_pt);
            if(dR1 < trg_Match_dR_cut && dPt1 < trg_Match_dP_cut)
              {
               if (verbose) cout<<"Matching L3 sucessfull muPt1 = " <<reco1_pt<<" trigPt = "<<hlt_pt<<" dR = "<<dR1<<endl;
               TempMomenta1.SetPtEtaPhiM(reco1_pt,reco1_eta,reco1_phi,reco1_mass);
                } //Found matching to muon 1
               else cout<<"Matching failed -> iTrig = "<< hlt_pt<<" eta = "<<hlt_eta<<" phi = "<< hlt_phi<<" dR ="<<dR1<<endl;
               } // checking best matching object with muon 1 
            } // matching to muon1 
           if (dR2<dR1)
            {
             if (dR2<dR2_minimum)
              {
             dR2_minimum = dR2;
             dPt2 = std::abs(reco2_pt - hlt_pt);
             if(dR2 < trg_Match_dR_cut && dPt2 < 0.3)
               {
               if (verbose) cout<<"Matching L3 sucessfull muPt2 = " <<reco2_pt<<" trigPt = "<<hlt_pt<<" dR = "<<dR2<<endl;
               TempMomenta2.SetPtEtaPhiM(reco2_pt,reco2_eta,reco2_phi,reco2_mass);
                }
             else cout<<"Matching failed -> iTrig = "<< hlt_pt<<" eta = "<<hlt_eta<<" phi = "<< hlt_phi<<" dR ="<<dR2<<endl;
                 } //checking best matching object  with muon 2
             } // matching to muon2
         } // loop over all HLT objects
        dR1 = dR1_minimum;
        dR2 = dR2_minimum;
        if( dR1 < trg_Match_dR_cut && dPt1 < trg_Match_dP_cut ){
          allRestTrigMuons.push_back(TempMomenta1);
           } // filling vector of matched muon 1
        if(dR2 < trg_Match_dR_cut && dPt2 < trg_Match_dP_cut ){
          allRestTrigMuons.push_back(TempMomenta2);
           } // filling vector of matched muon 2      
        if (verbose)
        {
        cout<<" All Rest Trigger matched muons size = "<<allRestTrigMuons.size()<<endl;
        cout<<" run/lumi/event : "<<run<<"/"<<lumi<<"/"<<event<<std::endl;
        if(allRestTrigMuons.size()>=1){
          cout<<"Atleast one of Rest muon Match with HLT Object "<<endl;
             }
          }
      if (allRestTrigMuons.size()>=1) return true;
      else return false;
     }
// ------------ method called for each event  ------------
void MuMuGammaRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

        using namespace trigger;
        using namespace pat;
	edm::Handle<pat::CompositeCandidateCollection> dimuons;
	iEvent.getByToken(dimuon_Label,dimuons);

	edm::Handle<pat::CompositeCandidateCollection> conversions;
	iEvent.getByToken(conversion_Label,conversions);

	edm::Handle< edm::View<pat::Muon> > muons;
	iEvent.getByToken(muon_Label, muons);

	edm::Handle<reco::VertexCollection> primaryVertices_handle;
	iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);
	reco::Vertex thePrimaryV;

	edm::Handle<reco::BeamSpot> theBeamSpot;
	iEvent.getByToken(bs_Label,theBeamSpot);
	reco::BeamSpot bs = *theBeamSpot;

	if ( primaryVertices_handle->begin() != primaryVertices_handle->end() ) {  
		thePrimaryV = reco::Vertex(*(primaryVertices_handle->begin()));
	}
	else {
		thePrimaryV = reco::Vertex(bs.position(), bs.covariance3D());
	}
	pv_x = thePrimaryV.x();
	pv_y = thePrimaryV.y();
	pv_z = thePrimaryV.z();

	/*edm::Handle<reco::ConversionCollection> conversionHandle;
	  iEvent.getByLabel(conversion_Label,conversionHandle);

	  edm::Handle<reco::PFCandidateCollection> pfcandidates;
	  iEvent.getByLabel("particleFlow",pfcandidates);
	  const reco::PFCandidateCollection pfphotons = selectPFPhotons(*pfcandidates);
	  */
         ++nEventsAnalyzed;

      if(int(iEvent.id().run())!=runNumber){
      runNumber = iEvent.id().run();
      if(verbose){
        std::cout<<" New run : "<<iEvent.id().run()<<std::endl;
       }
    const edm::Run * iRun_c = &iEvent.getRun();
    edm::Run * iRun = const_cast <edm::Run*> (iRun_c);
    beginRun(*iRun, iSetup);
      }

    edm::Handle<edm::TriggerResults> hltR;
    edm::Handle<trigger::TriggerEvent> hltE;
    iEvent.getByToken(triggerEventTok_,hltE);
    iEvent.getByToken(triggerResultsTok_, hltR);
    std::vector<std::string>  hlNames;
    triggerNames = iEvent.triggerNames(*hltR);
    hlNames=triggerNames.triggerNames();
    std::vector < std::string > triggersFoundToApply;
      bool theTriggerPassed = findTrigger(hltR, triggerList, triggersFoundToApply);
      if (theTriggerPassed){
         cout<<" theTriggerPassed = "<<theTriggerPassed<<" checkTrigger = "<<checkTrigger<<" trigFound = "<<triggersFoundToApply.size()<<endl;
        }
//
 
	if (!OnlyGen_) {
		numPrimaryVertices = -1;
		if (primaryVertices_handle.isValid()) 	numPrimaryVertices = (int) primaryVertices_handle->size();

		trigger = getTriggerBits(iEvent);
		run     = iEvent.id().run();
		lumi    = iEvent.id().luminosityBlock();
		event   = iEvent.id().event();
}	

	//if (run < 316569)		//a temporary run number selection, 316569 is the first run of 2018A prompt reco v3
	if(true)
	{

	dimuon_pdgId = 0;
	mother_pdgId = 0;
	irank = 0;

	v_mu1Charge.clear();
	v_mu2Charge.clear();
	v_mu1_d0.clear();
	v_mu1_d0err.clear();
	v_mu2_d0.clear();
	v_mu2_d0err.clear();
	v_mu1_dz.clear();
	v_mu1_dzerr.clear();
	v_mu2_dz.clear();
	v_mu2_dzerr.clear();
	v_mu1_vz.clear();
	v_mu2_vz.clear();
	v_mumufit_Mass.clear();
	v_mumufit_MassErr.clear();
	v_mumufit_VtxCL.clear();
	v_mumufit_VtxCL2.clear();
	v_mumufit_DecayVtxX.clear();
	v_mumufit_DecayVtxY.clear();
	v_mumufit_DecayVtxZ.clear();
	v_mumufit_DecayVtxXE.clear();
	v_mumufit_DecayVtxYE.clear();
	v_mumufit_DecayVtxZE.clear();
        mu1_filtersMatched.clear();
        mu2_filtersMatched.clear();
	dimuon_p4.SetPtEtaPhiM(0,0,0,0);
	mu1_p4.SetPtEtaPhiM(0,0,0,0);
	mu2_p4.SetPtEtaPhiM(0,0,0,0);
	mu1Charge = -10; 
	mu2Charge = -10; 
	mu1_d0 = -10; 
	mu1_d0err = -10; 
	mu2_d0 = -10; 
	mu2_d0err = -10; 
	mu1_dz = -1000;
	mu1_dzerr = -1000;
	mu2_dz = -1000;
	mu2_dzerr = -1000;
	mu1_vz = -1000;
	mu2_vz = -1000;
	mumufit_Mass = -10; 
	mumufit_MassErr = -10; 
	mumufit_VtxCL = -10; 
	mumufit_VtxCL2 = -10; 
	mumufit_DecayVtxX = -10; 
	mumufit_DecayVtxY = -10; 
	mumufit_DecayVtxZ = -10; 
	mumufit_DecayVtxXE = -10; 
	mumufit_DecayVtxYE = -10; 
	mumufit_DecayVtxZE = -10; 
	mumufit_p4.SetPtEtaPhiM(0,0,0,0);


	fourMuFit_Mass_allComb_mix.clear();
	fourMuFit_Mass_mix = -1;
	fourMuFit_MassErr_mix = -1;
	fourMuFit_VtxX_mix = -10;
	fourMuFit_VtxY_mix = -10;
	fourMuFit_VtxZ_mix = -100;
	fourMuFit_VtxProb_mix = -1;
	fourMuFit_Chi2_mix = -10;
	fourMuFit_ndof_mix = -1;
	fourMuFit_3plus1_mix = -1;
	fourMuFit_p4_mix.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu1p4_mix.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu2p4_mix.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu3p4_mix.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu4p4_mix.SetPtEtaPhiM(0,0,0,0);
	mu3_p4_mix.SetPtEtaPhiM(0,0,0,0);
	mu4_p4_mix.SetPtEtaPhiM(0,0,0,0);
	mu3_d0_mix = -10;
	mu3_d0err_mix = -10;
	mu4_d0_mix = -10;
	mu4_d0err_mix = -10;
	mu3_dz_mix = -100;
	mu3_dzerr_mix = -100;
	mu4_dz_mix = -100;
	mu4_dzerr_mix = -100;
	mu3Charge_mix = -10;
	mu4Charge_mix = -10;
	fourMuFit_Mass_mix3evts = -1;
	fourMuFit_VtxProb_mix3evts = -1;
	fourMuFit_p4_mix3evts.SetPtEtaPhiM(0,0,0,0);

	genbkg_mu1_Pt.clear();
	genbkg_mu1_Eta.clear();
	genbkg_mu1_Phi.clear();
	genbkg_mu1_Mass.clear();
	genbkg_mu2_Pt.clear();
	genbkg_mu2_Eta.clear();
	genbkg_mu2_Phi.clear();
	genbkg_mu2_Mass.clear();
	genbkg_mu3_Pt.clear();
	genbkg_mu3_Eta.clear();
	genbkg_mu3_Phi.clear();
	genbkg_mu3_Mass.clear();
	genbkg_mu4_Pt.clear();
	genbkg_mu4_Eta.clear();
	genbkg_mu4_Phi.clear();
	genbkg_mu4_Mass.clear();


	fourMuFit_Mass_allComb.clear();
	fourMuFit_Mass.clear();
	fourMuFit_MassErr.clear();
	fourMuFit_Pt.clear();
	fourMuFit_Eta.clear();
	fourMuFit_Phi.clear();
	fourMuFit_VtxX.clear();
	fourMuFit_VtxY.clear();
	fourMuFit_VtxZ.clear();
	fourMuFit_VtxProb.clear();
	fourMuFit_Chi2.clear();
	fourMuFit_ndof.clear();
	fourMuFit_mu1Pt.clear();
	fourMuFit_mu1Eta.clear();
	fourMuFit_mu1Phi.clear();
	fourMuFit_mu1E.clear();
	fourMuFit_mu2Pt.clear();
	fourMuFit_mu2Eta.clear();
	fourMuFit_mu2Phi.clear();
	fourMuFit_mu2E.clear();
	fourMuFit_mu3Pt.clear();
	fourMuFit_mu3Eta.clear();
	fourMuFit_mu3Phi.clear();
	fourMuFit_mu3E.clear();
	fourMuFit_mu4Pt.clear();
	fourMuFit_mu4Eta.clear();
	fourMuFit_mu4Phi.clear();
	fourMuFit_mu4E.clear();
	mu3_Pt.clear();
	mu3_Eta.clear();
	mu3_Phi.clear();
	mu3_E.clear();
	mu4_Pt.clear();
	mu4_Eta.clear();
	mu4_Phi.clear();
	mu4_E.clear();
	mu3_d0.clear();
	mu3_d0err.clear();
	mu4_d0.clear();
	mu4_d0err.clear();
	mu3_dz.clear();
	mu3_dzerr.clear();
	mu4_dz.clear();
	mu4_dzerr.clear();
	mu3Charge.clear();
	mu4Charge.clear();
	mu1_Tight.clear();
	mu2_Tight.clear();
	mu3_Tight.clear();
	mu4_Tight.clear();
	mu1_Medium.clear();
	mu2_Medium.clear();
	mu3_Medium.clear();
	mu4_Medium.clear();
	mu1_Loose.clear();
	mu2_Loose.clear();
	mu3_Loose.clear();
	mu4_Loose.clear();
	mu1_pdgID.clear();
	mu2_pdgID.clear();
	mu3_pdgID.clear();
	mu4_pdgID.clear();

	/*
		fourMuFit_Mass_allComb.clear();
		fourMuFit_Mass = -1;
		fourMuFit_MassErr = -1;
		fourMuFit_VtxX = -10;
		fourMuFit_VtxY = -10;
		fourMuFit_VtxZ = -100;
		fourMuFit_VtxProb = -1;
		fourMuFit_Chi2 = -10;
		fourMuFit_ndof = -1;
		fourMuFit_p4.SetPtEtaPhiM(0,0,0,0);
		fourMuFit_mu1p4.SetPtEtaPhiM(0,0,0,0);
		fourMuFit_mu2p4.SetPtEtaPhiM(0,0,0,0);
		fourMuFit_mu3p4.SetPtEtaPhiM(0,0,0,0);
		fourMuFit_mu4p4.SetPtEtaPhiM(0,0,0,0);
		mu3_p4.SetPtEtaPhiM(0,0,0,0);
		mu4_p4.SetPtEtaPhiM(0,0,0,0);
		mu3_d0 = -10;
		mu3_d0err = -10;
		mu4_d0 = -10;
		mu4_d0err = -10;
		mu3_dz = -100;
		mu3_dzerr = -100;
		mu4_dz = -100;
		mu4_dzerr = -100;
		mu3Charge = -10; 
		mu4Charge = -10; 
		mu1_Tight = -1;
		mu2_Tight = -1;
		mu3_Tight = -1;
		mu4_Tight = -1;
		mu3_pdgID = -1;
		mu4_pdgID = -1;
		theRestMuons.clear();
		upsilonMuons.clear();
		*/
	gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_mu1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_mu2_p4.SetPtEtaPhiM(0.,0.,0.,0.);


	dimuon_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu1_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu2_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu1Charge_bestYMass = -10; 
	mu2Charge_bestYMass = -10; 
	mu1_d0_bestYMass = -10; 
	mu1_d0err_bestYMass = -10; 
	mu2_d0_bestYMass = -10; 
	mu2_d0err_bestYMass = -10; 
	mu1_dz_bestYMass = -1000;
	mu1_dzerr_bestYMass = -1000;
	mu2_dz_bestYMass = -1000;
	mu2_dzerr_bestYMass = -1000;
	mumufit_Mass_bestYMass = -10; 
	mumufit_MassErr_bestYMass = -10; 
	mumufit_VtxCL_bestYMass = -10; 
	mumufit_VtxCL2_bestYMass = -10; 
	mumufit_DecayVtxX_bestYMass = -10; 
	mumufit_DecayVtxY_bestYMass = -10; 
	mumufit_DecayVtxZ_bestYMass = -10; 
	mumufit_DecayVtxXE_bestYMass = -10; 
	mumufit_DecayVtxYE_bestYMass = -10; 
	mumufit_DecayVtxZE_bestYMass = -10; 
	mumufit_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	bestVertex_and_bestYMass = -1;

	fourMuFit_Mass_allComb_mix_bestYMass.clear();
	fourMuFit_Mass_mix_bestYMass = -1;
	fourMuFit_MassErr_mix_bestYMass = -1;
	fourMuFit_VtxX_mix_bestYMass = -10;
	fourMuFit_VtxY_mix_bestYMass = -10;
	fourMuFit_VtxZ_mix_bestYMass = -100;
	fourMuFit_VtxProb_mix_bestYMass = -1;
	fourMuFit_Chi2_mix_bestYMass = -10;
	fourMuFit_ndof_mix_bestYMass = -1;
	fourMuFit_3plus1_mix_bestYMass = -1;
	fourMuFit_p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu1p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu2p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu3p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu4p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu3_p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu4_p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu3_d0_mix_bestYMass = -10;
	mu3_d0err_mix_bestYMass = -10;
	mu4_d0_mix_bestYMass = -10;
	mu4_d0err_mix_bestYMass = -10;
	mu3_dz_mix_bestYMass = -100;
	mu3_dzerr_mix_bestYMass = -100;
	mu4_dz_mix_bestYMass = -100;
	mu4_dzerr_mix_bestYMass = -100;
	mu3Charge_mix_bestYMass = -10;
	mu4Charge_mix_bestYMass = -10;

	fourMuFit_Mass_allComb_bestYMass.clear();
	fourMuFit_Mass_bestYMass = -1;
	fourMuFit_MassErr_bestYMass = -1;
	fourMuFit_VtxX_bestYMass = -10;
	fourMuFit_VtxY_bestYMass = -10;
	fourMuFit_VtxZ_bestYMass = -100;
	fourMuFit_VtxProb_bestYMass = -1;
	fourMuFit_Chi2_bestYMass = -10;
	fourMuFit_ndof_bestYMass = -1;
	fourMuFit_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu1p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu2p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu3p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu4p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu3_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu4_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu3_d0_bestYMass = -10;
	mu3_d0err_bestYMass = -10;
	mu4_d0_bestYMass = -10;
	mu4_d0err_bestYMass = -10;
	mu3_dz_bestYMass = -100;
	mu3_dzerr_bestYMass = -100;
	mu4_dz_bestYMass = -100;
	mu4_dzerr_bestYMass = -100;
	mu3Charge_bestYMass = -10;
	mu4Charge_bestYMass = -10;
	mu1_Tight_bestYMass = -1;
	mu2_Tight_bestYMass = -1;
	mu3_Tight_bestYMass = -1;
	mu4_Tight_bestYMass = -1;
	mu3_pdgID_bestYMass = -1;
	mu4_pdgID_bestYMass = -1;

// Trigger

  allL1TrigMuons.clear();
  allL2TrigMuons.clear();
  allL3TrigMuons.clear();
  allMuL1TriggerVectors.clear();
  allMuL2TriggerVectors.clear();
  allMuL3TriggerVectors_lowEff.clear();
  allMuL3TriggerVectors_highEff.clear();
  allMuHLTTriggerVectors.clear();
  if (verbose) cout<<"triggersFoundToApply.size()"<<triggersFoundToApply.size()<<endl;
  for(unsigned int iTrig=0;iTrig<triggersFoundToApply.size();++iTrig){
    lastTriggerModule = -1;
    if (verbose) cout<<"triggersFoundToApply.at(iTrig)"<<triggersFoundToApply.at(iTrig)<<endl;
    analyzeTrigger(hltR, hltE, triggersFoundToApply.at(iTrig));
                 }
         if (verbose)cout<<"Trigger analyzed Finished"<<endl;

	// Pruned particles are the one containing "important" stuff
	edm::Handle<reco::GenParticleCollection> pruned;
	iEvent.getByToken(genCands_, pruned);

	// Packed particles are all the status 1, so usable to remake jets
	// The navigation from status 1 to pruned is possible (the other direction should be made by hand)
	edm::Handle<pat::PackedGenParticleCollection> packed;
	iEvent.getByToken(packCands_,  packed);

	int nGoodGenCand = 0;

	TLorentzVector gen_bkg_mu1_p4,gen_bkg_mu2_p4,gen_bkg_mu12_p4;
	if ((isMC_ || OnlyGen_) && pruned.isValid()) {
		std::cout<<pruned->size()<<std::endl;
		for (size_t j=0; j<pruned->size(); j++) {
			const reco::Candidate * d1 = &(*pruned)[j];
			if (d1->status()!=1) continue;
			//std::cout<<d1->pdgId()<<std::endl;
			for (size_t k=j+1; k<pruned->size(); k++) {
				const reco::Candidate * d2 = &(*pruned)[k];
				if (d2->status()!=1) continue;
				if (abs(d1->pdgId()) == 13 && abs(d2->pdgId()) == 13 && d1->pdgId()+d2->pdgId() == 0) {
					gen_bkg_mu1_p4.SetPtEtaPhiM(d1->pt(),d1->eta(),d1->phi(),d1->mass());
					gen_bkg_mu2_p4.SetPtEtaPhiM(d2->pt(),d2->eta(),d2->phi(),d2->mass());
					gen_bkg_mu12_p4=gen_bkg_mu1_p4+gen_bkg_mu2_p4;
				}
				if (gen_bkg_mu12_p4.M()<9.2 || gen_bkg_mu12_p4.M()>9.7) continue;
				for (size_t l=0; l<pruned->size(); l++) {
					if  (l==j || l==k) continue;
					const reco::Candidate * d3 = &(*pruned)[l]; 
					if (d3->status()!=1) continue;
					for (size_t m=l+1; m<pruned->size(); m++) {
						if (m==j || m==k) continue;
						const reco::Candidate * d4 = &(*pruned)[k];
						if (d4->status()!=1) continue;
						if (abs(d3->pdgId()) == 13 && abs(d4->pdgId()) == 13 && d3->pdgId()+d4->pdgId() == 0) {
							nGoodGenCand++;
							genbkg_mu1_Pt.push_back(d1->pt());
							genbkg_mu1_Eta.push_back(d1->eta());
							genbkg_mu1_Phi.push_back(d1->phi());
							genbkg_mu1_Mass.push_back(d1->mass());
							genbkg_mu2_Pt.push_back(d2->pt());
							genbkg_mu2_Eta.push_back(d2->eta());
							genbkg_mu2_Phi.push_back(d2->phi());
							genbkg_mu2_Mass.push_back(d2->mass());
							genbkg_mu3_Pt.push_back(d3->pt());
							genbkg_mu3_Eta.push_back(d3->eta());
							genbkg_mu3_Phi.push_back(d3->phi());
							genbkg_mu3_Mass.push_back(d3->mass());
							genbkg_mu4_Pt.push_back(d4->pt());
							genbkg_mu4_Eta.push_back(d4->eta());
							genbkg_mu4_Phi.push_back(d4->phi());
							genbkg_mu4_Mass.push_back(d4->mass());
						}
					}
				}
			}
		}
	}
	if (isMC_ && nGoodGenCand>0) gen_tree->Fill();


	/*

		if ((isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid()) {
		dimuon_pdgId  = 0;
		gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
		int foundit   = 0;

		for (size_t i=0; i<pruned->size(); i++) {
		int p_id = abs((*pruned)[i].pdgId());
		const reco::Candidate *aonia = &(*pruned)[i];
		if (( p_id == pdgid_ ) && (aonia->status() == 2)) {
		dimuon_pdgId = p_id;
		foundit++;
		for (size_t j=0; j<packed->size(); j++) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
		const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
		const reco::Candidate * d = &(*packed)[j];
		if ( motherInPrunedCollection != nullptr && (d->pdgId() == 13 ) && isAncestor(aonia , motherInPrunedCollection) ){
		gen_mu2_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
		foundit++;
		} 
		if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aonia , motherInPrunedCollection) ) {
		gen_mu1_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
		foundit++;
		}
		if ( foundit == 3 ) break;               
		}
		if ( foundit == 3 ) {
		gen_dimuon_p4 = gen_mu2_p4 + gen_mu1_p4;   // this should take into account FSR
		mother_pdgId  = GetAncestor(aonia)->pdgId();
		break;
		} else {
		foundit = 0;
		dimuon_pdgId = 0;
		mother_pdgId = 0;
		gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
		}            
		}  // if ( p_id
		} // for (size

	// sanity check
	//if ( ! dimuon_pdgId ) std::cout << "MuMuGammaRootupler: does not found the given decay " << run << "," << event << std::endl;
	} */ // end if isMC

	float OniaMassMax_ = OniaMassCuts_[1];
	float OniaMassMin_ = OniaMassCuts_[0];

	// Kinematic fit

	edm::ESHandle<TransientTrackBuilder> theB; 
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 
	edm::ESHandle<MagneticField> bFieldHandle;
	iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

	int nGoodUpsilonCand = 0;
	float bestYMass = 1000;
	pat::CompositeCandidate DimuonCand_bestYMass;
	if ( ! OnlyGen_ 
			&& dimuons.isValid() && dimuons->size() > 0) {
		for(pat::CompositeCandidateCollection::const_iterator dimuonCand=dimuons->begin();dimuonCand!= dimuons->end(); ++dimuonCand)
		{
			if (dimuonCand->mass() < OniaMassMin_ || dimuonCand->mass() > OniaMassMax_) continue;
			if (dimuonCand->daughter("muon1")->charge() == dimuonCand->daughter("muon2")->charge() ) continue;
			if (dimuonCand->daughter("muon1")->pt()<2.0 || dimuonCand->daughter("muon2")->pt()<2.0 ) continue;
			if (fabs(dimuonCand->daughter("muon1")->eta())>2.4|| fabs(dimuonCand->daughter("muon2")->eta())>2.4) continue;

			//dimuon refit. 
			//Here we use the KinematicParticleVertexFitter with muon mass. But in the Onia2MuMu skim, it was just KalmanVertexFitter. 
			//Fitted Vertex from both methods are the same (dimuonCand->userData<reco::Vertex>("commonVertex")->z() == mumu_vFit_vertex_noMC->position().z()), 
			//but fitted p4 and mass are different. 
			reco::TrackRef JpsiTk[2]={  //this is from Chib code
				( dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon1") ) )->innerTrack(),
				( dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon2") ) )->innerTrack()
			};

			//std::vector<reco::TransientTrack> MuMuTT;
			//MuMuTT.push_back((*theB).build(&JpsiTk[0]));
			//MuMuTT.push_back((*theB).build(&JpsiTk[1]));	
			reco::TransientTrack muon1TT(JpsiTk[0], &(*bFieldHandle));
			reco::TransientTrack muon2TT(JpsiTk[1], &(*bFieldHandle));
			//std::cout<<muon1TT.isValid()<<" "<<muon1TT.track().pt()<<std::endl;
			reco::TrackTransientTrack muon1TTT(JpsiTk[0], &(*bFieldHandle));
			reco::TrackTransientTrack muon2TTT(JpsiTk[1], &(*bFieldHandle));

			KinematicParticleFactoryFromTransientTrack pmumuFactory;
			std::vector<RefCountedKinematicParticle> mumuParticles;
			mumuParticles.push_back(pmumuFactory.particle(muon1TT,muonMass,float(0),float(0),muonSigma));
			mumuParticles.push_back(pmumuFactory.particle(muon2TT,muonMass,float(0),float(0),muonSigma));

			KinematicParticleVertexFitter mumufitter;
			RefCountedKinematicTree mumuVertexFitTree;
			//try {mumuVertexFitTree = mumufitter.fit(mumuParticles);}
			//catch (...) {
			//std::cout<<"mumu fit: PerigeeKinematicState::kinematic state passed is not valid!"<<std::endl;
			//continue;
			//}   
			mumuVertexFitTree = mumufitter.fit(mumuParticles);
			RefCountedKinematicParticle mumu_vFit_noMC;
			RefCountedKinematicVertex mumu_vFit_vertex_noMC;

			if (!(mumuVertexFitTree->isValid())) continue;
			mumuVertexFitTree->movePointerToTheTop();
			mumu_vFit_noMC = mumuVertexFitTree->currentParticle();
			mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex();

			//if (mumu_vFit_noMC->currentState().mass() < 8 || mumu_vFit_noMC->currentState().mass() > 12) continue;
			if (fabs(mumu_vFit_noMC->currentState().mass()-upsilon_mass_) > (3*1.16*sqrt( mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) ))) continue; 
			if (ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom())) < 0.005) continue;
			//apply trigger cut: 36 for upsilon, 73 for Jpsi
			//if ((trigger&triggerCuts_)==0) continue;
			//if (dimuonCand->pt() < 7) ccontinue; //another method: using mumufit_ instead of dimuon			
			nGoodUpsilonCand++;
			pat::CompositeCandidate thisDimuonCand = *dimuonCand;
                        bool dimuon_trigger_matched = TriggerMatch(theTriggerPassed,thisDimuonCand);
                        if (!dimuon_trigger_matched) continue;
			fillUpsilonVector(mumuVertexFitTree,thisDimuonCand,bFieldHandle,bs);
			if (nGoodUpsilonCand==1) fillUpsilonBestVertex(mumuVertexFitTree,thisDimuonCand,bFieldHandle,bs);
			if (best4muonCand_ == false || (best4muonCand_ == true && nGoodUpsilonCand==1)) {
				/*		//4 muon mix
						if(muons_previousEvent.size()>3) { 
				//std::cout<<"-1 event: "<<muons_previousEvent.at(muons_previousEvent.size()-1).size()<<",  -2 event: "<<muons_previousEvent.at(muons_previousEvent.size()-2).size()<<",  -3 event: "<<muons_previousEvent.at(muons_previousEvent.size()-3).size()<<std::endl;
				if (muons_previousEvent.at(muons_previousEvent.size()-3).size() > 0) fourMuonMixFit(thisDimuonCand, muons, muons_previousEvent.at(muons_previousEvent.size()-3), bFieldHandle, bs, thePrimaryV);
				if (muons_previousEvent.at(muons_previousEvent.size()-3).size() > 0) fourMuonMixFit(thisDimuonCand, muons, muons_previousEvent.at(muons_previousEvent.size()-3), muons_previousEvent.at(muons_previousEvent.size()-2), bFieldHandle, bs, thePrimaryV);}
				*/		//4 muon fit
                                
				fourMuonFit(thisDimuonCand, muons, bFieldHandle, bs, thePrimaryV);
			}

			if (fabs(mumu_vFit_noMC->currentState().mass()-9.46)<bestYMass) {
				bestYMass=fabs(mumu_vFit_noMC->currentState().mass()-9.46);
				bestVertex_and_bestYMass = 0;
				if (nGoodUpsilonCand==1) bestVertex_and_bestYMass = 1;
				DimuonCand_bestYMass = *dimuonCand;
				fillUpsilonBestMass(mumuVertexFitTree,DimuonCand_bestYMass,bFieldHandle,bs);
			}
		} //end of Upsilon loop

		//if (nGoodUpsilonCand>0 && muons_previousEvent_bestYMass.size() > 0) fourMuonMixFit_bestYMass(DimuonCand_bestYMass, muons, muons_previousEvent_bestYMass, bFieldHandle, bs, thePrimaryV);

		//if (nGoodUpsilonCand>0) fourMuonFit_bestYMass(DimuonCand_bestYMass, muons, bFieldHandle, bs, thePrimaryV);
	}

	if (nGoodUpsilonCand>0) onia_tree->Fill();

	}		//end run number selection
}

void  MuMuGammaRootupler::fillUpsilonVector(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs) {
	mumuVertexFitTree->movePointerToTheTop();     
	RefCountedKinematicParticle mumu_vFit_noMC = mumuVertexFitTree->currentParticle();    
	RefCountedKinematicVertex mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex(); //fitted vertex is same as the commonVertex in the Onia2MuMu skim
	//KinematicParameters mymumupara=  mumu_vFit_noMC->currentState().kinematicParameters();
	//cout<<"mymumupara px="<<mymumupara.momentum().x()<<",py="<<mymumupara.momentum().y()<<", m="<<mumu_vFit_noMC->currentState().mass()<<endl;      
	//float mymumuonlyctau=GetcTau(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV,beamSpot);
	//float mymumuonlyctauerr=GetcTauErr(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV,beamSpot);     
	//cout<<"mymumuonlyctau="<<mymumuonlyctau<<endl;
	//v_mumufit_Ctau->push_back( mymumuonlyctau );
	//v_mumufit_Ctauerr->push_back( mymumuonlyctauerr );
	v_mumufit_Mass.push_back( mumu_vFit_noMC->currentState().mass());
	v_mumufit_MassErr.push_back( sqrt( mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) ) ) ; 
	v_mumufit_VtxCL.push_back( ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom())) );
	v_mumufit_VtxCL2.push_back( mumu_vFit_vertex_noMC->chiSquared() );
	v_mumufit_DecayVtxX.push_back( mumu_vFit_vertex_noMC->position().x() );
	v_mumufit_DecayVtxY.push_back( mumu_vFit_vertex_noMC->position().y() );
	v_mumufit_DecayVtxZ.push_back( mumu_vFit_vertex_noMC->position().z() );
	v_mumufit_DecayVtxXE.push_back( mumu_vFit_vertex_noMC->error().cxx() );
	v_mumufit_DecayVtxYE.push_back( mumu_vFit_vertex_noMC->error().cyy() );
	v_mumufit_DecayVtxZE.push_back( mumu_vFit_vertex_noMC->error().czz() );
	//v_mumufit_p4.SetXYZM( mumu_vFit_noMC->currentState().globalMomentum().x(), mumu_vFit_noMC->currentState().globalMomentum().y(), mumu_vFit_noMC->currentState().globalMomentum().z(), mumufit_Mass ); 
	//v_mumufit_mupIdx->push_back( std::distance(thePATMuonHandle->begin(), iMuon1));
	//v_mumufit_mumIdx->push_back( std::distance(thePATMuonHandle->begin(), iMuon2));
	//std::cout<<"mass0="<<dimuonCand->mass()<<", mass1="<<v_mumufit_Mass<<std::endl;
	//std::cout<<"vProb0="<<vProb<<", vProb1="<<v_mumufit_VtxCL<<std::endl;

	//raw dimuon and muon
	//dimuon_p4.SetPtEtaPhiM(dimuonCand.pt(), dimuonCand.eta(), dimuonCand.phi(), dimuonCand.mass());
	//reco::Candidate::LorentzVector vP = dimuonCand.daughter("muon1")->p4();
	//reco::Candidate::LorentzVector vM = dimuonCand.daughter("muon2")->p4();
	//std::cout<<"muon charge"<<dimuonCand->daughter("muon1")->charge()<<" "<<dimuonCand->daughter("muon2")->charge()<<std::endl;
	//if ( dimuonCand->daughter("muon1")->charge() < 0) {
	// vP = dimuonCand->daughter("muon2")->p4();
	// vM = dimuonCand->daughter("muon1")->p4();
	//}  
	//mu1_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
	//mu2_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
	v_mu1Charge.push_back( dimuonCand.daughter("muon1")->charge() );
	v_mu2Charge.push_back( dimuonCand.daughter("muon2")->charge() );

	reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
	reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
	reco::TrackTransientTrack muon1TTT(muTrack1_ref, &(*bFieldHandle));
	reco::TrackTransientTrack muon2TTT(muTrack2_ref, &(*bFieldHandle));
	v_mu1_d0.push_back( -muon1TTT.dxy(bs));
	v_mu1_d0err.push_back( muon1TTT.d0Error());
	v_mu1_dz.push_back( muon1TTT.dz());
	v_mu1_dzerr.push_back( muon1TTT.dzError());
	v_mu2_d0.push_back( -muon2TTT.dxy(bs));
	v_mu2_d0err.push_back( muon2TTT.d0Error());
	v_mu2_dz.push_back( muon2TTT.dz());
	v_mu2_dzerr.push_back( muon2TTT.dzError());
}


void  MuMuGammaRootupler::fillUpsilonBestVertex(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs) {
	mumuVertexFitTree->movePointerToTheTop();     
	RefCountedKinematicParticle mumu_vFit_noMC = mumuVertexFitTree->currentParticle();    
	RefCountedKinematicVertex mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex(); //fitted vertex is same as the commonVertex in the Onia2MuMu skim
	//KinematicParameters mymumupara=  mumu_vFit_noMC->currentState().kinematicParameters();
	//cout<<"mymumupara px="<<mymumupara.momentum().x()<<",py="<<mymumupara.momentum().y()<<", m="<<mumu_vFit_noMC->currentState().mass()<<endl;      
	//float mymumuonlyctau=GetcTau(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV,beamSpot);
	//float mymumuonlyctauerr=GetcTauErr(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV,beamSpot);     
	//cout<<"mymumuonlyctau="<<mymumuonlyctau<<endl;
	//mumufit_Ctau->push_back( mymumuonlyctau );
	//mumufit_Ctauerr->push_back( mymumuonlyctauerr );
	mumufit_Mass = mumu_vFit_noMC->currentState().mass();
	mumufit_MassErr = sqrt( mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) )  ;
	mumufit_VtxCL = ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom())) ;
	mumufit_VtxCL2 = mumu_vFit_vertex_noMC->chiSquared() ;
	mumufit_DecayVtxX = mumu_vFit_vertex_noMC->position().x() ;
	mumufit_DecayVtxY = mumu_vFit_vertex_noMC->position().y() ;
	mumufit_DecayVtxZ = mumu_vFit_vertex_noMC->position().z() ;
	mumufit_DecayVtxXE = mumu_vFit_vertex_noMC->error().cxx() ;
	mumufit_DecayVtxYE = mumu_vFit_vertex_noMC->error().cyy() ;
	mumufit_DecayVtxZE = mumu_vFit_vertex_noMC->error().czz() ;
	mumufit_p4.SetXYZM( mumu_vFit_noMC->currentState().globalMomentum().x(), mumu_vFit_noMC->currentState().globalMomentum().y(), mumu_vFit_noMC->currentState().globalMomentum().z(), mumufit_Mass ); 
	//mumufit_mupIdx->push_back( std::distance(thePATMuonHandle->begin(), iMuon1));
	//mumufit_mumIdx->push_back( std::distance(thePATMuonHandle->begin(), iMuon2));
	//std::cout<<"mass0="<<dimuonCand->mass()<<", mass1="<<mumufit_Mass<<std::endl;
	//std::cout<<"vProb0="<<vProb<<", vProb1="<<mumufit_VtxCL<<std::endl;

	//raw dimuon and muon
	dimuon_p4.SetPtEtaPhiM(dimuonCand.pt(), dimuonCand.eta(), dimuonCand.phi(), dimuonCand.mass());
	reco::Candidate::LorentzVector vP = dimuonCand.daughter("muon1")->p4();
	reco::Candidate::LorentzVector vM = dimuonCand.daughter("muon2")->p4();
	//std::cout<<"muon charge"<<dimuonCand->daughter("muon1")->charge()<<" "<<dimuonCand->daughter("muon2")->charge()<<std::endl;
	//if ( dimuonCand->daughter("muon1")->charge() < 0) {
	// vP = dimuonCand->daughter("muon2")->p4();
	// vM = dimuonCand->daughter("muon1")->p4();
	//}
	mu1_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
	mu2_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
	mu1Charge = dimuonCand.daughter("muon1")->charge();
	mu2Charge = dimuonCand.daughter("muon2")->charge();

	reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
	reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
	reco::TrackTransientTrack muon1TTT(muTrack1_ref, &(*bFieldHandle));
	reco::TrackTransientTrack muon2TTT(muTrack2_ref, &(*bFieldHandle));
	mu1_d0 = -muon1TTT.dxy(bs);
	mu1_d0err = muon1TTT.d0Error();
	mu1_dz = muon1TTT.dz();
	mu1_dzerr = muon1TTT.dzError();
	mu2_d0 = -muon2TTT.dxy(bs);
	mu2_d0err = muon2TTT.d0Error();
	mu2_dz = muon2TTT.dz();
	mu2_dzerr = muon2TTT.dzError();
}


void  MuMuGammaRootupler::fillUpsilonBestMass(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs) {
	mumuVertexFitTree->movePointerToTheTop();
	RefCountedKinematicParticle mumu_vFit_noMC = mumuVertexFitTree->currentParticle();
	RefCountedKinematicVertex mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex(); //fitted vertex is same as the commonVertex in the Onia2MuMu skim
	//KinematicParameters mymumupara=  mumu_vFit_noMC->currentState().kinematicParameters();
	//cout<<"mymumupara px="<<mymumupara.momentum().x()<<",py="<<mymumupara.momentum().y()<<", m="<<mumu_vFit_noMC->currentState().mass()<<endl;      
	//float mymumuonlyctau=GetcTau(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV,beamSpot);
	//float mymumuonlyctauerr=GetcTauErr(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV,beamSpot);     
	//cout<<"mymumuonlyctau="<<mymumuonlyctau<<endl;
	//mumufit_Ctau->push_back( mymumuonlyctau );
	//mumufit_Ctauerr->push_back( mymumuonlyctauerr );
	mumufit_Mass_bestYMass = mumu_vFit_noMC->currentState().mass();
	mumufit_MassErr_bestYMass = sqrt( mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) )  ;
	mumufit_VtxCL_bestYMass = ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom())) ;
	mumufit_VtxCL2_bestYMass = mumu_vFit_vertex_noMC->chiSquared() ;
	mumufit_DecayVtxX_bestYMass = mumu_vFit_vertex_noMC->position().x() ;
	mumufit_DecayVtxY_bestYMass = mumu_vFit_vertex_noMC->position().y() ;
	mumufit_DecayVtxZ_bestYMass = mumu_vFit_vertex_noMC->position().z() ;
	mumufit_DecayVtxXE_bestYMass = mumu_vFit_vertex_noMC->error().cxx() ;
	mumufit_DecayVtxYE_bestYMass = mumu_vFit_vertex_noMC->error().cyy() ;
	mumufit_DecayVtxZE_bestYMass = mumu_vFit_vertex_noMC->error().czz() ;
	mumufit_p4_bestYMass.SetXYZM( mumu_vFit_noMC->currentState().globalMomentum().x(), mumu_vFit_noMC->currentState().globalMomentum().y(), mumu_vFit_noMC->currentState().globalMomentum().z(), mumufit_Mass_bestYMass );
	//mumufit_mupIdx->push_back( std::distance(thePATMuonHandle->begin(), iMuon1));
	//mumufit_mumIdx->push_back( std::distance(thePATMuonHandle->begin(), iMuon2));
	//std::cout<<"mass0="<<dimuonCand->mass()<<", mass1="<<mumufit_Mass<<std::endl;
	//std::cout<<"vProb0="<<vProb<<", vProb1="<<mumufit_VtxCL<<std::endl;

	//raw dimuon and muon
	dimuon_p4_bestYMass.SetPtEtaPhiM(dimuonCand.pt(), dimuonCand.eta(), dimuonCand.phi(), dimuonCand.mass());
	reco::Candidate::LorentzVector vP = dimuonCand.daughter("muon1")->p4();
	reco::Candidate::LorentzVector vM = dimuonCand.daughter("muon2")->p4();
	//std::cout<<"muon charge"<<dimuonCand->daughter("muon1")->charge()<<" "<<dimuonCand->daughter("muon2")->charge()<<std::endl;
	//if ( dimuonCand->daughter("muon1")->charge() < 0) {
	// vP = dimuonCand->daughter("muon2")->p4();
	// vM = dimuonCand->daughter("muon1")->p4();
	//}
	mu1_p4_bestYMass.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
	mu2_p4_bestYMass.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
	mu1Charge_bestYMass = dimuonCand.daughter("muon1")->charge();
	mu2Charge_bestYMass = dimuonCand.daughter("muon2")->charge();

	reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
	reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
	reco::TrackTransientTrack muon1TTT(muTrack1_ref, &(*bFieldHandle));
	reco::TrackTransientTrack muon2TTT(muTrack2_ref, &(*bFieldHandle));
	mu1_d0_bestYMass = -muon1TTT.dxy(bs);
	mu1_d0err_bestYMass = muon1TTT.d0Error();
	mu1_dz_bestYMass = muon1TTT.dz();
	mu1_dzerr_bestYMass = muon1TTT.dzError();
	mu2_d0_bestYMass = -muon2TTT.dxy(bs);
	mu2_d0err_bestYMass = muon2TTT.d0Error();
	mu2_dz_bestYMass = muon2TTT.dz();
	mu2_dzerr_bestYMass = muon2TTT.dzError();
}


// ------------ method called once each job just before starting event loop  ------------
void MuMuGammaRootupler::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void MuMuGammaRootupler::endJob(const edm::Event & iEvent) {
}


// ------------ method called when ending the processing of a run  ------------
void MuMuGammaRootupler::beginRun(edm::Run & iRun, edm::EventSetup const& iSetup) {

  if(verbose){
    cout<<" New run..."<<endl;
  }
  //--- m_l1GtUtils.getL1GtRunCache(run, iSetup, true, false);
  bool hltConfigChanged;
  bool test = hltConfig_.init(iRun, iSetup, hltName_, hltConfigChanged);
  if (hltConfig_.init(iRun, iSetup, hltName_, hltConfigChanged)) {
    if(verbose){
      std::cout<<" hltConfig_.size() = "<<hltConfig_.size()<<std::endl;
    }
    // check if trigger name in (new) config
    if(verbose){
      if (triggerName_!="@") { // "@" means: analyze all triggers in config
        const unsigned int n(hltConfig_.size());
        const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
        cout<<" triggerIndex = "<<triggerIndex<<endl;
        if (triggerIndex>=n) {
          cout << "HLTEventAnalyzerAOD::beginRun:"
               << " TriggerName " << triggerName_
               << " not available in (new) config!" << endl;
          cout << "Available TriggerNames are: " << endl;
          hltConfig_.dump("Triggers");
        }
      }
      else{
        cout<<" Bad trigger name"<<endl;
      }
    }
  } else {
    cout << "beginRun:"
         << " HLT config extraction failure with process name "
         << hltName_ << endl;

  }     

}

// ------------ method called when starting to processes a run  ------------
 void MuMuGammaRootupler::endRun(edm::Run const &, edm::EventSetup const &) {
 }

// ------------ method called when starting to processes a luminosity block  ------------
void MuMuGammaRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void MuMuGammaRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuMuGammaRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

int MuMuGammaRootupler::looseMuon(edm::View<pat::Muon>::const_iterator rmu) {                                                        
	int goodLooseMuon=0;

	if(  muon::isLooseMuon(*rmu))
	goodLooseMuon = 1;
   return goodLooseMuon;
}

int MuMuGammaRootupler::mediumMuon(edm::View<pat::Muon>::const_iterator rmu) {
	int goodMediumMuon=0;

	if(  muon::isMediumMuon(*rmu))
	//bool goodGlob = rmu->isGlobalMuon() &&
	//	rmu->globalTrack()->normalizedChi2() < 3 &&
	//	rmu->combinedQuality().chi2LocalPosition < 12 &&
	//	rmu->combinedQuality().trkKink < 20;
	//if(  muon::isLooseMuon(*rmu) &&
	//		rmu->innerTrack()->validFraction() > 0.8 &&
	//		muon::segmentCompatibility(*rmu) > (goodGlob ? 0.303 : 0.451)
	//  ) 
	goodMediumMuon=1;

	return goodMediumMuon;
}

int MuMuGammaRootupler::tightMuon(edm::View<pat::Muon>::const_iterator rmu, reco::Vertex vertex) {
	int goodTightMuon=0;

	if( muon::isTightMuon(*rmu,vertex))
	//if( rmu->isGlobalMuon()
	//		&& rmu->isPFMuon()
	//		&& rmu->globalTrack()->normalizedChi2()<10.0
	//		&& rmu->globalTrack()->hitPattern().numberOfValidMuonHits()>0
	//		&& rmu->numberOfMatchedStations()>1
	//		&& fabs(rmu->muonBestTrack()->dxy( vertex.position() ))<0.2
	//		&& fabs(rmu->muonBestTrack()->dz( vertex.position() )) < 0.5
	//		&& rmu->innerTrack()->hitPattern().numberOfValidPixelHits()>0
	//		&& rmu->track()->hitPattern().trackerLayersWithMeasurement()>5
	//  ) 
	goodTightMuon = 1;

	return goodTightMuon;
}

int MuMuGammaRootupler::fourMuonMixFit(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, std::vector<pat::Muon> muons_previous1, std::vector<pat::Muon> muons_previous2, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	int nGoodFourMuonMix=0;
	fourMuFit_VtxProb_mix = -1;
	for (std::vector<pat::Muon>::iterator mu3 = muons_previous1.begin(), mu3end = muons_previous1.end(); mu3 != mu3end; ++mu3){
		for (std::vector<pat::Muon>::iterator mu4 = muons_previous2.begin(), mu4end = muons_previous2.end(); mu4 != mu4end; ++mu4){
			if (mu3->charge() == mu4->charge()) continue;
			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));
			KinematicConstrainedVertexFitter constVertexFitter;
			//fit w/ mass constraint
			//MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
			//RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles,upsilon_mtc);

			//fit w/o any mass constraint
			RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles);

			if(!fourMuTree->isEmpty())
			{
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();
				if (fitFourMu->currentState().isValid() &&
						ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb_mix)
				{ //Get chib         
					fourMuFit_Mass_mix3evts = fitFourMu->currentState().mass();
					fourMuFit_p4_mix3evts.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fourMuFit_Mass_mix3evts);
					fourMuFit_VtxProb_mix3evts = ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom()));
				}
			}
		}
	}
	return nGoodFourMuonMix;
}

int MuMuGammaRootupler::fourMuonMixFit(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, std::vector<pat::Muon> muons_previous, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	int nGoodFourMuonMix=0;
	for (std::vector<pat::Muon>::iterator mu3 = muons_previous.begin(), muend = muons_previous.end(); mu3 != muend; ++mu3){
		for(std::vector<pat::Muon>::iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){
			if (mu3->charge() == mu4->charge()) continue;

			TLorentzVector mu1p4,mu2p4,mu3p4,mu4p4,fourMup4;
			reco::Candidate::LorentzVector v1 = dimuonCand.daughter("muon1")->p4();
			reco::Candidate::LorentzVector v2 = dimuonCand.daughter("muon2")->p4();
			mu1p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
			mu2p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
			mu3p4.SetPtEtaPhiM(mu3->pt(), mu3->eta(), mu3->phi(), mu3->mass());
			mu4p4.SetPtEtaPhiM(mu4->pt(), mu4->eta(), mu4->phi(), mu4->mass());
			fourMup4 = mu1p4 + mu2p4 + mu3p4 + mu4p4;
			fourMuFit_Mass_allComb_mix.push_back(fourMup4.M());
		}
	}
	//std::cout<<"previousSize="<<muons_previous.size()<<", thisEventSize="<<muons->size()<<std::endl;
	int muons_previous_originalSize=muons_previous.size();

	for (edm::View<pat::Muon>::const_iterator mu = muons->begin(); mu !=  muons->end(); ++mu){
		if (mu->pt()<2.0 || fabs(mu->eta())>2.4) continue;
		if (mu-muons->begin() == dimuonCand.userInt("mu1Index"))  continue; 
		if (mu-muons->begin() == dimuonCand.userInt("mu2Index"))  continue;
		reco::GenParticleRef genMu;
		if (isMC_) genMu = mu->genParticleRef();
		if (!isMC_ || (isMC_ && !genMu.isNonnull())) {
			muons_previous.push_back(*mu);
		}
	}

	fourMuFit_VtxProb_mix = -1;
	//std::cout<<"combinedSize="<<muons_previous.size()<<std::endl;
	for (std::vector<pat::Muon>::iterator mu3 = muons_previous.begin(), muend = muons_previous.end(); mu3-muons_previous.begin() < muons_previous_originalSize; ++mu3){
		for(std::vector<pat::Muon>::iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){
			if (mu3->charge() == mu4->charge()) continue;

			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));
			KinematicConstrainedVertexFitter constVertexFitter;
			//fit w/ mass constraint
			//MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
			//RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles,upsilon_mtc);

			//fit w/o any mass constraint
			RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles);

			if(!fourMuTree->isEmpty())
			{
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();
				if (fitFourMu->currentState().isValid() &&
						ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb_mix)
				{ //Get chib         
					nGoodFourMuonMix = 1;
					fourMuFit_Mass_mix = fitFourMu->currentState().mass();
					fourMuFit_MassErr_mix = sqrt(fitFourMu->currentState().kinematicParametersError().matrix()(6,6));
					fourMuFit_p4_mix.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fourMuFit_Mass_mix);
					fourMuFit_VtxX_mix = FourMuDecayVertex->position().x();
					fourMuFit_VtxY_mix = FourMuDecayVertex->position().y();
					fourMuFit_VtxZ_mix = FourMuDecayVertex->position().z();
					fourMuFit_VtxProb_mix = ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom()));
					fourMuFit_Chi2_mix = FourMuDecayVertex->chiSquared();
					fourMuFit_ndof_mix = FourMuDecayVertex->degreesOfFreedom();
					if (mu4-muons_previous.begin() >= muons_previous_originalSize) fourMuFit_3plus1_mix = 1; 
					else fourMuFit_3plus1_mix = 0;

					//get first muon
					bool child = fourMuTree->movePointerToTheFirstChild();
					RefCountedKinematicParticle fitMu1 = fourMuTree->currentParticle();
					if(!child) break;
					double mu1M_fit_mix = fitMu1->currentState().mass();
					double mu1Px_fit_mix = fitMu1->currentState().kinematicParameters().momentum().x();
					double mu1Py_fit_mix = fitMu1->currentState().kinematicParameters().momentum().y();
					double mu1Pz_fit_mix = fitMu1->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu1p4_mix.SetXYZM( mu1Px_fit_mix, mu1Py_fit_mix, mu1Pz_fit_mix, mu1M_fit_mix );

					//get second muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu2 = fourMuTree->currentParticle();
					if(!child) break;
					double mu2M_fit_mix = fitMu2->currentState().mass();
					double mu2Px_fit_mix = fitMu2->currentState().kinematicParameters().momentum().x();
					double mu2Py_fit_mix = fitMu2->currentState().kinematicParameters().momentum().y();
					double mu2Pz_fit_mix = fitMu2->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu2p4_mix.SetXYZM( mu2Px_fit_mix, mu2Py_fit_mix, mu2Pz_fit_mix, mu2M_fit_mix );

					//get third muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu3 = fourMuTree->currentParticle();
					if(!child) break;
					double mu3M_fit_mix = fitMu3->currentState().mass();
					double mu3Px_fit_mix = fitMu3->currentState().kinematicParameters().momentum().x();
					double mu3Py_fit_mix = fitMu3->currentState().kinematicParameters().momentum().y();
					double mu3Pz_fit_mix = fitMu3->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu3p4_mix.SetXYZM( mu3Px_fit_mix, mu3Py_fit_mix, mu3Pz_fit_mix, mu3M_fit_mix );

					//get fourth muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu4 = fourMuTree->currentParticle();
					if(!child) break;
					double mu4M_fit_mix = fitMu4->currentState().mass();
					double mu4Px_fit_mix = fitMu4->currentState().kinematicParameters().momentum().x();
					double mu4Py_fit_mix = fitMu4->currentState().kinematicParameters().momentum().y();
					double mu4Pz_fit_mix = fitMu4->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu4p4_mix.SetXYZM( mu4Px_fit_mix, mu4Py_fit_mix, mu4Pz_fit_mix, mu4M_fit_mix );

					reco::Candidate::LorentzVector v3 = mu3->p4();
					reco::Candidate::LorentzVector v4 = mu4->p4();
					mu3_p4_mix.SetPtEtaPhiM(v3.pt(),v3.eta(),v3.phi(),v3.mass());
					mu4_p4_mix.SetPtEtaPhiM(v4.pt(),v4.eta(),v4.phi(),v4.mass());
					reco::TrackTransientTrack muon3TTT(muTrack3_ref, &(*bFieldHandle));
					reco::TrackTransientTrack muon4TTT(muTrack4_ref, &(*bFieldHandle));
					mu3_d0_mix = -muon3TTT.dxy(bs);
					mu3_d0err_mix = muon3TTT.d0Error();
					mu3_dz_mix = muon3TTT.dz();
					mu3_dzerr_mix = muon3TTT.dzError();
					mu4_d0_mix = -muon4TTT.dxy(bs);
					mu4_d0err_mix = muon4TTT.d0Error();
					mu4_dz_mix = muon4TTT.dz();
					mu4_dzerr_mix = muon4TTT.dzError();
					mu3Charge_mix = mu3->charge();
					mu4Charge_mix = mu4->charge();
				}
			}
		}
	}
	return nGoodFourMuonMix;
}


void MuMuGammaRootupler::fourMuonFit(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	std::vector<pat::Muon> theRestMuons;
	//fourMuFit_VtxProb = -1;
	//std::cout<<"mu1Index="<<mu1Index<<", mu2Index="<<mu2Index<<std::endl; 
	for (edm::View<pat::Muon>::const_iterator mu3 = muons->begin(), muend = muons->end(); mu3 != muend; ++mu3){
		if (mu3->pt()<2.0 || fabs(mu3->eta())>2.4)  continue;
		if (mu3-muons->begin() == dimuonCand.userInt("mu1Index"))  continue;
		if (mu3-muons->begin() == dimuonCand.userInt("mu2Index"))  continue;
		reco::GenParticleRef genMu3;
		if (isMC_) genMu3 = mu3->genParticleRef();
		if (!isMC_ || (isMC_ && !genMu3.isNonnull())) theRestMuons.push_back(*mu3);
		//std::cout<<"fill vector: "<<mu3->pt()<<std::endl;

		for(edm::View<pat::Muon>::const_iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){
			if (mu4->pt()<2.0 || fabs(mu4->eta())>2.4)  continue;
			if (mu4-muons->begin() == dimuonCand.userInt("mu1Index")) continue;
			if (mu4-muons->begin() == dimuonCand.userInt("mu2Index")) continue;
			reco::GenParticleRef genMu4;
			if (isMC_) genMu4 = mu4->genParticleRef();

			if (mu3->charge() == mu4->charge()) continue;
			/*if ( (tightMuon(muons->begin()+mu1Index)+
			  tightMuon(muons->begin()+mu2Index)+
			  tightMuon(mu3)+
			  tightMuon(mu4)) < 2 
			  ) continue;*/

			if (!isMC_ || (isMC_ && !genMu3.isNonnull() && !genMu4.isNonnull())) {
				TLorentzVector mu1p4,mu2p4,mu3p4,mu4p4,fourMup4;
				reco::Candidate::LorentzVector v1 = dimuonCand.daughter("muon1")->p4();
				reco::Candidate::LorentzVector v2 = dimuonCand.daughter("muon2")->p4();
				mu1p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
				mu2p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
				mu3p4.SetPtEtaPhiM(mu3->pt(), mu3->eta(), mu3->phi(), mu3->mass());
				mu4p4.SetPtEtaPhiM(mu4->pt(), mu4->eta(), mu4->phi(), mu4->mass());   
                                bool Rest_Muon_trigger_Matched = TriggerMatch_restMuons(mu3p4,mu4p4);
                                if (!Rest_Muon_trigger_Matched) continue;
				fourMup4 = mu1p4 + mu2p4 + mu3p4 + mu4p4;
				fourMuFit_Mass_allComb.push_back(fourMup4.M());
			}

			//std::cout<<"found good mu3mu4: "<<mu3->pt()<<" "<<mu4->pt()<<", mu1: "<<muon1TT.track().pt()<<", mu2: "<<muon2TT.track().pt()<<std::endl;
			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));

			KinematicConstrainedVertexFitter constVertexFitter;
			//fit w/ mass constraint
			//MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
			//RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles,upsilon_mtc);

			//fit w/o any mass constraint
			RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles);

			if(!fourMuTree->isEmpty())
			{
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();

				if (fitFourMu->currentState().isValid()
						// && ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb)
					)
					{ //Get chib         
						fourMuFit_Mass.push_back(fitFourMu->currentState().mass());
						fourMuFit_MassErr.push_back(sqrt(fitFourMu->currentState().kinematicParametersError().matrix()(6,6)));
						TLorentzVector p4;
						p4.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fitFourMu->currentState().mass());
						fourMuFit_Pt.push_back(p4.Pt());
						fourMuFit_Eta.push_back(p4.Eta());
						fourMuFit_Phi.push_back(p4.Phi());
						fourMuFit_VtxX.push_back(FourMuDecayVertex->position().x());
						fourMuFit_VtxY.push_back(FourMuDecayVertex->position().y());
						fourMuFit_VtxZ.push_back(FourMuDecayVertex->position().z());
						fourMuFit_VtxProb.push_back(ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())));
						fourMuFit_Chi2.push_back(FourMuDecayVertex->chiSquared());
						fourMuFit_ndof.push_back(FourMuDecayVertex->degreesOfFreedom());

						//get first muon
						bool child = fourMuTree->movePointerToTheFirstChild();
						RefCountedKinematicParticle fitMu1 = fourMuTree->currentParticle();
						if(!child) break;
						p4.SetXYZM( fitMu1->currentState().kinematicParameters().momentum().x(), fitMu1->currentState().kinematicParameters().momentum().y(), fitMu1->currentState().kinematicParameters().momentum().z(), fitMu1->currentState().mass() );
						fourMuFit_mu1Pt.push_back(p4.Pt());
						fourMuFit_mu1Eta.push_back(p4.Eta());
						fourMuFit_mu1Phi.push_back(p4.Phi());
						fourMuFit_mu1E.push_back(p4.E());

						//get second muon
						child = fourMuTree->movePointerToTheNextChild();
						RefCountedKinematicParticle fitMu2 = fourMuTree->currentParticle();
						if(!child) break;
						p4.SetXYZM( fitMu2->currentState().kinematicParameters().momentum().x(), fitMu2->currentState().kinematicParameters().momentum().y(), fitMu2->currentState().kinematicParameters().momentum().z(), fitMu2->currentState().mass() );
						fourMuFit_mu2Pt.push_back(p4.Pt());
						fourMuFit_mu2Eta.push_back(p4.Eta());
						fourMuFit_mu2Phi.push_back(p4.Phi());
						fourMuFit_mu2E.push_back(p4.E());

						//get third muon
						child = fourMuTree->movePointerToTheNextChild();
						RefCountedKinematicParticle fitMu3 = fourMuTree->currentParticle();
						if(!child) break;
						p4.SetXYZM( fitMu3->currentState().kinematicParameters().momentum().x(), fitMu3->currentState().kinematicParameters().momentum().y(), fitMu3->currentState().kinematicParameters().momentum().z(), fitMu3->currentState().mass() );
						fourMuFit_mu3Pt.push_back(p4.Pt());
						fourMuFit_mu3Eta.push_back(p4.Eta());
						fourMuFit_mu3Phi.push_back(p4.Phi());
						fourMuFit_mu3E.push_back(p4.E());

						//get fourth muon
						child = fourMuTree->movePointerToTheNextChild();
						RefCountedKinematicParticle fitMu4 = fourMuTree->currentParticle();
						if(!child) break;
						p4.SetXYZM( fitMu4->currentState().kinematicParameters().momentum().x(), fitMu4->currentState().kinematicParameters().momentum().y(), fitMu4->currentState().kinematicParameters().momentum().z(), fitMu4->currentState().mass() );
						fourMuFit_mu4Pt.push_back(p4.Pt());
						fourMuFit_mu4Eta.push_back(p4.Eta());
						fourMuFit_mu4Phi.push_back(p4.Phi());
						fourMuFit_mu4E.push_back(p4.E());

						mu3_Pt.push_back(mu3->pt());
						mu3_Eta.push_back(mu3->eta());
						mu3_Phi.push_back(mu3->phi());
						mu3_E.push_back(mu3->energy());
						mu4_Pt.push_back(mu4->pt());
						mu4_Eta.push_back(mu4->eta());
						mu4_Phi.push_back(mu4->phi());
						mu4_E.push_back(mu4->energy());
						reco::TrackTransientTrack muon3TTT(muTrack3_ref, &(*bFieldHandle));
						reco::TrackTransientTrack muon4TTT(muTrack4_ref, &(*bFieldHandle));
						mu3_d0.push_back(-muon3TTT.dxy(bs));
						mu3_d0err.push_back(muon3TTT.d0Error());
						mu3_dz.push_back(muon3TTT.dz());
						mu3_dzerr.push_back(muon3TTT.dzError());
						mu4_d0.push_back(-muon4TTT.dxy(bs));
						mu4_d0err.push_back(muon4TTT.d0Error());
						mu4_dz.push_back(muon4TTT.dz());
						mu4_dzerr.push_back(muon4TTT.dzError());
						mu3Charge.push_back(mu3->charge());
						mu4Charge.push_back(mu4->charge());
						mu1_Tight.push_back(tightMuon(muons->begin()+dimuonCand.userInt("mu1Index"),thePrimaryV));
						mu2_Tight.push_back(tightMuon(muons->begin()+dimuonCand.userInt("mu2Index"),thePrimaryV));
						mu3_Tight.push_back(tightMuon(mu3,thePrimaryV));
						mu4_Tight.push_back(tightMuon(mu4,thePrimaryV));
						mu1_Medium.push_back(mediumMuon(muons->begin()+dimuonCand.userInt("mu1Index")));
						mu2_Medium.push_back(mediumMuon(muons->begin()+dimuonCand.userInt("mu2Index")));
						mu3_Medium.push_back(mediumMuon(mu3));
						mu4_Medium.push_back(mediumMuon(mu4));
						mu1_Loose.push_back(looseMuon(muons->begin()+dimuonCand.userInt("mu1Index")));
						mu2_Loose.push_back(looseMuon(muons->begin()+dimuonCand.userInt("mu2Index")));
						mu3_Loose.push_back(looseMuon(mu3));
						mu4_Loose.push_back(looseMuon(mu4));

						if (isMC_) {
							reco::GenParticleRef genMu1 = (muons->begin()+dimuonCand.userInt("mu1Index"))->genParticleRef();
							reco::GenParticleRef genMu2 = (muons->begin()+dimuonCand.userInt("mu2Index"))->genParticleRef();
							if (genMu1.isNonnull() ) mu1_pdgID.push_back(genMu1->pdgId());
							else mu1_pdgID.push_back(0);
							if (genMu2.isNonnull() ) mu2_pdgID.push_back(genMu2->pdgId());
							else mu2_pdgID.push_back(0);
							//if (genMu1->motherRef()==genMu2->motherRef()) std::cout<<"genMu1->motherRef()->pdgId()="<<genMu1->motherRef()->pdgId()<<std::endl;
							//else std::cout<<"genMu1->motherRef()->pdgId()="<<genMu1->motherRef()->pdgId()<<", genMu2->motherRef()->pdgId()="<<genMu2->motherRef()->pdgId()
							// <<"genMu1->motherRef()->motherRef()->pdgId()="<<genMu1->motherRef()->motherRef()->pdgId()<<", genMu2->motherRef()->motherRef()->pdgId()="<<genMu2->motherRef()->motherRef()->pdgId()<<std::endl;

							//reco::GenParticleRef genMu3 = mu3->genParticleRef();
							//reco::GenParticleRef genMu4 = mu4->genParticleRef();
							if (genMu3.isNonnull() ){
								mu3_pdgID.push_back(genMu3->pdgId());
								/*
									if (genMu3->numberOfMothers()>0){ 
									reco::GenParticleRef mom3 = genMu3->motherRef();
									if (mom3.isNonnull()) { 
									std::cout<<""<<"mom pdgID= "<<mom3->pdgId()<<std::endl;
									if (mom3==genMu1->motherRef()) std::cout<<"same mother"<<std::endl;
									}    
									else std::cout<<"mom non"<<std::endl;
									}    
									else std::cout<<"# mom = 0"<<std::endl;
									*/
							}
							else mu3_pdgID.push_back(0);
							if (genMu4.isNonnull() ) mu4_pdgID.push_back(genMu4->pdgId());
							else mu4_pdgID.push_back(0);
						}
					}
			}
		}
	}
	muons_previousEvent.push_back(theRestMuons);
}


int MuMuGammaRootupler::fourMuonMixFit_bestYMass(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, std::vector<pat::Muon> muons_previous, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	int nGoodFourMuonMix=0;
	for (std::vector<pat::Muon>::iterator mu3 = muons_previous.begin(), muend = muons_previous.end(); mu3 != muend; ++mu3){
		for(std::vector<pat::Muon>::iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){
			if (mu3->charge() == mu4->charge()) continue;

			TLorentzVector mu1p4,mu2p4,mu3p4,mu4p4,fourMup4;
			reco::Candidate::LorentzVector v1 = dimuonCand.daughter("muon1")->p4();
			reco::Candidate::LorentzVector v2 = dimuonCand.daughter("muon2")->p4();
			mu1p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
			mu2p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
			mu3p4.SetPtEtaPhiM(mu3->pt(), mu3->eta(), mu3->phi(), mu3->mass());
			mu4p4.SetPtEtaPhiM(mu4->pt(), mu4->eta(), mu4->phi(), mu4->mass());
			fourMup4 = mu1p4 + mu2p4 + mu3p4 + mu4p4;
			fourMuFit_Mass_allComb_mix_bestYMass.push_back(fourMup4.M());
		}
	}
	//std::cout<<"previousSize="<<muons_previous.size()<<", thisEventSize="<<muons->size()<<std::endl;
	int muons_previous_originalSize=muons_previous.size();

	for (edm::View<pat::Muon>::const_iterator mu = muons->begin(); mu !=  muons->end(); ++mu){
		if (mu->pt()<2.0 || fabs(mu->eta())>2.4) continue;
		if (mu-muons->begin() == dimuonCand.userInt("mu1Index"))  continue;
		if (mu-muons->begin() == dimuonCand.userInt("mu2Index"))  continue;
		reco::GenParticleRef genMu;
		if (isMC_) genMu = mu->genParticleRef();
		if (!isMC_ || (isMC_ && !genMu.isNonnull())) {
			muons_previous.push_back(*mu);
		}
	}

	fourMuFit_VtxProb_mix_bestYMass = -1;
	//std::cout<<"combinedSize="<<muons_previous.size()<<std::endl;
	for (std::vector<pat::Muon>::iterator mu3 = muons_previous.begin(), muend = muons_previous.end(); mu3-muons_previous.begin() < muons_previous_originalSize; ++mu3){
		for(std::vector<pat::Muon>::iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){
			if (mu3->charge() == mu4->charge()) continue;

			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));
			KinematicConstrainedVertexFitter constVertexFitter;
			//fit w/ mass constraint
			//MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
			//RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles,upsilon_mtc);

			//fit w/o any mass constraint
			RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles);

			if(!fourMuTree->isEmpty())
			{
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();
				if (fitFourMu->currentState().isValid() &&
						ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb_mix_bestYMass)
				{ //Get chib         
					nGoodFourMuonMix = 1;
					fourMuFit_Mass_mix_bestYMass = fitFourMu->currentState().mass();
					fourMuFit_MassErr_mix_bestYMass = sqrt(fitFourMu->currentState().kinematicParametersError().matrix()(6,6));
					fourMuFit_p4_mix_bestYMass.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fourMuFit_Mass_mix_bestYMass);
					fourMuFit_VtxX_mix_bestYMass = FourMuDecayVertex->position().x();
					fourMuFit_VtxY_mix_bestYMass = FourMuDecayVertex->position().y();
					fourMuFit_VtxZ_mix_bestYMass = FourMuDecayVertex->position().z();
					fourMuFit_VtxProb_mix_bestYMass = ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom()));
					fourMuFit_Chi2_mix_bestYMass = FourMuDecayVertex->chiSquared();
					fourMuFit_ndof_mix_bestYMass = FourMuDecayVertex->degreesOfFreedom();
					if (mu4-muons_previous.begin() >= muons_previous_originalSize) fourMuFit_3plus1_mix_bestYMass = 1;
					else fourMuFit_3plus1_mix_bestYMass = 0;

					//get first muon
					bool child = fourMuTree->movePointerToTheFirstChild();
					RefCountedKinematicParticle fitMu1 = fourMuTree->currentParticle();
					if(!child) break;
					double mu1M_fit_mix = fitMu1->currentState().mass();
					double mu1Px_fit_mix = fitMu1->currentState().kinematicParameters().momentum().x();
					double mu1Py_fit_mix = fitMu1->currentState().kinematicParameters().momentum().y();
					double mu1Pz_fit_mix = fitMu1->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu1p4_mix_bestYMass.SetXYZM( mu1Px_fit_mix, mu1Py_fit_mix, mu1Pz_fit_mix, mu1M_fit_mix );
					//get second muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu2 = fourMuTree->currentParticle();
					if(!child) break;
					double mu2M_fit_mix = fitMu2->currentState().mass();
					double mu2Px_fit_mix = fitMu2->currentState().kinematicParameters().momentum().x();
					double mu2Py_fit_mix = fitMu2->currentState().kinematicParameters().momentum().y();
					double mu2Pz_fit_mix = fitMu2->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu2p4_mix_bestYMass.SetXYZM( mu2Px_fit_mix, mu2Py_fit_mix, mu2Pz_fit_mix, mu2M_fit_mix );

					//get third muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu3 = fourMuTree->currentParticle();
					if(!child) break;
					double mu3M_fit_mix = fitMu3->currentState().mass();
					double mu3Px_fit_mix = fitMu3->currentState().kinematicParameters().momentum().x();
					double mu3Py_fit_mix = fitMu3->currentState().kinematicParameters().momentum().y();
					double mu3Pz_fit_mix = fitMu3->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu3p4_mix_bestYMass.SetXYZM( mu3Px_fit_mix, mu3Py_fit_mix, mu3Pz_fit_mix, mu3M_fit_mix );

					//get fourth muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu4 = fourMuTree->currentParticle();
					if(!child) break;
					double mu4M_fit_mix = fitMu4->currentState().mass();
					double mu4Px_fit_mix = fitMu4->currentState().kinematicParameters().momentum().x();
					double mu4Py_fit_mix = fitMu4->currentState().kinematicParameters().momentum().y();
					double mu4Pz_fit_mix = fitMu4->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu4p4_mix_bestYMass.SetXYZM( mu4Px_fit_mix, mu4Py_fit_mix, mu4Pz_fit_mix, mu4M_fit_mix );

					reco::Candidate::LorentzVector v3 = mu3->p4();
					reco::Candidate::LorentzVector v4 = mu4->p4();
					mu3_p4_mix_bestYMass.SetPtEtaPhiM(v3.pt(),v3.eta(),v3.phi(),v3.mass());
					mu4_p4_mix_bestYMass.SetPtEtaPhiM(v4.pt(),v4.eta(),v4.phi(),v4.mass());
					reco::TrackTransientTrack muon3TTT(muTrack3_ref, &(*bFieldHandle));
					reco::TrackTransientTrack muon4TTT(muTrack4_ref, &(*bFieldHandle));
					mu3_d0_mix_bestYMass = -muon3TTT.dxy(bs);
					mu3_d0err_mix_bestYMass = muon3TTT.d0Error();
					mu3_dz_mix_bestYMass = muon3TTT.dz();
					mu3_dzerr_mix_bestYMass = muon3TTT.dzError();
					mu4_d0_mix_bestYMass = -muon4TTT.dxy(bs);
					mu4_d0err_mix_bestYMass = muon4TTT.d0Error();
					mu4_dz_mix_bestYMass = muon4TTT.dz();
					mu4_dzerr_mix_bestYMass = muon4TTT.dzError();
					mu3Charge_mix_bestYMass = mu3->charge();
					mu4Charge_mix_bestYMass = mu4->charge();
				}
			}
		}
	}
	return nGoodFourMuonMix;
}


void MuMuGammaRootupler::fourMuonFit_bestYMass(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	muons_previousEvent_bestYMass.clear();
	fourMuFit_VtxProb_bestYMass = -1;
	//std::cout<<"mu1Index="<<mu1Index<<", mu2Index="<<mu2Index<<std::endl; 
	for (edm::View<pat::Muon>::const_iterator mu3 = muons->begin(), muend = muons->end(); mu3 != muend; ++mu3){
		if (mu3->pt()<2.0 || fabs(mu3->eta())>2.4)  continue;
		if (mu3-muons->begin() == dimuonCand.userInt("mu1Index"))  continue;
		if (mu3-muons->begin() == dimuonCand.userInt("mu2Index"))  continue;
		reco::GenParticleRef genMu3;
		if (isMC_) genMu3 = mu3->genParticleRef();
		if (!isMC_ || (isMC_ && !genMu3.isNonnull())) muons_previousEvent_bestYMass.push_back(*mu3);
		//std::cout<<"fill vector: "<<mu3->pt()<<std::endl;

		for(edm::View<pat::Muon>::const_iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){
			if (mu4->pt()<2.0 || fabs(mu4->eta())>2.4)  continue;
			if (mu4-muons->begin() == dimuonCand.userInt("mu1Index")) continue;
			if (mu4-muons->begin() == dimuonCand.userInt("mu2Index")) continue;
			reco::GenParticleRef genMu4;
			if (isMC_) genMu4 = mu4->genParticleRef();

			if (mu3->charge() == mu4->charge()) continue;
			/*if ( (tightMuon(muons->begin()+mu1Index)+
			  tightMuon(muons->begin()+mu2Index)+
			  tightMuon(mu3)+
			  tightMuon(mu4)) < 2 
			  ) continue;*/

			if (!isMC_ || (isMC_ && !genMu3.isNonnull() && !genMu4.isNonnull())) {
				TLorentzVector mu1p4,mu2p4,mu3p4,mu4p4,fourMup4;
				reco::Candidate::LorentzVector v1 = dimuonCand.daughter("muon1")->p4();
				reco::Candidate::LorentzVector v2 = dimuonCand.daughter("muon2")->p4();
				mu1p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
				mu2p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
				mu3p4.SetPtEtaPhiM(mu3->pt(), mu3->eta(), mu3->phi(), mu3->mass());
				mu4p4.SetPtEtaPhiM(mu4->pt(), mu4->eta(), mu4->phi(), mu4->mass());
				fourMup4 = mu1p4 + mu2p4 + mu3p4 + mu4p4;
				fourMuFit_Mass_allComb_bestYMass.push_back(fourMup4.M());
			}

			//std::cout<<"found good mu3mu4: "<<mu3->pt()<<" "<<mu4->pt()<<", mu1: "<<muon1TT.track().pt()<<", mu2: "<<muon2TT.track().pt()<<std::endl;
			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));

			KinematicConstrainedVertexFitter constVertexFitter;
			//fit w/ mass constraint
			//MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
			//RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles,upsilon_mtc);

			//fit w/o any mass constraint
			RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles);

			if(!fourMuTree->isEmpty())
			{
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();

				if (fitFourMu->currentState().isValid() &&
						ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb_bestYMass)
				{ //Get chib         
					fourMuFit_Mass_bestYMass = fitFourMu->currentState().mass();
					fourMuFit_MassErr_bestYMass = sqrt(fitFourMu->currentState().kinematicParametersError().matrix()(6,6));
					fourMuFit_p4_bestYMass.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fourMuFit_Mass_bestYMass);
					fourMuFit_VtxX_bestYMass = FourMuDecayVertex->position().x();
					fourMuFit_VtxY_bestYMass = FourMuDecayVertex->position().y();
					fourMuFit_VtxZ_bestYMass = FourMuDecayVertex->position().z();
					fourMuFit_VtxProb_bestYMass = ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom()));
					fourMuFit_Chi2_bestYMass = FourMuDecayVertex->chiSquared();
					fourMuFit_ndof_bestYMass = FourMuDecayVertex->degreesOfFreedom();

					//get first muon
					bool child = fourMuTree->movePointerToTheFirstChild();
					RefCountedKinematicParticle fitMu1 = fourMuTree->currentParticle();
					if(!child) break;
					double mu1M_fit = fitMu1->currentState().mass();
					double mu1Px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
					double mu1Py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
					double mu1Pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu1p4_bestYMass.SetXYZM( mu1Px_fit, mu1Py_fit, mu1Pz_fit, mu1M_fit );

					//get second muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu2 = fourMuTree->currentParticle();
					if(!child) break;
					double mu2M_fit = fitMu2->currentState().mass();
					double mu2Px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
					double mu2Py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
					double mu2Pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu2p4_bestYMass.SetXYZM( mu2Px_fit, mu2Py_fit, mu2Pz_fit, mu2M_fit );

					//get third muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu3 = fourMuTree->currentParticle();
					if(!child) break;
					double mu3M_fit = fitMu3->currentState().mass();
					double mu3Px_fit = fitMu3->currentState().kinematicParameters().momentum().x();
					double mu3Py_fit = fitMu3->currentState().kinematicParameters().momentum().y();
					double mu3Pz_fit = fitMu3->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu3p4_bestYMass.SetXYZM( mu3Px_fit, mu3Py_fit, mu3Pz_fit, mu3M_fit );

					//get fourth muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu4 = fourMuTree->currentParticle();
					if(!child) break;
					double mu4M_fit = fitMu4->currentState().mass();
					double mu4Px_fit = fitMu4->currentState().kinematicParameters().momentum().x();
					double mu4Py_fit = fitMu4->currentState().kinematicParameters().momentum().y();
					double mu4Pz_fit = fitMu4->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu4p4_bestYMass.SetXYZM( mu4Px_fit, mu4Py_fit, mu4Pz_fit, mu4M_fit );

					//std::cout<<fourMuFit_Mass<<" "<<(fourMuFit_mu1p4+fourMuFit_mu2p4+fourMuFit_mu3p4+fourMuFit_mu4p4).M()<<std::endl;

					reco::Candidate::LorentzVector v3 = mu3->p4();
					reco::Candidate::LorentzVector v4 = mu4->p4();
					mu3_p4_bestYMass.SetPtEtaPhiM(v3.pt(),v3.eta(),v3.phi(),v3.mass());
					mu4_p4_bestYMass.SetPtEtaPhiM(v4.pt(),v4.eta(),v4.phi(),v4.mass());
					reco::TrackTransientTrack muon3TTT(muTrack3_ref, &(*bFieldHandle));
					reco::TrackTransientTrack muon4TTT(muTrack4_ref, &(*bFieldHandle));
					mu3_d0_bestYMass = -muon3TTT.dxy(bs);
					mu3_d0err_bestYMass = muon3TTT.d0Error();
					mu3_dz_bestYMass = muon3TTT.dz();
					mu3_dzerr_bestYMass = muon3TTT.dzError();
					mu4_d0_bestYMass = -muon4TTT.dxy(bs);
					mu4_d0err_bestYMass = muon4TTT.d0Error();
					mu4_dz_bestYMass = muon4TTT.dz();
					mu4_dzerr_bestYMass = muon4TTT.dzError();
					mu3Charge_bestYMass = mu3->charge();
					mu4Charge_bestYMass = mu4->charge();

					mu1_Tight_bestYMass = tightMuon(muons->begin()+dimuonCand.userInt("mu1Index"),thePrimaryV);
					mu2_Tight_bestYMass = tightMuon(muons->begin()+dimuonCand.userInt("mu2Index"),thePrimaryV);
					mu3_Tight_bestYMass = tightMuon(mu3,thePrimaryV);
					mu4_Tight_bestYMass = tightMuon(mu4,thePrimaryV);

					if (isMC_) {
						//reco::GenParticleRef genMu1 = (muons->begin()+mu1Index)->genParticleRef();
						//reco::GenParticleRef genMu2 = (muons->begin()+mu2Index)->genParticleRef();
						//if (genMu1->motherRef()==genMu2->motherRef()) std::cout<<"genMu1->motherRef()->pdgId()="<<genMu1->motherRef()->pdgId()<<std::endl;
						//else std::cout<<"genMu1->motherRef()->pdgId()="<<genMu1->motherRef()->pdgId()<<", genMu2->motherRef()->pdgId()="<<genMu2->motherRef()->pdgId()
						// <<"genMu1->motherRef()->motherRef()->pdgId()="<<genMu1->motherRef()->motherRef()->pdgId()<<", genMu2->motherRef()->motherRef()->pdgId()="<<genMu2->motherRef()->motherRef()->pdgId()<<std::endl;

						//reco::GenParticleRef genMu3 = mu3->genParticleRef();
						//reco::GenParticleRef genMu4 = mu4->genParticleRef();
						if (genMu3.isNonnull() ){
							mu3_pdgID_bestYMass = genMu3->pdgId();
							/*
								if (genMu3->numberOfMothers()>0){ 
								reco::GenParticleRef mom3 = genMu3->motherRef();
								if (mom3.isNonnull()) { 
								std::cout<<""<<"mom pdgID= "<<mom3->pdgId()<<std::endl;
								if (mom3==genMu1->motherRef()) std::cout<<"same mother"<<std::endl;
								}    
								else std::cout<<"mom non"<<std::endl;
								}    
								else std::cout<<"# mom = 0"<<std::endl;
								*/
						}
						else mu3_pdgID_bestYMass = 0;
						if (genMu4.isNonnull() ) mu4_pdgID_bestYMass = genMu4->pdgId();
						else mu4_pdgID_bestYMass = 0;
					}
				}
			}
		}
	}
}



//define this as a plug-in
DEFINE_FWK_MODULE(MuMuGammaRootupler);
