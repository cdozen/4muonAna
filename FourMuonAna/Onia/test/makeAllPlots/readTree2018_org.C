#define readTree2018_cxx
#include "readTree2018.h"
#include <TH2.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <TGraph.h>
#include "TLorentzVector.h"
#include <vector>

void readTree2018::Loop()
{
	const bool blind_signal = true;
	TString plot_format = "png";
	gStyle->SetOptStat(kFALSE);

	if (fChain == 0) return;

	std::vector<TLorentzVector> mu1_p4_vector;
	std::vector<TLorentzVector> mu2_p4_vector;
	std::vector<TLorentzVector> mu3_p4_vector;
	std::vector<TLorentzVector> mu4_p4_vector;
	std::vector<TLorentzVector> mu12_p4_vector;
	std::vector<TLorentzVector> mu34_p4_vector;
	std::vector<TLorentzVector> mu1boost_p4_vector;
	std::vector<TLorentzVector> mu2boost_p4_vector;
	std::vector<TLorentzVector> mu3boost_p4_vector;
	std::vector<TLorentzVector> mu4boost_p4_vector;
	std::vector<TLorentzVector> mu12boost_p4_vector;
	std::vector<TLorentzVector> mu34boost_p4_vector;


	TH1F *vProbMix = new TH1F("vProbMix","vProbMix",100,0,1);
	vProbMix->GetXaxis()->SetTitle("vertex fit probability");
	TH1F *vProbSig = new TH1F("vProbSig","vProbSig",100,0,1);

	TH2F *Ymumu2D_zerobias_mix = new TH2F("Ymumu2D_zerobias_mix","Ymumu2D_zerobias_mix", 120,-6,6,120,-6,6);
	Ymumu2D_zerobias_mix->GetXaxis()->SetTitle("#Upsilon #eta");
	Ymumu2D_zerobias_mix->GetYaxis()->SetTitle("#mu_{3}#mu_{4} #eta");
	TH2F *Ymumu2DBoost = new TH2F("Ymumu2DBoost","Ymumu2DBoost",120,-6,6,120,-6,6);
	Ymumu2DBoost->GetXaxis()->SetTitle("#Upsilon #eta Boost");
	Ymumu2DBoost->GetYaxis()->SetTitle("#mu_{3}#mu_{4} #eta Boost");
	TH2F *Ymumu2D = new TH2F("Ymumu2D","Ymumu2D",120,-6,6,120,-6,6);
	Ymumu2D->GetXaxis()->SetTitle("#Upsilon #eta");
	Ymumu2D->GetYaxis()->SetTitle("#mu_{3}#mu_{4} #eta");
	TH2F *Ymumu2DBoost_aftercut = new TH2F("Ymumu2DBoost_aftercut","Ymumu2DBoost_aftercut",60,-6,6,60,-6,6);
	Ymumu2DBoost_aftercut->GetXaxis()->SetTitle("#Upsilon #eta Boost");
	Ymumu2DBoost_aftercut->GetYaxis()->SetTitle("#mu_{3}#mu_{4} #eta Boost");
	TH2F *Ymumu2D_aftercut = new TH2F("Ymumu2D_aftercut","Ymumu2D_aftercut",60,-6,6,60,-6,6);
	Ymumu2D_aftercut->GetXaxis()->SetTitle("#Upsilon #eta");
	Ymumu2D_aftercut->GetYaxis()->SetTitle("#mu_{3}#mu_{4} #eta");
	TH2F *Ymumu2D_physkbg_mix = new TH2F("Ymumu2D_physkbg_mix","Ymumu2D_physkbg_mix",60,-6,6,60,-6,6);
	Ymumu2D_physkbg_mix->GetXaxis()->SetTitle("#Upsilon #eta");
	Ymumu2D_physkbg_mix->GetYaxis()->SetTitle("#mu_{3}#mu_{4} #eta");
	TH2F *Ymumu2DBoost_physkbg_mix = new TH2F("Ymumu2DBoost_physkbg_mix","Ymumu2DBoost_physkbg_mix",60,-6,6,60,-6,6);
	Ymumu2DBoost_physkbg_mix->GetXaxis()->SetTitle("#Upsilon #eta Boost");
	Ymumu2DBoost_physkbg_mix->GetYaxis()->SetTitle("#mu_{3}#mu_{4} #eta Boost");

	TH1F *mu12mass = new TH1F("mu12mass","mu12mass",60,8.5,11.5);
	mu12mass->GetXaxis()->SetTitle("#mu_{1}#mu_{2} mass (GeV)");
	TH1F *mu12massEBE = new TH1F("mu12massEBE","mu12massEBE",60,8.5,11.5);
	mu12massEBE->GetXaxis()->SetTitle("#mu_{1}#mu_{2} mass (GeV)");
	TH1F *mu34mass = new TH1F("mu34mass","mu34mass",60,8,11);
	mu34mass->GetXaxis()->SetTitle("#mu_{3}#mu_{4} mass (GeV)");
	TH1F *mu34massBkg = new TH1F("mu34massBkg","mu34massBlg",60,8,11);
	mu34massBkg->GetXaxis()->SetTitle("#mu_{3}#mu_{4} mass (GeV)");
	TH1F *mu34massBkgH = new TH1F("mu34massBkgH","mu34massBlgH",60,8,11);
	mu34massBkgH->GetXaxis()->SetTitle("#mu_{3}#mu_{4} mass (GeV)");
	TH1F *Ypt1 = new TH1F("Ypt1","Ypt1",20,0,40);
	Ypt1->GetXaxis()->SetTitle("#Upsilon p_{T} (GeV)");
	TH1F *Ypt2 = new TH1F("Ypt2","Ypt2",20,0,40);
	Ypt2->GetXaxis()->SetTitle("#Upsilon p_{T} (GeV)");
	TH1F *Ypt3 = new TH1F("Ypt3","Ypt3",20,0,40);
	Ypt3->GetXaxis()->SetTitle("#Upsilon p_{T} (GeV)");
	TH1F *Ypt4 = new TH1F("Ypt4","Ypt4",20,0,40);
	Ypt4->GetXaxis()->SetTitle("#Upsilon p_{T} (GeV)");
	TH2F *Y2Dpteta = new TH2F("Y2Dpteta","Y2Dpteta",20,-5,5,20,0,40);
	Y2Dpteta->GetXaxis()->SetTitle("#Upsilon #eta");
	Y2Dpteta->GetYaxis()->SetTitle("#Upsilon p_{T} (GeV)");
	TH2F *Y2Dptphi = new TH2F("Y2Dptphi","Y2Dptphi",20,-4,4,20,0,40);
	Y2Dptphi->GetXaxis()->SetTitle("#Upsilon #phi");
	Y2Dptphi->GetYaxis()->SetTitle("#Upsilon p_{T} (GeV)");
	TH2F *Y2Detaphi = new TH2F("Y2Detaphi","Y2Detaphi",20,-5,5,20,-4,4);
	Y2Detaphi->GetXaxis()->SetTitle("#Upsilon #eta");
	Y2Detaphi->GetYaxis()->SetTitle("#Upsilon #phi");


	TH1F *hnCand = new TH1F("hnCand","hnCand",9,1,10);
	hnCand->GetXaxis()->SetTitle("number of 4 muon candidates per event");
	TH1F *hnCand_aftercut = new TH1F("hnCand_aftercut","hnCand_aftercut",9,1,10);
	hnCand_aftercut->GetXaxis()->SetTitle("number of 4 muon candidates per event");

	TH1F *shnCand = new TH1F("shnCand","shnCand",9,1,10);
	shnCand->GetXaxis()->SetTitle("number of 4 muon candidates per event");
	TH1F *shnCand_aftercut = new TH1F("shnCand_aftercut","shnCand_aftercut",9,1,10);
	shnCand_aftercut->GetXaxis()->SetTitle("number of 4 muon candidates per event");


	TH1F *hfourMuMass_mix = new TH1F("hfourMuMass_mix","hfourMuMass_mix",100,0,100);
	hfourMuMass_mix->GetXaxis()->SetTitle("4 muon mass [GeV]");
	TH1F *hfourMuMass_mix_aftercut = new TH1F("hfourMuMass_mix_aftercut","hfourMuMass_mix_aftercut",50,0,100);
	hfourMuMass_mix_aftercut->GetXaxis()->SetTitle("4 muon mass [GeV]");
	TH1F *hfourMuMass_mix_highLumi = new TH1F("hfourMuMass_mix_highLumi","hfourMuMass_mix_highLumi",50,0,100);
	hfourMuMass_mix_highLumi->GetXaxis()->SetTitle("4 muon mass [GeV]");
	TH1F *hfourMuMass_mix_lowLumi = new TH1F("hfourMuMass_mix_lowLumi","hfourMuMass_mix_lowLumi",100,0,100);
	hfourMuMass_mix_lowLumi->GetXaxis()->SetTitle("4 muon mass [GeV]");
	TH1F *hfourMuMass_mix_midLumi = new TH1F("hfourMuMass_mix_midLumi","hfourMuMass_mix_midLumi",100,0,100);
	hfourMuMass_mix_midLumi->GetXaxis()->SetTitle("4 muon mass [GeV]");
	TH1F *hfourMuMass_mix_smallrange = new TH1F("hfourMuMass_mix_smallrange","hfourMuMass_mix_smallrange",30,17,20);
	hfourMuMass_mix_smallrange->GetYaxis()->SetTitle("Candidates / 0.1 GeV");
	hfourMuMass_mix_smallrange->GetXaxis()->SetTitle("4 muon mass [GeV]");
	TH1F *hfourMuMass_mix_highLumi_smallrange = new TH1F("hfourMuMass_mix_highLumi_smallrange","hfourMuMass_mix_highLumi_smallrange",52,13,26);
	hfourMuMass_mix_highLumi_smallrange->GetXaxis()->SetTitle("4 muon mass [GeV]");
	TH1F *hfourMuMass_mix_lowLumi_smallrange = new TH1F("hfourMuMass_mix_lowLumi_smallrange","hfourMuMass_mix_lowLumi_smallrange",52,13,26);
	hfourMuMass_mix_lowLumi_smallrange->GetXaxis()->SetTitle("4 muon mass [GeV]");
	TH1F *hfourMuMass_mix_midLumi_smallrange = new TH1F("hfourMuMass_mix_midLumi_smallrange","hfourMuMass_mix_midLumi_smallrange",52,13,26);
	hfourMuMass_mix_midLumi_smallrange->GetXaxis()->SetTitle("4 muon mass [GeV]");


	TH1F *h_number_of_Y = new TH1F("h_number_of_Y","h_number_of_Y",10,0,10);
	TH1F *h_fourMu_pt_order = new TH1F("h_fourMu_pt_order","h_fourMu_pt_order",6,0,6);
	TH1F *hfourMuMass = new TH1F("hfourMuMass","hfourMuMass",100,0,100);
	hfourMuMass->GetXaxis()->SetTitle("4 muon mass [GeV]");
	hfourMuMass->GetYaxis()->SetTitle("Candidates / GeV");
	hfourMuMass->GetYaxis()->SetLabelSize(0.03);
	TH1F *hfourMuMass_aftercut = new TH1F("hfourMuMass_aftercut","hfourMuMass_aftercut",50,0,100);
	hfourMuMass_aftercut->GetXaxis()->SetTitle("4 muon mass [GeV]");
	hfourMuMass_aftercut->GetYaxis()->SetTitle("Candidates / 2 GeV");
	TH1F *hfourMuMass_aftercut_smallrange = new TH1F("hfourMuMass_aftercut_smallrange","hfourMuMass_aftercut_smallrange",65,13,26);
	hfourMuMass_aftercut_smallrange->GetXaxis()->SetTitle("4 muon mass [GeV]");
	hfourMuMass_aftercut_smallrange->GetYaxis()->SetTitle("Candidates / 0.2 GeV");

	TH1F *hfourMuMass_physkbg_mix = new TH1F("hfourMuMass_physkbg_mix","hfourMuMass_physkbg_mix",50,0,100);
	hfourMuMass_physkbg_mix->GetXaxis()->SetTitle("4 muon mass [GeV]");
	TH1F *hfourMuMass_physkbg_mix_smallrange = new TH1F("hfourMuMass_physkbg_mix_smallrange","hfourMuMass_physkbg_mix_smallrange",65,13,26);
	hfourMuMass_physkbg_mix_smallrange->GetXaxis()->SetTitle("4 muon mass [GeV]");
	hfourMuMass_physkbg_mix_smallrange->GetXaxis()->SetTitle("Candidates / 0.2 GeV");
	TH1F *hfourMuMassBoost_physkbg_mix = new TH1F("hfourMuMassBoost_physkbg_mix","hfourMuMassBoost_physkbg_mix",50,0,100);
	hfourMuMassBoost_physkbg_mix->GetXaxis()->SetTitle("4 muon mass [GeV]");
	TH1F *hfourMuMassBoost_physkbg_mix_smallrange = new TH1F("hfourMuMassBoost_physkbg_mix_smallrange","hfourMuMassBoost_physkbg_mix_smallrange",65,13,26);
	hfourMuMassBoost_physkbg_mix_smallrange->GetXaxis()->SetTitle("4 muon mass [GeV]");


	TH1F *hfourMuMass_aftercut_highLumi = new TH1F("hfourMuMass_aftercut_highLumi","hfourMuMass_aftercut_highLumi",50,0,100);
	hfourMuMass_aftercut_highLumi->GetXaxis()->SetTitle("4 muon mass [GeV]");
	hfourMuMass_aftercut_highLumi->GetYaxis()->SetTitle("Candidates / 2 GeV");
	TH1F *hfourMuMass_aftercut_lowLumi = new TH1F("hfourMuMass_aftercut_lowLumi","hfourMuMass_aftercut_lowLumi",50,0,100);
	hfourMuMass_aftercut_lowLumi->GetXaxis()->SetTitle("4 muon mass [GeV]");
	hfourMuMass_aftercut_lowLumi->GetYaxis()->SetTitle("Candidates / 2 GeV");
	TH1F *hfourMuMass_aftercut_highLumi_smallrange = new TH1F("hfourMuMass_aftercut_highLumi_smallrange","hfourMuMass_aftercut_highLumi_smallrange",26,13,26);
	hfourMuMass_aftercut_highLumi_smallrange->GetXaxis()->SetTitle("4 muon mass [GeV]");
	TH1F *hfourMuMass_aftercut_lowLumi_smallrange = new TH1F("hfourMuMass_aftercut_lowLumi_smallrange","hfourMuMass_aftercut_lowLumi_smallrange",26,13,26);
	hfourMuMass_aftercut_lowLumi_smallrange->GetXaxis()->SetTitle("4 muon mass [GeV]");


	TH1F *nVertices_2012 = new TH1F("nVertices_2012","nVertices_2012",50,0,100);
	nVertices_2012->GetXaxis()->SetTitle("number of primary vertices per event");
	TH1F *nVertices_fourmuon = new TH1F("nVertices_fourmuon","nVertices_fourmuon",50,0,100);
	nVertices_fourmuon->GetXaxis()->SetTitle("number of primary vertices per event");
	TH1F *nVertices = new TH1F("nVertices","nVertices",50,0,100);
	nVertices->GetXaxis()->SetTitle("number of primary vertices per event");
	TH1F *nVertices_used = new TH1F("nVertices_used","nVertices_used",50,0,100);
	nVertices_used->GetXaxis()->SetTitle("number of primary vertices per event");
	TH1F *nVertices_fourmuon_aftercut = new TH1F("nVertices_fourmuon_aftercut","nVertices_fourmuon_aftercut",50,0,100);
	nVertices_fourmuon_aftercut->GetXaxis()->SetTitle("number of primary vertices per event");

	TH1F *hmu1Tight = new TH1F("hmu1Tight","hmu1Tight",2,0,2);
	TH1F *hmu2Tight = new TH1F("hmu2Tight","hmu2Tight",2,0,2);
	TH1F *hmu3Tight = new TH1F("hmu3Tight","hmu3Tight",2,0,2);
	TH1F *hmu4Tight = new TH1F("hmu4Tight","hmu4Tight",2,0,2);
	TH1F *hmu1Medium = new TH1F("hmu1Medium","hmu1Medium",2,0,2);
	TH1F *hmu2Medium = new TH1F("hmu2Medium","hmu2Medium",2,0,2);
	TH1F *hmu3Medium = new TH1F("hmu3Medium","hmu3Medium",2,0,2);
	TH1F *hmu4Medium = new TH1F("hmu4Medium","hmu4Medium",2,0,2);
	TH1F *hfourMuFitmu1Pt = new TH1F("hfourMuFitmu1Pt","hfourMuFitmu1Pt",40,0,20);
	TH1F *hfourMuFitmu2Pt = new TH1F("hfourMuFitmu2Pt","hfourMuFitmu2Pt",40,0,20);
	TH1F *hfourMuFitmu3Pt = new TH1F("hfourMuFitmu3Pt","hfourMuFitmu3Pt",40,0,20);
	TH1F *hfourMuFitmu4Pt = new TH1F("hfourMuFitmu4Pt","hfourMuFitmu4Pt",40,0,20);

	TH1F *sfourMuMass = new TH1F("sfourMuMass","sfourMuMass",60,17,20);
	sfourMuMass->GetXaxis()->SetTitle("4 muon mass [GeV]");
	sfourMuMass->GetYaxis()->SetTitle("Candidates / 50 MeV"); 
	TH1F *smu1Tight = new TH1F("smu1Tight","smu1Tight",2,0,2);
	TH1F *smu2Tight = new TH1F("smu2Tight","smu2Tight",2,0,2);
	TH1F *smu3Tight = new TH1F("smu3Tight","smu3Tight",2,0,2);
	TH1F *smu4Tight = new TH1F("smu4Tight","smu4Tight",2,0,2);
	TH1F *smu1Medium = new TH1F("smu1Medium","smu1Medium",2,0,2);
	TH1F *smu2Medium = new TH1F("smu2Medium","smu2Medium",2,0,2);
	TH1F *smu3Medium = new TH1F("smu3Medium","smu3Medium",2,0,2);
	TH1F *smu4Medium = new TH1F("smu4Medium","smu4Medium",2,0,2);
	TH1F *sfourMuFitmu1Pt = new TH1F("sfourMuFitmu1Pt","sfourMuFitmu1Pt",40,0,20);
	TH1F *sfourMuFitmu2Pt = new TH1F("sfourMuFitmu2Pt","sfourMuFitmu2Pt",40,0,20);
	TH1F *sfourMuFitmu3Pt = new TH1F("sfourMuFitmu3Pt","sfourMuFitmu3Pt",40,0,20);
	TH1F *sfourMuFitmu4Pt = new TH1F("sfourMuFitmu4Pt","sfourMuFitmu4Pt",40,0,20);


	Long64_t nentries = fChain->GetEntries();
	std::cout<<nentries<<std::endl;
	Long64_t nbytes = 0, nb = 0;

	float muonPtCut[5]={2,2.5,3,3.5,4};
	float nPassMuonPtCut[5]={0};
	float sigPassMuonPtCut[5]={0};
	float sigEffMuonPtCut[5]={0};
	float bkRejMuonPtCut[5]={0};
	float significanceMuonPtCut[5]={0};

	float vProbCut[14]={0,0.0001,0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
	float nPassVProbCut[14]={0};
	float sigPassVProbCut[14]={0};
	float sigEffVProbCut[14]={0};
	float bkRejVProbCut[14]={0};
	float significanceVProbCut[14]={0};

	float nPass2d[5][14]={0};
	float sigPass2d[5][14]={0};
	float sigEff2d[5][14]={0};
	float bkRej2d[5][14]={0};
	float significance2d[5][14]={0};

	float quaCut[6]={0,1,2,3,4,5};
	float nPassQuaCut[6]={0};
	float sigPassQuaCut[6]={0};
	float sigEffQuaCut[6]={0};
	float bkRejQuaCut[6]={0};
	float significanceQuaCut[6]={0};

	float quaCut_m[6]={0,1,2,3,4,5};
	float nPassQuaCut_m[6]={0};
	float sigPassQuaCut_m[6]={0};
	float sigEffQuaCut_m[6]={0};
	float bkRejQuaCut_m[6]={0};
	float significanceQuaCut_m[6]={0};

	std::ofstream EventNumber;
	EventNumber.open("EventNumber.txt");
	std::ofstream fourMumassoutput;
	fourMumassoutput.open("fourMumassoutput.txt");
	std::ofstream mu34massoutput;
	mu34massoutput.open("mu34massoutput.txt");
	std::ofstream mu34massBkgoutput;
	mu34massBkgoutput.open("mu34massBkgoutput.txt");
	std::ofstream mu34massBkgHoutput;
	mu34massBkgHoutput.open("mu34massBkgHoutput.txt");
	std::ofstream mu12Momoutput;
	mu12Momoutput.open("mu12Momoutput.txt");

	//for (Long64_t jentry=0; jentry<290310;jentry++) {
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		//Make plots
		nVertices->Fill(numPrimaryVertices);

		//if ((trigger&2)==0) continue; //2016 trigger
		//if ((trigger&128)==0) continue;
		//if ((trigger&256)==0) continue;			//HLT_Trimuon5_3p5_2_Upsilon_Muon
		//if ((trigger&512)==0) continue;		//HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon
		//if((trigger&768)==0) continue; 			// for both 256 and 512. triggers
		//std::cout << "v_mumufit_Mass->size() = " << v_mumufit_Mass->size() << std::endl;
		//std::cout << "fourMuFit_Mass->size() = " << fourMuFit_Mass->size() << std::endl;
		h_number_of_Y->Fill(v_mumufit_Mass->size());

		//orignial events in the MuOnia dataset
		int nCand = 0;
		int nCand_aftercut = 0;
		for (unsigned i=0;i<fourMuFit_Mass->size();++i) {
			//if (fabs(mumufit_Mass-9.4603) > (3*1.105*mumufit_MassErr)) continue;
			//if (fourMuFit_mu1Pt->at(i)<5 || fourMuFit_mu2Pt->at(i)<5 || fourMuFit_mu3Pt->at(i)<5 || fourMuFit_mu4Pt->at(i)<5) continue;
			TLorentzVector fourMuFit,mu1,mu2,mu3,mu4,mu12,mu34,mu1boost,mu2boost,mu3boost,mu4boost,mu12boost,mu34boost;

			fourMuFit.SetPtEtaPhiM(fourMuFit_Pt->at(i),fourMuFit_Eta->at(i),fourMuFit_Phi->at(i),fourMuFit_Mass->at(i));
			//std::cout<<fourMuFit_Mass->at(i)<<" "<<fourMuFit.M()<<" "<<fourMuFit_Pt->at(i)<<" "<<fourMuFit.Pt()<<" "<<fourMuFit.Pz()<<std::endl;
			mu1.SetPtEtaPhiE(fourMuFit_mu1Pt->at(i),fourMuFit_mu1Eta->at(i),fourMuFit_mu1Phi->at(i),fourMuFit_mu1E->at(i));
			//std::cout<<mu1.Pt()<<" "<<mu1.M()<<std::endl;
			mu2.SetPtEtaPhiE(fourMuFit_mu2Pt->at(i),fourMuFit_mu2Eta->at(i),fourMuFit_mu2Phi->at(i),fourMuFit_mu2E->at(i));
			mu3.SetPtEtaPhiE(fourMuFit_mu3Pt->at(i),fourMuFit_mu3Eta->at(i),fourMuFit_mu3Phi->at(i),fourMuFit_mu3E->at(i));
			mu4.SetPtEtaPhiE(fourMuFit_mu4Pt->at(i),fourMuFit_mu4Eta->at(i),fourMuFit_mu4Phi->at(i),fourMuFit_mu4E->at(i));
			mu12=mu1+mu2;
			mu34=mu3+mu4;
			EventNumber << event << "\n";
			hfourMuMass->Fill(fourMuFit_Mass->at(i));
			Ymumu2D->Fill(mu12.Eta(),mu34.Eta());

			//Boost to center of mass frame
			mu1boost.SetPtEtaPhiE(fourMuFit_mu1Pt->at(i),fourMuFit_mu1Eta->at(i),fourMuFit_mu1Phi->at(i),fourMuFit_mu1E->at(i));
			mu2boost.SetPtEtaPhiE(fourMuFit_mu2Pt->at(i),fourMuFit_mu2Eta->at(i),fourMuFit_mu2Phi->at(i),fourMuFit_mu2E->at(i));
			mu3boost.SetPtEtaPhiE(fourMuFit_mu3Pt->at(i),fourMuFit_mu3Eta->at(i),fourMuFit_mu3Phi->at(i),fourMuFit_mu3E->at(i));
			mu4boost.SetPtEtaPhiE(fourMuFit_mu4Pt->at(i),fourMuFit_mu4Eta->at(i),fourMuFit_mu4Phi->at(i),fourMuFit_mu4E->at(i));
			mu1boost.Boost(-fourMuFit.BoostVector());
			//std::cout<<mu1boost.Pt()<<" "<<mu1boost.M()<<std::endl;
			mu2boost.Boost(-fourMuFit.BoostVector());
			mu3boost.Boost(-fourMuFit.BoostVector());
			mu4boost.Boost(-fourMuFit.BoostVector());
			mu12boost=mu1boost+mu2boost;
			mu34boost=mu3boost+mu4boost;
			Ymumu2DBoost->Fill(mu12boost.Eta(),mu34boost.Eta());

			nCand++;

			//std::array<float, 4> s = {fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i), fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i)};
			//std::sort(s.begin(), s.end());		//sort from small to big
			//bool pass_trigger_pt = s.at(1) > 2 && s.at(2) > 3.5 && s.at(3) > 5;
			//if(pass_trigger_pt == false) std::cout << "mu1Pt = " << fourMuFit_mu1Pt->at(i)  << "mu2Pt = " << fourMuFit_mu2Pt->at(i)  << "mu3Pt = " << fourMuFit_mu3Pt->at(i)  << "mu4Pt = " << fourMuFit_mu4Pt->at(i) << std::endl;

			if (fourMuFit_Mass->at(i)>13 && fourMuFit_Mass->at(i)<23) nVertices_fourmuon->Fill(numPrimaryVertices);

			if (
					fourMuFit_VtxProb->at(i)>0.05
					&& (mu3_Medium->at(i) + mu4_Medium->at(i)) >= 2
					&& mu34.M()< 9.2
					&& fourMuFit_mu1Pt->at(i) >= muonPtCut[2] && fourMuFit_mu2Pt->at(i) >= muonPtCut[2] && fourMuFit_mu3Pt->at(i) >= muonPtCut[2] && fourMuFit_mu4Pt->at(i) >= muonPtCut[2]
					//&& pass_trigger_pt
					//&& fourMuFit_Mass->at(i)>13 && fourMuFit_Mass->at(i)<26
				) {
				mu12mass->Fill(mu12.M());
				/*if (fabs(mumufit_Mass-9.4603) > (3*1.105*mumufit_MassErr)) {
				  if (mumufit_Mass < 9.1 and mumufit_Mass > 8.7) {
				  mu34massBkg->Fill(mu34.M());
				  mu34massBkgoutput << mu34.M() << "\n"; 
				  }
				  if (mumufit_Mass < 10.2 and mumufit_Mass > 9.8) {
				  mu34massBkgH->Fill(mu34.M());
				  mu34massBkgHoutput << mu34.M() << "\n"; 
				  }
				  }
				  if (fabs(mumufit_Mass-9.4603) > (3*1.105*mumufit_MassErr)) continue;
				  mu12massEBE->Fill(mu12.M()); 
				  if (mumufit_Mass < (9.66) and mumufit_Mass > 9.26) mu34mass->Fill(mu34.M());
				  */
				mu34massoutput << mu34.M() << " "<< fourMuFit_mu1Pt->at(i) << " "<< fourMuFit_mu2Pt->at(i) << " "<< fourMuFit_mu3Pt->at(i) << " "<< fourMuFit_mu4Pt->at(i) << "\n";
				mu12Momoutput << mu12.Pt() << " "<< mu12.Eta()  << " " << mu12.Phi() << " " << mu12.M() << "\n";
				hfourMuMass_aftercut->Fill(fourMuFit_Mass->at(i));
				hfourMuMass_aftercut_smallrange->Fill(fourMuFit_Mass->at(i));
				fourMumassoutput << fourMuFit_Mass->at(i)<< "\n";
				Ymumu2D_aftercut->Fill(mu12.Eta(),mu34.Eta());
				Ymumu2DBoost_aftercut->Fill(mu12boost.Eta(),mu34boost.Eta());
				if (fabs(mu12.Eta())<0.2) Ypt1->Fill(mu12.Pt());
				if (fabs(mu12.Eta())<1.1 && fabs(mu12.Eta())>0.9) Ypt2->Fill(mu12.Pt());
				if (fabs(mu12.Eta())<2.1 && fabs(mu12.Eta())>1.9) Ypt3->Fill(mu12.Pt());
				if (fabs(mu12.Eta())<3.1 && fabs(mu12.Eta())>2.9) Ypt4->Fill(mu12.Pt());
				Y2Dpteta->Fill(mu12.Eta(),mu12.Pt());
				Y2Dptphi->Fill(mu12.Phi(),mu12.Pt()); 
				Y2Detaphi->Fill(mu12.Eta(),mu12.Phi()); 	

				bool fourMu_pt_fill = false;

				if (std::min(fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i)) > std::max(fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i)))
				{h_fourMu_pt_order->Fill(0); fourMu_pt_fill = true;}		//yymm
				if (std::max(fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i)) > std::max(fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i))
				&& std::max(fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i)) > std::min(fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i))
				&& std::min(fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i)) > std::min(fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i)))
				{h_fourMu_pt_order->Fill(1); fourMu_pt_fill = true;}		//ymym
				if (std::max(fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i)) > std::max(fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i))
				&& std::min(fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i)) > std::min(fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i)))
				{h_fourMu_pt_order->Fill(2); fourMu_pt_fill = true;}		//ymmy
				if (std::max(fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i)) > std::max(fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i))
				&& std::min(fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i)) > std::min(fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i)))
				{h_fourMu_pt_order->Fill(3); fourMu_pt_fill = true;}		//myym
				if (std::max(fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i)) > std::max(fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i))
				&& std::max(fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i)) > std::min(fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i))
				&& std::min(fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i)) > std::min(fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i)))
				{h_fourMu_pt_order->Fill(4); fourMu_pt_fill = true;}		//mymy
				if (std::min(fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i)) > std::max(fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i)))
				{h_fourMu_pt_order->Fill(5); fourMu_pt_fill = true;}		//mmyy

				if(fourMu_pt_fill == false) std::cout << "mu1Pt = " << fourMuFit_mu1Pt->at(i)  << "mu2Pt = " << fourMuFit_mu2Pt->at(i)  << "mu3Pt = " << fourMuFit_mu3Pt->at(i)  << "mu4Pt = " << fourMuFit_mu4Pt->at(i) << std::endl;

				nCand_aftercut++;

				if (numPrimaryVertices>30) {
					hfourMuMass_aftercut_highLumi->Fill(fourMuFit_Mass->at(i));
					hfourMuMass_aftercut_highLumi_smallrange->Fill(fourMuFit_Mass->at(i));
				}
				else {
					hfourMuMass_aftercut_lowLumi->Fill(fourMuFit_Mass->at(i));
					hfourMuMass_aftercut_lowLumi_smallrange->Fill(fourMuFit_Mass->at(i));
				}

				if (fourMuFit_Mass->at(i)>13 && fourMuFit_Mass->at(i)<23) nVertices_fourmuon_aftercut->Fill(numPrimaryVertices);


				//fill vectors for mixing later
				mu1_p4_vector.push_back(mu1);
				mu2_p4_vector.push_back(mu2);
				mu3_p4_vector.push_back(mu3);
				mu4_p4_vector.push_back(mu4);
				mu12_p4_vector.push_back(mu12);
				mu34_p4_vector.push_back(mu34);

				mu1boost_p4_vector.push_back(mu1boost);
				mu2boost_p4_vector.push_back(mu2boost);
				mu3boost_p4_vector.push_back(mu3boost);
				mu4boost_p4_vector.push_back(mu4boost);
				mu12boost_p4_vector.push_back(mu12boost);
				mu34boost_p4_vector.push_back(mu34boost);


			}
		}
		hnCand->Fill(nCand);
		hnCand_aftercut->Fill(nCand_aftercut);

	}


	//physics background mixing
	std::ofstream fourmuonMassMixBkg;
	fourmuonMassMixBkg.open("fourmuonMassMixBkg2.txt");
	//for (unsigned j=0; j<1000;j++){ 
	//	for (unsigned k=0; k<1000;k++){
	for (unsigned j=0; j<mu34_p4_vector.size();j++){ 
		for (unsigned k=0; k<mu12_p4_vector.size();k++){
			//for (unsigned k=j; k<j+1000 && k<mu12_p4_vector.size();k++){
			if (j==k) continue;
			//if ( fabs(mu12boost_p4_vector[k].Vect().DeltaR( mu12boost_p4_vector[j].Vect())) > 1.5) continue;
			if ( fabs(mu12_p4_vector[k].Vect().DeltaR( mu12_p4_vector[j].Vect())) > 0.3) continue;
			//if ( fabs(mu12_p4_vector[k].Vect().DeltaR( mu12_p4_vector[j].Vect())) > 0.05) continue;
			//if ( fabs(mu12_p4_vector[k].Phi() - mu12_p4_vector[j].Phi()) < 0.1) continue;
			//if ( fabs(mu12_p4_vector[k].Eta() - mu12_p4_vector[j].Eta()) > 0.05 ) continue;
			//if ( fabs(mu12_p4_vector[k].Phi() - mu12_p4_vector[j].Phi()) > 0.5 ) continue;
			//mu12_p4_vector[k].SetPtEtaPhiM(mu12_p4_vector[k].Pt(),mu12_p4_vector[k].Eta(),mu12_p4_vector[j].Phi(),mu12_p4_vector[k].M());
			//if ( fabs(mu12_p4_vector[k].Vect().DeltaR( mu12_p4_vector[j].Vect())) > 0.5 || fabs(mu12_p4_vector[k].Vect().DeltaR( mu12_p4_vector[j].Vect())) < 0.3) continue;
			//if ( mu12_p4_vector[k].Eta() * mu12_p4_vector[j].Eta() < 0 ) continue;
			//if (fabs( mu34_p4_vector[k].M()-mu34_p4_vector[j].M() )>0.5*mu34_p4_vector[k] ) continue;
			//if ( mu34.M() > 9.) continue;

			//std::array<double, 4> s = {mu1_p4_vector[k].Pt(), mu2_p4_vector[k].Pt(), mu3_p4_vector[j].Pt(), mu4_p4_vector[j].Pt()};
			//std::sort(s.begin(), s.end());		//sort from small to big
			//bool pass_trigger_pt = s.at(1) > 2 && s.at(2) > 3.5 && s.at(3) > 5;
			//if(pass_trigger_pt == false) continue;

			TLorentzVector mixFourMu;
			mixFourMu=mu12_p4_vector[k]+mu34_p4_vector[j];
			//std::cout<<mixFourMu.Pt()<<" "<<mixFourMu.Pz()<<" "<<mixFourMu.M()<<std::endl;
			//mixFourMu.Boost(fourMuFit.BoostVector());
			//std::cout<<mixFourMu.Pt()<<" "<<mixFourMu.Pz()<<" "<<mixFourMu.M()<<std::endl;
			hfourMuMass_physkbg_mix->Fill(mixFourMu.M());
			hfourMuMass_physkbg_mix_smallrange->Fill(mixFourMu.M());
			Ymumu2D_physkbg_mix->Fill(mu12_p4_vector[k].Eta(),mu34_p4_vector[j].Eta());
			if (mixFourMu.M()>13 && mixFourMu.M()<26) fourmuonMassMixBkg<<mixFourMu.M()<<"\n";

			TLorentzVector mixFourMuBoost;
			mixFourMuBoost=mu12boost_p4_vector[k]+mu34boost_p4_vector[j];
			//std::cout<<mixFourMu.Pt()<<" "<<mixFourMu.Pz()<<" "<<mixFourMu.M()<<std::endl;
			//mixFourMu.Boost(fourMuFit.BoostVector());
			//std::cout<<mixFourMu.Pt()<<" "<<mixFourMu.Pz()<<" "<<mixFourMu.M()<<std::endl;
			hfourMuMassBoost_physkbg_mix->Fill(mixFourMuBoost.M());
			hfourMuMassBoost_physkbg_mix_smallrange->Fill(mixFourMuBoost.M());
			Ymumu2DBoost_physkbg_mix->Fill(mu12boost_p4_vector[k].Eta(),mu34boost_p4_vector[j].Eta());

			/*if (mixFourMu.M()>18 && mixFourMu.M()<19) {
			  for (int p=0; p<5; p++) {
			  if (mu1_p4_vector[k].Pt() >= muonPtCut[p] && mu2_p4_vector[k].Pt() >= muonPtCut[p] && mu3_p4_vector[j].Pt() >= muonPtCut[p] && mu4_p4_vector[j].Pt() >= muonPtCut[p]) {
			  nPassMuonPtCut[p]++;
			  }
			  }
			  }*/
		}   
		}   
		fourmuonMassMixBkg.close();

		TCanvas *c1 = new TCanvas("c1","c1");
		sfourMuMass->Draw("e1");
		//c1->SaveAs("plots0p6/signalFourMuMass."+plot_format);


		//before cuts, mix and mix zoom in
		hfourMuMass_mix->SetMarkerStyle(24);
		hfourMuMass_mix->SetMarkerColor(kRed);
		hfourMuMass_mix->Draw("e1");
		//c1->SaveAs("plots0p6/hfourMuMass_mix."+plot_format);
		hfourMuMass_mix_smallrange->SetMarkerStyle(24);
		hfourMuMass_mix_smallrange->SetMarkerColor(kBlue);
		hfourMuMass_mix_smallrange->GetYaxis()->SetRangeUser(0,200);
		hfourMuMass_mix_smallrange->Draw("e1");
		//c1->SaveAs("plots0p6/hfourMuMass_mix_smallrange."+plot_format);

		//before cuts, compare original vs mix
		hfourMuMass->Draw("e1");
		//c1->SaveAs("plots0p6/hfourMuMass."+plot_format);
		hfourMuMass_mix->Draw("e1same");
		//c1->SaveAs("plots0p6/hfourMuMass_origVSmix."+plot_format);

		//after cuts, orignal
		hfourMuMass_aftercut->Draw("e1");
		//c1->SaveAs("plots0p6/hfourMuMass_afterCut."+plot_format);
		float scaleEntries = hfourMuMass_aftercut_smallrange->Integral();
		/*hfourMuMass_aftercut_smallrange->SetBinContent(23,0);
		  hfourMuMass_aftercut_smallrange->SetBinContent(24,0);
		  hfourMuMass_aftercut_smallrange->SetBinContent(25,0);
		  hfourMuMass_aftercut_smallrange->SetBinContent(26,0);
		  hfourMuMass_aftercut_smallrange->SetBinContent(27,0);
		  hfourMuMass_aftercut_smallrange->SetBinContent(28,0);
		  hfourMuMass_aftercut_smallrange->SetBinContent(29,0);
		  hfourMuMass_aftercut_smallrange->SetBinContent(30,0);
		  hfourMuMass_aftercut_smallrange->SetBinContent(31,0);
		  hfourMuMass_aftercut_smallrange->SetBinContent(32,0);*/
		hfourMuMass_aftercut_smallrange->Draw("e1");
		//c1->SaveAs("plots0p6/hfourMuMass_afterCut_smallrange."+plot_format);

		//after cuts, mix
		hfourMuMass_mix_aftercut->SetMarkerStyle(24);
		hfourMuMass_mix_aftercut->SetMarkerColor(kRed);
		hfourMuMass_mix_aftercut->Draw("e1");
		//c1->SaveAs("plots0p6/hfourMuMass_mix_afterCut."+plot_format);

		//after cuts, compare original vs mix
		hfourMuMass_aftercut->Draw("e1");
		hfourMuMass_mix_aftercut->Draw("e1same");
		//c1->SaveAs("plots0p6/hfourMuMass_afterCut_origVSmix."+plot_format);

		//after cuts, mix physics bkg
		hfourMuMass_physkbg_mix->SetMarkerStyle(22);
		hfourMuMass_physkbg_mix->SetMarkerColor(kBlue);
		hfourMuMass_physkbg_mix->Draw("e1");
		//c1->SaveAs("plots0p6/hfourMuMass_physkbg_mix."+plot_format);
		hfourMuMass_physkbg_mix_smallrange->SetMarkerStyle(22);
		hfourMuMass_physkbg_mix_smallrange->SetMarkerColor(kBlue);
		hfourMuMass_physkbg_mix_smallrange->Draw("e1");
		//c1->SaveAs("plots0p6/hfourMuMass_physkbg_mix_smallrange."+plot_format);

		//after cuts, compare original vs mixPhysPkg
		hfourMuMass_aftercut->SetMarkerStyle(20);
		hfourMuMass_aftercut->SetMarkerColor(kBlack);
		hfourMuMass_aftercut->SetLineColor(kBlack);
		hfourMuMass_aftercut->Draw("e1");
		std::cout << "hnCand_aftercut->Integral() = " << hnCand_aftercut->Integral() << std::endl;
		std::cout << "h_fourMu_pt_order->Integral() = " << h_fourMu_pt_order->Integral() << std::endl;
		std::cout << "hfourMuMass_aftercut->Integral() = " << hfourMuMass_aftercut->Integral() << std::endl;
		hfourMuMass_physkbg_mix->Scale(hfourMuMass_aftercut->Integral()/hfourMuMass_physkbg_mix->Integral());
		hfourMuMass_physkbg_mix->Draw("e1same");
		std::cout << "hfourMuMass_physkbg_mix->Integral() = " << hfourMuMass_physkbg_mix->Integral() << std::endl;
		TString total_data = "data (" + std::to_string((int)hfourMuMass_aftercut->Integral()) + " total)";
		TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
		leg->AddEntry(hfourMuMass_aftercut,total_data,"lep");
		leg->AddEntry(hfourMuMass_physkbg_mix,"mixed BG","lep");
		leg->Draw("same");
		c1->SaveAs("plots0p6/hfourMuMass_afterCut_origVSphyskbg."+plot_format);

		hfourMuMass_aftercut_smallrange->SetMarkerStyle(20);
		hfourMuMass_aftercut_smallrange->SetMarkerColor(kBlack);
		//hfourMuMass_aftercut_smallrange->Draw("e1");
		TH1F *hfourMuMass_aftercut_smallrange_copy = (TH1F*)hfourMuMass_aftercut_smallrange->Clone("hfourMuMass_aftercut_smallrange_copy");

		if(blind_signal)
		{
			for(int i=1; i <= hfourMuMass_aftercut_smallrange_copy->GetSize()-2; i++)
			{
				if(i>= hfourMuMass_aftercut_smallrange_copy->FindBin(17.5) && i<= hfourMuMass_aftercut_smallrange_copy->FindBin(19.5))
				hfourMuMass_aftercut_smallrange_copy->SetBinContent(i,0);
			}
		}
		hfourMuMass_aftercut_smallrange_copy->Draw("e1");

		hfourMuMass_physkbg_mix_smallrange->Scale(scaleEntries/hfourMuMass_physkbg_mix_smallrange->Integral());
		//hfourMuMass_physkbg_mix_smallrange->Scale(scaleEntries*5/2/hfourMuMass_physkbg_mix_smallrange->Integral());  //scale for different binning
		hfourMuMass_physkbg_mix_smallrange->Draw("e1same");
		leg->Draw("same");
		c1->SaveAs("plots0p6/hfourMuMass_afterCut_origVSphyskbg_smallrange."+plot_format);

		TFile* out_file = new TFile("plots0p6/fourMuMass.root","RECREATE");
		hfourMuMass_aftercut->Write();
		hfourMuMass_physkbg_mix->Write();
		hfourMuMass_aftercut_smallrange->Write();
		hfourMuMass_aftercut_smallrange_copy->Write();
		hfourMuMass_physkbg_mix_smallrange->Write();
		out_file->Close();

		//after cuts, compare original vs mixPhysPkg
		hfourMuMass_aftercut->Draw("e1");
		hfourMuMassBoost_physkbg_mix->SetMarkerStyle(23);
		hfourMuMassBoost_physkbg_mix->SetMarkerColor(kRed);
		hfourMuMassBoost_physkbg_mix->Scale(hfourMuMass_aftercut->Integral()/hfourMuMassBoost_physkbg_mix->Integral());
		hfourMuMassBoost_physkbg_mix->Draw("e1");
		hfourMuMass_aftercut->Draw("e1same");
		//c1->SaveAs("plots0p6/hfourMuMassBoost_afterCut_origVSphyskbg."+plot_format);
		hfourMuMass_aftercut_smallrange->Draw("e1");
		hfourMuMassBoost_physkbg_mix_smallrange->SetMarkerStyle(23);
		hfourMuMassBoost_physkbg_mix_smallrange->SetMarkerColor(kRed);
		hfourMuMassBoost_physkbg_mix_smallrange->Scale(scaleEntries*5/2/hfourMuMassBoost_physkbg_mix_smallrange->Integral());
		hfourMuMassBoost_physkbg_mix_smallrange->Draw("e1same");
		//c1->SaveAs("plots0p6/hfourMuMassBoost_afterCut_origVSphyskbg_smallrange."+plot_format);

		TCanvas *c9 = new TCanvas("c9","c9");
		Ymumu2D->Draw("colz");
		//c9->SaveAs("plots0p6/Ymumu2D."+plot_format);
		Ymumu2DBoost->Draw("colz");
		//c9->SaveAs("plots0p6/Ymumu2DBoost."+plot_format);
		Ymumu2D_aftercut->Draw("colz");
		//c9->SaveAs("plots0p6/Ymumu2D_aftercut."+plot_format);
		Ymumu2DBoost_aftercut->Draw("colz");
		//c9->SaveAs("plots0p6/Ymumu2DBoost_aftercut."+plot_format);
		Ymumu2D_zerobias_mix->Draw("colz");
		//c9->SaveAs("plots0p6/Ymumu2D_zerobias_mix."+plot_format);
		Ymumu2D_physkbg_mix->Draw("colz");
		//c9->SaveAs("plots0p6/Ymumu2D_physkbg_mix."+plot_format);
		Ymumu2DBoost_physkbg_mix->Draw("colz");
		//c9->SaveAs("plots0p6/Ymumu2DBoost_physkbg_mix."+plot_format);
		h_number_of_Y->Draw("e1");
		c9->SetLogy();
		c9->SaveAs("plots0p6/number_of_Y."+plot_format);

		TCanvas *c10 = new TCanvas("c10","c10");
		//Ypt1->Scale(1./Ypt1->Integral());
		Ypt1->Draw("e1");
		//c10->SaveAs("plots0p6/Ypt1."+plot_format);
		//Ypt2->Scale(Ypt1->Integral()/Ypt2->Integral());
		//Ypt2->SetMarkerColor(kRed);
		Ypt2->Draw("e1");
		//c10->SaveAs("plots0p6/Ypt2."+plot_format);
		//Ypt3->Scale(Ypt1->Integral()/Ypt3->Integral());
		//Ypt3->SetMarkerColor(kBlue);
		Ypt3->Draw("e1");
		//c10->SaveAs("plots0p6/Ypt3."+plot_format);
		//Ypt4->Scale(Ypt1->Integral()/Ypt4->Integral());
		//Ypt4->SetMarkerColor(kMagenta);
		Ypt4->Draw("e1");
		//c10->SaveAs("plots0p6/Ypt4."+plot_format);
		Y2Dpteta->Draw("colz");
		//c10->SaveAs("plots0p6/Y2Dpteta."+plot_format);
		Y2Dptphi->Draw("colz");
		//c10->SaveAs("plots0p6/Y2Dptphi."+plot_format);
		Y2Detaphi->Draw("colz");
		//c10->SaveAs("plots0p6/Y2Detaphi."+plot_format);
		mu12mass->Draw();
		mu12massEBE->SetLineColor(kRed);
		mu12massEBE->Draw("same");
		//c10->SaveAs("plots0p6/Ymass."+plot_format);
		mu34massBkg->Draw("e1");
		mu34massBkg->GetYaxis()->SetRangeUser(0,50);
		//c10->SaveAs("plots0p6/mu34massBkg."+plot_format);
		mu34massBkgH->Draw("e1");
		mu34massBkgH->GetYaxis()->SetRangeUser(0,50);
		//c10->SaveAs("plots0p6/mu34massBkgH."+plot_format);
		mu34mass->Draw("e1");
		mu34mass->GetYaxis()->SetRangeUser(0,50);
		//c10->SaveAs("plots0p6/mu34mass."+plot_format);
		mu34massBkg->SetMarkerColor(kRed);
		mu34massBkg->Draw("e1same");
		mu34massBkgH->SetMarkerColor(kBlue);
		mu34massBkgH->Draw("e1same");
		//c10->SaveAs("plots0p6/mu34massOverlay."+plot_format);
		h_fourMu_pt_order->Draw("e1");
		h_fourMu_pt_order->SetMinimum(0);
		c10->SaveAs("plots0p6/ourMu_pt_order."+plot_format);
		

		TCanvas *c8 = new TCanvas("c8","c8");
		hnCand_aftercut->SetMarkerColor(kRed);
		hnCand_aftercut->Draw("e1");
		//hnCand->Scale((float)hnCand_aftercut->Integral()/hnCand->Integral());
		//hnCand->SetMarkerSize(0.7);
		//hnCand->Draw("same");
		c8->SaveAs("plots0p6/nCand."+plot_format);
		shnCand_aftercut->SetMarkerColor(kRed);
		shnCand_aftercut->Draw("e1");
		shnCand->Scale((float)shnCand_aftercut->Integral()/shnCand->Integral());
		shnCand->SetMarkerSize(0.7);
		shnCand->Draw("same");
		//c8->SaveAs("plots0p6/nCand_signal."+plot_format);
		vProbMix->SetMarkerColor(kBlue);
		vProbMix->Draw("e1");
		c8->SetLogy();
		vProbSig->Scale((float)vProbMix->Integral()/vProbSig->Integral());
		vProbSig->SetMarkerColor(kRed);
		vProbSig->Draw("same");
		//c8->SaveAs("plots0p6/vProb_SandB."+plot_format);


		TCanvas *c7 = new TCanvas("c7","c7");	
		hfourMuMass_aftercut_lowLumi->SetMarkerColor(kRed);
		hfourMuMass_aftercut_lowLumi->Draw("e1");
		hfourMuMass_aftercut_highLumi->Scale((float)hfourMuMass_aftercut_lowLumi->Integral()/hfourMuMass_aftercut_highLumi->Integral());
		hfourMuMass_aftercut_highLumi->SetMarkerColor(kBlue);
		hfourMuMass_aftercut_highLumi->Draw("e1same");
		//c7->SaveAs("plots0p6/compareMassvsLumi."+plot_format);

		hfourMuMass_mix_lowLumi->SetMarkerColor(kRed);
		hfourMuMass_mix_lowLumi->Draw("e1");
		hfourMuMass_mix_midLumi->Scale((float)hfourMuMass_mix_lowLumi->Integral()/hfourMuMass_mix_midLumi->Integral());
		hfourMuMass_mix_midLumi->Draw("e1same");
		hfourMuMass_mix_highLumi->Scale((float)hfourMuMass_mix_lowLumi->Integral()/hfourMuMass_mix_highLumi->Integral());
		hfourMuMass_mix_highLumi->SetMarkerColor(kBlue);
		hfourMuMass_mix_highLumi->Draw("e1same");
		//c7->SaveAs("plots0p6/compareMassvsLumi_pileupbkg."+plot_format);

		hfourMuMass_aftercut_lowLumi_smallrange->SetMarkerColor(kRed);
		hfourMuMass_aftercut_lowLumi_smallrange->Draw("e1");
		hfourMuMass_aftercut_highLumi_smallrange->Scale((float)hfourMuMass_aftercut_lowLumi_smallrange->Integral()/hfourMuMass_aftercut_highLumi_smallrange->Integral());
		hfourMuMass_aftercut_highLumi_smallrange->SetMarkerColor(kBlue);
		hfourMuMass_aftercut_highLumi_smallrange->Draw("e1same");
		//c7->SaveAs("plots0p6/compareMassvsLumi_smallrange."+plot_format);

		hfourMuMass_mix_lowLumi_smallrange->SetMarkerColor(kRed);
		hfourMuMass_mix_lowLumi_smallrange->Draw("e1");
		hfourMuMass_mix_midLumi_smallrange->Scale((float)hfourMuMass_mix_lowLumi_smallrange->Integral()/hfourMuMass_mix_midLumi_smallrange->Integral());
		hfourMuMass_mix_midLumi_smallrange->Draw("e1same");	
		hfourMuMass_mix_highLumi_smallrange->Scale((float)hfourMuMass_mix_lowLumi_smallrange->Integral()/hfourMuMass_mix_highLumi_smallrange->Integral());
		hfourMuMass_mix_highLumi_smallrange->SetMarkerColor(kBlue);
		hfourMuMass_mix_highLumi_smallrange->Draw("e1same");
		//c7->SaveAs("plots0p6/compareMassvsLumi_pileupbkg_smallrange."+plot_format);



		TCanvas *c5 = new TCanvas("c5","c5");
		smu1Medium->GetXaxis()->SetTitle("Medium muon");
		smu1Medium->Draw("e1");
		smu2Medium->SetMarkerColor(kRed);
		smu2Medium->Draw("e1same");
		smu3Medium->SetMarkerColor(kBlue);
		smu3Medium->SetMarkerStyle(24);
		smu3Medium->Draw("e1same");
		smu4Medium->SetMarkerColor(kMagenta);
		smu4Medium->SetMarkerStyle(24);
		smu4Medium->Draw("e1same");
		//c5->SaveAs("plots0p6/signalMuonQua."+plot_format);

		sfourMuFitmu1Pt->GetXaxis()->SetTitle("muon p_{T} [GeV]");
		sfourMuFitmu1Pt->GetYaxis()->SetRangeUser(0,1000);
		sfourMuFitmu1Pt->Draw("e1");
		sfourMuFitmu2Pt->SetMarkerColor(kRed);
		sfourMuFitmu2Pt->Draw("e1same");
		sfourMuFitmu3Pt->SetMarkerColor(kBlue);
		sfourMuFitmu3Pt->SetMarkerStyle(24);
		sfourMuFitmu3Pt->Draw("e1same");
		sfourMuFitmu4Pt->SetMarkerColor(kMagenta);
		sfourMuFitmu4Pt->SetMarkerStyle(24);
		sfourMuFitmu4Pt->Draw("e1same");
		//c5->SaveAs("plots0p6/signalPt."+plot_format);

		hmu1Medium->GetXaxis()->SetTitle("Medium muon");
		hmu1Medium->GetYaxis()->SetRangeUser(0,700);
		hmu1Medium->Draw("e1");
		hmu2Medium->SetMarkerColor(kRed);
		hmu2Medium->Draw("e1same");
		hmu3Medium->SetMarkerColor(kBlue);
		hmu3Medium->SetMarkerStyle(24);
		hmu3Medium->Draw("e1same");
		hmu4Medium->SetMarkerColor(kMagenta);
		hmu4Medium->SetMarkerStyle(24);
		hmu4Medium->Draw("e1same");
		//c5->SaveAs("plots0p6/bkgMuonQua."+plot_format);

		hfourMuFitmu1Pt->GetXaxis()->SetTitle("muon p_{T} [GeV]");
		hfourMuFitmu1Pt->GetYaxis()->SetRangeUser(0,500);
		hfourMuFitmu1Pt->Draw("e1");
		hfourMuFitmu2Pt->SetMarkerColor(kRed);
		hfourMuFitmu2Pt->Draw("e1same");
		hfourMuFitmu3Pt->SetMarkerColor(kBlue);
		hfourMuFitmu3Pt->SetMarkerStyle(24);
		hfourMuFitmu3Pt->Draw("e1same");
		hfourMuFitmu4Pt->SetMarkerColor(kMagenta);
		hfourMuFitmu4Pt->SetMarkerStyle(24);
		hfourMuFitmu4Pt->Draw("e1same");
		//c5->SaveAs("plots0p6/bkgPt."+plot_format);

		TCanvas *c6 = new TCanvas("c6","c6");
		nVertices->Sumw2();
		nVertices->Draw();
		//c6->SaveAs("plots0p6/nVertices."+plot_format);
		nVertices_2012->SetMarkerColor(kRed);
		nVertices_2012->Scale((float)nVertices->Integral()/nVertices_2012->Integral());
		nVertices_2012->GetXaxis()->SetRangeUser(0,70);
		nVertices_2012->Draw();
		nVertices->Draw("same");
		//c6->SaveAs("plots0p6/nVerticesAnd2012."+plot_format);
		nVertices_fourmuon_aftercut->Draw("e1");
		//c6->SaveAs("plots0p6/nVertices_aftercuts."+plot_format);


		nVertices_fourmuon->SetMarkerColor(kRed);
		nVertices_fourmuon->Sumw2();
		nVertices_fourmuon->Divide(nVertices);
		nVertices_fourmuon->GetYaxis()->SetTitle("Normalized N_{4-muon candidates} (13-23 GeV)");
		nVertices_fourmuon->GetYaxis()->SetLabelSize(0.02);
		nVertices_fourmuon->GetYaxis()->SetRangeUser(-0.01,0.6);
		nVertices_fourmuon->GetXaxis()->SetRangeUser(0,90);
		nVertices_fourmuon->Draw("e1");
		nVertices_fourmuon_aftercut->Sumw2();
		nVertices_fourmuon_aftercut->SetMarkerStyle(24);
		nVertices_fourmuon_aftercut->SetMarkerColor(kBlue);
		nVertices_fourmuon_aftercut->Divide(nVertices);
		nVertices_fourmuon_aftercut->GetYaxis()->SetTitle("Normalized N_{4-muon candidates} (13-23 GeV)");
		nVertices_fourmuon_aftercut->Draw("e1same");
		//c6->SaveAs("plots0p6/4muonCandsPerEvt_compare."+plot_format);
		nVertices_fourmuon_aftercut->GetYaxis()->SetLabelSize(0.02);
		nVertices_fourmuon_aftercut->GetYaxis()->SetRangeUser(-0.002,0.04);
		nVertices_fourmuon_aftercut->GetXaxis()->SetRangeUser(0,90);
		nVertices_fourmuon_aftercut->Draw("e1");
		//c6->SaveAs("plots0p6/4muonCandsPerEvt_aftercut."+plot_format);


		/*	TCanvas *c2 = new TCanvas("c2","c2");
			TGraph *gSigEffMuonPtCut = new TGraph(5, muonPtCut, sigEffMuonPtCut);
			gSigEffMuonPtCut->SetMarkerColor(kRed);
			gSigEffMuonPtCut->SetLineColor(kRed);
			gSigEffMuonPtCut->SetMarkerStyle(20);
			gSigEffMuonPtCut->Draw("APL");
			gSigEffMuonPtCut->GetXaxis()->SetTitle("muon p_{T} [GeV]");
			gSigEffMuonPtCut->GetYaxis()->SetTitle("Efficiency");
			gSigEffMuonPtCut->GetYaxis()->SetRangeUser(0,1.2);
			TGraph *gbkRejMuonPtCut = new TGraph(5, muonPtCut, bkRejMuonPtCut);
			gbkRejMuonPtCut->SetMarkerColor(kBlue);
			gbkRejMuonPtCut->SetLineColor(kBlue);
			gbkRejMuonPtCut->SetMarkerStyle(22);
			gbkRejMuonPtCut->Draw("PL");
			c2->SaveAs("plots0p6/EffMuonPtCut."+plot_format);

			TGraph *gSignificanceMuonPtCut = new TGraph(5, muonPtCut, significanceMuonPtCut);
			gSignificanceMuonPtCut->Draw("APL");
			gSignificanceMuonPtCut->GetXaxis()->SetTitle("muon p_{T} [GeV]");
			gSignificanceMuonPtCut->GetYaxis()->SetTitle("S/sqrt(B)");
			c2->SaveAs("plots0p6/significanceMuonPtCut."+plot_format);
			*/
		TCanvas *c3 = new TCanvas("c3","c3");
		TGraph *gSigEffVProbCut = new TGraph(14, vProbCut, sigEffVProbCut);
		gSigEffVProbCut->SetMarkerColor(kRed);
		gSigEffVProbCut->SetLineColor(kRed);
		gSigEffVProbCut->SetMarkerStyle(20);
		gSigEffVProbCut->Draw("APL");
		gSigEffVProbCut->GetXaxis()->SetTitle("4 muons vertex fit probability [GeV]");
		gSigEffVProbCut->GetYaxis()->SetTitle("Efficiency");
		gSigEffVProbCut->GetYaxis()->SetRangeUser(0,1.2);

		TGraph *gbkRejVProbCut = new TGraph(14, vProbCut, bkRejVProbCut);
		gbkRejVProbCut->SetMarkerColor(kBlue);
		gbkRejVProbCut->SetLineColor(kBlue);
		gbkRejVProbCut->SetMarkerStyle(22);
		gbkRejVProbCut->Draw("PL");
		//c3->SaveAs("plots0p6/EffVProbCut."+plot_format);

		TGraph *gSignificanceVProbCut = new TGraph(14, vProbCut, significanceVProbCut);
		gSignificanceVProbCut->Draw("APL");
		gSignificanceVProbCut->GetXaxis()->SetTitle("4 muons vertex fit probability [GeV]");
		gSignificanceVProbCut->GetYaxis()->SetTitle("S/sqrt(B)");
		gSignificanceVProbCut->GetXaxis()->SetRangeUser(0.001,0.1);
		//c3->SaveAs("plots0p6/significanceVProbCut."+plot_format);

		TCanvas *c4 = new TCanvas("c4","c4");
		TGraph *gSigEffQuaCut = new TGraph(3, quaCut, sigEffQuaCut);
		gSigEffQuaCut->SetMarkerColor(kRed);
		gSigEffQuaCut->SetLineColor(kRed);
		gSigEffQuaCut->SetMarkerStyle(20);
		gSigEffQuaCut->Draw("APL");
		gSigEffQuaCut->GetXaxis()->SetTitle("number of muons passing ID cut");
		gSigEffQuaCut->GetYaxis()->SetTitle("Efficiency");
		gSigEffQuaCut->GetYaxis()->SetRangeUser(0,1.2);
		gSigEffQuaCut->GetXaxis()->SetRangeUser(0,3);
		gSigEffQuaCut->GetXaxis()->SetNdivisions(3);
		TGraph *gbkRejQuaCut = new TGraph(3, quaCut, bkRejQuaCut);
		gbkRejQuaCut->SetMarkerColor(kBlue);
		gbkRejQuaCut->SetLineColor(kBlue);
		gbkRejQuaCut->SetMarkerStyle(22);
		gbkRejQuaCut->Draw("PL");
		TGraph *gSigEffQuaCut_m = new TGraph(3, quaCut_m, sigEffQuaCut_m);
		gSigEffQuaCut_m->SetMarkerColor(kRed);
		gSigEffQuaCut_m->SetLineColor(kRed);
		gSigEffQuaCut_m->SetMarkerStyle(24);
		gSigEffQuaCut_m->Draw("PL");
		gSigEffQuaCut_m->GetYaxis()->SetRangeUser(0,1.2);
		gSigEffQuaCut_m->GetXaxis()->SetRangeUser(0,3);
		gSigEffQuaCut_m->GetXaxis()->SetNdivisions(3);
		TGraph *gbkRejQuaCut_m = new TGraph(3, quaCut_m, bkRejQuaCut_m);
		gbkRejQuaCut_m->SetMarkerColor(kBlue);
		gbkRejQuaCut_m->SetLineColor(kBlue);
		gbkRejQuaCut_m->SetMarkerStyle(26);
		gbkRejQuaCut_m->Draw("PL");
		//c4->SaveAs("plots0p6/EffQuaCut."+plot_format);

		TGraph *gSignificanceQuaCut = new TGraph(3, quaCut, significanceQuaCut);
		gSignificanceQuaCut->Draw("APL");
		gSignificanceQuaCut->GetXaxis()->SetTitle("number of muons passing ID cut");
		gSignificanceQuaCut->GetYaxis()->SetTitle("S/sqrt(B)");
		gSignificanceQuaCut->GetXaxis()->SetRangeUser(0,3);
		gSignificanceQuaCut->GetXaxis()->SetNdivisions(3);
		TGraph *gSignificanceQuaCut_m = new TGraph(3, quaCut_m, significanceQuaCut_m);
		gSignificanceQuaCut_m->SetMarkerStyle(24);
		gSignificanceQuaCut_m->Draw("PL");
		gSignificanceQuaCut_m->GetXaxis()->SetRangeUser(0,3);
		gSignificanceQuaCut_m->GetXaxis()->SetNdivisions(3);
		//c4->SaveAs("plots0p6/significanceQuaCut."+plot_format);

	}
