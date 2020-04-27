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
        //gROOT->ProcessLine(".x ~/setTDRStyle.C");
	const bool blind_signal = true;
	TString plot_format = "png";
	gStyle->SetOptStat(kFALSE);
        bool debug = false;
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
        TH1F *Ypt = new TH1F("Ypt","",20,0,40);
        TH1F *mix_Ypt = new TH1F("mix_Ypt","",20,0,40);
        TH1F *Yeta = new TH1F("Yeta","",20,-5.0,5.0);
        TH1F *mix_Yeta = new TH1F("mix_Yeta","",20,-5.0,5.0);
        TH1F *Yphi = new TH1F("Yphi","",20,-5.0,5.0);
        TH1F *mix_Yphi = new TH1F("mix_Yphi","",20,-5,5);
        Ypt->Sumw2();
        mix_Ypt->Sumw2();
        Yeta->Sumw2();
        mix_Yeta->Sumw2();
        Yphi->Sumw2();
        mix_Yphi->Sumw2();
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

	TH1F *hfourMuFitmu1Pt = new TH1F("hfourMuFitmu1Pt","",34,3,20);
	TH1F *hfourMuFitmu2Pt = new TH1F("hfourMuFitmu2Pt","",34,3,20);
	TH1F *hfourMuFitmu3Pt = new TH1F("hfourMuFitmu3Pt","",34,3,20);
	TH1F *hfourMuFitmu4Pt = new TH1F("hfourMuFitmu4Pt","",34,3,20);

        TH1F *hfourMuFitmu1Eta = new TH1F("hfourMuFitmu1Eta","",40,-2.5,2.5);
        TH1F *hfourMuFitmu2Eta = new TH1F("hfourMuFitmu2Eta","",40,-2.5,2.5);
        TH1F *hfourMuFitmu3Eta = new TH1F("hfourMuFitmu3Eta","",40,-2.5,2.5);
        TH1F *hfourMuFitmu4Eta = new TH1F("hfourMuFitmu4Eta","",40,-2.5,2.5);

        TH1F *hfourMuFitmu1Phi = new TH1F("hfourMuFitmu1Phi","",40,-2.5,2.5);
        TH1F *hfourMuFitmu2Phi = new TH1F("hfourMuFitmu2Phi","",40,-2.5,2.5);
        TH1F *hfourMuFitmu3Phi = new TH1F("hfourMuFitmu3Phi","",40,-2.5,2.5);
        TH1F *hfourMuFitmu4Phi = new TH1F("hfourMuFitmu4Phi","",40,-2.5,2.5);
        
        hfourMuFitmu1Pt->Sumw2();
        hfourMuFitmu2Pt->Sumw2();
        hfourMuFitmu3Pt->Sumw2();
        hfourMuFitmu4Pt->Sumw2();
        hfourMuFitmu1Eta->Sumw2();
        hfourMuFitmu2Eta->Sumw2();
        hfourMuFitmu3Eta->Sumw2();
        hfourMuFitmu4Eta->Sumw2();
        hfourMuFitmu1Phi->Sumw2();
        hfourMuFitmu2Phi->Sumw2();
        hfourMuFitmu3Phi->Sumw2();
        hfourMuFitmu4Phi->Sumw2();
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
	TH1F *h_mix_fourMuFitmu1Pt = new TH1F("h_mix_fourMuFitmu1Pt","",34,3,20);
	TH1F *h_mix_fourMuFitmu2Pt = new TH1F("h_mix_fourMuFitmu2Pt","",34,3,20);
	TH1F *h_mix_fourMuFitmu3Pt = new TH1F("h_mix_fourMuFitmu3Pt","",34,3,20);
	TH1F *h_mix_fourMuFitmu4Pt = new TH1F("h_mix_fourMuFitmu4Pt","",34,3,20);
 
        TH1F *h_mix_fourMuFitmu1Eta = new TH1F("h_mix_fourMuFitmu1Eta","",40,-2.5,2.5);
        TH1F *h_mix_fourMuFitmu2Eta = new TH1F("h_mix_fourMuFitmu2Eta","",40,-2.5,2.5);
        TH1F *h_mix_fourMuFitmu3Eta = new TH1F("h_mix_fourMuFitmu3Eta","",40,-2.5,2.5);
        TH1F *h_mix_fourMuFitmu4Eta = new TH1F("h_mix_fourMuFitmu4Eta","",40,-2.5,2.5);

        TH1F *h_mix_fourMuFitmu1Phi = new TH1F("h_mix_fourMuFitmu1Phi","",40,-2.5,2.5);
        TH1F *h_mix_fourMuFitmu2Phi = new TH1F("h_mix_fourMuFitmu2Phi","",40,-2.5,2.5);
        TH1F *h_mix_fourMuFitmu3Phi = new TH1F("h_mix_fourMuFitmu3Phi","",40,-2.5,2.5);
        TH1F *h_mix_fourMuFitmu4Phi = new TH1F("h_mix_fourMuFitmu4Phi","",40,-2.5,2.5);
              
        float h_mu_trg_dR_bins[19]= {0.0,0.001,0.003,0.005,0.007,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.25,0.3};
        Int_t  binnum = sizeof(h_mu_trg_dR_bins)/sizeof(Float_t) - 1;
        TH1F *h_mu_trg_dR = new TH1F("h_mu_trg_dR","h_mu_trg_dR", binnum, h_mu_trg_dR_bins);
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
        std::ofstream Muonsoutput;
        Muonsoutput.open("Muonsoutput.txt");

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
			std::array<float, 4> s = {fourMuFit_mu1Pt->at(i), fourMuFit_mu2Pt->at(i), fourMuFit_mu3Pt->at(i), fourMuFit_mu4Pt->at(i)};
			std::sort(s.begin(), s.end());		//sort from small to big
			bool pass_trigger_pt = s.at(1) > 2 && s.at(2) > 3.5 && s.at(3) > 5;
			//if(pass_trigger_pt == false) std::cout << "mu1Pt = " << fourMuFit_mu1Pt->at(i)  << "mu2Pt = " << fourMuFit_mu2Pt->at(i)  << "mu3Pt = " << fourMuFit_mu3Pt->at(i)  << "mu4Pt = " << fourMuFit_mu4Pt->at(i) << std::endl;
                        if (fourMuFit_Mass->at(i)>13 && fourMuFit_Mass->at(i)<23) nVertices_fourmuon->Fill(numPrimaryVertices);
			if (
					fourMuFit_VtxProb->at(i)>0.05
					&& (mu3_Medium->at(i) + mu4_Medium->at(i)) >= 2
					&& mu34.M()< 9.2
					&& fourMuFit_mu1Pt->at(i) >= muonPtCut[2] && fourMuFit_mu2Pt->at(i) >= muonPtCut[2] 
                                        && fourMuFit_mu3Pt->at(i) >= muonPtCut[2] && fourMuFit_mu4Pt->at(i) >= muonPtCut[2]
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
                                Muonsoutput << event << "\n";
                                Muonsoutput << fourMuFit_mu1Pt->at(i) << " "<< fourMuFit_mu2Pt->at(i) << " "<< fourMuFit_mu3Pt->at(i) << " "<< fourMuFit_mu4Pt->at(i) << "\n";
                                Muonsoutput << fourMuFit_mu1Eta->at(i) << " "<< fourMuFit_mu2Eta->at(i) << " "<< fourMuFit_mu3Eta->at(i) << " "<< fourMuFit_mu4Eta->at(i) << "\n";
                                Muonsoutput << fourMuFit_mu1Phi->at(i) << " "<< fourMuFit_mu2Phi->at(i) << " "<< fourMuFit_mu3Phi->at(i) << " "<< fourMuFit_mu4Phi->at(i) << "\n";
				hfourMuMass_aftercut->Fill(fourMuFit_Mass->at(i));
				hfourMuMass_aftercut_smallrange->Fill(fourMuFit_Mass->at(i));
				fourMumassoutput << fourMuFit_Mass->at(i)<< "\n";
				Ymumu2D_aftercut->Fill(mu12.Eta(),mu34.Eta());
				Ymumu2DBoost_aftercut->Fill(mu12boost.Eta(),mu34boost.Eta());
                                Ypt->Fill(mu12.Pt());
                                Yeta->Fill(mu12.Eta());
                                Yphi->Fill(mu12.Phi());
				if (fabs(mu12.Eta())<0.2) Ypt1->Fill(mu12.Pt());
				if (fabs(mu12.Eta())<1.1 && fabs(mu12.Eta())>0.9) Ypt2->Fill(mu12.Pt());
				if (fabs(mu12.Eta())<2.1 && fabs(mu12.Eta())>1.9) Ypt3->Fill(mu12.Pt());
				if (fabs(mu12.Eta())<3.1 && fabs(mu12.Eta())>2.9) Ypt4->Fill(mu12.Pt());
				Y2Dpteta->Fill(mu12.Eta(),mu12.Pt());
				Y2Dptphi->Fill(mu12.Phi(),mu12.Pt()); 
				Y2Detaphi->Fill(mu12.Eta(),mu12.Phi());
                                hfourMuFitmu1Pt->Fill(fourMuFit_mu1Pt->at(i));    
                                hfourMuFitmu2Pt->Fill(fourMuFit_mu2Pt->at(i));
                                hfourMuFitmu3Pt->Fill(fourMuFit_mu3Pt->at(i));
                                hfourMuFitmu4Pt->Fill(fourMuFit_mu4Pt->at(i));
                                hfourMuFitmu1Eta->Fill(fourMuFit_mu1Eta->at(i));
                                hfourMuFitmu2Eta->Fill(fourMuFit_mu2Eta->at(i));
                                hfourMuFitmu3Eta->Fill(fourMuFit_mu3Eta->at(i));
                                hfourMuFitmu4Eta->Fill(fourMuFit_mu4Eta->at(i));
                                hfourMuFitmu1Phi->Fill(fourMuFit_mu1Phi->at(i));
                                hfourMuFitmu2Phi->Fill(fourMuFit_mu2Phi->at(i));
                                hfourMuFitmu3Phi->Fill(fourMuFit_mu3Phi->at(i));
                                hfourMuFitmu4Phi->Fill(fourMuFit_mu4Phi->at(i));

                                float min_mu1_trg_dR = 111;
                                float min_mu2_trg_dR = 111;
                                float min_mu3_trg_dR = 111;
                                float min_mu4_trg_dR = 111;
                                for (unsigned l=0; l<mu1_trg_dR->size();l++)
                                {
                                   if (min_mu1_trg_dR > mu1_trg_dR->at(l)) min_mu1_trg_dR=mu1_trg_dR->at(l);
                                   } 
                                for (unsigned l=0; l<mu2_trg_dR->size();l++)  
                                 {
                                  if (min_mu2_trg_dR > mu2_trg_dR->at(l)) min_mu2_trg_dR=mu2_trg_dR->at(l);	
                                    }
                                for (unsigned l=0; l<fourMuFit_mu3_trg_dR->size();l++)
                                 {
                                  if (min_mu3_trg_dR > fourMuFit_mu3_trg_dR->at(l)) min_mu3_trg_dR=fourMuFit_mu3_trg_dR->at(l);
                                    }
                                for (unsigned l=0; l<fourMuFit_mu4_trg_dR->size();l++)
                                 {
                                  if (min_mu4_trg_dR > fourMuFit_mu4_trg_dR->at(l)) min_mu4_trg_dR=fourMuFit_mu4_trg_dR->at(l);
                                    }
                                h_mu_trg_dR->Fill(min_mu1_trg_dR);
                                h_mu_trg_dR->Fill(min_mu2_trg_dR);
                                if (min_mu4_trg_dR>min_mu3_trg_dR) h_mu_trg_dR->Fill(min_mu3_trg_dR);
                                if (min_mu4_trg_dR<min_mu3_trg_dR) h_mu_trg_dR->Fill(min_mu4_trg_dR);
                                 
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
        std::ofstream mix_Muonsoutput;
        mix_Muonsoutput.open("mix_Muonsoutput.txt");
	//for (unsigned j=0; j<1000;j++){ 
	//	for (unsigned k=0; k<1000;k++)
       cout<<"mu34_p4_vector.size()"<<mu34_p4_vector.size()<<endl;
       cout<<"mu12_p4_vector.size()"<<mu12_p4_vector.size()<<endl;
	for (unsigned j=0; j<mu34_p4_vector.size();j++){ 
		for (unsigned k=0; k<mu12_p4_vector.size();k++){
                         if (j==k) continue;
			
            
            Double_t DeltaPhi =(fabs(mu12_p4_vector[k].Phi() - mu12_p4_vector[j].Phi()));
            Double_t DeltaY =(fabs(mu12_p4_vector[k].Rapidity() - mu12_p4_vector[j].Rapidity()));
            Double_t DeltaR=sqrt((DeltaPhi*DeltaPhi)+ (DeltaY*DeltaY));
            if (DeltaR>0.3) continue;
            //for (unsigned k=j; k<j+1000 && k<mu12_p4_vector.si
			//if ( fabs(mu12boost_p4_vector[k].Vect().DeltaR( mu12boost_p4_vector[j].Vect())) > 1.5) continue;
			//if ( fabs(mu12_p4_vector[k].Vect().DeltaR( mu12_p4_vector[j].Vect())) > 0.3) continue;
			
			//if ( fabs(mu12_p4_vector[k].Vect().DeltaR( mu12_p4_vector[j].Vect())) > 0.6) continue;
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

                        mix_Muonsoutput << event << "\n";
                        mix_Muonsoutput << mu1_p4_vector[k].Pt() << " "<< mu2_p4_vector[k].Pt() << " "<< mu3_p4_vector[j].Pt() << " "<< mu4_p4_vector[j].Pt() << "\n";
                        mix_Muonsoutput << mu1_p4_vector[k].Eta() << " "<< mu2_p4_vector[k].Eta() << " "<< mu3_p4_vector[j].Eta() << " "<< mu4_p4_vector[j].Eta() << "\n";
                        mix_Muonsoutput << mu1_p4_vector[k].Phi() << " "<< mu2_p4_vector[k].Phi() << " "<< mu3_p4_vector[j].Phi() << " "<< mu4_p4_vector[j].Phi() << "\n";
			TLorentzVector mixFourMu;
			mixFourMu=mu12_p4_vector[k]+mu34_p4_vector[j];
			//std::cout<<mixFourMu.Pt()<<" "<<mixFourMu.Pz()<<" "<<mixFourMu.M()<<std::endl;
			//mixFourMu.Boost(fourMuFit.BoostVector());
			//std::cout<<mixFourMu.Pt()<<" "<<mixFourMu.Pz()<<" "<<mixFourMu.M()<<std::endl;
                        h_mix_fourMuFitmu1Pt->Fill(mu1_p4_vector[k].Pt());
                        h_mix_fourMuFitmu2Pt->Fill(mu2_p4_vector[k].Pt());
                        h_mix_fourMuFitmu3Pt->Fill(mu3_p4_vector[j].Pt());
                        h_mix_fourMuFitmu4Pt->Fill(mu4_p4_vector[j].Pt());

                        h_mix_fourMuFitmu1Eta->Fill(mu1_p4_vector[k].Eta());
                        h_mix_fourMuFitmu2Eta->Fill(mu2_p4_vector[k].Eta());
                        h_mix_fourMuFitmu3Eta->Fill(mu3_p4_vector[j].Eta());
                        h_mix_fourMuFitmu4Eta->Fill(mu4_p4_vector[j].Eta());

                        h_mix_fourMuFitmu1Phi->Fill(mu1_p4_vector[k].Phi());
                        h_mix_fourMuFitmu2Phi->Fill(mu2_p4_vector[k].Phi());
                        h_mix_fourMuFitmu3Phi->Fill(mu3_p4_vector[j].Phi());
                        h_mix_fourMuFitmu4Phi->Fill(mu4_p4_vector[j].Phi());
                        mix_Ypt->Fill(mu12_p4_vector[k].Pt());
                        mix_Yeta->Fill(mu12_p4_vector[k].Eta());
                        mix_Yphi->Fill(mu12_p4_vector[k].Phi());
			
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
		c1->SaveAs("plots0p3/signalFourMuMass."+plot_format);


		//before cuts, mix and mix zoom in
		hfourMuMass_mix->SetMarkerStyle(24);
		hfourMuMass_mix->SetMarkerColor(kRed);
		hfourMuMass_mix->Draw("e1");
		//c1->SaveAs("plots0p3/hfourMuMass_mix."+plot_format);
		hfourMuMass_mix_smallrange->SetMarkerStyle(24);
		hfourMuMass_mix_smallrange->SetMarkerColor(kBlue);
		hfourMuMass_mix_smallrange->GetYaxis()->SetRangeUser(0,200);
		hfourMuMass_mix_smallrange->Draw("e1");
		c1->SaveAs("plots0p3/hfourMuMass_mix_smallrange."+plot_format);

		//before cuts, compare original vs mix
		hfourMuMass->Draw("e1");
		c1->SaveAs("plots0p3/hfourMuMass."+plot_format);
		hfourMuMass_mix->Draw("e1same");
		c1->SaveAs("plots0p3/hfourMuMass_origVSmix."+plot_format);

		//after cuts, orignal
		hfourMuMass_aftercut->Draw("e1");
		c1->SaveAs("plots0p3/hfourMuMass_afterCut."+plot_format);
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
		c1->SaveAs("plots0p3/hfourMuMass_afterCut_smallrange."+plot_format);

		//after cuts, mix
		hfourMuMass_mix_aftercut->SetMarkerStyle(24);
		hfourMuMass_mix_aftercut->SetMarkerColor(kRed);
		hfourMuMass_mix_aftercut->Draw("e1");
		c1->SaveAs("plots0p3/hfourMuMass_mix_afterCut."+plot_format);

		//after cuts, compare original vs mix
		hfourMuMass_aftercut->Draw("e1");
		hfourMuMass_mix_aftercut->Draw("e1same");
		c1->SaveAs("plots0p3/hfourMuMass_afterCut_origVSmix."+plot_format);

		//after cuts, mix physics bkg
		hfourMuMass_physkbg_mix->SetMarkerStyle(22);
		hfourMuMass_physkbg_mix->SetMarkerColor(kBlue);
		hfourMuMass_physkbg_mix->Draw("e1");
		c1->SaveAs("plots0p3/hfourMuMass_physkbg_mix."+plot_format);
		hfourMuMass_physkbg_mix_smallrange->SetMarkerStyle(22);
		hfourMuMass_physkbg_mix_smallrange->SetMarkerColor(kBlue);
		hfourMuMass_physkbg_mix_smallrange->Draw("e1");
		c1->SaveAs("plots0p3/hfourMuMass_physkbg_mix_smallrange."+plot_format);

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
		c1->SaveAs("plots0p3/hfourMuMass_afterCut_origVSphyskbg."+plot_format);

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
		c1->SaveAs("plots0p3/hfourMuMass_afterCut_origVSphyskbg_smallrange."+plot_format);

		TFile* out_file = new TFile("plots0p3/fourMuMass.root","RECREATE");
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
		c1->SaveAs("plots0p3/hfourMuMassBoost_afterCut_origVSphyskbg."+plot_format);
		hfourMuMass_aftercut_smallrange->Draw("e1");
		hfourMuMassBoost_physkbg_mix_smallrange->SetMarkerStyle(23);
		hfourMuMassBoost_physkbg_mix_smallrange->SetMarkerColor(kRed);
		hfourMuMassBoost_physkbg_mix_smallrange->Scale(scaleEntries*5/2/hfourMuMassBoost_physkbg_mix_smallrange->Integral());
		hfourMuMassBoost_physkbg_mix_smallrange->Draw("e1same");
		c1->SaveAs("plots0p3/hfourMuMassBoost_afterCut_origVSphyskbg_smallrange."+plot_format);

		TCanvas *c2 = new TCanvas("c2","c2");
		Ymumu2D->Draw("colz");
		c2->SaveAs("plots0p3/Ymumu2D."+plot_format);
		Ymumu2DBoost->Draw("colz");
		c2->SaveAs("plots0p3/Ymumu2DBoost."+plot_format);
		Ymumu2D_aftercut->Draw("colz");
		c2->SaveAs("plots0p3/Ymumu2D_aftercut."+plot_format);
		Ymumu2DBoost_aftercut->Draw("colz");
		c2->SaveAs("plots0p3/Ymumu2DBoost_aftercut."+plot_format);
		Ymumu2D_zerobias_mix->Draw("colz");
		c2->SaveAs("plots0p3/Ymumu2D_zerobias_mix."+plot_format);
		Ymumu2D_physkbg_mix->Draw("colz");
		c2->SaveAs("plots0p3/Ymumu2D_physkbg_mix."+plot_format);
		Ymumu2DBoost_physkbg_mix->Draw("colz");
		c2->SaveAs("plots0p3/Ymumu2DBoost_physkbg_mix."+plot_format);
		h_number_of_Y->Draw("e1");
		c2->SetLogy();
		c2->SaveAs("plots0p3/number_of_Y."+plot_format);

		TCanvas *c3 = new TCanvas("c3","c3");
		//Ypt1->Scale(1./Ypt1->Integral());
		Ypt1->Draw("e1");
		//c3->SaveAs("plots0p3/Ypt1."+plot_format);
		//Ypt2->Scale(Ypt1->Integral()/Ypt2->Integral());
		//Ypt2->SetMarkerColor(kRed);
		Ypt2->Draw("e1");
		//c3->SaveAs("plots0p3/Ypt2."+plot_format);
		//Ypt3->Scale(Ypt1->Integral()/Ypt3->Integral());
		//Ypt3->SetMarkerColor(kBlue);
		Ypt3->Draw("e1");
		//c3->SaveAs("plots0p3/Ypt3."+plot_format);
		//Ypt4->Scale(Ypt1->Integral()/Ypt4->Integral());
		//Ypt4->SetMarkerColor(kMagenta);
		Ypt4->Draw("e1");
		c3->SaveAs("plots0p3/Ypt4."+plot_format);
		Y2Dpteta->Draw("colz");
		c3->SaveAs("plots0p3/Y2Dpteta."+plot_format);
		Y2Dptphi->Draw("colz");
		c3->SaveAs("plots0p3/Y2Dptphi."+plot_format);
		Y2Detaphi->Draw("colz");
		c3->SaveAs("plots0p3/Y2Detaphi."+plot_format);
		mu12mass->Draw();
		mu12massEBE->SetLineColor(kRed);
		mu12massEBE->Draw("same");
		c3->SaveAs("plots0p3/Ymass."+plot_format);
		mu34massBkg->Draw("e1");
		mu34massBkg->GetYaxis()->SetRangeUser(0,50);
		c3->SaveAs("plots0p3/mu34massBkg."+plot_format);
		mu34massBkgH->Draw("e1");
		mu34massBkgH->GetYaxis()->SetRangeUser(0,50);
		c3->SaveAs("plots0p3/mu34massBkgH."+plot_format);
		mu34mass->Draw("e1");
		mu34mass->GetYaxis()->SetRangeUser(0,50);
		c3->SaveAs("plots0p3/mu34mass."+plot_format);
		mu34massBkg->SetMarkerColor(kRed);
		mu34massBkg->Draw("e1same");
		mu34massBkgH->SetMarkerColor(kBlue);
		mu34massBkgH->Draw("e1same");
		c3->SaveAs("plots0p3/mu34massOverlay."+plot_format);
		h_fourMu_pt_order->Draw("e1");
		h_fourMu_pt_order->SetMinimum(0);
		//c3->SaveAs("plots0p3/ourMu_pt_order."+plot_format);
		

		TCanvas *c4 = new TCanvas("c4","c4");
		hnCand_aftercut->SetMarkerColor(kRed);
		hnCand_aftercut->Draw("e1");
		//hnCand->Scale((float)hnCand_aftercut->Integral()/hnCand->Integral());
		//hnCand->SetMarkerSize(0.7);
		//hnCand->Draw("same");
		c4->SaveAs("plots0p3/nCand."+plot_format);
		shnCand_aftercut->SetMarkerColor(kRed);
		shnCand_aftercut->Draw("e1");
		shnCand->Scale((float)shnCand_aftercut->Integral()/shnCand->Integral());
		shnCand->SetMarkerSize(0.7);
		shnCand->Draw("same");
		c4->SaveAs("plots0p3/nCand_signal."+plot_format);
		vProbMix->SetMarkerColor(kBlue);
		vProbMix->Draw("e1");
		c4->SetLogy();
		vProbSig->Scale((float)vProbMix->Integral()/vProbSig->Integral());
		vProbSig->SetMarkerColor(kRed);
		vProbSig->Draw("same");
		//c4->SaveAs("plots0p3/vProb_SandB."+plot_format);


		TCanvas *c5 = new TCanvas("c5","c5");	
		hfourMuMass_aftercut_lowLumi->SetMarkerColor(kRed);
		hfourMuMass_aftercut_lowLumi->Draw("e1");
		hfourMuMass_aftercut_highLumi->Scale((float)hfourMuMass_aftercut_lowLumi->Integral()/hfourMuMass_aftercut_highLumi->Integral());
		hfourMuMass_aftercut_highLumi->SetMarkerColor(kBlue);
		hfourMuMass_aftercut_highLumi->Draw("e1same");
		c5->SaveAs("plots0p3/compareMassvsLumi."+plot_format);

		hfourMuMass_mix_lowLumi->SetMarkerColor(kRed);
		hfourMuMass_mix_lowLumi->Draw("e1");
		hfourMuMass_mix_midLumi->Scale((float)hfourMuMass_mix_lowLumi->Integral()/hfourMuMass_mix_midLumi->Integral());
		hfourMuMass_mix_midLumi->Draw("e1same");
		hfourMuMass_mix_highLumi->Scale((float)hfourMuMass_mix_lowLumi->Integral()/hfourMuMass_mix_highLumi->Integral());
		hfourMuMass_mix_highLumi->SetMarkerColor(kBlue);
		hfourMuMass_mix_highLumi->Draw("e1same");
		c5->SaveAs("plots0p3/compareMassvsLumi_pileupbkg."+plot_format);

		hfourMuMass_aftercut_lowLumi_smallrange->SetMarkerColor(kRed);
		hfourMuMass_aftercut_lowLumi_smallrange->Draw("e1");
		hfourMuMass_aftercut_highLumi_smallrange->Scale((float)hfourMuMass_aftercut_lowLumi_smallrange->Integral()/hfourMuMass_aftercut_highLumi_smallrange->Integral());
		hfourMuMass_aftercut_highLumi_smallrange->SetMarkerColor(kBlue);
		hfourMuMass_aftercut_highLumi_smallrange->Draw("e1same");
		c5->SaveAs("plots0p3/compareMassvsLumi_smallrange."+plot_format);

		hfourMuMass_mix_lowLumi_smallrange->SetMarkerColor(kRed);
		hfourMuMass_mix_lowLumi_smallrange->Draw("e1");
		hfourMuMass_mix_midLumi_smallrange->Scale((float)hfourMuMass_mix_lowLumi_smallrange->Integral()/hfourMuMass_mix_midLumi_smallrange->Integral());
		hfourMuMass_mix_midLumi_smallrange->Draw("e1same");	
		hfourMuMass_mix_highLumi_smallrange->Scale((float)hfourMuMass_mix_lowLumi_smallrange->Integral()/hfourMuMass_mix_highLumi_smallrange->Integral());
		hfourMuMass_mix_highLumi_smallrange->SetMarkerColor(kBlue);
		hfourMuMass_mix_highLumi_smallrange->Draw("e1same");
		c5->SaveAs("plots0p3/compareMassvsLumi_pileupbkg_smallrange."+plot_format);

		TCanvas *c6 = new TCanvas("c6","c6");
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
		c6->SaveAs("plots0p3/signalMuonQua."+plot_format);

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
	    c6->SaveAs("plots0p3/bkgMuonQua."+plot_format);
                
		TCanvas *c7 = new TCanvas("c7","c7");
                float maximum_range;
                bool Ratio = true;
                hfourMuFitmu1Pt->GetXaxis()->SetTitle("Muon p_{T} [GeV]");
                hfourMuFitmu1Pt->Scale(1/hfourMuFitmu1Pt->Integral());
                h_mix_fourMuFitmu1Pt->Scale(1/h_mix_fourMuFitmu1Pt->Integral());
                maximum_range = max(hfourMuFitmu1Pt->GetMaximum(),h_mix_fourMuFitmu1Pt->GetMaximum());
                hfourMuFitmu1Pt->SetMarkerStyle(20);
                hfourMuFitmu1Pt->SetMarkerColor(kBlack);
                hfourMuFitmu1Pt->SetLineColor(kBlack);
                hfourMuFitmu1Pt->SetMinimum(0);
                hfourMuFitmu1Pt->SetMaximum(1.3*maximum_range);
                hfourMuFitmu1Pt->Draw("hist p");                 
                h_mix_fourMuFitmu1Pt->SetMarkerStyle(22);
                h_mix_fourMuFitmu1Pt->SetMarkerColor(kBlue);
                h_mix_fourMuFitmu1Pt->Draw("e1same");
                h_mix_fourMuFitmu1Pt->Draw("hist p same");
                TLegend* leg1 = new TLegend(0.7,0.7,0.9,0.9);
                leg1->AddEntry(hfourMuFitmu1Pt,"muons","lep");
                leg1->AddEntry(h_mix_fourMuFitmu1Pt,"Mixed muons","lep");
                leg1->Draw("same");
                if (Ratio) 
                {
                 c7->SetBottomMargin(0.3);
                 hfourMuFitmu1Pt->GetXaxis()->SetTitleSize(0);
                 hfourMuFitmu1Pt->GetXaxis()->SetLabelSize(0);
                 TH1F *Totratio = (TH1F*)h_mix_fourMuFitmu1Pt->Clone("Totratio");
                 Totratio->Divide(hfourMuFitmu1Pt);
                 Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                 Totratio->GetXaxis()->SetNoExponent(kTRUE);
                 TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 0.93, 0.93);
                 pad->SetTopMargin(0.7);
                 pad->SetRightMargin(0.03);
                 pad->SetFillColor(0);
                 pad->SetGridy(1);
                 pad->SetFillStyle(0);
                 pad->Draw();
                 pad->cd(0);
                 Totratio->GetXaxis()->SetTitle("p_{T}(#mu) [GeV]");
                 Totratio->GetYaxis()->SetTitle("Mixed(#mu)/#mu");
                 Totratio->GetYaxis()->SetTitleOffset(1.2);
                 Totratio->GetXaxis()->SetTitleOffset(1.2);
                 Totratio->GetYaxis()->SetLabelSize(0.03);
                 Totratio->SetMarkerStyle(20);
                 Totratio->SetMarkerSize(1.2);
                 Totratio->SetLineColor(1);
                 Totratio->SetMarkerColor(1);
                 Totratio->GetYaxis()->SetNdivisions(204);
                 Totratio->SetMinimum(0.0);
                 Totratio->SetMaximum(2.0);
                 Totratio->Draw("ep");
                 c7->Modified();
                 c7->Update();
                  }   
                c7->SaveAs("plots0p3/Muon_1_Pt."+plot_format);
                c7->Clear();
                hfourMuFitmu1Eta->GetXaxis()->SetTitle("Muon #eta");
                hfourMuFitmu1Eta->Scale(1/hfourMuFitmu1Eta->Integral());
                h_mix_fourMuFitmu1Eta->Scale(1/h_mix_fourMuFitmu1Eta->Integral());
                maximum_range = max(hfourMuFitmu1Eta->GetMaximum(),h_mix_fourMuFitmu1Eta->GetMaximum());
                hfourMuFitmu1Eta->SetMarkerStyle(20);
                hfourMuFitmu1Eta->SetMarkerColor(kBlack);
                hfourMuFitmu1Eta->SetLineColor(kBlack);
                hfourMuFitmu1Eta->SetMinimum(0);
                hfourMuFitmu1Eta->SetMaximum(1.7*maximum_range);
                hfourMuFitmu1Eta->Draw("hist p");
                h_mix_fourMuFitmu1Eta->SetMarkerStyle(22);
                h_mix_fourMuFitmu1Eta->SetMarkerColor(kBlue);
                h_mix_fourMuFitmu1Eta->Draw("hist p same");
                leg1->Draw("same");
                if (Ratio)
                {
                 c7->SetBottomMargin(0.3);
                 hfourMuFitmu1Eta->GetXaxis()->SetTitleSize(0);
                 hfourMuFitmu1Eta->GetXaxis()->SetLabelSize(0);
                 TH1F *Totratio = (TH1F*)h_mix_fourMuFitmu1Eta->Clone("Totratio");
                 Totratio->Divide(hfourMuFitmu1Eta);
                 Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                 Totratio->GetXaxis()->SetNoExponent(kTRUE);
                 TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.0, 0.93, 0.93);
                 pad1->SetTopMargin(0.7);
                 pad1->SetRightMargin(0.03);
                 pad1->SetFillColor(0);
                 pad1->SetGridy(1);
                 pad1->SetFillStyle(0);
                 pad1->Draw();
                 pad1->cd(0);
                 Totratio->GetXaxis()->SetTitle("#eta(#mu)");
                 Totratio->GetYaxis()->SetTitle("Mixed(#mu)/#mu");
                 Totratio->GetYaxis()->SetTitleOffset(1.2);
                 Totratio->GetXaxis()->SetTitleOffset(1.2);
                 Totratio->GetYaxis()->SetLabelSize(0.03);
                 Totratio->SetMarkerStyle(20);
                 Totratio->SetMarkerSize(1.2);
                 Totratio->SetLineColor(1);
                 Totratio->SetMarkerColor(1);
                 Totratio->GetYaxis()->SetNdivisions(204);
                 Totratio->SetMinimum(0.0);
                 Totratio->SetMaximum(2.0);
                 Totratio->Draw("ep");
                 c7->Modified();
                 c7->Update();
                  }
                c7->SaveAs("plots0p3/Muon_1_Eta."+plot_format);
                c7->Clear();
                hfourMuFitmu1Phi->GetXaxis()->SetTitle("Muon #phi");
                hfourMuFitmu1Phi->Scale(1/hfourMuFitmu1Phi->Integral());
                h_mix_fourMuFitmu1Phi->Scale(1/h_mix_fourMuFitmu1Phi->Integral());
                maximum_range = max(hfourMuFitmu1Phi->GetMaximum(),h_mix_fourMuFitmu1Phi->GetMaximum());
                hfourMuFitmu1Phi->SetMarkerStyle(20);
                hfourMuFitmu1Phi->SetMarkerColor(kBlack);
                hfourMuFitmu1Phi->SetLineColor(kBlack);
                hfourMuFitmu1Phi->SetMinimum(0);
                hfourMuFitmu1Phi->SetMaximum(1.7*maximum_range);
                hfourMuFitmu1Phi->Draw("hist p");
                h_mix_fourMuFitmu1Phi->SetMarkerStyle(22);
                h_mix_fourMuFitmu1Phi->SetMarkerColor(kBlue);
                h_mix_fourMuFitmu1Phi->Draw("hist p same");
                leg1->Draw("same");
                if (Ratio)
                {
                 c7->SetBottomMargin(0.3);
                 hfourMuFitmu1Phi->GetXaxis()->SetTitleSize(0);
                 hfourMuFitmu1Phi->GetXaxis()->SetLabelSize(0);
                 TH1F *Totratio = (TH1F*)h_mix_fourMuFitmu1Phi->Clone("Totratio");
                 Totratio->Divide(hfourMuFitmu1Phi);
                 Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                 Totratio->GetXaxis()->SetNoExponent(kTRUE);
                 TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 0.93, 0.93);
                 pad2->SetTopMargin(0.7);
                 pad2->SetRightMargin(0.03);
                 pad2->SetFillColor(0);
                 pad2->SetGridy(1);
                 pad2->SetFillStyle(0);
                 pad2->Draw();
                 pad2->cd(0);
                 Totratio->GetXaxis()->SetTitle("#Phi(#mu)");
                 Totratio->GetYaxis()->SetTitle("Mixed(#mu)/#mu");
                 Totratio->GetYaxis()->SetTitleOffset(1.2);
                 Totratio->GetXaxis()->SetTitleOffset(1.2);
                 Totratio->GetYaxis()->SetLabelSize(0.03);
                 Totratio->SetMarkerStyle(20);
                 Totratio->SetMarkerSize(1.2);
                 Totratio->SetLineColor(1);
                 Totratio->SetMarkerColor(1);
                 Totratio->GetYaxis()->SetNdivisions(204);
                 Totratio->SetMinimum(0.0);
                 Totratio->SetMaximum(2.0);
                 Totratio->Draw("ep");
                 c7->Modified();
                 c7->Update();
                  }
               
                c7->SaveAs("plots0p3/Muon_1_Phi."+plot_format);
                c7->Clear();
                hfourMuFitmu2Pt->GetXaxis()->SetTitle("Muon p_{T} [GeV]");
                hfourMuFitmu2Pt->Scale(1/hfourMuFitmu2Pt->Integral());
                h_mix_fourMuFitmu2Pt->Scale(1/h_mix_fourMuFitmu2Pt->Integral());
                maximum_range = max(hfourMuFitmu2Pt->GetMaximum(),h_mix_fourMuFitmu2Pt->GetMaximum());
                hfourMuFitmu2Pt->SetMarkerStyle(20);
                hfourMuFitmu2Pt->SetMarkerColor(kBlack);
                hfourMuFitmu2Pt->SetLineColor(kBlack);
                hfourMuFitmu2Pt->SetMinimum(0);
                hfourMuFitmu2Pt->SetMaximum(1.3*maximum_range);
                hfourMuFitmu2Pt->Draw("hist p");
                h_mix_fourMuFitmu2Pt->SetMarkerStyle(22);
                h_mix_fourMuFitmu2Pt->SetMarkerColor(kBlue);
                h_mix_fourMuFitmu2Pt->Draw("hist p same");
                leg1->Draw("same");
                if (Ratio)
                {
                 c7->SetBottomMargin(0.3);
                 hfourMuFitmu2Pt->GetXaxis()->SetTitleSize(0);
                 hfourMuFitmu2Pt->GetXaxis()->SetLabelSize(0);
                 TH1F *Totratio = (TH1F*)h_mix_fourMuFitmu2Pt->Clone("Totratio");
                 Totratio->Divide(hfourMuFitmu2Pt);
                 Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                 Totratio->GetXaxis()->SetNoExponent(kTRUE);
                 TPad *pad3 = new TPad("pad3", "pad3", 0.0, 0.0, 0.93, 0.93);
                 pad3->SetTopMargin(0.7);
                 pad3->SetRightMargin(0.03);
                 pad3->SetFillColor(0);
                 pad3->SetGridy(1);
                 pad3->SetFillStyle(0);
                 pad3->Draw();
                 pad3->cd(0);
                 Totratio->GetXaxis()->SetTitle("p_{T}(#mu) [GeV]");
                 Totratio->GetYaxis()->SetTitle("Mixed(#mu)/#mu");
                 Totratio->GetYaxis()->SetTitleOffset(1.2);
                 Totratio->GetXaxis()->SetTitleOffset(1.2);
                 Totratio->GetYaxis()->SetLabelSize(0.03);
                 Totratio->SetMarkerStyle(20);
                 Totratio->SetMarkerSize(1.2);
                 Totratio->SetLineColor(1);
                 Totratio->SetMarkerColor(1);
                 Totratio->GetYaxis()->SetNdivisions(204);
                 Totratio->SetMinimum(0.0);
                 Totratio->SetMaximum(2.0);
                 Totratio->Draw("ep");
                 c7->Modified();
                 c7->Update();
                  }
                c7->SaveAs("plots0p3/Muon_2_Pt."+plot_format);
                c7->Clear();
                hfourMuFitmu2Eta->GetXaxis()->SetTitle("Muon #eta");
                hfourMuFitmu2Eta->Scale(1/hfourMuFitmu2Eta->Integral());
                h_mix_fourMuFitmu2Eta->Scale(1/h_mix_fourMuFitmu2Eta->Integral());
                maximum_range = max(hfourMuFitmu2Eta->GetMaximum(),h_mix_fourMuFitmu2Eta->GetMaximum());
                hfourMuFitmu2Eta->SetMarkerStyle(20);
                hfourMuFitmu2Eta->SetMarkerColor(kBlack);
                hfourMuFitmu2Eta->SetLineColor(kBlack);
                hfourMuFitmu2Eta->SetMinimum(0);
                hfourMuFitmu2Eta->SetMaximum(1.7*maximum_range);
                hfourMuFitmu2Eta->Draw("hist p"); 
                h_mix_fourMuFitmu2Eta->SetMarkerStyle(22);
                h_mix_fourMuFitmu2Eta->SetMarkerColor(kBlue);
                h_mix_fourMuFitmu2Eta->Draw("hist p same");
                leg1->Draw("same");
                if (Ratio)
                {
                 c7->SetBottomMargin(0.3);
                 hfourMuFitmu2Eta->GetXaxis()->SetTitleSize(0);
                 hfourMuFitmu2Eta->GetXaxis()->SetLabelSize(0);
                 TH1F *Totratio = (TH1F*)h_mix_fourMuFitmu2Eta->Clone("Totratio");
                 Totratio->Divide(hfourMuFitmu2Eta);
                 Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                 Totratio->GetXaxis()->SetNoExponent(kTRUE);
                 TPad *pad4 = new TPad("pad4", "pad4", 0.0, 0.0, 0.93, 0.93);
                 pad4->SetTopMargin(0.7);
                 pad4->SetRightMargin(0.03);
                 pad4->SetFillColor(0);
                 pad4->SetGridy(1);
                 pad4->SetFillStyle(0);
                 pad4->Draw();
                 pad4->cd(0);
                 Totratio->GetXaxis()->SetTitle("#eta(#mu)");
                 Totratio->GetYaxis()->SetTitle("Mixed(#mu)/#mu");
                 Totratio->GetYaxis()->SetTitleOffset(1.2);
                 Totratio->GetXaxis()->SetTitleOffset(1.2);
                 Totratio->GetYaxis()->SetLabelSize(0.03);
                 Totratio->SetMarkerStyle(20);
                 Totratio->SetMarkerSize(1.2);
                 Totratio->SetLineColor(1);
                 Totratio->SetMarkerColor(1);
                 Totratio->GetYaxis()->SetNdivisions(204);
                 Totratio->SetMinimum(0.0);
                 Totratio->SetMaximum(2.0);
                 Totratio->Draw("ep");
                 c7->Modified();
                 c7->Update();
                  }
                c7->SaveAs("plots0p3/Muon_2_Eta."+plot_format);
                c7->Clear();
                hfourMuFitmu2Phi->GetXaxis()->SetTitle("Muon #phi");
                hfourMuFitmu2Phi->Scale(1/hfourMuFitmu2Phi->Integral());
                h_mix_fourMuFitmu2Phi->Scale(1/h_mix_fourMuFitmu2Phi->Integral());
                maximum_range = max(hfourMuFitmu2Phi->GetMaximum(),h_mix_fourMuFitmu2Phi->GetMaximum());
                hfourMuFitmu2Phi->SetMarkerStyle(20);
                hfourMuFitmu2Phi->SetMarkerColor(kBlack);
                hfourMuFitmu2Phi->SetLineColor(kBlack);
                hfourMuFitmu2Phi->SetMinimum(0);
                hfourMuFitmu2Phi->SetMaximum(1.7*maximum_range);
                hfourMuFitmu2Phi->Draw("hist p");
                h_mix_fourMuFitmu2Phi->SetMarkerStyle(22);
                h_mix_fourMuFitmu2Phi->SetMarkerColor(kBlue);
                h_mix_fourMuFitmu2Phi->Draw("hist p same");
                leg1->Draw("same");
                if (Ratio)
                {
                 c7->SetBottomMargin(0.3);
                 hfourMuFitmu2Phi->GetXaxis()->SetTitleSize(0);
                 hfourMuFitmu2Phi->GetXaxis()->SetLabelSize(0);
                 TH1F *Totratio = (TH1F*)h_mix_fourMuFitmu2Phi->Clone("Totratio");
                 Totratio->Divide(hfourMuFitmu2Phi);
                 Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                 Totratio->GetXaxis()->SetNoExponent(kTRUE);
                 TPad *pad5 = new TPad("pad5", "pad5", 0.0, 0.0, 0.93, 0.93);
                 pad5->SetTopMargin(0.7);
                 pad5->SetRightMargin(0.03);
                 pad5->SetFillColor(0);
                 pad5->SetGridy(1);
                 pad5->SetFillStyle(0);
                 pad5->Draw();
                 pad5->cd(0);
                 Totratio->GetXaxis()->SetTitle("#eta(#mu)");
                 Totratio->GetYaxis()->SetTitle("Mixed(#mu)/#mu");
                 Totratio->GetYaxis()->SetTitleOffset(1.2);
                 Totratio->GetXaxis()->SetTitleOffset(1.2);
                 Totratio->GetYaxis()->SetLabelSize(0.03);
                 Totratio->SetMarkerStyle(20);
                 Totratio->SetMarkerSize(1.2);
                 Totratio->SetLineColor(1);
                 Totratio->SetMarkerColor(1);
                 Totratio->GetYaxis()->SetNdivisions(204);
                 Totratio->SetMinimum(0.0);
                 Totratio->SetMaximum(2.0);
                 Totratio->Draw("ep");
                 c7->Modified();
                 c7->Update();
                  }
                c7->SaveAs("plots0p3/Muon_2_Phi."+plot_format);
                c7->Clear();
                hfourMuFitmu3Pt->GetXaxis()->SetTitle("Muon p_{T} [GeV]");
                hfourMuFitmu3Pt->Scale(1/hfourMuFitmu3Pt->Integral());
                h_mix_fourMuFitmu3Pt->Scale(1/h_mix_fourMuFitmu3Pt->Integral());
                maximum_range = max(hfourMuFitmu3Pt->GetMaximum(),h_mix_fourMuFitmu3Pt->GetMaximum());
                hfourMuFitmu3Pt->SetMarkerStyle(20);
                hfourMuFitmu3Pt->SetMarkerColor(kBlack);
                hfourMuFitmu3Pt->SetLineColor(kBlack);
                hfourMuFitmu3Pt->SetMinimum(0);
                hfourMuFitmu3Pt->SetMaximum(1.3*maximum_range);
                hfourMuFitmu3Pt->Draw("hist p");
                h_mix_fourMuFitmu3Pt->SetMarkerStyle(22);
                h_mix_fourMuFitmu3Pt->SetMarkerColor(kBlue);
                h_mix_fourMuFitmu3Pt->Draw("hist p same");
                leg1->Draw("same");
                if (Ratio)
                {
                 c7->SetBottomMargin(0.3);
                 hfourMuFitmu3Pt->GetXaxis()->SetTitleSize(0);
                 hfourMuFitmu3Pt->GetXaxis()->SetLabelSize(0);
                 TH1F *Totratio = (TH1F*)h_mix_fourMuFitmu3Pt->Clone("Totratio");
                 Totratio->Divide(hfourMuFitmu3Pt);
                 Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                 Totratio->GetXaxis()->SetNoExponent(kTRUE);
                 TPad *pad6 = new TPad("pad6", "pad6", 0.0, 0.0, 0.93, 0.93);
                 pad6->SetTopMargin(0.7);
                 pad6->SetRightMargin(0.03);
                 pad6->SetFillColor(0);
                 pad6->SetGridy(1);
                 pad6->SetFillStyle(0);
                 pad6->Draw();
                 pad6->cd(0);
                 Totratio->GetXaxis()->SetTitle("p_{T}(#mu) [GeV]");
                 Totratio->GetYaxis()->SetTitle("Mixed(#mu)/#mu");
                 Totratio->GetYaxis()->SetTitleOffset(1.2);
                 Totratio->GetXaxis()->SetTitleOffset(1.2);
                 Totratio->GetYaxis()->SetLabelSize(0.03);
                 Totratio->SetMarkerStyle(20);
                 Totratio->SetMarkerSize(1.2);
                 Totratio->SetLineColor(1);
                 Totratio->SetMarkerColor(1);
                 Totratio->GetYaxis()->SetNdivisions(204);
                 Totratio->SetMinimum(0.0);
                 Totratio->SetMaximum(2.0);
                 Totratio->Draw("ep");
                 c7->Modified();
                 c7->Update();
                  }
                c7->SaveAs("plots0p3/Muon_3_Pt."+plot_format);
                c7->Clear();

                hfourMuFitmu3Eta->GetXaxis()->SetTitle("Muon #eta");
                hfourMuFitmu3Eta->Scale(1/hfourMuFitmu3Eta->Integral());
                h_mix_fourMuFitmu3Eta->Scale(1/h_mix_fourMuFitmu3Eta->Integral());
                maximum_range = max(hfourMuFitmu3Eta->GetMaximum(),h_mix_fourMuFitmu3Eta->GetMaximum());
                hfourMuFitmu3Eta->SetMarkerStyle(20);
                hfourMuFitmu3Eta->SetMarkerColor(kBlack);
                hfourMuFitmu3Eta->SetLineColor(kBlack);
                hfourMuFitmu3Eta->SetMinimum(0);
                hfourMuFitmu3Eta->SetMaximum(1.7*maximum_range);
                hfourMuFitmu3Eta->Draw("hist p");
                h_mix_fourMuFitmu3Eta->SetMarkerStyle(22);
                h_mix_fourMuFitmu3Eta->SetMarkerColor(kBlue);
                 h_mix_fourMuFitmu3Eta->Draw("hist p same");
                leg1->Draw("same");
                if (Ratio)
                {
                 c7->SetBottomMargin(0.3);
                 hfourMuFitmu3Eta->GetXaxis()->SetTitleSize(0);
                 hfourMuFitmu3Eta->GetXaxis()->SetLabelSize(0);
                 TH1F *Totratio = (TH1F*)h_mix_fourMuFitmu3Eta->Clone("Totratio");
                 Totratio->Divide(hfourMuFitmu3Eta);
                 Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                 Totratio->GetXaxis()->SetNoExponent(kTRUE);
                 TPad *pad7 = new TPad("pad7", "pad7", 0.0, 0.0, 0.93, 0.93);
                 pad7->SetTopMargin(0.7);
                 pad7->SetRightMargin(0.03);
                 pad7->SetFillColor(0);
                 pad7->SetGridy(1);
                 pad7->SetFillStyle(0);
                 pad7->Draw();
                 pad7->cd(0);
                 Totratio->GetXaxis()->SetTitle("#eta(#mu)");
                 Totratio->GetYaxis()->SetTitle("Mixed(#mu)/#mu");
                 Totratio->GetYaxis()->SetTitleOffset(1.2);
                 Totratio->GetXaxis()->SetTitleOffset(1.2);
                 Totratio->GetYaxis()->SetLabelSize(0.03);
                 Totratio->SetMarkerStyle(20);
                 Totratio->SetMarkerSize(1.2);
                 Totratio->SetLineColor(1);
                 Totratio->SetMarkerColor(1);
                 Totratio->GetYaxis()->SetNdivisions(204);
                 Totratio->SetMinimum(0.0);
                 Totratio->SetMaximum(2.0);
                 Totratio->Draw("ep");
                 c7->Modified();
                 c7->Update();
                  }

                c7->SaveAs("plots0p3/Muon_3_Eta."+plot_format);
                c7->Clear();
                
                hfourMuFitmu3Phi->GetXaxis()->SetTitle("Muon #phi");
                hfourMuFitmu3Phi->Scale(1/hfourMuFitmu3Phi->Integral());
                h_mix_fourMuFitmu3Phi->Scale(1/h_mix_fourMuFitmu3Phi->Integral());
                maximum_range = max(hfourMuFitmu3Phi->GetMaximum(),h_mix_fourMuFitmu3Phi->GetMaximum());
                hfourMuFitmu3Phi->SetMarkerStyle(20);
                hfourMuFitmu3Phi->SetMarkerColor(kBlack);
                hfourMuFitmu3Phi->SetLineColor(kBlack);
                hfourMuFitmu3Phi->SetMinimum(0);
                hfourMuFitmu3Phi->SetMaximum(1.7*maximum_range);
                hfourMuFitmu3Phi->Draw("hist p");
                h_mix_fourMuFitmu3Phi->SetMarkerStyle(22);
                h_mix_fourMuFitmu3Phi->SetMarkerColor(kBlue);
                h_mix_fourMuFitmu3Phi->Draw("hist p same");
                leg1->Draw("same");
                if (Ratio)
                {
                 c7->SetBottomMargin(0.3);
                 hfourMuFitmu3Phi->GetXaxis()->SetTitleSize(0);
                 hfourMuFitmu3Phi->GetXaxis()->SetLabelSize(0);
                 TH1F *Totratio = (TH1F*)h_mix_fourMuFitmu3Phi->Clone("Totratio");
                 Totratio->Divide(hfourMuFitmu3Phi);
                 Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                 Totratio->GetXaxis()->SetNoExponent(kTRUE);
                 TPad *pad8 = new TPad("pad8", "pad8", 0.0, 0.0, 0.93, 0.93);
                 pad8->SetTopMargin(0.7);
                 pad8->SetRightMargin(0.03);
                 pad8->SetFillColor(0);
                 pad8->SetGridy(1);
                 pad8->SetFillStyle(0);
                 pad8->Draw();
                 pad8->cd(0);
                 Totratio->GetXaxis()->SetTitle("#phi(#mu)");
                 Totratio->GetYaxis()->SetTitle("Mixed(#mu)/#mu");
                 Totratio->GetYaxis()->SetTitleOffset(1.2);
                 Totratio->GetXaxis()->SetTitleOffset(1.2);
                 Totratio->GetYaxis()->SetLabelSize(0.03);
                 Totratio->SetMarkerStyle(20);
                 Totratio->SetMarkerSize(1.2);
                 Totratio->SetLineColor(1);
                 Totratio->SetMarkerColor(1);
                 Totratio->GetYaxis()->SetNdivisions(204);
                 Totratio->SetMinimum(0.0);
                 Totratio->SetMaximum(2.0);
                 Totratio->Draw("ep");
                 c7->Modified();
                 c7->Update();
                  }
                c7->SaveAs("plots0p3/Muon_3_Phi."+plot_format);
                c7->Clear();
                hfourMuFitmu4Pt->GetXaxis()->SetTitle("Muon p_{T} [GeV]");
                hfourMuFitmu4Pt->Scale(1/hfourMuFitmu4Pt->Integral());
                h_mix_fourMuFitmu4Pt->Scale(1/h_mix_fourMuFitmu4Pt->Integral());
                maximum_range = max(hfourMuFitmu4Pt->GetMaximum(),h_mix_fourMuFitmu4Pt->GetMaximum());
                hfourMuFitmu4Pt->SetMarkerStyle(20);
                hfourMuFitmu4Pt->SetMarkerColor(kBlack);
                hfourMuFitmu4Pt->SetLineColor(kBlack);
                hfourMuFitmu4Pt->SetMinimum(0);
                hfourMuFitmu4Pt->SetMaximum(1.3*maximum_range);
                hfourMuFitmu4Pt->Draw("hist p");
                h_mix_fourMuFitmu4Pt->SetMarkerStyle(22);
                h_mix_fourMuFitmu4Pt->SetMarkerColor(kBlue);
                h_mix_fourMuFitmu4Pt->Draw("hist p same");
                leg1->Draw("same");
                if (Ratio)
                {
                 c7->SetBottomMargin(0.3);
                 hfourMuFitmu4Pt->GetXaxis()->SetTitleSize(0);
                 hfourMuFitmu4Pt->GetXaxis()->SetLabelSize(0);
                 TH1F *Totratio = (TH1F*)h_mix_fourMuFitmu4Pt->Clone("Totratio");
                 Totratio->Divide(hfourMuFitmu4Pt);
                 Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                 Totratio->GetXaxis()->SetNoExponent(kTRUE);
                 TPad *pad9 = new TPad("pad9", "pad9", 0.0, 0.0, 0.93, 0.93);
                 pad9->SetTopMargin(0.7);
                 pad9->SetRightMargin(0.03);
                 pad9->SetFillColor(0);
                 pad9->SetGridy(1);
                 pad9->SetFillStyle(0);
                 pad9->Draw();
                 pad9->cd(0);
                 Totratio->GetXaxis()->SetTitle("p_{T}(#mu) [GeV]");
                 Totratio->GetYaxis()->SetTitle("Mixed(#mu)/#mu");
                 Totratio->GetYaxis()->SetTitleOffset(1.2);
                 Totratio->GetXaxis()->SetTitleOffset(1.2);
                 Totratio->GetYaxis()->SetLabelSize(0.03);
                 Totratio->SetMarkerStyle(20);
                 Totratio->SetMarkerSize(1.2);
                 Totratio->SetLineColor(1);
                 Totratio->SetMarkerColor(1);
                 Totratio->GetYaxis()->SetNdivisions(204);
                 Totratio->SetMinimum(0.0);
                 Totratio->SetMaximum(2.0);
                 Totratio->Draw("ep");
                 c7->Modified();
                 c7->Update();
                  }
                c7->SaveAs("plots0p3/Muon_4_Pt."+plot_format);
                c7->Clear();
                hfourMuFitmu4Eta->GetXaxis()->SetTitle("Muon #eta");
                hfourMuFitmu4Eta->Scale(1/hfourMuFitmu4Eta->Integral());
                h_mix_fourMuFitmu4Eta->Scale(1/h_mix_fourMuFitmu4Eta->Integral());
                maximum_range = max(hfourMuFitmu4Eta->GetMaximum(),h_mix_fourMuFitmu4Eta->GetMaximum());
                hfourMuFitmu4Eta->SetMarkerStyle(20);
                hfourMuFitmu4Eta->SetMarkerColor(kBlack);
                hfourMuFitmu4Eta->SetLineColor(kBlack);
                hfourMuFitmu4Eta->SetMinimum(0);
                hfourMuFitmu4Eta->SetMaximum(1.7*maximum_range);
                hfourMuFitmu4Eta->Draw("hist p");
                h_mix_fourMuFitmu4Eta->SetMarkerStyle(22);
                h_mix_fourMuFitmu4Eta->SetMarkerColor(kBlue);
                h_mix_fourMuFitmu4Eta->Draw("hist p same");
                leg1->Draw("same");
                if (Ratio)
                {
                 c7->SetBottomMargin(0.3);
                 hfourMuFitmu4Eta->GetXaxis()->SetTitleSize(0);
                 hfourMuFitmu4Eta->GetXaxis()->SetLabelSize(0);
                 TH1F *Totratio = (TH1F*)h_mix_fourMuFitmu4Eta->Clone("Totratio");
                 Totratio->Divide(hfourMuFitmu4Eta);
                 Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                 Totratio->GetXaxis()->SetNoExponent(kTRUE);
                 TPad *pad10 = new TPad("pad10", "pad10", 0.0, 0.0, 0.93, 0.93);
                 pad10->SetTopMargin(0.7);
                 pad10->SetRightMargin(0.03);
                 pad10->SetFillColor(0);
                 pad10->SetGridy(1);
                 pad10->SetFillStyle(0);
                 pad10->Draw();
                 pad10->cd(0);
                 Totratio->GetXaxis()->SetTitle("#eta(#mu)");
                 Totratio->GetYaxis()->SetTitle("Mixed(#mu)/#mu");
                 Totratio->GetYaxis()->SetTitleOffset(1.2);
                 Totratio->GetXaxis()->SetTitleOffset(1.2);
                 Totratio->GetYaxis()->SetLabelSize(0.03);
                 Totratio->SetMarkerStyle(20);
                 Totratio->SetMarkerSize(1.2);
                 Totratio->SetLineColor(1);
                 Totratio->SetMarkerColor(1);
                 Totratio->GetYaxis()->SetNdivisions(204);
                 Totratio->SetMinimum(0.0);
                 Totratio->SetMaximum(2.0);
                 Totratio->Draw("ep");
                 c7->Modified();
                 c7->Update();
                  }

                c7->SaveAs("plots0p3/Muon_4_Eta."+plot_format);
                c7->Clear();
                
                hfourMuFitmu4Phi->GetXaxis()->SetTitle("Muon #phi");
                hfourMuFitmu4Phi->Scale(1/hfourMuFitmu4Phi->Integral());
                h_mix_fourMuFitmu4Phi->Scale(1/h_mix_fourMuFitmu4Phi->Integral());
                maximum_range = max(hfourMuFitmu4Phi->GetMaximum(),h_mix_fourMuFitmu4Phi->GetMaximum());
                hfourMuFitmu4Phi->SetMarkerStyle(20);
                hfourMuFitmu4Phi->SetMarkerColor(kBlack);
                hfourMuFitmu4Phi->SetLineColor(kBlack);
                hfourMuFitmu4Phi->SetMinimum(0);
                hfourMuFitmu4Phi->SetMaximum(1.7*maximum_range);
                hfourMuFitmu4Phi->Draw("hist p");
                h_mix_fourMuFitmu4Phi->SetMarkerStyle(22);
                h_mix_fourMuFitmu4Phi->SetMarkerColor(kBlue);
                h_mix_fourMuFitmu4Phi->Draw("hist p same");
                leg1->Draw("same");
                if (Ratio)
                {
                 c7->SetBottomMargin(0.3);
                 hfourMuFitmu4Phi->GetXaxis()->SetTitleSize(0);
                 hfourMuFitmu4Phi->GetXaxis()->SetLabelSize(0);
                 TH1F *Totratio = (TH1F*)h_mix_fourMuFitmu4Phi->Clone("Totratio");
                 Totratio->Divide(hfourMuFitmu4Phi);
                 Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                 Totratio->GetXaxis()->SetNoExponent(kTRUE);
                 TPad *pad11 = new TPad("pad11", "pad11", 0.0, 0.0, 0.93, 0.93);
                 pad11->SetTopMargin(0.7);
                 pad11->SetRightMargin(0.03);
                 pad11->SetFillColor(0);
                 pad11->SetGridy(1);
                 pad11->SetFillStyle(0);
                 pad11->Draw();
                 pad11->cd(0);
                 Totratio->GetXaxis()->SetTitle("#phi(#mu)");
                 Totratio->GetYaxis()->SetTitle("Mixed(#mu)/#mu");
                 Totratio->GetYaxis()->SetTitleOffset(1.2);
                 Totratio->GetXaxis()->SetTitleOffset(1.2);
                 Totratio->GetYaxis()->SetLabelSize(0.03);
                 Totratio->SetMarkerStyle(20);
                 Totratio->SetMarkerSize(1.2);
                 Totratio->SetLineColor(1);
                 Totratio->SetMarkerColor(1);
                 Totratio->GetYaxis()->SetNdivisions(204);
                 Totratio->SetMinimum(0.0);
                 Totratio->SetMaximum(2.0);
                 Totratio->Draw("ep");
                 c7->Modified();
                 c7->Update();
                  }
                c7->SaveAs("plots0p3/Muon_4_Phi."+plot_format);
                c7->Clear();



		Ypt->GetXaxis()->SetTitle("Upsilon p_{T} [GeV]"); 
		Ypt->Scale(1/Ypt->Integral());
		mix_Ypt->Scale(1/mix_Ypt->Integral());
		maximum_range = max(Ypt->GetMaximum(),mix_Ypt->GetMaximum());
		Ypt->SetMinimum(0);
		Ypt->SetMaximum(1.3*maximum_range);
		Ypt->SetMarkerStyle(20);
		Ypt->SetMarkerColor(kBlack);
		Ypt->SetLineColor(kBlack);
        Ypt->Draw("hist p");
		TLegend* leg3 = new TLegend(0.7,0.7,0.9,0.9);
		leg3->AddEntry(Ypt,"#Upsilon","lep");
		leg3->AddEntry(mix_Ypt,"mixed #Upsilon","lep");
		leg3->Draw("same");
		mix_Ypt->SetMarkerStyle(22);
		mix_Ypt->SetMarkerColor(kBlue);
        mix_Ypt->Draw("hist p same");
                  if (Ratio)
                             {
                		c7->SetBottomMargin(0.3);
		                Ypt->GetXaxis()->SetTitleSize(0);
                		Ypt->GetXaxis()->SetLabelSize(0);
                 		TH1F *Totratio = (TH1F*)mix_Ypt->Clone("Totratio");
          		        Totratio->Divide(Ypt);
       			        Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
               			Totratio->GetXaxis()->SetNoExponent(kTRUE);
                                TPad *pad12 = new TPad("pad12", "pad12", 0.0, 0.0, 0.93, 0.93);
                                pad12->SetTopMargin(0.7);
                                pad12->SetRightMargin(0.03);
                                pad12->SetFillColor(0);
                                pad12->SetGridy(1);
	                        pad12->SetFillStyle(0);
               	                pad12->Draw();
                 	        pad12->cd(0);
              		        Totratio->GetXaxis()->SetTitle("p_{T}(#Upsilon) [GeV]");
                 		Totratio->GetYaxis()->SetTitle("Mixed(#Upsilon)/#Upsilon");
                 		Totratio->GetYaxis()->SetTitleOffset(1.2);
                 		Totratio->GetXaxis()->SetTitleOffset(1.2);
                 		Totratio->GetYaxis()->SetLabelSize(0.03);
                 		Totratio->SetMarkerStyle(20);
                 		Totratio->SetMarkerSize(1.2);
                 		Totratio->SetLineColor(1);
                 		Totratio->SetMarkerColor(1);
                 		Totratio->GetYaxis()->SetNdivisions(204);
                 		Totratio->SetMinimum(0.0);
                 		Totratio->SetMaximum(2.0);
                 		Totratio->Draw("ep");
                 		//c7->Modified();
                 		//c7->Update();
                  		}
		c7->SaveAs("plots0p3/Upsilon_Pt."+plot_format);
                c7->Clear();

                Yeta->GetXaxis()->SetTitle("#Upsilon #eta ");
                Yeta->Scale(1/Yeta->Integral());
                mix_Yeta->Scale(1/mix_Yeta->Integral());
                maximum_range = max(Yeta->GetMaximum(),mix_Yeta->GetMaximum());
                Yeta->SetMinimum(0);
                Yeta->SetMaximum(1.7*maximum_range);
                Yeta->SetMarkerStyle(20);
                Yeta->SetMarkerColor(kBlack);
                Yeta->SetLineColor(kBlack);
                Yeta->Draw("hist p");
                leg3->Draw("same");
                mix_Yeta->SetMarkerStyle(22);
                mix_Yeta->SetMarkerColor(kBlue);
             
             mix_Yeta->Draw("hist p same");
                  if (Ratio)
                              {
                                c7->SetBottomMargin(0.3);
                                Yeta->GetXaxis()->SetTitleSize(0);
                                Yeta->GetXaxis()->SetLabelSize(0);
                                TH1F *Totratio = (TH1F*)mix_Yeta->Clone("Totratio");
                                Totratio->Divide(Yeta);
                                Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                                Totratio->GetXaxis()->SetNoExponent(kTRUE);
                                TPad *pad13 = new TPad("pad13", "pad13", 0.0, 0.0, 0.93, 0.93);
                                pad13->SetTopMargin(0.7);
                                pad13->SetRightMargin(0.03);
                                pad13->SetFillColor(0);
                                pad13->SetGridy(1);
                                pad13->SetFillStyle(0);
                                pad13->Draw();
                                pad13->cd(0);
                                Totratio->GetXaxis()->SetTitle("#eta (#Upsilon)");
                                Totratio->GetYaxis()->SetTitle("Mixed(#Upsilon)/#Upsilon");
                                Totratio->GetYaxis()->SetTitleOffset(1.2);
                                Totratio->GetXaxis()->SetTitleOffset(1.2);
                                Totratio->GetYaxis()->SetLabelSize(0.03);
                                Totratio->SetMarkerStyle(20);
                                Totratio->SetMarkerSize(1.2);
                                Totratio->SetLineColor(1);
                                Totratio->SetMarkerColor(1);
                                Totratio->GetYaxis()->SetNdivisions(204);
                                Totratio->SetMinimum(0.0);
                                Totratio->SetMaximum(2.0);
                                Totratio->Draw("ep");
                                c7->Modified();
                                c7->Update();
                                }
                c7->SaveAs("plots0p3/Upsilon_Eta."+plot_format);
                c7->Clear();
                Yphi->GetXaxis()->SetTitle("Upsilon #phi ");
                Yphi->Scale(1/Yphi->Integral());
                mix_Yphi->Scale(1/mix_Yphi->Integral());
                maximum_range = max(Yphi->GetMaximum(),mix_Yphi->GetMaximum());
                Yphi->SetMinimum(0);
                Yphi->SetMaximum(1.7*maximum_range);
                Yphi->SetMarkerStyle(20);
                Yphi->SetMarkerColor(kBlack);
                Yphi->SetLineColor(kBlack);
                Yphi->Draw("hist p");
                leg3->Draw("same");
                mix_Yphi->SetMarkerStyle(22);
                mix_Yphi->SetMarkerColor(kBlue);
                mix_Yphi->Draw("hist p same");
                 if (Ratio)
                              {
                                c7->SetBottomMargin(0.3);
                                Yphi->GetXaxis()->SetTitleSize(0);
                                Yphi->GetXaxis()->SetLabelSize(0);
                                TH1F *Totratio = (TH1F*)mix_Yphi->Clone("Totratio");
                                Totratio->Divide(Yphi);
                                Totratio->GetXaxis()->SetMoreLogLabels(kTRUE);
                                Totratio->GetXaxis()->SetNoExponent(kTRUE);
                                TPad *pad14 = new TPad("pad14", "pad14", 0.0, 0.0, 0.93, 0.93);
                                pad14->SetTopMargin(0.7);
                                pad14->SetRightMargin(0.03);
                                pad14->SetFillColor(0);
                                pad14->SetGridy(1);
                                pad14->SetFillStyle(0);
                                pad14->Draw();
                                pad14->cd(0);
                                Totratio->GetXaxis()->SetTitle("#phi (#Upsilon)");
                                Totratio->GetYaxis()->SetTitle("Mixed(#Upsilon)/#Upsilon");
                                Totratio->GetYaxis()->SetTitleOffset(1.2);
                                Totratio->GetXaxis()->SetTitleOffset(1.2);
                                Totratio->GetYaxis()->SetLabelSize(0.03);
                                Totratio->SetMarkerStyle(20);
                                Totratio->SetMarkerSize(1.2);
                                Totratio->SetLineColor(1);
                                Totratio->SetMarkerColor(1);
                                Totratio->GetYaxis()->SetNdivisions(204);
                                Totratio->SetMinimum(0.0);
                                Totratio->SetMaximum(2.0);
                                Totratio->Draw("ep");
                                c7->Modified();
                                c7->Update();
                                }
                c7->SaveAs("plots0p3/Upsilon_Phi."+plot_format);
                c7->Clear();
		
		TCanvas *c8 = new TCanvas("c8","c8");
        h_mu_trg_dR->GetXaxis()->SetTitle("dR(offline muon,trigger obj)");
		h_mu_trg_dR->GetXaxis()->SetRangeUser(0.00001,0.3);
		h_mu_trg_dR->SetMarkerColor(kRed);
		h_mu_trg_dR->Draw();
		c8->SetLogy();
		//c8->SetLogx();
		c8->SaveAs("plots0p3/Trig_Match_dR."+plot_format); 

		TCanvas *c9 = new TCanvas("c9","c9");
		nVertices->Sumw2();
		nVertices->Draw();
		c9->SaveAs("plots0p3/nVertices."+plot_format);
		nVertices_2012->SetMarkerColor(kRed);
		nVertices_2012->Scale((float)nVertices->Integral()/nVertices_2012->Integral());
		nVertices_2012->GetXaxis()->SetRangeUser(0,70);
		nVertices_2012->Draw();
		nVertices->Draw("same");
		c9->SaveAs("plots0p3/nVerticesAnd2012."+plot_format);
		nVertices_fourmuon_aftercut->Draw("e1");
		c9->SaveAs("plots0p3/nVertices_aftercuts."+plot_format);
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
		c9->SaveAs("plots0p3/4muonCandsPerEvt_compare."+plot_format);
		nVertices_fourmuon_aftercut->GetYaxis()->SetLabelSize(0.02);
		nVertices_fourmuon_aftercut->GetYaxis()->SetRangeUser(-0.002,0.04);
		nVertices_fourmuon_aftercut->GetXaxis()->SetRangeUser(0,90);
		nVertices_fourmuon_aftercut->Draw("e1");
		c9->SaveAs("plots0p3/4muonCandsPerEvt_aftercut."+plot_format);


		TCanvas *c10 = new TCanvas("c10","c10");
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
		//c10->SaveAs("plots0p3/EffMuonPtCut."+plot_format);

		TGraph *gSignificanceMuonPtCut = new TGraph(5, muonPtCut, significanceMuonPtCut);
		gSignificanceMuonPtCut->Draw("APL");
		gSignificanceMuonPtCut->GetXaxis()->SetTitle("muon p_{T} [GeV]");
		gSignificanceMuonPtCut->GetYaxis()->SetTitle("S/sqrt(B)");
		//c10->SaveAs("plots0p3/significanceMuonPtCut."+plot_format);
			
		TCanvas *c11 = new TCanvas("c11","c11");
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
		//c11->SaveAs("plots0p3/EffVProbCut."+plot_format);

		TGraph *gSignificanceVProbCut = new TGraph(14, vProbCut, significanceVProbCut);
		gSignificanceVProbCut->Draw("APL");
		gSignificanceVProbCut->GetXaxis()->SetTitle("4 muons vertex fit probability [GeV]");
		gSignificanceVProbCut->GetYaxis()->SetTitle("S/sqrt(B)");
		gSignificanceVProbCut->GetXaxis()->SetRangeUser(0.001,0.1);
		//c11->SaveAs("plots0p3/significanceVProbCut."+plot_format);

		TCanvas *c12 = new TCanvas("c12","c12");
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
		//c12->SaveAs("plots0p3/EffQuaCut."+plot_format);

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
		//c12->SaveAs("plots0p3/significanceQuaCut."+plot_format);

	}
