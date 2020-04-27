/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SetHistograms.h
 * Author: Candan
 *
 * Created on November 20, 2019, 11:02 PM
 */

#ifndef SETHISTOGRAMS_H
#define SETHISTOGRAMS_H
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
 
        ///*******************************//
        /// Define Histograms
        ///*******************************//
        TFile *fout = new TFile("Plots.root","recreate");  
        TH1F *mix_Ypt        = new TH1F("mix_Ypt","",20,0,40);
        TH1F *mix_Yeta       = new TH1F("mix_Yeta","",20,-5.0,5.0);
        TH1F *mix_Yrapid     = new TH1F("mix_Yrapid","",20,-3.0,3.0);
        TH1F *mix_Yphi       = new TH1F("mix_Yphi","",20,-5,5);
        TH2F *mix_Y2Dpteta   = new TH2F("mix_Y2Dpteta","mix_Y2Dpteta",20,-5,5,20,0,40);
        TH2F *mix_Y2Dptrapid = new TH2F("mix_Y2Dptrapid","mix_Y2Dptrapid",20,-5,5,20,0,40);
        TH2F *mix_Y2Detarapid = new TH2F("mix_Y2Detarapid","mix_Y2Detarapid",20,-5,5,20,-5,5);
               ///
        TH1F *Ypt    = new TH1F("Ypt","",20,0,40);        
        TH1F *Yeta   = new TH1F("Yeta","",20,-5.0,5.0);
        TH1F *Yrapid = new TH1F("Yrapid","",20,-3.0,3.0);
        TH1F *Yphi   = new TH1F("Yphi","",20,-5.0,5.0);
        TH2F *Y2Dpteta = new TH2F("Y2Dpteta","Y2Dpteta",20,-5,5,20,0,40);
        TH2F *Y2Dptrapid = new TH2F("Y2Dptrapid","Y2Dptrapid",20,-5,5,20,0,40);
        TH2F *Y2Detarapid = new TH2F("Y2Detarapid","Y2Detarapid",20,-5,5,20,-5,5);
        TH2F *Y2Dptphi = new TH2F("Y2Dptphi","Y2Dptphi",20,-4,4,20,0,40);
        TH2F *Y2Detaphi = new TH2F("Y2Detaphi","Y2Detaphi",20,-5,5,20,-4,4);
        

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
        TH1F *hfourMuFitmu1Rap = new TH1F("hfourMuFitmu1Rap","",40,-2.5,2.5);
        TH1F *hfourMuFitmu2Rap = new TH1F("hfourMuFitmu2Rap","",40,-2.5,2.5);
       TH1F *hfourMuFitmu3Rap = new TH1F("hfourMuFitmu3Rap","",40,-2.5,2.5);
        TH1F *hfourMuFitmu4Rap = new TH1F("hfourMuFitmu4Rap","",40,-2.5,2.5);
                       
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


        TH1F *h_mix_fourMuFitmu1Rap = new TH1F("h_mix_fourMuFitmu1Rap","",40,-2.5,2.5);
        TH1F *h_mix_fourMuFitmu2Rap = new TH1F("h_mix_fourMuFitmu2Rap","",40,-2.5,2.5);
        TH1F *h_mix_fourMuFitmu3Rap = new TH1F("h_mix_fourMuFitmu3Rap","",40,-2.5,2.5);
        TH1F *h_mix_fourMuFitmu4Rap = new TH1F("h_mix_fourMuFitmu4Rap","",40,-2.5,2.5);


        //fourMuMass histograms
        TH1F *hfourMuMass = new TH1F("hfourMuMass","hfourMuMass",100,0,100);
        TH1F *hfourMuMass_aftercut = new TH1F("hfourMuMass_aftercut","hfourMuMass_aftercut",50,0,100);
        TH1F *hfourMuMass_aftercut_smallrange = new TH1F("hfourMuMass_aftercut_smallrange","hfourMuMass_aftercut_smallrange",65,13,26);
        TH1F *hfourMuMass_mix = new TH1F("hfourMuMass_mix","hfourMuMass_mix",100,0,100);
        TH1F *hfourMuMass_mix_aftercut = new TH1F("hfourMuMass_mix_aftercut","hfourMuMass_mix_aftercut",50,0,100);
        TH1F *hfourMuMass_mix_smallrange = new TH1F("hfourMuMass_mix_smallrange","hfourMuMass_mix_smallrange",30,17,20);
        
        TH1F *hfourMuMass_physkbg_mix = new TH1F("hfourMuMass_physkbg_mix","hfourMuMass_physkbg_mix",50,0,100);
	    TH1F *hfourMuMass_physkbg_mix_smallrange = new TH1F("hfourMuMass_physkbg_mix_smallrange","hfourMuMass_physkbg_mix_smallrange",65,13,26);
        TH2F *Ymumu2D_physkbg_mix = new TH2F("Ymumu2D_physkbg_mix","Ymumu2D_physkbg_mix",60,-6,6,60,-6,6);
	

        ///*******************************//
	// Clear the histos
        ///*******************************//  
        void clearHistograms(){
        mix_Ypt->Clear();
        mix_Yeta->Clear()       ;
        mix_Yrapid->Clear()     ;
        mix_Yphi->Clear()      ;
        mix_Y2Dpteta->Clear()   ;
        mix_Y2Dptrapid->Clear();
        
       ///
        Ypt->Clear();        
        Yeta->Clear();
        Yrapid->Clear();
        Yphi->Clear();
        Y2Dpteta->Clear()   ;
        Y2Dptrapid->Clear();
       //// 
       hfourMuFitmu1Pt->Clear(); 
       hfourMuFitmu2Pt->Clear();
       hfourMuFitmu3Pt->Clear();
       hfourMuFitmu4Pt->Clear();
       hfourMuFitmu1Eta->Clear();
       hfourMuFitmu2Eta->Clear();
       hfourMuFitmu3Eta->Clear();
       hfourMuFitmu4Eta->Clear();
       hfourMuFitmu1Phi->Clear();
       hfourMuFitmu2Phi->Clear();
       hfourMuFitmu3Phi->Clear();
       hfourMuFitmu4Phi-> Clear();
       hfourMuFitmu1Rap->Clear();
       hfourMuFitmu2Rap->Clear();
       hfourMuFitmu3Rap->Clear();
       hfourMuFitmu4Rap-> Clear();
       h_mix_fourMuFitmu1Pt-> Clear();
       h_mix_fourMuFitmu2Pt-> Clear();
       h_mix_fourMuFitmu3Pt-> Clear();
       h_mix_fourMuFitmu4Pt-> Clear();
       h_mix_fourMuFitmu1Eta-> Clear();
       h_mix_fourMuFitmu2Eta-> Clear();
       h_mix_fourMuFitmu3Eta-> Clear();
       h_mix_fourMuFitmu4Eta-> Clear();
       h_mix_fourMuFitmu1Phi->Clear();
       h_mix_fourMuFitmu2Phi->Clear();
       h_mix_fourMuFitmu3Phi->Clear();
       h_mix_fourMuFitmu4Phi->Clear();
       h_mix_fourMuFitmu1Rap->Clear();
       h_mix_fourMuFitmu2Rap->Clear();
       h_mix_fourMuFitmu3Rap->Clear();
       h_mix_fourMuFitmu4Rap->Clear();

        
        
        
        
        ////
        hfourMuMass->Clear(); 
        hfourMuMass_aftercut->Clear();
        hfourMuMass_aftercut_smallrange->Clear();
        hfourMuMass_mix->Clear();
        hfourMuMass_mix_aftercut->Clear();
        hfourMuMass_mix_smallrange->Clear();

        hfourMuMass_physkbg_mix->Clear();
	    hfourMuMass_physkbg_mix_smallrange->Clear() ;
        Ymumu2D_physkbg_mix->Clear() ;
        
        }

        ///*******************************//
	// Set initial values
        ///*******************************//
        void define_and_setHistograms(){
        mix_Ypt->Sumw2();
        mix_Yrapid->Sumw2();
        mix_Yeta->Sumw2();
        mix_Yphi->Sumw2();
        Ypt->Sumw2();
        Yeta->Sumw2();
        Yrapid->Sumw2();     
        Yphi->Sumw2();
        
        Ypt->GetXaxis()->SetTitle("Upsilon p_{T} [GeV]"); 
        Yeta->GetXaxis()->SetTitle("Upsilon #eta"); 
        Yrapid->GetXaxis()->SetTitle("Upsilon Rapidity"); 
        Yphi->GetXaxis()->SetTitle("Upsilon #phi"); 
        
        mix_Y2Dpteta->GetXaxis()->SetTitle("#Upsilon_mix #eta");
	    mix_Y2Dpteta->GetYaxis()->SetTitle("#Upsilon_mix p_{T} (GeV)");
	    mix_Y2Dptrapid->GetXaxis()->SetTitle("#Upsilon_mix Rapidity");
	    mix_Y2Dptrapid->GetYaxis()->SetTitle("#Upsilon_mix p_{T} (GeV)");

        hfourMuFitmu1Phi->GetXaxis()->SetTitle("Muon1 #phi");
        hfourMuFitmu2Phi->GetXaxis()->SetTitle("Muon2 #phi");
        hfourMuFitmu3Phi->GetXaxis()->SetTitle("Muon3 #phi");
        hfourMuFitmu4Phi->GetXaxis()->SetTitle("Muon4 #phi");
        hfourMuFitmu1Eta->GetXaxis()->SetTitle("Muon1 #eta");
        hfourMuFitmu2Eta->GetXaxis()->SetTitle("Muon2 #eta");
        hfourMuFitmu3Eta->GetXaxis()->SetTitle("Muon3 #eta");
        hfourMuFitmu4Eta->GetXaxis()->SetTitle("Muon4 #eta");
        hfourMuFitmu1Pt->GetXaxis()->SetTitle("Muon1 p_{T} [GeV]");
        hfourMuFitmu2Pt->GetXaxis()->SetTitle("Muon2 p_{T} [GeV]");
        hfourMuFitmu3Pt->GetXaxis()->SetTitle("Muon3 p_{T} [GeV]");
        hfourMuFitmu4Pt->GetXaxis()->SetTitle("Muon4 p_{T} [GeV]");
        hfourMuFitmu1Rap->GetXaxis()->SetTitle("Muon1 y");
        hfourMuFitmu2Rap->GetXaxis()->SetTitle("Muon2 y");
        hfourMuFitmu3Rap->GetXaxis()->SetTitle("Muon3 y");
        hfourMuFitmu4Rap->GetXaxis()->SetTitle("Muon4 y");
        
        
        
        
        
        
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
        hfourMuFitmu1Rap->Sumw2();
        hfourMuFitmu2Rap->Sumw2();
        hfourMuFitmu3Rap->Sumw2();
        hfourMuFitmu4Rap->Sumw2();

       
       
       
       
       hfourMuMass->SetStats(0);
        hfourMuMass_aftercut->SetStats(0);
        hfourMuMass_mix->SetStats(0);
        hfourMuMass_mix_smallrange->SetStats(0);
        hfourMuMass_physkbg_mix->SetStats(0);
        hfourMuMass_physkbg_mix_smallrange->SetStats(0);

        hfourMuMass->GetXaxis()->SetTitle("4 muon mass [GeV]");
        hfourMuMass->GetYaxis()->SetTitle("Candidates / GeV");
        hfourMuMass->GetYaxis()->SetLabelSize(0.03);
        hfourMuMass_aftercut->GetXaxis()->SetTitle("4 muon mass [GeV]");
        hfourMuMass_aftercut->GetYaxis()->SetTitle("Candidates / 2 GeV");
        hfourMuMass_aftercut_smallrange->GetXaxis()->SetTitle("4 muon mass [GeV]");
        hfourMuMass_aftercut_smallrange->GetYaxis()->SetTitle("Candidates / 0.2 GeV");
        hfourMuMass_mix->GetXaxis()->SetTitle("4 muon mass [GeV]");
        hfourMuMass_mix_aftercut->GetXaxis()->SetTitle("4 muon mass [GeV]");
        hfourMuMass_mix_smallrange->GetYaxis()->SetTitle("Candidates / 0.1 GeV");
        hfourMuMass_mix_smallrange->GetXaxis()->SetTitle("4 muon mass [GeV]");
        
        
        hfourMuMass_physkbg_mix->GetXaxis()->SetTitle("4 muon mass [GeV]");
	    hfourMuMass_physkbg_mix_smallrange->GetXaxis()->SetTitle("4 muon mass [GeV]");
	    hfourMuMass_physkbg_mix_smallrange->GetYaxis()->SetTitle("Candidates / 0.2 GeV");
	    Ymumu2D_physkbg_mix->GetXaxis()->SetTitle("#Upsilon #eta");
	    Ymumu2D_physkbg_mix->GetYaxis()->SetTitle("#mu_{3}#mu_{4} #eta");
        }
        
        
        //******************************************//
        // Plot Function
        //******************************************//        
        void ScaledPlot(
        
                    TString CanvasName,
                        TH1F *hist1,
                        TH1F *hist2,
                        TString Xtitle,
                        TString Ytitle,
                        bool Ratio=true
                    )
        {
                
        TCanvas *c1 = new TCanvas(CanvasName.Data(),CanvasName.Data(),800,600);
		hist1->SetStats(0);
		hist2->SetStats(0);
        hist1->SetTitle(CanvasName.Data()); 
		hist1->Scale(1/hist1->Integral());
		hist2->Scale(1/hist2->Integral());
		double maximum_range = max(hist1->GetMaximum(),hist2->GetMaximum());
		hist1->SetMinimum(0);
		hist1->SetMaximum(1.3*maximum_range);
		hist1->SetMarkerStyle(20);
		hist1->SetMarkerColor(kBlack);
		hist1->SetLineColor(kBlack);
        hist1->Draw("hist p");
		TLegend* leg3 = new TLegend(0.7,0.7,0.9,0.9);
		leg3->AddEntry(hist1,hist1->GetName(),"lep");
		leg3->AddEntry(hist2,hist2->GetName(),"lep");
		leg3->Draw("same");
		hist2->SetMarkerStyle(22);
		hist2->SetMarkerColor(kBlue);
        hist2->Draw("hist p same");
                  if (Ratio)
                             {
                		c1->SetBottomMargin(0.3);
		                hist1->GetXaxis()->SetTitleSize(0);
                		hist1->GetXaxis()->SetLabelSize(0);
                 		TH1F *Totratio = (TH1F*)hist2->Clone("Totratio");
          		        Totratio->Divide(hist1);
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
              		    Totratio->GetXaxis()->SetTitle(Xtitle); //"p_{T}(#Upsilon) [GeV]""
                 		Totratio->GetYaxis()->SetTitle(Ytitle);       //"Mixed(#Upsilon)/#Upsilon"
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
                 		c1->Modified();
                 		c1->Update();
                  		}
	    c1->SaveAs(CanvasName+".png");
        c1->Clear();
        
        }
#endif /* SETHISTOGRAMS_H */

