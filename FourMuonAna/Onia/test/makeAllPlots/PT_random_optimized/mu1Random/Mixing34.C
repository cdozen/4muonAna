/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   tester.h
 * Mentor: Zack 
 * Author: Candan
 *
 * Created on November 20, 2019, 8:38 PM
 */

#ifndef TESTER_H
#define TESTER_H

/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
///  'ADDITION AND CONVOLUTION' RooFit tutorial macro #208
///
///  One-dimensional numeric convolution
///  (require ROOT to be compiled with --enable-fftw3)
///
///  pdf = landau(t) (x) gauss(t)
///
/// \macro_image
/// \macro_output
/// \macro_code
/// \author 07/2008 - Wouter Verkerke 

#include "TROOT.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TCut.h"
#include "TLorentzVector.h"
#include "FourMu_tree.C"
#include "SetHistograms.h"
#include "RooKeysPdf.h"
#include "TF1.h"
#include "RooGenericPdf.h"
#include "TFormula.h"
using namespace RooFit;

void ProduceRandom(TLorentzVector *inputvector, // Input TLorentzVector
                   vector<TLorentzVector> &outputVector, // Output vector of TLorentzVector
                   TF1 *func, // Convoluted fit PDF
                   int i,// Event counter
                   double bin1,
                   double bin2
                   )
{



    int count=0;
    for (int k = 0; k < 10000; k++)
    {
        double val=0;
    

            val = func->GetRandom();
         //  cout<<val<<" ---- "<<i<<" \n";
        //cout<<" Produced data ["<<i<<"] :"<<pdfObs->getRealValue("Pt_upsilon")<<"\n";
        TLorentzVector temp_vector_mix;
        temp_vector_mix.SetPtEtaPhiM(val, inputvector->Eta(), inputvector->Phi(), inputvector->M());
        
        if( !(abs(temp_vector_mix.Eta())>bin1 && abs(inputvector->Eta()) < bin2)) continue;
        
        outputVector.push_back(temp_vector_mix);
        count++;
        if (count==1000) break;
    }

}

  


   



//*************************************//
// Produce Random Distribution
//*************************************//

void Mixing34()
{   
    gROOT->SetBatch();
    define_and_setHistograms();
//*************************************//
//  Eta  -- Eta bins, 
//*************************************//    
 
    bool blind_signal=true; 
//*************************************//
    vector<double> ybins = {0, 0.35, 0.65, 0.9, 1.15, 1.4, 1.65, 2.0, 2.5};

        int y_size = ybins.size();



// S e t u p   c o m p o n e n t   p d f s 

//TFile *f = new TFile("fourMuMass_tree.root");  //2018 with mu34 mas cut
    TFile *f = new TFile("/uscms_data/d3/cdozen/CMSSW_10_2_5/src/FourMuonAna/Onia/test/makeAllPlots/plots0p3/fourMuMass_tree.root");  //2018 with mu34 mas cut
    TTree *fTree = (TTree*) f->Get("FourMu_tree");

    FourMu_tree *ReadTree = new FourMu_tree(fTree);

    // Construct observable
    
    RooRealVar *var_mu3_pt = new RooRealVar("mu3_pt", "mu3_pt", 0.0, 50, "GeV");
    RooRealVar *var_mu3_eta = new RooRealVar("mu3_eta", "mu3_eta", -5.0, 5.0, "GeV");
    RooRealVar *var_mu4_pt = new RooRealVar("mu4_pt", "mu4_pt", 0.0, 50, "GeV");
    RooRealVar *var_mu4_eta = new RooRealVar("mu4_eta", "mu4_eta", -5.0, 5.0, "GeV");
    RooArgSet *observables3 = new RooArgSet("Observables");
    RooArgSet *observables4 = new RooArgSet("Observables");
    observables3->add(*var_mu3_pt);
    observables3->add(*var_mu3_eta);
    observables4->add(*var_mu4_pt);
    observables4->add(*var_mu4_eta);
    
// Construct landau(t,ml,sl) ;

    RooRealVar ml("ml", "mean landau", 5., -0.0, 40);
    RooRealVar sl("sl", "sigma landau", 1, 0.1, 20);
    RooLandau landau("lx", "lx", *var_mu3_pt, ml, sl);
    RooLandau landau4("lx4", "lx4", *var_mu4_pt, ml, sl);

    
// Construct gauss(t,mg,sg)

    RooRealVar mg("mg", "mg", 0);
    RooRealVar sg("sg", "sg", 2, 0.01, 10);
    RooGaussian gauss("gauss", "gauss", *var_mu3_pt, mg, sg);
    RooGaussian gauss4("gauss4", "gauss4", *var_mu4_pt, mg, sg);

// ---------------------------------------
// C o n s t r u c t   ERF    p d f 
// ---------------------------------------

    RooRealVar* m0shift = new RooRealVar("m0shift", "m0shift", 1, 0., 5.);
    RooRealVar* width = new RooRealVar("width", "width", 0.5, 0., 5.);
    RooRealVar* par3 = new RooRealVar("par3", "par3", 5, 0., 15.);
    RooRealVar* par4 = new RooRealVar("par4", "par4", 5, 0., 15.);


//**************************************//
// Loop on Every eta bin
//**************************************//

    TString GraphTitle3;
    TString GraphTitle4;
    for (int n = 0; n < y_size - 1; n++)
    {
        int Ntotal = fTree->GetEntries();
        TCut mycuty3; // okudugun tree yi bu cut lara gore al.
        TCut mycuty4; // okudugun tree yi bu cut lara gore al.
            mycuty3 = Form("(abs(%s) > %.2f) && (abs(%s) < %.2f)", "mu3_eta", ybins[n], "mu3_eta", ybins[n + 1]);
            mycuty4 = Form("(abs(%s) > %.2f) && (abs(%s) < %.2f)", "mu4_eta", ybins[n], "mu4_eta", ybins[n + 1]);
        
        RooDataSet *dataset3 = new RooDataSet("DataSet", "DataSet", *observables3, RooFit::Import(*fTree), RooFit::Cut(mycuty3));
        TString totalbin3;
            totalbin3 = Form("%s_%.2f_%.2f", "mu3_eta", ybins[n], ybins[n + 1]);
        RooDataSet *dataset4 = new RooDataSet("DataSet", "DataSet", *observables4, RooFit::Import(*fTree), RooFit::Cut(mycuty4));
        TString totalbin4;
            totalbin4 = Form("%s_%.2f_%.2f", "mu4_eta", ybins[n], ybins[n + 1]);
        
        
//*******************************************************************//
///PDF Models
// Activate only 1 fit model
//*******************************************************************//      
        // Fitno1-) Construct landau (x) gauss

        RooFFTConvPdf *lxg3=new RooFFTConvPdf ("lxg3","landau (X) gauss",*var_mu3_pt,landau,gauss) ;
        RooFFTConvPdf *lxg4=new RooFFTConvPdf ("lxg4","landau (X) gauss",*var_mu4_pt,landau4,gauss4) ;
        
        //*******************************************************************//      
        
        //Fitno2-) Construct Keys Pdf
        
         //RooKeysPdf *lxg= new RooKeysPdf("KDE","KDE", *var_mu3_pt, *dataset, RooKeysPdf::MirrorBoth, 1);         
        
        //*******************************************************************//      
        
        //Fitno3-) Construct ERF
        
        //RooGenericPdf* lxg = new RooGenericPdf("ERF", "LikeSignPdf", "exp(-@0/par3)*(TMath::Erf((@0-m0shift)/width)+1)", RooArgList(*var_mu3_pt, *m0shift, *width, *par3));
        
        //*******************************************************************//      
        //*******************************************************************//      



       //*******************************************************************//
        /// Fitter ...
        //*******************************************************************//
        cout << "!!******************************!" << endl;
        cout << "Fitting : " << lxg3->GetName() << "\n";
        cout << "!!******************************!" << endl;
        
        // ----------------------------------------------------------------------
        // S a m p l e ,   f i t   a n d   p l o t   c o n v o l u t e d   p d f 
        // ----------------------------------------------------------------------
        // Fit gxlx to data
            lxg3->fitTo(*dataset3, RooFit::NumCPU(4, RooFit::BulkPartition));
            TF1 *func3 ;
            func3 = lxg3->asTF(RooArgList(*var_mu3_pt));
            // Plot data, landau pdf, landau (X) gauss pdf
            GraphTitle3 = lxg3->GetName();
            GraphTitle3 += "_mu3_EtaBin";
        
            RooPlot* frame3 = var_mu3_pt->frame(Title(GraphTitle3 + totalbin3));
            dataset3->plotOn(frame3);
            lxg3->plotOn(frame3);
    
        cout << "!!******************************!" << endl;
        cout << "Fitting : " << lxg4->GetName() << "\n";
        cout << "!!******************************!" << endl;

            lxg4->fitTo(*dataset4, RooFit::NumCPU(4, RooFit::BulkPartition));
            TF1 *func4 ;
            func4 = lxg4->asTF(RooArgList(*var_mu4_pt));
            GraphTitle4 = lxg4->GetName();
            GraphTitle4 += "_mu4_EtaBin";
        
            RooPlot* frame4 = var_mu4_pt->frame(Title(GraphTitle4 + totalbin4));
            dataset4->plotOn(frame4);
            lxg4->plotOn(frame4);
//****************************//
//Dataset Generation
//****************************//

        cout << "!!***************************************!!\n";
        cout << " Readed   :" << Ntotal << "\n";
        cout << "!!***************************************!!\n";

    vector<TLorentzVector> mu3_MixingVect;
    vector<TLorentzVector> mu4_MixingVect;

        for (int i = 0; i < Ntotal; i++)
        {


            ReadTree->GetEntry(i);
            
            
            if (!((abs(ReadTree->mu3_eta) > ybins[n] && abs(ReadTree->mu3_eta) < ybins[n + 1]) || (abs(ReadTree->mu4_eta) > ybins[n] && abs(ReadTree->mu4_eta) < ybins[n + 1]))) continue;
          
            cout << "\r" << "Entry :" << i << "/" << Ntotal << flush;
          ProduceRandom(ReadTree->mu3, mu3_MixingVect, func3, i, ybins[n], ybins[n + 1]);
           //cout <<"mu3 mixing vector size "<< mu3_MixingVect.size() <<endl; 
            
            //if (!(abs(ReadTree->mu4_eta) > ybins[n] && abs(ReadTree->mu4_eta) < ybins[n + 1])) continue;
            cout << "\r" << "Entry :" << i << "/" << Ntotal << flush;
          ProduceRandom(ReadTree->mu4, mu4_MixingVect, func4, i, ybins[n], ybins[n + 1]);
           //cout <<"mu4 mixing vector size "<< mu4_MixingVect.size() <<endl; 
            


            Ypt->Fill(ReadTree->mu12->Pt());
            Yeta->Fill(ReadTree->mu12->Eta());
            Yrapid->Fill(ReadTree->mu12->Rapidity());
            Yphi->Fill(ReadTree->mu12->Phi());
            hfourMuFitmu4Pt->Fill(ReadTree->mu34->Pt());
            hfourMuFitmu4Eta->Fill(ReadTree->mu34->Eta());
            hfourMuFitmu4Rap->Fill(ReadTree->mu34->Rapidity());
            Y2Dpteta->Fill(ReadTree->mu12->Eta(), ReadTree->mu12->Pt());
            Y2Detarapid->Fill(ReadTree->mu12->Eta(), ReadTree->mu12->Rapidity());
            Y2Dptrapid->Fill(ReadTree->mu12->Rapidity(), ReadTree->mu12->Pt());
            hfourMuMass_aftercut->Fill(ReadTree->Mass4mu);
            hfourMuMass_aftercut_smallrange->Fill(ReadTree->Mass4mu);
 
 

        int Nm4 =mu4_MixingVect.size();
        int Nm3 =mu3_MixingVect.size();
        int runner =0;
        if(Nm3>Nm4 ) runner =Nm4;
        else runner =Nm3;
        
        
        for(int j= 0; j< runner; j++ )

        {
                
                TLorentzVector mu34Mixing = mu4_MixingVect[j] + mu3_MixingVect[j];


//if( abs(mu34Mixing.Eta())- abs((*ReadTree->mu34).Eta()) > 0.1) continue;
//if( abs(mu34Mixing.Rapidity())- abs((*ReadTree->mu34).Rapidity()) > 0.1) continue;

            TLorentzVector mixFourMu = mu34Mixing + *ReadTree->mu12 ;

//if (!( mu34Mixing.M()< 9.2 && 
//    (*ReadTree->mu1).Pt() >= muonPtCut[2] &&  (*ReadTree->mu2).Pt() >= muonPtCut[2]&& 
//    mu3_MixingVect[j].Pt() >= muonPtCut[2] && mu4_MixingVect[j].Pt() >= muonPtCut[2]))continue;



                cout<<"\r Run :"<< i<<"/"<<Ntotal<<flush;  
                




                mix_Ypt->Fill(mu34Mixing.Pt());
                mix_Yeta->Fill(mu34Mixing.Eta());
                mix_Yrapid->Fill(mu34Mixing.Rapidity());
                mix_Yphi->Fill(mu34Mixing.Phi());
                mix_Y2Dpteta->Fill(mu34Mixing.Eta(), mu34Mixing.Pt());
                mix_Y2Detarapid->Fill(mu34Mixing.Eta(), mu34Mixing.Rapidity());
                mix_Y2Dptrapid->Fill(mu34Mixing.Rapidity(), mu34Mixing.Pt());
                
                hfourMuMass_physkbg_mix->Fill(mixFourMu.M());
                hfourMuMass_physkbg_mix_smallrange->Fill(mixFourMu.M());
                
                //Ymumu2D_physkbg_mix->Fill(mu34Mixing.Eta(), ReadTree->mu12->Eta());
                //h_mix_fourMuFitmu4Pt->Fill(mu4_MixingVect[s].Pt());
                //h_mix_fourMuFitmu4Eta->Fill(mu4_MixingVect[s].Eta());
                //h_mix_fourMuFitmu4Rap->Fill(mu4_MixingVect[s].Rapidity());
            
                 
        }//mixing size loop

        

        mu4_MixingVect.clear();
        mu3_MixingVect.clear();
        //mu34_MixingVect.clear();
  }//event loop
    
TString MainfCanv = "fitplot_" + GraphTitle3 + totalbin3;
TCanvas *c1 = new TCanvas(MainfCanv, MainfCanv, 600, 600);
gPad->SetLeftMargin(0.15);
frame3->GetYaxis()->SetTitleOffset(1.4);
frame3->Draw();
c1->SaveAs(MainfCanv + ".png");

TString MainfCanv4 = "fitplot_" + GraphTitle4 + totalbin4;
TCanvas *c111 = new TCanvas(MainfCanv4, MainfCanv4, 600, 600);
gPad->SetLeftMargin(0.15);
frame4->GetYaxis()->SetTitleOffset(1.4);
frame4->Draw();
c111->SaveAs(MainfCanv4 + ".png");
    
} //end loop eta bin
//***********************************************************//  
// Draw frame on canvas
//***********************************************************//  
    TString Cname = "Graph_total_"+ GraphTitle3;
    TString Cname4 = "Graph_total_"+ GraphTitle4;
    
          
    ScaledPlot(Cname, Ypt, mix_Ypt, "p_{T}(#Upsilon) [GeV]", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot(Cname4, Ypt, mix_Ypt, "p_{T}(#Upsilon) [GeV]", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot("Y_Etavsmix_YEta", Yeta, mix_Yeta, "#Eta (#Upsilon)", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot("Y_Rapidityvsmix_YRapidty", Yrapid, mix_Yrapid, "Rapidity (#Upsilon) [GeV]", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot("muon_pt vs mixed_muon_pt",hfourMuFitmu4Pt, h_mix_fourMuFitmu4Pt, "p_{T}(#muon) [GeV]", "Mixed(#muon)/#muon");
    ScaledPlot("muon_eta vs mixed_muons_eta",hfourMuFitmu4Eta, h_mix_fourMuFitmu4Eta, "#Eta (#muon) [GeV]", "Mixed(#muon)/#muon");
    ScaledPlot("muon_y vs mixed_muon_y",hfourMuFitmu4Rap, h_mix_fourMuFitmu4Rap, "Rapidty (#muon) [GeV]", "Mixed(#muon)/#muon");
    
    ///draw Histograms:
    
    TCanvas *c2 = new TCanvas("c2","c2",800,600);

    Y2Dpteta->Draw("colz");
    c2->SaveAs("plots2D/Y2Dpteta.png");
    mix_Y2Dpteta->Draw("colz");
    c2->SaveAs("plots2D/Y2Dpteta_mix.png");
    Y2Dptrapid->Draw("colz");
    c2->SaveAs("plots2D/Y2Dptrapidity.png");
    Y2Detarapid->Draw("colz");
    c2->SaveAs("plots2D/Y2Detarapidity.png");
    mix_Y2Dptrapid->Draw("colz");
    c2->SaveAs("plots2D/Y2Dptrapidity_mix.png");
    mix_Y2Detarapid->Draw("colz");
    //Y2Dptphi->Draw("colz");
    //c2->SaveAs("plots0p3/Y2Dptphi."+plot_format);
    //Y2Detaphi->Draw("colz");
    //c2->SaveAs("plots0p3/Y2Detaphi."+plot_format);




    //before cuts, mix and mix zoom in
    hfourMuMass_mix->SetMarkerStyle(24);
    hfourMuMass_mix->SetMarkerColor(kRed);
    hfourMuMass_mix->Draw("e1");
    //c2->SaveAs("hfourMuMass_mix.png");
    hfourMuMass_mix_smallrange->SetMarkerStyle(24);
    hfourMuMass_mix_smallrange->SetMarkerColor(kBlue);
    hfourMuMass_mix_smallrange->GetYaxis()->SetRangeUser(0,200);
    hfourMuMass_mix_smallrange->Draw("e1");
    //c2->SaveAs("hfourMuMass_mix_smallrange.png");

    //before cuts, compare original vs mix
    hfourMuMass->Draw("e1");
    //c2->SaveAs("hfourMuMass.png");
    hfourMuMass_mix->Draw("e1same");
    c2->SaveAs("plots/hfourMuMass_origvsmix.png");
  
    //after cuts, orignal
    hfourMuMass_aftercut->Draw("e1");
    //c2->SaveAs("hfourMuMass_afterCut.png");
    float scaleEntries = hfourMuMass_aftercut_smallrange->Integral();
    hfourMuMass_aftercut_smallrange->SetStats(0);
    hfourMuMass_aftercut_smallrange->Draw("e1");
    //c2->SaveAs("hfourMuMass_afterCut_smallrange.png");
  
    //after cuts, mix
    hfourMuMass_mix_aftercut->SetMarkerStyle(24);
    hfourMuMass_mix_aftercut->SetMarkerColor(kRed);
    hfourMuMass_mix_aftercut->Draw("e1");  
    //c2->SaveAs("hfourMuMass_mix_afterCut.png");

    //after cuts, compare original vs mix
    hfourMuMass_aftercut->Draw("e1");
    hfourMuMass_mix_aftercut->Draw("e1same");
    c2->SaveAs("plots/hfourMuMass_afterCut_origvsmix.png");
 
     //after cuts, mix physics bkg an zoom in
    hfourMuMass_physkbg_mix->SetMarkerStyle(22);
    hfourMuMass_physkbg_mix->SetMarkerColor(kBlue);
    hfourMuMass_physkbg_mix->Draw("e1");
    c2->SaveAs("plots/hfourMuMass_physkbg_mix.png");
    hfourMuMass_physkbg_mix_smallrange->SetMarkerStyle(22);
    hfourMuMass_physkbg_mix_smallrange->SetMarkerColor(kBlue);
    hfourMuMass_physkbg_mix_smallrange->SetMarkerStyle(22);
    hfourMuMass_physkbg_mix_smallrange->SetMarkerColor(kBlue);
    hfourMuMass_physkbg_mix_smallrange->Draw("e1");
    c2->SaveAs("plots/hfourMuMass_physkbg_mix_smallrange.png");
    //c2->SaveAs(MixPhysBG_aftercuts_small".png");


    //after cuts, compare original vs mixPhysPkg
    TCanvas *c3 = new TCanvas("c3","c3",800,600);
    hfourMuMass_aftercut->SetMarkerStyle(20);
    hfourMuMass_aftercut->SetMarkerColor(kBlack);
    hfourMuMass_aftercut->SetLineColor(kBlack);
    hfourMuMass_aftercut->Draw("e1");
    std::cout << "hfourMuMass_aftercut->Integral() = " << hfourMuMass_aftercut->Integral() << std::endl;
    hfourMuMass_physkbg_mix->Scale(hfourMuMass_aftercut->Integral()/hfourMuMass_physkbg_mix->Integral());
    hfourMuMass_physkbg_mix->Draw("e1same");
    std::cout << "hfourMuMass_physkbg_mix->Integral() = " << hfourMuMass_physkbg_mix->Integral() << std::endl;


    TString total_data = "data (" + std::to_string((int)hfourMuMass_aftercut->Integral()) + " total)";
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(hfourMuMass_aftercut,total_data,"lep");
    leg->AddEntry(hfourMuMass_physkbg_mix,"mixed BG","lep");
    leg->Draw("same");
    c3->SaveAs("plots/hfourMuMass_afterCut_origVSphyskbg.png");
    hfourMuMass_aftercut_smallrange->SetMarkerStyle(20);
    hfourMuMass_aftercut_smallrange->SetMarkerColor(kBlack);
    TH1F *hfourMuMass_aftercut_smallrange_copy = (TH1F*)hfourMuMass_aftercut_smallrange->Clone("hfourMuMass_aftercut_smallrange_copy");
    if(blind_signal)
        { for(int i=1; i <= hfourMuMass_aftercut_smallrange_copy->GetSize()-2; i++)
            {
                if(i>= hfourMuMass_aftercut_smallrange_copy->FindBin(17.5) && i<= hfourMuMass_aftercut_smallrange_copy->FindBin(19.5))
                     hfourMuMass_aftercut_smallrange_copy->SetBinContent(i,0);
   
            }    
         }
    hfourMuMass_aftercut_smallrange_copy->Draw("e1");
    hfourMuMass_physkbg_mix_smallrange->Scale(scaleEntries/hfourMuMass_physkbg_mix_smallrange->Integral());
    hfourMuMass_physkbg_mix_smallrange->Draw("e1same");
    leg->Draw("same");
    c3->SaveAs("plots/hfourMuMass_afterCut_origVSphyskbg_smallrange.png");
   
   
    clearHistograms();
 
fout->Write();
f->Close();
fout->Close();
    
}


#endif /* TESTER_H */



