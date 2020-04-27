/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
//mu3 and mu4 choosen randomly in the same event and then mixed.
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

void ProduceRandom(TLorentzVector *inputvector1, // Input TLorentzVector
                   TLorentzVector *inputvector2,
                   vector<TLorentzVector> &outputVector, // Output vector of TLorentzVector
                   TF1 *func1, // Convoluted fit PDF
                   TF1 *func2, // Convoluted fit PDF
                   int i,// Event counter
                   double bin1,
                   double bin2,
                   bool Set_Eta
                   )
{



    int count=0;
    for (int k = 0; k < 10000; k++)
    {
        double val1=0;
        double val2=0;
    

            val1 = func1->GetRandom();
            val2 = func2->GetRandom();
         //  cout<<val<<" ---- "<<i<<" \n";
        //cout<<" Produced data ["<<i<<"] :"<<pdfObs->getRealValue("Pt_upsilon")<<"\n";
        TLorentzVector temp_vector_mix12;
        TLorentzVector temp_vector_mix1;
        TLorentzVector temp_vector_mix2;
        temp_vector_mix1.SetPtEtaPhiM(val1, inputvector1->Eta(), inputvector1->Phi(), inputvector1->M());
        temp_vector_mix2.SetPtEtaPhiM(val2, inputvector2->Eta(), inputvector2->Phi(), inputvector2->M());
        if(Set_Eta) if( !(abs(temp_vector_mix1.Eta())>bin1 && abs(temp_vector_mix1.Eta()) < bin2)) continue;
        if(Set_Eta) if( !(abs(temp_vector_mix2.Eta())>bin1 && abs(temp_vector_mix2.Eta()) < bin2)) continue;
        if(!Set_Eta)if( !(abs(temp_vector_mix1.Rapidity())>bin1 && abs(temp_vector_mix1.Rapidity()) < bin2)) continue;
        if(!Set_Eta)if( !(abs(temp_vector_mix2.Rapidity())>bin1 && abs(temp_vector_mix2.Rapidity()) < bin2)) continue;
        
        temp_vector_mix12 = temp_vector_mix1 + temp_vector_mix2; 
        if (Set_Eta) if( !(abs(temp_vector_mix12.Eta())>bin1 && abs(temp_vector_mix12.Eta()) < bin2)) continue;
        if(!Set_Eta)if( !(abs(temp_vector_mix12.Rapidity())>bin1 && abs(temp_vector_mix12.Rapidity()) < bin2)) continue;
        
        outputVector.push_back(temp_vector_mix12);
        count++;
        if (count==1000) break;
    }

}


//*************************************//
// Produce Random Distribution
//*************************************//

void mu3_mu4ptrandom()
{   
    gROOT->SetBatch();
    define_and_setHistograms();
    //*************************************//
    //  ----- Eta(true) or Rapidity (false)bins ------ 
    //*************************************//    
    bool Set_Eta =true; 
    bool blind_signal=true; 
    //*************************************//
     vector<double> ybins;
    if(!Set_Eta)
    {
        ybins = {0, 0.5, 0.85, 1.10, 1.3, 1.5, 1.7, 1.9, 2.5};
    }else
    {
    ybins = {0, 0.35, 0.65, 0.9, 1.15, 1.4, 1.65, 2.0, 2.5};
    //ybins = {0,0.35, 0.7,0.9, 1.2,1.4, 1.6, 1.9, 2.2, 2.5, 2.9, 5.0};
    //vector<double> ybins = {0, 2.5};
    //vector<double> ybins = {0., 0.7, 1.2, 1.6, 1.9, 2.2, 2.5};
    }
    int y_size = ybins.size();


    // S e t u p   c o m p o n e n t   p d f s 

    //TFile *f = new TFile("fourMuMass_tree.root");  //2018 with mu34 mas cut
    TFile *f = new TFile("/uscms_data/d3/cdozen/CMSSW_10_2_5/src/FourMuonAna/Onia/test/makeAllPlots/plots0p3/fourMuMass_tree.root");  //2018 with mu34 mas cut
    TTree *fTree = (TTree*) f->Get("FourMu_tree");

    FourMu_tree *ReadTree = new FourMu_tree(fTree);

    // Construct observable
    
    RooRealVar *var_mu3_pt = new RooRealVar("mu3_pt", "mu3_pt", 0.0, 50, "GeV");
    RooRealVar *var_mu3_eta = new RooRealVar("mu3_eta", "mu3_eta", -5.0, 5.0, "GeV");
    RooRealVar *var_mu3_y = new RooRealVar("mu3_y", "mu3_y", -5.0, 5.0, "GeV");
    RooRealVar *var_mu4_pt = new RooRealVar("mu4_pt", "mu4_pt", 0.0, 50, "GeV");
    RooRealVar *var_mu4_eta = new RooRealVar("mu4_eta", "mu4_eta", -5.0, 5.0, "GeV");
    RooRealVar *var_mu4_y = new RooRealVar("mu4_y", "mu4_y", -5.0, 5.0, "GeV");
    RooArgSet *observables = new RooArgSet("Observables");
    observables->add(*var_mu3_pt);
    observables->add(*var_mu3_eta);
    observables->add(*var_mu4_eta);
    observables->add(*var_mu3_y);
    observables->add(*var_mu4_pt);
    observables->add(*var_mu4_y);
    
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
          if(Set_Eta){
          mycuty3 = Form("(abs(%s) > %.2f) && (abs(%s) < %.2f)", "mu3_eta", ybins[n], "mu3_eta", ybins[n + 1]);
            mycuty4 = Form("(abs(%s) > %.2f) && (abs(%s) < %.2f)", "mu4_eta", ybins[n], "mu4_eta", ybins[n + 1]);
        }else
        {
            mycuty3 = Form("(abs(%s) > %.2f) && (abs(%s) < %.2f)", "mu3_y", ybins[n], "mu3_eta", ybins[n + 1]);
            mycuty4 = Form("(abs(%s) > %.2f) && (abs(%s) < %.2f)", "mu4_y", ybins[n], "mu4_eta", ybins[n + 1]);
        }
        RooDataSet *dataset = new RooDataSet("DataSet", "DataSet", *observables, RooFit::Import(*fTree), RooFit::Cut(mycuty3+mycuty4));
        TString totalbin3,totalbin4;
        if(Set_Eta){
        totalbin3 = Form("%s_%.2f_%.2f", "mu3_eta", ybins[n], ybins[n + 1]);
        totalbin4 = Form("%s_%.2f_%.2f", "mu4_eta", ybins[n], ybins[n + 1]);
        }else
        {
        totalbin3 = Form("%s_%.2f_%.2f", "mu3_y", ybins[n], ybins[n + 1]);
        totalbin4 = Form("%s_%.2f_%.2f", "mu4_y", ybins[n], ybins[n + 1]);
        }
        //RooDataSet *dataset3 = new RooDataSet("DataSet", "DataSet", *observables, RooFit::Import(*fTree), RooFit::Cut(mycuty3));
        //TString totalbin3;
          //  totalbin3 = Form("%s_%.2f_%.2f", "mu3_eta", ybins[n], ybins[n + 1]);
        //RooDataSet *dataset4 = new RooDataSet("DataSet", "DataSet", *observables, RooFit::Import(*fTree), RooFit::Cut(mycuty4));
        //TString totalbin4;
          //  totalbin4 = Form("%s_%.2f_%.2f", "mu4_eta", ybins[n], ybins[n + 1]);
        
        
//*******************************************************************//
///PDF Models
// Activate only 1 fit model
//*******************************************************************//      
        // Fitno1-) Construct landau (x) gauss

        //RooFFTConvPdf *lxg3=new RooFFTConvPdf ("lxg3","landau (X) gauss",*var_mu3_pt,landau,gauss) ;
        //RooFFTConvPdf *lxg4=new RooFFTConvPdf ("lxg4","landau (X) gauss",*var_mu4_pt,landau4,gauss4) ;
        
        //*******************************************************************//      
        
        //Fitno2-) Construct Keys Pdf
        
         RooKeysPdf *lxg3= new RooKeysPdf("KDE","KDE", *var_mu3_pt, *dataset, RooKeysPdf::MirrorBoth, 1);         
         RooKeysPdf *lxg4= new RooKeysPdf("KDE","KDE", *var_mu4_pt, *dataset, RooKeysPdf::MirrorBoth, 1);         
        
        //*******************************************************************//      
        
        //Fitno3-) Construct ERF
        
        //RooGenericPdf* lxg3 = new RooGenericPdf("ERF", "LikeSignPdf", "exp(-@0/par3)*(TMath::Erf((@0-m0shift)/width)+1)", RooArgList(*var_mu3_pt, *m0shift, *width, *par3));
        //RooGenericPdf* lxg4 = new RooGenericPdf("ERF", "LikeSignPdf", "exp(-@0/par4)*(TMath::Erf((@0-m0shift)/width)+1)", RooArgList(*var_mu4_pt, *m0shift, *width, *par4));
        
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
            lxg3->fitTo(*dataset, RooFit::NumCPU(4, RooFit::BulkPartition));
            TF1 *func1 ;
            func1 = lxg3->asTF(RooArgList(*var_mu3_pt));
            // Plot data, landau pdf, landau (X) gauss pdf
            GraphTitle3 = lxg3->GetName();
            if(Set_Eta) {
                GraphTitle3 += "_mu3_EtaBin";
            }else
            {
                GraphTitle3 += "_mu3_rapidityBin";
                }
            
            RooPlot* frame3 = var_mu3_pt->frame(Title(GraphTitle3 + totalbin3));
            dataset->plotOn(frame3);
            lxg3->plotOn(frame3);
    
        cout << "!!******************************!" << endl;
        cout << "Fitting : " << lxg4->GetName() << "\n";
        cout << "!!******************************!" << endl;

            lxg4->fitTo(*dataset, RooFit::NumCPU(4, RooFit::BulkPartition));
            TF1 *func2 ;
            func2= lxg4->asTF(RooArgList(*var_mu4_pt));
            GraphTitle4 = lxg4->GetName();
            
            if(Set_Eta) {
                GraphTitle4 += "_mu4_EtaBin";
            }else
            {
                GraphTitle3 += "_mu4_rapidityBin";
                }
            RooPlot* frame4 = var_mu4_pt->frame(Title(GraphTitle4 + totalbin4));
            dataset->plotOn(frame4);
            lxg4->plotOn(frame4);

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
 


//****************************//
//Dataset Generation
//****************************//

        cout << "!!***************************************!!\n";
        cout << " Readed   :" << Ntotal << "\n";
        cout << "!!***************************************!!\n";

    vector<TLorentzVector> mu34_MixingVect;
    //vector<TLorentzVector> mu4_MixingVect;

        for (int i = 0; i < Ntotal; i++)
        {


            ReadTree->GetEntry(i);
            
            
          if (Set_Eta) if (!((abs(ReadTree->mu3_eta) > ybins[n] && abs(ReadTree->mu3_eta) < ybins[n + 1]) && (abs(ReadTree->mu4_eta) > ybins[n] && abs(ReadTree->mu4_eta) < ybins[n + 1])))continue;
          if (!Set_Eta) if (!((abs(ReadTree->mu3_y) > ybins[n] && abs(ReadTree->mu3_y) < ybins[n + 1]) && (abs(ReadTree->mu4_y) > ybins[n] && abs(ReadTree->mu4_y) < ybins[n + 1])))continue;
            cout << "\r" << "Entry :" << i << "/" << Ntotal << flush;
          
          ProduceRandom(ReadTree->mu3,ReadTree->mu4, mu34_MixingVect, func1,func2, i, ybins[n], ybins[n + 1],Set_Eta);
            

            Ypt->Fill(ReadTree->mu12->Pt());
            Yeta->Fill(ReadTree->mu12->Eta());
            Yrapid->Fill(ReadTree->mu12->Rapidity());
            Yphi->Fill(ReadTree->mu12->Phi());
            hfourMuFitmu4Pt->Fill(ReadTree->mu34->Pt());
            hfourMuFitmu4Eta->Fill(ReadTree->mu34->Eta());
            hfourMuFitmu4Rap->Fill(ReadTree->mu34->Rapidity());
            //Y2Dpteta->Fill(ReadTree->mu12->Eta(), ReadTree->mu12->Pt());
            //Y2Detarapid->Fill(ReadTree->mu12->Eta(), ReadTree->mu12->Rapidity());
            //Y2Dptrapid->Fill(ReadTree->mu12->Rapidity(), ReadTree->mu12->Pt());
            //Ypt->Fill(ReadTree->mu34->Pt());
            //Yeta->Fill(ReadTree->mu34->Eta());
            //Yrapid->Fill(ReadTree->mu34->Rapidity());
            //Yphi->Fill(ReadTree->mu34->Phi());
            Y2Dpteta->Fill(ReadTree->mu34->Eta(), ReadTree->mu34->Pt());
            Y2Detarapid->Fill(ReadTree->mu34->Eta(), ReadTree->mu34->Rapidity());
            Y2Dptrapid->Fill(ReadTree->mu34->Rapidity(), ReadTree->mu34->Pt());
            hfourMuMass_aftercut->Fill(ReadTree->Mass4mu);
            hfourMuMass_aftercut_smallrange->Fill(ReadTree->Mass4mu);
 
 

           
        
        int Nm34 =mu34_MixingVect.size();
        for(int j= 0; j<Nm34; j++ )

        {
                
         if(Set_Eta)  if (!(abs(mu34_MixingVect[j].Eta()) > ybins[n] && abs(mu34_MixingVect[j].Eta()) < ybins[n + 1]))continue;
         if(!Set_Eta)  if (!(abs(mu34_MixingVect[j].Rapidity()) > ybins[n] && abs(mu34_MixingVect[j].Rapidity()) < ybins[n + 1]))continue;
        
        //Double_t DeltaPhi=(fabs(mu34_MixingVect[j].Phi())- abs((*ReadTree->mu34).Phi()));
        //Double_t DeltaY=(fabs(mu34_MixingVect[j].Rapidity())- abs((*ReadTree->mu34).Rapidity()));
        //Double_t DeltaR=sqrt(DeltaPhi*DeltaPhi)+(DeltaY*DeltaY);
        //if(Set_Eta)if(DeltaR>0.3)continue;
        if(Set_Eta)if(! (abs(mu34_MixingVect[j].Eta())- abs((*ReadTree->mu12).Eta()) < 0.1)) continue;
        //if (Set_Eta)if( abs(mu34_MixingVect[j].Phi())- abs((*ReadTree->mu34).Phi()) >0.6) continue;
        //if( abs(mu34_MixingVect[j].Rapidity())- abs((*ReadTree->mu34).Rapidity()) > 0.05) continue;
//if( abs(mu34Mixing.Rapidity())- abs((*ReadTree->mu34).Rapidity()) > 0.1) continue;
//if (fabs( ReadTree->mu34->M()-mu34_MixingVect[j].M())>0.5*ReadTree->mu34->M() ) continue;
//if (mu34_MixingVect[j].M()> 9.2) continue;
            
            TLorentzVector mixFourMu = mu34_MixingVect[j] + *ReadTree->mu12 ;

                cout<<"\r Run :"<< i<<"/"<<Ntotal<<flush;  
                

                mix_Ypt->Fill(mu34_MixingVect[j].Pt());
                mix_Yeta->Fill(mu34_MixingVect[j].Eta());
                mix_Yrapid->Fill(mu34_MixingVect[j].Rapidity());
                mix_Yphi->Fill(mu34_MixingVect[j].Phi());
                mix_Y2Dpteta->Fill(mu34_MixingVect[j].Eta(), mu34_MixingVect[j].Pt());
                mix_Y2Detarapid->Fill(mu34_MixingVect[j].Eta(), mu34_MixingVect[j].Rapidity());
                mix_Y2Dptrapid->Fill(mu34_MixingVect[j].Rapidity(), mu34_MixingVect[j].Pt());
                
                hfourMuMass_physkbg_mix->Fill(mixFourMu.M());
                hfourMuMass_physkbg_mix_smallrange->Fill(mixFourMu.M());
                
                Ymumu2D_physkbg_mix->Fill(mu34_MixingVect[j].Eta(), ReadTree->mu12->Eta());
                h_mix_fourMuFitmu4Pt->Fill(mu34_MixingVect[j].Pt());
                h_mix_fourMuFitmu4Eta->Fill(mu34_MixingVect[j].Eta());
                h_mix_fourMuFitmu4Rap->Fill(mu34_MixingVect[j].Rapidity());
            
                 
        }//mixing size loop
        

        mu34_MixingVect.clear();
      //  mu3_MixingVect.clear();
        //mu34_MixingVect.clear();
  }//event loop
    
} //end loop eta bin
//***********************************************************//  
// Draw frame on canvas
//***********************************************************//  
    TString Cname = "Graph_total_"+ GraphTitle3;
    TString Cname4 = "Graph_total_"+ GraphTitle4;
    
          
    ScaledPlot(Cname, Ypt, mix_Ypt, "p_{T}(#Upsilon) [GeV]", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot(Cname4, Ypt, mix_Ypt, "p_{T}(#Upsilon) [GeV]", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot("Yptvsmu34mixpt", Ypt, mix_Ypt, "p_{T}(#Upsilon) [GeV]", "Mixed(#Upsilon)/#Upsilon");
    //ScaledPlot("Y_Etavsmix_YEta", Yeta, mix_Yeta, "#Eta (#Upsilon)", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot("Y_Rapidityvsmix_YRapidty", Yrapid, mix_Yrapid, "Rapidity (#Upsilon) [GeV]", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot("muon34_pt vs mixed_muon34_pt",hfourMuFitmu4Pt, h_mix_fourMuFitmu4Pt, "p_{T}(#muon) [GeV]", "Mixed(#muon)/#muon");
    ScaledPlot("muon34_eta vs mixed_muon34_eta",hfourMuFitmu4Eta, h_mix_fourMuFitmu4Eta, "#Eta (#muon) [GeV]", "Mixed(#muon)/#muon");
    ScaledPlot("muon34_y vs mixed_muon34_y",hfourMuFitmu4Rap, h_mix_fourMuFitmu4Rap, "Rapidty (#muon) [GeV]", "Mixed(#muon)/#muon");
    
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


