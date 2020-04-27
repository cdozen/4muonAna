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
        // cout<<" Produced data ["<<i<<"] :"<<pdfObs->getRealValue("Pt_upsilon")<<"\n";
        TLorentzVector temp_vector_mix;
        temp_vector_mix.SetPtEtaPhiM(val, inputvector->Eta(), inputvector->Phi(), inputvector->M());
        if(!(abs(temp_vector_mix.Eta())> bin1 && abs(temp_vector_mix.Eta())< bin2)) continue; 
        outputVector.push_back(temp_vector_mix);
        count++;
        if (count==1000) break;
    }

}
//*************************************//
// Produce Random Distribution
//*************************************//

void mu4Random()
{   
    gROOT->SetBatch();
    define_and_setHistograms();

//*************************************//
//  Eta  -- Eta bins, 
//*************************************//    
    bool fullRange = false;
    const bool blind_signal = true;
    
//*************************************//
    vector<double> ybins;
        if (fullRange) ybins = {0, 2.5};
        if (!fullRange) ybins = {0, 0.35, 0.65, 0.9, 1.15, 1.4, 1.65, 2.0, 2.5};
        //if (!fullRange) ybins = {0, 0.7, 1.2, 1.6, 1.9, 2.2, 2.5, 2.9, 5.0};

        int y_size = ybins.size();

// S e t u p   c o m p o n e n t   p d f s 
// ---------------------------------------

    //TFile *f = new TFile("fourMuMass_tree.root");  //2018 with mu34 mas cut
    TFile *f = new TFile("/uscms_data/d3/cdozen/CMSSW_10_2_5/src/FourMuonAna/Onia/test/makeAllPlots/plots0p3/fourMuMass_tree.root");  //2018 with mu34 mas cut
    TTree *fTree = (TTree*) f->Get("FourMu_tree");

    FourMu_tree *ReadTree = new FourMu_tree(fTree);
   
// Construct observable
    
    RooRealVar *var_mu4_pt = new RooRealVar("mu4_pt", "mu4_pt", 0.0, 50, "GeV");
    RooRealVar *var_mu4_eta = new RooRealVar("mu4_eta", "mu4_eta", -5.0, 5.0, "GeV");
    RooArgSet *observables = new RooArgSet("Observables");
    observables->add(*var_mu4_pt);
    observables->add(*var_mu4_eta);
    
// Construct landau(t,ml,sl) ;

    RooRealVar ml("ml", "mean landau", 5., -0.0, 40);
    RooRealVar sl("sl", "sigma landau", 1, 0.1, 20);
    RooLandau landau("lx", "lx", *var_mu4_pt, ml, sl);

    
// Construct gauss(t,mg,sg)

    RooRealVar mg("mg", "mg", 0);
    RooRealVar sg("sg", "sg", 2, 0.01, 10);
    RooGaussian gauss("gauss", "gauss", *var_mu4_pt, mg, sg);

// ---------------------------------------
// C o n s t r u c t   ERF    p d f 
// ---------------------------------------

    RooRealVar* m0shift = new RooRealVar("m0shift", "m0shift", 1, 0., 5.);
    RooRealVar* width = new RooRealVar("width", "width", 0.5, 0., 5.);
    RooRealVar* par3 = new RooRealVar("par3", "par3", 5, 0., 15.);


//**************************************//
// Loop on Every rapidity bin
//**************************************//

    TString GraphTitle;
    for (int n = 0; n < y_size - 1; n++)
    {
        TCut mycuty;
            mycuty = Form("(abs(%s) > %.2f) && (abs(%s) < %.2f)", "mu4_eta", ybins[n], "mu4_eta", ybins[n + 1]);
        
        int Ntotal = fTree->GetEntries();
        RooDataSet *dataset = new RooDataSet("DataSet", "DataSet", *observables, RooFit::Import(*fTree), RooFit::Cut(mycuty));
        
        TString totalbin;
            totalbin = Form("%s_%.2f_%.2f", "mu4_eta", ybins[n], ybins[n + 1]);
        
        
//*******************************************************************//
///PDF Models
// Activate only 1 fit model
//*******************************************************************//      
        // Fitno1-) Construct landau (x) gauss

        //RooFFTConvPdf *lxg=new RooFFTConvPdf ("lxg","landau (X) gauss",*var_mu4_pt,landau,gauss) ;
        
        //*******************************************************************//      
        
        //Fitno2-) Construct Keys Pdf
        
        RooKeysPdf *lxg= new RooKeysPdf("KDE","KDE", *var_mu4_pt, *dataset, RooKeysPdf::MirrorBoth, 1);         
        
        //*******************************************************************//      
        
        //Fitno3-) Construct ERF
        
        //RooGenericPdf* lxg = new RooGenericPdf("ERF", "LikeSignPdf", "exp(-@0/par3)*(TMath::Erf((@0-m0shift)/width)+1)", RooArgList(*var_mu4_pt, *m0shift, *width, *par3));
        
        //*******************************************************************//      
        //*******************************************************************//      



       //*******************************************************************//
        /// Fitter ...
        //*******************************************************************//
        cout << "!!******************************!" << endl;
        cout << "Fitting : " << lxg->GetName() << "\n";
        cout << "!!******************************!" << endl;

        // ----------------------------------------------------------------------
        // S a m p l e ,   f i t   a n d   p l o t   c o n v o l u t e d   p d f 
        // ----------------------------------------------------------------------
        // Fit gxlx to data
        lxg->fitTo(*dataset, RooFit::NumCPU(4, RooFit::BulkPartition));

         // Set func for error from Muhammad
          TF1 *func ;
          func = lxg->asTF(RooArgList(*var_mu4_pt));
        


        // Plot data, landau pdf, landau (X) gauss pdf
            GraphTitle = lxg->GetName();
            GraphTitle += "_mu4_EtaBin";
        
        RooPlot* frame = var_mu4_pt->frame(Title(GraphTitle + totalbin));
        dataset->plotOn(frame);
        lxg->plotOn(frame);

//****************************//
//Dataset Generation
//****************************//

        cout << "!!***************************************!!\n";
        cout << " Readed   :" << Ntotal << "\n";
        cout << "!!***************************************!!\n";

        vector<TLorentzVector> mu12_MixFourMuVect;
        vector<TLorentzVector> mu12_MixingVect;
        vector<TLorentzVector> mu3_MixingVect; //for mu3 mixing
        vector<TLorentzVector> mu4_MixingVect; //for mu4 mixing
        vector<TLorentzVector*> mu12_OriginalVect;

        for (int i = 0; i < Ntotal; i++)
        {


            ReadTree->GetEntry(i);
            //  TCut mycutytree  = Form("(abs(%.2f) > %.2f) && (abs(%.2f) < %.2f)", ReadTree->Y_upsilon, ybins[n], ReadTree->Y_upsilon, ybins[n + 1]);
            //for eta binning
            
            if (!(abs(ReadTree->mu4_eta) > ybins[n] && abs(ReadTree->mu4_eta) < ybins[n + 1]))continue;
            
            cout << "\r" << "Entry :" << i << "/" << Ntotal << flush;
            
            //ProduceRandom(ReadTree->mu12,RandomYvector,datasetgen,lxg,i);
            ProduceRandom(ReadTree->mu4, mu4_MixingVect, func, i, ybins[n], ybins[n + 1]);
            // cout<<"finito.."<<endl;


///*******************************//
/// Fill Histograms
///*******************************//
            Ypt->Fill(ReadTree->mu12->Pt());
            Yeta->Fill(ReadTree->mu12->Eta());
            Yrapid->Fill(ReadTree->mu12->Rapidity());
            Yphi->Fill(ReadTree->mu12->Phi());
            hfourMuFitmu4Pt->Fill(ReadTree->mu4->Pt());
            hfourMuFitmu4Eta->Fill(ReadTree->mu4->Eta());
            hfourMuFitmu4Rap->Fill(ReadTree->mu4->Rapidity());
            Y2Dpteta->Fill(ReadTree->mu12->Eta(), ReadTree->mu12->Pt());
            Y2Detarapid->Fill(ReadTree->mu12->Eta(), ReadTree->mu12->Rapidity());
            Y2Dptrapid->Fill(ReadTree->mu12->Rapidity(), ReadTree->mu12->Pt());
            TLorentzVector mixFourMu;

            hfourMuMass_aftercut->Fill(ReadTree->Mass4mu);
            hfourMuMass_aftercut_smallrange->Fill(ReadTree->Mass4mu);
            
            
            
            //mix fourMuon
            for (auto s = 0; s < mu4_MixingVect.size(); s++)
                
            {
            if (!(abs(mu4_MixingVect[s].Eta()) > ybins[n] && abs(mu4_MixingVect[s].Eta()) < ybins[n + 1]))continue;
            if (!(abs((*ReadTree->mu3).Eta()) > ybins[n] && abs((*ReadTree->mu3).Eta()) < ybins[n + 1]))continue;
            if(!(abs(mu4_MixingVect[s].Eta())- abs((*ReadTree->mu3).Eta())< 0.1)) continue;
            
            
            TLorentzVector mu4Mixing= *ReadTree->mu3 + mu4_MixingVect[s];
                
            //if (!(abs(mu4Mixing.Eta()) > ybins[n] && abs(mu4Mixing.Eta()) < ybins[n + 1]))continue;
            //if (!(abs(mu4Mixing.Eta())- abs((*ReadTree->mu12).Eta())<0.5))continue;
                
                mixFourMu = *ReadTree->mu12 + mu4Mixing;
                mix_Ypt->Fill(mu4Mixing.Pt());
                mix_Yeta->Fill(mu4Mixing.Eta());
                mix_Yrapid->Fill(mu4Mixing.Rapidity());
                mix_Yphi->Fill(mu4Mixing.Phi());
                mix_Y2Dpteta->Fill(mu4Mixing.Eta(), mu4Mixing.Pt());
                mix_Y2Detarapid->Fill(mu4Mixing.Eta(), mu4Mixing.Rapidity());
                mix_Y2Dptrapid->Fill(mu4Mixing.Rapidity(), mu4Mixing.Pt());
                hfourMuMass_physkbg_mix->Fill(mixFourMu.M());
                hfourMuMass_physkbg_mix_smallrange->Fill(mixFourMu.M());
                //Ymumu2D_physkbg_mix->Fill(mu3Mixing.Eta(), ReadTree->mu34->Eta());
                h_mix_fourMuFitmu4Pt->Fill(mu4_MixingVect[s].Pt());
                h_mix_fourMuFitmu4Eta->Fill(mu4_MixingVect[s].Eta());
                h_mix_fourMuFitmu4Rap->Fill(mu4_MixingVect[s].Rapidity());
            }
            //

            //************************************//   
            //clear vect
            //************************************//   
            mu3_MixingVect.clear();
            mu4_MixingVect.clear();
            mu12_MixingVect.clear();
            mu12_OriginalVect.clear();


        }// Event loop


        //tata..



//************************************//   
//clear vect
//************************************//   
        TString MainfCanv = "fitplot_" + GraphTitle + totalbin;
        TCanvas *c1 = new TCanvas(MainfCanv, MainfCanv, 600, 600);
        gPad->SetLeftMargin(0.15);
        frame->GetYaxis()->SetTitleOffset(1.4);
        frame->Draw();
        c1->SaveAs(MainfCanv + ".png");

    } // End Loop on eta  bin
  
//***********************************************************//  
// Draw frame on canvas
//***********************************************************//  
    TString Cname = "Graph_total_" + GraphTitle;
    ;
    if (fullRange) Cname += "_FullRange_";
    if (!fullRange) Cname += "_BinnedTotal_";
          
    ScaledPlot(Cname, Ypt, mix_Ypt, "p_{T}(#Upsilon) [GeV]", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot("Y_Etavsmix_YEta", Yeta, mix_Yeta, "#Eta (#Upsilon)", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot("Y_Rapidityvsmix_YRapidty", Yrapid, mix_Yrapid, "Rapidity (#Upsilon) [GeV]", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot("muon_pt vs mixed_muon_pt",hfourMuFitmu4Pt, h_mix_fourMuFitmu4Pt, "p_{T}(#muon) [GeV]", "Mixed(#muon)/#muon");
    ScaledPlot("muon_eta vs mixed_muons_eta",hfourMuFitmu4Eta, h_mix_fourMuFitmu4Eta, "#Eta (#muon) [GeV]", "Mixed(#muon)/#muon");
    ScaledPlot("muon_y vs mixed_muon_y",hfourMuFitmu4Rap, h_mix_fourMuFitmu4Rap, "Rapidty (#muon) [GeV]", "Mixed(#muon)/#muon");
    
    ///draw Histograms:
    
    TCanvas *c2 = new TCanvas("c2","c2",800,600);

    Y2Dpteta->Draw("colz");
    c2->SaveAs("plots2D_mu4/Y2Dpteta.png");
    mix_Y2Dpteta->Draw("colz");
    c2->SaveAs("plots2D_mu4/Y2Dpteta_mix.png");
    Y2Dptrapid->Draw("colz");
    c2->SaveAs("plots2D_mu4/Y2Dptrapidity.png");
    Y2Detarapid->Draw("colz");
    c2->SaveAs("plots2D_mu4/Y2Detarapidity.png");
    mix_Y2Dptrapid->Draw("colz");
    c2->SaveAs("plots2D_mu4/Y2Dptrapidity_mix.png");
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
    c2->SaveAs("plots_mu4/hfourMuMass_origvsmix.png");
  
    //after cuts, orignal
    hfourMuMass_aftercut->Draw("e1");
    //c2->SaveAs("hfourMuMass_afterCut.png");
    float scaleEntries = hfourMuMass_aftercut_smallrange->Integral();
    hfourMuMass_aftercut_smallrange->SetStats(0);
    hfourMuMass_aftercut_smallrange->Draw("e1");
    //c2->SaveAs("hfourMuMass_afterCut_smallrange.png");
  /*hfourMuMass_aftercut_smallrange->SetBinContent(23,0);
   * hfourMuMass_aftercut_smallrange->SetBinContent(24,0);
   * hfourMuMass_aftercut_smallrange->SetBinContent(25,0);
   * hfourMuMass_aftercut_smallrange->SetBinContent(26,0);
   * hfourMuMass_aftercut_smallrange->SetBinContent(27,0);
   * hfourMuMass_aftercut_smallrange->SetBinContent(28,0);
   * hfourMuMass_aftercut_smallrange->SetBinContent(29,0);
   * hfourMuMass_aftercut_smallrange->SetBinContent(30,0);
   * hfourMuMass_aftercut_smallrange->SetBinContent(31,0);
   * hfourMuMass_aftercut_smallrange->SetBinContent(32,0);*/
  
    //after cuts, mix
    hfourMuMass_mix_aftercut->SetMarkerStyle(24);
    hfourMuMass_mix_aftercut->SetMarkerColor(kRed);
    hfourMuMass_mix_aftercut->Draw("e1");  
    //c2->SaveAs("hfourMuMass_mix_afterCut.png");

    //after cuts, compare original vs mix
    hfourMuMass_aftercut->Draw("e1");
    hfourMuMass_mix_aftercut->Draw("e1same");
    c2->SaveAs("plots_mu4/hfourMuMass_afterCut_origvsmix.png");
 
     //after cuts, mix physics bkg an zoom in
    hfourMuMass_physkbg_mix->SetMarkerStyle(22);
    hfourMuMass_physkbg_mix->SetMarkerColor(kBlue);
    hfourMuMass_physkbg_mix->Draw("e1");
    c2->SaveAs("plots_mu4/hfourMuMass_physkbg_mix.png");
    hfourMuMass_physkbg_mix_smallrange->SetMarkerStyle(22);
    hfourMuMass_physkbg_mix_smallrange->SetMarkerColor(kBlue);
    hfourMuMass_physkbg_mix_smallrange->SetMarkerStyle(22);
    hfourMuMass_physkbg_mix_smallrange->SetMarkerColor(kBlue);
    hfourMuMass_physkbg_mix_smallrange->Draw("e1");
    c2->SaveAs("plots_mu4/hfourMuMass_physkbg_mix_smallrange.png");
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
    c3->SaveAs("plots_mu4/hfourMuMass_afterCut_origVSphyskbg.png");
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
    c3->SaveAs("plots_mu4/hfourMuMass_afterCut_origVSphyskbg_smallrange.png");
   
   
    clearHistograms();
    fout->Write();
    f->Close();
    fout->Close();

}//end Testrandomproduction function



#endif /* TESTER_H */


