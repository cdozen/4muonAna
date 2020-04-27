/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   tester.h
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

void PlotHistograms()
{


}
//
/*
int countBin1=0;
int countBin2=0;
int countBin3=0;
int countBin4=0;
int countBin5=0;
int countBin6=0;
int countBin7=0;
//
int countInBin1=0;
int countInBin2=0;
int countInBin3=0;
int countInBin4=0;
int countInBin5=0;
int countInBin6=0;
int countInBin7=0;
*/ 
//
TF1 * Janjan(int n)
{
   // double ups_rap_value = mu12_p4.Rapidity();

    //double pT_rndm = 0;
    TF1 *f1;

    //if (abs(ups_eta_value)>0.0&&abs(ups_eta_value)<0.7)
    //if (abs(ups_rap_value) > 0.0 && abs(ups_rap_value) < 0.6)
    if (n==0)
    {  //  countInBin1++;
        //TF1 *f1 = new TF1("f1","exp(-x/9.20832286467)*(TMath::Erf((x-4.99999970055)/3.07316560256)+1)",0,50);
    //f1 = new TF1("f1", "exp(-x/7.43632004328)*(TMath::Erf((x-4.76905277365)/4.64398418226)+1)", 0, 50);
    //Ups. rapidy binning error func. 
    f1 = new TF1("f1", "exp(-x/7.43632004328)*(TMath::Erf((x-4.76905277365)/4.64398418226)+1)", 0, 50);
       
    }else
    //if (abs(ups_eta_value)>0.7&&abs(ups_eta_value)<1.2)
    //if (abs(ups_rap_value) > 0.6 && abs(ups_rap_value) < 0.9)
        if (n==1)
    {
       // countInBin2++;
        //TF1 *f1 = new TF1("f1","exp(-x/7.79191986237)*(TMath::Erf((x-4.99999951034)/2.48208577327)+1)",0,50);
        //f1 = new TF1("f1", "exp(-x/7.31342937311)*(TMath::Erf((x-1.54878888981)/2.25132915199)+1)", 0, 50);
    //Ups. rapidy binning error func. 
        f1 = new TF1("f1", "exp(-x/7.31342937311)*(TMath::Erf((x-1.54878888981)/2.25132915199)+1)", 0, 50);
       // pT_rndm = f1->GetRandom();
    }else
    //if (abs(ups_eta_value)>1.2&&abs(ups_eta_value)<1.6)
    //if (abs(ups_rap_value) > 0.9 && abs(ups_rap_value) < 1.1)
        if (n==2)
    {
       // countInBin3++;
        //TF1 *f1 = new TF1("f1","exp(-x/6.51155087776)*(TMath::Erf((x-4.84495922421)/2.74353432729)+1)",0,50);
        //f1 = new TF1("f1", "exp(-x/6.61877801248)*(TMath::Erf((x-1.47229454565)/1.59635908955)+1)", 0, 50);
    //Ups. rapidy binning error func. 
        f1 = new TF1("f1", "exp(-x/6.61877801248)*(TMath::Erf((x-1.47229454565)/1.59635908955)+1)", 0, 50);
       // pT_rndm = f1->GetRandom();
    }else
    //if (abs(ups_eta_value)>1.6&&abs(ups_eta_value)<1.9)
    //if (abs(ups_rap_value) > 1.1 && abs(ups_rap_value) < 1.3)
        if (n==3)
    { //countInBin4++;
        //TF1 *f1 = new TF1("f1","exp(-x/7.49015940087)*(TMath::Erf((x-2.63184664589)/1.42873868111)+1)",0,50);
        //f1 = new TF1("f1", "exp(-x/5.92246543461)*(TMath::Erf((x-1.34362745427)/0.872836572706)+1)", 0, 50);
    //Ups. rapidy binning error func. 
        f1 = new TF1("f1", "exp(-x/5.92246543461)*(TMath::Erf((x-1.34362745427)/0.872836572706)+1)", 0, 50);
      //  pT_rndm = f1->GetRandom();
    }else
    //if (abs(ups_eta_value)>1.9&&abs(ups_eta_value)<2.2)
    //if (abs(ups_rap_value) > 1.3 && abs(ups_rap_value) < 1.5)
        if (n==4)
    {
       // countInBin5++;
        //TF1 *f1 = new TF1("f1","exp(-x/5.03883937761)*(TMath::Erf((x-2.6210839199)/1.23349333644)+1)",0,50);
        //f1 = new TF1("f1", "exp(-x/5.84513508567)*(TMath::Erf((x-1.24642564857)/1.03897957692)+1)", 0, 50);
    //Ups. rapidy binning error func. 
        f1 = new TF1("f1", "exp(-x/5.84513508567)*(TMath::Erf((x-1.24642564857)/1.03897957692)+1)", 0, 50);
       // pT_rndm = f1->GetRandom();
    }else
    //if (abs(ups_eta_value)>2.2&&abs(ups_eta_value)<2.5)
    //if (abs(ups_rap_value) > 1.5 && abs(ups_rap_value) < 1.8)
        if (n==5)
    {   
        //countInBin6++;
        //TF1 *f1 = new TF1("f1","exp(-x/3.06152702812)*(TMath::Erf((x-2.57072754024)/1.14786510452)+1)",0,50);
        //f1 = new TF1("f1", "exp(-x/5.38129850368)*(TMath::Erf((x-1.51008066104)/0.986344895293)+1)", 0, 50);
    //Ups. rapidy binning error func. 
        f1 = new TF1("f1", "exp(-x/5.38129850368)*(TMath::Erf((x-1.51008066104)/0.986344895293)+1)", 0, 50);
      //  pT_rndm = f1->GetRandom();
    }else
    //if (abs(ups_eta_value)>2.5&&abs(ups_eta_value)<2.9)
    //if (abs(ups_rap_value) > 1.8 && abs(ups_rap_value) < 2.4)
        if (n==6)
    {   //countInBin7++;
        //  TF1 *f1 = new TF1("f1","exp(-x/1.51109342873)*(TMath::Erf((x-2.80111484471)/1.278065495)+1)",0,50);
        //f1 = new TF1("f1", "exp(-x/3.9449553223)*(TMath::Erf((x-2.10744727896)/1.45974398534)+1)", 0, 50);
    //Ups. rapidy binning error func. 
        f1 = new TF1("f1", "exp(-x/3.9449553223)*(TMath::Erf((x-2.10744727896)/1.45974398534)+1)", 0, 50);
      //  pT_rndm = f1->GetRandom();
    }
   //  pT_rndm = f1->GetRandom();
   // f1->Delete();
    return f1;

    
    //
}

template <class T>
void ProduceRandom(TLorentzVector *mu34, // Input TLorentzVector
                   TLorentzVector *outputVector, // Output vector of TLorentzVector
                   RooDataSet *datasetgen, //
                   //RooArgSet *observables,                // Observable
                   T lxg, // Convoluted fit PDF
                   //RooKeysPdf* lxg,                     // Convoluted fit PDF
                   int i, // Event counter
                   double bin1,
                   double bin2
                   )
{


    RooArgSet* pdfObs = lxg->getObservables(*datasetgen);

    *pdfObs = *datasetgen->get(i);
    // cout<<" Produced data ["<<i<<"] :"<<pdfObs->getRealValue("Pt_upsilon")<<"\n";
    // TLorentzVector mu12_p4_vector_mix;
    outputVector->SetPtEtaPhiM(pdfObs->getRealValue("Pt34_upsilon"), mu34->Eta(), mu34->Phi(), mu34->M());
    // outputVector.push_back(mu12_p4_vector_mix);


}
//template <class T> 

void ProduceRandom(TLorentzVector *mu34, // Input TLorentzVector
                   vector<TLorentzVector> &outputVector, // Output vector of TLorentzVector
                   //RooDataSet *datasetgen,                //
                   //RooArgSet *observables,                // Observable
                   TF1 *func, // Convoluted fit PDF
                   //RooKeysPdf* lxg,                     // Convoluted fit PDF
                   int i // Event counter
                   )
{


    //RooArgSet* pdfObs = lxg->getObservables(*datasetgen) ;

    // *pdfObs = *datasetgen->get(i);
    for (int k = 0; k < 1000; k++)
    {
        double val=0;
    

            val = func->GetRandom();
         //  cout<<val<<" ---- "<<i<<" \n";
        // cout<<" Produced data ["<<i<<"] :"<<pdfObs->getRealValue("Pt_upsilon")<<"\n";
        TLorentzVector mu34_p4_vector_mix;
        mu34_p4_vector_mix.SetPtEtaPhiM(val, mu34->Eta(), mu34->Phi(), mu34->M());
        outputVector.push_back(mu34_p4_vector_mix);
    }

}
//*************************************//
// Produce Random Distribution
//*************************************//

void testRandomProd_v2()
{
    gROOT->SetBatch();
    define_and_setHistograms();
    //*************************************//
    //  Set_Eta  -- true Eta bins, false Rapidity bins
    //*************************************//    
    bool Set_Eta = false;
    bool fullRange = false;
    bool test_janjan= false; // Muhammad error fuc±..
    const bool blind_signal = true;
    
    //*************************************//
    vector<double> ybins;
    // Define rapidity bins ..
    // double ybins[]= { 0,0.7,1.2,1.6,1.9,2.2,2.5};
    if (!Set_Eta)
    {
        if (fullRange)ybins = {0, 2.5}; // rapidity bins
        //if (!fullRange)ybins = {0, 0.6, 0.9, 1.1, 1.3, 1.5, 1.8, 2.5}; // rapidity bins
        if (!fullRange)ybins = {0, 0.5, 0.85, 1.1, 1.3, 1.5, 1.7, 1.9, 2.5}; // rapidity bins
    }
    else
    {
        if (fullRange) ybins = {0, 5.0};
        if (!fullRange) ybins = {0, 0.7, 1.2, 1.6, 1.9, 2.2, 2.5, 2.9, 5.0};
        //if (!fullRange) ybins = {0, 0.5, 0.8, 1.0, 1.3, 1.6, 1.9, 2.4, 5.0};
    }


    // double ybins[] = {0.0, 2.5};
    int y_size = ybins.size();

    // S e t u p   c o m p o n e n t   p d f s 
    // ---------------------------------------

    //TFile *f = new TFile("fourMuMass_tree34.root");
    TFile *f = new TFile("fourMuMass_treemu34on.root"); //with mu34 mass cut
    TTree *fTree = (TTree*) f->Get("FourMu_tree");

    FourMu_tree *ReadTree = new FourMu_tree(fTree);
    // Construct observable
    RooRealVar *var_Pt34_upsilon = new RooRealVar("Pt34_upsilon", "Pt34_upsilon", 0.0, 50, "GeV");
    RooRealVar *var_Y34_upsilon = new RooRealVar("Y34_upsilon", "Y34_upsilon", -5.0, 5.0, "GeV");
    RooRealVar *var_Eta34_upsilon = new RooRealVar("Eta34_upsilon", "Eta34_upsilon", -5.0, 5.0, "GeV");
    RooArgSet *observables = new RooArgSet("Observables");
    observables->add(*var_Pt34_upsilon);
    if (!Set_Eta) observables->add(*var_Y34_upsilon);
    if (Set_Eta) observables->add(*var_Eta34_upsilon);

    // Construct landau(t,ml,sl) ;
    RooRealVar ml("ml", "mean landau", 5., -0.0, 40);
    RooRealVar sl("sl", "sigma landau", 1, 0.1, 20);
    RooLandau landau("lx", "lx", *var_Pt34_upsilon, ml, sl);

    // Construct gauss(t,mg,sg)
    RooRealVar mg("mg", "mg", 0);
    RooRealVar sg("sg", "sg", 2, 0.01, 10);
    RooGaussian gauss("gauss", "gauss", *var_Pt34_upsilon, mg, sg);

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
        if (!Set_Eta)
        {
            //mycuty = Form("(abs(%s) > %.2f) && (abs(%s) < %.2f)", "Y_upsilon", ybins[n], "Y_upsilon", ybins[n + 1]);
            mycuty = Form("(abs(%s) > %.2f) && (abs(%s) <= %.2f)", "Y34_upsilon", ybins[n], "Y34_upsilon", ybins[n + 1]);
        }
        else
        {
            mycuty = Form("(abs(%s) > %.2f) && (abs(%s) < %.2f)", "Eta34_upsilon", ybins[n], "Eta34_upsilon", ybins[n + 1]);
        }
        int Ntotal = fTree->GetEntries();
        RooDataSet *dataset = new RooDataSet("DataSet", "DataSet", *observables, RooFit::Import(*fTree), RooFit::Cut(mycuty));
        TString totalbin;
        if (!Set_Eta)
        {
            totalbin = Form("%s_%.2f_%.2f", "Y34_rapidity", ybins[n], ybins[n + 1]);
        }
        else
        {
            totalbin = Form("%s_%.2f_%.2f", "Y34_Eta", ybins[n], ybins[n + 1]);
        }
        //*******************************************************************//
        ///PDF Models
        // Activate only 1 fit model
        //*******************************************************************//      
        // Fitno1-) Construct landau (x) gauss
         // RooFFTConvPdf *lxg=new RooFFTConvPdf ("lxg","landau (X) gauss",*var_Pt34_upsilon,landau,gauss) ;
        //*******************************************************************//      
        //Fitno2-) Construct Keys Pdf
          //RooKeysPdf *lxg= new RooKeysPdf("KDE","KDE", *var_Pt_upsilon, *dataset, RooKeysPdf::MirrorBoth, 1);         
        //*******************************************************************//      
        //Fitno3-) Construct ERF
         RooGenericPdf* lxg = new RooGenericPdf("ERF", "LikeSignPdf", "exp(-@0/par3)*(TMath::Erf((@0-m0shift)/width)+1)", RooArgList(*var_Pt34_upsilon, *m0shift, *width, *par3));
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

         // Fonksiyon set ettik
          TF1 *func ;
          if(test_janjan)
          { 
              func = Janjan(n);
          }
          else{
          func = lxg->asTF(RooArgList(*var_Pt34_upsilon));
          }
        


        // Plot data, landau pdf, landau (X) gauss pdf
        GraphTitle = lxg->GetName();
        if (!Set_Eta)
        {
            GraphTitle += "_Y_Pt34_RapidityBin";
        }
        else
        {
            GraphTitle += "_Y_Pt34_EtaBin";
        }
        RooPlot* frame = var_Pt34_upsilon->frame(Title(GraphTitle + totalbin));
        dataset->plotOn(frame);
        lxg->plotOn(frame);
        //keys_pdfs->plotOn(frame, RooFit::LineColor(kGreen));

        //****************************//
        //Dataset Generation
        //****************************//

        // RooDataSet *datasetgen= lxg->generate(*observables,Ntotal);
        //  datasetgen->plotOn(frame,LineStyle(kDashed),MarkerColor(kRed));
        cout << "!!***************************************!!\n";
        cout << " Readed   :" << Ntotal << "\n";
        //cout<<" Produced :"<<datasetgen->numEntries()<<"\n";
        cout << "!!***************************************!!\n";

        vector<TLorentzVector> mu34_MixFourMuVect;
        vector<TLorentzVector> mu34_MixingVect;
        vector<TLorentzVector*> mu12_OriginalVect;
        vector<TLorentzVector*> mu34_OriginalVect;
        vector<TLorentzVector*> mu1_OriginalVect;
        vector<TLorentzVector*> mu2_OriginalVect;
        vector<TLorentzVector*> mu3_OriginalVect;
        vector<TLorentzVector*> mu4_OriginalVect;
        //TLorentzVector *inputYvector;
        //TLorentzVector *RandomYvector=new TLorentzVector();

        //   TLorentzVector *mixFourMu= new TLorentzVector();   
        for (int i = 0; i < Ntotal; i++)
        {


            ReadTree->GetEntry(i);
            //  TCut mycutytree  = Form("(abs(%.2f) > %.2f) && (abs(%.2f) < %.2f)", ReadTree->Y_upsilon, ybins[n], ReadTree->Y_upsilon, ybins[n + 1]);
            //for rpaidiy binning 
            if(!Set_Eta) if (!(abs(ReadTree->Y34_upsilon) > ybins[n] && abs(ReadTree->Y34_upsilon) < ybins[n + 1]))continue;
            if(Set_Eta) if (!(abs(ReadTree->Eta34_upsilon) > ybins[n] && abs(ReadTree->Eta34_upsilon) < ybins[n + 1]))continue;
            cout << "\r" << "Entry :" << i << "/" << Ntotal << flush;
            //ProduceRandom(ReadTree->mu12,RandomYvector,datasetgen,lxg,i);
            ProduceRandom(ReadTree->mu34, mu34_MixingVect, func, i);
            // cout<<"finito.."<<endl;

            mu34_OriginalVect.push_back(ReadTree->mu34);
            // mu12_MixingVect.push_back(RandomYvector);

            //  cout<<"Pt :"<<ReadTree->fourMuFit->Pt()<<" Eta :"<<ReadTree->fourMuFit->Eta()<<" Rapidity: "<<ReadTree->fourMuFit->Rapidity()<<" Mass:"<<ReadTree->fourMuFit->M()<<endl;
            //  cout<<"Pt :"<<mixFourMu->Pt()          <<" Eta :"<<mixFourMu->Eta()          <<" Rapidity: "<<mixFourMu->Rapidity()          <<" Mass:"<<mixFourMu->M()<<endl;        

            ///*******************************//
            /// Fill Histograms
            ///*******************************//
            Ypt->Fill(ReadTree->mu34->Pt());
            Yeta->Fill(ReadTree->mu34->Eta());
            Yrapid->Fill(ReadTree->mu34->Rapidity());
            Yphi->Fill(ReadTree->mu34->Phi());
            TLorentzVector mixFourMu;

            hfourMuMass_aftercut->Fill(ReadTree->Mass4mu);
            hfourMuMass_aftercut_smallrange->Fill(ReadTree->Mass4mu);
            
            
            
            //mix
            for (auto s = 0; s < mu34_MixingVect.size(); s++)
            {
                mixFourMu = *ReadTree->mu12 + mu34_MixingVect[s];
                // mu12_MixFourMuVect.push_back(mixFourMu);
                mix_Ypt->Fill(mu34_MixingVect[s].Pt());
                mix_Yeta->Fill(mu34_MixingVect[s].Eta());
                mix_Yrapid->Fill(mu34_MixingVect[s].Rapidity());
                mix_Yphi->Fill(mu34_MixingVect[s].Phi());
                mix_Y2Dpteta->Fill(mu34_MixingVect[s].Pt(), mu34_MixingVect[s].Eta());
                mix_Y2Dptrapid->Fill(mu34_MixingVect[s].Pt(), mu34_MixingVect[s].Rapidity());
                hfourMuMass_physkbg_mix->Fill(mixFourMu.M());
                hfourMuMass_physkbg_mix_smallrange->Fill(mixFourMu.M());
                Ymumu2D_physkbg_mix->Fill(mu34_MixingVect[s].Eta(), ReadTree->mu12->Eta());
                //mixFourMu->delete;

            }
            //
            //  hfourMuMass_physkbg_mix_smallrange->Fill(mixFourMu->M());

            //************************************//   
            //clear vect
            //************************************//   

            mu34_MixingVect.clear();
            mu34_OriginalVect.clear();


        }// Event loop

        //TLorentzVector *RandomYvector=new TLorentzVector();

        //tata..



        //************************************//   
        //clear vect
        //************************************//   
        /*
        mu12_MixFourMuVect.clear();
        mu34_OriginalVect.clear();
        mu1_OriginalVect.clear();
        mu2_OriginalVect.clear();
        mu3_OriginalVect.clear();
        mu4_OriginalVect.clear();*/
        //TLorentzVector *inputYvector;
        TString MainfCanv = "fitplot_" + GraphTitle + totalbin;
        TCanvas *c1 = new TCanvas(MainfCanv, MainfCanv, 600, 600);
        gPad->SetLeftMargin(0.15);
        frame->GetYaxis()->SetTitleOffset(1.4);
        frame->Draw();
        c1->SaveAs(MainfCanv + ".png");

    } // End Loop on Rapidity bin
    //***********************************************************//  
    // Draw frame on canvas
    //***********************************************************//  
    TString Cname = "Graph_total_" + GraphTitle;
    ;
    if (fullRange) Cname += "_FullRange_";
    if (!fullRange) Cname += "_BinnedTotal_";
    if(test_janjan) Cname += "Test_Formulas";
          
    ScaledPlot(Cname, Ypt, mix_Ypt, "p_{T}(#Upsilon) [GeV]", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot("Y_Etavsmix_YEta", Yeta, mix_Yeta, "#Eta (#Upsilon)", "Mixed(#Upsilon)/#Upsilon");
    ScaledPlot("Y_Rapidityvsmix_YRapidty", Yrapid, mix_Yrapid, "Rapidity (#Upsilon) [GeV]", "Mixed(#Upsilon)/#Upsilon"); 

    
    ///draw Histograms:
    TCanvas *c2 = new TCanvas("c2","c2",800,600);
   
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
    c2->SaveAs("hfourMuMass_origvsmix.png");
  
    //after cuts, orignal
    hfourMuMass_aftercut->Draw("e1");
    //c2->SaveAs("hfourMuMass_afterCut.png");
    float scaleEntries = hfourMuMass_aftercut_smallrange->Integral();
    hfourMuMass_aftercut_smallrange->SetStats(0);
    hfourMuMass_aftercut_smallrange->Draw("e1");
    //c2->SaveAs("hfourMuMass_afterCut_smallrange.png");
    //
  
  
    //after cuts, mix
    hfourMuMass_mix_aftercut->SetMarkerStyle(24);
    hfourMuMass_mix_aftercut->SetMarkerColor(kRed);
    hfourMuMass_mix_aftercut->Draw("e1");
    //c2->SaveAs("hfourMuMass_mix_afterCut.png");

    //after cuts, compare original vs mix
    hfourMuMass_aftercut->Draw("e1");
    hfourMuMass_mix_aftercut->Draw("e1same");
    c2->SaveAs("hfourMuMass_afterCut_origvsmix.png");
 
     //after cuts, mix physics bkg an zoom in
    hfourMuMass_physkbg_mix->SetMarkerStyle(22);
    hfourMuMass_physkbg_mix->SetMarkerColor(kBlue);
    hfourMuMass_physkbg_mix->Draw("e1");
    c2->SaveAs("hfourMuMass_physkbg_mix.png");
    hfourMuMass_physkbg_mix_smallrange->SetMarkerStyle(22);
    hfourMuMass_physkbg_mix_smallrange->SetMarkerColor(kBlue);
    hfourMuMass_physkbg_mix_smallrange->SetMarkerStyle(22);
    hfourMuMass_physkbg_mix_smallrange->SetMarkerColor(kBlue);
    hfourMuMass_physkbg_mix_smallrange->Draw("e1");
    c2->SaveAs("hfourMuMass_physkbg_mix_smallrange.png");
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
    c3->SaveAs("hfourMuMass_afterCut_origVSphyskbg.png");
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
    c3->SaveAs("hfourMuMass_afterCut_origVSphyskbg_smallrange.png");
   
   
   
    clearHistograms();
    fout->Write();
    f->Close();
    fout->Close();

}//end Testrandomproduction function



#endif /* TESTER_H */


