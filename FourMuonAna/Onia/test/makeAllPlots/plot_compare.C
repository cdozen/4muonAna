int plot_compare()
{
//TString histo = "hfourMuMass_aftercut";
//TString histo = "hfourMuMass_physkbg_mix";
//TString histo = "hfourMuMass_aftercut_smallrange";
//TString histo = "hfourMuMass_aftercut_smallrange_copy";
TString histo = "hfourMuMass_physkbg_mix_smallrange";
TCanvas* mycanvas = new TCanvas();
TFile *f1 = new TFile("plots_2017_rereco/fourMuMass.root");
TH1D *h1 = (TH1D*)f1->Get(histo);
h1->Scale(1/h1->Integral());
std::cout << "h1->Integral() = " << h1->Integral() << std::endl;
h1->SetLineColor(kRed);
h1->SetMarkerColor(kRed);
h1->Draw("e1");

TFile *f2 = new TFile("plots_2018_promptReco/fourMuMass.root");
//f2->ls();
TH1D *h2 = (TH1D*)f2->Get(histo);
h2->Scale(1/h2->Integral());
std::cout << "h2->Integral() = " << h2->Integral() << std::endl;
h2->SetLineColor(kBlue);
h2->SetMarkerColor(kBlue);
h2->Draw("same");

TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
leg->AddEntry(h1,"2017 mixed BG","lep");
leg->AddEntry(h2,"2018 mixed BG","lep");
leg->SetTextSize(0.04);
leg->Draw("same");

mycanvas->SaveAs("plots/2017_vs_2018_" + histo + ".png");
return 0;
}
