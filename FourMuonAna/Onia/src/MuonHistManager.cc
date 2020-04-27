#include "FourMuonAna/Onia/src/MuonHistManager.h"
#include <algorithm>

using namespace std;


  MuonHistManager::MuonHistManager(){

    cout << "Initializing Histogram Manager..." << endl;

  }  


  MuonHistManager::~MuonHistManager(){

  }


  void MuonHistManager::writeHists(TFile* theFile){

    vector<string> theFolders;
    vector<string>::iterator fit;
    theFile->cd();

    map<string,pair<TH1*,string> >::const_iterator mapit;
    for (mapit = theMap.begin(); mapit != theMap.end(); mapit++){
      string folder = (*mapit).second.second.c_str();
      fit = find(theFolders.begin(), theFolders.end(), folder);
      if (fit == theFolders.end()){
        theFolders.push_back(folder);
        theFile->mkdir(folder.c_str());
      }
      theFile->cd((*mapit).second.second.c_str());
      (*mapit).second.first->Write();
      theFile->cd();
    }

    theFile->cd();
    map<string,pair<TH2*,string> >::const_iterator mapit_2D;
    for (mapit_2D = theMap_2D.begin(); mapit_2D != theMap_2D.end(); ++mapit_2D){
      string folder = (*mapit_2D).second.second.c_str();
      fit = find(theFolders.begin(), theFolders.end(), folder);
      if (fit == theFolders.end()){
        theFolders.push_back(folder);
        theFile->mkdir(folder.c_str());
      }
      theFile->cd((*mapit_2D).second.second.c_str());
      (*mapit_2D).second.first->Write();
      theFile->cd();
    }
  }


  void MuonHistManager::insertPlot(TH1* thePlot, string name, string folder){

    theMap[name] = pair<TH1*,string>(thePlot, folder);

  }
  


  void MuonHistManager::fill1DHist(double x, double w , string name,  string title ,
                               int bins , double xmin , double xmax , string folder ){

    map<string,pair<TH1*,string> >::iterator it;
    it = theMap.find(name);
    if (it == theMap.end()){
      theMap[name] = pair<TH1*,string>(new TH1D(name.c_str(),title.c_str(),bins,xmin,xmax), folder);
    }


    theMap[name].first->Fill(x,w);

  }


  void MuonHistManager::fill2DHist(double x, double y,  double w ,string name , string title ,
                               int binsx , double xmin , double xmax ,
                               int binsy , double ymin , double ymax , string folder ){

    map<string,pair<TH2*,string> >::iterator it;
    it = theMap_2D.find(name);
    if (it == theMap_2D.end()){
      theMap_2D[name] = pair<TH2*,string>(new TH2D(name.c_str(),title.c_str(),binsx,xmin,xmax,binsy,ymin,ymax),folder);
    }

    theMap_2D[name].first->Fill(x,y,w);

  }



  void MuonHistManager::fillProfile(float x, float y, string name, string title,
                                int binsx, float xmin, float xmax,
                                float ymin, float ymax, string folder){

    map<string,pair<TH1*,string> >::iterator it;
    it = theMap.find(name);
    if (it == theMap.end()){
      theMap[name] = pair<TProfile*,string>(new TProfile(name.c_str(),title.c_str(),binsx,xmin,xmax,ymin,ymax), folder);
    }

    theMap[name].first->Fill(x,y);

  }
