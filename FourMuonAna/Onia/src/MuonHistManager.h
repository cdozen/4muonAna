#ifndef FourmuonAnalysis_MuonHistManager4_H
#define FourmuonAnalysis_MuonHistManager4_H



// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"


class MuonHistManager{

  public:

  // constructor
  MuonHistManager();

  // destructor
  ~MuonHistManager();

    // write histograms the theFile
  void writeHists(TFile* theFile);

  // insert any TH1 into the big map
  void insertPlot(TH1* thePlot, std::string name, std::string folder);

  // fill 1D histogram 
  void fill1DHist(double x, double w , std::string name, std::string title,
                  int bins, double xmin, double xmax, std::string folder);

  // fill 2D histogram
  void fill2DHist(double x, double y, double w , 
                  std::string name, std::string title,
                  int binsx, double xmin, double xmax,
                  int binsy, double ymin, double ymax, std::string folder);


  // make a profile histogram
  void fillProfile(float x, float y, std::string name, std::string title, 
                   int binsx, float xmin, float xmax,
                   float ymin, float ymax, std::string folder);

  protected:

  private:

  // map to hold histograms
  std::map<std::string,std::pair<TH1*,std::string> > theMap;
  std::map<std::string,std::pair<TH2*,std::string> > theMap_2D;
  //std::map<std::string,std::pair<TH2D*,std::string> > theMap_2D;
  //std::map<std::string,std::pair<TH1D*,std::string> > theMap_1D;


};

#endif   
