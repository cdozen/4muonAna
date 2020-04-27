#include "readTree2018.C"
void runReadTree2018()
{
   TChain *chain = new TChain("oniaTree","");
	//2018 Old files 
/*	chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2018Av1-MuMuGammaRootuple-2018Av1_promptReco.root");
	chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2018Av2-MuMuGammaRootuple-2018Av2_promptReco.root");
	chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2018Av3-MuMuGammaRootuple-2018Av3_promptReco.root");
	chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2018Bv1-MuMuGammaRootuple-2018B_promptReco.root");
	chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2018Bv2-MuMuGammaRootuple-2018B_promptReco.root");
	chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2018Cv1-MuMuGammaRootuple-2018C_promptReco.root");
	chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2018Cv2-MuMuGammaRootuple-2018C_promptReco.root");
	chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2018Cv3-MuMuGammaRootuple-2018C_promptReco.root");
	chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2018Dv2-MuMuGammaRootuple-2018D_promptReco.root");
*/	//new 2018 MuOnia 
	//chain->Add("root://cmsxrootd.fnal.gov//store/user/muahmad/FourMuon_Analysis/NTuples/Data_2018_17Sept/MuOnia/Run2018_NoDuplicates.root");
chain->Add("root://cmsxrootd.fnal.gov//store/user/muahmad/FourMuon_Analysis/NTuples/Data_2018_17Sept_v3/MuOnia/Run2018A_NoDuplicates.root");
chain->Add("root://cmsxrootd.fnal.gov//store/user/muahmad/FourMuon_Analysis/NTuples/Data_2018_17Sept_v3/MuOnia/Run2018B_NoDuplicates.root");
chain->Add("root://cmsxrootd.fnal.gov//store/user/muahmad/FourMuon_Analysis/NTuples/Data_2018_17Sept_v3/MuOnia/Run2018C_NoDuplicates.root");
chain->Add("root://cmsxrootd.fnal.gov//store/user/muahmad/FourMuon_Analysis/NTuples/Data_2018_17Sept_v3/MuOnia/Run2018D_NoDuplicates.root");
	//chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2017Bv1-MuMuGammaRootuple-2017_Rereco.root");
	//chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2017Cv1-MuMuGammaRootuple-2017_Rereco.root");
	//chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2017Dv1-MuMuGammaRootuple-2017_Rereco.root");
	//chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2017Ev1-MuMuGammaRootuple-2017_Rereco.root");
	//chain->Add("root://cmsxrootd.fnal.gov//store/user/huiwang/MuOnia/Run2017Fv1-MuMuGammaRootuple-2017_Rereco.root");

   readTree2018 a(chain);
   a.Loop();
}
