# FourMuonAna
CMSSW version for 2016 : CMSSW_8_0_X 
#####
cmsrel CMSSW_8_0_29
####
cd CMSSW_8_0_29/src
###
cmsenv
#####
CMSSW version for 2017: CMSSW_9_4_X 
###
cmsrel CMSSW_9_4_2
####
cd CMSSW_9_4_2/src
###
cmsenv
####

CMSSW  version for 2018:  CMSSW_10_2_X
###
cmsrel CMSSW_10_2_5
#####
cd CMSSW_10_2_5/src


git clone https://github.com/zhenhu/FourMuonAna.git

#####

COMPILE:

scram b -j28

####

# Global Tags :
2016-->80X_dataRun2_2016SeptRepro_v7


2017-->94X_dataRun2_ReReco_EOY17_v6

# 2018 Prompt Reco

2018D-->102X_dataRun2_Prompt_v4

# 2018 ReReco

2018ABC--->102X_dataRun2_Sep2018ABC_v2










