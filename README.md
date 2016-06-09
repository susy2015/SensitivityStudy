# SensitivityStudy

cmsrel CMSSW_8_0_10
cd CMSSW_8_0_10/src/
cmsenv
git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git
git clone -b git@github.com:susy2015/SusyAnaTools.git
git clone git@github.com:susy2015/SensitivityStudy.git
scram b -j9
