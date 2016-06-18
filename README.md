#SensitivityStudy

cmsrel CMSSW_8_0_10

cd CMSSW_8_0_10/src/

cmsenv

git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git

git clone -b Ana_Prod_merged_June17_2016_fix_top_projection_bug_data_topoff git@github.com:susy2015/SusyAnaTools.git

git clone git@github.com:susy2015/SensitivityStudy.git

scram b -j9

cd SensitivityStudy/SensitivityStudy

xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/hua/Skimmed_2015Nov15/Sensitivity_MC_v6/signalScan_SMS-T1tttt_forHua.root SignalScanBeforeBaseline/

xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/hua/Skimmed_2015Nov15/Sensitivity_MC_v6/signalScan_SMS-T2tt_forHua.root SignalScanBeforeBaseline/

0.To produce the SSTree for a quick study:

cd SSTreeMaker

Make

Then run the type of MC you want to get SSTree(In priciple all of them!)

And then hadd and move them into EOS, change the runList files in the SensitivityStudy/SensitivityStudy directory

1.To study Signal/MC in designed search bin:

./SS SSAllMC runList_Sensitivity_MC_SSSkimmed_v6_BG.txt runList_Sensitivity_MC_SSSkimmed_v6_SG.txt runList_Sensitivity_MC_SSSkimmed_v6_MuCS.txt

2.To study the CS in designed search bin:

./SS SSCS runList_Sensitivity_MC_SSSkimmed_v6_BG.txt runList_Sensitivity_MC_SSSkimmed_v6_SG.txt runList_Sensitivity_MC_SSSkimmed_v6_MuCS.txt

3.To generate Signal Data Card

./SS SignalCardT2tt runList_Sensitivity_MC_SSSkimmed_v6_BG.txt runList_Sensitivity_MC_SSSkimmed_v6_SG.txt runList_Sensitivity_MC_SSSkimmed_v6_MuCS.txt

./SS SignalCardT1tttt runList_Sensitivity_MC_SSSkimmed_v6_BG.txt runList_Sensitivity_MC_SSSkimmed_v6_SG.txt runList_Sensitivity_MC_SSSkimmed_v6_MuCS.txt

3.To make 2D plots in designed search bin:

./SSPlots SSPlots20160504

and 1D plots to show the signal and BG distribution:

./SSAUX1DPlots SSAUX1DPlots20160517

4.To compare the real data card and fake data card

./SSDataCardCompare

5.To test the SBGeometry.h(Need to be improved!)

g++ -std=c++11 SSBinIDTest.cc

./a.out

