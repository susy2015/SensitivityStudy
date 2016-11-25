#SensitivityStudy

```
cmsrel CMSSW_8_0_23
cd CMSSW_8_0_23/src/
cmsenv
git cms-init
git cms-merge-topic -u kpedro88:METfix8022
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate
git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_80X_V2
git clone -b Moriond2017 git@github.com:susy2015/SusyAnaTools.git
git clone git@github.com:susy2015/SensitivityStudy.git
scram b -j9
cd SensitivityStudy/SensitivityStudy
xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/hua/Skimmed_2015Nov15/Sensitivity_MC_v6/signalScan_SMS-T1tttt_forHua.root SignalScanBeforeBaseline/
xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/hua/Skimmed_2015Nov15/Sensitivity_MC_v6/signalScan_SMS-T2tt_forHua.root SignalScanBeforeBaseline/
source rmsetup.csh
source setup.csh
```

To Checkout TopTagger Code:
```
## Checkout OpenCV
cd $CMSSW_BASE/src
git clone git@github.com:susy2015/opencv.git
cd opencv
git checkout 3.1.0_StopBugFix
cmake .
make -j 8
## Checkout Tagtagger
cd $CMSSW_BASE/src
git clone git@github.com:susy2015/TopTagger.git
scram b -j 8
cd TopTagger/TopTagger/test/
make -j 8
```

You can then compile the SUSYAnaTools
```
cd $CMSSW_BASE/src/SusyAnaTools/Tools/
make
```

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

