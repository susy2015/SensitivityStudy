#SensitivityStudy

```
cmsrel CMSSW_8_0_23
cd CMSSW_8_0_23/src/
cmsenv
```
Download source code from github and compile plugins:

SusyAnaTools:
```
git cms-init
git cms-merge-topic -u kpedro88:METfix8022
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate
git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_80X_V2
git clone -b Ana_Dec9_2016_Moriond2017_updatesForSBStudies_v2 git@github.com:susy2015/SusyAnaTools.git
```

TopTagger:
```
## Checkout and build OpenCV library for Toptagger
cd $CMSSW_BASE/src
git clone git@github.com:susy2015/opencv.git
cd $CMSSW_BASE/src/opencv
git checkout 3.1.0_StopBugFix
cmake .
make -j 8
## Checkout Tagtagger
cd $CMSSW_BASE/src
git clone -b HadStopAnaDevel_v4_Moriond2017_Dec17_2016 git@github.com:susy2015/TopTagger.git
```

CMS Build application:
```
cd $CMSSW_BASE/src
scram b -j 10
```

Build SusyAnaTools and TopTagger library:
```
cd $CMSSW_BASE/src/TopTagger/TopTagger/test
make -j 8
cd $CMSSW_BASE/src/SusyAnaTools/Tools
make
```
Please make sure compile the TopTagger first then SusyAnaTools/Tools! Since baselineDef.cc is rely on the toptagger library!

Sensitivity Study:
```
cd $CMSSW_BASE/src
git clone git@github.com:susy2015/SensitivityStudy.git
```

0.To produce the SSTree for a quick study:

```
cd SensitivityStudy/SSTreeMaker
source reset.csh
source $CMSSW_BASE/src/SusyAnaTools/Tools/setup.csh
$CMSSW_BASE/src/TopTagger/Tools/getTaggerCfg.sh -t MVAAK8_Medium_v1.1.0 -d /uscms_data/d3/hwei/stop
$CMSSW_BASE/src/TopTagger/Tools/getTaggerCfg.sh -t Legacy_AK4Only_v0.0.2 -f Legacy_TopTagger.cfg -d /uscms_data/d3/hwei/stop
make
```

Then run the type of MC you want to get SSTree(In priciple all of them!)

And then hadd and move them into EOS, change the runList files in the SensitivityStudy/SensitivityStudy directory

1.To study the CS in designed search bin:

```
cd SensitivityStudy/SensitivityStudy
./SS SSCS runList_Sensitivity_MC_SSSkimmed_v11p3b_BG.txt runList_Sensitivity_MC_SSSkimmed_v11p3b_SG.txt runList_Sensitivity_MC_SSSkimmed_v11p3b_MuCS.txt
```

2.To generate data card for all background:

```
cd SensitivityStudy/SensitivityStudy
./SS SSAllMC runList_Sensitivity_MC_SSSkimmed_v11p2_BG.txt runList_Sensitivity_MC_SSSkimmed_v11p2_SG.txt runList_Sensitivity_MC_SSSkimmed_v11p2_MuCS.txt
./SS SSAllMC runList_Sensitivity_MC_SSSkimmed_v11p2a_BG.txt runList_Sensitivity_MC_SSSkimmed_v11p2a_SG.txt runList_Sensitivity_MC_SSSkimmed_v11p2a_MuCS.txt
./SS SSAllMC runList_Sensitivity_MC_SSSkimmed_v11p2b_BG.txt runList_Sensitivity_MC_SSSkimmed_v11p2b_SG.txt runList_Sensitivity_MC_SSSkimmed_v11p2b_MuCS.txt
```

3.To generate data card for all signal points:

```
cd SensitivityStudy/SensitivityStudy
./SS SignalCardT2tt runList_Sensitivity_MC_SSSkimmed_v11p2_BG.txt runList_Sensitivity_MC_SSSkimmed_v11p2_SG.txt runList_Sensitivity_MC_SSSkimmed_v11p2_MuCS.txt
./SS SignalCardT1tttt runList_Sensitivity_MC_SSSkimmed_v11p2_BG.txt runList_Sensitivity_MC_SSSkimmed_v11p2_SG.txt runList_Sensitivity_MC_SSSkimmed_v11p2_MuCS.txt
```

4.To compare the real data card and fake data card

```
cd SensitivityStudy/SensitivityStudy
./SSDataCardCompare
```
