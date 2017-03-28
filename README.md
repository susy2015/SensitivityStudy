#SensitivityStudy

```
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src/
cmsenv
```
Download source code from github and compile plugins:

SusyAnaTools:
```
git cms-init
git cms-merge-topic cms-met:METRecipe_8020 -u
git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git
git clone git@github.com:susy2015/JetToolbox.git JMEAnalysis/JetToolbox -b fix_NoLep_jetToolbox_80X_V3
git clone -b ana_reMINIAOD_Mar06_2017 git@github.com:susy2015/SusyAnaTools.git
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
git clone -b HadStopAnaDevel_v8_Moriond2017_Mar27_2017 git@github.com:susy2015/TopTagger.git
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
Note, setup.csh not updated with new cfg file, change in local with new cfg tag

```
cd SensitivityStudy/SSTreeMaker
source reset.csh
source $CMSSW_BASE/src/SusyAnaTools/Tools/setup.csh
sh cache_all.sh
tar --exclude-caches-all -zcf ${CMSSW_VERSION}.tar.gz -C ${CMSSW_BASE}/.. ${CMSSW_VERSION}
tar -zcf SensitivityTxt.tar.gz SensitivityTxt
make
```

Manually script:
```
$CMSSW_BASE/src/TopTagger/Tools/getTaggerCfg.sh -t MVAAK8_Tight_noQGL_binaryCSV_v1.0.2 -d /uscms_data/d3/hwei/stop
$CMSSW_BASE/src/TopTagger/Tools/getTaggerCfg.sh -t Legacy_AK4Only_v0.0.2 -f Legacy_TopTagger.cfg -d /uscms_data/d3/hwei/stop
```

Then run the type of MC you want to get SSTree(In priciple all of them!)

And then hadd and move them into EOS, change the runList files in the SensitivityStudy/SensitivityStudy directory

1.To study the CS in designed search bin:

```
cd SensitivityStudy/SensitivityStudy
./SS SSCS runList_Sensitivity_MC_SSSkimmed_v11p5d_BG.txt runList_Sensitivity_MC_SSSkimmed_v11p5d_SG.txt runList_Sensitivity_MC_SSSkimmed_v11p5d_MuCS.txt
```

2.To generate data card for all background:

```
cd SensitivityStudy/SensitivityStudy
./SS SSAllMC runList_Sensitivity_MC_SSSkimmed_v11p5d_BG.txt runList_Sensitivity_MC_SSSkimmed_v11p5d_SG.txt runList_Sensitivity_MC_SSSkimmed_v11p5d_MuCS.txt
```

3.To generate data card for all signal points:

```
cd SensitivityStudy/SensitivityStudy
./SS SignalCardT2tt runList_Sensitivity_MC_SSSkimmed_v11p5d_BG.txt runList_Sensitivity_MC_SSSkimmed_v11p5d_SG.txt runList_Sensitivity_MC_SSSkimmed_v11p5d_MuCS.txt
./SS SignalCardT1tttt runList_Sensitivity_MC_SSSkimmed_v11p5d_BG.txt runList_Sensitivity_MC_SSSkimmed_v11p5d_SG.txt runList_Sensitivity_MC_SSSkimmed_v11p5d_MuCS.txt
```

4.To compare the real data card and fake data card

```
cd SensitivityStudy/SensitivityStudy
./SSDataCardCompare
```
