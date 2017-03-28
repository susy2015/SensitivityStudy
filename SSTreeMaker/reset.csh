## the following line can be commented out once the bug fixed in the SusyAnaTools/Tools/setup.csh
setenv LD_LIBRARY_PATH ./:${CMSSW_BASE}/src/opencv/lib/:${CMSSW_BASE}/src/TopTagger/TopTagger/test/:${CMSSW_BASE}/src/SusyAnaTools/Tools/obj/:${LD_LIBRARY_PATH}

## Unlink ISR Reweighting files
if (-f allINone_ISRJets.root) then
  unlink allINone_ISRJets.root
endif

if (-f ISRWeights.root) then
  unlink ISRWeights.root
endif

## Unlink btagging factor files
if (-f CSVv2_Moriond17_B_H.csv) then
  unlink CSVv2_Moriond17_B_H.csv
endif

if (-f allINone_bTagEff.root) then
  unlink allINone_bTagEff.root
endif

## Unlink lepton factor file
if (-f allINone_leptonSF_Moriond17.root) then
  unlink allINone_leptonSF_Moriond17.root
endif

## Unlink Pileup Reweighting file
if (-f PileupHistograms_0121_69p2mb_pm4p6.root) then
  unlink PileupHistograms_0121_69p2mb_pm4p6.root
endif

## Unlink W softdrop mass correction file
if (-f puppiCorr.root) then
  unlink puppiCorr.root
endif

## Unlink Top tagger related files
if (-f TopTagger.cfg) then
  unlink TopTagger.cfg
endif

if (-f TrainingOutput_dR20_pt30_depth12_500tree_noQGL_binaryCSV_2017_Mar24.model) then
  unlink TrainingOutput_dR20_pt30_depth12_500tree_noQGL_binaryCSV_2017_Mar24.model
endif

if (-f Legacy_TopTagger.cfg) then
  unlink Legacy_TopTagger.cfg
endif
