## setenv LD_LIBRARY_PATH ./:${CMSSW_BASE}/src/opencv/lib/:${CMSSW_BASE}/src/TopTagger/TopTagger/test/:${CMSSW_BASE}/src/SusyAnaTools/Tools/obj/:${LD_LIBRARY_PATH}
## Delete the btagging file
if (-f CSVv2_ichep.csv) then
  unlink CSVv2_ichep.csv
endif

if (-f TTbarNoHad_bTagEff.root) then
  unlink TTbarNoHad_bTagEff.root
endif

## Delete Pileup Reweighting file
if (-f PileupHistograms_Nov17.root) then
  unlink PileupHistograms_Nov17.root
endif

## Delete the Top tagger files
if (-f TopTagger.cfg) then
  unlink TopTagger.cfg
endif

if (-f TrainingOutput_dR20_pt30_depth14_2016_Dec2.model) then
  unlink TrainingOutput_dR20_pt30_depth14_2016_Dec2.model
endif

if (-f Legacy_TopTagger.cfg) then
  unlink Legacy_TopTagger.cfg
endif
