## Delete the btagging file
if (-f CSVv2_ichep.csv) then
  rm CSVv2_ichep.csv
endif

if (-f TTbarNoHad_bTagEff.root) then
  rm TTbarNoHad_bTagEff.root
endif

## Pileup Reweighting
if (-f PileupHistograms_Nov17.root) then
  rm PileupHistograms_Nov17.root
endif

## Delete the Top tagger files
if (-f ICHEPTaggerConfig.cfg) then
  rm ICHEPTaggerConfig.cfg
endif

if (-f TopTaggerConfig.cfg) then
  rm TopTaggerConfig.cfg
endif

if (-f TrainingOutput.model) then
  rm TrainingOutput.model
endif
