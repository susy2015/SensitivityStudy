## Delete the btagging file
if (-f CSVv2_ichep.csv) then
  rm CSVv2_ichep.csv
endif

if (-f TTbarNoHad_bTagEff.root) then
  rm TTbarNoHad_bTagEff.root
endif

## Delete Pileup Reweighting file
if (-f PileupHistograms_Nov17.root) then
  rm PileupHistograms_Nov17.root
endif

## Delete the Top tagger files
if (-f Example_TopTagger.cfg) then
  rm Example_TopTagger.cfg
endif

if (-f Example_Legacy_TopTagger.cfg) then
  rm Example_Legacy_TopTagger.cfg
endif

if (-f TrainingOutput.model) then
  rm TrainingOutput.model
endif
