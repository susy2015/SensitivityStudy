#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$1/src/SensitivityStudy/SSTreeMaker/obj:$1/src/TopTagger/TopTagger/test:$1/src/opencv/lib

cd $1/src
eval `scramv1 runtime -sh`

cd ${_CONDOR_SCRATCH_DIR}
python NTuple_SSSignalScan.py $2

#find . -name "*.root" -exec xrdcp {} "root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/Skimmed_2015Nov15/" \;
#find . -name "*.root" -exec rm {} \;
