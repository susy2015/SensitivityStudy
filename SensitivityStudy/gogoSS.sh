#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

cd $5/src
eval `scramv1 runtime -sh`

cd ${_CONDOR_SCRATCH_DIR}

#find . -name "*.root" -exec rm {} \;
#find . -name "*.h" -exec rm {} \;
#find . -name "_*.*" -exec rm {} \;

#xrdcp root://cmsxrootd-site.fnal.gov//store/user/lpcsusyhad/PHYS14_720_Mar14_2014_v2/rootlist_$1.txt .
./SS $1 $2 $3 $4

#find . -name "*.root" -exec xrdcp {} "root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/AnaOut_QCD/" \;
#find . -name "*.h" -exec xrdcp {} "root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/AnaOut_QCD/" \;
#find . -name "_*.*" -exec xrdcp {} "root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/AnaOut_QCD/" \;
