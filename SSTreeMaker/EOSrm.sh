#!/bin/bash
export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

#cd $1/src
eval `scramv1 runtime -sh`

#for i in `xrdfs root://cmseos.fnal.gov/ ls /store/group/lpcsusyhad/hua/Skimmed_2015Nov15 | grep -E 'SSSignalScan_.*root'`
for i in `xrdfs root://cmseos.fnal.gov/ ls /store/group/lpcsusyhad/hua/Skimmed_2015Nov15 | grep -E 'SSTrimAndSlimmed_.*root'` 
#for i in `xrdfs root://cmseos.fnal.gov/ ls /store/group/lpcsusyhad/hua/Skimmed_2015Nov15/Sensitivity_MC_v11_p3a | grep -E 'SSTrimAndSlimmed_.*root'`
#for i in `xrdfs root://cmseos.fnal.gov/ ls /store/group/lpcsusyhad/hua/Skimmed_2015Nov15 | grep -E 'SFCheckTrimmed_.*root'`
do 
  eos root://cmseos.fnal.gov/ rm $i .
done
