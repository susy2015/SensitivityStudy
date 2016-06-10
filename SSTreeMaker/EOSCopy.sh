#!/bin/bash
export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

#cd $1/src
eval `scramv1 runtime -sh`

for i in `xrdfs root://cmseos.fnal.gov/ ls /store/group/lpcsusyhad/Spring15_74X_Dec_2015_Ntp_v4X | grep 'txt'` 
do 
  xrdcp root://cmseos.fnal.gov/$i .
done

#xrdfs root://cmseos.fnal.gov/ ls /store/group/lpcsusyhad/hua/Skimmed_2015Nov15 | grep -E 'SSTrimmed_Spring15_74X_Feb_2016_Ntp_v6p0_forMoriond_TTJets_DiLeptstopFlatNtuples_.*root'
