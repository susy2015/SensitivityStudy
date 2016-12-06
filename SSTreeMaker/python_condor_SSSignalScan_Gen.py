import sys
import glob, os

def PrintCondorHeaderLine():
  print("universe = vanilla")
  print("request_disk   = 50 GB")
  print("request_memory = 1.5 GB")
  print("executable = $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/goSSSignalScan.sh")
  print("should_transfer_files = YES")
  print("#when_to_transfer_output = ON_EXIT")
  print ""

def PrintTransferFileLine(directory, sampletype, isfirst, islast):
  if(isfirst):
    sys.stdout.write('transfer_input_files = $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SSSignalScan, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/NTuple_SSSignalScan.py, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/goSSSignalScan.sh, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/CSVv2_ichep.csv, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/TTbarNoHad_bTagEff.root, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/PileupHistograms_Nov17.root, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/Example_Legacy_TopTagger.cfg, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/Example_TopTagger.cfg, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/TrainingOutput.model, ')
  for dirname, dirnames, filenames in os.walk(directory):
    for filename in filenames:
      if ( sampletype in filename ):
        sys.stdout.write('$ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SensitivityTxt/' + filename + ', ')
      else:
        continue
  if(islast):
    print ""
    print ""

def PrintCondorLogLine():
  print ("Output = $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SensitivityTxt/res/Trim_$(Process).stdout")
  print ("Error = $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SensitivityTxt/res/Trim_$(Process).stderr")
  print ("Log = $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SensitivityTxt/res/Trim_$(Process).log")
  print ("notify_user = hua.wei@cern.ch")
  print ""

def PrintCondorSubmitLine(directory, sampletype):
  print("#### "+ sampletype +" ####")
  print ""
  for dirname, dirnames, filenames in os.walk(directory):
    for filename in filenames:
      if ( sampletype in filename ):
        print ("arguments = $ENV(CMSSW_BASE)  " + filename)
        print ("Queue")
        print ""
      else:
        continue

d = "/uscms_data/d3/hwei/stop/SS/CMSSW_8_0_23/src/SensitivityStudy/SSTreeMaker/SensitivityTxt"
runtype = sys.argv[1]
print ("#The valid run types for SS are Signal, Background! While the current run type is : " + runtype)

if(runtype == "Signal"):
  PrintCondorHeaderLine()
  print("##transfer file list for " + runtype + " samples")
  PrintTransferFileLine(d, "T1tttt", True, False)
  PrintTransferFileLine(d, "T2tt", False, False)
  PrintTransferFileLine(d, "T5ttcc", False, True)
  PrintCondorLogLine()
  PrintCondorSubmitLine(d, "T1tttt")
  PrintCondorSubmitLine(d, "T2tt")
  PrintCondorSubmitLine(d, "T5ttcc")

else:
  print ("#Invalid run type for SSSignalScan! What the fuck is going on ??!!")
