import sys
import glob, os

def PrintCondorHeaderLine():
  print("universe = vanilla")
  print("request_disk   = 50 GB")
  print("request_memory = 1.5 GB")
  print("executable = $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/goSSTrimAndSlim.sh")
  print("should_transfer_files = YES")
  print("#when_to_transfer_output = ON_EXIT")
  print ""

def PrintTransferFileLine(directory, sampletype, isfirst, islast):
  if(isfirst):
    sys.stdout.write('transfer_input_files = $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SSTrimAndSlim_LLHadTau, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SSTrimAndSlim_ZinvQCD, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SSTrimAndSlim_TTZ, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SSTrimAndSlim_Signal, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/NTuple_SSTrimAndSlim.py, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/goSSTrimAndSlim.sh, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/CSVv2_ichep.csv, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/TTbarNoHad_bTagEff.root, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/PileupHistograms_Nov17.root, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/Example_Legacy_TopTagger.cfg, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/Example_TopTagger.cfg, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/TrainingOutput.model, ')
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

def PrintCondorSubmitLine(directory, runoption, sampletype):
  print("#### "+ sampletype +" ####")
  print ""
  for dirname, dirnames, filenames in os.walk(directory):
    for filename in filenames:
      if ( sampletype in filename ):
        print ("arguments = $ENV(CMSSW_BASE) " + runoption + " " + filename)
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
  PrintCondorSubmitLine(d, "Signal", "T1tttt")
  PrintCondorSubmitLine(d, "Signal", "T2tt")
  PrintCondorSubmitLine(d, "Signal", "T5ttcc")

elif(runtype == "Background"):
  PrintCondorHeaderLine()
  print("##transfer file list for " + runtype + " samples")
  PrintTransferFileLine(d, "TTJets_", True, False)
  PrintTransferFileLine(d, "WJetsToLNu_HT-", False, False)
  PrintTransferFileLine(d, "ST_tW_", False, False)
  PrintTransferFileLine(d, "ZJetsToNuNu_HT-", False, False)
  PrintTransferFileLine(d, "QCD_HT", False, False)
  PrintTransferFileLine(d, "TTWJets", False, False)
  PrintTransferFileLine(d, "TTZ", False, True)
  #PrintTransferFileLine(d, "WJetsToLNu_HT-200To400", True, False)
  #PrintTransferFileLine(d, "ZJetsToNuNu_HT-200To400", False, False)
  #PrintTransferFileLine(d, "QCD_HT300to500", False, True)

  PrintCondorLogLine()
  PrintCondorSubmitLine(d, "LLHadTau", "TTJets_")
  PrintCondorSubmitLine(d, "LLHadTau", "WJetsToLNu_HT-")
  PrintCondorSubmitLine(d, "LLHadTau", "ST_tW_")
  PrintCondorSubmitLine(d, "ZinvQCD", "ZJetsToNuNu_HT-")
  PrintCondorSubmitLine(d, "ZinvQCD", "QCD_HT")
  PrintCondorSubmitLine(d, "TTZ", "TTWJets")
  PrintCondorSubmitLine(d, "TTZ", "TTZ")
  #PrintCondorSubmitLine(d, "LLHadTau", "WJetsToLNu_HT-200To400")
  #PrintCondorSubmitLine(d, "ZinvQCD", "ZJetsToNuNu_HT-200To400")
  #PrintCondorSubmitLine(d, "ZinvQCD", "QCD_HT300to500")
else:
  print ("#Invalid run type for SSTrimAndSlim! What the fuck is going on ??!!")
