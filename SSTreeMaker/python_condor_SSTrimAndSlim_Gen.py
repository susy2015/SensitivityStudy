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

def PrintTransferFileLine(sampletype, isfirst, islast):
  if(isfirst):
    sys.stdout.write('transfer_input_files = $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/$ENV(CMSSW_VERSION).tar.gz, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SSTrimAndSlim_LLHadTau, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SSTrimAndSlim_ZinvQCD, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SSTrimAndSlim_TTZRare, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SSTrimAndSlim_Signal, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/NTuple_SSTrimAndSlim.py, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/goSSTrimAndSlim.sh, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/allINone_ISRJets.root, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/ISRWeights.root, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/CSVv2_Moriond17_B_H.csv, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/allINone_bTagEff.root, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/allINone_leptonSF_Moriond17.root, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/PileupHistograms_0121_69p2mb_pm4p6.root, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/puppiCorr.root, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/TopTagger.cfg, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/TrainingOutput_dR20_pt30_depth12_500tree_noQGL_binaryCSV_2017_Mar24.model, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/Legacy_TopTagger.cfg, $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SensitivityTxt.tar.gz')
  if(islast):
    print ""
    print ""

def PrintCondorLogLine(runtype):
  print ("Output = $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SensitivityTxt/res/SSTrimAndSlim_" + runtype + "_$(Process).stdout")
  print ("Error = $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SensitivityTxt/res/SSTrimAndSlim_" + runtype + "_$(Process).stderr")
  print ("Log = $ENV(CMSSW_BASE)/src/SensitivityStudy/SSTreeMaker/SensitivityTxt/res/SSTrimAndSlim_" + runtype + "_$(Process).log")
  print ("notify_user = hua.wei@cern.ch")
  print ""

def PrintCondorSubmitLine(directory, sampletype, command):
  print("#### "+ sampletype +" ####")
  print ""
  for dirname, dirnames, filenames in os.walk(directory):
    for filename in filenames:
      if ( sampletype in filename ):
        print ("arguments = $ENV(CMSSW_VERSION) $ENV(SCRAM_ARCH) " + command + " SensitivityTxt/" + filename)
        print ("Queue")
        print ""
      else:
        continue

d = os.environ.get('CMSSW_BASE') + "/src/SensitivityStudy/SSTreeMaker/SensitivityTxt"
runtype = sys.argv[1]
print ("#The valid run types for SS are Signal, Background! While the current run type is : " + runtype)

if(runtype == "Signal"):
  PrintCondorHeaderLine()
  print("##transfer file list for " + runtype + " samples")
  PrintTransferFileLine("T1tttt",  True, False)
  PrintTransferFileLine("T1ttbb", False, False)
  PrintTransferFileLine("T2tt"  , False, False)
  PrintTransferFileLine("T5tttt", False, False)
  PrintTransferFileLine("T5ttcc", False,  True)
  PrintCondorLogLine(runtype)
  PrintCondorSubmitLine(d, "T1tttt", runtype)
  PrintCondorSubmitLine(d, "T1ttbb", runtype)
  PrintCondorSubmitLine(d, "T2tt"  , runtype)
  PrintCondorSubmitLine(d, "T5tttt", runtype)
  PrintCondorSubmitLine(d, "T5ttcc", runtype)
elif(runtype == "Background"):
  PrintCondorHeaderLine()
  print("##transfer file list for " + runtype + " samples")
  PrintTransferFileLine("TTJets_"        , True , False)
  PrintTransferFileLine("WJetsToLNu_HT-" , False, False)
  PrintTransferFileLine("ST_tW_"         , False, False)
  PrintTransferFileLine("ZJetsToNuNu_HT-", False, False)
  PrintTransferFileLine("QCD_HT"         , False, False)
  PrintTransferFileLine("TTWJets"        , False, False)
  PrintTransferFileLine("TTZ"            , False, False)
  PrintTransferFileLine("WW", False, False)
  PrintTransferFileLine("WZ", False, False)
  PrintTransferFileLine("ZZ", False, True)
  PrintCondorLogLine(runtype)
  PrintCondorSubmitLine(d, "TTJets_"        , "LLHadTau")
  PrintCondorSubmitLine(d, "WJetsToLNu_HT-" , "LLHadTau")
  PrintCondorSubmitLine(d, "ST_tW_"         , "LLHadTau")
  PrintCondorSubmitLine(d, "ZJetsToNuNu_HT-", "ZinvQCD" )
  PrintCondorSubmitLine(d, "QCD_HT"         , "ZinvQCD" )
  PrintCondorSubmitLine(d, "TTWJets"        , "TTZRare" )
  PrintCondorSubmitLine(d, "TTZ"            , "TTZRare" )
  PrintCondorSubmitLine(d, "WW", "TTZRare")   
  PrintCondorSubmitLine(d, "WZ", "TTZRare")   
  PrintCondorSubmitLine(d, "ZZ", "TTZRare")   
else:
  print ("#Invalid run type for SSTrimAndSlim! What the fuck is going on ??!!")
