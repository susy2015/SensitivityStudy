import glob, os, sys

#d= "root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/Skimmed_2015Nov15"
d_hadd_source = "/store/group/lpcsusyhad/hua/Skimmed_2015Nov15"
d_final_target = "/store/group/lpcsusyhad/hua/Skimmed_2015Nov15/Sensitivity_MC_v6"
MGM = "root://cmseos.fnal.gov/"
#example of sample tag: SSTrimmed_SMS-T2tt_mStop-850_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8stopFlatNtuples_, always end with "_"
sample_tag = sys.argv[1]
#rootfile_tag = file_name[:-4]
#print(rootfile_tag) 

#hadd myTarget.root `xrdfsls -u /store/group/lpcsusyhad/hua/Skimmed_2015Nov15 | grep 'SSTrimmed_SMS-T2tt_mStop-850_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8stopFlatNtuples_'`
#hadd
cmd = 'hadd ' + sample_tag[:-1] + '.root `xrdfsls -u ' + d_hadd_source + ' | grep \'' + sample_tag + '\'`'
print(cmd)
#os.system(cmd)
cmd = 'xrdcp ' + sample_tag[:-1] + '.root ' + MGM + d_final_target
print(cmd)

#xrdcp root://cmseos.fnal.gov//store/user/jjesus/EOSFile.txt \? root://cmseos.fnal.gov//store/user/jjesus/EOSFile1.txt
cmd = 'xrdcp ' + MGM + d_hadd_source + '/' + sample_tag[:-1] + '.root ' + MGM + d_final_target + '/' + sample_tag[:-1] + '.root' 
print(cmd)
cmd = 'eosrm ' + d_hadd_source + '/' + sample_tag[:-1] + '.root'
print(cmd)

