#!/bin/bash
#where UHH2 code installed
pathGL_code=/nfs/dust/cms/user/deleokse/analysis/CMSSW_10_2_10/src/UHH2/
#where (NOT MERGED) trees after preselection stored
path_data=/nfs/dust/cms/user/deleokse/analysis/ZPrime_2017/Analysis_2017v2_2502_ApplyDNN/muon/uhh2.AnalysisModuleRunner. 

mkdir $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_2017v2_AfterSelection_2502
cd $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_2017v2_AfterSelection_2502


# #MC

for sample_name in TTToSemiLeptonic_2017v2 TTToHadronic_2017v2 TTTo2L2Nu_2017v2 WW_2017v2 WZ_2017v2 ZZ_2017v2 ST_t-channel_antitop_2017v2 ST_t-channel_top_2017v2 ST_tW_antitop_5f_inclusiveDecays_2017v2 ST_tW_top_5f_inclusiveDecays_2017v2 DYJetsToLL_M-50_HT-100to200_2017v2 DYJetsToLL_M-50_HT-1200to2500_2017v2 DYJetsToLL_M-50_HT-200to400_2017v2 DYJetsToLL_M-50_HT-2500toInf_2017v2 DYJetsToLL_M-50_HT-400to600_2017v2 DYJetsToLL_M-50_HT-600to800_2017v2 DYJetsToLL_M-50_HT-800to1200_2017v2 WJetsToLNu_HT-70To100_2017v2 WJetsToLNu_HT-100To200_2017v2 WJetsToLNu_HT-200To400_2017v2 WJetsToLNu_HT-400To600_2017v2 WJetsToLNu_HT-600To800_2017v2 WJetsToLNu_HT-800To1200_2017v2 WJetsToLNu_HT-1200To2500_2017v2 WJetsToLNu_HT-2500ToInf_2017v2 WJetsToLNu_HT-70To100_2017v2_B WJetsToLNu_HT-100To200_2017v2_B WJetsToLNu_HT-200To400_2017v2_B WJetsToLNu_HT-400To600_2017v2_B WJetsToLNu_HT-600To800_2017v2_B WJetsToLNu_HT-800To1200_2017v2_B WJetsToLNu_HT-1200To2500_2017v2_B WJetsToLNu_HT-2500ToInf_2017v2_B WJetsToLNu_HT-70To100_2017v2_C WJetsToLNu_HT-100To200_2017v2_C WJetsToLNu_HT-200To400_2017v2_C WJetsToLNu_HT-400To600_2017v2_C WJetsToLNu_HT-600To800_2017v2_C WJetsToLNu_HT-800To1200_2017v2_C WJetsToLNu_HT-1200To2500_2017v2_C WJetsToLNu_HT-2500ToInf_2017v2_C WJetsToLNu_HT-70To100_2017v2_L WJetsToLNu_HT-100To200_2017v2_L WJetsToLNu_HT-200To400_2017v2_L WJetsToLNu_HT-400To600_2017v2_L WJetsToLNu_HT-600To800_2017v2_L WJetsToLNu_HT-800To1200_2017v2_L WJetsToLNu_HT-1200To2500_2017v2_L WJetsToLNu_HT-2500ToInf_2017v2_L  

do
    echo $sample_name

       $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"MC."${sample_name}"*.root" MC_$sample_name.xml
       python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 MC_$sample_name.xml True
done

## # #DATA
#for sample_name in DATA_SingleMuon_Run2017v2_C DATA_SingleMuon_Run2017v2_D DATA_SingleMuon_Run2017v2_E DATA_SingleMuon_Run2017v2_F DATA_SingleElectron_Run2017v2_C DATA_SingleElectron_Run2017v2_D DATA_SingleElectron_Run2017v2_E DATA_SingleElectron_Run2017v2_F 
#
#do
#    echo $sample_name 
#    $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"DATA."${sample_name}"*.root" DATA_$sample_name.xml
#    python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 DATA_$sample_name.xml True
#
#done
pwd
cd $pathGL_code/ZprimeSemiLeptonic/macros