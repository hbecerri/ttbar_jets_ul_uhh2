#!/bin/bash
#where UHH2 code installed
pathGL_code=/nfs/dust/cms/user/hugobg/UHH2_UL/CMSSW_10_6_28/src/UHH2/
#where (NOT MERGED) trees after preselection stored
path_data=/nfs/dust/cms/user/hugobg/output_UL/Preselection_UL16postVFP/workdir_Preselection_UL16postVFP/

mkdir $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_2016UL
cd $pathGL_code/ZprimeSemiLeptonic/data/Skimming_datasets_2016UL

 #MC
#for sample_name in DYJetsToLL_M-50_HT-600to800_UL16postVFP DYJetsToLL_M-50_HT-600to800_UL16postVFP DYJetsToLL_M-50_HT-100to200_UL16postVFP DYJetsToLL_M-50_HT-1200to2500_UL16postVFP DYJetsToLL_M-50_HT-200to400_UL16postVFP DYJetsToLL_M-50_HT-2500toInf_UL16postVFP DYJetsToLL_M-50_HT-400to600_UL16postVFP DYJetsToLL_M-50_HT-70to100_UL16postVFP DYJetsToLL_M-50_HT-800to1200_UL16postVFP QCD_HT1000to1500_UL16postVFP QCD_HT100to200_UL16postVFP QCD_HT1500to2000_UL16postVFP QCD_HT2000toInf_UL16postVFP QCD_HT200to300_UL16postVFP QCD_HT300to500_UL16postVFP QCD_HT500to700_UL16postVFP QCD_HT50to100_UL16postVFP QCD_HT700to1000_UL16postVFP ST_s-channel_4f_leptonDecays_UL16postVFP ST_t-channel_antitop_4f_InclusiveDecays_UL16postVFP ST_t-channel_top_4f_InclusiveDecays_UL16postVFP ST_tW_antitop_5f_NoFullyHadronicDecays_UL16postVFP ST_tW_top_5f_NoFullyHadronicDecays_UL16postVFP TTTo2L2Nu_UL16postVFP TTToHadronic_UL16postVFP TTToSemiLeptonic_UL16postVFP WJetsToLNu_HT-100To200_UL16postVFP WJetsToLNu_HT-1200To2500_UL16postVFP WJetsToLNu_HT-200To400_UL16postVFP WJetsToLNu_HT-2500ToInf_UL16postVFP WJetsToLNu_HT-400To600_UL16postVFP WJetsToLNu_HT-600To800_UL16postVFP WJetsToLNu_HT-70To100_UL16postVFP WJetsToLNu_HT-800To1200_UL16postVFP WW_UL16postVFP WZ_UL16postVFP ZZ_UL16postVFP 

#do
#    echo $sample_name
#
#       $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"uhh2.AnalysisModuleRunner.MC."${sample_name}"*.root" MC_$sample_name.xml
#       python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 MC_$sample_name.xml True
#done

# # #DATA
for sample_name in DATA_SingleMuon_RunG DATA_SingleMuon_RunF # DATA_SingleMuon_RunG DATA_SingleMuon_RunH DATA_SingleMuon_RunF_UL16postVFP

do
    echo $sample_name
    $pathGL_code/scripts/create-dataset-xmlfile ${path_data}"uhh2.AnalysisModuleRunner.DATA."${sample_name}"*.root" DATA_$sample_name.xml
    python $pathGL_code/scripts/crab/readaMCatNloEntries.py 10 DATA_$sample_name.xml True

done
pwd
cd $pathGL_code/ZprimeSemiLeptonic/macros
