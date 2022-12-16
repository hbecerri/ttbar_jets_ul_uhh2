#!/bin/bash

declare -a StringArray=("M_tt M_th M_tl")

for var in ${StringArray[@]}; do

python makePlots_fromAnalysisTree.py ${var}

done


