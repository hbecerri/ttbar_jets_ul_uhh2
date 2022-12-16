#!/bin/bash

declare -a StringArray=("M_th" "M_tl" "M_tt")

for var in ${StringArray[@]}; do

rm *.pdf 
python makePlots_fromAnalysisTree.py ${var}

done


