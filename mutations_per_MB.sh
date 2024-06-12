#!/bin/bash

# This script calculates mutation rates from MAF files
# found in DERMALTAS directories 'all', 'independent'
# and 'onePerPatient'. 

echo -e "\nThis script calculates TMB from MAF files."
echo -e "Directories 'all', 'independent', 'onePerPatient' are expected in the current working directory."
echo -e "TMB will be calculated using MAFs found in these directories.\n\n"

if [[ ! -d "all" && ! -d "independent" && ! -d "onePerPatient" ]]; then
	echo -e "None of these directories were found. Exiting.\n\n"
	exit 1
fi

for h in all independent onePerPatient; do
	if [[ -d $h ]]; then
        cd $h
        echo $h
        for f in keep keep_matched keep_vaf_filt keep_vaf_filt_matched keep_vaf_size_filt keep_vaf_size_filt_matched; do
            for g in `cut -f 11 ${f}_*.maf |grep PD |sort -u`; do 
                echo -ne "$g\t"; muts=`grep $g ${f}_*maf | cut -f 4,5 |sort -u|wc -l`; echo $muts/48.225157 | bc -l;
            done > plots_${f}/mutations_per_Mb.tsv;
        done
        cd ..
        echo $h
    fi
done
