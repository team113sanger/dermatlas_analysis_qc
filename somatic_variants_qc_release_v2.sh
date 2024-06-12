#!/bin/bash

vcflist=$1
maf_file=$2
PROJECTDIR=$3
which=$4

#RSCRIPT=/software/team113/dermatlas/R/R-4.2.2/bin/Rscript
#RSCRIPT2="module load R/4.1.3; Rscript"
#RSCRIPT2="/software/R-4.1.3/bin/Rscript"

if [[ -z "$PROJECTDIR" || -z $vcflist ]]; then
	echo -e "\nUsage: somatic_variants_qc.sh vcflist output_file /path/to/projectdir/ [filter1|filter2]\n"
	echo -e "where 'filter1' indicates that the calls should be filtered with VAF >= 0.1 and 2filter2' indicates filtering by REF and ALT length and 'ONP' variants (e.g indel calls)."
	echo -e "Leave the 4th option blank if no filtering by VAF is required."
	echo -e "Example: somatic_variants_qc.sh vcfs.list out.maf /lustre/scratch125/casm/team113da/projects/dermatlas_pu3_project_dir/6684_3062_DERMATLAS_Porocarcinoma_WES/ filter\n"

	exit 0
fi

if [[ ! -z $which && ! ($which == "filter1" || $which == "filter2") ]]; then
	echo "Leave the 4th option blank or use 'filter1' to indicate filtering by VAF >= 0.1 or 'filter2' to filter by 'ONP', VAF and REF and ALT lengths (indels)"
fi

# Check that directories and files exist

echo "Looking for script directory"

SCRIPTDIR=$PROJECTDIR/scripts/MAF

if [[ ! -d $SCRIPTDIR ]]; then
	echo "Directory $SCRIPTDIR does not exist"
	exit
elif [[ ! -e $SCRIPTDIR/reformat_vcf2maf.pl ]]; then
	echo "File $SCRIPTDIR/reformat_vcf2maf.pl does not exist"
	exit
elif [[ ! -e $SCRIPTDIR/plot_vaf_vs_depth_from_maf.v3.R ]]; then
	echo "File $SCRIPTDIR/plot_vaf_vs_depth_from_maf.v3.R does not exist"
	exit
elif [[ ! -e $SCRIPTDIR/maketileplot_from_maf.v4.R ]]; then
	echo "$SCRIPTDIR/maketileplot_from_maf.v4.R does not exist"
	exit
fi

echo "Output directory is current directory: $PWD"

# Reformat VCFs to MAF
#
#	Usage: reformatVCF2TSV.pl --vcflist [file] 
#	
#	where file is a file with a list of VCF files (full or relative paths). File formats
#	parsed are MuTect (v1), Strelka2, cgpCaVEMan and cgpPindel and may be a combination
#	of any of these. VCF type will be determined by header information.
#
#	Optional:
#
#		--af_col    Column name with population AFs to be used when filtering variants to 'keep'
#		            and  variants of interest (along with variant consequences. Population AF 
#		            cutoffs 0.01 (common SNPs) and 0.001 used for 'keep' and 'voi', respectively).
#		--build     NCBI Genome build [default: GRCh38]
#
#	Output options:
#		--pass      Print out PASS calls only (based only on PASS filter)
#		--voi_only  Print out variants of interest only (not compatible with --keep)
#		--keep      Print out 'keep' and 'keep-PA' variants only (not compatible with --voi_only or --keepPA)
#		--keepPA    Print out 'keep-PA' variants only (not compatible with --voi_only or --keep)
#
#		--tum_vaf   Print out variants with tumour VAF above this value [default: not used].

#		            Works with any of the 3 output options above.


# Make QC plots

#  Usage: Rscript plot.R [options]
#
#  Options:
#    --file, -f character  Variants file [required]
#    --genefile, -g character List of genes to include. 
#    --colname, -c character Column in variants file with gene name matching the gene list. Required with --genelist
#    --suffix, -s character Suffix to be added to output files when using a custom gene list. Required with --genelist
#    --noindels            Ingore indels (only SNVs and MNVs will be considered)
#    --indels_only, -i     Plot indels only
#    --protein_alt, -p     Exclude variants with main consequence synonymous, stop/start retained, splice_region
#
#  Options for AF vs Depth plots:
#    --width, -w numeric   Width (inches) of AF vs Depth by sample plots
#    --height numeric      Height (inches) of AF vs Depth by sample plots
#    --ncol integer        Number of columns for AF vs Depth plot
# 
#

$SCRIPTDIR/reformat_vcf2maf.pl --vcflist $vcflist --build GRCh38 --pass --af_col gnomAD_AF --sample_list sample_list.tsv > pass_${maf_file}
$SCRIPTDIR/reformat_vcf2maf.pl --vcflist $vcflist --build GRCh38 --voi_only --af_col gnomAD_AF > voi_${maf_file}

head -n1 pass_${maf_file} > keep_${maf_file}
awk 'BEGIN{FS="\t"}{if($28~/keep/){print}}' pass_${maf_file} >> keep_${maf_file}

head -n1 keep_${maf_file} > keepPA_${maf_file}
awk 'BEGIN{FS="\t"}{if($28~/keep-PA/){print}}' keep_${maf_file} >> keepPA_${maf_file}

# Remove tumours run with in silico normal

grep -v PDv38is_wes_v2 pass_${maf_file} > pass_matched_${maf_file}
grep -v PDv38is_wes_v2 keep_${maf_file} > keep_matched_${maf_file}
grep -v PDv38is_wes_v2 keepPA_${maf_file} > keepPA_matched_${maf_file}
grep -v PDv38is_wes_v2 voi_${maf_file} > voi_matched_${maf_file}

grep -v wes_v2 sample_list.tsv > sample_list_matched.tsv

cut -f 1 sample_list.tsv > sample_list_tum.tsv
cut -f 1 sample_list_matched.tsv > sample_list_matched_tum.tsv

# Filter VCFs
##--vcflist test.list  --build GRCh38  --keepPA --af_col gnomAD_AF --tum_vaf 0.1 --vaf_filter_type indel --indel_filter

if [[ $which == "filter1" || $which == "filter2" ]]; then
	$SCRIPTDIR/reformat_vcf2maf.pl --vcflist $vcflist --build GRCh38 --pass --af_col gnomAD_AF --tum_vaf 0.1 --vaf_filter_type indel > pass_vaf_filt_${maf_file}
	$SCRIPTDIR/reformat_vcf2maf.pl --vcflist $vcflist --build GRCh38 --voi_only --af_col gnomAD_AF --tum_vaf 0.1 --vaf_filter_type indel > voi_vaf_filt_${maf_file}
	head -n1 pass_vaf_filt_${maf_file} > keep_vaf_filt_${maf_file}
	awk 'BEGIN{FS="\t"}{if($28~/keep/){print}}' pass_vaf_filt_${maf_file} >> keep_vaf_filt_${maf_file}

	head -n1 keep_vaf_filt_${maf_file} > keepPA_vaf_filt_${maf_file}
	awk 'BEGIN{FS="\t"}{if($28~/keep-PA/){print}}' keep_vaf_filt_${maf_file} >> keepPA_vaf_filt_${maf_file}

	grep -v PDv38is_wes_v2 pass_vaf_filt_${maf_file} > pass_vaf_filt_matched_${maf_file}
	grep -v PDv38is_wes_v2 keep_vaf_filt_${maf_file} > keep_vaf_filt_matched_${maf_file}
	grep -v PDv38is_wes_v2 keepPA_vaf_filt_${maf_file} > keepPA_vaf_filt_matched_${maf_file}
	grep -v PDv38is_wes_v2 voi_vaf_filt_${maf_file} > voi_vaf_filt_matched_${maf_file}

fi

if [[ $which == "filter2" ]]; then
	$SCRIPTDIR/reformat_vcf2maf.pl --vcflist $vcflist --build GRCh38 --pass --af_col gnomAD_AF --tum_vaf 0.1 --vaf_filter_type indel --indel_filter > pass_vaf_size_filt_${maf_file}
	$SCRIPTDIR/reformat_vcf2maf.pl --vcflist $vcflist --build GRCh38 --voi_only --af_col gnomAD_AF --tum_vaf 0.1 --vaf_filter_type indel --indel_filter > voi_vaf_size_filt_${maf_file}

	head -n1 pass_vaf_size_filt_${maf_file} > keep_vaf_size_filt_${maf_file}
	awk 'BEGIN{FS="\t"}{if($28~/keep/){print}}' pass_vaf_size_filt_${maf_file} >> keep_vaf_size_filt_${maf_file}

	head -n1 keep_vaf_size_filt_${maf_file} > keepPA_vaf_size_filt_${maf_file}
	awk 'BEGIN{FS="\t"}{if($28~/keep-PA/){print}}' keep_vaf_size_filt_${maf_file} >> keepPA_vaf_size_filt_${maf_file}

	grep -v PDv38is_wes_v2 pass_vaf_size_filt_${maf_file} > pass_vaf_size_filt_matched_${maf_file}
	grep -v PDv38is_wes_v2 keep_vaf_size_filt_${maf_file} > keep_vaf_size_filt_matched_${maf_file}
	grep -v PDv38is_wes_v2 keepPA_vaf_size_filt_${maf_file} > keepPA_vaf_size_filt_matched_${maf_file}
	grep -v PDv38is_wes_v2 voi_vaf_size_filt_${maf_file} > voi_vaf_size_filt_matched_${maf_file}

fi

# Make plots

echo "Loading R/4.1.3"
module load R/4.1.3
which Rscript
RSCRIPT=Rscript

for type in keep keepPA voi; do
	mkdir -p plots_${type}
	cd plots_${type}
	$RSCRIPT $SCRIPTDIR/plot_vaf_vs_depth_from_maf.v3.R --file ../${type}_${maf_file} --width 8 --height 10 -ncol 7
	cut -f 3 top_recurrently_mutated_genes.tsv |sort -u | grep -v Hugo_ > top_genes.list
	echo "Plotting tile plot $type"
	$RSCRIPT $SCRIPTDIR/maketileplot_from_maf.v4.R -a ../${type}_${maf_file} -s ../sample_list_tum.tsv  -g top_genes.list --sortbyfrequency -w 8 -t 5 
	echo "$RSCRIPT $SCRIPTDIR/maketileplot_from_maf.v4.R -a ../${type}_${maf_file} -s ../sample_list_tum.tsv  -g top_genes.list --sortbyfrequency -w 8 -t 5"
	cd ..
done

for type in keep_matched keepPA_matched voi_matched; do
	mkdir -p plots_${type}
	cd plots_${type}
	$RSCRIPT $SCRIPTDIR/plot_vaf_vs_depth_from_maf.v3.R --file ../${type}_${maf_file} --width 8 --height 10 -ncol 7
	cut -f 3 top_recurrently_mutated_genes.tsv | sort -u | grep -v Hugo_ > top_genes.list
	echo "Plotting tile plot $type"
	$RSCRIPT $SCRIPTDIR/maketileplot_from_maf.v4.R -a ../${type}_${maf_file} -s ../sample_list_matched_tum.tsv -g top_genes.list --sortbyfrequency -w 8 -t 5
	echo "$RSCRIPT $SCRIPTDIR/maketileplot_from_maf.v4.R -a ../${type}_${maf_file} -s ../sample_list_matched_tum.tsv -g top_genes.list --sortbyfrequency -w 8 -t 5"
	cd ..
done

if [[ $which == "filter1" || $which == "filter2" ]]; then
	type_array=("keep_vaf_filt" "keepPA_vaf_filt" "voi_vaf_filt")
	if [[ $which == "filter2" ]]; then
		type_array+=("keep_vaf_size_filt" "keepPA_vaf_size_filt" "voi_vaf_size_filt")
	fi
	for type in ${type_array[*]}; do
		mkdir -p plots_${type}
		cd plots_${type}
		$RSCRIPT $SCRIPTDIR/plot_vaf_vs_depth_from_maf.v3.R --file ../${type}_${maf_file} --width 8 --height 10 -ncol 7
		cut -f 3 top_recurrently_mutated_genes.tsv |sort -u | grep -v Hugo_ > top_genes.list
		echo "Plotting tile plot $type"
		$RSCRIPT $SCRIPTDIR/maketileplot_from_maf.v4.R -a ../${type}_${maf_file} -s ../sample_list_tum.tsv  -g top_genes.list --sortbyfrequency -w 8 -t 5 
		echo "$RSCRIPT $SCRIPTDIR/maketileplot_from_maf.v4.R -a ../${type}_${maf_file} -s ../sample_list_tum.tsv  -g top_genes.list --sortbyfrequency -w 8 -t 5"
		cd ..
	done

	type_matched_array=("keep_vaf_filt_matched" "keepPA_vaf_filt_matched" "voi_vaf_filt_matched")
	if [[ $which == "filter2" ]]; then
		type_matched_array+=("keep_vaf_size_filt_matched" "keepPA_vaf_size_filt_matched" "voi_vaf_size_filt_matched")
	fi
	for type in ${type_matched_array[*]}; do
		mkdir -p plots_${type}
		cd plots_${type}
		$RSCRIPT $SCRIPTDIR/plot_vaf_vs_depth_from_maf.v3.R --file ../${type}_${maf_file} --width 8 --height 10 -ncol 7
		cut -f 3 top_recurrently_mutated_genes.tsv | sort -u | grep -v Hugo_ > top_genes.list
		echo "Plotting tile plot $type"
		$RSCRIPT $SCRIPTDIR/maketileplot_from_maf.v4.R -a ../${type}_${maf_file} -s ../sample_list_matched_tum.tsv -g top_genes.list --sortbyfrequency -w 8 -t 5
		echo "$RSCRIPT $SCRIPTDIR/maketileplot_from_maf.v4.R -a ../${type}_${maf_file} -s ../sample_list_matched_tum.tsv -g top_genes.list --sortbyfrequency -w 8 -t 5"
		cd ..
	done
fi

