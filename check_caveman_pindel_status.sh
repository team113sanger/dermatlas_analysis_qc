#!/bin/bash

PROJECTDIR=$1
projID=$2
getfiles=$3

if [[ -z "$PROJECTDIR" || -z "$projID" ]]; then
	echo -e "\nUsage: check_caveman_pindel_status.sh /path/to/projectdir/ project_ID [getfiles]\n"
	echo -e "\twhere project_ID is the Canapps project number and getfiles is a flag to rsync vcfs from nst_links.\n"
	echo -e "Example: bash check_caveman_pindel_status.sh /lustre/scratch125/casm/team113da/projects/dermatlas_pu3_project_dir/6684_3062_DERMATLAS_Porocarcinoma_WES/ 3062 getfiles\n"
	exit 0
fi

# Check that directories and files exist

if [[ -d $PROJECTDIR ]]; then
	cd $PROJECTDIR
	if [[ -d "./metadata" ]]; then
		cd metadata
	else
		echo "Directory $PROJECTDIR/metadata does not exist"
		exit
	fi
else
	echo "Directory $PROJECTDIR does not exist"
	exit
fi


# Check for submitted samples from the biosample manifest

manifest=( `dir -1 *manifest*completed* | grep -v xls` )

if [[ ${#manifest[@]} > 1 ]]; then
	echo "Multiple manifests found: ${manifest[@]}"
	exit
elif [[ ${#manifest[@]} == 0 ]]; then
	echo "No manifest file found in $PROJECTDIR/metadata"
	exit
else 
	manifest=${manifest[0]}
fi

# Check for nst_links directory

if [[ ! -d /nfs/cancer_ref01/nst_links/live/$projID ]]; then
	echo "Directory /nfs/cancer_ref01/nst_links/live/$projID does not exist"
	exit
fi

# Check for pindel and Caveman outputs:

pindel=0
caveman=0
totsamples=0

echo "##############################"
echo -e "Checking manifest $PROJECTDIR/metadata/$manifest for submitted samples\n"

for sample in `awk 'BEGIN{FS="\t"}{if($8=="Y"){print $4}}' $manifest`; do
	let "totsamples=totsamples+1"
	pindelfile=/nfs/cancer_ref01/nst_links/live/$projID/$sample/$sample.pindel.flagged.vcf.gz
	if [[ ! -e $pindelfile ]]; then
		echo "$sample Pindel not complete"
	else
		let "pindel=pindel+1"
	fi
	cavemanfile=/nfs/cancer_ref01/nst_links/live/$projID/$sample/$sample.smartphase.vep.vcf.gz
	if [[ ! -e $cavemanfile ]]; then
		echo "$sample Caveman/smartphase/vep not complete"
	else
		let "caveman=caveman+1"
	fi
done

echo "Found $caveman of $totsamples *smartphase.vep.vcf.gz files in /nfs/cancer_ref01/nst_links/live/$projID/"
echo -e "Found $pindel of $totsamples *pindel.flagged.vcf.gz files in /nfs/cancer_ref01/nst_links/live/$projID/\n"

# Copy files

pindel=0
caveman=0
totsamples=0

if [[ ! -z "$getfiles" ]]; then

	# Check if output directory exists and check for files
	caveman_dir=$PROJECTDIR/analysis/caveman_files
	pindel_dir=$PROJECTDIR/analysis/pindel_files

	echo "##############################"
	echo -e "Copying Caveman and Pindel VCFs to $caveman_dir and $pindel_dir\n"
	if [[ -d $caveman_dir ]]; then
		echo "$caveman_dir already exists. Exiting."
		exit	
	else
		mkdir $caveman_dir
	fi
	if [[ -d $pindel_dir ]]; then
		echo "$pindel_dir already exists. Exiting."
		exit	
	else
		mkdir $pindel_dir
	fi

	for sample in `awk 'BEGIN{FS="\t"}{if($8=="Y"){print $4}}' $manifest`; do
		let "totsamples=totsamples+1"
		if [[ -e $pindelfile ]]; then
			rsync -avL /nfs/cancer_ref01/nst_links/live/$projID/$sample/*pindel*vcf* $pindel_dir/$sample/
			let "pindel=pindel+1"
		fi
		if [[ -e $cavemanfile ]]; then
			rsync -avL /nfs/cancer_ref01/nst_links/live/$projID/$sample/*caveman*vcf* $caveman_dir/$sample/
			rsync -avL /nfs/cancer_ref01/nst_links/live/$projID/$sample/*smartphase*vcf* $caveman_dir/$sample/
			let "caveman=caveman+1"
		fi
	done
	echo -e "\nrsync'ed Caveman output files for $caveman of $totsamples samples in /nfs/cancer_ref01/nst_links/live/$projID/"
	echo -e "rsync'ed Pindel output files for $pindel of $totsamples sample in /nfs/cancer_ref01/nst_links/live/$projID/\n"
fi
			



