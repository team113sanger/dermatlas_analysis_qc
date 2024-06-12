#!/bin/bash


usage () {
	echo -e "\nAdd dbSNP155 common SNP tag to a VCF INFO column.\n\nUsage $0: [options]
	Required options:
	-p   Full path to project directory 
	-o   Relative or absolute path to output directory
	-v   VCF input file
	'snpflagged' will be added to the VCF file name, e.g. vcfin.snpflagged.vcf.gz"
	exit 1;
}


while getopts ":p:o:v:h" flag; do
	case "${flag}" in
		p) PROJECTDIR=$OPTARG
		   ;;
		o) OUTDIR=$OPTARG
		   ;;
		v) VCF_IN=$OPTARG
		   ;;
		: )
		   echo "Invalid option: $OPTARG requires an argument" 1>&2
		   ;;	
		h | *) 
		   usage
		   exit 0
		   ;;
    esac
done
shift $((OPTIND -1))

# Check the output and project directories

if [[ -z "$PROJECTDIR" || -z $OUTDIR ]]; then
	usage
	exit 0
fi


if [[ ! -d $OUTDIR ]]; then
	echo "Directory $OUTDIR does not exist"
	exit
fi
#else
#	cd $OUTDIR
#	echo "Current directory is $OUTDIR"
#fi


# Check that the resources directories and files exist

echo "Looking for dbsnp resources directory"

SNPDIR=$PROJECTDIR/resources/dbsnp/

if [[ ! -d $SNPDIR ]]; then
	echo "Directory $SNPDIR does not exist"
	exit
elif [[ ! -e $SNPDIR/dbSNP155_common.tsv.gz ]]; then
	echo "File $SNPDIR/dbSNP155_common.tsv does not exist"
	exit
elif [[ ! -e $SNPDIR/addheader.txt ]]; then
	echo "$SNPDIR/addheader.txt does not exist"
	exit
fi

echo "Found dbsnp directory"

# Check for bcftools

if which bcftools &> /dev/null; then
	echo "Found bcftools"
    bcftools --version
else 
    echo "bcftools is not available."
    exit 1
fi

# Generate the output file name

vcf_prefix=`basename $VCF_IN .vcf.gz`
VCF_OUT=$OUTDIR/$vcf_prefix.snpflagged.vcf.gz

echo "Writing to file $VCF_OUT"

# Run bcftools annotate to add COMMON_SNP tag to the INFO column in a VCF

bcftools annotate -a $PROJECTDIR/resources/dbsnp/dbSNP155_common.tsv.gz -h $PROJECTDIR/resources/dbsnp/addheader.txt -c CHROM,POS,-,REF,ALT,dbSNP -O z -o $VCF_OUT $VCF_IN
tabix -p vcf $VCF_OUT




