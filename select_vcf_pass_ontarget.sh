#!/bin/bash

vcf=$1
bedfile=$2
vcfout=${vcf/.vcf/.filt.vcf}

if [[ -z $vcf || -z $bedfile ]]; then
	echo -e "\nThis script takes a VCF and uses bcftools to create a VCF with\nvariants in the DERMATLAS bait regions BED file. The output\nfile is *.filt.vcf*\n"
	echo -e "\tUsage: $0 [VCF_file]\n"
	exit 0
elif [[ ! -f $bedfile ]]; then
	echo "File $bedfile does not exist."
	exit 1
fi

echo -e "Outfile $vcfout"

#bedfile=/lustre/scratch125/casm/team113da/projects/DERMATLAS/metadata/references/baitset/DNA/GRCh38_WES5_canonical_pad100.merged.bed

if which tabix &> /dev/null; then
    echo "tabix is available."
    tabix --version
else 
    echo "tabix is not available."
    exit 1
fi

if which bcftools &> /dev/null; then
    echo "bcftools is available."
    bcftools --version
else 
    echo "bcftools is not available."
    exit 1
fi


bcftools view -f PASS -O z -o $vcfout -T $bedfile $vcf
tabix -p vcf $vcfout
