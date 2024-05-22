#!/bin/bash

vcf=$1
vcfout=${vcf/.vcf/.filt.vcf}

if [ -z "$vcf" ]; then
	echo "Input a VCF file"
	exit
fi

echo -e "Outfile $vcfout"

bedfile=/lustre/scratch125/casm/team113da/projects/DERMATLAS/metadata/references/baitset/DNA/GRCh38_WES5_canonical_pad100.merged.bed

module load bcftools/1.9 tabix/1.9

bcftools view -f PASS -O z -o $vcfout -T $bedfile $vcf
tabix -p vcf $vcfout
