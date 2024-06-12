#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $infile;
my $prefix;
my $rejected;

my $check = GetOptions ("infile=s" => \$infile,
            "prefix=s" => \$prefix,
            "reject=s" => \$rejected);

my $help = <<END;
	
	Usage: $0 --infile biosample_manifest.tsv --prefix \${STUDY}_\${PROJ} --reject rejected_samples.list

	where 
	--infile is a completed biosample manifest with samples selected [Required]
	--prefix is a prefix for the output files [Suggested: Sequencescape studyID and CanApps ID] [Required]
	--reject is a list of samples to exclude/ignore [Optional]

END

if (!$check) {
	print STDERR $help;
	exit;
}

# Read in a completed biosample manifest and print out t/n pairs,

if (!$infile || !$prefix) {
	print STDERR $help;
	exit;
}

my %sample;
my %reject;

if ($rejected) {
	print STDERR "Reading list of rejected samples\n";
	open R, "<$rejected" or die  "Cannot open $rejected\n";
	while (<R>) {
		if (/(\S+)/) {
			$reject{$1} = 1;
			print STDERR "$1\n";
		}
	}
	close R;
}

open ALL, ">$prefix-analysed_all.tsv" or die "Cannot write to file $prefix-analysed_all.tsv\n";
open ALL_UNMATCHED, ">$prefix-analysed_unmatched.tsv" or die "Cannot write to file $prefix-analysed_samples_unmatched.tsv\n";
open ALL_MATCHED, ">$prefix-analysed_matched.tsv" or die "Cannot write to file $prefix-analysed_samples_matched.tsv\n";
open F, "<$infile" or die "Can't open $infile\n";

# Get the column numbers for specific columns

my $tum_col=`awk -v RS='\t' '/Tumour Sample\$/{print NR; exit}' $infile`;
my $norm_col=`awk -v RS='\t' '/Normal Sample\$/{print NR; exit}' $infile`;
my $submit_col=`awk -v RS='\t' '/Submit to/{print NR; exit}' $infile`;
my $analyse_col=`awk -v RS='\t' '/Use for analysis/{print NR; exit}' $infile`;

chomp($tum_col, $norm_col, $submit_col, $analyse_col);

# Column number to array index

$tum_col--;
$norm_col--;
$submit_col--;
$analyse_col--;

foreach my $header ($tum_col, $norm_col, $submit_col, $analyse_col) {
	if (!$header || ($header && $header < 0)) {
		print STDERR "\nOne or more columns not present: Tumour sample, Normal sample, Submit to Canpipe, Use for analysis.\n\n";
		exit;
	}
}

print STDERR "Tumour column: $tum_col, Normal column: $norm_col, Submit column: $submit_col, Analyse column: $analyse_col\n";

# Iterate through the biosample manifest

while (<F>) {
	if (/Biosample/) {
		next;
	}
	chomp;
	my $line = $_;
	my @col = split(/\t/, $line);
	my $norm = $col[$norm_col];
	my $tum = $col[$tum_col];

	print STDERR "TUM NORM: $tum $norm\n";

	# Check list of rejected samples
	if ($reject{$tum} && $reject{$tum} == 1) {
		print STDERR "Skipping $tum\n";
		next;
	} elsif ($reject{$norm} && $reject{$norm} == 1 && !$reject{$tum}) {
		print STDERR "Warning: $norm is on the reject list and $tum is not. Not excluding\n";
	}

	# Check 'use for analysis' column
	if ($col[$analyse_col] =~ /Y|y/) {
		my $patient = $1 if $tum =~ /(\S+)\S$/;
		$sample{$patient}{$tum}{line} = $line;
		$sample{$patient}{$tum}{norm} = $norm;
	}
	# Check 'Submit for analysis' column; all samples, including duplicates
	if ($col[$submit_col] =~ /Y|y/) {
		my $fh;
		if ($norm eq 'PDv38is_wes_v2') {
			#print ALL_UNMATCHED "$tum\t$norm\n";
			$fh = *ALL_UNMATCHED;
		} else {
			#print ALL_MATCHED "$tum\t$norm\n";
			$fh = *ALL_MATCHED;
		}
		print $fh "$tum\t$norm\n";
		print ALL "$tum\t$norm\n";
	}
}
close F;
close ALL;
close ALL_UNMATCHED;
close ALL_MATCHED;

# Print out all independent and one-per-patient tumours

open ONE, ">$prefix-one_tumour_per_patient_all.tsv" or die "Cannot write to file $prefix-one_tumour_per_patient_all.tsv\n";
open INDEP, ">$prefix-independent_tumours_all.tsv" or die "Cannot write to file $prefix-independent_tumours_all.tsv\n";
open ONE_UNMATCHED, ">$prefix-one_tumour_per_patient_unmatched.tsv" or die "Cannot write to file $prefix-one_tumour_per_patient_unmatched.tsv\n";
open INDEP_UNMATCHED, ">$prefix-independent_tumours_unmatched.tsv" or die "Cannot write to file $prefix-independent_tumours_unmatched.tsv\n";
open ONE_MATCHED, ">$prefix-one_tumour_per_patient_matched.tsv" or die "Cannot write to file $prefix-one_tumour_per_patient_matched.tsv\n";
open INDEP_MATCHED, ">$prefix-independent_tumours_matched.tsv" or die "Cannot write to file $prefix-independent_tumours_matched.tsv\n";

foreach my $patient (sort keys %sample) {
	my @tums = (sort keys %{$sample{$patient}});
	# Print one sample per patient, sort by alphabetical order
	my $norm = $sample{$patient}{$tums[0]}{norm};
	my $fh1;
	if ($norm eq 'PDv38is_wes_v2') {
		#print ONE_UNMATCHED "$tums[0]\t$sample{$patient}{$tums[0]}{norm}\n";
		$fh1 = *ONE_UNMATCHED;
	} else {
		#print ONE_MATCHED "$tums[0]\t$sample{$patient}{$tums[0]}{norm}\n";
		$fh1 = *ONE_MATCHED;
	}
	print $fh1 "$tums[0]\t$sample{$patient}{$tums[0]}{norm}\n";
	print ONE "$tums[0]\t$sample{$patient}{$tums[0]}{norm}\n";
	##print $sample{$patient}{$tums[0]}[0] . "\n";
	foreach my $i (1..$#tums) {
		print STDERR "Excluding lesion from same patient: $tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n$sample{$patient}{$tums[$i]}{line}\n";
	}

	# Print all independent tumours
	foreach my $i (0..$#tums) {
		my $fh2;
		my $norm = $sample{$patient}{$tums[$i]}{norm};
		if ($norm eq 'PDv38is_wes_v2') {
			#print INDEP_UNMATCHED "$tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n";
			$fh2 = *INDEP_UNMATCHED;
		} else {
			#print INDEP_MATCHED "$tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n";
			$fh2 = *INDEP_MATCHED;
		}
		print $fh2 "$tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n";
		print INDEP "$tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n";
	}
}

close ONE;
close INDEP;
close ONE_UNMATCHED;
close INDEP_UNMATCHED;
close ONE_MATCHED;
close INDEP_MATCHED;

# Create files with just the tumours

foreach my $file ('one_tumour_per_patient', 'independent_tumours', 'analysed') {
	foreach my $type ('matched', 'unmatched', 'all') {
		`cut -f 1 $prefix-${file}_${type}.tsv | sort > $prefix-${file}_${type}_tum.txt`;
	}
}

# Create a readme file explaining the files

open READ, ">readme_sample_lists.txt" or die "Cannot write to readme_sample_lists.txt\n";

print READ "The following files were derived from parsing file $infile\n";
if ($rejected) {
	print READ "\nSamples from $rejected were excluded:\n\n";
	print READ `cat $rejected`;
}

print READ <<END1;

********************************************************************
NOTE: THESE LISTS ARE AUTOMATICALLY GENERATED. NOT EVERY COHORT
WILL HAVE DUPLICATES OR MULTIPLE INDEPENDENT TUMOURS. PLEASE COMPARE
SAMPLE LISTS BEFLOW TO DETERMINE WHAT IS REQUIRED FOR YOUR ANALYSES.
********************************************************************

END1

print READ "\nFILES:\n\n";

foreach my $file ('one_tumour_per_patient', 'independent_tumours', 'analysed') {
	foreach my $type ('matched', 'unmatched', 'all') {
		print READ "$prefix-${file}_${type}.tsv\n";
	}
	foreach my $type ('matched', 'unmatched', 'all') {
		print READ "$prefix-${file}_${type}_tum.txt\n";
	}
}

print READ <<END2;

----------------------
File naming convention
----------------------

one_tumour_per_patient
----------------------

One tumour per patient is selected fromi the QC-pass/non-rejected
list, based on alphabetical order

independent_tumours
-------------------

All tumours from each patient from the QC-pass/non-rejected list,
annotated as independent tumours

analysed_tumours
----------------

All tumours from each patient from the QC-pass/non-rejected list,
including samples from the same tumour, that were run through
Caveman/Pindel

matched
-------

Tumour samples run with matched normals in Caveman/Pindel

unmatched
---------

Unmatched tumour samples run with in silico normal in Caveman/Pindel

all
---

Tumours with and without matched normal

*_tum.txt
---------

Same as above, but contains a single column with tumour sample
name only




END2

close READ;

