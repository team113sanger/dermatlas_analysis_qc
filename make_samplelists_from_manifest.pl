#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $infile;
my $prefix;
my $metadata_file;
my $rejected;

my $check = GetOptions ("infile=s" => \$infile,
            "prefix=s" => \$prefix,
            "metadata=s" => \$metadata_file,
            "reject=s" => \$rejected);

my $help = <<END;
	
	Usage: $0 --infile biosample_manifest.tsv --prefix \${STUDY}_\${PROJ} --reject rejected_samples.list --metadata metadata.tsv

	where 
	--infile is a completed biosample manifest with samples selected (TSV) [Required]
	--prefix is a prefix for the output files [Suggested: Sequencescape studyID and CanApps ID] [Required]
	--metadata is a metadata file with a 'Tumour type' column (TSV) [optional]
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

my %valid_types = map { $_ => 1 } ("P", "R", "M", "U", "-");
my %tum_types;

# Make a lookup table for tumour type
if (defined($metadata_file)) {
	my $tum_type_col=`awk -v RS='\t' '/Tumour type/{print NR; exit}' $metadata_file`;
	my $sanger_id_col=`awk -v RS='\t' '/Sanger_DNA_ID/{print NR; exit}' $metadata_file`;
	# Column number to array index
	$tum_type_col--;
	$sanger_id_col--;


	open MET, "<$metadata_file" or die "Can't open $metadata_file\n";
	while (<MET>) {
		if (/Sanger/) {
			next;
		}
		chomp;
		my $line = $_;
		my @col = split(/\t/, $line);
		my $sanger_id = $col[$sanger_id_col];
		my $type = $col[$tum_type_col];

		if (exists $valid_types{$type}) {
			$tum_types{$sanger_id} = $type;
		} else {
			die("Error: found tumour type $type for sample $sanger_id; must be one of " . join(",", keys(%valid_types)) . "\n");
		}
	}
	close MET;
}

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

		if (defined($metadata_file)) {
			if (exists($tum_types{$tum})) {
				$sample{$patient}{$tum}{type} = $tum_types{$tum};
			} else {
				die("Error: tumour type not found for $tum\n")
			}
		}
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

if (defined($metadata_file)) {
	open RELATED, ">$prefix-related_tumours_all.tsv" or die "Cannot write to file $prefix-related_tumours_all.tsv\n";
	open RELATED_UNMATCHED, ">$prefix-related_tumours_unmatched.tsv" or die "Cannot write to file $prefix-related_tumours_unmatched.tsv\n";
	open RELATED_MATCHED, ">$prefix-related_tumours_matched.tsv" or die "Cannot write to file $prefix-related_tumours_matched.tsv\n";
}

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
    # Optionally, check tumour types
	# One-per-patient tumour priority: P, R, M, U, "-"
	my $found = undef;
	if (defined($metadata_file)) {
		foreach my $valid_type (keys(%valid_types)) {
			foreach my $i (0..$#tums) {
				my $sample_tum_type = $sample{$patient}{$tums[$i]}{type};
				if ($sample_tum_type eq $valid_type) {
					print $fh1 "$tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n";
					print ONE "$tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n";
					$found = $i;
					last;
				}
			}
			last if defined($found);
		}
	} else {
		print $fh1 "$tums[0]\t$sample{$patient}{$tums[0]}{norm}\n";
		print ONE "$tums[0]\t$sample{$patient}{$tums[0]}{norm}\n";
		$found = 0;
	}

	foreach my $i (0..$#tums) {
		next if $i == $found;
		print STDERR "Excluding lesion from same patient (one-per-patient list): $tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n$sample{$patient}{$tums[$i]}{line}\n";
	}

	# Print all independent and related tumours
	foreach my $i (0..$#tums) {
		my $fh2;
		my $norm = $sample{$patient}{$tums[$i]}{norm};
		if (!defined($metadata_file) || (defined($metadata_file) && $sample{$patient}{$tums[$i]}{type} eq "P")) {
			if ($norm eq 'PDv38is_wes_v2') {
				#print INDEP_UNMATCHED "$tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n";
				$fh2 = *INDEP_UNMATCHED;
			} else {
				#print INDEP_MATCHED "$tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n";
				$fh2 = *INDEP_MATCHED;
			}
			print INDEP "$tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n";
		} else {
			if ($norm eq 'PDv38is_wes_v2') {
				$fh2 = *RELATED_UNMATCHED;
			} else {
				$fh2 = *RELATED_MATCHED;
			}
			print RELATED "$tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n";
		}
		print $fh2 "$tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n";
	}
}

close ONE;
close INDEP;
close ONE_UNMATCHED;
close INDEP_UNMATCHED;
close ONE_MATCHED;
close INDEP_MATCHED;

if (defined($metadata_file)) {
	close RELATED;
	close RELATED_MATCHED;
	close RELATED_UNMATCHED;
}

# Create files with just the tumours

foreach my $file ('one_tumour_per_patient', 'independent_tumours', 'analysed', 'related_tumours') {
	last if $file eq "related_tumours" && !defined($metadata_file);
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

foreach my $file ('one_tumour_per_patient', 'independent_tumours', 'analysed', 'related') {
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

One tumour per patient is selected from the QC-pass/non-rejected
list, based on alphabetical order. If a metadata file with tumour
type was provided, tumour type is further prioritised in order:
P, R, M, U, "-"

independent_tumours
-------------------

All tumours from each patient from the QC-pass/non-rejected list,
annotated as originating different tumours, not samples from the same patient.
If a metadata file was provided, this list will only include primary tumours
and all other tumours will be classified as 'related tumours'.

analysed_tumours
----------------

All tumours from each patient from the QC-pass/non-rejected list,
including samples from the same tumour, that were run through
Caveman/Pindel

related_tumours
---------------
If a metadata file with tumour type was provided for sample list generation,
non-primary tumours will be separated from primary tumours in the independent
tumours list, and placed in the related tumours matched/unmatched lists

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

