#!/usr/bin/env perl

# kw10
# This script parses tsv versions of DERMATLAS
# metadata files and list T/N pairs (biosamples),
# whether samples should be analysed, phenotype, 
# diagnosis and comments

use strict;
use warnings;
use List::MoreUtils qw(first_index indexes uniq);
use Data::Dumper;
use Getopt::Std;

# Get options

my %opts = ();
getopts("hi:", \%opts);

if (! keys(%opts) || $opts{h} ) {
	print STDERR ("Usage: $0 -i metadata_files.tsv\n");
	exit;
}

# Columns to parse

my @get_col = (
	"Diagnosis",
	"Phenotype",
	"T/N linked",
	"Sanger DNA ID",
	"OK_to_analyse_DNA?",
	"Sanger RNA ID",
	"OK_to_analyse_RNA?",
	"Clinical notes",
	"Comments",
	"Collaborator"
);

# Outout headers

my @headers = qw(
	Tumour_sample
	Normal_sample
	OK_to_analyse_tum
	OK_to_analyse_norm
	Phenotype_tum
	Phenotype_norm
	Diagnosis_tum
	Diagnosis_norm
	Clinical_notes
	Comments
	Collaborator
);

print join("\t", @headers) . "\n";

my $headerfound = 0;
my %idx;
my %biosample;
my $infile = $opts{i};

# Read and parse the file

open INFILE, "<$infile" or die "Can't open file $infile\n";

while (<INFILE>) {
	my $line = $_;
	chomp $line;
	my @col = split(/\t/, $line);
	if (!$headerfound && $line =~ /^LvdW ID/) {
		foreach my $col (@get_col) {
			$idx{$col} = first_index { $_ eq $col } @col;
		}
		$headerfound = 1;
#		print Dumper %idx;
		next;
	}	
	push @{$biosample{ $col[ $idx{"T/N linked"} ] }{"Diagnosis"} }, $col[ $idx{"Diagnosis"} ];
	push @{$biosample{ $col[ $idx{"T/N linked"} ] }{"Phenotype"} }, $col[ $idx{"Phenotype"} ];
	push @{$biosample{ $col[ $idx{"T/N linked"} ] }{"Sanger DNA ID"} }, $col[ $idx{"Sanger DNA ID"} ];
	push @{$biosample{ $col[ $idx{"T/N linked"} ] }{"OK_to_analyse_DNA"} }, $col[ $idx{"OK_to_analyse_DNA?"} ];
	push @{$biosample{ $col[ $idx{"T/N linked"} ] }{"Sanger RNA ID"} }, $col[ $idx{"Sanger RNA ID"} ];
	push @{$biosample{ $col[ $idx{"T/N linked"} ] }{"OK_to_analyse_RNA"} }, $col[ $idx{"OK_to_analyse_RNA?"} ];
	push @{$biosample{ $col[ $idx{"T/N linked"} ] }{"Clinical notes"} }, $col[ $idx{"Clinical notes"} ];
	push @{$biosample{ $col[ $idx{"T/N linked"} ] }{"Comments"} }, $col[ $idx{"Comments"} ];
	push @{$biosample{ $col[ $idx{"T/N linked"} ] }{"Collaborator"} }, $col[ $idx{"Collaborator"} ];
}

close INFILE;

# Get all T/N pairs

#print Dumper %biosample;

my @printlines;

for my $group (sort keys %biosample) {
	my @tum_idx = indexes { $_ =~ /T/ } @{ $biosample{$group}{"Phenotype"} };
	my @norm_idx = indexes { $_ =~ /N/ } @{ $biosample{$group}{"Phenotype"} };
	foreach my $tum_idx (@tum_idx) {
		my $tum = $biosample{$group}{"Sanger DNA ID"}[$tum_idx];
		next if $tum eq '-';
		my $ok_tum_dna = $biosample{$group}{"OK_to_analyse_DNA"}[$tum_idx];
		my $pheno_tum = $biosample{$group}{"Phenotype"}[$tum_idx];
		my $diagnosis_tum = $biosample{$group}{"Diagnosis"}[$tum_idx];
		$biosample{$group}{"Clinical notes"}[$tum_idx] = "-" if !$biosample{$group}{"Clinical notes"}[$tum_idx];
		$biosample{$group}{"Comments"}[$tum_idx] = "-" if !$biosample{$group}{"Comments"}[$tum_idx];
		if (@norm_idx < 1) {
			my $clinical = $biosample{$group}{"Clinical notes"}[$tum_idx];
			my $comments = $biosample{$group}{"Comments"}[$tum_idx];
			my $collab = $biosample{$group}{"Collaborator"}[$tum_idx];
			push @printlines, join("\t", $tum, "-", $ok_tum_dna, "-", $pheno_tum, "-", $diagnosis_tum, "-", $clinical, $comments, $collab) . "\n";
			next;
		}
		foreach my $norm_idx (@norm_idx) {
			my $norm = $biosample{$group}{"Sanger DNA ID"}[$norm_idx];
			my $ok_norm_dna = $biosample{$group}{"OK_to_analyse_DNA"}[$norm_idx];
			my $pheno_norm = $biosample{$group}{"Phenotype"}[$norm_idx];
			my $diagnosis_norm = $biosample{$group}{"Diagnosis"}[$norm_idx];
			$biosample{$group}{"Clinical notes"}[$norm_idx] = "-" if !$biosample{$group}{"Clinical notes"}[$norm_idx];
			$biosample{$group}{"Comments"}[$norm_idx] = "-" if !$biosample{$group}{"Comments"}[$norm_idx];
			my $clinical = join(";", uniq( $biosample{$group}{"Clinical notes"}[$tum_idx], $biosample{$group}{"Clinical notes"}[$norm_idx] ));
			my $comments = join(";", uniq( $biosample{$group}{"Comments"}[$tum_idx], $biosample{$group}{"Comments"}[$norm_idx] ));
			my $collab = join(";", uniq( $biosample{$group}{"Collaborator"}[$tum_idx], $biosample{$group}{"Collaborator"}[$norm_idx] ));
			push @printlines, join("\t", $tum, $norm, $ok_tum_dna, $ok_norm_dna, $pheno_tum, $pheno_norm, $diagnosis_tum, $diagnosis_norm, $clinical, $comments, $collab) . "\n";
		}
	}
}

my @printlines_rna;

for my $group (sort keys %biosample) {
	foreach my $i (0..scalar( @{ $biosample{$group}{"OK_to_analyse_RNA"} } )  - 1 ) {
		my $rna = $biosample{$group}{"Sanger RNA ID"}[$i];
		next if $rna eq '-';
		my $ok_rna = $biosample{$group}{"OK_to_analyse_RNA"}[$i];
		my $pheno_rna = $biosample{$group}{"Phenotype"}[$i];
		my $diagnosis_rna = $biosample{$group}{"Diagnosis"}[$i];
		$biosample{$group}{"Clinical notes"}[$i] = "-" if !$biosample{$group}{"Clinical notes"}[$i];
		$biosample{$group}{"Comments"}[$i] = "-" if !$biosample{$group}{"Comments"}[$i];
		my $clinical = $biosample{$group}{"Clinical notes"}[$i];
		my $comments = $biosample{$group}{"Comments"}[$i];
		my $collab = $biosample{$group}{"Collaborator"}[$i];
		push @printlines_rna, join("\t", $rna, "-", $ok_rna, "-", $pheno_rna, "-", $diagnosis_rna, "-", $clinical, $comments, $collab) . "\n";
	}
}

print sort(@printlines);
print sort(@printlines_rna);
