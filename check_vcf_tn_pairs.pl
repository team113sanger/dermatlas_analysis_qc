#!/usr/bin/env perl

use strict;
use warnings;

my $manifest = shift @ARGV;
my $projdir = shift @ARGV;

# Check manifest and projdir

if (! -e $manifest) {
	die "File $manifest does not exist\n";
} elsif (! -d $projdir) {
	die "Directory $projdir does not exist\n";
}

my %pairs;

open M, "<$manifest" or die "Can't open $manifest";
while (<M>) {
	my $line = $_;
	my @cols = split(/\t/, $line);
	chomp @cols;
	if ($line =~ /PDID/) {
		if ($cols[3] !~ /Tumour Sample/ || $cols[4]  !~ /Normal Sample/) {
			die "Manifest file is not the correct format\n";
		}
		next;
	}
	if ($cols[7] =~ /y|yes/i) {
		print STDERR "Found $cols[7] $cols[3] $cols[4]\n";
		$pairs{$cols[3]} = $cols[4];
		$pairs{$cols[3]} = $cols[4];
	}
}
close M;

my @cavemanlist = `dir -1 $projdir/analysis/caveman_files/*/*smartphase.flag.vcf.gz`;
my @pindellist = `dir -1 $projdir/analysis/pindel_files/*/*.flagged.vcf.gz`;
chomp @cavemanlist;
chomp @pindellist;

my %vcfpairs;

foreach my $file (@cavemanlist, @pindellist) {
	open F, "zcat $file | grep \"#\" |" or die "Can't open file $file\n";
	my $tum;
	my $norm;
	my $type;
	while (<F>) {
		if (/^##SAMPLE/ && /NORMAL/) {
			if (/SampleName=(\w+)/) {
				$norm = $1;
			} else {
				$norm = "NOT_FOUND";
			}
		} elsif (/^##SAMPLE/ && /TUMOUR/) {
			if (/SampleName=(\w+)/) {
				$tum = $1;
			} else {
				$tum = "NOT_FOUND";
			}
		} elsif (/vcfProcessLog.+pindel/i) {
			$type = "pindel";
		} elsif (/vcfProcessLog.+caveman/i) {
			$type = "caveman";
		}
	}
	close F;
	print STDERR "Found $type $tum $norm\n";
	$vcfpairs{$type}{$tum} = $norm;
}

# Check all tn pairs in the manifest 

my @tumchecked;

foreach my $type (sort keys %vcfpairs) {
	foreach my $tum (sort keys %pairs) {
		push @tumchecked, $tum;
		if (exists($vcfpairs{$type}{$tum}) && $vcfpairs{$type}{$tum} eq $pairs{$tum}) {
			print "$type\t$tum\t$pairs{$tum}\t$tum\t$vcfpairs{$type}{$tum}\tFOUND\n";
		} elsif (exists($vcfpairs{$type}{$tum}) && $vcfpairs{$type}{$tum} ne $pairs{$tum}) {
			print "$type\t$tum\t$pairs{$tum}\t$$tum\t$vcfpairs{$type}{$tum}\tDIFFERENT_NORMAL\n";
		} else {
			print "$type\t$tum\t$pairs{$tum}\tNA\tNA\tDIFFERENT_NORMAL\n";
		}
	}
}

# Check if there any extra tn pairs in the vcfs

foreach my $vcftype (sort keys %vcfpairs) {
	foreach my $tum (sort keys %{$vcfpairs{$vcftype}}) {
		next if grep( /^$tum$/, @tumchecked);
		print "$vcftype\t$tum\t$vcfpairs{$vcftype}{$tum}\tNOT_IN_MANIFEST\n";
	}
}


