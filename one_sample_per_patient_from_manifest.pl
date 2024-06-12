#!/usr/bin/env perl

use strict;
use warnings;

# Read in a completed biosample manifest and print out t/n pairs,
# one per patient

my $infile = shift @ARGV or die "Input a completed biosample manifest and optional file with a list of samples to exclude\n";
my $rejected = shift @ARGV;

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

open F, "<$infile" or die "Can't open $infile\n";
while (<F>) {
	if (/Biosample/) {
		next;
	}
	chomp;
	my $line = $_;
	my @col = split(/\t/, $line);
	if ($col[8] =~ /Y|y/) {
		my $norm = $col[4];
		my $tum = $col[3];
		if ($reject{$tum} && $reject{$tum} == 1) {
			print STDERR "Skipping $tum\n";
			next;
		} elsif ($reject{$norm} && $reject{$norm} == 1 && !$reject{$tum}) {
			print STDERR "Warning: $norm is on the reject list and $tum is not. Not excluding\n";
		}
		my $patient = $1 if $tum =~ /(\S+)\S$/;
		$sample{$patient}{$tum}{line} = $line;
		$sample{$patient}{$tum}{norm} = $norm;
	}
}
close F;

foreach my $patient (sort keys %sample) {
	my @tums = (sort keys %{$sample{$patient}});
	print "$tums[0]\t$sample{$patient}{$tums[0]}{norm}\n";
	##print $sample{$patient}{$tums[0]}[0] . "\n";
	foreach my $i (1..$#tums) {
		print STDERR "Excluding lesion from same patient: $tums[$i]\t$sample{$patient}{$tums[$i]}{norm}\n$sample{$patient}{$tums[$i]}{line}\n";
	}
}



