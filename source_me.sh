#!/bin/bash

# R and R libraries 

_R_BIN_PATH="/software/team113/dermatlas/R/R-4.2.2/bin"
export PATH="${_R_BIN_PATH:?empty-path-variable}:${PATH}"
export R_LIBS="/software/team113/dermatlas/R/R-4.2.2/lib/R/library"

# bcftools (note: loading this bcftools causes 'perl' to point to perl v5.36)

module load bcftools-1.9/python-3.11.6

# Perl and Perl libraries

module load perl-5.38.0
export PERL5LIB=/software/team113/dermatlas/perl5/5.38.0/perl5/lib/perl5

