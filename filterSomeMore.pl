#!/usr/bin/perl

use warnings;
use strict;

# filter vcf files


### stringency variables, edits as desired
my $maxCoverage =  8433; # maximum depth to avoid repeats, mean + 3sd == gs
#my $maxCoverage =  8142; # maximum depth to avoid repeats, mean + 3sd == gus

my $in = shift(@ARGV);
open (IN, $in) or die "Could not read the infile = $in\n";
$in =~ m/^([a-zA-Z0-9_]+\.vcf)$/ or die "Failed to match the variant file\n";
open (OUT, "> morefilter_$1") or die "Could not write the outfile\n";

my $flag = 0;
my $cnt = 0;

while (<IN>){
	chomp;
	$flag = 1;
	print "\n";
	if (m/^\#/){ ## header row, always write
		$flag = 1;
	}
	elsif (m/^Sc/){ ## this is a sequence line, you migh need to edit this reg. expr.
		$flag = 1;
		m/DP=(\d+)/ or die "Syntax error, DP not found\n";
		if ($1 > $maxCoverage){
			$flag = 0;
			print "fail DP\n";
		}
		if ($flag == 1){
			$cnt++; ## this is a good SNV
		}
	}
	else{
		print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
		$flag = 0;
	}
	if ($flag == 1){
		print OUT "$_\n";
	}
}
close (IN);
close (OUT);

print "Finished filtering $in\nRetained $cnt variable loci\n";
