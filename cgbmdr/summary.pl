#!/usr/bin/perl
use strict;
use warnings;

print "Hello, World...\n";
open(R,">summary.txt");
for (my $i = 0; $i < 100; $i++) {
	my $file="1loci_result".$i.".txt";
	open(F, "<$file") or die;
	while(<F>) {
		chop();
		if($_=~/Best estimation/) {
			print R $_, "\n";
			last;
		}
	}
}
