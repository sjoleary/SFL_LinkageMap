#/usr/bin/perl -w

use strict;

my $map = shift;
my $fasta = shift;

open(MAP, "<", $map) or die $!;
open(FAS, "<", $fasta) or die $!;
open(OUT, ">", 'filter.fasta') or die $!;

my %mapped;
<MAP>;
while(<MAP>) {
	chomp;
	my ($locus, $lg, $pos) = split;
	$mapped{$locus} = 1;
}


my $keep = 0;
while(<FAS>) {
	if ($_ =~ /^>([\w\d]+)/) { # a header line
		my $locus = $1;
		if ($mapped{$locus}) {
			print OUT $_;
			$keep = 1;
			next;
		} else {
			next;
		}
	} else { # a sequence line
		if ($keep == 1) {
			print OUT $_;
			$keep = 0;
			next;
		} else {
			next;
		}
	}
}

# Chris Hollenbeck