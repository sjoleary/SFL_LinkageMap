#!/usr/bin/perl -w

# Chris Hollenbeck

use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use List::MoreUtils qw/uniq/;
use Bio::PopGen::Population;
use Bio::PopGen::Individual;

open(DUMP, ">", 'dumper.out');

pod2usage(-verbose => 1) if @ARGV == 0;

my $infile = '';
my $informat = '';
my $outfile = '';
my $outformat = '';
my $sex_specific = '';
my $female_id = '';
my $male_id = '';
my $loc_cutoff = 0;
my $ind_cutoff = 0;
my $mapfile = '';

my @ind_blacklist;



GetOptions(		'infile|i=s' => \$infile,
			'informat|n=s' => \$informat,
			'outfile|o=s' => \$outfile,
			'outformat|u=s' => \$outformat,
			'sex_specific|s=s' => \$sex_specific,
			'female_id|F=s' => \$female_id,
			'male_id|M=s' => \$male_id,
			'loc_cutoff|l=s' => \$loc_cutoff,
			'ind_cutoff|c=s' => \$ind_cutoff,
			'map|a=s' => \$mapfile,

	);

if (! $female_id || ! $male_id) {
	die "Female and male parent IDs are required (-F and -M flags)\n";
}

my $pop; # The Bio::PopGen::Population object for the cross
my @progeny; # An array of Bio::PopGen::Individual objects for progeny
my @prog_names; # An array of progeny names

if ($informat eq 'csv') {
	my $ind_names;
	($pop, $ind_names) = read_csv($infile);
	foreach my $ind ($pop->get_Individuals) {
		if ($ind->unique_id() eq $female_id || $ind->unique_id() eq $male_id) {
			next;
		}
		push @progeny, $ind;
	}
	# Loop through one more time to preserve the order of the individuals in the original file
	foreach my $ind (@$ind_names) {
		if ($ind eq $female_id || $ind eq $male_id) {
			next;
		}
		push @prog_names, $ind;
	}
}

if ($informat eq 'tsv') {
	my $ind_names;
	($pop, $ind_names) = read_tsv($infile, $female_id, $male_id);
	foreach my $ind ($pop->get_Individuals) {
		if ($ind->unique_id() eq $female_id || $ind->unique_id() eq $male_id) {
			next;
		}
		push @progeny, $ind;
	}
	# Loop through one more time to preserve the order of the individuals in the original file
	foreach my $ind (@{$ind_names}) {
		if ($ind eq $female_id || $ind eq $male_id) {
			next;
		}
		push @prog_names, $ind;
	}
}

if ($informat eq 'joinmap') {
	my $ind_names;
	($pop, $ind_names) = read_joinmap($infile, $female_id, $male_id);
	foreach my $ind ($pop->get_Individuals) {
		if ($ind->unique_id() eq $female_id || $ind->unique_id() eq $male_id) {
			next;
		}
		push @progeny, $ind;
	}
	# Loop through one more time to preserve the order of the individuals in the original file
	foreach my $ind (@$ind_names) {
		if ($ind eq $female_id || $ind eq $male_id) {
			next;
		}
		push @prog_names, $ind;
	}
}


if ($outformat eq 'mm') {
	if (! $sex_specific) {
		die "Need to specify parent with mm output (use -s flag)\n";
	}
}

# Sort the progeny based on the order in the original file

my @sorted_progeny;
foreach my $ind (@prog_names) {
	foreach my $ind_obj (@progeny) {
		if ($ind_obj->unique_id() eq $ind) {
			push @sorted_progeny, $ind_obj;
			last;
		}
	}
}

@progeny = @sorted_progeny;


#print DUMP Dumper(\$pop);

#Initialize a hash to keep track of various stats

my %stats = (
	'both_hom' => [0, [ ]],
	'id_hets' => [0, [ ]],
	'miss_f' => [0, [ ]],
	'miss_m' => [0, [ ]],
	'strange_seg' => [0, [ ]],
	'total_loci' => [0, [ ]],
	'added_loci' => [0, [ ]],
	'excluded_loci' => [0, [ ]],
	'total_inds' => [0, [ ]],
	'added_inds' => [0, [ ]],
	'excluded_inds' => [0, [ ]],
);

my %ind_missing;
my %loc_missing;

my @whitelist_progeny;
foreach my $ind (@progeny) {
	$stats{'total_inds'}[0]++;
	push @{$stats{'total_inds'}[1]}, $ind->unique_id();
	my $ind_id = $ind->unique_id();
	if (grep(/$ind_id/, @ind_blacklist)) {
		$stats{'excluded_inds'}[0]++;
		push @{$stats{'excluded_inds'}[1]}, $ind->unique_id();
	}
	push @whitelist_progeny, $ind;
	$stats{'added_inds'}[0]++;
	push @{$stats{'added_inds'}[1]}, $ind->unique_id();
}

@progeny = @whitelist_progeny;

my %progeny_codes;
my %progeny_codes_f;
my %progeny_codes_m;
my @lines;
my @lines1;
my @lines2;
open(ERR, ">", 'errors.txt') or die $!;
#open(TEMP, ">", 'temp.out');

# Cycle through each locus
foreach my $locus ($pop->get_marker_names) {
	my @male_ind = $pop->get_Individuals(-unique_id => $male_id);
	my @female_ind = $pop->get_Individuals(-unique_id => $female_id);
	my @m_geno = $male_ind[0]->get_Genotypes(-marker => $locus);
	my @f_geno = $female_ind[0]->get_Genotypes(-marker => $locus);
	my ($m0, $m1) = $m_geno[0]->get_Alleles;
	my ($f0, $f1) = $f_geno[0]->get_Alleles;
	#print join("\t", $locus, $m0, $m1, $f0, $f1), "\n";

	$stats{'total_loci'}[0]++;
	push @{$stats{'total_loci'}[1]}, $locus;

	if (! $f0 || ! $f1) {	# Missing maternal genotype
		$stats{'miss_f'}[0]++;
		push @{$stats{'miss_f'}[1]}, $locus;
		$stats{'excluded_loci'}[0]++;
		push @{$stats{'excluded_loci'}[1]}, $locus;
		#print join("-", $locus, $f0, $f1), "\n";
		next;
	}

	if (! $m0 || ! $m1) {	# Missing paternal genotype
		$stats{'miss_m'}[0]++;
		push @{$stats{'miss_m'}[1]}, $locus;
		$stats{'excluded_loci'}[0]++;
		push @{$stats{'excluded_loci'}[1]}, $locus;
		next;
	}

	if ($f0 eq $f1 && $m0 eq $m1) {  # Both homozygotes
		$stats{'both_hom'}[0]++;
		push @{$stats{'both_hom'}[1]}, $locus;
		$stats{'excluded_loci'}[0]++;
		push @{$stats{'excluded_loci'}[1]}, $locus;
		next;
	}
	if ($f0 eq $m0 && $f1 eq $m1) {	# Identical heterozygotes
		$stats{'id_hets'}[0]++;
		push @{$stats{'id_hets'}[1]}, $locus;
		$stats{'excluded_loci'}[0]++;
		push @{$stats{'excluded_loci'}[1]}, $locus;
		next;
	}


	my @alleles = ($f0, $f1, $m0, $m1);
	my $unique = scalar(uniq(@alleles));

	# Determine the segregation type of the locus, based on Table 4 in the JoinMap4 documentation (with slight modification). In brief:
	#
	# Type 'A': <abxcd> - both parents heterozygotes, four alleles
	# Type 'B': <efxeg> - both parents heterozygotes, three alleles
	# Type 'C': <hkxhk> - both parents heterozygotes, two alleles
	# Type 'D': <lmxll> - first parent heterozygous
	# Type 'E': <nnxnp> - second parent heterozygous
	#
	# The modification of this in this program is that type 'C' loci are excluded (because the alleles can't be tracked) and a different case
	# of segregation <abxcc> (or alternatively <aaxbc>) is called type 'C' segregation initially. Because this is effectively the same as Type 'D' or 'E',
	# depending on the heterozygous parent, type 'C' is later converted to either 'D' or 'E'.
	#

	my $seg_type;
	if ($unique == 2) {
		if ($f0 eq $f1) {
			$seg_type = 'E';
		} else {
			$seg_type = 'D';
		}
	} elsif ($unique == 3) {
		if ($f0 ne $f1 && $m0 ne $m1) {
			$seg_type = 'B';
		} else {
			$seg_type = 'C';
		}
	} elsif ($unique == 4) {
		$seg_type = 'A';
	}

	if (! $seg_type) {
		$stats{'strange_seg'}[0]++;
		push @{$stats{'strange_seg'}[1]}, $locus;
		$stats{'excluded_loci'}[0]++;
		push @{$stats{'excluded_loci'}[1]}, $locus;
		next;
	}


	# Determine the proper coding for each individual given the genotype and segregation type
	$progeny_codes{$locus} = [];

	if ($outformat eq 'joinmap') {
		my $type;
		my @progeny_codes;
		if ($seg_type eq 'A') {
			if ($sex_specific eq 'F') {
				$type = '<lmxll>';
			} elsif ($sex_specific eq 'M') {
				$type = '<nnxnp>';
			} else {
				$type = '<abxcd>';
			}
			foreach my $ind (@progeny) {

				my @ind_geno = $ind->get_Genotypes(-marker => $locus);
				my @ind_alleles = $ind_geno[0]->get_Alleles();
				my $geno_string;
				if ($ind_alleles[0] eq '' || $ind_alleles[1] eq '') {
					$geno_string = '';
				} else {
					$geno_string = $ind_alleles[0] . $ind_alleles[1];
				}

				my $a = $f0;
				my $b = $f1;
				my $c = $m0;
				my $d = $m1;
				if ($geno_string =~ /$a$c/i || $geno_string =~ /$c$a/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'lm';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'np';
					} else {
						push @progeny_codes, 'ac';
					}
				} elsif ($geno_string =~ /$a$d/i || $geno_string =~ /$d$a/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'lm';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'nn';
					} else {
						push @progeny_codes, 'ad';
					}
				} elsif ($geno_string =~ /$b$c/i || $geno_string =~ /$c$b/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'll';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'np';
					} else {
						push @progeny_codes, 'bc';
					}
				} elsif ($geno_string =~ /$b$d/i || $geno_string =~ /$d$b/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'll';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'nn';
					} else {
						push @progeny_codes, 'bd';
					}
				} elsif ($geno_string eq '') {
					push @progeny_codes, '--';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
				} else {
					push @progeny_codes, '--';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
					print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
				}
			}
		} elsif ($seg_type eq 'B') {
			if ($sex_specific eq 'F') {
				$type = '<lmxll>';
			} elsif ($sex_specific eq 'M') {
				$type = '<nnxnp>';
			} else {
				$type = '<efxeg>';
			}
			foreach my $ind (@progeny) {

				my @ind_geno = $ind->get_Genotypes(-marker => $locus);
				my @ind_alleles = $ind_geno[0]->get_Alleles();
				my $geno_string;
				if ($ind_alleles[0] eq '' || $ind_alleles[1] eq '') {
					$geno_string = '';
				} else {
					$geno_string = $ind_alleles[0] . $ind_alleles[1];
				}

				my $e = $f0;
				my $f = $f1;
				my $g;
				if ($m0 eq $e) {
					 $g = $m1;
				} elsif ($m1 eq $e) {
					$g = $m0;
				} else {
					$e = $f1;
					$f = $f0;
					if ($m0 eq $e) {
						$g = $m1;
					} else {
						$g = $m0;
					}
				}
				if ($geno_string =~ /$e$e/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'lm';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'np';
					} else {
						push @progeny_codes, 'ee';
					}
				} elsif ($geno_string =~ /$e$f/i || $geno_string =~ /$f$e/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'll';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'np';
					} else {
						push @progeny_codes, 'ef';
					}
				} elsif ($geno_string =~ /$e$g/i || $geno_string =~ /$g$e/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'lm';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'nn';
					} else {
						push @progeny_codes, 'eg';
					}
				} elsif ($geno_string =~ /$f$g/i || $geno_string =~ /$g$f/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'll';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'nn';
					} else {
						push @progeny_codes, 'fg';
					}
				} elsif ($geno_string eq '') {
					push @progeny_codes, '--';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
				} else {
					push @progeny_codes, '--';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
					print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
				}
			}
		} elsif ($seg_type eq 'C') {

			my $x;
			my $n;
			my $p;
			my $l;
			my $m;

			if ($f0 eq $f1) {
				$type = '<nnxnp>';
				if ($sex_specific eq 'F') {
					$stats{'excluded_loci'}[0]++;
					push @{$stats{'excluded_loci'}[1]}, $locus;
					next;
				}
				$x = $f0;
				$n = $m0;
				$p = $m1;

			} else {
				$type = '<lmxll>';
				if ($sex_specific eq 'M') {
					$stats{'excluded_loci'}[0]++;
					push @{$stats{'excluded_loci'}[1]}, $locus;
					next;
				}
				$l = $f0;
				$m = $f1;
				$x = $m0;
			}


			foreach my $ind (@progeny) {

				my @ind_geno = $ind->get_Genotypes(-marker => $locus);
				my @ind_alleles = $ind_geno[0]->get_Alleles();
				my $geno_string;
				if ($ind_alleles[0] eq '' || $ind_alleles[1] eq '') {
					$geno_string = '';
				} else {
					$geno_string = $ind_alleles[0] . $ind_alleles[1];
				}

				if ($type eq '<nnxnp>') {
					if ($geno_string =~ /$n$x/i || $geno_string =~ /$x$n/i) {
						push @progeny_codes, 'nn';
					} elsif ($geno_string =~ /$x$p/i || $geno_string =~ /$p$x/i) {
						push @progeny_codes, 'np';
					} elsif ($geno_string eq '') {
						push @progeny_codes, '--';
						$ind_missing{$ind->unique_id()}++;
						$loc_missing{$locus}++;
					} else {
						push @progeny_codes, '--';
						$ind_missing{$ind->unique_id()}++;
						$loc_missing{$locus}++;
						print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
					}
				} else {
					if ($geno_string =~ /$l$x/i || $geno_string =~ /$x$l/i) {
						push @progeny_codes, 'll';
					} elsif ($geno_string =~ /$m$x/i || $geno_string =~ /$x$m/i) {
						push @progeny_codes, 'lm';
					} elsif ($geno_string eq '') {
						push @progeny_codes, '--';
						$ind_missing{$ind->unique_id()}++;
						$loc_missing{$locus}++;
					} else {
						push @progeny_codes, '--';
						$ind_missing{$ind->unique_id()}++;
						$loc_missing{$locus}++;
						print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
					}
				}

			}
		} elsif ($seg_type eq 'D') {
			$type = '<lmxll>';
			if ($sex_specific eq 'M') {
				$stats{'excluded_loci'}[0]++;
				push @{$stats{'excluded_loci'}[1]}, $locus;
				next;
			}
			foreach my $ind (@progeny) {

				my @ind_geno = $ind->get_Genotypes(-marker => $locus);
				my @ind_alleles = $ind_geno[0]->get_Alleles();
				my $geno_string;
				if ($ind_alleles[0] eq '' || $ind_alleles[1] eq '') {
					$geno_string = '';
				} else {
					$geno_string = $ind_alleles[0] . $ind_alleles[1];
				}

				my $l;
				my $m;
				if ($f0 eq $m0) {
					$l = $f0;
					$m = $f1;
				} else {
					$l = $f1;
					$m = $f0;
				}

				if ($geno_string =~ /$l$l/i || $geno_string =~ /$m$m/i) {
					push @progeny_codes, 'll';
				} elsif ($geno_string =~ /$l$m/i || $geno_string =~ /$m$l/i) {
					push @progeny_codes, 'lm';
				} elsif ($geno_string eq '') {
					push @progeny_codes, '--';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
				} else {
					push @progeny_codes, '--';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
					print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
				}
			}
		} elsif ($seg_type eq 'E') {
			$type = '<nnxnp>';
			if ($sex_specific eq 'F') {
				$stats{'excluded_loci'}[0]++;
				push @{$stats{'excluded_loci'}[1]}, $locus;
				next;
			}
			foreach my $ind (@progeny) {

				my @ind_geno = $ind->get_Genotypes(-marker => $locus);
				my @ind_alleles = $ind_geno[0]->get_Alleles();
				my $geno_string;
				if ($ind_alleles[0] eq '' || $ind_alleles[1] eq '') {
					$geno_string = '';
				} else {
					$geno_string = $ind_alleles[0] . $ind_alleles[1];
				}

				my $n = $f0;
				my $p;
				if ($m0 eq $n) {
					$p = $m1;
				} else {
					$p = $m0;
				}
				if ($geno_string =~ /$n$n/i) {
					push @progeny_codes, 'nn';
				} elsif ($geno_string =~ /$n$p/i || $geno_string =~ /$p$n/i) {
					push @progeny_codes, 'np';
				} elsif ($geno_string eq '') {
					push @progeny_codes, '--';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
				} else {
					push @progeny_codes, '--';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
					print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
				}
			}
		} else {
			$stats{'strange_seg'}[0]++;
			push @{$stats{'strange_seg'}[1]}, $locus;
			$stats{'excluded_loci'}[0]++;
			push @{$stats{'excluded_loci'}[1]}, $locus;
			next;
		}

		$progeny_codes{$locus} = [$type, \@progeny_codes];

		#my $line = join("\t", $locus, $type, @progeny_codes);
		#push @lines1, $line;
		#print join("\t", $locus, $type, @progeny_codes), "\n";

		$stats{'added_loci'}[0]++;
		push @{$stats{'added_loci'}[1]}, $locus;


	}

	if ($outformat eq 'mm') {
		my $missing_data;
		my @progeny_codes_f;
		my @progeny_codes_m;
		if ($seg_type eq 'A') {
			foreach my $ind (@progeny) {
				my @genotypes = $ind->get_Genotypes(-marker => $locus);
				my @alleles = $genotypes[0]->get_Alleles;
				my $geno;
				if ($alleles[0] eq '' || $alleles[1] eq '') {
					$geno = '';
				} else {
					$geno = $alleles[0] . $alleles[1];
				}
				#next if grep(/$ind_names[$ind_counter]/, @no_include);
				if ($geno =~ /$f0$m0/i || $geno =~ /$m0$f0/i) {
					push @progeny_codes_f, 'A';
					push @progeny_codes_m, 'A';
				} elsif ($geno =~ /$f0$m1/i || $geno =~ /$m1$f0/i) {
					push @progeny_codes_f, 'A';
					push @progeny_codes_m, 'H';
				} elsif ($geno =~ /$f1$m0/i || $geno =~ /$m0$f1/i) {
					push @progeny_codes_f, 'H';
					push @progeny_codes_m, 'A';
				} elsif ($geno =~ /$f1$m1/i || $geno =~ /$m1$f1/i) {
					push @progeny_codes_f, 'H';
					push @progeny_codes_m, 'H';
				} elsif ($geno eq '') {
					push @progeny_codes_f, '-';
					push @progeny_codes_m, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
				} else {
					push @progeny_codes_f, '-';
					push @progeny_codes_m, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
					print ERR "Error: Check locus $locus, genotype ", $geno, "\n";
				}
			}
		} elsif ($seg_type eq 'B') {
			foreach my $ind (@progeny) {
				my @genotypes = $ind->get_Genotypes(-marker => $locus);
				my @alleles = $genotypes[0]->get_Alleles;
				my $geno;
				if ($alleles[0] eq '' || $alleles[1] eq '') {
					$geno = '';
				} else {
					$geno = $alleles[0] . $alleles[1];
				}
				#next if grep(/$ind_names[$ind_counter]/, @no_include);
				if ($geno =~ /$f0$m0/i || $geno =~ /$m0$f0/i) {
					push @progeny_codes_f, 'A';
					push @progeny_codes_m, 'A';
				} elsif ($geno =~ /$f0$m1/i || $geno =~ /$m1$f0/i) {
					push @progeny_codes_f, 'A';
					push @progeny_codes_m, 'H';
				} elsif ($geno =~ /$m0$f1/i || $geno =~ /$f1$m0/i) {
					push @progeny_codes_f, 'H';
					push @progeny_codes_m, 'A';
				} elsif ($geno =~ /$f1$m1/i || $geno =~ /$m1$f1/i) {
					push @progeny_codes_f, 'H';
					push @progeny_codes_m, 'H';
				} elsif ($geno eq '') {
					push @progeny_codes_f, '-';
					push @progeny_codes_m, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
				} else {
					push @progeny_codes_f, '-';
					push @progeny_codes_m, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
					print ERR "Error: Check locus $locus, genotype ", $geno, "\n";
				}
			}

		} elsif ($seg_type eq 'C') {
			my $type;
			my $x;
			my $n;
			my $p;
			my $l;
			my $m;

			if ($f0 eq $f1) {
				$type = '<nnxnp>';
				if ($sex_specific eq 'F') {
					$stats{'excluded_loci'}[0]++;
					push @{$stats{'excluded_loci'}[1]}, $locus;
					next;
				}
				$x = $f0;
				$n = $m0;
				$p = $m1;

			} else {
				$type = '<lmxll>';
				if ($sex_specific eq 'M') {
					$stats{'excluded_loci'}[0]++;
					push @{$stats{'excluded_loci'}[1]}, $locus;
					next;
				}
				$l = $f0;
				$m = $f1;
				$x = $m0;
			}

			foreach my $ind (@progeny) {

				my @ind_geno = $ind->get_Genotypes(-marker => $locus);
				my @ind_alleles = $ind_geno[0]->get_Alleles();
				my $geno_string;
				if ($ind_alleles[0] eq '' || $ind_alleles[1] eq '') {
					$geno_string = '';
				} else {
					$geno_string = $ind_alleles[0] . $ind_alleles[1];
				}

				if ($type eq '<nnxnp>') {
					if ($geno_string =~ /$n$x/i || $geno_string =~ /$x$n/i) {
						push @progeny_codes_m, 'A';
					} elsif ($geno_string =~ /$x$p/i || $geno_string =~ /$p$x/i) {
						push @progeny_codes_m, 'H';
					} elsif ($geno_string eq '') {
						push @progeny_codes_m, '-';
						$ind_missing{$ind->unique_id()}++;
						$loc_missing{$locus}++;
					} else {
						push @progeny_codes_m, '-';
						$ind_missing{$ind->unique_id()}++;
						$loc_missing{$locus}++;
						print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
					}
				} else {
					if ($geno_string =~ /$l$x/i || $geno_string =~ /$x$l/i) {
						push @progeny_codes_f, 'A';
					} elsif ($geno_string =~ /$m$x/i || $geno_string =~ /$x$m/i) {
						push @progeny_codes_f, 'H';
					} elsif ($geno_string eq '') {
						push @progeny_codes_f, '-';
						$ind_missing{$ind->unique_id()}++;
						$loc_missing{$locus}++;
					} else {
						push @progeny_codes_f, '-';
						$ind_missing{$ind->unique_id()}++;
						$loc_missing{$locus}++;
						print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
					}
				}

			}


		} elsif ($seg_type eq 'D') {

			if ($sex_specific eq 'M') {
				$stats{'excluded_loci'}[0]++;
				push @{$stats{'excluded_loci'}[1]}, $locus;
				next;
			}

			foreach my $ind (@progeny) {
				my @genotypes = $ind->get_Genotypes(-marker => $locus);
				my @alleles = $genotypes[0]->get_Alleles;
				my $geno;
				if ($alleles[0] eq '' || $alleles[1] eq '') {
					$geno = '';
				} else {
					$geno = $alleles[0] . $alleles[1];
				}
				#next if grep(/$ind_names[$ind_counter]/, @no_include);
				if ($geno =~ /$f0$m0/i || $geno =~ /$m0$f0/i) {
					push @progeny_codes_f, 'A';
				} elsif ($geno =~ /$m0$f1/i || $geno =~ /$f1$m0/i) {
					push @progeny_codes_f, 'H';
				} elsif ($geno eq '') {
					push @progeny_codes_f, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
				} else {
					push @progeny_codes_f, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
					print ERR "Error: Check locus $locus, genotype ", $geno, "\n";
				}
			}
		} elsif ($seg_type eq 'E') {

			if ($sex_specific eq 'F') {
				$stats{'excluded_loci'}[0]++;
				push @{$stats{'excluded_loci'}[1]}, $locus;
				next;
			}

			foreach my $ind (@progeny) {
				my @genotypes = $ind->get_Genotypes(-marker => $locus);
				my @alleles = $genotypes[0]->get_Alleles;
				my $geno;
				if ($alleles[0] eq '' || $alleles[1] eq '') {
					$geno = '';
				} else {
					$geno = $alleles[0] . $alleles[1];
				}
				#next if grep(/$ind_names[$ind_counter]/, @no_include);
				if ($geno =~ /$f0$m0/i || $geno =~ /$m0$f0/i) {
					push @progeny_codes_m, 'A';
				} elsif ($geno =~ /$f0$m1/i || $geno =~ /$m1$f0/i) {
					push @progeny_codes_m, 'H';
				} elsif ($geno eq '') {
					push @progeny_codes_m, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
				} else {
					push @progeny_codes_m, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
					print ERR "Error: Check locus $locus, genotype ", $geno, "\n";
				}
			}
		} else {
			$stats{'strange_seg'}[0]++;
			push @{$stats{'strange_seg'}[1]}, $locus;
			$stats{'excluded_loci'}[0]++;
			push @{$stats{'excluded_loci'}[1]}, $locus;
			next;
		}

		$progeny_codes_f{$locus} = [$seg_type, \@progeny_codes_f];
		$progeny_codes_m{$locus} = [$seg_type, \@progeny_codes_m];


		$stats{'added_loci'}[0]++;
		push @{$stats{'added_loci'}[1]}, $locus;


	}

	if ($outformat eq 'onemap') {
		my $type;
		my @progeny_codes;
		if ($seg_type eq 'A') {
			if ($sex_specific eq 'F') {
				$type = 'D1.10';
			} elsif ($sex_specific eq 'M') {
				$type = 'D2.15';
			} else {
				$type = 'A.1';
			}
			foreach my $ind (@progeny) {

				my @ind_geno = $ind->get_Genotypes(-marker => $locus);
				my @ind_alleles = $ind_geno[0]->get_Alleles();
				my $geno_string;
				if ($ind_alleles[0] eq '' || $ind_alleles[1] eq '') {
					$geno_string = '';
				} else {
					$geno_string = $ind_alleles[0] . $ind_alleles[1];
				}

				my $a = $f0;
				my $b = $f1;
				my $c = $m0;
				my $d = $m1;
				if ($geno_string =~ /$a$c/i || $geno_string =~ /$c$a/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'ab';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'ab';
					} else {
						push @progeny_codes, 'ac';
					}
				} elsif ($geno_string =~ /$a$d/i || $geno_string =~ /$d$a/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'ab';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'a';
					} else {
						push @progeny_codes, 'ad';
					}
				} elsif ($geno_string =~ /$b$c/i || $geno_string =~ /$c$b/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'a';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'ab';
					} else {
						push @progeny_codes, 'bc';
					}
				} elsif ($geno_string =~ /$b$d/i || $geno_string =~ /$d$b/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'a';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'a';
					} else {
						push @progeny_codes, 'bd';
					}
				} elsif ($geno_string eq '') {
					push @progeny_codes, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
				} else {
					push @progeny_codes, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
					print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
				}
			}
		} elsif ($seg_type eq 'B') {
			if ($sex_specific eq 'F') {
				$type = 'D1.10';
			} elsif ($sex_specific eq 'M') {
				$type = 'D2.15';
			} else {
				$type = 'A.2';
			}
			foreach my $ind (@progeny) {

				my @ind_geno = $ind->get_Genotypes(-marker => $locus);
				my @ind_alleles = $ind_geno[0]->get_Alleles();
				my $geno_string;
				if ($ind_alleles[0] eq '' || $ind_alleles[1] eq '') {
					$geno_string = '';
				} else {
					$geno_string = $ind_alleles[0] . $ind_alleles[1];
				}

				my $e = $f0;
				my $f = $f1;
				my $g;
				if ($m0 eq $e) {
					 $g = $m1;
				} elsif ($m1 eq $e) {
					$g = $m0;
				} else {
					$e = $f1;
					$f = $f0;
					if ($m0 eq $e) {
						$g = $m1;
					} else {
						$g = $m0;
					}
				}
				if ($geno_string =~ /$e$e/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'ab';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'ab';
					} else {
						push @progeny_codes, 'a';
					}
				} elsif ($geno_string =~ /$e$f/i || $geno_string =~ /$f$e/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'a';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'ab';
					} else {
						push @progeny_codes, 'ba';
					}
				} elsif ($geno_string =~ /$e$g/i || $geno_string =~ /$g$e/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'ab';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'a';
					} else {
						push @progeny_codes, 'ac';
					}
				} elsif ($geno_string =~ /$f$g/i || $geno_string =~ /$g$f/i) {
					if ($sex_specific eq 'F') {
						push @progeny_codes, 'a';
					} elsif ($sex_specific eq 'M') {
						push @progeny_codes, 'a';
					} else {
						push @progeny_codes, 'bc';
					}
				} elsif ($geno_string eq '') {
					push @progeny_codes, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
				} else {
					push @progeny_codes, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
					print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
				}
			}
		} elsif ($seg_type eq 'C') {

			my $x;
			my $n;
			my $p;
			my $l;
			my $m;

			if ($f0 eq $f1) {
				$type = 'D2.15';
				if ($sex_specific eq 'F') {
					$stats{'excluded_loci'}[0]++;
					push @{$stats{'excluded_loci'}[1]}, $locus;
					next;
				}
				$x = $f0;
				$n = $m0;
				$p = $m1;

			} else {
				$type = 'D1.10';
				if ($sex_specific eq 'M') {
					$stats{'excluded_loci'}[0]++;
					push @{$stats{'excluded_loci'}[1]}, $locus;
					next;
				}
				$l = $f0;
				$m = $f1;
				$x = $m0;
			}


			foreach my $ind (@progeny) {

				my @ind_geno = $ind->get_Genotypes(-marker => $locus);
				my @ind_alleles = $ind_geno[0]->get_Alleles();
				my $geno_string;
				if ($ind_alleles[0] eq '' || $ind_alleles[1] eq '') {
					$geno_string = '';
				} else {
					$geno_string = $ind_alleles[0] . $ind_alleles[1];
				}

				if ($type eq 'D2.15') {
					if ($geno_string =~ /$n$x/i || $geno_string =~ /$x$n/i) {
						push @progeny_codes, 'a';
					} elsif ($geno_string =~ /$x$p/i || $geno_string =~ /$p$x/i) {
						push @progeny_codes, 'ab';
					} elsif ($geno_string eq '') {
						push @progeny_codes, '-';
						$ind_missing{$ind->unique_id()}++;
						$loc_missing{$locus}++;
					} else {
						push @progeny_codes, '-';
						$ind_missing{$ind->unique_id()}++;
						$loc_missing{$locus}++;
						print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
					}
				} else {
					if ($geno_string =~ /$l$x/i || $geno_string =~ /$x$l/i) {
						push @progeny_codes, 'a';
					} elsif ($geno_string =~ /$m$x/i || $geno_string =~ /$x$m/i) {
						push @progeny_codes, 'ab';
					} elsif ($geno_string eq '') {
						push @progeny_codes, '-';
						$ind_missing{$ind->unique_id()}++;
						$loc_missing{$locus}++;
					} else {
						push @progeny_codes, '-';
						$ind_missing{$ind->unique_id()}++;
						$loc_missing{$locus}++;
						print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
					}
				}

			}
		} elsif ($seg_type eq 'D') {
			$type = 'D1.10';
			if ($sex_specific eq 'M') {
				$stats{'excluded_loci'}[0]++;
				push @{$stats{'excluded_loci'}[1]}, $locus;
				next;
			}
			foreach my $ind (@progeny) {

				my @ind_geno = $ind->get_Genotypes(-marker => $locus);
				my @ind_alleles = $ind_geno[0]->get_Alleles();
				my $geno_string;
				if ($ind_alleles[0] eq '' || $ind_alleles[1] eq '') {
					$geno_string = '';
				} else {
					$geno_string = $ind_alleles[0] . $ind_alleles[1];
				}

				my $l;
				my $m;
				if ($f0 eq $m0) {
					$l = $f0;
					$m = $f1;
				} else {
					$l = $f1;
					$m = $f0;
				}

				if ($geno_string =~ /$l$l/i || $geno_string =~ /$m$m/i) {
					push @progeny_codes, 'a';
				} elsif ($geno_string =~ /$l$m/i || $geno_string =~ /$m$l/i) {
					push @progeny_codes, 'ab';
				} elsif ($geno_string eq '') {
					push @progeny_codes, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
				} else {
					push @progeny_codes, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
					print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
				}
			}
		} elsif ($seg_type eq 'E') {
			$type = 'D2.15';
			if ($sex_specific eq 'F') {
				$stats{'excluded_loci'}[0]++;
				push @{$stats{'excluded_loci'}[1]}, $locus;
				next;
			}
			foreach my $ind (@progeny) {

				my @ind_geno = $ind->get_Genotypes(-marker => $locus);
				my @ind_alleles = $ind_geno[0]->get_Alleles();
				my $geno_string;
				if ($ind_alleles[0] eq '' || $ind_alleles[1] eq '') {
					$geno_string = '';
				} else {
					$geno_string = $ind_alleles[0] . $ind_alleles[1];
				}

				my $n = $f0;
				my $p;
				if ($m0 eq $n) {
					$p = $m1;
				} else {
					$p = $m0;
				}
				if ($geno_string =~ /$n$n/i) {
					push @progeny_codes, 'a';
				} elsif ($geno_string =~ /$n$p/i || $geno_string =~ /$p$n/i) {
					push @progeny_codes, 'ab';
				} elsif ($geno_string eq '') {
					push @progeny_codes, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
				} else {
					push @progeny_codes, '-';
					$ind_missing{$ind->unique_id()}++;
					$loc_missing{$locus}++;
					print ERR "Error: Check locus $locus, genotype $geno_string - ", $ind->unique_id(), "\n";
				}
			}
		} else {
			$stats{'strange_seg'}[0]++;
			push @{$stats{'strange_seg'}[1]}, $locus;
			$stats{'excluded_loci'}[0]++;
			push @{$stats{'excluded_loci'}[1]}, $locus;
			next;
		}

		$progeny_codes{$locus} = [$type, \@progeny_codes];

		#my $line = join("\t", $locus, $type, @progeny_codes);
		#push @lines1, $line;
		#print join("\t", $locus, $type, @progeny_codes), "\n";

		$stats{'added_loci'}[0]++;
		push @{$stats{'added_loci'}[1]}, $locus;


	}

}

print DUMP Dumper(\%ind_missing);
print DUMP Dumper(\%loc_missing);
print DUMP Dumper(\%progeny_codes_m);

close ERR;

# Build headers and write files

if ($outformat eq 'joinmap') {

	if ($outfile) {
		$outfile = $outfile . '.loc';
	} elsif (! $outfile && $infile) {
		$outfile = $infile . '.loc';
	} else {
		$outfile = 'joinmap.loc';
	}


	my %ind_discard = filter_inds(\@progeny, \%ind_missing, $ind_cutoff);
	my %loc_discard = filter_loci(\%progeny_codes, \%loc_missing, $loc_cutoff);

	print DUMP Dumper(\%progeny_codes);

	open(OUT, ">", $outfile) or die $!;

	print OUT 'name = ', $pop->name, "\n";
	print OUT 'popt = ', 'CP', "\n";
	print OUT 'nloc = ', $stats{'added_loci'}[0] - $loc_discard{'total'}, "\n";
	print OUT 'nind = ', $stats{'added_inds'}[0] - $ind_discard{'total'}, "\n";

	foreach my $locus (sort keys %progeny_codes) {
		next unless $progeny_codes{$locus}[0];
		if ($loc_discard{$locus} == 0) {
			print OUT join("\t", $locus, $progeny_codes{$locus}[0]), "\t";
			for (my $i = 0; $i < scalar(@{$progeny_codes{$locus}[1]}); $i++) {
				if ($ind_discard{$progeny[$i]->unique_id()} == 0) {
					print OUT $progeny_codes{$locus}[1][$i], "\t";
				}
			}
			print OUT "\n";
		}

	}

	print OUT "\n";
	print OUT "individual names:\n";
	foreach my $ind (@progeny) {
		print OUT $ind->unique_id(), "\n" unless $ind_discard{$ind->unique_id()} == 1;
	}
	close OUT;

}

if ($outformat eq 'mm') {

	if ($outfile) {
		$outfile = $outfile . '.raw';
	} elsif (! $outfile && $infile) {
		$outfile = $infile . '.raw';
	} else {
		$outfile = 'mm.raw';
	}

	if ($sex_specific eq 'F') {
		%progeny_codes = %progeny_codes_f;
	} elsif ($sex_specific eq 'M') {
		%progeny_codes = %progeny_codes_m;
	}

	my %ind_discard = filter_inds(\@progeny, \%ind_missing, $ind_cutoff);
	my %loc_discard = filter_loci(\%progeny_codes, \%loc_missing, $loc_cutoff);

	open(OUT, ">", $outfile) or die $!;

	print OUT 'data type f2 backcross', "\n";
	print OUT join(' ', scalar(@progeny) - $ind_discard{'total'}, $stats{'added_loci'}[0] - $loc_discard{'total'}, 1), "\n";

	foreach my $locus (sort keys %progeny_codes) {
		if ($loc_discard{$locus} == 0) {
			print OUT "*$locus ";
			for (my $i = 0; $i < scalar(@{$progeny_codes{$locus}[1]}); $i++) {
				if ($ind_discard{$progeny[$i]->unique_id()} == 0) {
					print OUT $progeny_codes{$locus}[1][$i];
				}
			}
			print OUT "\n";
		}
	}

	close OUT;

	my $base;
	if ($outfile =~ /(\w+)\.(\w+)/) {
		$base = $1;
	}
	open(IND, ">", $base . '.ind') or die $!;
	foreach my $ind (@progeny) {
		print IND $ind->unique_id(), "\n" unless $ind_discard{$ind->unique_id()} == 1;
	}
	close IND;


	if ($mapfile) {
		my %map = read_map($mapfile);
		my %new_map;
		foreach my $locus (sort keys %progeny_codes) {
			if ($loc_discard{$locus} == 0) {
				my $mapped = 0;
				foreach my $lg (keys %map) {
					if ($map{$lg}{$locus}) {
						$new_map{$lg}{$locus} = $map{$lg}{$locus};
						$mapped = 1;
						last;
					}
				}
				if ($mapped == 0) {
					$new_map{'0'}{$locus} = 0;
				}
			}
		}
		open(MAP, ">", $base . '.map');
		#print join("\n", sort keys %new_map);
		foreach my $lg (sort {$a <=> $b} keys %new_map) {
			foreach my $locus (keys %{$new_map{$lg}}) {
				print MAP join("\t", $lg, "*$locus"), "\n";
			}
		}
		close MAP;
	}

}


if ($outformat eq 'onemap') {

	if ($outfile) {
		$outfile = $outfile . '.txt';
	} elsif (! $outfile && $infile) {
		$outfile = $infile . '.txt';
	} else {
		$outfile = 'onemap.txt';
	}


	my %ind_discard = filter_inds(\@progeny, \%ind_missing, $ind_cutoff);
	my %loc_discard = filter_loci(\%progeny_codes, \%loc_missing, $loc_cutoff);

	print DUMP Dumper(\%progeny_codes);

	open(OUT, ">", $outfile) or die $!;

	print OUT join("\t", $stats{'added_inds'}[0] - $ind_discard{'total'}, $stats{'added_loci'}[0] - $loc_discard{'total'}), "\n";

	foreach my $locus (sort keys %progeny_codes) {
		next unless $progeny_codes{$locus}[0];
		if ($loc_discard{$locus} == 0) {
			my $code_string = '';
			print OUT join("\t", "*$locus", $progeny_codes{$locus}[0]), "\t";
			for (my $i = 0; $i < scalar(@{$progeny_codes{$locus}[1]}); $i++) {
				if ($ind_discard{$progeny[$i]->unique_id()} == 0) {
					$code_string = $code_string .  $progeny_codes{$locus}[1][$i] . ',';
				}
			}
			$code_string =~ s/\,$//;
			print OUT $code_string, "\n";
		}

	}

	open(INDS, ">", $outfile . '.inds') or die $!;
	foreach my $ind (@progeny) {
		print INDS $ind->unique_id(), "\n" unless $ind_discard{$ind->unique_id()} == 1;
	}
	close INDS;

}


#print DUMP Dumper(\%stats);

##### Subroutines #####

sub read_csv {

	# Read in a csv file of genotypes and store data in BioPerl objects

	my $file = $_[0];


	open(CSV, "<", $file) or die $!;

	my $pop = Bio::PopGen::Population->new();
	$pop->name($file);

	my @loci;
	my @prog;

	while(<CSV>) {
		last if $_ =~ /^,/;
		if ($_ =~ /^Sample/) {
			my @cols = split(',', $_);
			for (my $i = 1; $i < scalar(@cols); $i += 2) {
				$cols[$i] =~ s/\n//g;
				push @loci, $cols[$i];
			}
			#print join("\n", @loci), "\n";
		} else {
			my @cols = split(',', $_);
			my $samp_name = splice(@cols, 0, 1);
			my $ind = Bio::PopGen::Individual->new(-unique_id => $samp_name);
			foreach my $locus (@loci) {
				my @alleles = splice(@cols, 0, 2);
				foreach my $allele (@alleles) {
					$allele =~ s/[\n\r]//g;
					$allele = sprintf("%03s", $allele);
					if ($allele eq '000') {
						$allele = '';
					}
				}
				my $genotype = Bio::PopGen::Genotype->new(-marker_name => $locus, -alleles => \@alleles);
				$ind->add_Genotype($genotype);
			}
			$pop->add_Individual($ind);
			push @prog, $samp_name;
		}
	}

	return $pop, \@prog;
}

sub read_tsv {

	my $file = $_[0];
	my $f_id = $_[1];
	my $m_id = $_[2];

	open(TSV, "<", $file) or die $!;

	my $pop = Bio::PopGen::Population->new();
	$pop->name($file);

	my @loci;
	my $header = <TSV>;
	my @fields = split(/\s/, $header);
	chomp(my @ind_names = @fields[4..scalar(@fields)-1]);

	my %inds;
	while(<TSV>) {
		last if $_ =~ /^\s/;
		my ($locus, $par_gen, $count, $f, @genotypes) = split;
		if ($loc_cutoff) {
			my $prop = $count / scalar(@ind_names);
			next if $prop < $loc_cutoff;
		}
		my ($f_gen, $m_gen) = split('/', $par_gen);
		my @f_alleles = (substr($f_gen, 0, 1), substr($f_gen, 1, 1));
		my @m_alleles = (substr($m_gen, 0, 1), substr($m_gen, 1, 1));
		$inds{$f_id}{$locus} = \@f_alleles;
		$inds{$m_id}{$locus} = \@m_alleles;
		for (my $i = 0; $i < scalar(@genotypes); $i++) {
			my @alleles = (substr($genotypes[$i], 0, 1), substr($genotypes[$i], 1, 1));
			$inds{$ind_names[$i]}{$locus} = \@alleles;
		}
	}
	foreach my $samp_name (@ind_names, $f_id, $m_id) {
		my $ind = Bio::PopGen::Individual->new(-unique_id => $samp_name);
		foreach my $locus (keys %{$inds{$samp_name}}) {
			my $genotype = Bio::PopGen::Genotype->new(-marker_name => $locus, -alleles => $inds{$samp_name}{$locus});
			$ind->add_Genotype($genotype);
		}
		$pop->add_Individual($ind);

	}
	return $pop, \@ind_names;

}

sub read_joinmap {
	my $file = $_[0];
	my $f_id = $_[1];
	my $m_id = $_[2];

	open(LOC, "<", $file) or die $!;

	my $pop = Bio::PopGen::Population->new();
	$pop->name($file);

	my %data;
	my @loci;
	my $name_line = <LOC>;
	my $pop_type_line = <LOC>;
	my $num_loc_line = <LOC>;
	my $num_ind_line = <LOC>;
	# foreach my $line ($name_line, $pop_type_line, $num_loc_line, $num_ind_line) {
		# my ($left, $right) = split('=', $line);
		# $right =~ s/\s//g;
		# $line = $right;
	# }
	# $data{'name'} = $name_line;
	# $data{'pop_type'} = $pop_type_line;
	# $data{'num_loci'} = $num_loc_line;
	# $data{'num_inds'} = $num_ind_line;

	my $ind_parse = 0;
	my @inds;
	while(<LOC>) {
		next if $_ =~ /^\s/;
		if ($_ =~ /^individual/) {
			$ind_parse = 1;
			next;
		}
		if ($ind_parse == 1) {
			chomp;
			push @inds, $_;
			next;
		}
		chomp;
		my @fields = split;
		my $locus = shift @fields;
		my $seg_type = shift @fields;
		$data{$locus} = [$seg_type, \@fields];
		#print join("\t", $locus, $seg_type), "\n";
	}
	$data{'inds'} = \@inds;
	print DUMP Dumper(\@inds);

	for (my $i = 0; $i < scalar(@{$data{'inds'}}); $i++) {
		my $ind = Bio::PopGen::Individual->new(-unique_id => $data{'inds'}[$i]);
		foreach my $locus (keys %data) {
			next if $locus eq 'inds';
			my $geno_string = $data{$locus}[1][$i];
			my @alleles = split('', $geno_string);
			foreach my $allele (@alleles) {
				if ($allele eq '-') {
					$allele = '';
				}
			}
			my $genotype = Bio::PopGen::Genotype->new(-marker_name => $locus, -alleles => \@alleles);
			$ind->add_Genotype($genotype);
			#print join("\t", $ind->unique_id(), $locus, @alleles), "\n";
		}
		$pop->add_Individual($ind);
	}
	my $f_ind = Bio::PopGen::Individual->new(-unique_id => $f_id);
	my $m_ind = Bio::PopGen::Individual->new(-unique_id => $m_id);
	$Bio::PopGen::Genotype::BlankAlleles = '[\s-?]';

	foreach my $locus (keys %data) {
		next if $locus eq 'inds';
		my ($f_geno, $m_geno) = split('x', $data{$locus}[0]);
		my @f_alleles = (substr($f_geno, 1, 1), substr($f_geno, 2, 1));
		my @m_alleles = (substr($m_geno, 0, 1), substr($m_geno, 1, 1));
		my $f_gen = Bio::PopGen::Genotype->new(-marker_name => $locus, -alleles => \@f_alleles);
		my $m_gen = Bio::PopGen::Genotype->new(-marker_name => $locus, -alleles => \@m_alleles);
		$f_ind->add_Genotype($f_gen);
		$m_ind->add_Genotype($m_gen);
		#print join("\t", $f_ind->unique_id(), $locus, @f_alleles), "\n";
		#print join("\t", $m_ind->unique_id(), $locus, @m_alleles), "\n";

	}
	$pop->add_Individual($f_ind);
	$pop->add_Individual($m_ind);
	#print DUMP Dumper($pop);
	#foreach my $ind ($pop->get_Individuals(-unique_id => 'FamA_F')) {
	#	foreach my $gen ($ind->get_Genotypes) {
	#		my @alleles = $gen->get_Alleles();
	#		print join("\t", $ind->unique_id(), $gen->marker_name, @alleles), "\n";
	#	}
	#}
	return $pop, \@inds;
}

sub optimize {
	my $file = $_[0];
	my @inds = @{$_[1]};
	open(OPT, "<", $file) or die $!;
	my %ind_miss;
	my %loc_miss;
	my %lines;
	my $total_loci;
	while(<OPT>) {
		next unless $_ =~ /^\*/;
		chomp;
		my ($locus, $geno_string) = split;
		my $loc_miss;
		$locus =~ s/\*//;
		$lines{$locus} = $geno_string;
		for (my $i = 0; $i < length($geno_string); $i++) {
			my $score = substr($geno_string, $i, 1);
			if ($score eq '-') {
				$loc_miss{$locus}++;
				$ind_miss{$inds[$i]->unique_id()}++;
			}
		}
		$total_loci++;
	}
	my @keep;
	my @discarded_inds;
	foreach my $ind (@inds) {
		my $prop = $ind_miss{$ind->unique_id()} / $total_loci;
		if ($prop > (1 - $ind_cutoff)) {
			push @keep, 0;
			push @discarded_inds, $ind->unique_id(); # Add names to the array, not objects
		} else {
			push @keep, 1;
		}
		print join("\t", $ind->unique_id(), $ind_miss{$ind->unique_id()}, $total_loci, $prop), "\n";
	}
	open(OUT, ">", "opt.$file") or die $!;
	print OUT 'data type f2 backcross', "\n";
	print OUT join(' ', scalar(@inds) - scalar(@discarded_inds), $total_loci, 1), "\n";
	foreach my $locus (keys %lines) {
		my $prop = $loc_miss{$locus} / length($lines{$locus});
		next if $prop > (1 - $loc_cutoff);
		print OUT "*$locus ";
		for (my $i = 0; $i < length($lines{$locus}); $i++) {
			print OUT substr($lines{$locus}, $i, 1) if $keep[$i] == 1;
		}
		print OUT "\n";
	}

	my $base;
	if ($file =~ /(\w+)\.(\w+)/) {
		$base = $1;
	}
	open(IND, ">", $base . '.opt.ind') or die $!;
	foreach my $ind (@inds) {
		my $name = $ind->unique_id();
		if (grep(/$name/, @discarded_inds)) {
			next;
		} else {
			print IND $name, "\n";
		}
	}
	close IND;
}

sub filter_inds {
	my @inds = @{$_[0]};
	my %missing_inds = %{$_[1]};
	my $cutoff = $_[2];

	my %discard;
	my $discarded = 0;
	foreach my $ind (@inds) {
		my $prop_missing;
		if (! $missing_inds{$ind->unique_id()}) {
			$prop_missing = 0;
		} else {
			$prop_missing = $missing_inds{$ind->unique_id()} / $stats{'added_loci'}[0];
		}
		#print $ind->unique_id(), "\t", $prop_missing, "\n";
		if ($prop_missing > (1 - $cutoff)) {
			$discard{$ind->unique_id()} = 1;
			$discarded++;
		} else {
			$discard{$ind->unique_id()} = 0;
		}
	}
	$discard{'total'} = $discarded;
	return %discard;
}

sub filter_loci {
	my %progeny_codes = %{$_[0]};
	my %missing_loc = %{$_[1]};
	my $cutoff = $_[2];

	my %discard;
	my $discarded = 0;
	foreach my $locus (keys %progeny_codes) {
		my $prop_missing;
		if (! $missing_loc{$locus}) {
			$prop_missing = 0;
		} else {
			$prop_missing = $missing_loc{$locus} / scalar(@{$progeny_codes{$locus}[1]});
		}
		if ($prop_missing > (1 - $cutoff)) {
			$discard{$locus} = 1;
			$discarded++;
		} else {
			$discard{$locus} = 0;
		}
	}
	$discard{'total'} = $discarded;
	return %discard;
}

sub read_map {
	my $mapfile = shift;
	open(MAP, "<", $mapfile) or die $!;
	my %lgs;
	my $group;
	while(<MAP>) {
		next if $_ =~ /^\s/;
		if ($_ =~ /^group/) {
			if ($_ =~ /^group (\w+)/) {
				$group = $1;
			}
			next;
		}
		my ($locus, $pos) = split;
		chomp($pos);
		$lgs{$group}{$locus} = $pos;

	}
	close MAP;

	return %lgs;
}

__END__

=head1 NAME

map_convert.pl

=head1 SYNOPSIS

perl map_convert.pl

Options:
     -i     in_file
     -n     in_format
     -o		out_file
     -u		out_format
     -s		sex_specific_parent
     -F		female_id
     -M		male_id
     -l		locus_cutoff
	 -c		ind_cutoff

=head1 OPTIONS

=over 8

=item B<-i>

Name of input file to be converted. Currently supports:

csv (Comma separated Excel genotypes)
tsv (Stacks raw tab-separated genotype file *_genotypes.tsv)
joinmap (Joinmap format '.loc', with individual names included)

=item B<-n>

See above. Currently supported options are 'csv', 'tsv', 'joinmap'

=item B<-o>

Name of outfile - only requires a basename, not a file extension.

=item B<-u>

Format of outfile. Currently supports:

joinmap (Joinmap format '.loc', with individual names included)
mm (mapmaker '.raw' backcross output, suitable for input into MapMaker, R/qtl, or CarthaGene)
onemap (OneMap format)

=item B<-s>

Write sex-specific output for the parent specified. Accepts 'F' or 'M'.

=item B<-F>

ID of female parent. This is required for all input formats.

=item B<-M>

ID of male parent. This is required for all input formats.

=item B<-l>

Cutoff for proportion of scored individuals to keep a locus. For example, to only keep loci with 95% of the progeny scored, use 0.95.

=item B<-c>

Cutoff for proportion of scored loci to keep an individual. For example, to only keep individuals with 95% of their loci genotyped, use 0.95.

=back

=head1 DESCRIPTION

B<map_convert.pl> converts linkage mapping data for outcrossed species between various formats. It can also do some limited data filtering, if you ask nicely.

=cut
