#!/usr/bin/perl

######################################################################
########################################################################
## Copyright 2016 Fundação Oswaldo Cruz
## Author: Laura Rabelo Leite
## parse_gene.pl is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
##
## META_pipeline.pl is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY;
## without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with get_kmers.pl (file: COPYING).
## If not, see <http://www.gnu.org/licenses/>.
##
########################################################################
########################################################################

use Bio::SeqIO;
use File::Basename;

my $in = shift or die "Must provide an input file\n";
my $out = shift or die "Must provide an output name \n";

open(OUT, ">".$out);
open(FILE, $in) or die "Can't open $in\n";
	while (my $line =<FILE>){
	chomp ($line);
	$line =~ s/,/\t/g;
	my @config = split(/\t/, $line);
	my $new_line="$config[0]\t$config[1]\t$config[2]\t$config[3]\t$config[18]\t$config[19]\n";
	$new_line=~ s/gene_id=/gene_id_/g;
	$new_line=~ s/length=//g;
	print OUT $new_line;
	}
close(FILE);
close(OUT);

