#!/usr/bin/perl

##########################################################################################################
## Filter BLAST tabular output by minimum coverage 
## You may distribute this module under the same terms as perl itself
## Laura Rabelo Leite
## GGBC Fiocruz-MG - 11_05_2015
##########################################################################################################

use strict;
use Getopt::Long;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use Data::Dumper;
use diagnostics;

my $help;
my $blast_m8;
my $fasta_file;
my $file_output;
my $cut_off;


GetOptions(
	'i=s'    => \$blast_m8,
        'f=s'    => \$fasta_file,
	'o=s'    => \$file_output,
	'c=s'    => \$cut_off,
	'help|?' => \$help
);


main();



########################
# MAIN 
########################

sub main{
	validate_parameters();
	print "Reading file $blast_m8 ...\n";
	read_fasta();
	read_blast();	
	process_m8();
	print "\nFile generated: $file_output\n";
	print "Process finished!\n";
}




########################
## VALIDATION
#########################


sub validate_parameters {
	my $allExists  = 1;
	my $fileExists = 1;

	if ( defined $help ) {
		print usage();
		exit 0;
	}

	unless ( defined $blast_m8 ) {
		$allExists = 0;
	}
	
	unless ( defined $fasta_file ) {
		$allExists = 0;
	}
	
	unless ( defined $cut_off ) {
		$allExists = 0;
	}

        unless ( defined $file_output ) {
                $allExists = 0;
        }

        unless ($allExists) {
                print usage();
                exit 0;
        }

	unless ($fileExists) {
		print STDERR "Program execution aborted.\n";
		exit 0;
	}

}


sub usage {
	my $usage = <<FOO;
Usage:
	perl $0 -i blast_m8_file -f fasta_file -o output_name -c coverage_cut_off
	blast_m8_file		Blast output in tabular format
	fasta_file		Fasta file correspondent to blast output
	output_name		Name of output file
	coverage_cut_off	Percentage of query bases covering in the reference file
FOO
	return $usage;

}



########################
### READ INPUT FILES
##########################

my %fasta_hash=();
my %blast_hash=();


sub read_fasta{
	
	my $seqio_object = Bio::SeqIO->new(-file => $fasta_file);
	while(my $seq_object = $seqio_object->next_seq){	
		$fasta_hash{$seq_object->id}=length($seq_object->seq);
	}

#	print Dumper(\%fasta_hash);

}



sub read_blast{
        open(IN2, $blast_m8);

	while(<IN2>){
#		print $_."\n";
		chomp;
        	my @line = split(/\t/,$_);
#		print "$line[3]\n";
		my @value_blast_hash;
		push(@value_blast_hash, $line[3]);
		my $ref_line = \@line;
		push(@value_blast_hash, $ref_line);
#		print "@value_blast_hash\n";
#		print $line[0]."\n";
		my $ref_value = \@value_blast_hash;
		$blast_hash{$line[0]}=$ref_value;
        	}
        close(IN2);
#       print Dumper(\%blast_hash);
}


########################
#### FILTER BLAST FILE
###########################


sub process_m8{
	open(OUT, ">", $file_output);
	my @value;
	my @length;
	foreach my $key1 (sort keys %fasta_hash){
			foreach my $key2 (sort keys %blast_hash){
				if ($key1 eq $key2){
					my $length_seq = $fasta_hash{$key1};
					my $alignment = $blast_hash{$key2}[0];
					my $percentage = ($alignment*100)/($length_seq+1);
#					print "$percentage\=\($alignment\*100\)\/\($length_seq\+1\)";
					if($percentage >= $cut_off){
						my @reference = @{$blast_hash{$key2}[1]};
						foreach(@reference){
							print OUT "$_\t";
						}
						print OUT "\n";
					}
				}
			}	
	}
	close(OUT);
}


