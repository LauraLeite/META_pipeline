#!/usr/bin/perl
######################################################################
#######################################################################
# Copyright 2011 Fundação Oswaldo Cruz
# Author: Eric Aguiar and Laura Rabelo
# template.pl is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
#
# template.pl is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with template.pl (file: COPYING).
# If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################
#######################################################################


use Getopt::Long;
use Bio::SeqIO;

my $usage = "

$0 -i <input file> 
$0 -h

-i <input file>		: Fasta Input file
-h			: Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.

my $inputFile;
my $outputFile;
my %hash;
GetOptions ("i=s" => \$inputFile,
			"h!" => \$help);

if ($help) {
	die $usage;
}

if (not(defined($inputFile))) {
	die "\nGive an input file name",$usage;
}
	

our $sum=0;
my $div;
our $nContigs;
my $nMaior100b=0;
my $totBaseMaior100b=0;
my $nMaior200b=0;
my $totBaseMaior200b=0;
my $nMaior500b=0;
my $totBaseMaior500b=0;
my $nMaior1k=0;
my $totBaseMaior1k=0;
my $nMaior2k=0;
my $totBaseMaior1k=0;
my $nMaior2k=0;
my $totBaseMaior1k=0;
my $nMaior3k=0;
my $totBaseMaior3k=0;
my $nMaior4k=0;
my $totBaseMaior4k=0;
my $nMaior5k=0;
my $totBaseMaior5k=0;
my $nMaior6k=0;
my $totBaseMaior6k=0;
my $nMaior7k=0;
my $totBaseMaior7k=0;
my $nMaior8k=0;
my $totBaseMaior8k=0;
my $nMaior9k=0;
my $totBaseMaior9k=0;
my $nMaior10k=0;
my $totBaseMaior10k=0;

my @v;
my $j=0;
readFile();
our $average = $sum / $nContigs;
calc();

sub readFile{
	
	#print "\nLoading Contigs... \n";
	my $seqio_object = Bio::SeqIO->new(-file => $inputFile);
	while($seq_object = $seqio_object->next_seq){	
		
		$nContigs++;
		$hash{$seq_object->id}=length($seq_object->seq);
		$totBase+=length($seq_object->seq);
		$hash{$seq_object->id}=length($seq_object->seq);
		$v[$j]=length($seq_object->seq);
		$j++;
		
		if (length($seq_object->seq)> 100){
			$nMaior100b++;	
			$totBaseMaior100b+=length($seq_object->seq);
				if (length($seq_object->seq)> 200){
                        		$nMaior200b++;
                        		$totBaseMaior200b+=length($seq_object->seq);
                                        if (length($seq_object->seq)> 500){
                        		$nMaior500b++;
                        		$totBaseMaior500b+=length($seq_object->seq);
						if (length($seq_object->seq)> 1000){
                                        	$nMaior1k++;
                                        	$totBaseMaior1k+=length($seq_object->seq);
							if (length($seq_object->seq)> 2000){
							$nMaior2k++;
							$totBaseMaior2k+=length($seq_object->seq);
								if (length($seq_object->seq)> 3000){
								$nMaior3k++;
								$totBaseMaior3k+=length($seq_object->seq);
									if (length($seq_object->seq)> 4000){
									$nMaior4k++;
									$totBaseMaior4k+=length($seq_object->seq);
										if (length($seq_object->seq)> 5000){
										$nMaior5k++;
										$totBaseMaior5k+=length($seq_object->seq);
											if (length($seq_object->seq)> 6000){
											$nMaior6k++;
											$totBaseMaior6k+=length($seq_object->seq);
												if (length($seq_object->seq)> 7000){
												$nMaior7k++;
												$totBaseMaior7k+=length($seq_object->seq);
													if (length($seq_object->seq)> 8000){
													$nMaior8k++;
													$totBaseMaior8k+=length($seq_object->seq);
														if (length($seq_object->seq)> 9000){
														$nMaior9k++;
														$totBaseMaior9k+=length($seq_object->seq);
															if (length($seq_object->seq)> 10000){
															$nMaior10k++;
															$totBaseMaior10k+=length($seq_object->seq);
															}
														}
													}
												}
											}
										}
									}
								}
							}
                                               }
					}
				}
		}	
		$sum+=length($seq_object->seq);
	}
	#print "Close file\n";
	$div = $sum /2;
	$div80 = $sum/3.2;
	#print "Sum: $sum\n Sum /2 :$div\n";
}
sub calcMediana{
	@v2= sort{$a <=> $b} @v;
	$pm=$nContigs / 2;
	$parImpar=$nContigs %2;
#	print "@v2\n $pm -$nContigs -  $parImpar\n";
	if ($nContigs%2 == 0 ){
#		print "mediana:".$v2[$pm]."\n"; 
		return $v2[$pm];
	}else{
		$pm = int($pm) +1;
		return $v2[$pm];
	}	
	

}
sub hashValueAscendingNum {
   $hash{$a} <=> $hash{$b};
}
sub hashValueDescendingNum {
   $hash{$b} <=> $hash{$a};
}

sub calc{
	
#print "\nGRADES IN DESCENDING NUMERIC ORDER:\n";
my $value=0;
$average = $sum / $nContigs;
my $largest;
foreach $key (sort hashValueDescendingNum (keys(%hash))) {
  #print "\t$hash{$key} \t\t $key\n";
   if ($value==0){
      $largest = $hash{$key};
   }
  # print "$hash{$key} \n";
   $value+=$hash{$key};
   if ($value > $div){
	   my $perc100b=($nMaior100b*100)/$nContigs;
	   $perc100b=sprintf("%.2f", $perc100b);
	   my $perc200b=($nMaior200b*100)/$nContigs;
           $perc200b=sprintf("%.2f", $perc200b);
           my $perc500b=($nMaior500b*100)/$nContigs;
           $perc100b=sprintf("%.2f", $perc100b);
	   my $perc1k=($nMaior1k*100)/$nContigs;
           $perc1k=sprintf("%.2f", $perc1k);
           my $perc2k=($nMaior2k*100)/$nContigs;
           $perc2k=sprintf("%.2f", $perc2k);
           my $perc3k=($nMaior3k*100)/$nContigs;
           $perc3k=sprintf("%.2f", $perc3k);
           my $perc4k=($nMaior4k*100)/$nContigs;
           $perc4k=sprintf("%.2f", $perc4k);
           my $perc5k=($nMaior5k*100)/$nContigs;
           $perc5k=sprintf("%.2f", $perc5k);
           my $perc6k=($nMaior6k*100)/$nContigs;
           $perc6k=sprintf("%.2f", $perc6k);
           my $perc7k=($nMaior7k*100)/$nContigs;
           $perc7k=sprintf("%.2f", $perc7k);
           my $perc8k=($nMaior8k*100)/$nContigs;
           $perc8k=sprintf("%.2f", $perc8k);
           my $perc9k=($nMaior9k*100)/$nContigs;
           $perc9k=sprintf("%.2f", $perc9k);
           my $perc10k=($nMaior10k*100)/$nContigs;
           $perc10k=sprintf("%.2f", $perc10k);
		
	
	   print "$inputFile\n";

	   print "Number of contigs\t".$nContigs."\nNumber of bases in all contigs\t".$totBase."\nN50\t".$hash{$key}."\nLongest contig\t".$largest."\nMedian contig size\t".calcMediana()."\nNumber of contigs > 100 b\t".$nMaior100b."\nTotal of bases into contigs > 100 b\t".$totBaseMaior100b."\n% contigs > 100 b\t".$perc100b."\nNumber of contigs > 200 b\t".$nMaior200b."\nTotal of bases into contigs > 200 b\t".$totBaseMaior200b."\n% contigs > 200 b\t".$perc200b."\nNumber of contigs > 500 b\t".$nMaior500b."\nTotal of bases into contigs > 500 b\t".$totBaseMaior500b."\n% contigs > 500 b\t".$perc500b."\nNumber of contigs > 1 kb\t".$nMaior1k."\nTotal of bases into contigs > 1 kb\t".$totBaseMaior1k."\n% contigs > 1 kb\t".$perc1k."\nNumber of contigs > 2 kb\t".$nMaior2k."\nTotal of bases into contigs > 2 kb\t".$totBaseMaior2k."\n% contigs > 2 kb\t".$perc2k."\nNumber of contigs > 3 kb\t".$nMaior3k."\nTotal of bases into contigs > 3 kb\t".$totBaseMaior3k."\n% contigs > 3 kb\t".$perc3k."\nNumber of contigs > 4 kb\t".$nMaior4k."\nTotal of bases into contigs > 4 kb\t".$totBaseMaior4k."\n% contigs > 4 kb\t".$perc4k."\nNumber of contigs > 5 kb\t".$nMaior5k."\nTotal of bases into contigs > 5 kb\t".$totBaseMaior5k."\n% contigs > 5 kb\t".$perc5k."\nNumber of contigs > 6 kb\t".$nMaior6k."\nTotal of bases into contigs > 6 kb\t".$totBaseMaior6k."\n% contigs > 6 kb\t".$perc6k."\nNumber of contigs > 7 kb\t".$nMaior7k."\nTotal of bases into contigs > 7 kb\t".$totBaseMaior7k."\n% contigs > 7 kb\t".$perc7k."\nNumber of contigs > 8 kb\t".$nMaior8k."\nTotal of bases into contigs > 8 kb\t".$totBaseMaior8k."\n% contigs > 8 kb\t".$perc8k."\nNumber of contigs > 9 kb\t".$nMaior9k."\nTotal of bases into contigs > 9 kb\t".$totBaseMaior9k."\n% contigs > 9 kb\t".$perc9k."\nNumber of contigs > 10 kb\t".$nMaior10k."\nTotal of bases into contigs > 10 kb\t".$totBaseMaior10k."\n% contig9s > 10 kb\t".$perc10k."\n";
   		die;
   }
   
  
}

}
