#!/usr/bin/perl

######################################################################
#######################################################################
# Copyright 2016 Fundação Oswaldo Cruz
# Author: Laura Rabelo Leite
# META_pipeline.pl is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
#
# META_pipeline.pl is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with META_pipeline.pl (file: COPYING).
# If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################
#######################################################################

#TUCUNARE

use Bio::SeqIO;
use File::Basename;

my $in = shift or die "Must provide an input directory name\nusage: META_pipeline.pl <input directory> <config file> <output directory>\n";
my $param = shift or die "Must provide an config file that contain the sequence type (PE for paired-end or SE for single-end) and program parameters (one program per line. Settable programs: trimmomatic, spades, gmhmm, diamond, kanalyze, megablast, blastm8_filter, bowtie2. eg. diamond: -e 0.1 --id 50)\nusage: META_pipeline.pl <input directory> <config file> <output directory>\n";
my $out = shift or die "Must provide an output directory name\nusage: META_pipeline.pl <input directory> <config file> <output directory>\n";
my @list = <$in>;
@list = glob("*.fastq");


##Quality control (QC) and assembly

##Getting information from config file
my %config_p = ();
open(FILE, $param) or die "Can't open $param\n";
	while (my $line =<FILE>){
	chomp ($line);
	my @config = split(/\t/, $line);
	my $key=$config[0];
	my $val=$config[1];
	$config_p{$key}=$val;
	}
close(FILE);

#while (($key, $value) = each (%config_p))
#{
#  print "$key\t$config_p{$key}\n";
#}


foreach (@list) {

	if (exists($config_p{'trimmomatic'})){
	#QC with configuration
		my $config_line = $config_p{'trimmomatic'};
		if (grep { $_ =~ 'PE' } values %config_p ) {
			$config_line =~ s/PE//g;
			#QC for paired-end sequences (PE)
			print scalar localtime()."\n";
			print "trimmomatic-0.32.jar PE -threads 40 $list[0] $list[1] $out/output_1.fastq $out/output_1.unp.fastq $out/output_2.fastq $out/output_2.unp.fastq $config_line\n";
			`trimmomatic-0.32.jar PE -threads 40 $list[0] $list[1] $out/output_1.fastq $out/output_1.unp.fastq $out/output_2.fastq $out/output_2.unp.fastq $config_line &> trimmomatic.log`;
			if (exists($config_p{'spades'})){
				my $config_line2 = $config_p{'spades'};
				##Metagenome assembly with configuration
				print scalar localtime()."\n";
				print "spades.py -1 $out/output_1.fastq -2 $out/output_2.fastq -s  $out/output_all.unp.fastq -o  $out/spades $config_line2\n";
				`spades.py -1 $out/output_1.fastq -2 $out/output_2.fastq -s  $out/output_all.unp.fastq -o  $out/spades $config_line2`;
				#fastq to fasta conversion
				print scalar localtime()."\n";
				print "seqtk fq2fa $out/output_1.fastq >> $out/reads.fasta\nseqtk fq2fa $out/output_2.fastq >> $out/reads.fasta\nseqtk fq2fa $out/output_all.unp.fastq >> $out/reads.fasta\n";
				`seqtk fq2fa $out/output_1.fastq >> $out/reads.fasta`;
				`seqtk fq2fa $out/output_2.fastq >> $out/reads.fasta`;
				`seqtk fq2fa $out/output_all.unp.fastq >> $out/reads.fasta`;
			}else{
				##Metagenome assembly  without configuration
				print scalar localtime()."\n";
				print "spades.py -t 40 --cov-cutoff auto --phred-offset 33 -1 $out/output_1.fastq -2 $out/output_2.fastq -s  $out/output_all.unp.fastq --careful -k 105,115,127 -o  $out/spades\n";
				`spades.py -t 40 --cov-cutoff auto --phred-offset 33 -1 $out/output_1.fastq -2 $out/output_2.fastq -s  $out/output_all.unp.fastq --careful -k 105,115,127 -o  $out/spades`;
				#fastq to fasta conversion
				print scalar localtime()."\n";
				print "seqtk fq2fa $out/output_1.fastq >> $out/reads.fasta\nseqtk fq2fa $out/output_2.fastq >> $out/reads.fasta\nseqtk fq2fa $out/output_all.unp.fastq >> $out/reads.fasta\n";
				`seqtk fq2fa $out/output_1.fastq >> $out/reads.fasta`;
				`seqtk fq2fa $out/output_2.fastq >> $out/reads.fasta`;
				`seqtk fq2fa $out/output_all.unp.fastq >> $out/reads.fasta`;
			}
		}
		else{
			#QC for single-end sequenced (SE)
			$config_line =~ s/SE//g;
			print scalar localtime()."\n";
			print "trimmomatic-0.32.jar SE -threads 40 $list[0] $out/output_SE.fastq $config_line\n";		
			`trimmomatic-0.32.jar SE -threads 40 $list[0] $out/output_SE.fastq $config_line  &> trimmomatic.log`;
			if (exists($config_p{'spades'})){
				my $config_line2 = $config_p{'spades'};
				print scalar localtime()."\n";
				print "spades.py -s $out/output_SE.fastq -o $out/spades $config_line2\n";
				#Assembly for SE  with configuration
				`spades.py -s $out/output_SE.fastq -o $out/spades $config_line2`;
				#fastq to fasta conversion
				`seqtk fq2fa $out/output_SE.fastq >> $out/reads.fasta`;
			}else{
				#Assembly for SE  without configuration
				print scalar localtime()."\n";
				print "spades.py -s $out/output_SE.fastq -o $out/spades -t 40 -k 105,115,127\n";
				`spades.py -s $out/output_SE.fastq -o $out/spades -t 40 -k 105,115,127`;
				#fastq to fasta conversion
				print scalar localtime()."\n";
				print "seqtk fq2fa $out/output_SE.fastq >> $out/reads.fasta\n";
				`seqtk fq2fa $out/output_SE.fastq >> $out/reads.fasta`;
			}
		}
	}else {
	#QC without configuration
		if (grep { $_ == 'PE' } values %config_p) {
			#QC for paired-end sequences (PE)
			print scalar localtime()."\n";
			print "trimmomatic-0.32.jar PE -threads 40 -phred33 -trimlog trimmomatic.log  $list[0] $list[1] $out/output_1.fastq $out/output_1.unp.fastq $out/output_2.fastq $out/output_2.unp.fastq SLIDINGWINDOW:4:15 MINLEN:100\n";
			`trimmomatic-0.32.jar PE -threads 40 -phred33 -trimlog trimmomatic.log  $list[0] $list[1] $out/output_1.fastq $out/output_1.unp.fastq $out/output_2.fastq $out/output_2.unp.fastq SLIDINGWINDOW:4:15 MINLEN:100 &> trimmomatic.log`;
			`cat  $out/output_1.unp.fastq  $out/output_2.unp.fastq >  $out/output_all.unp.fastq`;
			##Metagenome assembly
			print scalar localtime()."\n";
			print "spades.py -t 40 --cov-cutoff auto --phred-offset 33 -1 $out/output_1.fastq -2 $out/output_2.fastq -s  $out/output_all.unp.fastq --careful -k 105,115,127 -o  $out/spades\n";
			`spades.py -t 40 --cov-cutoff auto --phred-offset 33 -1 $out/output_1.fastq -2 $out/output_2.fastq -s  $out/output_all.unp.fastq --careful -k 105,115,127 -o  $out/spades`;
			#fastq to fasta conversion
			print scalar localtime()."\n";
			print "seqtk fq2fa $out/output_1.fastq >> $out/reads.fasta\nseqtk fq2fa $out/output_2.fastq >> $out/reads.fasta\nseqtk fq2fa $out/output_all.unp.fastq >> $out/reads.fasta\n";
			`seqtk fq2fa $out/output_1.fastq >> $out/reads.fasta`;
			`seqtk fq2fa $out/output_2.fastq >> $out/reads.fasta`;
			`seqtk fq2fa $out/output_all.unp.fastq >> $out/reads.fasta`;
		}elsif (grep { $_ == 'SE' } values %config_p){
			#QC for single-end sequenced (SE)
			print scalar localtime()."\n";
			print "trimmomatic-0.32.jar SE -phred33 $list[0] $out/output_SE.fastq SLIDINGWINDOW:4:15  MINLEN:100\n";
			`trimmomatic-0.32.jar SE -phred33 $list[0] $out/output_SE.fastq SLIDINGWINDOW:4:15  MINLEN:100 &> trimmomatic.log`;
			#Assembly for SE
			print scalar localtime()."\n";
			print "spades.py -s $out/output_SE.fastq -o $out/spades -t 40 -k 105,115,127\n";
			`spades.py -s $out/output_SE.fastq -o $out/spades -t 40 -k 105,115,127`;
			#fastq to fasta conversion
			print scalar localtime()."\n";
			print "seqtk fq2fa $out/output_SE.fastq >> $out/reads.fasta\n";
			`seqtk fq2fa $out/output_SE.fastq >> $out/reads.fasta`;
		}
	}

	#CDS prediction
	#create CDS directory
	`mkdir CDS_annot`;
	#scaffolds statistic
	`perl calcN50.pl -i $out/spades/scaffolds.fasta > $out/spades/scaffolds.stat`;
	if (exists($config_p{'gmhmmp'})){
		print scalar localtime()."\n";
		print "gmhmmp $out/spades/scaffolds.fasta -m MetaGeneMark_v1.mod -o $out/CDS_annot/gmhmm -f G -a -d $config_line3\n";
		my $config_line3 = $config_p{'gmhmmp'};
		`gmhmmp $out/spades/scaffolds.fasta -m MetaGeneMark_v1.mod -o $out/CDS_annot/gmhmm -f G -a -d $config_line3`;
		`perl aa_from_gff.pl < $out/CDS_annot/gmhmm > $out/CDS_annot/gmhmm.faa`;

	}
	else{
		print scalar localtime()."\n";
		print "gmhmmp $out/spades/scaffolds.fasta -m MetaGeneMark_v1.mod -o $out/CDS_annot/gmhmm -f G -a -d\n";
		`gmhmmp $out/spades/scaffolds.fasta -m MetaGeneMark_v1.mod -o $out/CDS_annot/gmhmm -f G -a -d`;
		`perl aa_from_gff.pl < $out/CDS_annot/gmhmm > $out/CDS_annot/gmhmm.faa`;
	}


	#Functional annotation (based on the predicted protein set)
	if (exists($config_p{'diamond'})){
		my $config_line4 = $config_p{'diamond'};
		print scalar localtime()."\n";
		print "diamond blastp -d ueko -q $out/CDS_annot/gmhmm.faa -a $out/CDS_annot/gmhmm.kegg.daa -t $out/CDS_annot $config_line4\n";
		`diamond blastp -d ueko -q $out/CDS_annot/gmhmm.faa -a $out/CDS_annot/gmhmm.kegg.daa -t $out/CDS_annot $config_line4`;

		`diamond view -a $out/CDS_annot/gmhmm.kegg.daa -o $out/CDS_annot/gmhmm.kegg.m8`;

		#Filtration
		`cat $out/CDS_annot/gmhmm.kegg.m8 | perl evfilter -e 1e-5 | perl alnfilter -l 50 > $out/CDS_annot/filt1_gmhmm.kegg.m8`;
		`python BlastBestHit.py $out/CDS_annot/filt1_gmhmm.kegg.m8 > $out/CDS_annot/filt_gmhmm.kegg.m8`;
		`cut -f 1 -d "|" $out/CDS_annot/filt_gmhmm.kegg.m8 > $out/CDS_annot/interprot.m8`;


	}
	else{
		print scalar localtime()."\n";
		print "diamond blastp -p 40 -d ueko -q $out/CDS_annot/gmhmm.faa -e 0.00001 --id 30 -a $out/CDS_annot/gmhmm.kegg.daa -t $out/CDS_annot\n";
		`diamond blastp -p 40 -d ueko -q $out/CDS_annot/gmhmm.faa -e 0.00001 --id 30 -a $out/CDS_annot/gmhmm.kegg.daa -t $out/CDS_annot`;
		`diamond view -a $out/CDS_annot/gmhmm.kegg.daa -o $out/CDS_annot/gmhmm.kegg.m8`;

		`cat $out/CDS_annot/gmhmm.kegg.m8 | perl evfilter -e 1e-5 | perl alnfilter -l 50 > $out/CDS_annot/filt1_gmhmm.kegg.m8`;
		`python BlastBestHit.py $out/CDS_annot/filt1_gmhmm.kegg.m8 > $out/CDS_annot/filt_gmhmm.kegg.m8`;
		`cut -f 1 -d "|" $out/CDS_annot/filt_gmhmm.kegg.m8 > $out/CDS_annot/interprot.m8`;

	}



	#Taxonomic annotation (similarity based and LCA)
	#MEGABLAST
	`mkdir $out/TAX_annot`;
	if (exists($config_p{'megablast'})){
		my $config_line5 = $config_p{'megablast'};
		print scalar localtime()."\n";
		print "megablast -i $out/reads.fasta -d nt -o $out/TAX_annot/MEGABLAST_NT.txt $config_line5\n";
		`megablast -i $out/reads.fasta -d nt -o $out/TAX_annot/MEGABLAST_NT.txt $config_line5`;
	}
	else{
		print scalar localtime()."\n";
		print "megablast -i $out/reads.fasta -d nt -e 1e05 -p 90 -a 40 -m 8 -b 10 -o $out/TAX_annot/MEGABLAST_NT.txt\n";
		`megablast -i $out/reads.fasta -d nt -e 1e05 -p 90 -a 40 -m 8 -b 10 -o $out/TAX_annot/MEGABLAST_NT.txt`;
	}
	#BLAST filter
	if (exists($config_p{'blastm8_filter'})){
		print scalar localtime()."\n";
		print "perl coverage_blastm8_filter.pl -i $out/TAX_annot/MEGABLAST_NT.txt -f $out/reads.fasta -o $out/TAX_annot/FILT_MEGABLAST_NT.txt  $config_line6\n";
		my $config_line6 = $config_p{'blastm8_filter'};
		`perl coverage_blastm8_filter.pl -i $out/TAX_annot/MEGABLAST_NT.txt -f $out/reads.fasta -o $out/TAX_annot/FILT_MEGABLAST_NT.txt  $config_line6`;
	}
	else{
		print scalar localtime()."\n";
		print "perl coverage_blastm8_filter.pl -i $out/TAX_annot/MEGABLAST_NT.txt -f $out/reads.fasta -o $out/TAX_annot/FILT_MEGABLAST_NT.txt  -c 90\n";
		`perl coverage_blastm8_filter.pl -i $out/TAX_annot/MEGABLAST_NT.txt -f $out/reads.fasta -o $out/TAX_annot/FILT_MEGABLAST_NT.txt  -c 90`;
	}
	print scalar localtime()."\n";
	print "blast2lca -i $out/TAX_annot/FILT_MEGABLAST_NT.txt -a2t /sto3data-2/vale/test_META2/nucl_acc2tax-Nov2016.abin --accessionTags -f BlastTab -m BlastN -o $out/TAX_annot/FILT_MEGABLAST_NT-taxonomy.txt";
	`blast2lca -i $out/TAX_annot/FILT_MEGABLAST_NT.txt -a2t /sto3data-2/vale/test_META2/nucl_acc2tax-Nov2016.abin --accessionTags -f BlastTab -m BlastN -o $out/TAX_annot/FILT_MEGABLAST_NT-taxonomy.txt`;
#	versão para nt que ainda tem gi
#print "blast2lca -dict /sto3data-2/vale/databases/gi_taxid_nucl.bin -levels=superkingdom:phylum:class:family:genus -nodes /usr/local/bioinformatics/Blast2lca-master/nodes.dmp -names /usr/local/bioinformatics/Blast2lca-master/names.dmp $out/TAX_annot/FILT_MEGABLAST_NT.txt >  $out/TAX_annot/FILT_MEGABLAST_NT.taxonomy\n";
#`blast2lca -dict /sto3data-2/vale/databases/gi_taxid_nucl.bin -levels=superkingdom:phylum:class:family:genus -nodes /usr/local/bioinformatics/Blast2lca-master/nodes.dmp -names /usr/local/bioinformatics/Blast2lca-master/names.dmp $out/TAX_annot/FILT_MEGABLAST_NT.txt >  $out/TAX_annot/FILT_MEGABLAST_NT.taxonomy`;

	#Linking taxonomy and function
	#Map reads to contigs with bowtie
	print scalar localtime()."\n";
	print "bowtie2-build $out/spades/scaffolds.fasta $out/spades/scaffolds\n";
	`bowtie2-build $out/spades/scaffolds.fasta $out/spades/scaffolds`;

	if (exists($config_p{'bowtie2'})){
		print scalar localtime()."\n";
		print "bowtie2-align-s -x $out/spades/scaffolds -1 $out/output_1.fastq -2 $out/output_2.fastq -U $out/output_all.unp.fastq -S $out/TAX_annot/bowtie2_scaffolds.sam $config_line7\n";
		my $config_line7 = $config_p{'bowtie2'};
		`bowtie2-align-s -x $out/spades/scaffolds -1 $out/output_1.fastq -2 $out/output_2.fastq -U $out/output_all.unp.fastq -S $out/TAX_annot/bowtie2_scaffolds.sam $config_line7`;
	}
	else{
		print scalar localtime()."\n";
		print "bowtie2-align-s --phred33 -p 40 -x $out/spades/scaffolds -U $out/reads.fasta -S $out/TAX_annot/bowtie2_scaffolds.sam\n";
		`bowtie2-align --phred33 -p 40 -x $out/spades/scaffolds -U $out/reads.fasta -S $out/TAX_annot/bowtie2_scaffolds.sam`;
	}
	# Parse the alignments with samtools and bedtools
	print scalar localtime()."\n";
	print "samtools view -bS $out/TAX_annot/bowtie2_scaffolds.sam > $out/TAX_annot/bowtie2_scaffolds.bam\n";
	`samtools view -bS $out/TAX_annot/bowtie2_scaffolds.sam > $out/TAX_annot/bowtie2_scaffolds.bam`;
	print scalar localtime()."\n";
	print "samtools index $out/TAX_annot/bowtie2_scaffolds.bam\n";
	`samtools index $out/TAX_annot/bowtie2_scaffolds.bam`;
	print scalar localtime()."\n";
	print "samtools sort $out/TAX_annot/bowtie2_scaffolds.bam $out/TAX_annot/bowtie2_scaffolds.sorted\n";
	`samtools sort $out/TAX_annot/bowtie2_scaffolds.bam $out/TAX_annot/bowtie2_scaffolds.sorted`;
	print scalar localtime()."\n";
	print "samtools index $out/TAX_annot/bowtie2_scaffolds.sorted.bam\n";
	`samtools index $out/TAX_annot/bowtie2_scaffolds.sorted.bam`;
	print scalar localtime()."\n";
	print "bedtools intersect -abam $out/TAX_annot/bowtie2_scaffolds.sorted.bam -b $out/CDS_annot/gmhmm -bed > $out/TAX_annot/gmhmm.filt.bed\n";
	`bedtools intersect -abam $out/TAX_annot/bowtie2_scaffolds.sorted.bam -b $out/CDS_annot/gmhmm -bed > $out/TAX_annot/gmhmm.filt.bed`;
	`bedtools intersect -b $out/CDS_annot/gmhmm -a $out/TAX_annot/gmhmm.filt.bed -wb > $out/TAX_annot/semifinal_bed.txt `;
	`perl parse_gene.pl $out/TAX_annot/semifinal_bed.txt $out/TAX_annot/final.bed`;

	#Genomes from metagenome
	#create genomes directory
#	`mkdir genomes`;
#	print scalar localtime()."\n";
#	print "perl run_MaxBin.pl -contig scaffolds.fasta -reads $out/output_1.fastq -out $out/genomes/out maxbin\n";
#	`perl run_MaxBin.pl -contig scaffolds.fasta -reads $out/output_1.fastq -out $out/genomes/out maxbin`;
#olhar com Fausto passar da dourado para a tucunare	`checkm taxonomy_wf domain Bacteria data -x fasta $out/genomes/out maxbin`;

	#Creating the MYSQL database
	`mkdir $out/MYSQL`;
	`cp $out/TAX_annot/final.bed $out/MYSQL/gene.txt`;
	`cp $out/TAX_annot/FILT_MEGABLAST_NT-taxonomy.txt $out/MYSQL/blastn.txt`;
	`cut -f 1,4,5,6 KEGG_hier.2016.01.05 > $out/MYSQL/KEGG_hier.txt`;
	`cp  $out/CDS_annot/interprot.m8 $out/MYSQL/ko.txt`;
#	`cd $out/MYSQL`;
#	`cat createdb.sql | mysql -u laura -pcebioabcd`;
}
