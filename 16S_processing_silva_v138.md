## Sequence processing


```R
#A list of files to assembly is generated.
ls *.fastq.gz | perl -pe 's/_R.*.fastq.gz//g' | sort | uniq >lista
```


```R
#Use assembly sh script
The "assembly.sh" scrip generates the works to assembly files from list.
```


```R
# assembly.sh

#!/bin/bash

SEQS=$(pwd)
SALIDAS=$(pwd)
BIN=/usr/bin
BIN2=/usr/bin
COUNT=0
for FAA in `cat lista`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo "zcat $SEQS/$FAA"_R1.fastq.gz" | $BIN2/fastx_trimmer -l 250 > $SEQS/$FAA"tr_R1.fastq" &" >>$*.$COUNT.scr mo
echo "zcat $SEQS/$FAA"_R2.fastq.gz" | $BIN2/fastx_trimmer -l 250 > $SEQS/$FAA"tr_R2.fastq"" >>$*.$COUNT.scr
echo "$BIN/pandaseq -B -f $SEQS/$FAA"tr_R1.fastq" -r $SEQS/$FAA"tr_R2.fastq" -t 0.95 -l 250 -L 470 -o 15 -w $SALIDAS/$FAA"_$*.fasta" -G $SALIDAS/$FAA"-$*.log.bz2"" >>$*.$COUNT.scr

done
```


```R
Preceding script utilices “fastx_trimmer” that crop sequence to a size that is defined with -l. 
The obtained file is assembled with “pandaseq”, -t can take values for 0 to 1 and alignments with lower values 
are discarded; -l is the minimum sequence length -L maximum sequence length; -o minimum overlap between assembled 
sequences.
```


```R
#Run script to generate works
bash assembly.sh ensamble
```


```R
#Running works in the cluster:
for N in `ls *.scr`; do qsub $N; done
```


```R
#Sequences names are changed with “header.fasta.numbers.pl”.

# Luis David Alcaraz 2013-04-11

my $prefix = $ARGV[0]; chomp $prefix;
my $f =  1;

my $fasta_file = $ARGV [1]; chomp $fasta_file;

my $fh;
open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";
open(OUT, ">$fasta_file.numbered.fas") || die "can't open $fasta_file.numbered.f
as\n";

my %sequence_data;
while (read_fasta_sequence($fh, \%sequence_data)) {
   print OUT ">$sequence_data{header}\n$sequence_data{seq}\n";
}

close $fh;
close OUT;

sub read_fasta_sequence {
   my ($fh, $seq_info) = @_;

$seq_info->{seq} = undef; # clear out previous sequence

   # put the header into place
   $seq_info->{header} = $seq_info->{next_header} if $seq_info->{next_header};

   my $file_not_empty = 0;
   while (<$fh>) {
      $file_not_empty = 1;
      next if /^\s*$/;  # skip blank lines
      chomp;    

      if (/^>/) { # fasta header line
         my $h = $_;    
         $h =~ s/>/$prefix\_$f\ /;
     $f++;
         if ($seq_info->{header}) {
            $seq_info->{next_header} = $h;
            return $seq_info;
       
         }              
         else { # first time through only
            $seq_info->{header} = $h;
         }              
      }         
      else {    
         s/\s+//;  # remove any white space
         s/\n\n/\n/;
         $seq_info->{seq} .= $_;
      }         
   }    

   if ($file_not_empty) {
      return $seq_info;
   }    
   else {
      # clean everything up
      $seq_info->{header} = $seq_info->{seq} = $seq_info->{next_header} = undef;

      return;   
   }    
}
```


```R
#Run the script for each sample
perl header.fasta.numbers.pl substrate substrate_ensamble.fasta
perl header.fasta.numbers.pl rhizosphere rhizosphere_ensamble.fasta
```


```R
#Concatenate all samples
cat *.numbered.fas > biofert.fas
```


```R
#Edit sequence name to leave only the first part, which refers to sample name. 
perl -i.bak -pe "s/\ .*//g" biofert.fas
```


```R
#Sequence clustering of OTUs at 97%  of identity was done with cd-hit-est, -c indicates the identity for clustering.

cd-hit-est -i biofert.fas -c 0.97 -o biofert.fasout -T 20 -M 0
```


```R
#The clustering file is converted in a file that can be read by QIIME.
perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g;s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_ \n"; exit}' biofert.fasout.clstr > biofert.otu
sed -i '1d' biofert.otu
```


```R
#Extract representative OTUs, -i indicates input file, -f is the fasta file with the sequences to be extracted. 

pick_rep_set.py -i biofert.otu -f biofert.fas -o biofert.rep.fna
```


```R
#Filter sequences that are not 16S by blasting against a smaller DB (OTUs 70% id).


parallel_assign_taxonomy_blast.py -i biofert.rep.fna -o no16S_screen -r /qiime/gg_otus-13_8-release/rep_set/70_otus.fasta -t /qiime/gg_otus-13_8-release/taxonomy/70_otu_taxonomy.txt

cat no16S_screen/biofert.rep_tax_assignments.txt | grep -c "No blast hit"
cat no16S_screen/biofert.rep_tax_assignments.txt | grep -v "No blast hit" | cut -f1 >ids_screened.txt
cat no16S_screen/biofert.rep_tax_assignments.txt | grep "No blast hit" | cut -f1 >ids_REMOVE_biom.txt

#Extract 16S sequences and make a new representative file
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ids_screened.txt biofert.rep.fna > biofert.screened.fna
```


```R
#Taxonomic assignment was done with QIIME using Silva database as reference.

parallel_assign_taxonomy_blast.py -i biofert.screened.fna -o taxonomy -r SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta -t SILVA_128_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_7_levels.txt

cat taxonomy/biofert.screened_tax_assignments.txt | grep "No_blast_hit" | cut -f1 >>ids_REMOVE_biom.txt
cat taxonomy/biofert.screened_tax_assignments.txt | grep -i "mitochondria" | cut -f1 >>ids_REMOVE_biom.txt
cat taxonomy/biofert.screened_tax_assignments.txt | grep -i "chloroplast" | cut -f1 >>ids_REMOVE_biom.txt
```


```R
#Update taxonomy with silva v138.1
$R
library(dada2)
library(Biostrings)
biofert <- readDNAStringSet("biofert.screened.fna") 
seqs <- getSequences(biofert)  

taxa <- assignTaxonomy(seqs, silva_nr99_v138_train_set.fa.gz, multithread = 4); save.image()
asv_seqs <- colnames(biofert)
asv_headers <- names(biofert)


asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "Silva_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
savehistory()
q()

awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' Silva_taxonomy.tsv > Silva_taxonomy_clean.tsv
sed -i '1d' Silva_taxonomy_clean.tsv
```


```R
#Make otu table
sed -i 's/ /_/g' taxonomy/biofert.screened_tax_assignments.txt

make_otu_table.py -i biofert.otu -t taxonomy/biofert.screened_tax_assignments.txt -o biofert.biom 

#Remove singletons and non 16S sequences
filter_otus_from_otu_table.py -i biofert.biom -e ids_REMOVE_biom.txt -o biofert_screened.biom -n2 ; mv biofert_screened.biom biofert.biom
```


```R
#Remove chimeras

#Align sequences
parallel_align_seqs_pynast.py -i biofert.screened.fna -o chimalign -X biofert

#Identify chimeric sequences
parallel_identify_chimeric_seqs.py -m blast_fragments -i biofert.screened.fna -a chimalign/biofert.screened_aligned.fasta -o biofert.chimera.txt -X biofertblast --id_to_taxonomy_fp 97_otu_taxonomy.txt -r 97_otus.fasta

#Filter OTUs
filter_otus_from_otu_table.py -i biofert.biom -e biofert.chimera.txt -o biofert_chimera.biom; mv biofert_chimera.biom biofert.biom

#Make list of non-chimeric sequences
awk '{print $1}' biofert.chimera.txt > chim_list

cut -f1  Silva_taxonomy_clean.tsv > all_otus

cat all_otus chim_list | sort | uniq -c | grep -w '1' | awk '{print$2}' > no_chim_list

```


```R
#Transform biom table to tabular format
biom convert --to-tsv -i biofert.biom -o biofert.biom.tsv --header-key=taxonomy

```


```R
#Split OTU abundance for sample from OTU table: 
cut -f 1-3 biofert.biom.tsv > biofert_16S_otu.tsv
sed -i -e "1d" biofert_16S_otu.tsv
```
