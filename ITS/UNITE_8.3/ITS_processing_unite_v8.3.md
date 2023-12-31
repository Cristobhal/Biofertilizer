# Sequence processing

### Fungal sequences


```R
#A list of files to assembly is generated.
ls *.fastq.gz | perl -pe 's/_R.*.fastq.gz//g' | sort | uniq > lista
```


```R
#Use assembly sh script
The "assembly.sh" scrip generates the works to assembly files from list.
```


```R
#!/bin/bash
# Use: bash assemblyCASPER.sh NOMBRE_TRABAJO

SEQS=$(pwd)
SALIDAS=$(pwd)
BIN=/usr/local/bin
COUNT=0

for FAA in `cat lista`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr

echo "$BIN/casper <(cat $SEQS/$FAA"_R1.fastq") <(cat($SEQS/$FAA"_R2.fastq") -o $FAA.assembly.fastq -o $FAA"_assembly"" >>$*.$COUNT.scr

done
```


```R
#Quality filter, remove sequence with Quality score lower than 20. 
#This step is made to remove low quality assemblies

#!/bin/bash
SEQS=$(pwd)
SALIDAS=$(pwd)
for FAA in `cat lista`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr

echo fastq_quality_filter -q 20 -p 80 -i "$SEQS/$FAA".fastq -o "$SEQS/$FAA"_tr.fastq >>$*.$COUNT.scr

done
```


```R
#Transform filtered sequences to fasta
for x in `cat lista`; do sed -n '1~4s/^@/>/p;2~4p' "$x"_tr.fastq > "$x"_ensamble.fasta; done
```


```R
Sequences names are changed with “header.fasta.numbers.pl”.

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
## Extract representative OTUs, -i indicates input file, -f is the fasta file with the sequences to be extracted. 

pick_rep_set.py -i biofert.otu -f biofert.fas -o biofert.rep.fna
```


```R
#Assign taxonomy against most recent UNITE database. We use the complete database to accuretly assign non fungal taxa

parallel_assign_taxonomy_blast.py -i biofert.rep.fna -o taxonomy -r /unite_8.3/sh_qiime_release_all_10.05.2021/sh_refs_qiime_ver8_97_all_10.05.2021.fasta -t /unite_8.3/sh_qiime_release_all_10.05.2021/sh_taxonomy_qiime_ver8_97_all_10.05.2021.txt

#Create a list of fungal ITS sequences
cat taxonomy/biofert.rep_tax_assignments.txt | grep "k__Fungi" | cut -f1 >ids_screened.txt

#Create a list of sequences that are not fungi
cat taxonomy/biofert.rep_tax_assignments.txt | grep -v "k__Fungi" | cut -f1 >ids_REMOVE_biom.txt

#Extract ITS sequences from ITS list and make a new representative file
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ids_screened.txt biofert.rep.fna > biofert.screened.fna
```


```R
#Make biom object

make_otu_table.py -i biofert.otu -t taxonomy/biofert.rep_tax_assignments.txt -o biofert.biom 

#Remove singletons and non fungal sequences. Chimeras are removed later
filter_otus_from_otu_table.py -i biofert.biom -e ids_REMOVE_biom.txt -o biofert_screened.biom -n2 
```


```R
#Remove chimeras from representative sequences

#Identify chimera
parallel_identify_chimeric_seqs.py -m blast_fragments -i biofert.screened.fna -o biofert.chimera.txt -X biofertblast --id_to_taxonomy_fp /home/cristobal/DB/unite_8.3/sh_qiime_release_all_10.05.2021/sh_taxonomy_qiime_ver8_97_all_10.05.2021.txt -r /home/cristobal/DB/unite_8.3/sh_qiime_release_all_10.05.2021/sh_refs_qiime_ver8_97_all_10.05.2021.fasta


#Filter chimeric sequences from biom table
filter_otus_from_otu_table.py -i biofert_screened.biom -e biofert.chimera.txt -o biofert_chimera.biom
```


```R
#Create tables
biom convert --to-tsv -i biofert_chimera.biom -o biofert.biom.tsv --table-type "Taxon table" --header-key=taxonomy

#Remove OTUs labeled as "None"
grep -v 'None' biofert.biom.tsv > biofert_filt.biom.tsv; mv biofert_filt.biom.tsv biofert.biom.tsv
```


```R
#Split taxonomy table from OTU table.
perl -pe 's/\; /\;/g' biofert.biom.tsv | awk '{print $1,"\t",$NF}' | perl -pe 's/\;/\t/g' > biofert_ITS_tax.tsv

#Split OTU abundance for sample from OTU table: 
cut -f 1-3 biofert.biom.tsv > biofert_ITS_otu.tsv
```

### All eukaryotic sequences


```R
#All ITS sequences are used for further processing

#Remove no ITS sequences
cat taxonomy/biofert.rep_tax_assignments.txt | grep -v "No blast hit" | cut -f1 > all_ids_screened.txt

#Create a list of sequences that are not fungi
cat taxonomy/biofert.rep_tax_assignments.txt | grep "No blast hit" | cut -f1 > all_ids_REMOVE_biom.txt

#Extract ITS sequences from ITS list and make a new representative file
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' all_ids_screened.txt biofert.rep.fna > all_biofert.screened.fna
```


```R
#Remove singletons and non ITS sequences. Chimeras are removed later

filter_otus_from_otu_table.py -i biofert.biom -e all_ids_REMOVE_biom.txt -o all_biofert_screened.biom -n2 
```


```R
#Remove chimeras from representative sequences

#Identify chimera

parallel_identify_chimeric_seqs.py -m blast_fragments -i all_biofert.screened.fna -o all_biofert.chimera.txt -X biofertblast --id_to_taxonomy_fp /home/cristobal/DB/unite_8.3/sh_qiime_release_all_10.05.2021/sh_taxonomy_qiime_ver8_97_all_10.05.2021.txt -r /home/cristobal/DB/unite_8.3/sh_qiime_release_all_10.05.2021/sh_refs_qiime_ver8_97_all_10.05.2021.fasta


#Filter chimeric sequences from biom table

filter_otus_from_otu_table.py -i all_biofert_screened.biom -e all_biofert.chimera.txt -o all_biofert_chimera.biom
```


```R
#Create tables
biom convert --to-tsv -i all_biofert_chimera.biom -o all_biofert.biom.tsv --table-type "Taxon table" --header-key=taxonomy

#Remove OTUs labeled as "None"
grep -v 'None' all_biofert.biom.tsv > all_biofert_filt.biom.tsv; mv all_biofert_filt.biom.tsv all_biofert.biom.tsv
```


```R
#Split taxonomy table from OTU table.
perl -pe 's/\; /\;/g' all_biofert.biom.tsv | awk '{print $1,"\t",$NF}' | perl -pe 's/\;/\t/g' > all_biofert_ITS_tax.tsv

#Split OTU abundance for sample from OTU table: 
cut -f 1-3 all_biofert.biom.tsv > all_biofert_ITS_otu.tsv
```
