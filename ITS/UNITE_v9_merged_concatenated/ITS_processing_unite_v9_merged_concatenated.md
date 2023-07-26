# Sequence processing

```R
##Concatenate merged and unmerged sequences
#!/bin/bash
#LDA 2023
#Usage: ./merge_reads.sh FILE1 FILE2 OUTPUT
#files are fastq or fastq.gz files
# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <forward_file> <reverse_file> <output_prefix>"
    exit 1
fi

# Input file paths
forward_file="$1"
reverse_file="$2"

# Output prefix
output_prefix="$3"

# Output file paths
concat_fastq="${output_prefix}_concat.fastq"
merged_fastq="${output_prefix}_merged.fastq"
concat_fasta="${output_prefix}_concat.fas"
merged_fasta="${output_prefix}_merged.fas"
concat_and_merged_fasta="${output_prefix}_concat_and_merged.fas"

# Run pandaseq to merge the reads
pandaseq -f "$forward_file" -r "$reverse_file" -O 10 -w "$concat_fastq" -U "$merged_fastq"

# Convert concatenated fastq to fasta format
seqtk seq -A "$concat_fastq" > "$concat_fasta"

# Convert merged fastq to fasta format
seqtk seq -A "$merged_fastq" > "$merged_fasta"

echo "Total assembled sequences:"
grep -c ">" $concat_fasta

echo "Total merged -- sequences:"
grep -c ">" $merged_fasta

# Concatenate the fasta files and compress using pbzip
cat "$concat_fasta" "$merged_fasta" | pbzip2 -p2 > "$concat_and_merged_fasta.bz2"

# Remove temporary files
rm "$concat_fastq" "$merged_fastq" "$concat_fasta" "$merged_fasta"

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
perl header.fasta.numbers.pl rhizosphere rhizoLDA_concat_and_merged.fas
perl header.fasta.numbers.pl substrate substrateLDA_concat_and_merged.fas
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

parallel_assign_taxonomy_blast.py -i biofert.rep.fna -o taxonomy -r /home/cristobal/DB/UNITE_9_all/sh_refs_qiime_ver9_97_all_29.11.2022.fasta -t /home/cristobal/DB/UNITE_9_all/sh_taxonomy_qiime_ver9_97_all_29.11.2022.txt

#Create a list of fungal ITS sequences, this step removes sequences without hits
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

parallel_identify_chimeric_seqs.py -m blast_fragments -i biofert.screened.fna -o biofert.chimera.txt -X biofertblast --id_to_taxonomy_fp /home/cristobal/DB/UNITE_9_all/sh_taxonomy_qiime_ver9_97_all_29.11.2022.txt -r /home/cristobal/DB/UNITE_9_all/sh_refs_qiime_ver9_97_all_29.11.2022.fasta

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
