cat AF254382.1.fa KP761313.1.fa > 18s.fa
#putting these two files into a single file

muscle -in 18s.fa -out 18s.aligned.fa
#aligning the two sequences

alan 18s.aligned.fa
#view the aligned to see conserved/ nonconserved regions

t_coffee -other_pg seq_reformat -in 18s.aligned.fa -output sim

can calculate the percent identify between two aligned sequences

#comment out the genome you don't want.

genome="/home/ubuntu/refseq/invertebrate/GCF_001417965.1/GCF_001417965.1_Aiptasia_genome_1.1_genomic"

#genome="/home/ubuntu/refseq/invertebrate/GCF_000209225.1/GCF_000209225.1_ASM20922v1_genomic"

 for geneid in 110253432

 do

 for isoform in X2

 do

echo $geneid

echo $isoform

 grep GeneID:${geneid} ${genome}.gff > ${geneid}.gff

 grep ${isoform} ${geneid}.gff > ${geneid}.${isoform}.gff

 grep CDS ${geneid}.${isoform}.gff > ${geneid}.${isoform}.cds.gff

#sometimes, there are multiple RNAs associated with an isoform. Let's just look at the first one of these

#by sorting and then selecting only unique star coordinates for the exons in the CDS

 cat ${geneid}.${isoform}.cds.gff | cut -f1-8 | sort | uniq > ${geneid}.${isoform}.cds.uniq.gff

 bedtools getfasta -s -fi ${genome}.fna -fo ${geneid}.${isoform}.cds.fa -bed ${geneid}.${isoform}.cds.uniq.gff 

  union ${geneid}.${isoform}.cds.fa -outseq ${geneid}.${isoform}.cds.union.fa

#let's change the name in the fasta file to the gene ID

sed -i "1 s/^.*$/>${geneid}/" ${geneid}.${isoform}.cds.union.fa

 

 transeq ${geneid}.${isoform}.cds.union.fa -outseq ${geneid}.${isoform}.aa.fa

#note: currently this outputs info for all isoforms and transcripts of this gene. 

 gt gff3 -sort -tidy ${geneid}.gff | gt stat -genelengthdistri -exonlengthdistri -intronlengthdistri -cdslengthdistri -addintrons -force -o ${geneid}.counts.txt 

 done

 done

#finding orthologs and their information

cat NEMVEDRAFT_v1g189134.cds.union.fa  110253432.X2.cds.union.fa > NucRec.cds.fa
#putting all the orthologs into one file.  

muscle -in NucRec.cds.fa -out NucRec.cds.aligned.fa
#aligning CDS in muscle

alan NucRec.cds.aligned.fa
#view the alignment in alan

t_coffee -other_pg seq_reformat -in NucRec.cds.aligned.fa -output sim
#calculation for percent identity

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < NucRec.cds.aligned.fa | while read line; do echo $line | grep -v '>' | grep -o "-" | sort | uniq -c; echo $line | grep '>' ; done 
