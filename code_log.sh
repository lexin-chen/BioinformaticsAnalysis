######### lab 3
ls *.gz
#finding files ending with gz

#Ctrl C cancel command -> fresh new prompt
#Ctrl R reverse search through command history
#Ctrl L clear command will clear screen
#b go backward/ g go to beginning/ G go to the end/ q to quit

cp GCF_000209225.1_ASM20922v1_genomic.gff GCF_000209225.1_ASM20922v1_genomic.copy.gff
# to make a copy of document

mv GCF_000209225.1_ASM20922v1_genomic.copy.gff backup
# to move backup file to backup 

chmod -w GCF_000209225.1_ASM20922v1_genomic.backup.gff
chmod -x GCF_000209225.1_ASM20922v1_genomic.backup.gff
# change permission (-) subtract permission (+) add permission
# w means write x means execute

rm -r backup
#delete everythng in the directory

grep -i NNNNNNNNNN GCF_000209225.1_ASM20922v1_genomic.fna | less
# finding a sequence in a file

 for geneid in 5509978 5519844 5500059 5511186 5515112 5499603 5504018 5497282 5497102 5502183  5517249 5512933 5497772 5514511 5515994 5495809 5521394 5518406 5519811 5500343 5518755 5513869

 do

echo $geneid

 grep GeneID:${geneid} /home/ubuntu/refseq/invertebrate/GCF_000209225.1/GCF_000209225.1_ASM20922v1_genomic.gff > ${geneid}.gff

 grep exon ${geneid}.gff > ${geneid}.mrna.gff

 grep CDS ${geneid}.gff > ${geneid}.cds.gff

 bedtools getfasta -s -fi /home/ubuntu/refseq/invertebrate/GCF_000209225.1/GCF_000209225.1_ASM20922v1_genomic.fna -fo ${geneid}.cds.fa -bed ${geneid}.cds.gff 

 bedtools getfasta -s -fi /home/ubuntu/refseq/invertebrate/GCF_000209225.1/GCF_000209225.1_ASM20922v1_genomic.fna -fo ${geneid}.mrna.fa -bed ${geneid}.mrna.gff 

  union ${geneid}.cds.fa -outseq ${geneid}.cds.union.fa

 transeq ${geneid}.cds.union.fa -outseq ${geneid}.aa.fa

 gt gff3 -sort -tidy ${geneid}.gff | gt stat -genelengthdistri -exonlengthdistri -intronlengthdistri -cdslengthdistri -addintrons -force -o ${geneid}.counts.txt 

 done

#using a bash script to figure out gene length distribution/ intron exon lengths/ finding count of these

bash analyze-gene-script.sh

cat counts.txt 
#to view specifics of that gene

######### lab 4
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

###### lab 5

ncbi-acc-download -F fasta -m protein XP_001630613
#download PROTEIN from genbank

blastp -db  ~/blastdbaa/allprotein.fas  -query ~/lab5/XP_001630613.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -out ~/lab5/XP_001630613.blastp.detail.out
#perform a blast search using this query protein. (remember to change cd location and the file title. last part is the output location 

blastp -db  ~/blastdbaa/allprotein.fas  -query ~/lab5/XP_001630613.fa -outfmt 0 | less -S
#you will be able to see the alignment of the high scoring pairs. similar to blast results

awk '{if ($3>=207 && $7 >50 && $6<0.0000000001)print $0 }' ~/lab5/XP_001630613.blastp.detail.out > ~/lab5/XP_001630613.blastp.detail.filtered.out
#207 means 207 amino acids. It is from the original amount/2. 50 is matches that alignemnets with a length that is equivalent to 50% of the length of the query protein. evalue is 1x10^-10

grep Nematostella ~/lab5/XP_001630613.blastp.detail.filtered.out | awk '!_[$1]++' | less -S
grep Exaiptasia ~/lab5/XP_001630613.blastp.detail.filtered.out | awk '!_[$1]++' | less -S
grep Homo ~/lab5/XP_001630613.blastp.detail.filtered.out | awk '!_[$1]++' | less -S
grep Orbicella ~/lab5/XP_001630613.blastp.detail.filtered.out | awk '!_[$1]++' | less -S
grep Strongylocentrotus ~/lab5/XP_001630613.blastp.detail.filtered.out | awk '!_[$1]++' | less -S
#use these grep tools to find the top hits per category

for proteinid in XP_001626384.1 XP_001623938.1 XP_001618979.1 XP_001640902.1 XP_028517720.1 XP_020913254.1 XP_020903999.1 XP_024307277.1 NP_787106.4 NP_001340387.1 NP_001073878.2 NP_060770.3 XP_020621543.1 XP_020608251.1 XP_020602351.1 XP_020602425.1 XP_011679460.1 XP_011665765.1  XP_011684022.1 XP_011675659.1 XP_011669479.1
do
echo $proteinid
ncbi-acc-download -F fasta -m protein ${proteinid}
sleep 1
done
#run this script and it will download all the homologs (when doing this do not redownload the sequence beforehand only the homologs not previously downloaded or else multiple .fas)

cat *.fa >allhomologs.fa
#will put all the sequences in one file. Be wary that if there is anything else with *.fa extension, this will also put into the file (so remember to ls -a before doing this)

muscle -in allhomologs.fa -out allhomologs.aligned.fa

t_coffee -other_pg seq_reformat -in  allhomologs.aligned.fa -output sim
#average percent identity between the query sequence and all other sequences is different from average percent identity among all sequences.former means look at your file and find the corresponding average. Latter means find the total average. 

alv -kli --majority allhomologs.aligned.fa | less -RS
#best viewer for alignments

sed -i '1 s/^.*$/>5509158/' NEMVEDRAFT_v1g189134.cds.union.fa
#rename the definition line with its gene ID


##### lab 6

iqtree -s ~/lab5/allhomologs.aligned.fa
#making a phylogenetic tree in IQ-TREE

less allhomologs.aligned.fa.iqtree
#look at what's up

gotree reroot midpoint -i allhomologs.aligned.fa.treefile -o allhomologs.aligned.fa.midpoint.treefile

#reroot and get midpoint


