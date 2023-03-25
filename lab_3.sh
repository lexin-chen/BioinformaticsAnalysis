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
