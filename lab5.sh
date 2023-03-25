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

grep Nematostella ~/lab5/XP_001630613.blastp.detail.filtered.out | awk '!_[$1]++' | cat 
grep Exaiptasia ~/lab5/XP_001630613.blastp.detail.filtered.out | awk '!_[$1]++' | cat 
grep Homo ~/lab5/XP_001630613.blastp.detail.filtered.out | awk '!_[$1]++' | cat
grep Orbicella ~/lab5/XP_001630613.blastp.detail.filtered.out | awk '!_[$1]++' | cat
grep Strongylocentrotus ~/lab5/XP_001630613.blastp.detail.filtered.out | awk '!_[$1]++' | cat

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
#do not download original one again

muscle -in allhomologs.fa -out allhomologs.aligned.fa

t_coffee -other_pg seq_reformat -in  allhomologs.aligned.fa -output sim
#average percent identity between the query sequence and all other sequences is different from average percent identity among all sequences.former means look at your file and find the corresponding average. Latter means find the total average. 

alv -k allhomologs.aligned.fa | less -R
alv -kl allhomologs.aligned.fa | less -RS
alv -kli --majority allhomologs.aligned.fa | less -RS
#last one is best viewer for alignments because columns are only colored if majority of the residues is in >=50% of the sequence

sed -i '1 s/^.*$/>5509158/' NEMVEDRAFT_v1g189134.cds.union.fa
#rename the definition line with its gene ID

t_coffee -other_pg seq_reformat -in  allhomologs.aligned.fa -action +rm_gap 50 -out   allhomologs.aligned.r50.fa 
#improve any column that contains greater than 50% gapped residued 

alv -kli --majority allhomologs.aligned.r50.fa | less -RS
#you will be able to calculate number of columns removed from alignment width.


