#make sure to use large instance when doing the domain stuff.
 ~/my_interproscan/interproscan-5.36-75.0/interproscan.sh -appl Pfam -i ~/lab5/allhomologs.fa -f TSV,GFF3,HTML,SVG  -goterms -iprlookup
#Run interproscan

grep XP_001630613  allhomologs.fa.tsv |less -S
#pull out the anootations for one protein

cut -f 1 allhomologs.fa.tsv  | sort | uniq -c
to compare the number of annotations among proteins
--
awk 'BEGIN{FS="\t"} {print $1"\t"$3"\t"$7"@"$8"@"$5}' allhomologs.fa.tsv | datamash -sW --group=1,2 collapse 3 | sed 's/,/\t/g' |  sed 's/@/,/g' > allhomologs.fa.rear.tsv
