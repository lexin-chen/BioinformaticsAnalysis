#running bootstrap

iqtree -s allhomologs.aligned.r50.fa -b 100 
#will take literally forever
#screen
#ctrl A+D: will run in the background and you can do other stuff
#screen -r running thing will come back

gotree reroot midpoint -i allhomologs.aligned.r50.fa.treefile -o allhomologs.aligned.r50.fa.midpoint.treefile
#midpoint resulting tree to that file

#need to make speciestree.tre before this step and make sure the names of the species match up
java -jar ~/Notung-2.9.1.2/Notung-2.9.1.2.jar -g ~/lab6/allhomologs.aligned.r50.fa.midpoint.treefile -s speciestree.tre --speciestag prefix --reconcile --savepng --treestats --events --homologtabletabs

#events file:see duplication+losses 
#reconciled get png of tree/# of losses and duplications displayed
sftp -i C:\Users\aestr\Downloads\lexin.pem ubuntu@3.229.208.174

get allhomologs.aligned.r50.fa.midpoint.treefile.reconciled.xml C:\Users\aestr\Downloads
put (file) to ___
#last two are instruction to transfer files from terminal to local device

python /home/ubuntu/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g allhomologs.aligned.r50.fa.midpoint.treefile.reconciled --include.species
#looking at the species tree on RecPhyloXML
#http://phylariane.univ-lyon1.fr/recphyloxml/recphylovisu
#https://cloudconvert.com/svg-to-pdf

#8)According to this tree, what is the relationship of XP_020913254.1 to XP_001640902.1 with respect to the earliest divergence among Cnidarians (n93)? Because it is with respect to earliest divergence or speciation from before. 

#to find out which node has the highest number of losses. 
use allhomologs.aligned.r50.fa.midpoint.treefile.reconciled.events.txt
