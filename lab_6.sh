iqtree -s ~/lab5/allhomologs.aligned.fa
#making a phylogenetic tree in IQ-TREE

less allhomologs.aligned.fa.iqtree
#look at what's up

gotree reroot midpoint -i allhomologs.aligned.fa.treefile -o allhomologs.aligned.fa.midpoint.treefile

#reroot and get midpoint


