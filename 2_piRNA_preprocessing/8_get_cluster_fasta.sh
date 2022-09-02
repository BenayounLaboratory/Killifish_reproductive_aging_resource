#get_cluster_fasta.sh

#Obtain fasta sequences of clusters

#this script uses bedtools v2.27.1

bedtools getfasta -fi Masked_Genome.fa -bed merged_cluster_coordinates.CLEAN.bed -name > clusters.fasta