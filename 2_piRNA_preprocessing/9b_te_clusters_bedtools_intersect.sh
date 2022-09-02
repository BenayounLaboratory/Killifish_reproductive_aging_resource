#te_clusters_bedtools_intersect.sh

#this script uses bedtools v2.27.1

bedtools intersect -wo -a merged_cluster_coordinates.CLEAN.bed -b repeatmasker_bedfile.bed > intersect.txt