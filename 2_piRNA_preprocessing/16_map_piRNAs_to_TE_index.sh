#map_piRNAs_to_TE_index.sh

# this script uses bowtie v1.2.3
# this script uses samtools v1.10

for f in $(find "." -name '*_piRNA.fq.gz')
do
        f2=$(basename "${f}");
        of=$(basename "${f}" | sed 's/_piRNA\.fq\.gz/_mapped/g');
        of_bam="${of}.bam"
        of_sorted="${of}_sorted"
        of_sorted_bam="${of_sorted}.bam"
        of_sorted_bambai="${of_sorted}.bam.bai"
        of_sorted_bai="${of_sorted}.bai"

        bowtie -v 3 -a --best --strata -S killi_TE_bowtie_ref $f2 $of
        samtools view -bS -o $of_bam $of
        rm $of
        samtools sort $of_bam -o $of_sorted_bam
        rm $of_bam
        samtools index $of_sorted_bam
        mv $of_sorted_bambai $of_sorted_bai

done