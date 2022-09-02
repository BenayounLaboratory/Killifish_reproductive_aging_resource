#make_TE_only_bowtie_index.sh

# this script generates a TE-only bowtie index that will be used
# for te locus level ping-pong analysis

# this script uses bowtie v1.2.3

bowtie-build Nothobranchius_furzeri_TEs.fa ./killi_TE_bowtie_ref