#!/bin/bash

module load python
source activate rnaseq_basic

fa=$1
piRNA_ta_index=$2
piRNA_bt_index=$3
genome_bt_index=$4
tmpdir=$5

[ ! -d $tmpdir ] && mkdir -p $tmpdir

cd $tmpdir

N_reads=$(cat $fa | grep ">" | awk -F'\t' '{OFS="\t"; split($1,a,":"); print a[2]}' | addCols stdin)
echo -e "Number of starting reads = $N_reads"

uniq_fasta_to_uniq_fastq.py $fa > tmp.fq 

# for some reason tailor_sam_to_bed wont report tails longer than like 9 nt so use samtools sam2bed for this, grep GTGGAATTTGAGGAAAC for an example using both commands
#tailor_v11 map -i tmp.fq -p $index -n 4 | tailor_sam_to_bed | grep "TAGTGAGAGTGTTTTCACTTGAGTGCTGTG:957"

# find "tails" longer than 10 starting with G/A
#tailor_v11 map -i tmp.fq -p $piRNA_ta_index -n 4 |\
#    sam2bed - |\
#    awk '($2<=2 && $6=="+")' |\
#    awk -F'\t' '{
#        OFS="\t"
#        
#        split($15,a,":")
#        
#        split($4,b,":")
#        
#        if ( (length(a[3])>=17 && length(a[3])<=19) && (substr(a[3],1,1)=="G" || substr(a[3],1,1)=="A" )) {
#            print $1,$2,$3,b[1],b[2],$6, substr(b[1], 1, length(b[1])-length(a[3])), a[3]
#            }
#        }' \
#    > candidates.bed

tailor_v11 map -i tmp.fq -p $piRNA_ta_index -n 4 |\
    sam2bed - |\
    awk '($2<=2 && $6=="+")' |\
    awk -F'\t' '{
        OFS="\t"
        
        split($15,a,":")
        
        split($4,b,":")
        
        if ( length(a[3])>=10 ) {
            print $1,$2,$3,b[1],b[2],$6, substr(b[1], 1, length(b[1])-length(a[3])), a[3]
            }
        }' \
    > candidates.bed

N_reads_tail_gt10=$(cat candidates.bed | cut -f5 | sort | uniq | addCols stdin )
echo -e "Number of reads with tails greater than 10nt long = $N_reads_tail_gt10"

# find which "tails" map to the reference
cat candidates.bed | awk -F'\t' '{ print ">" $4":"$8":"$5 "\n" $8 }' > anti.fa

# align anti to genome, keep reads that only 1X
bowtie \
    -x $genome_bt_index \
    -f anti.fa \
    -p 4 \
    -m 1 \
    -v 0 \
    -S anti.genome.sam

cat anti.genome.sam | sam2bed - | cut -f4 | awk -F'\t' '{ split($1,a,":"); print ">"$1"\n"a[2] }' > anti.uniMapper.fa

N_uniMapper_genome=$(cat anti.uniMapper.fa | grep ">" | awk -F'\t' '{OFS="\t"; split($1,a,":"); print a[3]}' | addCols stdin)
echo -e "Number of reads uniquely map to genome = $N_uniMapper_genome"

# align with bowtie to the reference 
bowtie \
    -x $piRNA_bt_index \
    -f anti.uniMapper.fa \
    -p 4 \
    -a \
    --best \
    --strata \
    -m 1 \
    -v 0 \
    -S anti.sam
    
sam2bed < anti.sam | awk '($6=="-")' | cut -f 1-6 > anti.bed 

N_reads_map_antipiRNA=$(cat anti.bed | cut -f4 | awk -F'\t' '{ split($1,a,":"); print a[3] }' | addCols stdin)

echo -e "Number of reads mapping antisense to a piRNA = $N_reads_map_antipiRNA"
