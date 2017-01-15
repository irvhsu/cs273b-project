# This script accepts a bed file and creates a bed file with sequence
# level and/or base level activity scores / importance scores

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

IN=$1
FASTA=$2 # must be absolute path
MODEL=$3 # must be absolute path
OUT=$3

############################################
mkdir $OUT
cp $IN $OUT/in.bed
cd $OUT

############################################

python $HERE/fix_bed_lengths.py in.bed corrected_lengths.bed

############################################

bedtools getfasta -fi $FASTA -bed corrected_lengths.bed > in.fa

############################################

python score_regions.py in.fa > sequence_level.bed

############################################

python score_nucleotides.py in.fa > nucleotide_level.bed
