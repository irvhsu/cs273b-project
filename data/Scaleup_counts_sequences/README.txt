This directory contains files on the sequences tested in the scale-up experiments and the raw DNA and mRNA counts as summarized below.

Files:

The 145-bp sequences tested in the scale-up experiments are specified in the files
*ScaleUpDesign1.sequences.txt
*ScaleupDesign2.sequences.txt  

Each file corresponds to one of the two array designs.
The format of the files are two columns. The first column contains an ID for the sequence and the second the actual 145-bp DNA sequence.
IDs are of form CELLTYPE_STATE_REGIONIDINSTATE_TILEPOS_CHR_CENTER where: 
-CELLTYPE and STATE are the cell type (H1hesc, Hepg2, K562, or Huvec) and chromatin state ID (1..25) respectively 
based on which the corresponding regulatory region of the sequence was selected
-REGIONIDINSTATE is an integer ID for the region that the sequence corresponds to which is 
unique among region selections from the same CELLTYPE and STATE combination 
-TILEPOS is the tile position of the sequence within its region and ranges from 0 to 30 where TILEPOS=0 is the first tile and 30 the last tile. 
Tile positions are spaced in consecutive 5-bp intervals. 
-CHR and CENTER are the chromosome and coordinate of the center base (hg19) respectively of the corresponding regulatory region being tiled. 


=============================================================================================================================================


The coordinates of the 295-bp regions tiled in the scale-up are specified in the files
*coords_ScaleUpDesign1_hg19.txt
*coords_ScaleUpDesign2_hg19.txt 

The format of the files are four columns. The first column contains an ID for the region. The second, third, and fourth columns
specify the chromosome, start, and end coordinates of the region respectively.
Coordinate intervals are 0-based with first coordinate inclusive and second exclusive.
IDs are of the form CELLTYPE_STATE_REGIONIDINSTATE_CHR_CENTER where CELLTYPE, STATE, REGIONIDINSTATE, CHR, CENTER are defined
as above. There are no TILEPOS in these IDs as there is only one entry for each region.


=============================================================================================================================================

The DNA plasmid counts are in the following files:
*DNACOUNTS/ScaleUpDesign1_SV40P_Plasmid.counts
*DNACOUNTS/ScaleUpDesign1_minP_Plasmid.counts
*DNACOUNTS/ScaleUpDesign2_SV40P_Plasmid.counts
*DNACOUNTS/ScaleUpDesign2_minP_Plasmid.counts

The file name indicates which array design and for promoter type, minimal (minP) or SV40 (SV40P), experiments
the counts correspond. 
The first line is a header line indicating only one barcode per sequence.
The remaining lines have two columns with the first column containing the ID of the sequence as in the *.sequences.txt file
and the second column the corresponding count.


=============================================================================================================================================


The mRNA counts are in the following files:
*HEPG2/HepG2_ScaleUpDesign1_SV40P_mRNA_Rep1.counts
*HEPG2/HepG2_ScaleUpDesign1_SV40P_mRNA_Rep2.counts
*HEPG2/HepG2_ScaleUpDesign1_minP_mRNA_Rep1.counts
*HEPG2/HepG2_ScaleUpDesign1_minP_mRNA_Rep2.counts
*HEPG2/HepG2_ScaleUpDesign2_SV40P_mRNA_Rep1.counts
*HEPG2/HepG2_ScaleUpDesign2_SV40P_mRNA_Rep2.counts
*HEPG2/HepG2_ScaleUpDesign2_minP_mRNA_Rep1.counts
*HEPG2/HepG2_ScaleUpDesign2_minP_mRNA_Rep2.counts
*K562/K562_ScaleUpDesign1_SV40P_mRNA_Rep1.counts
*K562/K562_ScaleUpDesign1_SV40P_mRNA_Rep2.counts
*K562/K562_ScaleUpDesign1_minP_mRNA_Rep1.counts
*K562/K562_ScaleUpDesign1_minP_mRNA_Rep2.counts
*K562/K562_ScaleUpDesign2_SV40P_mRNA_Rep1.counts
*K562/K562_ScaleUpDesign2_SV40P_mRNA_Rep2.counts
*K562/K562_ScaleUpDesign2_minP_mRNA_Rep1.counts
*K562/K562_ScaleUpDesign2_minP_mRNA_Rep2.counts


The directory and file name indicates if the experiments were conducted in HepG2 or K562 cells.
The file names also specify for which array design (ScaleUpDesign1 or ScaleUpDesign2), 
whether a minimal promoter (minP) or SV40 promoter (SV40P) 
was used, and which of the two replicates for the combination (Rep1 or Rep2).
The format of the files are the same for the DNA count files above.