This directory contains files on the sequences tested in the pilot experiments and the raw DNA and mRNA counts as summarized below.

Files:


The 145-bp sequences tested in the pilot experiments are specified in the file
*PilotDesign.sequences.txt 

The format of the files are two columns. The first column contains an ID for the sequence and the second the actual 145-bp DNA sequence.
IDs are of the form GROUP_REGIONIDINGROUP_TILEPOS where: 
- GROUP is from the set {high,range,controlhigh} where 
'high' denotes a sequence from regions selected based on being in candidate enhancer chromatin states in HepG2 with the strongest H3K27ac dip scores 
'range' denotes a sequence from regions selected for having a range of dip scores 
'controlhigh' denotes regions selected from being in a candidate enhancer with strong H3K27ac dip scores in K562, 
but were in low activity states in HepG2.
- REGIONIDINGROUP is an ID among selections in the GROUP in order of their selection. 
- TILEPOS is the tile position and ranges from 0 to 8 where TILEPOS=0 is the first tile and 8 the last tile. 
Tile positions are spaced in consecutive 30-bp intervals. 

=============================================================================================================================

The coordinates of the 145bp sequences tested in the pilot designs are in the files
*coords_PilotDesign_hg18.txt (original hg18 coordinates)
*coords_PilotDesign_hg19.txt (hg19 liftover coordinates)

The files contain four columns.
The first column contains an ID for the sequence and is defined as above. 
The second, third, and fourth columns specify the chromosome, start, and end coordinates of the region respectively.
Coordinate intervals are 0-based with first coordinate inclusive and second exclusive.


=============================================================================================================================
The DNA plasmid counts are in the following files:
*DNACOUNTS/PilotDesign_SV40P_Plasmid_Rep1.counts 
*DNACOUNTS/PilotDesign_SV40P_Plasmid_Rep2.counts

The file name indicates which replicate (Rep1 or Rep2) to which the counts correspond.
In the files the first line is a header liner with each column header corresponding to one of the 24-barcodes for the sequence tested.
The remaining lines first have an entry containing the ID of the sequence as specified in the above files
and the remaining columns contain a count for each barcode with which the sequence was tested.


=============================================================================================================================

The mRNA counts are in the following files:
*HEPG2/HepG2_PilotDesign_SV40P_mRNA_Rep1.counts
*HEPG2/HepG2_PilotDesign_SV40P_mRNA_Rep2.counts
*K562/K562_PilotDesign_SV40P_mRNA_Rep1.counts
*K562/K562_PilotDesign_SV40P_mRNA_Rep2.counts

The directory and file name indicates if the experiments were conducted in HepG2 or K562 cells.
The file name also indicates which of the two replicates (Rep1 or Rep2) the experiments correspond.
The format of the files are the same for the DNA count files above.