This directory contains files with normalized data from the pilot. 

Files are:
HEPG2/tablenorm_recenterends_HepG2_Rep1_20.txt
HEPG2/tablenorm_recenterends_HepG2_Rep2_20.txt
K562/tablenorm_recenterends_K562_Rep1_20.txt
K562/tablenorm_recenterends_K562_Rep2_20.txt

The directory and file name indicates if the experiments were conducted in HepG2 or K562 cells.
The file name also indicates which of the two replicates (Rep1 or Rep2) the experiments correspond.

Each row corresponds to one of the regions tested. 
First column is an ID for the region.
IDs are of the form GROUP_REGIONIDINGROUP where: 
- GROUP is from the set {high,range,controlhigh} where 
'high' denotes a sequence from regions selected based on being in candidate enhancer chromatin states in HepG2 with the strongest H3K27ac dip scores 
'range' denotes a sequence from regions selected for having a range of dip scores 
'controlhigh' denotes regions selected from being in a candidate enhancer with strong H3K27ac dip scores in K562, 
but were in low activity states in HepG2.
- REGIONIDINGROUP is an ID among selections in the GROUP in order of their selection. 

The remaining columns correspond to the normalized values for each
of the nine tile positions in order within the region.

