"""
Given a file formatted as

chrom, start, end, sequence

create a file with deeplift scores for
"""
from math import ceil

SEQ_LEN = 145
STRIDE  = 29

def fragment_seq(seq):
    center = ceil(len(seq) / 2)
    shifts = int((center - SEQ_LEN / 2) / STRIDE)
    seqs = []
    for i in range(-shifts, shifts + 1):
        start = int(center - i*STRIDE - SEQ_LEN / 2)
        end   = int(center - i*STRIDE + SEQ_LEN / 2 + 1)
        
        seqs += [seq[start:end].upper()]
    return seqs


bases = ['A', 'T', 'C', 'G']

def one_hot_encode_seq(seq):
    result = np.zeros((len(bases), len(seq)))
    for i, base in enumerate(seq):
        result[bases.index(base), i] = 1
    return result
    
def seqs_to_encoded_matrix(seqs):
    # Wrangle the data into a shape that Dragonn wants.
    result = np.concatenate(
        map(one_hot_encode_seq, seqs)
    ).reshape(
        len(seqs), 1, len(bases), len(seqs[0])
    )
    
    # Check we actually did the encoding right.
    for i in range(len(seqs)):
        for j in range(len(seqs[0])):
            assert sum(result[i, 0, :, j]) == 1
            
    return result

seqs, coords = [], []

with open(sys.argv[1]) as dnase:
    for line in dnase:
        chrom, start, end, seq = line.strip().split()
        start, end = int(start), int(end)

        seqs += fragment_seq(seq)
        coords += [(chrom, start, end, len(seqs))]

X = seqs_to_encoded_matrix(seqs)
deeplift = model.deeplift(X)
np.save(sys.argv[1] + '.deeplift.npy', deeplift)

# Map deeplift output back into genomic coordinates
seq_idx = 0
for chrom, start, end, num in coords:
    for i in range(num):
        scores = deeplift[i]

    seq_idx += num
        
