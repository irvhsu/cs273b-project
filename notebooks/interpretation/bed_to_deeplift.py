from math import ceil
import numpy as np
from dragonn.models import SequenceDNN
import gzip
import sys

SEQ_LEN = 145
STRIDE  = 29

def fragment_seq(seq):
    center = ceil(len(seq) / 2)
    shifts = int((center - SEQ_LEN / 2) / STRIDE)
    seqs = []
    for i in range(-shifts, shifts + 1):
        start = int(center + i*STRIDE - SEQ_LEN / 2)
        end   = int(center + i*STRIDE + SEQ_LEN / 2 + 1)
        
        seqs += [seq[start:end].upper()]
    return seqs

bases = ['A', 'T', 'C', 'G']

def one_hot_encode_seq(seq):
    if "N" in seq:
        seq = seq.replace("N", "A")
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
    return result

######################################################################

import deeplift
from deeplift.conversion import keras_conversion as kc

model = SequenceDNN.load('../models/models/100n1_100n2_8w1_15w2.arch.json', '../models/models/100n1_100n2_8w1_15w2.weights.h5')

deeplift_model = kc.convert_sequential_model(
                    model.model, nonlinear_mxts_mode=deeplift.blobs.NonlinearMxtsMode.DeepLIFT)
deeplift_contribs_func = deeplift_model.get_target_contribs_func(find_scores_layer_idx=0, target_layer_idx=-1)

#######################################################################

def batch(coords, seqs, out, deepift_func):
    X = seqs_to_encoded_matrix(seqs)
    scores = np.array(deeplift_contribs_func(task_idx=0, input_data_list=[X], batch_size=128, progress_update=1000))
    
    # Map deeplift output back into genomic coordinates
    seq_idx = 0
    activity = {}
    for chrom, start, end, num in coords:
        center = ceil((end - start) / 2)
        shifts = int((center - SEQ_LEN / 2) / STRIDE)
        begin = int(center - shifts*STRIDE - SEQ_LEN / 2)
        for i in range(num):
            s = scores[i+seq_idx].reshape(4, SEQ_LEN).T
            for j, nt in enumerate(s):
                if abs(min(nt)) > max(nt):
                    act = min(nt)
                else:
                    act = max(nt)
                pos = start + begin + i * STRIDE + j
                
                if chrom not in activity:
                    activity[chrom] = {}
                if pos not in activity[chrom]:
                    activity[chrom][pos] = []
                activity[chrom][pos] += [act]
        seq_idx += num

    for chrom, positions in activity.items():
        for pos, score in positions.items():
            out.write('\t'.join([chrom, str(pos), str(np.mean(score))]) + '\n')

#######################################################################

with gzip.open("../../data/dnase/{}.dnase.fa.gz".format(sys.argv[1])) as dnase:
    chrom = None
    site_num = 0
    seqs, coords = [], []
    for line in dnase:
        site_num += 1
        if line[0] == '>':
            if chrom:
                frags = fragment_seq(seq)
                seqs += frags
                coords += [(chrom, start, end, len(frags))]
            count, rest = line[1:].strip().split('::')
            chrom, rest = rest.split(':')
            start, end = map(int, rest.split('-'))
            seq = ''
        else:
            seq += line.strip()

        if site_num % 20000 == 0:
            if site_num > 50000:
                with gzip.open("../../data/dnase/{}_deeplift_{}.tsv.gz".format(sys.argv[1], site_num), 'w') as out:
                    batch(coords, seqs, out, deeplift_contribs_func)
            seqs, coords = [], []

