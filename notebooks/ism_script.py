from collections import *
import math
import numpy as np
from subprocess import Popen, PIPE

letterindex = {'A': 0, 'a': 0, 'T': 1, 't': 1, 'C': 2, 'c': 2, 'G': 3, 'g': 3, 'N': -1, 'n': -1}

def bases(chrom, start, end):
    seq_count = int(math.ceil((float(end - start) / 60.0)))
    sum_seq = ""
    for i in xrange(seq_count - 1):
        p = Popen(['samtools', 'faidx', '../../Genomes/hg19.fa', 'chr' + str(chrom) + ':' + str(start + 1 + i * 60) + '-' + str(start + 1 + (i + 1) * 60)], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate()
        sum_seq = sum_seq + output.split('\n')[1]
    p = Popen(['samtools', 'faidx', '../../Genomes/hg19.fa', 'chr' + str(chrom) + ':' + str(start + 1 + (seq_count - 1) * 60) + '-' + str(end)], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output, err = p.communicate()
    sum_seq = sum_seq + output.split('\n')[1]
    return sum_seq

import numpy as np
from collections import defaultdict
import json
from dragonn import models

model = models.SequenceDNN_Regression.load("model.arch.json", "model.weights.h5")

f = open("../../id_dict_gen/id_dict.txt", 'r')
id_to_seq = json.loads(f.readlines()[0])

experiments = [("minP", "HepG2"), ("minP", "K562"), ("SV40P", "HepG2"), ("SV40P", "K562")]
ism = {}
base_to_row = {'A': 0, 'T': 1, 'C': 2, 'G': 3}

for name in id_to_seq.keys()[0:2000]:
    print name
    sequence, coords = str(id_to_seq[name][0]), id_to_seq[name][1]
    chrom, start, end = str(coords[0]), int(coords[1]), int(coords[2])
    for i in xrange(-4, 0):
        model_input = np.zeros((1, 1, 4, 145))
        subseq = bases(chrom, start + (i * 29), start + (i * 29) + 145).upper().replace("N", "A")
        for j in xrange(145):
            model_input[0][0][base_to_row[subseq[j]]][j] = 1
        ISM = model.in_silico_mutagenesis(model_input)
        for j in xrange(145):
            # we are looking at position: start + (i * 29) + j
            if (start + (i * 29) + j) in range(start, end):
                for k in xrange(len(experiments)):
                    if (experiments[k], chrom, (start + (i * 29) + j)) in ism:
                        ism[(experiments[k], chrom, (start + (i * 29) + j))].append(np.amax(np.abs(ISM[k][0][0][:,j])))
                    else:
                        ism[(experiments[k], chrom, (start + (i * 29) + j))] = [np.amax(np.abs(ISM[k][0][0][:,j]))]
    for i in xrange(0, 6):
        model_input = np.zeros((1, 1, 4, 145))
        subseq = sequence[(i * 29) : (i * 29) + 145].upper().replace("N", "A")
        for j in xrange(145):
            model_input[0][0][base_to_row[subseq[j]]][j] = 1
        ISM = model.in_silico_mutagenesis(model_input)
        for j in xrange(145):
            # we are looking at position: start + (i * 29) + j
            for k in xrange(len(experiments)):
                if (experiments[k], chrom, (start + (i * 29) + j)) in ism:
                    ism[(experiments[k], chrom, (start + (i * 29) + j))].append(np.amax(np.abs(ISM[k][0][0][:,j])))
                else:
                    ism[(experiments[k], chrom, (start + (i * 29) + j))] = [np.amax(np.abs(ISM[k][0][0][:,j]))]
    for i in xrange(6, 10):
        model_input = np.zeros((1, 1, 4, 145))
        subseq = bases(chrom, start + (i * 29), start + (i * 29) + 145).upper().replace("N", "A")
        for j in xrange(145):
            model_input[0][0][base_to_row[subseq[j]]][j] = 1
        ISM = model.in_silico_mutagenesis(model_input)
        for j in xrange(145):
            # we are looking at position: start + (i * 29) + j
            if (start + (i * 29) + j) in range(start, end):
                for k in xrange(len(experiments)):
                    if (experiments[k], chrom, (start + (i * 29) + j)) in ism:
                        ism[(experiments[k], chrom, (start + (i * 29) + j))].append(np.amax(np.abs(ISM[k][0][0][:,j])))
                    else:
                        ism[(experiments[k], chrom, (start + (i * 29) + j))] = [np.amax(np.abs(ISM[k][0][0][:,j]))]

import pickle
pickle.dump(ism, open("ism_test_2k.p", 'wb'))
