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
from in_silico_mutagenesis import in_silico_mutagenesis

model = models.SequenceDNN_Regression.load("models/models/145_weighted.arch.json", "models/models/145_weighted.weights.h5")

f = open("../../id_dict_gen/id_dict.txt", 'r')
id_to_seq = json.loads(f.readlines()[0])

experiments = [("minP", "HepG2"), ("minP", "K562"), ("SV40P", "HepG2"), ("SV40P", "K562")]

for name in id_to_seq.keys():
    sequence, coords = str(id_to_seq[name][0]), id_to_seq[name][1]
    chrom, start, end = str(coords[0]), int(coords[1]), int(coords[2])
    big_seq = bases(chrom, start - 72, end + 72).upper().replace('N', 'A')
    for i in xrange(295):
        middle_seq = big_seq[i : i + 145]
        ISM = in_silico_mutagenesis(model, middle_seq)
        for j in xrange(4):
            print '\t'.join(map(str, [chrom, start + i, start + i + 1, j, ISM[j], '+']))
