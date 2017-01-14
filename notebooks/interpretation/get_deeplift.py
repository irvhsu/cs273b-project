# outputs dict "deep", keys: ("HepG2", "minP"), ("K562", "minP"), ("HepG2", "SV40P"), ("K562", "SV40P"), values chroms
# chroms is a dict keys: chrom (string), values positions
# positions is a dict keys: position (int), values deeplift score

import sys
import json
import numpy as np
from dragonn import models

begin = int(sys.argv[1])
end = int(sys.argv[2])
out = open(sys.argv[3], 'w')

model = models.SequenceDNN_Regression.load("model.arch.json", "model.weights.h5")

f = open("../../id_dict_gen/id_dict.txt", 'r')
id_to_seq = json.loads(f.readlines()[0])
f.close()

base_to_row = {'A': 0, 'T': 1, 'C': 2, 'G': 3}

for name in id_to_seq.keys()[int(sys.argv[1]):int(sys.argv[2])]:
    sequence, coords = str(id_to_seq[name][0]), id_to_seq[name][1]
    chrom, start, end = str(coords[0]), coords[1], coords[2]
    for i in xrange(31):
        model_input = np.zeros((1, 1, 4, 145))
        subseq = sequence[5 * i : 145 + 5 * i]
        for j in xrange(145):
            # coordinates chrom, start + 5 * i + j
            model_input[0][0][base_to_row[subseq[j]]][j] = 1
        D = model.deeplift(model_input)
        entry = []
        for task in range(4):
            for pos in range(145):
                scores = D[task][0][0][:, pos]
                if max(scores) != 0:
                    entry += [max(scores)]
                else:
                    entry += [min(scores)]
        out.write(','.join([chrom, str(start + 5 * i)] + map(str, entry)) + '\n')
out.close()
