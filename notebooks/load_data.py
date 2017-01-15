import numpy as np
from math import log
from sklearn.model_selection import train_test_split

data_dir = '~/cs273b-project/data/Scaleup_counts_sequences'
promoters = ['minP', 'SV40P']
cell_types = ['HEPG2', 'K562']
designs = ['1', '2']
bases = ['A', 'T', 'C', 'G']

def normalized_scores(dna, rna1, rna2, labels):
    rep1, rep2, dna_count = [], [], []
    total_r_rep1, total_r_rep2, total_d = 0, 0, 0
    for i, lines in enumerate(zip(dna, rna1, rna2)[1:]):
        d, r1, r2 = lines
        d_name, d_val   =  d.strip().split()
        r1_name, r1_val = r1.strip().split()
        r2_name, r2_val = r2.strip().split()
        d_val, r1_val, r2_val = map(float, [d_val, r1_val, r2_val])
        
        assert d_name == r1_name == r2_name == labels[i]
        
        if d_val < 20:
            rep1      += [0]
            rep2      += [0]
            dna_count += [0]
        else:
            rep1 += [log(r1_val + 1, 2) - log(d_val + 1, 2)]
            rep2 += [log(r2_val + 1, 2) - log(d_val + 1, 2)]
            dna_count += [d_val]
            
            total_r_rep1 += r1_val + 1
            total_r_rep2 += r2_val + 1
            total_d += d_val + 1

    rep1 = map(lambda x: x + log(total_d, 2) - log(total_r_rep1, 2), rep1)
    rep2 = map(lambda x: x + log(total_d, 2) - log(total_r_rep2, 2), rep2)
    return rep1, rep2, dna_count

def get_weights(dna_count, rep1, rep2):
    # Fit transform and regressor on valid data
    valid_avg       = [(r1+r2)    / 2 for r1, r2, dna in zip(rep1, rep2, dna_count) if dna > 19]
    valid_variance  = [(r1-r2)**2 / 2 for r1, r2, dna in zip(rep1, rep2, dna_count) if dna > 19]
    valid_dna_count = [dna for dna in dna_count if dna > 19]
    
    scaler = StandardScaler()
    valid_X = scaler.fit_transform(np.array([valid_dna_count, valid_avg]).T)
    regressor = KNeighborsRegressor(n_neighbors = 200).fit(valid_X, valid_variance)

    # Run on complete dataset and then mask nonvalid values
    avg = [(r1+r2) / 2 for r1, r2 in zip(rep1, rep2)]
    X = scaler.transform(np.array([dna_count, avg]).T)
    
    return [1 / var if valid else 0 for var, valid in zip(regressor.predict(X), dna_count)]

def get_labels():
    labels = {}
    for design in designs:
        with open("{}/DNACOUNTS/ScaleUpDesign{}_{}_Plasmid.counts".format(data_dir, design, promoters[0])) as _dna:
            _dna.readline()
            labels[design] = [line.strip().split()[0] for line in _dna]

def get_activities():
    y = []
    w = []
    for promoter in promoters:
        for cell_type in cell_types:
            merged_y, merged_w = [], []
            for design in ['1', '2']:
                dna = open("{}/DNACOUNTS/ScaleUpDesign{}_{}_Plasmid.counts".format(data_dir, design, promoter))
                rna1 = open("{}/{}/{}_ScaleUpDesign{}_{}_mRNA_Rep1.counts".format(
                    data_dir, cell_type, cell_type, design, promoter))
                rna2 = open("{}/{}/{}_ScaleUpDesign{}_{}_mRNA_Rep2.counts".format(
                    data_dir, cell_type, cell_type, design, promoter))
                rep1, rep2, dna_count = normalized_scores(dna, rna1, rna2, labels[design])
                
                merged_y += [(r1 + r2) / 2 for r1, r2 in zip(rep1, rep2)]
                merged_w += get_weights(dna_count, rep1, rep2)
        
                dna.close()
                rna1.close()
                rna2.close()
            scaler = StandardScaler().fit([val for val, weight in zip(merged_y, merged_w) if weight != 0])
            y += [list(scaler.transform(np.array(merged_y)))]
            w += [merged_w]
    y = np.array(y).T
    w = np.array(w).T
    return y, w
        
##############################################################################
########  Get sequences from files ##########################################

def get_seqs(labels):
    seqs = []
    for design in designs:
        with open("../data/Scaleup_counts_sequences/ScaleUpDesign{}.sequences.txt".format(design)) as f:
            for i, line in enumerate(f):
                key, seq = line.strip().split()
                assert key == labels[design][i]
                if 'N' in seq: seq = seq.replace('N', 'A')
                seqs += [seq]
    return seqs

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
    return result

def get_one_hot_seqs(labels):
    return seqs_to_encoded_matrix(get_seqs(labels))

##########################################################################
######## The whole process ##############################################

def get_data():
    labels = get_labels()
    activities = get_activities(labels)
    dna_counts = get_dna_counts(labels)
    seqs = get_one_hot_seqs(labels)

    return labels, activities, seqs, dna_counts
