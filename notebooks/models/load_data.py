from models import SequenceDNN_Regression
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.model_selection import train_test_split
from collections import OrderedDict
from pprint import pprint
from warnings import warn
import numpy as np
import math

def get_seqs():
    key_to_seq = OrderedDict()
    seq_len = 145
    with open("../../data/Scaleup_counts_sequences/ScaleUpDesign1.sequences.txt") as f:
        for line in f:
            key, seq = line.strip().split()
            if "N" in seq:
                seq = seq.replace("N", "A")
            assert key not in key_to_seq
            key_to_seq[key] = seq

    with open("../../data/Scaleup_counts_sequences/ScaleUpDesign2.sequences.txt") as f:
        for line in f:
            key, seq = line.strip().split()
            if "N" in seq:
                seq = seq.replace("N", "A")
            assert key not in key_to_seq
            key_to_seq[key] = seq
    return key_to_seq

def get_expression():
    data = {}
    sample_weights = {}
    cell_types =  ["HepG2", "K562"]
    promoters = ["SV40P", "minP"]
    design_names = ["ScaleUpDesign1", "ScaleUpDesign2"]

    for cell_type in cell_types:
        for promoter in promoters:
            experiment_key = (cell_type, promoter)
            data[experiment_key] = {}
            sample_weights[experiment_key] = {}
        
            for design_name in design_names:
                with open("../../data/Scaleup_normalized/{}_{}_{}_mRNA_Rep1.normalized".format(cell_type, design_name, promoter)) as f:
                    for line in f:
                        parts = line.strip().split()
                        key = parts[0]
                        val = float(parts[1])
                        if parts[2] == "1":
                            data[experiment_key][key] = val
                        
                with open("../../data/Scaleup_normalized/{}_{}_{}_mRNA_Rep2.normalized".format(cell_type, design_name, promoter)) as f:
                    for line in f:
                        parts = line.strip().split()
                        key = parts[0]
                        val = float(parts[1])
                        if parts[2] == "1" and key in data[experiment_key]:
                            dot_prod = (val + data[experiment_key][key])
                            norm = math.sqrt(val**2 + data[experiment_key][key]**2)
                            cos = dot_prod/(math.sqrt(2) * norm)
                            sin = math.sqrt(1-cos**2)
                            w = sin * norm
                        
                            sample_weights[experiment_key][key] = w
                            data[experiment_key][key] = (val + data[experiment_key][key]) / 2.0
    return data, sample_weights
                        
# One hot encode DNA sequences the standard way.
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

def filter_scale(data, sample_weights):
    valid_keys = list(reduce(
        lambda acc, d: acc.intersection(d.keys()),
        data.values()[1:],
        set(data.values()[0].keys())
    ))

    scaler = StandardScaler()
    experiment_labels = []
    weights = []
    for experiment_key, key_to_normalized in data.items():
        filtered_normalized = np.array([key_to_normalized[key] for key in valid_keys]).reshape(-1, 1)
        filtered_weights = np.array([sample_weights[experiment_key][key] for key in valid_keys]).reshape(-1, 1)
        
        scaled = scaler.fit_transform(filtered_normalized)
    
        experiment_labels.append(scaled)
        weights.append(filtered_weights)

    y = np.hstack(experiment_labels)
    weights = np.hstack(weights).mean(axis=1) # What?

    scaler = MinMaxScaler(feature_range=(0,1))
    weights = scaler.fit_transform(-weights)
    return y, weights, valid_keys

def get_data():
    key_to_seq = get_seqs()
    data, sample_weights = get_expression()
    y, weights, valid_keys = filter_scale(data, sample_weights)
    X = seqs_to_encoded_matrix([key_to_seq[key] for key in valid_keys])

    X_temp, X_test, y_temp, y_test, weights_temp, weights_test = train_test_split(
        X, y, weights, test_size=0.2, random_state=42
    )
    X_train, X_valid, y_train, y_valid, weights_train, weights_valid = train_test_split(
        X_temp, y_temp, weights_temp, test_size=0.125, random_state=42
    )
    return X_train, X_valid, X_test, y_train, y_valid, y_test, weights_train, weights_valid, weights_test

