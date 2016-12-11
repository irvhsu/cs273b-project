from warnings import warn
import numpy as np
from collections import OrderedDict

class MrpaData:
    cell_types =  ['HepG2', 'K562']
    promoters = ['SV40P', 'minP']
    design_names = ['ScaleUpDesign1', 'ScaleUpDesign2']
    bases = ['A', 'T', 'C', 'G']
    
    def __init__(self):
        self.split_data = self._load_data()
        self.data = self._merge_data()
        self.valid_keys = self._get_valid_keys()
        self.seqs = self._get_seqs()
        self.one_hot_seqs = self._one_hot_encode_seqs()
        
    def y_multitask(self):
        """
        Returns a N x 4 np.array of experimental data averaged between reps
        Data follows the order given by self.valid keys
        row 0: cell_type[0], promoters[0]
        row 1: cell_type[0], promoters[1]
        ...
        """
        return np.array([
                [self.data[experiment_key][key] for key in self.valid_keys]
                for experiment_key in self._experiment_keys()
                ]).T

    def y_merged_promoters(self):
        """
        Returns a N x 2 np.array of experimental data averaged accross reps and
        promoters
        """
        vector = self.y_multitask()
        hep_g2 = (vector[0, :] + vector[1, :]) / 2.0
        k562   = (vector[2, :] + vector[3, :]) / 2.0
        return np.array([hep_g2, k562]).T

    def X_one_hot(self):
        return self.one_hot_seqs

    def get_data(self):
        """
        Returns raw data in a dictionary of (cell_type, promoter) -> mean score
        """
        return {key: val for key, val in self.data.items()}
        
    def get_rep_data(self):
        """
        Returns raw data in a dictionary of (cell_type, promoter) -> (rep1, rep2)
        """
        return {key: val for key, val in self.split_data.items()}        

    def _experiment_keys(self):
        return [(cell_type, promoter) for cell_type in self.cell_types for promoter in self.promoters]

    def _load_data(self):
        split_data = OrderedDict()
        for cell_type in self.cell_types:
            for promoter in self.promoters:
                experiment_key = (cell_type, promoter)
                split_data[experiment_key] = {}
                for design_name in self.design_names:
                    with open("../data/Scaleup_normalized/{}_{}_{}_mRNA_Rep1.normalized".format(cell_type, design_name, promoter)) as f:
                        for line in f:
                            parts = line.strip().split()
                            key = parts[0]
                            val = float(parts[1])
                            if parts[2] == "1":
                                assert key not in split_data[experiment_key]
                                split_data[experiment_key][key] = (val, 0)

                    with open("../data/Scaleup_normalized/{}_{}_{}_mRNA_Rep2.normalized".format(cell_type, design_name, promoter)) as f:
                        for line in f:
                            parts = line.strip().split()
                            key = parts[0]
                            val = float(parts[1])
                            if parts[2] == "1" and key in split_data[experiment_key]:
                                assert split_data[experiment_key][key][1] == 0
                                split_data[experiment_key][key] = (split_data[experiment_key][key][0], val)
        return split_data

    def _merge_data(self):
        data = OrderedDict()
        for experiment_key, key_to_split in self.split_data.items():
            data[experiment_key] = {}
            for key, val in key_to_split.items():
                data[experiment_key][key] = sum(list(self.split_data[experiment_key][key])) / 2.0
        return data

    def _get_valid_keys(self):
        return list(reduce(
                lambda acc, d: acc.intersection(d.keys()),
                self.data.values()[1:], 
                set(self.data.values()[0].keys())
                ))

    def _get_seqs(self):
        key_to_seq = {}

        with open("../data/Scaleup_counts_sequences/ScaleUpDesign1.sequences.txt") as f:
            for line in f:
                key, seq = line.strip().split()
                if "N" in seq:
                    warn("Replacing 'N' bases in seq with 'A' in seq {}.".format(seq))
                    seq = seq.replace("N", "A")
                assert key not in key_to_seq
                key_to_seq[key] = seq
        
        with open("../data/Scaleup_counts_sequences/ScaleUpDesign2.sequences.txt") as f:
            for line in f:
                key, seq = line.strip().split()
                if "N" in seq:
                    warn("Replacing 'N' bases in seq with 'A' in seq {}.".format(seq))
                    seq = seq.replace("N", "A")
                assert key not in key_to_seq
                key_to_seq[key] = seq
        return key_to_seq
    
    def _one_hot_encode_seq(self, seq):
        result = np.zeros((len(self.bases), len(seq)))
        for i, base in enumerate(seq):
            if base == 'N': base = 'A'
            result[self.bases.index(base), i] = 1
        return result

    def _one_hot_encode_seqs(self):
        shape = (len(self.valid_keys), 1, len(self.bases), len(self.seqs.values()[0]))
        return np.concatenate([self._one_hot_encode_seq(self.seqs[key]) for key in self.valid_keys]).reshape(*shape)
