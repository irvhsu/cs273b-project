def one_hot(seq):
    

def main(model, fasta):
    # Create input vector and record chromsomal positions
    seqs, labels = []
    for line in fasta:
        if line[0] == '>':
            if sequence:
                seqs += one_hot(sequence)
            labels += line.strip()
            sequence = ''
        else:
            sequence += line.strip()

            X = np.array(seqs)
    # Predict
    model.predict(X)
    # Link predictions with labels
