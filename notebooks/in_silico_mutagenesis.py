bases = ['A', 'C', 'G', 'T']
length = 145
middle = length / 2

def in_silico_mutagenesis(model, seq):
    assert len(seq) == length
    model_input = np.zeros((4, 1, 4, length))
    for i in range(4):
        model_input[i, 0, i, middle] = 1
        for j in range(middle) + range(middle+1, length):
            model_input[i, 0, bases.index(seq[j]), j] = 1

    activities = model.predict(model_input) # first axis is bases, second is experiment
    out = []
    baseline = np.array([activities[bases.index(seq[middle]), :] for i in range(4)])
    deltas = activities - baseline
    for i in range(4):
        d = deltas[:, i]
        out += [max(d)] if abs(max(d)) > abs(min(d)) else [min(d)]
    return out
