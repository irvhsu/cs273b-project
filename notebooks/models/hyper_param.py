from load_data import get_data
from sklearn.metrics import mean_squared_error
from models import SequenceDNN_Regression

X_train, X_valid, X_test, y_train, y_valid, y_test, w_train, w_valid, w_test = get_data()
conv1_n = (100,)
conv2_n = (100,)
conv1_w = (8,)
conv2_w = (15,)

for n1 in conv1_n:
    for n2 in conv2_n:
        for w1 in conv1_w:
            for w2 in conv2_w:
                name = "{}n1_{}n2_{}w1_{}w2_act".format(n1, n2, w1, w2)
                fn = "models/{}".format(name)
                try:
                    model = SequenceDNN_Regression.load(fn + '.arch.json', fn + '.weights.h5')
                except IOError:
                    model = SequenceDNN_Regression(
                        seq_length=X_train.shape[3],
                        num_filters=(n1, n2),
                        conv_width=(w1, w2),
                        pool_width=X_train.shape[3] - w1 - w2,
                        num_tasks=y_train.shape[1],
                        dropout=0.1
                    )
                    SequenceDNN_Regression.save(model, fn)
                    model.train(X_train, y_train, (X_valid, y_valid),
                                train_sample_weight=w_train, valid_sample_weight=w_valid)

                model.plot_architecture(fn + '.png')
                SequenceDNN_Regression.save(model, fn)

                print name
                print mean_squared_error(y_train, model.predict(X_train), sample_weight = w_train)
                print mean_squared_error(y_valid, model.predict(X_valid), sample_weight = w_valid)
