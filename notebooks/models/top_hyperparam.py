from load_data import get_data
from sklearn.metrics import mean_squared_error
from models import Basset, SequenceDNN_Regression

X_train, X_valid, X_test, y_train, y_valid, y_test, w_train, w_valid, w_test = get_data()

# Try deeper top of network

models = [Basset(
    seq_length = X_train.shape[3],
    num_filters = (100, 100, 100),
    conv_width  = (8, 15, 8),
    pool_width  = (1, 4, 24),
    num_task    = 4,
    dropout = 0.1,
    L1 = (0, 0, 0)),
          Basset(
    seq_length  = X_train.shape[3],
    num_filters = (100, 100, 100),
    conv_width  = (8, 15, 8),
    pool_width = (1, 24, 4),
    num_task   = 4,
    dropout = 0.1,
    L1 = (0, 0, 0))               
    ]

for i, model in enumerate(models):
    if i:
        name = 'thin'
    else:
        name = 'wide'
    fn = 'models/' + name
    model.train(X_train, y_train, (X_valid, y_valid),
                train_sample_weight=w_train, valid_sample_weight=w_valid)

    model.plot_architecture(fn + '.png')
    SequenceDNN_Regression.save(model, fn)

    print name
    print mean_squared_error(y_train, model.predict(X_train), sample_weight = w_train)
    print mean_squared_error(y_valid, model.predict(X_valid), sample_weight = w_valid)

