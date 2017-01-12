from __future__ import absolute_import, division, print_function
import numpy as np
from collections import OrderedDict
from sklearn.metrics import auc, log_loss, precision_recall_curve, roc_auc_score, mean_squared_error, mean_absolute_error, median_absolute_error, r2_score
from prg.prg import create_prg_curve, calc_auprg


def loss(labels, predictions):
    return log_loss(labels, predictions)


def positive_accuracy(labels, predictions, threshold=0.5):
    return 100 * (predictions[labels] > threshold).mean()


def negative_accuracy(labels, predictions, threshold=0.5):
    return 100 * (predictions[~labels] < threshold).mean()


def balanced_accuracy(labels, predictions, threshold=0.5):
    return (positive_accuracy(labels, predictions, threshold) +
            negative_accuracy(labels, predictions, threshold)) / 2


def auROC(labels, predictions):
    return roc_auc_score(labels, predictions)


def auPRC(labels, predictions):
    precision, recall = precision_recall_curve(labels, predictions)[:2]
    return auc(recall, precision)


def auPRG(labels, predictions):
    return calc_auprg(create_prg_curve(labels, predictions))


def recall_at_precision_threshold(labels, predictions, precision_threshold):
    precision, recall = precision_recall_curve(labels, predictions)[:2]
    return 100 * recall[np.searchsorted(precision - precision_threshold, 0)]

class RegressionResult(object):

    def __init__(self, labels, predictions, sample_weight=None, task_names=None):
        self.results = [OrderedDict((
            ('Mean Squared Error', mean_squared_error(task_labels, task_predictions, sample_weight=sample_weight)),
            ('Mean Absolute Error', mean_absolute_error(task_labels, task_predictions, sample_weight=sample_weight)),
            ('Median Absolute Error', median_absolute_error(task_labels, task_predictions)),
            ('R2 Score', r2_score(task_labels, task_predictions, sample_weight=sample_weight)),
        )) for task_labels, task_predictions in zip(labels.T, predictions.T)]
        self.task_names = task_names
        self.multitask = labels.shape[1] > 1

    def __str__(self):
        return '\n'.join(
            '{}Mean Squared Error: {:.4f}\tMean Absolute Error: {:.4f}\t '
            'Median Absolute Error: {:.4f}\t R2 Score: {:.4f}'.format(
                '{}: '.format('Task {}'.format(
                    self.task_names[task_index]
                    if self.task_names is not None else task_index))
                if self.multitask else '', *results.values())
            for task_index, results in enumerate(self.results))

    def __getitem__(self, item):
        return np.array([task_results[item] for task_results in self.results])
