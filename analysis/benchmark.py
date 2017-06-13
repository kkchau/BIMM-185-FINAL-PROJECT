"""
    Functions to provide benchmarking for the given predictive model
"""


import pymysql
import getpass
import matplotlib.pyplot as plt
import numpy as np
from inference_model import control
from inference_model import posterior


def bench(cur, threshold):
    """
        Given the specified threshold, for each variant, determine whether
        the model gives a true or false positive or negative
        Returns counts of 'tp', 'tn', 'fp', 'fn'
    """
    confusion_dist = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0}

    # posterior function
    pos_splice, _ = control(cur)
    neg_splice, _ = control(cur, positive=False)
    post_func = posterior(pos_splice, neg_splice)

    # get all variants and scores into memory
    cur.execute(
        "SELECT PrimaryPhenotype,SpliceScore FROM scored_denovo_db_alt;"
    )
    
    # count values of the confusion matrix
    for record in cur.fetchall():
        prob = post_func(float(record[1]))
        if record[0] == 'autism':
            if prob < threshold:
                confusion_dist['fn'] += 1
            else:
                confusion_dist['tp'] += 1
        else:
            if prob < threshold:
                confusion_dist['tn'] += 1
            else:
                confusion_dist['fp'] += 1

    return confusion_dist


def bench_thresh(cur):
    """
        Get a set of benchmarking values
        Returns a list of dictionaries
    """
    ret_array = []
    for thresh in np.arange(0, 1, 0.05):
        print(thresh)
        ret_array.append(bench(cur, thresh)) 

    return ret_array


def svs(bench_vals):
    """
        Sensitivity and Specificity value calculations
        Returns an array for sensitivity and specificity each per threshold
    """
    sens = []
    spec = []

    for vals in bench_vals:
        if vals['tp'] > 0 or vals['fp'] > 0:
            sens.append(
                vals['tp'] / (vals['tp'] + vals['fn'])
            )
        else:
            sens.append(0)
        if vals['tn'] > 0 or vals['fn'] > 0:
            spec.append(
                vals['tn'] / (vals['tn'] + vals['fp'])
            )
        else:
            spec.append(0)

    return sens, spec


def roc(bench_vals):
    """
    Create receiver operator curve values from an array of benchmark values
    Uses table: tus
    Returns two arrays: True Positive Rate and False Positive Rate
    """
    tpr = []
    fpr = []
    for value in bench_vals:
        if value['tp'] + value['fn'] != 0:
            tpr.append(value['tp'] / (value['tp'] + value['fn']))
        else:
            tpr.append(0)
        if value['tn'] + value['fp'] != 0:
            fpr.append(value['fp'] / (value['tn'] + value['fp']))
        else:
            fpr.append(0)

    return tpr, fpr


def acc(bench_vals):
    """
    Return an array of accuracy values for each set of benchmark values
    """
    return [
        (val['tp'] + val['tn']) / (
            val['tp'] + val['tn'] + val['fp'] + val['fn']
        ) for val in bench_vals
    ]


if __name__ == '__main__':

    connection = pymysql.connect(
        host='localhost',
        user='root',
        db='bimm185',
        passwd=getpass.getpass("Input password: ")
    )

    print("Benchmarking")
    benchmarking_values = bench_thresh(connection.cursor())

    print("SVS")
    sensitivity, specificity = svs(benchmarking_values)
    print("ROC")
    true_pos, false_pos = roc(benchmarking_values)
    print("ACC")
    accuracy = acc(benchmarking_values)

    thresh_array = np.arange(0, 1, 0.05)
    plt.xlim([0, 1])
    svs_fig = plt.figure()
    svs_fig.suptitle("Sensitivity vs Specificity")
    ax = svs_fig.add_subplot(111)
    ax.set_xlabel("Posterior Threshold")
    ax.set_ylabel("Performance")
    ax.plot(thresh_array, sensitivity)
    ax.plot(thresh_array, specificity)
    ax.legend(['Sensitivity', 'Specificity'])
    svs_fig.savefig("svs.png")

    roc_fig = plt.figure()
    roc_fig.suptitle("Receiver Operator Curve")
    ax = roc_fig.add_subplot(111)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.plot(false_pos, true_pos)
    roc_fig.savefig("roc.png")

    acc_fig = plt.figure()
    acc_fig.suptitle("Accuracy per Threshold Value")
    ax = acc_fig.add_subplot(111)
    ax.set_xlabel("Posterior Threshold")
    ax.set_ylabel("Accuracy")
    ax.plot(thresh_array, accuracy)
    acc_fig.savefig("acc.png")

    connection.close()
