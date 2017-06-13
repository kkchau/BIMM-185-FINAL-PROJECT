import pymysql
import getpass
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.stats import gaussian_kde


def control(cur):
    """
        Return distributions
    """
    cur.execute(
        "SELECT SpliceScore FROM scored_denovo_db"
        + " WHERE PrimaryPhenotype='autism'"
        + " AND PathogenScore>=65;"
    )
    results = cur.fetchall()
    autism_scores = [r[0] for r in results]
    cur.execute(
        "SELECT SpliceScore FROM scored_denovo_db"
        + " WHERE PrimaryPhenotype='control'"
        + " AND PathogenScore<65;"
    )
    results = cur.fetchall()
    control_scores = [r[0] for r in results]

    return autism_scores, control_scores


if __name__ == '__main__':
    connection = pymysql.connect(
        host='localhost',
        user='root',
        db='bimm185',
        passwd=getpass.getpass("Input password: ")
    )

    aut, con = control(connection.cursor())

    # positive and control KDE distributions
    plt.clf()
    plt.xlabel("Risk Score")
    plt.ylabel("Frequency")
    plt.title("Distribution of Splice Probabilities")
    p_kde = gaussian_kde(aut)
    n_kde = gaussian_kde(con)
    x_val = np.arange(-50, 10, 0.01)
    plt.xlim((-50, 10))
    #plt.plot(x_val, p_kde(x_val))
    #plt.plot(x_val, n_kde(x_val))
    plt.hist(aut)
    plt.hist(con)
    plt.legend(["Autism", "Control"])
    plt.show()

