"""
    Functions for calculations of the inference model
"""


import pymysql
import getpass
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde


def control(cur, positive=True):
    """
        Takes in a pymysql connection cursor
        Read MySQL database and calculate distribution of values for control
        Return: array of splice score
                array of pathogen*risk score for each period
    """
    if positive:
        phen = 'autism'
    else:
        phen = 'control'
    cur.execute(
        "SELECT SpliceScore,PathogenScore FROM scored_denovo_db"
        + " WHERE PrimaryPhenotype='{}';".format(phen)
    )
    result = cur.fetchall()
    splice_scores = [spl[0] / 100 for spl in result]
    pathog_scores = [pat[1] / 100 for pat in result]
    pathog_per_period = []
    for per in range(8):
        cur.execute(
            "SELECT P{}Risk FROM scored_denovo_db".format(per + 1)
            + " WHERE PrimaryPhenotype='{}';".format(phen)
        )
        p_result = cur.fetchall()
        pathog_per_period.append(
            [p*g for p,g in zip(
                pathog_scores, [risk[0] for risk in p_result]
            )]
        )

    return splice_scores, pathog_per_period


if __name__ == '__main__':
    connection = pymysql.connect(
        host='localhost',
        user='root',
        db='bimm185',
        passwd=getpass.getpass("Input password: ")
    )
    cursor = connection.cursor()
    pos_splice, pos_pathog = control(cursor)
    neg_splice, neg_pathog = control(cursor, positive=False)

    # positive
    pos_fig = plt.figure()
    pos_fig.suptitle("Splice Probability vs. Pathogenicity in ASD SNVs")

    # negative
    neg_fig = plt.figure()
    neg_fig.suptitle("Splice Probability vs. Pathogenicity in Control SNVs")

    # scatter plots for splice score vs pathogenicity
    for i in range(8):

        # positive
        ax = pos_fig.add_subplot(241 + i)
        # ax.set_xlabel("Splice Score")
        # ax.set_ylabel("Pathogenicity Score")
        ax.scatter(pos_splice, pos_pathog[i], s=2)

        # negative
        an = neg_fig.add_subplot(241 + i)
        # an.set_xlabel("Splice Score")
        # an.set_ylabel("Pathogenicity Score")
        an.scatter(neg_splice, neg_pathog[i], s=2)


    pos_fig.tight_layout()
    neg_fig.tight_layout()

    pos_fig.subplots_adjust(top=0.9)
    neg_fig.subplots_adjust(top=0.9)

    pos_fig.savefig("pos_corr.png")
    neg_fig.savefig("neg_corr.png")

    # positive and control KDE distributions
    plt.clf()
    plt.xlabel("Splice Score")
    plt.ylabel("Frequency")
    plt.title("Distribution of Splice Probabilities")
    p_kde = gaussian_kde(pos_splice)
    n_kde = gaussian_kde(neg_splice)
    x_val = np.arange(-1, 0.5, 0.01)
    plt.xlim((-1, 0.5))
    plt.plot(x_val, p_kde(x_val))
    plt.plot(x_val, n_kde(x_val))
    plt.legend(["Autism", "Control"])
    plt.savefig("splice_distributions")
    plt.show()
