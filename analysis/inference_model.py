"""
    Functions for calculations of the inference model
"""


import pymysql
import getpass
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from scipy.stats import linregress


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
        "SELECT SpliceScore,PathogenScore FROM scored_denovo_db_alt"
        + " WHERE PrimaryPhenotype='{}';".format(phen)
    )
    all_result = cur.fetchall()
    result = []
    for index in np.random.randint(0, len(all_result), len(all_result)//2):
        result.append(all_result[index])
    splice_scores = [spl[0] for spl in result]
    pathog_scores = [pat[1] for pat in result]
    pathog_per_period = []
    for per in range(8):
        cur.execute(
            "SELECT P{}Risk FROM scored_denovo_db_alt".format(per + 1)
            + " WHERE PrimaryPhenotype='{}';".format(phen)
        )
        p_result = cur.fetchall()
        pathog_per_period.append(
            [p for p,g in zip(
                pathog_scores, [risk[0] for risk in p_result]
            )]
        )

    return splice_scores, pathog_per_period


def posterior(pos_control, neg_control):
    p_kde = gaussian_kde(pos_control)
    #p_kde.covariance_factor = lambda: 0.2
    #p_kde._compute_covariance()
    n_kde = gaussian_kde(neg_control)
    #n_kde.covariance_factor = lambda: 0.2
    #n_kde._compute_covariance()
    
    def _probability(x):
        return(p_kde(x)/(p_kde(x) + n_kde(x)))

    return _probability

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

    """
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
        slope, intercept, rval, pval, stderr = linregress(pos_splice, pos_pathog[i])
        ax.scatter(pos_splice, pos_pathog[i], s=2)

        # negative
        an = neg_fig.add_subplot(241 + i)
        an.scatter(neg_splice, neg_pathog[i], s=2)


    pos_fig.tight_layout()
    neg_fig.tight_layout()

    pos_fig.subplots_adjust(top=0.9)
    neg_fig.subplots_adjust(top=0.9)

    #pos_fig.savefig("pos_corr.png")
    #neg_fig.savefig("neg_corr.png")
    """

    # positive and control KDE distributions
    plt.clf()
    plt.xlabel("Splice Score")
    plt.ylabel("Frequency")
    plt.title("Distribution of Splice Probabilities")
    p_kde = gaussian_kde(pos_splice)
    p_kde.covariance_factor = lambda: 0.1
    p_kde._compute_covariance()
    n_kde = gaussian_kde(neg_splice)
    n_kde.covariance_factor = lambda: 0.1
    n_kde._compute_covariance()
    x_val = np.arange(-2, 0.5, 0.05)
    plt.xlim((-2, 0.5))
    plt.plot(x_val, p_kde(x_val))
    plt.plot(x_val, n_kde(x_val))
    plt.legend(["Autism", "Control"])
    plt.savefig("splice_distributions")
    plt.show()

    # posterior
    plt.clf()
    plt.title("Unweighted Posterior Probability")
    plt.xlabel("Splice Score")
    plt.ylabel("Probability")
    post = posterior(pos_splice, neg_splice)
    plt.plot(np.arange(-2, 0.5, 0.01), post(np.arange(-2, 0.5, 0.01)))
    plt.savefig("posterior.png")
    plt.show()
