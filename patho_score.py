"""
    Add pathogenicity scores to list of variants
    From umd_predictions.tsv --> splice_score_denovo.tsv
    Averages pathogenicity scores for all transcripts derived from the variant
"""


from collections import defaultdict


# map prediction scores to variants
pred_scores = defaultdict(list)
with open('umd_predictions.tsv', 'r') as pred:
    next(pred)      # skip header
    for line in pred:
        line = line.strip().split()
        identifier = tuple(line[0:4])
        score = line[-2] if 'Probab' not in line[-2] else line[-3]
        if score == 'NA':
            score = 0
        pred_scores[identifier].append(float(score))

# read splice_scored variants
with open('splice_score_denovo.tsv', 'r') as splice:

    # output file
    with open('pathogen_score.tsv', 'w') as patho:

        for record in splice:
            record = record.strip().split()
            ref, var = record[5].strip().split('>')
            variant = (
                'chr{}'.format(record[3]), record[4], ref, var
            )
            if variant in pred_scores:
                record.append(str(
                    sum(pred_scores[variant]) / len(pred_scores[variant])
                ))
            else:
                record.append('0')

            patho.write('\t'.join(record) + '\n')

