"""
    Append the risk for each variant
    Risk is calculated as num_coex_partners (per period) / max_partners
"""


# map gene to risk per period
periods = [{} for _ in range(8)]

for period in range(1,9):
    with open('gene_networks/bSpan_{}.txt'.format(period), 'r') as network:
        all_values = []
        for record in network:
            record = record.strip().split()
            periods[period - 1][record[0]] = float(record[1])
            all_values.append(int(record[1]))
        max_partners = max(all_values)
        for gene in periods[period - 1]:
            periods[period - 1][gene] = periods[period - 1][gene] / max_partners

# append values to pathogen file
with open('pathogen_score.tsv', 'r') as pathogens:
    with open('risk_score.tsv', 'w') as risks:

        for line in pathogens:
            line = line.strip().split()
            for i in range(8):
                if line[1] in periods[i]:
                    line.append(periods[i][line[1]])
                else:
                    line.append(0)

            risks.write('\t'.join([str(x) for x in line]) + '\n')
