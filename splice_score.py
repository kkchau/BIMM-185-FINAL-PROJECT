"""
    Score mutations for splicing
"""


import subprocess


# input and output files
with open('filtered_denovo_db.tsv', 'r') as db:
    with open('splice_score_denovo.tsv', 'w') as splice:

        # iterate through input file
        for line in db:
            line = line.strip().split()

            # relevant variables
            chromosome = line[3]
            position = line[4]
            ref, var = line[5].strip().split('>')

            # get splicing scores
            output = subprocess.check_output([
                'tabix',
                'spidex_public_noncommercial_v1.0/spidex_public_noncommercial_v1_0.tab.gz',
                'chr{}:{}-{}'.format(chromosome, position, position)
            ])

            variant = None

            # get the variant if exists
            if output:
                output = output.splitlines()
                for record in output:
                    record = record.strip().split()
                    if len(record) == 0:
                        continue
                    if record[2] == ref and record[3] == var:
                        variant = record
                        break

            # if the variant exists, record the record with splice score
            if variant:
                line.append(variant[4])
                splice.write('\t'.join(line) + '\n')
