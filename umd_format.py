"""
    Convert variants to uploadable file to UMD-predict
"""


written_records = []


with open('splice_score_denovo.tsv', 'r') as collapsed:
    with open('umd_upload.tsv', 'w') as upload:
        for line in collapsed:
            line = line.strip().split()
            ref, var = line[5].strip().split('>')

            if (line[3], line[4], ref, var) in written_records:
                continue
            else:
                upload.write('\t'.join([
                    "chr{}".format(line[3]),
                    line[4],
                    ref,
                    var
                ]))
                upload.write('\n')

                written_records.append((line[3], line[4], ref, var))
