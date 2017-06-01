"""
    Collapse transcripts together that have same gene/variant
"""


found_variants = []


with open('splice_score_denovo.tsv', 'r') as uncollapsed:
    with open('collapsed_splice_denovo.tsv', 'w') as collapsed:
        for line in uncollapsed:
            line = line.strip().split()
            if (line[0], line[1], line[4], line[5]) not in found_variants:
                collapsed.write('\t'.join(line) + '\n')
                found_variants.append((line[0], line[1], line[4], line[5]))

