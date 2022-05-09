#######################################################################
# Output the number of different SNP genotype between two MTB samples #
#######################################################################

import sys
import re
import numpy as np
import pandas as pd


def compare_two_seq(seq1, seq2):
    return sum([1 if i != j else 0 for i, j in zip(seq1, seq2)])


out_dir = sys.argv[1]
data_dir = sys.argv[2]
# ploidy = 1 or 2
ploidy = sys.argv[3]

# input lineage representative list
try:
    with open(data_dir + "/lineage_sample_list.txt", 'r') as f:
        lineage_samples = f.read().strip().split('\n')
except IOError:
    lineage_samples = []
################################################################

f = open(out_dir + "/all_samples_snp_genotype.txt", 'r')

buf = f.readline().strip()  # header
header = buf.split('\t')
sample_num = len(header) - 4 + 1  # exclude header column CHROM, POS, REF, ALT, samples include REF
seq_list = [''] * sample_num
seq_name = []
for i in range(sample_num):
    if i == 0:
        seq_name = ['REF']
    else:
        # if a strain sample is in lineage representative strain list
        # convert the strain sample name format to something like "L1.1.1_SRR5067315"
        # this name format will be better viewed in a phylogenetic tree
        seqname = re.findall(r'\[[0-9]+](.*):GT', header[i+3])
        lineage_sample_matched = [i for i in lineage_samples if seqname[0] in i]
        if lineage_sample_matched:
            seq_name += [i.replace('lineage_representative_strain/', '').replace('/', '_')
                         for i in lineage_sample_matched]
        else:
            seq_name += seqname

buf = f.readline().strip()
while buf:
    items = buf.split('\t')
    for i in range(sample_num):
        if i == 0:
            seq_list[i] += items[2]
        else:
            if ploidy == "1":
                seq_list[i] += items[i+3]
            elif ploidy == "2":
                hom_alt_list = ["%s/%s" % (alt, alt) for alt in items[3].split(',')]
                if items[i+3] in hom_alt_list:
                    seq_list[i] += (items[3])[0]
                else:
                    seq_list[i] += items[2]
    buf = f.readline().strip()
f.close()

with open(out_dir + "/all_samples_snp_genotype.fasta", 'w') as o:
    print(seq_name, file=sys.stderr)
    for i in range(sample_num):
        o.write('>' + seq_name[i] + '\n')
        o.write(seq_list[i] + '\n')

# output the number of different SNP between two samples
output_matrix = []
for i in range(sample_num):
    diff = []
    for j in range(sample_num):
        diff.append(compare_two_seq(seq_list[i], seq_list[j]))
    output_matrix.append(diff)

snp_diff_df = pd.DataFrame(np.array(output_matrix, dtype=int), index=seq_name, columns=seq_name)
print('\n' * 3, file=sys.stderr)
print("### Matrix - Different SNP numbers between each sample pair ###\n", file=sys.stderr)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
print(snp_diff_df.to_string(header=False), file=sys.stderr)
print('\n' * 3, file=sys.stderr)
snp_diff_df.to_csv(out_dir + "/all_samples_snp_diff.txt", sep='\t')
