#!/home/lyw/software/anaconda3/bin/python

########################################################################
# Extract the percentage of reads mapping to MTBC from a kraken report #
########################################################################

import sys

kraken_report = sys.argv[1]
output = sys.argv[2]
sample_name = sys.argv[3]

f_out = open(output, 'w')
cov = 0

f_in = open(kraken_report, 'r')
i = 0
cat = ""
for lines in f_in:
    fields = lines.rstrip("\r\n").split("\t")
    if i == 10:
        cat = fields[5].strip()
    fields = lines.rstrip("\r\n").split("\t")
    if fields[5].find("Mycobacterium tuberculosis complex") != -1:
        cov += float(fields[0])
    i += 1
f_in.close()

f_out.write("%s\t%s\t%s\n" % (sample_name, cov, cat))
f_out.close()
